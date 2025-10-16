import os
from pathlib import Path
import shutil
import json
from collections import Counter
from zipfile import ZipFile

import pandas as pd
from tqdm import tqdm

try:
    import pyarrow.parquet as pq  # schema introspection without loading data
except ImportError:
    pq = None

# ----------------------------
# User settings (edit freely)
# ----------------------------
PLINDER_ITERATION = os.environ.get("PLINDER_ITERATION", "v2")  # unused but kept for logging
OUTPUT_DIR        = Path("./plinder_patterns")                 # where to write patterns
MAX_PATTERNS      = None                                       # e.g. 50 to test, or None for all

# Dedup strategy:
#   "uniqueness"  -> group by 'uniqueness' from splits/split.parquet and keep best-quality system per group
#   "nonredundant"-> simply use index/annotation_table_nonredundant.parquet (skip grouping)
DEDUP_STRATEGY    = "uniqueness"

# Light QC thresholds (tweak as you like)
REQUIRE_SINGLE_LIGAND = True        # 1 ligand chain per system
MIN_INTERACTIONS      = 1           # >=1 proper PL interactions (if column exists)
MAX_RESOLUTION_A      = 3.2         # enforce if numeric; keep NaN (NMR/cryo-EM)
MAX_CLASHSCORE        = 15.0        # enforce if numeric; column may not exist

# If you have many SDFs per system and want to relax strictness:
STRICT_SINGLE_SDF     = True        # if False, pick the largest .sdf by file size

# -------------------------------------------------
# Paths (point to your local PLINDER copy)
# -------------------------------------------------
PL_ROOT   = Path("plinder/v2")
IDX_DIR   = PL_ROOT / "index"
SPLIT_PARQ= PL_ROOT / "splits" / "split.parquet"
SYS_DIR   = PL_ROOT / "systems"     # contains either zip chunks or unpacked system dirs

# ----------------------------
# Helpers
# ----------------------------
def safe_float(x):
    try:
        return float(x)
    except Exception:
        return None

def quality_score(row, clash_col):
    """
    Higher is better.
    Combines low resolution, low clashscore, high interactions.
    Missing values are downweighted but not fatal.
    """
    res   = safe_float(row.get("entry_resolution"))
    clash = safe_float(row.get(clash_col)) if clash_col else None
    ints  = safe_float(row.get("system_proper_num_interactions"))

    # Normalize to 0..1 with gentle caps
    res_s   = None if res   is None else max(0.0, min(1.0, 1.0 - (res - 1.0) / 3.0))  # 1Å->1.0, 4Å->0
    clash_s = None if clash is None else max(0.0, min(1.0, 1.0 - (clash / 50.0)))     # 0->1.0, 50->0
    ints_s  = None if ints  is None else max(0.0, min(1.0, (ints / 20.0)))            # >=20 -> 1.0

    parts = []
    if res_s   is not None: parts.append(0.35 * res_s)
    if clash_s is not None: parts.append(0.25 * clash_s)
    if ints_s  is not None: parts.append(0.40 * ints_s)
    return sum(parts) if parts else 0.0

def ensure_split_column(df):
    """Normalize any split-like column to 'split', or create default 'all'."""
    candidates = ["split", "split_name", "split_set", "dataset_split"]
    found = None
    for c in candidates:
        if c in df.columns:
            found = c
            break
    if found is None:
        df["split"] = "all"
    elif found != "split":
        df = df.rename(columns={found: "split"})
    return df

def index_system_dirs(systems_root: Path) -> dict:
    """Map system_id -> system_dir by scanning systems/*/<system_id>."""
    sys_map = {}
    if not systems_root.exists():
        return sys_map
    for chunk in systems_root.iterdir():
        if not chunk.is_dir():
            continue
        # Some releases flatten systems directly under systems_root
        lig_dir = chunk / "ligand_files"
        if lig_dir.exists():
            sys_map[chunk.name] = chunk
            continue
        for sdir in chunk.iterdir():
            if sdir.is_dir():
                sys_map[sdir.name] = sdir
    return sys_map

def ensure_system_dirs_unzipped(systems_root: Path) -> dict:
    """
    Make sure system directories exist.
    If only zip chunks are present, extract them in-place and re-index.
    """
    zips = sorted(p for p in systems_root.glob("*.zip") if p.is_file())
    to_extract: list[tuple[Path, set[str]]] = []
    for zpath in zips:
        try:
            with ZipFile(zpath) as zf:
                members = {Path(name).parts[0] for name in zf.namelist() if name and not name.endswith("/")}
        except Exception as exc:
            raise SystemExit(f"Failed to inspect {zpath.name}: {exc}") from exc
        if not members:
            continue
        if all((systems_root / member).exists() for member in members):
            continue
        to_extract.append((zpath, members))

    if to_extract:
        print(f"[info] extracting {len(to_extract)} zip chunks into {systems_root} ...")
        for zpath, _ in tqdm(to_extract, desc="Extracting system zips", unit="zip"):
            try:
                with ZipFile(zpath) as zf:
                    zf.extractall(systems_root)
            except Exception as exc:
                raise SystemExit(f"Failed to extract {zpath.name}: {exc}") from exc

    return index_system_dirs(systems_root)

def copy_pattern(system_row: dict, out_root: Path, sdir: Path) -> tuple[bool, str]:
    """
    Copy receptor.pdb -> pocket.pdb (or receptor.cif -> pocket.cif)
    and the single ligand SDF -> ligand.sdf.
    Returns (ok, reason).
    """
    sid = system_row["system_id"]

    receptor_pdb = sdir / "receptor.pdb"
    receptor_cif = sdir / "receptor.cif"
    lig_dir      = sdir / "ligand_files"

    if not lig_dir.is_dir():
        return False, "ligand_dir_missing"

    sdf_files = list(lig_dir.glob("*.sdf"))
    if STRICT_SINGLE_SDF:
        if len(sdf_files) != 1:
            return False, f"sdf_count_{len(sdf_files)}"
        sdf_pick = sdf_files[0]
    else:
        if not sdf_files:
            return False, "sdf_count_0"
        # pick largest file
        sdf_pick = max(sdf_files, key=lambda p: p.stat().st_size)

    if not receptor_pdb.exists() and not receptor_cif.exists():
        return False, "no_receptor_pdb_or_cif"

    patt_dir = out_root / sid
    patt_dir.mkdir(parents=True, exist_ok=True)

    if receptor_pdb.exists():
        shutil.copy2(receptor_pdb, patt_dir / "pocket.pdb")
    else:
        shutil.copy2(receptor_cif, patt_dir / "pocket.cif")

    shutil.copy2(sdf_pick, patt_dir / "ligand.sdf")

    # Minimal metadata
    meta = {
        "system_id": sid,
        "entry_pdb_id": system_row.get("entry_pdb_id"),
        "entry_resolution": system_row.get("entry_resolution"),
        "entry_validation_clashscore": system_row.get("entry_validation_clashscore"),
        "system_proper_num_interactions": system_row.get("system_proper_num_interactions"),
        "uniqueness": system_row.get("uniqueness"),
        "split": system_row.get("split", "all"),
    }
    with open(patt_dir / "meta.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)

    return True, "ok"

# ----------------------------
# Load metadata (LOCAL parquet)
# ----------------------------
print(f"Using PLINDER at: {PL_ROOT}")
if not PL_ROOT.exists():
    raise SystemExit("PLINDER directory not found.")

if not SPLIT_PARQ.exists():
    raise SystemExit(f"Missing splits file: {SPLIT_PARQ}")

# Choose local annotation_table*.parquet
cand = sorted(IDX_DIR.glob("annotation_table*.parquet"))
if not cand:
    raise SystemExit(f"No annotation_table*.parquet found in {IDX_DIR}")
INDEX_PARQ = cand[0]
print(f"[info] using index parquet: {INDEX_PARQ.name}")

# Load only columns that exist
if pq is not None:
    have_idx = set(pq.ParquetFile(INDEX_PARQ).schema.names)
else:
    have_idx = set(pd.read_parquet(INDEX_PARQ).columns)
need_idx = {"system_id", "entry_pdb_id", "entry_resolution", "entry_validation_clashscore", "system_num_ligand_chains"}
use_cols = [c for c in need_idx if c in have_idx]
idx_df = pd.read_parquet(INDEX_PARQ, columns=use_cols).dropna(subset=["system_id","entry_pdb_id"])

split_df = pd.read_parquet(SPLIT_PARQ)
# some releases keep interactions & uniqueness here
if "system_proper_num_interactions" not in split_df.columns:
    split_df["system_proper_num_interactions"] = pd.NA

# Merge & normalize split column
df = pd.merge(idx_df, split_df, on="system_id", how="inner")
df = df.dropna(subset=["system_id"]).copy()
df = ensure_split_column(df)

print("after merge:", len(df))

# ----------------------------
# Filter & QC
# ----------------------------
if REQUIRE_SINGLE_LIGAND:
    if "system_num_ligand_chains" in df.columns:
        df = df[df["system_num_ligand_chains"] == 1]
    elif "system_proper_num_ligand_chains" in df.columns:
        df = df[df["system_proper_num_ligand_chains"] == 1]
print("after single-ligand:", len(df))

if MIN_INTERACTIONS is not None and "system_proper_num_interactions" in df.columns:
    df = df[df["system_proper_num_interactions"].fillna(0) >= MIN_INTERACTIONS]
print("after min interactions:", len(df))

if MAX_RESOLUTION_A is not None and "entry_resolution" in df.columns:
    df = df[df["entry_resolution"].apply(safe_float).isna() | (df["entry_resolution"] <= MAX_RESOLUTION_A)]

clash_col = "entry_validation_clashscore" if "entry_validation_clashscore" in df.columns else None
if MAX_CLASHSCORE is not None and clash_col:
    df = df[df[clash_col].apply(safe_float).isna() | (df[clash_col] <= MAX_CLASHSCORE)]
print("after QC gates:", len(df))

# ----------------------------
# Deduplicate
# ----------------------------
if DEDUP_STRATEGY == "nonredundant":
    nonred_path = IDX_DIR / "annotation_table_nonredundant.parquet"
    if not nonred_path.exists():
        raise SystemExit(f"Requested nonredundant strategy, but file not found: {nonred_path}")
    nonred = pd.read_parquet(nonred_path, columns=["system_id"])
    df = df.merge(nonred, on="system_id", how="inner")
else:
    if "uniqueness" in df.columns:
        df["__q"] = df.apply(lambda r: quality_score(r, clash_col), axis=1)
        df = df.sort_values("__q", ascending=False)
        df = df.drop_duplicates(subset=["uniqueness"], keep="first").drop(columns="__q")
    else:
        print("[warn] 'uniqueness' not in splits — skipping dedup by uniqueness.")
print("after dedup:", len(df))

# Optional: cap number of patterns for a quick dry-run
if isinstance(MAX_PATTERNS, int) and MAX_PATTERNS > 0:
    df = df.head(MAX_PATTERNS)
    print(f"[info] capped to first {MAX_PATTERNS} rows for dry-run.")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ----------------------------
# Build system_id -> path map
# ----------------------------
sys_map = ensure_system_dirs_unzipped(SYS_DIR)
print(f"[info] system dirs indexed: {len(sys_map)}")

if not sys_map:
    zips = list(SYS_DIR.glob("*.zip"))
    if zips:
        raise SystemExit("No system directories found even after extracting zip chunks. Check archive layout.")
    raise SystemExit("No system directories found under systems/. Check your data layout.")

# ----------------------------
# Materialize patterns
# ----------------------------
kept = skipped = 0
reasons = Counter()

rows = df.to_dict(orient="records")
for row in tqdm(rows, desc="Creating patterns"):
    sid = row["system_id"]
    sdir = sys_map.get(sid)
    if sdir is None:
        skipped += 1
        reasons["system_dir_not_found"] += 1
        continue
    try:
        ok, why = copy_pattern(row, OUTPUT_DIR, sdir)
    except Exception:
        ok, why = False, "exception"
    if ok:
        kept += 1
    else:
        skipped += 1
        reasons[why] += 1

# Save inventory (only existing cols)
inv_cols_wanted = [
    "system_id", "entry_pdb_id", "split", "uniqueness",
    "entry_resolution", "entry_validation_clashscore",
    "system_proper_num_interactions",
]
inv_cols = [c for c in inv_cols_wanted if c in df.columns]
(df[inv_cols]
   .drop_duplicates("system_id", keep="first")
   .to_csv(OUTPUT_DIR / "inventory.csv", index=False))

print(f"\nDone. Patterns created: {kept}, skipped: {skipped}")
print("Skip reasons:", dict(reasons))
print(f"Output directory: {OUTPUT_DIR.resolve()}")
