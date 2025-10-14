import os
import pickle
import shutil
import sys

if len(sys.argv) < 3:
    print("Usage: python move_patterns_from_index.py <index.pkl> <destination_root>")
    sys.exit(1)

index_path = sys.argv[1]
dest_root = sys.argv[2]

# Root dataset folder (same directory as index.pkl)
base_dir = os.path.dirname(index_path)
dataset_root = os.path.join(base_dir, "crossdocked_pocket10")

os.makedirs(dest_root, exist_ok=True)

print(f"[INFO] Loading index file: {index_path}")
with open(index_path, "rb") as f:
    data = pickle.load(f)

print(f"[INFO] Loaded {len(data)} entries.")
print(f"[INFO] Dataset root assumed as: {dataset_root}")

def safe_copy(src, dst_folder):
    """Copy a file if it exists."""
    if not src or not isinstance(src, str):
        print(f"[WARN] Invalid path (None or not str): {src}")
        return False
    if not os.path.exists(src):
        print(f"[WARN] Missing file: {src}")
        return False
    shutil.copy2(src, os.path.join(dst_folder, os.path.basename(src)))
    return True

skipped = 0
copied = 0

for i, entry in enumerate(data):
    if not isinstance(entry, (list, tuple)) or len(entry) < 2:
        print(f"[WARN] Skipping invalid entry #{i}: {entry}")
        skipped += 1
        continue

    rel_pocket_pdb, rel_ligand_sdf = entry[:2]

    # skip if any is None or empty
    if not rel_pocket_pdb or not rel_ligand_sdf:
        print(f"[WARN] Entry #{i} missing pocket or ligand — skipped.")
        skipped += 1
        continue

    # prepend dataset root to make full paths
    pocket_pdb = os.path.join(dataset_root, rel_pocket_pdb)
    ligand_sdf = os.path.join(dataset_root, rel_ligand_sdf)

    # folder name = ligand file name (without .sdf)
    ligand_name = os.path.splitext(os.path.basename(ligand_sdf))[0]
    target_dir = os.path.join(dest_root, ligand_name)
    os.makedirs(target_dir, exist_ok=True)

    copied_flag = False
    if safe_copy(pocket_pdb, target_dir):
        copied_flag = True
    if safe_copy(ligand_sdf, target_dir):
        copied_flag = True

    if copied_flag:
        copied += 1
        print(f"[{i+1}/{len(data)}] Copied to {target_dir}")
    else:
        skipped += 1

print(f"\n✅ Done! {copied} patterns copied. {skipped} skipped due to missing info.")
