import os.path
from pathlib import Path
import urllib.request
import os

input_dir = Path("example")
output_dir = Path("example")

target = "example (3rfm)" #@param ["example (3rfm)", "upload structure"]

if target == "example (3rfm)":
  pdbfile = Path(input_dir, '3rfm.pdb')
  urllib.request.urlretrieve('http://files.rcsb.org/download/3rfm.pdb', pdbfile)

elif target == "upload structure":
  uploaded = files.upload()
  fn = list(uploaded.keys())[0]
  pdbfile = Path(input_dir, fn)
  Path(fn).rename(pdbfile)

import sys
sys.path.append("/usr/local/lib/python3.9/site-packages")
import py3Dmol

view = py3Dmol.view(js='https://3dmol.org/build/3Dmol.js',)
view.addModel(open(Path("example", "3rfm.pdb"), 'r').read(), 'pdb')
view.setStyle({'model': -1}, {'cartoon': {'color': 'lime'}})
# view.addSurface(py3Dmol.VDW, {'opacity': 0.4, 'color': 'lime'})
view.addModelsAsFrames(open(Path("example", "3rfm_mol.sdf"), 'r').read())
view.setStyle({'model': -1}, {'stick': {}})
view.zoomTo({'model': -1})
view.zoom(0.5)
if target == "example (3rfm)":
  view.rotate(90, 'y')
view.animate({'loop': "forward", 'interval': 1000})
out = Path("structure_view.html")
out.write_text(view._make_html())   # tạo file HTML tự chứa
print(f"Saved viewer to: {out.resolve()}\nOpen it in a browser.")  
