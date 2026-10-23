"""Retrieve and decode the CC BY 3.0 Stanford MURaM intensity sequence."""
from pathlib import Path
import requests,tarfile,numpy as np,hashlib,json
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'data/active-region';OUT.mkdir(parents=True,exist_ok=True)
URL='https://stacks.stanford.edu/file/dv883vb9686/I_out_every200steps.tar'
SHA='90c7afc88f93d4092f63152c5aaa3019e51e9b275db14edc2152e5f5bcb0482e'
p=OUT/'I_out_every200steps.tar'
if not p.exists() or hashlib.sha256(p.read_bytes()).hexdigest()!=SHA:
 with requests.get(URL,stream=True,timeout=(20,150)) as r:
  r.raise_for_status()
  with p.open('wb') as f:
   for c in r.iter_content(1024*1024):f.write(c)
if hashlib.sha256(p.read_bytes()).hexdigest()!=SHA:raise ValueError('Archive changed; review provenance')
a=[];times=[];iters=[]
with tarfile.open(p,'r:') as tar:
 for m in tar:
  if not m.isfile():continue
  data=np.frombuffer(tar.extractfile(m).read(),dtype='<f4');nx,ny=map(int,data[1:3]);assert len(data)==4+nx*ny
  a.append(data[4:].reshape(ny,nx));times.append(float(data[3]));iters.append(int(m.name.split('.')[-1]))
order=np.argsort(times)
np.savez_compressed(OUT/'intensity-sequence.npz',intensity_cgs=np.stack(a)[order],time_s=np.array(times)[order],iteration=np.array(iters)[order],dx_Mm=.192,dy_Mm=.192)
(OUT/'intensity-provenance.json').write_text(json.dumps(dict(url=URL,sha256=SHA,license='CC BY 3.0',citation='Cheung et al., Nature Astronomy 3,160–166 (2019); online 2018',doi='10.1038/s41550-018-0629-3',archive='https://purl.stanford.edu/dv883vb9686',native_grid_km=192,dimensions_Mm=[98.304,49.152],modifications='Float32 arrays decoded in documented native layout; only sorted by actual stored timestamps; no generated morphology or dynamics'),indent=2)+'\n')
print('Verified and decoded',len(a),'native MURaM intensity maps')
