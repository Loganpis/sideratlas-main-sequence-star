import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
import numpy as np
from continuum import formal,ROOT
from geometry import linear_to_srgb
from PIL import Image
mu=np.array([.1,.3,1.]);az=np.array([0.,1.]);z=np.linspace(0,1,75)
a=np.full((75,8,8),2.,np.float32);s=np.full_like(a,3.)
x=formal(a,s,z,1.,mu,az,1)
assert np.max(abs(x-3))<1e-5
s=1+np.broadcast_to(z[:,None,None],a.shape).astype('f4')
x=formal(a,s,z,1.,mu,az,1)
expected=2-(1-np.exp(-2/mu))*mu/2
err=float(np.max(abs(x-expected[:,None,None,None])))
assert err<2e-6
p=ROOT/'data/continuum/385-17-2.npz';d=np.load(p);rgb=d['rgb'];assert np.isfinite(rgb).all() and (rgb>0).all()
yrow=np.linalg.inv(np.array([[3.2404542,-1.5371385,-.4985314],[-.969266,1.8760108,.041556],[.0556434,-.2040259,1.0572252]]))[1]
y=rgb@yrow
center=rgb[-1,0];norm=float((center@yrow).mean());im=1-np.exp(-center/norm)
Image.fromarray(np.uint8(np.clip(linear_to_srgb(im),0,1)*255)).resize((756,756)).save(ROOT/'validation/first-continuum.png')
result=dict(formal_solver_isothermal_max_error=float(np.max(abs(formal(a,np.full_like(a,3.),z,1.,mu,az,1)-3))),linear_source_max_error=err,limb_darkening=[float(v) for v in y.mean(axis=(1,2,3))/y[-1].mean()],finite_positive_RGB=True,disk_center_broadband_rms_contrast=float(y[-1,0].std()/y[-1,0].mean()),checks='Analytic constant and linearly varying source functions; no MHD convergence claim')
(ROOT/'validation/transfer.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
