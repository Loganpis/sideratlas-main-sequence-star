import sys,json,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'src'))
from continuum import *

def relative(a,b):
 return dict(mean_relative_difference=float(abs(a.mean()/b.mean()-1)),rms_relative_to_mean=float(np.sqrt(np.mean((a-b)**2))/b.mean()),max_relative_to_mean=float(np.max(abs(a-b))/b.mean()))
result={}
# Wavelength convergence in the spectral hybrid; same input atmosphere/rays.
synthesize(385,33,2,True)
a=np.load(ROOT/'data/continuum/385-17-2-solar.npz')['rgb'];b=np.load(ROOT/'data/continuum/385-33-2-solar.npz')['rgb'];result['17_vs_33_wavelengths']=relative(a,b)
# Integrating the same piecewise-linear atmosphere with half layer spacing.
T,ne,nh,z,dx,t,teff=inputs(385);w=550.
bf,ff=coefficients(w,T);chi,S=radiation_arrays(T,ne,nh,ff.astype('f4'),bf,w)
mu=np.array([.1,.3,.55,.78,1.]);az=np.arange(4)*np.pi/2
coarse=formal(chi,S,z,dx,mu,az,4)
z2=np.empty(2*len(z)-1);z2[::2]=z;z2[1::2]=.5*(z[:-1]+z[1:])
ch2=np.empty((len(z2),)+chi.shape[1:],'f4');s2=np.empty_like(ch2);ch2[::2]=chi;s2[::2]=S;ch2[1::2]=.5*(chi[:-1]+chi[1:]);s2[1::2]=.5*(S[:-1]+S[1:])
fine=formal(ch2,s2,z2,dx,mu,az,4);result['vertical_formal_integration_refinement']=relative(coarse,fine)
# Direct intermediate oblique ray vs the preview's angular interpolation.
true=formal(chi,S,z,dx,np.array([.665]),np.array([np.pi/4]),4)[0,0]
interp=.25*(coarse[2,0]+coarse[2,1]+coarse[3,0]+coarse[3,1]);result['oblique_angular_interpolation']=relative(interp,true)
# Continuum output stride does not coarsen the source volume or transfer ray.
full=formal(chi,S,z,dx,np.array([1.]),np.array([0.]),1)[0,0]
result['output_grid_point_consistency']=relative(coarse[-1,0],full[::4,::4])
while not all((ROOT/f'data/raw/{v}-387.json').exists() for v in ['lgtg','lgne','lgn1']):time.sleep(5)
synthesize(387,17,2,True)
while not (ROOT/'data/continuum/389-17-2-solar.npz').exists():time.sleep(5)
a=np.load(ROOT/'data/continuum/385-17-2-solar.npz')['rgb'];b=np.load(ROOT/'data/continuum/389-17-2-solar.npz')['rgb'];m=np.load(ROOT/'data/continuum/387-17-2-solar.npz')['rgb'];result['40_second_temporal_midpoint']=relative(.5*(a+b),m)
result['scope']='Post-processing convergence only, not independent radiation-MHD grid convergence; maximum point errors may be larger in magnetic lanes'
(ROOT/'validation/convergence.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
