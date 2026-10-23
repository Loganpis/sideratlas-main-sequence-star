"""Validate the actual MURaM radiation/field data and its rendered sampling."""
from pathlib import Path
import sys,json
import numpy as np
from scipy.ndimage import gaussian_filter,label,center_of_mass
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'src'))
from active_region import ActiveRegion,active_sample,pixel
from solar_surface import Surface,YROW
cfg=json.loads((ROOT/'config-active.json').read_text());q=Surface(cfg);ar=ActiveRegion(cfg,q)
a=ar.intensity;times=ar.times;b=np.load(ROOT/'data/active-region/magnetic/270000.npz');bz=abs(b['bz_G']);bt=np.sqrt(b['bx_G']**2+b['by_G']**2+b['bz_G']**2);r=a[0]/ar.quiet
ny,nx=r.shape;yy,xx=np.mgrid[:ny,:nx];x=xx*.192;y=yy*.192;window=(x-73.75)**2+(y-27.65)**2<16**2
umbra=(r<.5)&window;pen=(r>=.5)&(r<.9)&(bt>500)&window;fac=(r>1.1)&(bz>500)&window
stats={}
for name,m in [('dark_cores',umbra),('magnetic_intermediate_intensity',pen),('magnetic_bright_points',fac)]:
 stats[name]=dict(native_pixel_count=int(m.sum()),median_vertical_field_G=float(np.median(bz[m])),median_total_field_G=float(np.median(bt[m])),mean_physical_intensity_over_quiet=float(np.mean(r[m])))
mid=[]
for j in range(1,len(a)-1,2):
 f=(times[j]-times[j-1])/(times[j+1]-times[j-1]);estimate=(1-f)*a[j-1]+f*a[j+1];mid.append(float(np.sqrt(np.mean((estimate-a[j])**2))/ar.quiet))
# Exact unsmoothed sampling when display contrast exponent is one.
y0,x0=np.array([27.65,73.75])/.192
expected=float(pixel(a,0,x0,y0)/ar.quiet)
col,w=active_sample(ar.center[None,:],np.zeros(1),float(times[0]),ar.center,ar.east,ar.north,a,ar.large_scale,times,np.ones(1),q.mu,q.mean_mu_y,ar.colors,ar.tempgrid,15.,np.array([73.75,27.65]),ar.quiet,1.)
tb0=5772*expected**.25;correction=np.array([np.interp(tb0,ar.tempgrid,ar.colors[:,c]) for c in range(3)])
error=abs(float(col[0]@YROW)/(expected*(correction@YROW))-1);assert error<1e-6 and w[0]==1
# Native source temporal evolution, no camera movement or imposed animation.
k=int(np.argmin(abs(times-times[0]-320)))
corr=lambda m:float(np.corrcoef(a[0][m],a[k][m])[0,1])
result=dict(dataset='https://purl.stanford.edu/dv883vb9686',license='CC BY 3.0',native_map_shape=list(a.shape[1:]),horizontal_spacing_km=192,source_domain_Mm=[98.304,49.152],timestamps_s=[float(times[0]),float(times[-1])],number_source_frames=len(times),source_cadence_s=np.percentile(np.diff(times),[0,50,100]).tolist(),field_plane='tau=0.1; native reader planes 5,6,7 times sqrt(4 pi) -> bx,bz,by in gauss',fixed_quiet_intensity_cgs=ar.quiet,quiet_reference='mean first-frame outgoing intensity where total field <100 G; no per-frame normalization',physical_intensity_classes=stats,temporal_holdout_rms_relative_to_quiet=dict(mean=float(np.mean(mid)),maximum=float(max(mid))),radiance_sampling_relative_error=error,physical_source_correlation_after_320s=dict(all_patch=corr(window),weak_field_patch=corr(window&(bt<100))),new_MHD_solver_run=False,limitations=['Gray disk-center radiation with a thermal optical-band correction, not local line-by-line spectral synthesis','Quiet-star limb law is used for active region; Wilson depression and directional facular hot-wall contrast are not newly solved','192 km grid cannot resolve the finest penumbral filaments','Patch is blended into independently simulated Bifrost background; not a globally coupled MHD solution','One 40:1 clock, separate source epochs, fixed 15x magnification; no fabricated plasma advection'])
(ROOT/'validation/active-region.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
