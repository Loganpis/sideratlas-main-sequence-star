import json,sys
from pathlib import Path
import numpy as np
from scipy.ndimage import label,gaussian_filter
from scipy.fft import fft2,ifft2
from scipy.optimize import curve_fit
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'src'))
from continuum import inputs
from solar_surface import YROW

def field(v,s,z0):
 d=np.load(ROOT/f'data/raw/{v}-{s}.npz');k=np.argmin(abs(d['z_Mm']-z0));a=d['values'][k];return (10**a if v.startswith('lg') else a),float(d['z_Mm'][k])
T,ne,nh,z,dx,t,teff=inputs(385)
i=np.argmin(abs(z));temp=T[i];uz,zh=field('uz',385,.061);bz,_=field('bz',385,.061);rho,_=field('lgr',385,.061)
d=np.load(ROOT/'data/continuum/385-17-2.npz');rgb=d['rgb'][-1,0];y=rgb@YROW;contrast=float(y.std()/y.mean())
# Granule segmentation is diagnostic, not a generator. Connected bright
# basins are biased by threshold and fragmentation; report its definition.
ys=gaussian_filter(y,1.0);mask=ys>np.median(ys);lab,n=label(mask);sizes=np.bincount(lab.ravel())[1:];sizes=sizes[sizes>=9]
diam=np.sqrt(4*sizes/np.pi)*dx*2/1e3
spectrum=abs(fft2(y/y.mean()-1))**2;freq=np.fft.fftfreq(y.shape[0],d=dx*2/1e6);kk=np.hypot(freq[:,None],freq[None,:]);bins=np.linspace(.1,3,60);which=np.digitize(kk.ravel(),bins);shell=np.array([spectrum.ravel()[which==k].sum() for k in range(1,len(bins))]);freqpeak=.5*(bins[:-1]+bins[1:])[np.argmax(shell)]
# Correlate at the actual archived mean optical-depth-unity reference height.
u=uz[::2,::2];b=abs(bz[::2,::2])*1e4
corr=lambda a,b:float(np.corrcoef(a.ravel(),b.ravel())[0,1])
result=dict(source='Carlsson et al. 2016 Bifrost en024048_hion',snapshot=385,source_teff_K=teff,solar_g_m_s2=274.,grid_horizontal_km=dx/1000,grid_photosphere_vertical_km=float(np.median(np.diff(z)))/1000,selected_height_Mm=zh,velocity_Bunit='archived SI: m/s and tesla',vertical_speed_rms_km_s=float(np.std(uz)/1000),vertical_speed_percentiles_km_s=(np.percentile(uz,[1,50,99])/1000).tolist(),mean_unsigned_vertical_B_G=float(np.mean(abs(bz))*1e4),vertical_B_abs_99_G=float(np.percentile(abs(bz),99)*1e4),continuum_upflow_correlation=corr(y,u),continuum_contrast=contrast,granule_threshold_diameter_median_km=float(np.median(diam)),granule_threshold_diameter_10_90_km=np.percentile(diam,[10,90]).tolist(),intensity_power_shell_peak_scale_Mm=float(1/freqpeak),segmentation_definition='connected components above median of sigma=1-pixel smoothed native continuum map; >=9 samples; diagnostic only',no_sunspot_model=True,energy_validation_scope='Published relaxed radiation-MHD model near solar Teff; no independent full MHD resolution or energy-closure test from these cropped snapshots')
paths=sorted((ROOT/'data/continuum').glob('*-17-2-solar.npz'));cor=[];time=[]
for p in paths:
 a=np.load(p);yy=a['rgb'][-1,0]@YROW;time.append(float(a['time_s'])-t);cor.append(corr(y,yy))
result['same_position_temporal_correlation']={'lag_s':time,'correlation':cor}
(ROOT/'validation/source-statistics.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
