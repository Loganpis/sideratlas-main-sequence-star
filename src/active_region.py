"""Time-dependent MURaM photospheric active-region patch (CC BY 3.0).
Geometry, dark/bright fine structure and their evolution come from computed
outgoing radiation; no spots, filaments or bright-point masks are painted.
"""
from pathlib import Path
import numpy as np
from numba import njit,prange
from continuum import ROOT,H,C,K,XYZ_RGB
from solar_surface import YROW
from geometry import orientation

@njit(cache=True)
def pixel(a,it,x,y):
    ix=int(np.floor(x));iy=int(np.floor(y));fx=x-ix;fy=y-iy
    return (1-fy)*((1-fx)*a[it,iy,ix]+fx*a[it,iy,ix+1])+fy*((1-fx)*a[it,iy+1,ix]+fx*a[it,iy+1,ix+1])

@njit(parallel=True,cache=True)
def active_sample(n,delay,time,center,east,north,fields,large_scale,times,mu,ld_mu,ld_y,colors,tempgrid,mag,sourcexy,quiet,contrast):
    out=np.zeros((len(n),3),np.float32);weight=np.zeros(len(n),np.float32)
    for j in prange(len(n)):
        front=np.dot(n[j],center)
        if front<=0:continue
        u=np.arctan2(np.dot(n[j],east),front)*695.7/mag
        v=np.arctan2(np.dot(n[j],north),front)*695.7/mag
        r=np.hypot(u,v)
        if r>=20:continue
        w=1. if r<=16 else .5*(1+np.cos(np.pi*(r-16)/4))
        x=(u+sourcexy[0])/.192;y=(v+sourcexy[1])/.192
        if x<0 or y<0 or x>=fields.shape[2]-1 or y>=fields.shape[1]-1:continue
        t=time-delay[j]
        it=max(0,min(len(times)-2,np.searchsorted(times,t)-1));ft=(t-times[it])/(times[it+1]-times[it])
        ratio=((1-ft)*pixel(fields,it,x,y)+ft*pixel(fields,it+1,x,y))/quiet
        baseline=((1-ft)*pixel(large_scale,it,x,y)+ft*pixel(large_scale,it+1,x,y))/quiet
        # Gray outgoing intensity is not a measured Teff map. The fourth-root
        # brightness-temperature proxy supplies a documented thermal color
        # approximation; actual MURaM outgoing intensity supplies luminance.
        tb=5772*ratio**.25;k=max(0,min(len(tempgrid)-2,np.searchsorted(tempgrid,tb)-1));f=min(1.,max(0.,(tb-tempgrid[k])/(tempgrid[k+1]-tempgrid[k])))
        chrom=(1-f)*colors[k]+f*colors[k+1]
        ld=np.interp(mu[j],ld_mu,ld_y)
        # Same fixed contrast exponent as the quiet photosphere. This acts
        # on already-calculated radiation, never on a synthetic spot shape.
        out[j]=ld*baseline*max(0.,ratio/baseline)**contrast*chrom
        weight[j]=w
    return out,weight

class ActiveRegion:
    def __init__(self,cfg,quiet_surface):
        self.cfg=cfg['active_region'];self.clock=cfg['physical_seconds_per_display_second'];d=np.load(ROOT/'data/active-region/intensity-sequence.npz')
        take=d['time_s']<=d['time_s'][0]+1320
        self.intensity=d['intensity_cgs'][take]
        from scipy.ndimage import gaussian_filter
        self.large_scale=gaussian_filter(self.intensity,sigma=(0,self.cfg.get('baseline_sigma_native_pixels',2),self.cfg.get('baseline_sigma_native_pixels',2)),mode='wrap')
        self.times=d['time_s'][take];self.start=float(self.times[0])+2.4
        b=np.load(ROOT/'data/active-region/magnetic/270000.npz');btotal=np.sqrt(b['bx_G']**2+b['by_G']**2+b['bz_G']**2)
        quiet=btotal<100
        self.quiet=float(np.mean(self.intensity[0][quiet]))
        xy=np.array(self.cfg['center_camera_at_start']);n=np.r_[xy,np.sqrt(1-np.dot(xy,xy))]
        from geometry import camera_direction
        inc,az=camera_direction(0,cfg)
        self.center=orientation(n,inc,az);self.east=np.cross(np.array([0.,0.,1.]),self.center);self.east/=np.linalg.norm(self.east);self.north=np.cross(self.center,self.east)
        self.surface=quiet_surface
        self.tempgrid=np.linspace(2500,9000,651);wave=np.linspace(380,800,421);cmf=np.loadtxt(ROOT/'data/CIE_xyz_1931_2deg.csv',delimiter=',');xyz=np.stack([np.interp(wave,cmf[:,0],cmf[:,j]) for j in [1,2,3]],-1)
        import h5py
        with h5py.File(ROOT/'data/tsis1-hsrs-v2-1nm.nc') as f:ref=np.interp(wave,np.array(f['Vacuum Wavelength']),np.array(f['SSI']))
        lam=wave*1e-9;refplanck=1/np.expm1(H*C/(lam*K*5772))
        flux=ref[None,:]/np.expm1(H*C/(lam[None,:]*K*self.tempgrid[:,None]))/refplanck[None,:]
        colors=np.trapz(flux[:,:,None]*xyz[None,:,:],wave,axis=1)@XYZ_RGB.T
        ref_y=float(np.interp(5772,self.tempgrid,colors@YROW))
        # Divide by the gray bolometric ratio because active_sample already
        # multiplies by that radiation. This retains the optical fraction
        # of the thermal spectrum instead of treating gray and optical
        # intensity as interchangeable.
        self.colors=colors/(ref_y*(self.tempgrid/5772)**4)[:,None]
    def sample(self,samples,t,inc,az):
        if not len(samples.mu):return np.zeros((0,3),np.float32),np.zeros(0,np.float32)
        physical=t*self.clock;rotation=physical/(25.38*86400)*360
        n=orientation(samples.camera,inc,az-rotation)
        if self.start+physical>self.times[-1]:raise ValueError('MURaM sequence too short')
        return active_sample(n,samples.delay,self.start+physical,self.center,self.east,self.north,self.intensity,self.large_scale,self.times,samples.mu,self.surface.mu,self.surface.mean_mu_y,self.colors,self.tempgrid,self.cfg['spatial_magnification'],np.array(self.cfg['source_center_Mm']),self.quiet,self.cfg['contrast_exponent'])
