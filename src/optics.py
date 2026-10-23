"""Reused normalized Airy/jitter/Moffat imaging from SiderAtlas No.003.
Solar spectrum replaces DA atmosphere. Halo gain is a documented display aid.
"""
import numpy as np,h5py
from scipy.fft import rfft2,irfft2,rfftfreq,fftfreq
from continuum import ROOT

def airy_otf(frequency,first_zero):
    rho=np.minimum(1,frequency*(np.pi/3.8317059702075125)*first_zero)
    return 2/np.pi*(np.arccos(rho)-rho*np.sqrt(np.maximum(0,1-rho*rho)))

class Optics:
    def __init__(self,shape,cfg):
        self.h,self.w=shape;self.cfg=cfg;scale=self.w/1080;self.pad=round(300*scale)
        self.shape=(self.h+2*self.pad,self.w+2*self.pad)
        freq=np.hypot(fftfreq(self.shape[0])[:,None],rfftfreq(self.shape[1])[None,:])
        if 'monochromatic_vacuum_nm' in cfg:
            wave=np.array([cfg['monochromatic_vacuum_nm']]);weights=np.ones(1)
        else:
            with h5py.File(ROOT/'data/tsis1-hsrs-v2-1nm.nc') as f:
                w=np.array(f['Vacuum Wavelength']);s=np.array(f['SSI'])
            wave=np.linspace(380,800,65);cmf=np.loadtxt(ROOT/'data/CIE_xyz_1931_2deg.csv',delimiter=',')
            weights=np.interp(wave,w,s)*np.interp(wave,cmf[:,0],cmf[:,2]);weights[[0,-1]]*=.5;weights/=weights.sum()
        core=np.zeros(freq.shape)
        for lam,wt in zip(wave,weights):core+=wt*airy_otf(freq,cfg['airy_first_zero_550nm_at1080']*scale*lam/550)
        self.core=core*np.exp(-2*np.pi**2*(cfg['jitter_sigma_at1080']*scale)**2*freq**2)
        q=2*np.pi*cfg['moffat_alpha_at1080']*scale*freq;self.wing=(1+q)*np.exp(-q)
        self.flux_errors=[];self.crop_fractions=[]
        self.display_halo_mask=None
        if cfg.get('extra_halo_outside_disk_only',False):
            yy,xx=np.mgrid[:self.shape[0],:self.shape[1]]
            cx,cy=np.array(cfg['halo_mask_center_fraction'])*[self.w,self.h]+self.pad
            r=cfg['halo_mask_radius_fraction_of_width']*self.w
            f=np.clip((np.hypot(xx-cx,yy-cy)/r-.85)/.25,0,1)
            self.display_halo_mask=f*f*(3-2*f)

    def convolve(self,plane,otf):
        a=np.pad(plane,self.pad);return irfft2(rfft2(a,workers=2)*otf,s=self.shape,workers=2)
    def process(self,rgb):
        eta=self.cfg['scatter_fraction'];core=np.stack([self.convolve(rgb[...,c],self.core) for c in range(3)],-1)
        halo=np.stack([self.convolve(rgb[...,c],self.wing) for c in range(3)],-1)
        physical=(1-eta)*core+eta*halo
        self.flux_errors.append(float(abs(physical.sum()/rgb.sum()-1)))
        p=self.pad;self.crop_fractions.append(float(physical[p:p+self.h,p:p+self.w].sum()/physical.sum()))
        if self.cfg.get('extra_halo_outside_disk_only',False):
            # Explicit display-only glare suppression: retain the normalized
            # physical PSF everywhere; add enhanced glow beyond the disk.
            hdr=physical+eta*(self.cfg['halo_display_gain']-1)*halo*self.display_halo_mask[...,None]
        else:
            hdr=(1-eta)*core+eta*self.cfg['halo_display_gain']*halo
        hdr=np.maximum(hdr[p:p+self.h,p:p+self.w],0)
        if self.cfg.get('tone_mapping')=='smooth-max-rgb':
            # A smooth common-channel shoulder avoids the abrupt gamut cap
            # that otherwise flattens a saturated golden photosphere.
            a=self.cfg['exposure']*hdr
            return a/(1+a.max(axis=-1))[...,None]
        if self.cfg.get('tone_mapping')=='luminance-preserving':
            from solar_surface import YROW
            y=hdr@YROW;l=-np.expm1(-self.cfg['exposure']*y)
            out=hdr*np.divide(l,y,out=np.zeros_like(y),where=y>1e-15)[...,None]
            # A common scale limits rare out-of-gamut highlights without
            # independent channel clipping or a yellow color multiplier.
            return out/np.maximum(1.,out.max(axis=-1))[...,None]
        return -np.expm1(-self.cfg['exposure']*hdr)
