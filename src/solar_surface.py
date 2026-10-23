"""Box-in-star reconstruction of actual evolving radiation-MHD intensities.
Six overlapping angular charts, periodic MHD domain, no procedural dynamics.
Spatial and temporal sampling assumptions are explicit in config and docs.
"""
from pathlib import Path
import numpy as np
from numba import njit,prange
from continuum import ROOT,MU,AZ
from geometry import orientation
YROW=np.linalg.inv(np.array([[3.2404542,-1.5371385,-.4985314],[-.969266,1.8760108,.041556],[.0556434,-.2040259,1.0572252]]))[1]

@njit(cache=True)
def interp_map(fields,it,im,ia,x,y,ch):
    ny,nx=fields.shape[3:5];ix=int(np.floor(x));iy=int(np.floor(y));fx=x-ix;fy=y-iy
    i0=ix%nx;i1=(ix+1)%nx;j0=iy%ny;j1=(iy+1)%ny
    return (1-fy)*((1-fx)*fields[it,im,ia,j0,i0,ch]+fx*fields[it,im,ia,j0,i1,ch])+fy*((1-fx)*fields[it,im,ia,j1,i0,ch]+fx*fields[it,im,ia,j1,i1,ch])

@njit(parallel=True,cache=True)
def sample_atlas_reference(normals,observer,mu,delay,fields,times,physical_time,magnification,mu_grid,center_mu_y,contrast):
    npts=len(normals);out=np.zeros((npts,3),np.float32);scale=695.7/magnification/(24/fields.shape[4])
    for j in prange(npts):
        n=normals[j];weights=np.abs(n)**24;weights/=weights.sum()
        for ax in range(3):
            if weights[ax]<1e-4:continue
            p=(ax+1)%3;q=(ax+2)%3;sign=1. if n[ax]>=0 else -1.
            face=2*ax+int(sign<0);u=np.arctan2(n[p],abs(n[ax]));v=np.arctan2(n[q],abs(n[ax]));rot=face%4
            uu=u;vv=v
            if rot==1:uu=-v;vv=u
            elif rot==2:uu=-u;vv=-v
            elif rot==3:uu=v;vv=-u
            x=uu*scale+face*37.;y=vv*scale+face*83.
            # Local orthonormal tangent basis. Chart metric distortion is a
            # documented spherical reconstruction approximation.
            eu=np.zeros(3);eu[p]=1.;eu-=n*n[p];eu/=np.sqrt((eu*eu).sum())
            ev=np.array([n[1]*eu[2]-n[2]*eu[1],n[2]*eu[0]-n[0]*eu[2],n[0]*eu[1]-n[1]*eu[0]])
            if ev[q]<0:ev=-ev
            ph=np.arctan2((ev*observer).sum(),(eu*observer).sum())+rot*np.pi/2
            af=(ph%(2*np.pi))/(np.pi/2);ia=int(af)%4;fa=af-np.floor(af)
            im=max(0,min(len(mu_grid)-2,np.searchsorted(mu_grid,mu[j])-1));fm=min(1.,max(0.,(mu[j]-mu_grid[im])/(mu_grid[im+1]-mu_grid[im])))
            # Same 40:1 physical clock everywhere; independent fixed local
            # phases avoid globally synchronizing the source box's p-modes.
            t=physical_time+face*60.-delay[j]
            if t<times[0] or t>times[-1]:
                # Guard is also checked in Python with the full time range.
                t=min(max(t,times[0]),times[-1])
            it=max(0,min(len(times)-2,np.searchsorted(times,t)-1));ft=(t-times[it])/(times[it+1]-times[it])
            rgb=np.zeros(3)
            for dt in range(2):
                wt=(1-ft) if dt==0 else ft
                for dm in range(2):
                    wm=(1-fm) if dm==0 else fm
                    for da in range(2):
                        wa=(1-fa) if da==0 else fa
                        for ch in range(3):rgb[ch]+=wt*wm*wa*interp_map(fields,it+dt,im+dm,(ia+da)%4,x,y,ch)
            # Fixed local contrast enhancement about a fixed mean limb law.
            # It preserves the calculated spatial/temporal phase and spectrum;
            # its nonlinear radiometry is explicitly a display operation.
            yy=np.dot(rgb,YROW);mean=np.interp(mu[j],mu_grid,center_mu_y)
            if yy>0:rgb*=max(0.,yy/mean)**(contrast-1)
            out[j]+=weights[ax]*rgb
    return out

@njit(cache=True,inline='always')
def _pixel_fast(fields,it,im,ia,i0,i1,j0,j1,fx,fy,ch):
    return (1-fy)*((1-fx)*fields[it,im,ia,j0,i0,ch]+fx*fields[it,im,ia,j0,i1,ch])+fy*((1-fx)*fields[it,im,ia,j1,i0,ch]+fx*fields[it,im,ia,j1,i1,ch])

@njit(parallel=True,cache=True)
def sample_atlas(normals,observer,mu,delay,fields,times,physical_time,magnification,mu_grid,center_mu_y,contrast):
    out=np.zeros((len(normals),3),np.float32);ny,nx=fields.shape[3:5];scale=695.7/magnification/(24/nx)
    for j in prange(len(normals)):
        n=normals[j];wx=abs(n[0])**24;wy=abs(n[1])**24;wz=abs(n[2])**24;ws=wx+wy+wz
        nod=n[0]*observer[0]+n[1]*observer[1]+n[2]*observer[2]
        for ax in range(3):
            weight=(wx if ax==0 else wy if ax==1 else wz)/ws
            if weight<1e-4:continue
            p=(ax+1)%3;q=(ax+2)%3;sgn=1. if n[ax]>=0 else -1.;face=2*ax+int(sgn<0)
            u=np.arctan2(n[p],abs(n[ax]));v=np.arctan2(n[q],abs(n[ax]));rot=face%4;uu=u;vv=v
            if rot==1:uu=-v;vv=u
            elif rot==2:uu=-u;vv=-v
            elif rot==3:uu=v;vv=-u
            x=uu*scale+face*37.;y=vv*scale+face*83.
            ix=int(np.floor(x));iy=int(np.floor(y));fx=x-ix;fy=y-iy;i0=ix%nx;i1=(ix+1)%nx;j0=iy%ny;j1=(iy+1)%ny
            norm=np.sqrt(max(1e-30,1-n[p]*n[p]))
            ou=(observer[p]-n[p]*nod)/norm
            ov=sgn*(n[ax]*observer[q]-n[q]*observer[ax])/norm
            ph=np.arctan2(ov,ou)+rot*np.pi/2
            af=(ph%(2*np.pi))/(np.pi/2);ia=int(af)%4;fa=af-np.floor(af)
            im=max(0,min(len(mu_grid)-2,np.searchsorted(mu_grid,mu[j])-1));fm=min(1.,max(0.,(mu[j]-mu_grid[im])/(mu_grid[im+1]-mu_grid[im])))
            t=min(max(physical_time+face*60.-delay[j],times[0]),times[-1]);it=max(0,min(len(times)-2,np.searchsorted(times,t)-1));ft=(t-times[it])/(times[it+1]-times[it])
            rr=0.;gg=0.;bb=0.
            for dt in range(2):
                wt=(1-ft) if dt==0 else ft
                for dm in range(2):
                    wm=(1-fm) if dm==0 else fm
                    for da in range(2):
                        wa=(1-fa) if da==0 else fa;f=wt*wm*wa;ind=(ia+da)%4
                        rr+=f*_pixel_fast(fields,it+dt,im+dm,ind,i0,i1,j0,j1,fx,fy,0)
                        gg+=f*_pixel_fast(fields,it+dt,im+dm,ind,i0,i1,j0,j1,fx,fy,1)
                        bb+=f*_pixel_fast(fields,it+dt,im+dm,ind,i0,i1,j0,j1,fx,fy,2)
            yy=rr*YROW[0]+gg*YROW[1]+bb*YROW[2];mean=np.interp(mu[j],mu_grid,center_mu_y)
            gain=max(0.,yy/mean)**(contrast-1) if yy>0 else 1.
            out[j,0]+=weight*gain*rr;out[j,1]+=weight*gain*gg;out[j,2]+=weight*gain*bb
    return out


class Surface:
    def __init__(self,cfg):
        self.cfg=cfg;self.files=sorted(p for p in (ROOT/'data/continuum').glob('*-17-2-solar.npz') if p.with_suffix('.json').exists())
        if len(self.files)<2:raise RuntimeError('Continuum sequence not ready')
        a=[np.load(p) for p in self.files];self.fields=np.stack([d['rgb'] for d in a]);self.times=np.array([float(d['time_s']) for d in a])
        self.mu=a[0]['mu'];self.mean_mu_y=(self.fields[0]@YROW).mean(axis=(1,2,3));self.normalization=self.mean_mu_y[-1]
        self.mean_mu_y/=self.normalization;self.fields/=self.normalization
        self.min_time=float(self.times[0]);self.max_time=float(self.times[-1]);self.star_redshift=(1-2*6.67430e-11*1.98847e30/(6.957e8*299792458.**2))**.5
    def ensure(self,maximum_time):
        """Allow the renderer to consume completed transfer maps as they arrive.
        Initial normalization stays fixed; no atmosphere or physical evolution
        is changed by this I/O optimization.
        """
        if maximum_time<=self.max_time:return
        import time
        start=time.monotonic();print('Waiting for source time',round(maximum_time,2),'available',self.max_time,flush=True)
        while maximum_time>self.max_time:
            files=sorted(p for p in (ROOT/'data/continuum').glob('*-17-2-solar.npz') if p.with_suffix('.json').exists())
            if files!=self.files:
                a=[np.load(p) for p in files]
                self.fields=np.stack([d['rgb'] for d in a])/self.normalization
                self.times=np.array([float(d['time_s']) for d in a]);self.files=files;self.min_time=float(self.times[0]);self.max_time=float(self.times[-1])
            if time.monotonic()-start>1800:raise RuntimeError('No required completed radiation map after 30 minutes')
            if maximum_time>self.max_time:time.sleep(3)
        print('Source available through',self.max_time,flush=True)
    def sample(self,samples,t,inc,az):
        if not len(samples.mu):return np.zeros((0,3),np.float32)
        pt=self.cfg['source_start_s']+t*self.cfg['physical_seconds_per_display_second']
        if pt-samples.delay.max()<self.min_time-1e-6 or pt+300>self.max_time+1e-6:raise ValueError('Out of released time sequence')
        rot=t*self.cfg['physical_seconds_per_display_second']/(25.38*86400)*360
        n=orientation(samples.camera,inc,az-rot);o=orientation(np.array([0.,0.,1.]),inc,az-rot)
        return sample_atlas(n,o,samples.mu,samples.delay,self.fields,self.times,pt,self.cfg['spatial_magnification'],self.mu,self.mean_mu_y,self.cfg['granulation_contrast_exponent'])*self.star_redshift**4
