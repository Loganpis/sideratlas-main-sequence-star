"""Non-synchronous box-in-star radiation atlas; no procedural fluid evolution.

Static tile phases/translations/quarter turns sample the existing local MHD
sequence. Compact overlap windows soften joins, with no claimed global MHD
continuity. The complete magnetic region uses a separate coherent source map.
"""
import sys
from pathlib import Path
import numpy as np
from numba import njit,prange
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'src'))
from solar_surface import Surface,YROW,_pixel_fast
from geometry import orientation
from active_region import ActiveRegion

@njit(cache=True)
def _window(x):
    x=abs(x)
    if x<=7.2:return 1.
    if x>=12:return 0.
    return .5*(1+np.cos(np.pi*(x-7.2)/4.8))

@njit(cache=True)
def _hash(face,x,y):
    # A deterministic label for a rigid tile transform, never evolving noise.
    h=((x+32768)*73856093)^((y+32768)*19349663)^((face+41004)*83492791)
    return h&2147483647

@njit(parallel=True,cache=True)
def sample_patches(normals,observer,mu,delay,fields,times,physical_time,mag,mu_grid,mean_y,contrast):
    out=np.zeros((len(normals),3),np.float32);ny,nx=fields.shape[3:5];dx=24/nx
    for j in prange(len(normals)):
        n=normals[j];weights=np.abs(n)**24;weights/=weights.sum()
        rgb=np.zeros(3)
        for ax in range(3):
            if weights[ax]<1e-4:continue
            p=(ax+1)%3;q=(ax+2)%3;sgn=1. if n[ax]>=0 else -1.;face=2*ax+int(sgn<0)
            u=np.arctan2(n[p],abs(n[ax]))*695.7/mag;v=np.arctan2(n[q],abs(n[ax]))*695.7/mag
            norm=np.sqrt(max(1e-30,1-n[p]*n[p]));nod=np.dot(n,observer)
            ou=(observer[p]-n[p]*nod)/norm;ov=sgn*(n[ax]*observer[q]-n[q]*observer[ax])/norm
            ph=np.arctan2(ov,ou)
            cellx=int(np.floor(u/19.2));celly=int(np.floor(v/19.2));acc=np.zeros(3);total=0.
            im=max(0,min(len(mu_grid)-2,np.searchsorted(mu_grid,mu[j])-1));fm=min(1.,max(0.,(mu[j]-mu_grid[im])/(mu_grid[im+1]-mu_grid[im])))
            for tx in range(cellx,cellx+2):
                du=u-tx*19.2;wu=_window(du)
                if wu<=0:continue
                for ty in range(celly,celly+2):
                    dv=v-ty*19.2;ww=wu*_window(dv)
                    if ww<=0:continue
                    h=_hash(face,tx,ty);rot=h%4;uu=du;vv=dv
                    if rot==1:uu=-dv;vv=du
                    elif rot==2:uu=-du;vv=-dv
                    elif rot==3:uu=dv;vv=-du
                    x=(uu+12)/dx+((h//4)%nx);y=(vv+12)/dx+((h//1024)%ny)
                    ix=int(np.floor(x));iy=int(np.floor(y));fx=x-ix;fy=y-iy;i0=ix%nx;i1=(ix+1)%nx;j0=iy%ny;j1=(iy+1)%ny
                    t=physical_time+((h//1048576)%6)*60.-delay[j]
                    it=max(0,min(len(times)-2,np.searchsorted(times,t)-1));ft=(t-times[it])/(times[it+1]-times[it])
                    af=((ph+rot*np.pi/2)%(2*np.pi))/(np.pi/2);ia=int(af)%4;fa=af-np.floor(af)
                    for dt in range(2):
                        wt=(1-ft) if dt==0 else ft
                        for dm in range(2):
                            wm=(1-fm) if dm==0 else fm
                            for da in range(2):
                                wa=(1-fa) if da==0 else fa;f=ww*wt*wm*wa;ind=(ia+da)%4
                                for ch in range(3):acc[ch]+=f*_pixel_fast(fields,it+dt,im+dm,ind,i0,i1,j0,j1,fx,fy,ch)
                    total+=ww
            if total>0:rgb+=weights[ax]*acc/total
        yy=np.dot(rgb,YROW);mean=np.interp(mu[j],mu_grid,mean_y)
        gain=max(0.,yy/mean)**(contrast-1) if yy>0 else 1.
        out[j]=rgb*gain
    return out

class ReconstructedSurface(Surface):
    def sample(self,samples,t,inc,az):
        if not len(samples.mu):return np.zeros((0,3),np.float32)
        pt=self.cfg['source_start_s']+t*self.cfg['physical_seconds_per_display_second']
        if pt-samples.delay.max()<self.min_time or pt+300>self.max_time:raise ValueError('Missing required physical source time')
        rot=t*self.cfg['physical_seconds_per_display_second']/(25.38*86400)*360
        n=orientation(samples.camera,inc,az-rot);o=orientation(np.array([0.,0.,1.]),inc,az-rot)
        return sample_patches(n,o,samples.mu,samples.delay,self.fields,self.times,pt,self.cfg['spatial_magnification'],self.mu,self.mean_mu_y,self.cfg['granulation_contrast_exponent'])*self.star_redshift**4

@njit(parallel=True,cache=True)
def active_rectangle(n,delay,time,center,east,north,fields,large_scale,times,mu,ld_mu,ld_y,colors,tempgrid,mag,quiet,contrast):
    out=np.zeros((len(n),3),np.float32);weight=np.zeros(len(n),np.float32)
    from_x=49.152;from_y=24.576
    for j in prange(len(n)):
        front=np.dot(n[j],center)
        if front<=0:continue
        u=np.arctan2(np.dot(n[j],east),front)*695.7/mag;v=np.arctan2(np.dot(n[j],north),front)*695.7/mag
        if abs(u)>=48.960 or abs(v)>=24.384:continue
        wx=1. if abs(u)<=39.552 else .5*(1+np.cos(np.pi*(abs(u)-39.552)/(48.960-39.552)))
        wy=1. if abs(v)<=17.472 else .5*(1+np.cos(np.pi*(abs(v)-17.472)/(24.384-17.472)))
        x=(u+from_x)/.192;y=(v+from_y)/.192
        ix=int(np.floor(x));iy=int(np.floor(y));fx=x-ix;fy=y-iy
        t=time-delay[j];it=max(0,min(len(times)-2,np.searchsorted(times,t)-1));ft=(t-times[it])/(times[it+1]-times[it])
        vals=np.zeros(2)
        for k in range(2):
            array=fields if k==0 else large_scale
            for dt in range(2):
                a=array[it+dt];val=(1-fy)*((1-fx)*a[iy,ix]+fx*a[iy,ix+1])+fy*((1-fx)*a[iy+1,ix]+fx*a[iy+1,ix+1])
                vals[k]+=((1-ft) if dt==0 else ft)*val/quiet
        ratio=vals[0];baseline=vals[1];tb=5772*ratio**.25
        k=max(0,min(len(tempgrid)-2,np.searchsorted(tempgrid,tb)-1));f=min(1.,max(0.,(tb-tempgrid[k])/(tempgrid[k+1]-tempgrid[k])))
        chrom=(1-f)*colors[k]+f*colors[k+1]
        out[j]=np.interp(mu[j],ld_mu,ld_y)*baseline*max(0.,ratio/baseline)**contrast*chrom
        weight[j]=wx*wy
    return out,weight

class CoherentActiveRegion(ActiveRegion):
    def sample(self,samples,t,inc,az):
        if not len(samples.mu):return np.zeros((0,3),np.float32),np.zeros(0,np.float32)
        physical=t*self.clock;rotation=physical/(25.38*86400)*360
        n=orientation(samples.camera,inc,az-rotation)
        if self.start+physical>self.times[-1]:raise ValueError('MURaM sequence too short')
        return active_rectangle(n,samples.delay,self.start+physical,self.center,self.east,self.north,self.intensity,self.large_scale,self.times,samples.mu,self.surface.mu,self.surface.mean_mu_y,self.colors,self.tempgrid,self.cfg['spatial_magnification'],self.quiet,self.cfg['contrast_exponent'])
