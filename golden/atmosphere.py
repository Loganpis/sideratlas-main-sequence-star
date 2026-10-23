"""Optically thin projection of the actual 2.5D GOL radiation-MHD sequence.

The invariant third dimension is given a disclosed finite 6 Mm window. This
is a 2.5D-volume display approximation, NOT newly calculated 3D jet structure.
Positive conservative voxel projection integrates j dV / projected pixel area.
An opaque sphere occults the far side. No particles or field lines are drawn.
"""
from pathlib import Path
import json
import numpy as np
from numba import njit
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'src'))
from geometry import orientation,camera_direction

@njit(cache=True)
def project(fields,times,physical_time,x,z,dz,centers,easts,norths,basis,mag,depth,ns,cx,cy,rpix,width,height):
    image=np.zeros((height,width),np.float64);R=695.7
    dx=x[1]-x[0];dy=depth/ns;pixarea=(R/rpix)**2;roll=np.deg2rad(-20)
    cr=np.cos(roll);sr=np.sin(roll)
    for patch in range(len(centers)):
        center=centers[patch];east=easts[patch];north=norths[patch]
        for k in range(ns):
            y=(k+.5)*dy-depth/2
            wy=1. if abs(y)<depth/2-1 else .5*(1+np.cos(np.pi*(abs(y)-(depth/2-1))))
            v=y*mag/R;tv=np.tan(v);sv=1+tv*tv
            for ix in range(len(x)):
                u=(x[ix]-(x[0]+x[-1])/2)*mag/R;tu=np.tan(u);su=1+tu*tu;dn=np.sqrt(1+tu*tu+tv*tv)
                wx=min(1.,max(0.,min(x[ix]-x[0],x[-1]-x[ix])/2.));wx=.5-.5*np.cos(np.pi*wx)
                # Inverse camera basis; source frame is fixed to the star.
                nx=(center[0]+tu*east[0]+tv*north[0])/dn
                ny=(center[1]+tu*east[1]+tv*north[1])/dn
                nz=(center[2]+tu*east[2]+tv*north[2])/dn
                a=nx*basis[0,0]+ny*basis[1,0]+nz*basis[2,0]
                b=nx*basis[0,1]+ny*basis[1,1]+nz*basis[2,1]
                c=nx*basis[0,2]+ny*basis[1,2]+nz*basis[2,2]
                sx=cr*a-sr*b;sy=sr*a+cr*b
                jac=su*sv/(dn**3)
                for iz in range(len(z)):
                    rr=R+mag*z[iz];px=sx*rr;py=sy*rr;pz=c*rr
                    if pz<0 and px*px+py*py<R*R:continue
                    emitted_time=physical_time-(R-pz)*1e6/299792458.
                    if emitted_time<times[0] or emitted_time>times[-1]:
                        raise ValueError('Visible emission falls outside available source times')
                    it=max(0,min(len(times)-2,np.searchsorted(times,emitted_time)-1))
                    ft=(emitted_time-times[it])/(times[it+1]-times[it])
                    em=(1-ft)*fields[it,iz,ix]+ft*fields[it+1,iz,ix]
                    if em<=0:continue
                    xx=cx+px/R*rpix-.5;yy=cy-py/R*rpix-.5
                    xp=int(np.floor(xx));yp=int(np.floor(yy));fx=xx-xp;fy=yy-yp
                    if xp<0 or yp<0 or xp>=width-1 or yp>=height-1:continue
                    top=1. if z[iz]<12 else max(0.,.5+.5*np.cos(np.pi*(z[iz]-12)/2))
                    val=em*mag**2*(rr/R)**2*jac*dx*dy*dz[iz]*1e8/pixarea*wx*wy*top
                    image[yp,xp]+=val*(1-fx)*(1-fy);image[yp,xp+1]+=val*fx*(1-fy)
                    image[yp+1,xp]+=val*(1-fx)*fy;image[yp+1,xp+1]+=val*fx*fy
    return image

class Atmosphere:
    def __init__(self,cfg):
        self.cfg=cfg;self.a=cfg['atmosphere'];self.paths=sorted((HERE/'data/spicules-uv').glob('*.npz'))
        if not self.paths:raise RuntimeError('Fetch the actual spicule sequence first')
        self.times=[];self.fields=[]
        for p in self.paths:
            with np.load(p) as d:
                self.times.append(float(d['time_s']));self.fields.append(d['jUV']);self.x=d['x_Mm'].astype(float);self.z=d['z_Mm'].astype(float);self.dz=np.diff(d['z_edges_Mm']).astype(float)
        self.times=np.array(self.times);self.fields=np.stack(self.fields)
        self.centers=[];self.easts=[];self.norths=[]
        inc,az=camera_direction(0,cfg)
        for p in self.a['patch_camera_vectors_at_start']:
            center=orientation(np.array(p)/np.linalg.norm(p),inc,az)
            # Source x spans the limb, source invariant y initially along LOS.
            east=orientation(np.array([-p[1],p[0],0.])/np.hypot(p[0],p[1]),inc,az)
            north=np.cross(center,east)
            self.centers.append(center);self.easts.append(east);self.norths.append(north)
        self.centers=np.array(self.centers);self.easts=np.array(self.easts);self.norths=np.array(self.norths)
    def image(self,t,depth_samples=None):
        a=self.a;c=self.cfg;pt=a['source_start_s']+t*c['physical_seconds_per_display_second']
        if pt<self.times[0] or pt>self.times[-1]:raise ValueError('Missing GOL source time')
        it=max(0,min(len(self.times)-2,np.searchsorted(self.times,pt)-1));f=(pt-self.times[it])/(self.times[it+1]-self.times[it])
        inc,az=camera_direction(t/c['full_reel_seconds'],c);rot=t*c['physical_seconds_per_display_second']/(25.38*86400)*360
        basis=orientation(np.eye(3),inc,az-rot).T
        # orientation(I) returns transformed row vectors; transpose gives the
        # usual column transform. To obtain camera coordinates use dot columns.
        return project(self.fields,self.times,pt,self.x,self.z,self.dz,self.centers,self.easts,self.norths,basis,a['spatial_magnification'],a['invariant_depth_Mm'],depth_samples or a['depth_samples'],c['width']*c['center_fraction'][0],c['height']*c['center_fraction'][1],c['width']*c['radius_fraction_of_width'],c['width'],c['height'])
