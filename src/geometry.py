"""SiderAtlas No.003 geometry reused; solar physics is supplied separately."""
import numpy as np
from pathlib import Path
def srgb_to_linear(x):return np.where(x<=.04045,x/12.92,((x+.055)/1.055)**2.4)

def linear_to_srgb(x):return np.where(x<=.0031308,12.92*x,1.055*x**(1/2.4)-.055)

def camera_direction(progress,config):
    if config['camera'].get('motion')=='uniform-orbit':
        # Keep unwrapped phase: four complete orbits with no easing or resets.
        # Physical stellar rotation belongs in the atmosphere, never here.
        p=np.asarray(progress)
        return (np.full_like(p,config['camera']['inclination_deg'],dtype=float),
                config['camera']['azimuth_start_deg']+360*config['camera']['revolutions']*p)
    f=(1-np.cos(np.pi*np.asarray(progress)))/2
    inc=config['camera']['inclination_deg'];az=config['camera']['azimuth_deg']
    return inc[0]+f*(inc[1]-inc[0]),az[0]+f*(az[1]-az[0])

def orientation(camera,inc,az):
    a,b,c=np.moveaxis(camera,-1,0);i,p=np.deg2rad([inc,az])
    x=np.cos(i)*a-np.sin(i)*c;z=np.sin(i)*a+np.cos(i)*c
    return np.stack((np.cos(p)*x-np.sin(p)*b,np.sin(p)*x+np.cos(p)*b,z),axis=-1)

def circle_rectangle_area(x0,x1,y0,y1,radius):
    """Analytic circle/square overlap, avoiding quantized limb coverage."""
    def primitive(x,y):
        sign=np.sign(x)*np.sign(y);a=np.minimum(abs(x),radius);b=np.minimum(abs(y),radius)
        t=np.sqrt(np.maximum(0,radius*radius-b*b));lo=np.minimum(a,t)
        def integral(s):
            return .5*(s*np.sqrt(np.maximum(0,radius*radius-s*s))+radius*radius*np.arcsin(s/radius))
        return sign*(b*lo+integral(a)-integral(lo))
    area=primitive(x1,y1)-primitive(x0,y1)-primitive(x1,y0)+primitive(x0,y0)
    return np.clip(area,0,(x1-x0)*(y1-y0))


from numpy.polynomial.legendre import leggauss
C=299792458.
def ray_map(mu_em, compactness, quadrature=768):
    """Direct Schwarzschild rays: cos(psi), delay/(R/c), in the static frame.

    b/R=sqrt(1-mu_em²)/sqrt(1-u). Substitution x=R/r=1-s² removes
    the tangent-ray endpoint singularity. The delay is relative to the
    radial ray; a common propagation distance and Shapiro offset cancel.
    """
    mu = np.asarray(mu_em, float)
    if np.any((mu < 0) | (mu > 1)) or not 0 <= compactness < 0.1:
        raise ValueError('Weakly compact star and outward rays required')
    nodes, weights = leggauss(quadrature)
    s = (nodes + 1) / 2
    w = weights / 2
    x = 1 - s * s
    b2 = (1 - mu * mu) / (1 - compactness)
    f = mu[..., None] ** 2 + b2[..., None] * s * s * (1 + x - compactness * (1 + x + x * x))
    if np.any(f <= 0):
        raise ValueError('Invalid direct-ray radicand')
    root = np.sqrt(f)
    psi = np.sqrt(b2) * np.sum(2 * s * w / root, axis=-1)
    delay = np.sum(2 * s * w * b2[..., None] / (root * (1 + root)), axis=-1)
    return (np.cos(psi), delay)

class SolarStar:
    radius_m=6.957e8
    mass_kg=1.98847e30
    compactness=2*6.67430e-11*mass_kg/(radius_m*C*C)
    redshift_factor=np.sqrt(1-compactness)

class Samples:
    def __init__(self,x,y,star,ray_grid):
        q2=(x*x+y*y)*star.redshift_factor**2
        self.mask=q2<=1
        x=x[self.mask];y=y[self.mask]
        self.mu=np.sqrt(np.maximum(0,1-q2[self.mask]))
        grid,cp,dt=ray_grid
        cp=np.clip(np.interp(self.mu,grid,cp),-1,1)
        self.delay=np.interp(self.mu,grid,dt)*star.radius_m/C
        sp=np.sqrt(np.maximum(0,1-cp*cp));rho=np.hypot(x,y)
        nx=np.divide(x,rho,out=np.zeros_like(x),where=rho!=0)*sp
        ny=np.divide(y,rho,out=np.zeros_like(y),where=rho!=0)*sp
        roll=np.deg2rad(-20)
        self.camera=np.stack((np.cos(roll)*nx+np.sin(roll)*ny,-np.sin(roll)*nx+np.cos(roll)*ny,cp),axis=-1)

class Tile:
    def __init__(self,row0,row1,config,star,ray_grid,samples_class=Samples):
        self.row0=row0;self.row1=row1;w=config['width'];ss=config['supersample']
        self.ss=ss;self.limb_ss=config['limb_supersample'];cx,cy=np.array(config['center_fraction'])*[w,config['height']]
        radius=config['radius_fraction_of_width']*w
        yy,xx=np.mgrid[row0*ss:row1*ss,:w*ss]
        x=((xx+.5)/ss-cx)/radius;y=-((yy+.5)/ss-cy)/radius
        self.regular=samples_class(x,y,star,ray_grid)
        # Only pixels intersecting the circular limb use an 8x8 subgrid.
        py,px=np.mgrid[row0:row1,:w]
        rho=np.hypot(px+.5-cx,py+.5-cy)
        self.edge=np.abs(rho-radius/star.redshift_factor)<=np.sqrt(.5)
        ey,ex=np.nonzero(self.edge);n=self.limb_ss
        apparent_radius=radius/star.redshift_factor
        self.edge_area=circle_rectangle_area(ex-cx,ex+1-cx,ey+row0-cy,ey+row0+1-cy,apparent_radius)
        dy,dx=np.mgrid[:n,:n]
        x=(ex[:,None]+(dx.ravel()+.5)/n-cx)/radius
        y=-(ey[:,None]+row0+(dy.ravel()+.5)/n-cy)/radius
        self.edge_samples=samples_class(x,y,star,ray_grid)
        self.edge_counts=self.edge_samples.mask.sum(axis=1)
        self.fallback=(self.edge_counts==0)&(self.edge_area>1e-10)
        # Extremely thin slivers may contain none of the 8x8 intensity nodes.
        # Pick an interior point halfway from the nearest pixel corner to
        # its intersection with the circle, retaining the exact overlap area.
        fx=ex[self.fallback];fy=ey[self.fallback]+row0
        p0=np.stack((np.clip(cx,fx,fx+1)-cx,np.clip(cy,fy,fy+1)-cy),axis=-1)
        direction=np.stack((fx+.5-cx,fy+.5-cy),axis=-1)-p0
        aa=np.sum(direction*direction,axis=1);bb=2*np.sum(p0*direction,axis=1)
        cc=np.sum(p0*p0,axis=1)-apparent_radius**2
        root=(-bb+np.sqrt(np.maximum(0,bb*bb-4*aa*cc)))/(2*aa)
        point=p0+.5*root[:,None]*direction
        self.fallback_samples=samples_class(point[:,0]/radius,-point[:,1]/radius,star,ray_grid)

    def project(self,regular_colors,edge_colors,fallback_colors):
        ss=self.ss;rows=self.row1-self.row0;w=self.edge.shape[1]
        image=np.zeros(self.regular.mask.shape+(3,),'f4')
        image[self.regular.mask]=regular_colors
        image=image.reshape(rows,ss,w,ss,3).mean(axis=(1,3))
        sparse=np.zeros(self.edge_samples.mask.shape+(3,),'f4')
        sparse[self.edge_samples.mask]=edge_colors
        mean=np.divide(sparse.sum(axis=1),self.edge_counts[:,None],out=np.zeros((len(self.edge_counts),3),'f4'),where=self.edge_counts[:,None]>0)
        mean[self.fallback]=fallback_colors
        image[self.edge]=mean*self.edge_area[:,None]
        return image
