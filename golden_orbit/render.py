"""Final golden, four-orbit display of the existing MHD/radiation calculations.

All colour operations are fixed presentation transforms, never new emission.
World geometry, source clocks and spectral intensities are inherited unchanged.
"""
from pathlib import Path
import sys,json,time,argparse
import numpy as np
from PIL import Image
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path[:0]=[str(ROOT/'golden'),str(ROOT/'photosphere'),str(ROOT),str(ROOT/'src')]
import surface,stellar_renderer as base,active_region
base.Surface=surface.ReconstructedSurface
active_region.ActiveRegion=surface.CoherentActiveRegion
from solar_surface import YROW
from geometry import srgb_to_linear,linear_to_srgb
from video_export import encoder,metadata
from atmosphere import Atmosphere

def grade_photosphere(physical,reference,cfg):
    """Continuous fixed colour LUT with exactly preserved pre-tone CIE Y."""
    y=physical@YROW
    relative=np.divide(y,reference,out=np.zeros_like(y),where=reference>0)
    knots=np.array(cfg['relative_intensity_knots'])
    gains=np.array(cfg['linear_rgb_gains_at_knots'])
    grade=np.stack([np.interp(relative,knots,gains[:,ch]) for ch in range(3)],axis=-1)
    rgb=physical*grade;gy=rgb@YROW
    rgb*=np.divide(y,gy,out=np.zeros_like(y),where=gy>0)[...,None]
    return rgb

class Renderer(base.Renderer):
    def __init__(self,cfg):
        super().__init__(cfg);self.atmosphere=Atmosphere(cfg);self.metrics=[]
        yy,xx=np.mgrid[:cfg['height'],:cfg['width']]
        cx,cy=np.array(cfg['center_fraction'])*[cfg['width'],cfg['height']]
        self.rho=np.hypot(xx+.5-cx,yy+.5-cy)/(cfg['width']*cfg['radius_fraction_of_width'])
        mu=np.sqrt(np.maximum(0,1-(self.rho*self.star.redshift_factor)**2))
        # An initial, fixed calculated limb law, not a per-frame normalization.
        self.grade_reference=np.interp(mu,self.surface.mu,self.surface.mean_mu_y)*self.surface.star_redshift**4
        a=cfg['atmosphere']
        self.uv_visibility=np.minimum(a['radial_display_gain_cap'],np.exp(np.maximum(0,self.rho-1)/a['radial_display_e_fold_R']))
    def title(self):
        from PIL import ImageDraw,ImageFont
        c=self.cfg;s=c['width']/1080;ss=3
        im=Image.new('L',(c['width']*ss,c['height']*ss));d=ImageDraw.Draw(im)
        for x,y,v,size,name,anchor in [(84,250,'@sideratlas',34,'lmroman10-regular.otf','lt'),(996,253,'NO.004',25,'lmroman10-regular.otf','rt'),(540,320,'THE MAIN-SEQUENCE STAR',58,'lmroman12-bold.otf','mt')]:
            f=ImageFont.truetype(str(ROOT/'fonts'/name),round(size*s*ss))
            d.text((x*s*ss,y*s*ss),v,font=f,fill=255,anchor=anchor)
        return np.asarray(im.resize((c['width'],c['height']),Image.Resampling.LANCZOS),'f4')/255
    def graded_source(self,t):
        physical=super().source(t)
        graded=grade_photosphere(physical,self.grade_reference,self.cfg['color_grade'])
        error=float(np.max(abs(graded@YROW-physical@YROW)))
        uv=self.atmosphere.image(t);a=self.cfg['atmosphere']
        channel=a['display_gain']*np.arcsinh(uv*self.uv_visibility/a['display_asinh_width_erg_cm2_s_sr'])/np.arcsinh(a['intensity_reference_erg_cm2_s_sr']/a['display_asinh_width_erg_cm2_s_sr'])
        color=np.array(a['false_color_linear_rgb']);color/=color@YROW
        self.metrics.append(dict(time_s=t,grade_max_absolute_Y_error=error,physical_UV_sum=float(uv.sum()),physical_UV_max=float(uv.max()),display_UV_max=float(channel.max())))
        return graded+channel[...,None]*color
    def frame(self,index):
        t=self.cfg.get('time_offset_s',0)+(index+.5)/self.cfg['fps']
        im=self.optics.process(self.graded_source(t))
        im=im*(1-self.alpha)+self.alpha*float(srgb_to_linear(np.array(244/255)))
        assert np.isfinite(im).all() and (im>=0).all() and (im<=1).all()
        return np.uint16(np.clip(linear_to_srgb(im),0,1)*65535+.5)

def main():
    p=argparse.ArgumentParser();p.add_argument('--config',default='config-preview.json')
    p.add_argument('--still',action='store_true');p.add_argument('--time',type=float,default=0);args=p.parse_args()
    c=json.loads((HERE/args.config).read_text());r=Renderer(c);start=time.monotonic()
    stem='sideratlas-no004-golden-orbit-32s' if c['seconds']==32 else 'sideratlas-no004-orbit-preview'
    (HERE/'output').mkdir(exist_ok=True);(HERE/'validation').mkdir(exist_ok=True)
    if args.still:
        im=r.frame(round(args.time*c['fps']))
        Image.fromarray((im//257).astype('u1')).save(HERE/f'output/still-{args.time:g}.png')
        print('Still',time.monotonic()-start,r.metrics[-1],flush=True);return
    out=HERE/f'output/{stem}.mp4';n=round(c['seconds']*c['fps']);enc=encoder(out,c['width'],c['height'],c['fps'])
    for i in range(n):
        frame=r.frame(i);enc.stdin.write(frame.tobytes())
        if i in [0,n//2,n-1]:Image.fromarray((frame//257).astype('u1')).save(HERE/f'output/frame-{i:03}.png')
        if i%15==0:print(i,'/',n,'elapsed',round(time.monotonic()-start,1),'s',flush=True)
    enc.stdin.close();assert enc.wait()==0
    Image.open(HERE/'output/frame-000.png').save(HERE/'output/cover.png')
    (HERE/f'validation/{stem}-render.json').write_text(json.dumps(dict(config=c,frames=n,elapsed_s=time.monotonic()-start,max_normalized_PSF_flux_error=max(r.optics.flux_errors),minimum_crop_fraction=min(r.optics.crop_fractions),emission_metrics=r.metrics,export=metadata(),GOL_source_times_s=r.atmosphere.times.tolist()),indent=2)+'\n')
    print('Completed',out,flush=True)

if __name__=='__main__':main()
