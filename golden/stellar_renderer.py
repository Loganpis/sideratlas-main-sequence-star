"""Shared solar optical renderer; final assembly is in render.py."""
import sys,json,time,argparse
from pathlib import Path
import numpy as np
from PIL import Image,ImageFont,ImageDraw
ROOT=Path(__file__).resolve().parent.parent;sys.path.insert(0,str(ROOT/'src'))
from geometry import SolarStar,Samples,Tile,ray_map,camera_direction,linear_to_srgb,srgb_to_linear
from solar_surface import Surface
from optics import Optics
from video_export import encoder,metadata

class Renderer:
    def __init__(self,cfg):
        self.cfg=cfg;self.star=SolarStar();self.surface=Surface(cfg)
        self.active=None
        if 'active_region' in cfg:
            from active_region import ActiveRegion
            self.active=ActiveRegion(cfg,self.surface)
        grid=np.linspace(0,1,4097);cp,dt=ray_map(grid,self.star.compactness);rays=(grid,cp,dt)
        self.tiles=[Tile(y,min(y+cfg['tile_rows'],cfg['height']),cfg,self.star,rays) for y in range(0,cfg['height'],cfg['tile_rows'])]
        for tile in self.tiles:
            if not np.all(tile.fallback_samples.mask):
                ey,ex=np.nonzero(tile.edge);cx,cy=np.array(cfg['center_fraction'])*[cfg['width'],cfg['height']];rad=cfg['radius_fraction_of_width']*cfg['width']
                x=(ex[tile.fallback]+.5-cx)/rad;y=-(ey[tile.fallback]+tile.row0+.5-cy)/rad
                f=(1-1e-12)/(self.star.redshift_factor*np.hypot(x,y));tile.fallback_samples=Samples(x*f,y*f,self.star,rays)
        self.optics=Optics((cfg['height'],cfg['width']),cfg['optics']);self.alpha=self.title()[...,None]
    def title(self):
        c=self.cfg;s=c['width']/1080;ss=3;im=Image.new('L',(c['width']*ss,c['height']*ss));d=ImageDraw.Draw(im)
        def text(x,y,value,size,bold=False,anchor='lt'):
            name='lmroman12-bold.otf' if bold else 'lmroman10-regular.otf'
            font=ImageFont.truetype(str(ROOT/'fonts'/name),round(size*s*ss))
            while value=='THE MAIN-SEQUENCE STAR' and d.textlength(value,font=font)>936*s*ss:
                size-=1;font=ImageFont.truetype(str(ROOT/'fonts'/name),round(size*s*ss))
            d.text((x*s*ss,y*s*ss),value,font=font,fill=255,anchor=anchor)
        text(72,270,'@sideratlas',34);text(1008,275,'NO.004',25,anchor='rt');text(540,335,'THE MAIN-SEQUENCE STAR',70,True,'mt')
        return np.asarray(im.resize((c['width'],c['height']),Image.Resampling.LANCZOS),'f4')/255
    def source(self,t):
        c=self.cfg
        if c.get('stream_source_maps',False):self.surface.ensure(c['source_start_s']+t*c['physical_seconds_per_display_second']+300)
        inc,az=camera_direction(t/c['full_reel_seconds'],c);image=np.zeros((c['height'],c['width'],3),np.float32)
        for tile in self.tiles:
            colors=[]
            for samples in [tile.regular,tile.edge_samples,tile.fallback_samples]:
                c=self.surface.sample(samples,t,inc,az)
                if self.active is not None:
                    a,w=self.active.sample(samples,t,inc,az);c=c*(1-w[:,None])+a*w[:,None]
                colors.append(c)
            image[tile.row0:tile.row1]=tile.project(*colors)
        return image
    def frame(self,index):
        t=(index+.5)/self.cfg['fps'];source=self.source(t);im=self.optics.process(source)
        im=im*(1-self.alpha)+self.alpha*float(srgb_to_linear(np.array(244/255)))
        assert np.isfinite(im).all() and (im>=0).all() and (im<=1).all()
        return np.uint16(np.clip(linear_to_srgb(im),0,1)*65535+.5)
