"""Decode every frame; verify the final delivery and its single-generation export."""
from pathlib import Path
import sys,json,subprocess,hashlib,argparse
import numpy as np
from PIL import Image
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path.insert(0,str(ROOT/'src'))
from video_export import DECODE_FILTER

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--preview',action='store_true');a=ap.parse_args()
    c=json.loads((HERE/('config-preview.json' if a.preview else 'config-final.json')).read_text());stem='sideratlas-no004-orbit-preview' if a.preview else 'sideratlas-no004-golden-orbit-32s';movie=HERE/f'output/{stem}.mp4'
    p=json.loads(subprocess.check_output(['ffprobe','-v','error','-show_streams','-show_format','-of','json',str(movie)]));s=p['streams'][0];n=round(c['fps']*c['seconds'])
    assert int(s['nb_frames'])==n and s['width']==c['width'] and s['height']==c['height']
    assert abs(float(p['format']['duration'])-c['seconds'])<.001
    assert s['color_space']==s['color_transfer']==s['color_primaries']=='bt709' and s['color_range']=='tv' and s['pix_fmt']=='yuv420p'
    q=subprocess.Popen(['ffmpeg','-v','error','-xerror','-i',str(movie),'-vf',DECODE_FILTER,'-f','rawvideo','-pix_fmt','rgb24','-'],stdout=subprocess.PIPE)
    size=c['width']*c['height']*3;last=None;diff=[];frames=0
    while True:
        b=q.stdout.read(size)
        if not b:break
        assert len(b)==size
        image=np.frombuffer(b,'u1').reshape(c['height'],c['width'],3).astype(float)
        if frames==0:
            ref=np.asarray(Image.open(HERE/'output/frame-000.png')).astype(float);err=image-ref;mae=float(abs(err).mean());psnr=float(10*np.log10(255**2/np.mean(err**2)))
        if last is not None:diff.append(float(abs(image-last).mean()))
        last=image;frames+=1
    assert q.wait()==0 and frames==n and min(diff)>0
    render=json.loads((HERE/f'validation/{stem}-render.json').read_text());assert render['max_normalized_PSF_flux_error']<1e-4
    rec=dict(verified=True,file=movie.name,sha256=hashlib.sha256(movie.read_bytes()).hexdigest(),size_bytes=movie.stat().st_size,width=s['width'],height=s['height'],frames=n,duration_s=float(p['format']['duration']),complete_decode=True,first_frame_MAE_8bit=mae,first_frame_PSNR_dB=psnr,minimum_adjacent_frame_MAE_8bit=min(diff),max_normalized_PSF_flux_error=render['max_normalized_PSF_flux_error'],minimum_crop_fraction=render['minimum_crop_fraction'],only_labels=['THE MAIN-SEQUENCE STAR','@sideratlas','NO.004'])
    (HERE/f'validation/{stem}-export.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))

if __name__=='__main__':main()
