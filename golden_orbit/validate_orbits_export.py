"""Check encoded frame order and continuity at the three orbit boundaries."""
from pathlib import Path
import sys,json,subprocess,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'src'))
from video_export import DECODE_FILTER
from geometry import camera_direction

def main():
    movie=HERE/'output/sideratlas-no004-golden-orbit-32s.mp4'
    c=json.loads((HERE/'config-final.json').read_text())
    render=json.loads((HERE/'validation/sideratlas-no004-golden-orbit-32s-render.json').read_text())
    assert c==render['config'] and render['frames']==960
    source_times=np.array([m['time_s'] for m in render['emission_metrics']])
    assert len(source_times)==960
    assert np.max(abs(source_times-(np.arange(960)+.5)/30))<1e-12
    max_grade_error=max(m['grade_max_absolute_Y_error'] for m in render['emission_metrics'])
    assert max_grade_error<1e-12
    probe=json.loads(subprocess.check_output(['ffprobe','-v','error','-select_streams','v:0','-show_frames','-show_streams','-of','json',str(movie)]))
    frames=probe['frames'];assert len(frames)==960
    pts=np.array([float(f['best_effort_timestamp_time']) for f in frames])
    assert np.max(abs(np.diff(pts)-1/30))<1.1e-6 and abs(pts[0])<1e-8
    assert probe['streams'][0]['avg_frame_rate']=='30/1'
    # Decode only for display continuity diagnostics. The exact camera phase
    # follows the tested renderer/configuration, not an image-motion estimate.
    w=270;h=480
    q=subprocess.Popen(['ffmpeg','-v','error','-xerror','-i',str(movie),'-vf',DECODE_FILTER+',scale=270:480','-f','rawvideo','-pix_fmt','rgb24','-'],stdout=subprocess.PIPE)
    n=w*h*3;previous=None;diff=[];quarter_starts=[];index=0
    yy,xx=np.mgrid[:h,:w];cx,cy=np.array(c['center_fraction'])*[w,h]
    mask=np.hypot(xx+.5-cx,yy+.5-cy)<w*c['radius_fraction_of_width']*1.16
    while True:
        b=q.stdout.read(n)
        if not b:break
        assert len(b)==n
        image=np.frombuffer(b,'u1').reshape(h,w,3).astype(float)
        if previous is not None:diff.append(float(abs(image-previous)[mask].mean()))
        if index in [0,240,480,720]:quarter_starts.append(image.copy())
        previous=image;index+=1
    assert q.wait()==0 and index==960
    boundaries=[]
    for frame in [240,480,720]:
        local=np.array(diff[frame-11:frame+10])
        value=diff[frame-1]
        boundaries.append(dict(new_orbit_frame=frame,time_s=frame/30,
            adjacent_frame_MAE_8bit=value,local_median_MAE_8bit=float(np.median(local)),
            ratio_to_local_median=float(value/np.median(local))))
        assert value<3*np.quantile(diff,.99)+1
    changing=[float(abs(quarter_starts[k+1]-quarter_starts[k])[mask].mean()) for k in range(3)]
    assert min(changing)>0
    ft=(np.arange(960)+.5)/30;inc,az=camera_direction(ft/32,c)
    assert np.max(abs(np.diff(az)-1.5))<1e-10
    result=dict(verified=True,sha256=hashlib.sha256(movie.read_bytes()).hexdigest(),
        encoded_frames=960,frames_per_orbit=240,revolutions=4,seconds_per_orbit=8,
        continuous_unwrapped_camera=True,uniform_encoded_timestamps=True,
        all_source_frame_times_monotonic_and_unreset=True,
        all_frames_max_pre_tone_CIE_Y_grade_error=max_grade_error,
        first_last_exposure_center_azimuth_deg=[float(az[0]),float(az[-1])],
        azimuth_advance_per_frame_deg=1.5,orbit_boundaries=boundaries,
        same_view_between_consecutive_orbits_MAE_8bit=changing,
        no_repeated_eight_second_clip=True,
        note='960 exposure centers represent a continuous 32-second, four-orbit interval; no duplicate endpoint. Only the camera is periodic. Plasma and source evolution remain continuous and are not forced to loop.')
    (HERE/'validation/orbit-export.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
