"""High-depth input, single-generation Instagram H.264 export."""
from pathlib import Path
import shutil,subprocess

FILTER=(
    'format=gbrp16le,'
    'zscale=primariesin=709:transferin=iec61966-2-1:matrixin=gbr:rangein=full:'
    'primaries=709:transfer=709:matrix=709:range=limited:'
    'chromal=left:filter=spline36:agamma=0:dither=error_diffusion,'
    'format=yuv420p'
)
DECODE_FILTER=(
    'zscale=primariesin=709:transferin=709:matrixin=709:rangein=limited:chromalin=left:'
    'primaries=709:transfer=iec61966-2-1:matrix=gbr:range=full:'
    'filter=spline36:agamma=0:dither=none,format=gbrp16le,format=rgb24'
)


def encoder(path,width,height,fps):
    ffmpeg=shutil.which('ffmpeg')
    if not ffmpeg:raise RuntimeError('FFmpeg with libx264 and zscale required')
    command=[ffmpeg,'-hide_banner','-loglevel','error','-y','-f','rawvideo','-pix_fmt','rgb48le',
        '-s',f'{width}x{height}','-framerate',str(fps),'-i','-','-an','-vf',FILTER,
        '-c:v','libx264','-threads','4','-preset','slow','-crf','14','-maxrate','20M','-bufsize','40M',
        '-profile:v','high','-pix_fmt','yuv420p','-r',str(fps),'-fps_mode','cfr',
        '-g',str(2*fps),'-keyint_min',str(2*fps),'-bf','0','-sc_threshold','0','-flags','+cgop',
        '-color_primaries','bt709','-color_trc','bt709','-colorspace','bt709','-color_range','tv',
        '-chroma_sample_location','left','-bsf:v',
        'h264_metadata=video_full_range_flag=0:colour_primaries=1:transfer_characteristics=1:matrix_coefficients=1',
        '-movflags','+faststart','-use_editlist','0',str(path)]
    return subprocess.Popen(command,stdin=subprocess.PIPE)


def metadata():
    return dict(input='16-bit-per-channel sRGB RGB48LE',output='H.264 High / yuv420p / Rec.709 limited / MP4',
        crf=14,maxrate_bps=20000000,preset='slow',encoder_threads=4,audio=False,faststart=True,
        quantization='zscale error diffusion at final 8-bit conversion, not astrophysical texture or added decorative grain',
        filter=FILTER,ffmpeg_version=subprocess.check_output(['ffmpeg','-version'],text=True).splitlines()[0])
