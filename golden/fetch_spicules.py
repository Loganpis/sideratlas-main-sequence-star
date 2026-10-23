"""Public Bifrost GOL spicule sequence -> conserved optically-thin C IV emissivity.

The source is 2.5D, not a 3D jet catalogue. No jet trajectories are synthesized.
Only compact emissivity maps and immutable byte hashes are retained, not the
large temporary thermodynamic slabs. Requires the recorded AIA response table.
"""
from pathlib import Path
import sys,json,hashlib,concurrent.futures,time
import numpy as np,requests
from scipy.io import readsav
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'src'))
from archive import header
BASE='https://sdc.uio.no/vol/simulations/en096014_gol/atmos/'
OUT=HERE/'data/spicules-uv';OUT.mkdir(parents=True,exist_ok=True)
from atomic import coefficient
logT,G=coefficient();logR=np.log10(np.maximum(G,1e-100))
rpath=HERE/'data/atomic/CIV-coefficient.npz'

def field(var,snap):
    url=BASE+f'BIFROST_en096014_gol_{var}_{snap}.fits'
    sess=requests.Session();sess.mount('https://',HTTPAdapter(max_retries=Retry(total=4,backoff_factor=1)))
    def get(a,b):
        x=sess.get(url,headers={'Range':f'bytes={a}-{b}','Accept-Encoding':'identity'},timeout=(30,160));x.raise_for_status()
        if x.status_code!=206 or len(x.content)!=b-a+1:raise ValueError('Incomplete range')
        return x.content
    head=get(0,28799);h,offset=header(head);nx,ny,nz=[h[f'NAXIS{i}'] for i in [1,2,3]];assert ny==1
    ext=offset+((nx*nz*4+2879)//2880)*2880;tail=get(ext,ext+11519);zh,zo=header(tail)
    z=np.frombuffer(tail[zo:zo+nz*4],'>f4').astype(float)
    take=np.flatnonzero((z>=.6)&(z<=14));lo,hi=int(take[0]),int(take[-1]+1)
    a=offset+lo*nx*4;b=offset+hi*nx*4-1;payload=get(a,b)
    arr=np.frombuffer(payload,'>f4').reshape(hi-lo,nx).astype('f4')
    rec=dict(url=url,variable=var,snapshot=snap,time_s=float(h['ELAPSED']),source_shape=[nz,ny,nx],range=[a,b],sha256=hashlib.sha256(payload).hexdigest(),header_sha256=hashlib.sha256(head).hexdigest(),units=h['BUNIT'],dx_Mm=float(h['CDELT1']),source_header=dict(h))
    return arr,z[lo:hi],rec

def fetch(snap):
    dest=OUT/f'{snap}.npz'
    if dest.exists() and dest.with_suffix('.json').exists():return
    begin=time.monotonic();a,z,ra=field('lgtg',snap);b,zb,rb=field('lgne',snap)
    assert np.array_equal(z,zb) and ra['time_s']==rb['time_s']
    # ne is SI m^-3. AIA response uses ne^2 cm^-6 and LOS cm.
    T=10.**a.astype(float);ne=10.**(b.astype(float)-6)
    j=ne*ne*10.**np.interp(np.log10(T),logT,logR,left=-100,right=-100)
    j[T<1e4]=0
    # Conservation-weighted block average of emissivity itself (never T or ne).
    fz,fx=4,4;nz=(len(z)//fz)*fz;nx=(a.shape[1]//fx)*fx
    edges=np.r_[z[0]-(z[1]-z[0])/2,(z[1:]+z[:-1])/2,z[-1]+(z[-1]-z[-2])/2]
    dz=np.diff(edges);weights=dz[:nz].reshape(-1,fz)
    out=(j[:nz,:nx].reshape(nz//fz,fz,nx//fx,fx)*weights[:,:,None,None]).sum(axis=(1,3))/(fx*weights.sum(axis=1)[:,None])
    zz=(z[:nz].reshape(-1,fz)*weights).sum(axis=1)/weights.sum(axis=1)
    xx=(np.arange(nx//fx)*fx+(fx-1)/2)*ra['dx_Mm']
    total=float((j[:nz,:nx]*dz[:nz,None]).sum()*ra['dx_Mm']);reduced=float((out*weights.sum(axis=1)[:,None]).sum()*ra['dx_Mm']*fx)
    np.savez_compressed(dest,jUV=out.astype('f4'),x_Mm=xx.astype('f4'),z_Mm=zz.astype('f4'),time_s=ra['time_s'],z_edges_Mm=edges[:nz+1:4],physical_emission_integral=total)
    dest.with_suffix('.json').write_text(json.dumps(dict(source_fields=[ra,rb],time_s=ra['time_s'],lines_A=[1548.187,1550.772],response_sha256=hashlib.sha256(rpath.read_bytes()).hexdigest(),shape=list(out.shape),reduction_integral_error=abs(reduced/total-1),scope='Actual 2.5D spicule-producing MHD thermodynamics; equilibrium optically thin C IV doublet emission with radiative cascades; excludes Halpha and non-equilibrium carbon kinetics',elapsed_s=time.monotonic()-begin),indent=2,default=str)+'\n')
    print('Ready',snap,ra['time_s'],out.shape,round(time.monotonic()-begin,1),'s',flush=True)

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--test',action='store_true');args=ap.parse_args()
    selected=[308] if args.test else list(range(300,371))
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:list(pool.map(fetch,selected))
