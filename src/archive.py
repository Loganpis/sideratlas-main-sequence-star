"""Small, explicit byte-range subsets of the published Bifrost FITS cubes."""
from pathlib import Path
import json,hashlib,time
import numpy as np
import requests
from astropy.io import fits
ROOT=Path(__file__).resolve().parents[1]
BASE='https://sdc.uio.no/vol/simulations/en024048_hion/atmos/'

def byte_range(session,url,start,stop):
 r=session.get(url,headers={'Range':f'bytes={start}-{stop}','Accept-Encoding':'identity'},timeout=(20,150));r.raise_for_status()
 if r.status_code!=206:raise RuntimeError('Server did not honor uncompressed byte range')
 return r.content

def header(data):
 end=next(i+80 for i in range(0,len(data),80) if data[i:i+8]==b'END     ')
 return fits.Header.fromstring(data[:end].decode('ascii'),sep=''),(end+2879)//2880*2880

def fetch(variable,snapshot,zmin=-.6,zmax=.85):
 out=ROOT/'data/raw'/f'{variable}-{snapshot}.npz';meta=out.with_suffix('.json')
 out.parent.mkdir(parents=True,exist_ok=True)
 if out.exists() and meta.exists():return str(out)
 url=BASE+f'BIFROST_en024048_hion_{variable}_{snapshot}.fits'
 s=requests.Session();s.headers['Accept-Encoding']='identity';begin=time.monotonic()
 head=byte_range(s,url,0,28799);h,offset=header(head);nx,ny,nz=[h[f'NAXIS{i}'] for i in [1,2,3]]
 size=nx*ny*nz*4;ext=offset+(size+2879)//2880*2880;end=byte_range(s,url,ext,ext+5759);eh,eo=header(end)
 z=np.frombuffer(end[eo:eo+nz*4],dtype='>f4').astype('f8');take=np.flatnonzero((z>=zmin)&(z<=zmax));lo,hi=int(take[0]),int(take[-1]+1)
 start=offset+lo*nx*ny*4;stop=offset+hi*nx*ny*4-1
 data=byte_range(s,url,start,stop);assert len(data)==(hi-lo)*nx*ny*4
 a=np.frombuffer(data,dtype='>f4').reshape(hi-lo,ny,nx).astype('f4')
 np.savez_compressed(out,values=a,z_Mm=z[lo:hi],dx_Mm=h['CDELT1'],time_s=h['ELAPSED'],source_teff_K=h['TEFF'])
 record=dict(source=url,variable=variable,snapshot=snapshot,shape=list(a.shape),source_shape=[nz,ny,nx],source_header=dict(h),byte_range=[start,stop],download_sha256=hashlib.sha256(data).hexdigest(),file_sha256=hashlib.sha256(out.read_bytes()).hexdigest(),elapsed_s=time.monotonic()-begin,scope='Unmodified numeric log/linear field subset; no interpolated dynamics')
 meta.write_text(json.dumps(record,indent=2,default=str)+'\n');print('Fetched',variable,snapshot,a.shape,round(record['elapsed_s'],1),'s',flush=True)
 return str(out)
