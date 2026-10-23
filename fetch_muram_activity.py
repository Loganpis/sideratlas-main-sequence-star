"""Fetch matching magnetic slices from the CC BY 3.0 MURaM archive.
Only exact HTTP ranges of public numeric fields; no motion is synthesized.
"""
from pathlib import Path
import requests,numpy as np,tarfile,json,hashlib,time
from concurrent.futures import ThreadPoolExecutor
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
ROOT=Path(__file__).resolve().parent
URL='https://stacks.stanford.edu/file/dv883vb9686/tau_slices_every200steps.tar'
OUT=ROOT/'data/active-region/magnetic';OUT.mkdir(parents=True,exist_ok=True)
STRIDE=5243904;PLANE=512*256*4

def fetch(it):
 p=OUT/f'{it}.npz'
 if p.exists():return
 index=(int(it)-270000)//200;offset=index*STRIDE
 s=requests.Session();s.mount('https://',HTTPAdapter(max_retries=Retry(total=4,backoff_factor=1,status_forcelist=[429,500,502,503,504])))
 def rng(a,b):
  r=s.get(URL,headers={'Range':f'bytes={a}-{b}','Accept-Encoding':'identity'},timeout=(20,100));r.raise_for_status()
  if r.status_code!=206 or len(r.content)!=b-a+1:raise ValueError('Wrong range response')
  return r.content
 head=rng(offset,offset+1023);h=tarfile.TarInfo.frombuf(head[:512],'utf-8','surrogateescape');assert int(h.name.split('.')[-1])==it
 info=np.frombuffer(head[512:528],'<f4');assert list(info[:3])==[10,512,256]
 start=offset+512+16+5*PLANE;end=start+3*PLANE-1;data=rng(start,end)
 a=np.frombuffer(data,'<f4').reshape(3,256,512)*np.sqrt(4*np.pi)
 np.savez_compressed(p,bx_G=a[0],by_G=a[2],bz_G=a[1],time_s=info[3],iteration=it)
 p.with_suffix('.json').write_text(json.dumps(dict(url=URL,member=h.name,byte_range=[start,end],sha256=hashlib.sha256(data).hexdigest(),license='CC BY 3.0',time_s=float(info[3]),field_order='bx, vertical bz, horizontal by in archived planes 5,6,7; multiplied by sqrt(4 pi) to gauss'),indent=2)+'\n')
 print('Fetched magnetic',it,float(info[3]),flush=True)

if __name__=='__main__':
 d=np.load(ROOT/'data/active-region/intensity-sequence.npz');take=d['time_s']<=float(d['time_s'][0])+1320
 with ThreadPoolExecutor(max_workers=3) as ex:list(ex.map(fetch,d['iteration'][take]))
