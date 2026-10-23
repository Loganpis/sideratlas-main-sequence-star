"""Reproduce this final visualization from public byte-verified numerical data.

No new hydrodynamic solver is run. Requires substantial source downloads and
space for the original optical atmosphere subsets. Cached inputs are reused.
"""
from pathlib import Path
import subprocess,sys,os,json,hashlib
from concurrent.futures import ThreadPoolExecutor
ROOT=Path(__file__).resolve().parent;sys.path[:0]=[str(ROOT/'src'),str(ROOT/'golden')]
os.environ.setdefault('NUMBA_NUM_THREADS','3');os.environ.setdefault('OPENBLAS_NUM_THREADS','1')

def run(name):subprocess.run([sys.executable,str(ROOT/name)],check=True,cwd=ROOT)

def main():
    for d in ['data/raw','data/continuum','data/provenance','golden/data','golden/output','golden/validation']:(ROOT/d).mkdir(parents=True,exist_ok=True)
    run('fetch_auxiliary.py')
    from archive import fetch
    snaps=sorted(list(range(385,542,4))+[387,544])
    with ThreadPoolExecutor(max_workers=3) as ex:
        jobs=[ex.submit(fetch,v,s) for s in snaps for v in ['lgtg','lgne','lgn1']]
        for j in jobs:j.result()
    # Assert source identity against the retained manifest rather than silently
    # accepting changed institutional inputs. Download transport is not physics.
    manifest=json.loads((ROOT/'data/source-manifest.json').read_text())
    for record in manifest['records']:
        current=json.loads((ROOT/f"data/raw/{record['variable']}-{record['snapshot']}.json").read_text())
        if current['download_sha256']!=record['download_sha256']:raise ValueError('Optical source payload changed; review provenance')
    from continuum import synthesize
    synthesize(385,17,2,False)
    for s in snaps:synthesize(s,17,2,True)
    run('prepare_muram_intensity.py');run('fetch_muram_activity.py')
    run('golden/fetch_spicules.py')
    expected=json.loads((ROOT/'golden/data/source-manifest.json').read_text())
    import numpy as np
    for record in expected['source_states']:
        path=ROOT/'golden/data/spicules-uv'/record['file']
        current=json.loads(path.with_suffix('.json').read_text())
        if [f['sha256'] for f in current['source_fields']]!=[f['sha256'] for f in record['source_fields']]:raise ValueError('Spicule source payload changed; review provenance')
        d=np.load(path)
        actual=hashlib.sha256(d['jUV'].tobytes()).hexdigest()
        if actual!=record['processed_emissivity_sha256']:
            values=d['jUV'].astype(float);moments=[values.sum(),np.sum(values*values),values.max()]
            if not np.allclose(moments,record['processed_moments'],rtol=5e-6,atol=0):raise ValueError('Reproduced UV emissivity changed beyond numerical tolerance; review provenance')
            print('Platform rounding differs within 5e-6 moment tolerance:',record['file'])
    print('Public model inputs and calculated radiation are ready; render with golden_orbit/render.py --config config-final.json')

if __name__=='__main__':main()
