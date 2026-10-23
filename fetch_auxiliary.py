"""Fetch the institutional solar reference spectrum, retaining its checksum."""
from pathlib import Path
import hashlib,requests,json
ROOT=Path(__file__).resolve().parent
p=ROOT/'data/tsis1-hsrs-v2-1nm.nc'
record=json.loads((ROOT/'data/source-manifest.json').read_text())['reference_spectrum']
if not p.exists():
 r=requests.get(record['source'],timeout=(20,120));r.raise_for_status();p.write_bytes(r.content)
actual=hashlib.sha256(p.read_bytes()).hexdigest()
if actual!=record['sha256']:raise ValueError('Reference spectrum changed; review provenance before using it')
print('Verified official TSIS-1 HSRS v2:',actual)
