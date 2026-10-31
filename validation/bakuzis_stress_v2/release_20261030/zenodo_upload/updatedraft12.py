"""updatedraft12.py: refresh the UNPUBLISHED version 1.2.0 draft (deposition 23814198) to catalog v8.
Syncs files to files_to_upload.txt by md5 and PUTs metadata. No new version, never publishes."""
import requests, hashlib, os, json
T = open(os.path.expanduser('~/.config/zenodo/token')).read().strip(); H = {'Authorization': 'Bearer ' + T}
A = 'https://zenodo.org/api/deposit/depositions'; DID = 23814198
d = requests.get(f'{A}/{DID}', headers=H, timeout=30).json()
print('draft', DID, 'state', d['state'], 'submitted', d.get('submitted'), 'files', len(d['files']), 'version', d['metadata'].get('version'))
if d.get('submitted') or d['state'] != 'unsubmitted': raise SystemExit('draft is not an unpublished draft; stopping')
bucket = d['links']['bucket']
srv = {f['filename']: (f['checksum'].replace('md5:', ''), f['id']) for f in d['files']}
local = {}
for ln in open('files_to_upload.txt'):
    ln = ln.split('#', 1)[0].strip()
    if ln: local[os.path.basename(ln)] = ln
removed = put = 0
for b, (md5, fid) in srv.items():
    if b not in local:
        r = requests.delete(f'{A}/{DID}/files/{fid}', headers=H, timeout=60); print('removed', b, r.status_code); removed += 1
for b, p in local.items():
    m = hashlib.md5(open(p, 'rb').read()).hexdigest()
    if srv.get(b, ('', ''))[0] == m: continue
    with open(p, 'rb') as fh: r = requests.put(f'{bucket}/{b}', data=fh, headers=H, timeout=600)
    ok = r.status_code in (200, 201) and r.json().get('checksum') == 'md5:' + m
    print('put', b, r.status_code, ok); put += 1
    if not ok: print(r.text[:200])
m = d['metadata']; m.pop('doi', None); m = {k: v for k, v in m.items() if v not in ('', [], None)}
new = json.load(open('zenodo_metadata.json'))['metadata']; m.update(new)
r = requests.put(f'{A}/{DID}', headers={**H, 'Content-Type': 'application/json'}, data=json.dumps({'metadata': m}), timeout=60)
print('metadata PUT', r.status_code, r.text[:300] if r.status_code >= 400 else 'ok')
d = requests.get(f'{A}/{DID}', headers=H, timeout=30).json(); srv = {f['filename']: f['checksum'].replace('md5:', '') for f in d['files']}
print('removed', removed, 'uploaded', put)
print('files', len(srv), 'manifest', len(local), 'missing', [b for b in local if b not in srv], 'extra', [b for b in srv if b not in local],
      'mismatch', [b for b, p in local.items() if hashlib.md5(open(p, 'rb').read()).hexdigest() != srv.get(b)],
      'state', d['state'], 'submitted', d.get('submitted'), 'version', d['metadata'].get('version'), 'reserved', d['metadata'].get('prereserve_doi'))
