"""newversion12.py: stage version 1.2.0 of the fvs-modern catalog record as an UNPUBLISHED draft.
Creates a new version from 23788727 (v1.1.0), syncs files to files_to_upload.txt by md5, PUTs metadata.
Never calls the publish action. Token from ~/.config/zenodo/token, never printed."""
import requests, hashlib, os, json
T = open(os.path.expanduser('~/.config/zenodo/token')).read().strip(); H = {'Authorization': 'Bearer ' + T}
A = 'https://zenodo.org/api/deposit/depositions'
r = requests.post(f'{A}/23788727/actions/newversion', headers=H, timeout=60)
print('newversion', r.status_code, r.text[:300] if r.status_code >= 400 else '')
if r.status_code >= 400:
    raise SystemExit('newversion failed; inspect before retrying')
draft = r.json()['links']['latest_draft']
d = requests.get(draft, headers=H, timeout=30).json(); did = d['id']; bucket = d['links']['bucket']
print('draft id', did, 'state', d['state'], 'submitted', d.get('submitted'), 'inherited files', len(d['files']))
srv = {f['filename']: (f['checksum'].replace('md5:', ''), f['id']) for f in d['files']}
local = {}
for ln in open('files_to_upload.txt'):
    ln = ln.split('#', 1)[0].strip()
    if ln: local[os.path.basename(ln)] = ln
for b, (md5, fid) in srv.items():
    if b not in local:
        r = requests.delete(f'{A}/{did}/files/{fid}', headers=H, timeout=60); print('removed', b, r.status_code)
for b, p in local.items():
    m = hashlib.md5(open(p, 'rb').read()).hexdigest()
    if srv.get(b, ('', ''))[0] == m: continue
    with open(p, 'rb') as fh: r = requests.put(f'{bucket}/{b}', data=fh, headers=H, timeout=600)
    ok = r.status_code in (200, 201) and r.json().get('checksum') == 'md5:' + m
    print('put', b, r.status_code, ok)
    if not ok: print(r.text[:200])
m = d['metadata']; m.pop('prereserve_doi', None); m.pop('doi', None); m = {k: v for k, v in m.items() if v not in ('', [], None)}
new = json.load(open('zenodo_metadata.json'))['metadata']; m.update(new)
r = requests.put(f'{A}/{did}', headers={**H, 'Content-Type': 'application/json'}, data=json.dumps({'metadata': m}), timeout=60)
print('metadata PUT', r.status_code, r.text[:300] if r.status_code >= 400 else 'ok')
d = requests.get(f'{A}/{did}', headers=H, timeout=30).json(); srv = {f['filename']: f['checksum'].replace('md5:', '') for f in d['files']}
print('files', len(srv), 'manifest', len(local), 'missing', [b for b in local if b not in srv], 'extra', [b for b in srv if b not in local],
      'mismatch', [b for b, p in local.items() if hashlib.md5(open(p, 'rb').read()).hexdigest() != srv.get(b)],
      'state', d['state'], 'submitted', d.get('submitted'), 'version', d['metadata'].get('version'), 'reserved', d['metadata'].get('prereserve_doi'))
print('edit URL', d['links'].get('html'))
