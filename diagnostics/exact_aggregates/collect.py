"""Collect exact aggregate packets; never open site clinical caches."""
from pathlib import Path
import csv,json,hashlib,sys,subprocess
args=sys.argv[1:]
if len(args)<2:raise SystemExit('Usage: collect.py NEW_OUTPUT SITE_SHARE [SITE_SHARE ...]')
out=Path(args[0]);shares=[Path(p) for p in args[1:]]
if out.exists():raise ValueError('Choose a fresh collection')
sites=[];consensus=[];tables={}
for share in shares:
 m=json.loads((share/'site_return_manifest.json').read_text())
 assert m['schema']=='exact-aggregate-return-v1' and m['status']=='passed'
 site=m['site'];assert site not in sites;sites.append(site);consensus.append(m['primary_consensus_id'])
 actual={str(p.relative_to(share)) for p in share.rglob('*') if p.is_file() and p.name!='site_return_manifest.json'}
 assert actual==set(m['files'])
 for rel,expected in m['files'].items():
  p=Path(rel);assert not p.is_absolute() and '..' not in p.parts
  file=share/p;assert hashlib.md5(file.read_bytes()).hexdigest()==expected
  if p.suffix!='.csv':continue
  with file.open() as f:
   reader=csv.DictReader(f);cols=reader.fieldnames
   assert not {'patient_id','hospitalization_id','t0','event_time','recorded_dttm','sequence','marker','threshold'}.intersection(cols)
   rows=list(reader)
  for r in rows:
   assert r['site']==site
   for base in ['n','events','patients','observed','missing','at_risk','classified','case','control']:
    lo,hi=r.get(base+'_lower',''),r.get(base+'_upper','')
    if lo and hi:assert float(lo)==float(hi)
   assert r.get('cell_status') not in {'suppressed_partition','suppressed_cluster'}
  tables.setdefault(rel,[]).extend(rows)
assert len(set(consensus))==1
out.mkdir(parents=True)
for rel,rows in tables.items():
 if not rows:continue
 fields=list(dict.fromkeys(k for r in rows for k in r));p=out/rel;p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
subprocess.run(['Rscript','diagnostics/descriptives/collect.R',str(out/'publication_tables'),*[str(p) for p in shares]],check=True)
(out/'collection_manifest.json').write_text(json.dumps({'status':'passed','sites':sites,'primary_consensus_id':consensus[0],'aggregate_policy':'exact counts; no disclosure suppression','patient_level_inputs_read':False,'counts_restored_from_local_caches_not_old_bands':True,'site_distinct_patient_counts_not_summed':True,'clinical_model_refits':False},indent=2))
print('Exact aggregate collection passed')
