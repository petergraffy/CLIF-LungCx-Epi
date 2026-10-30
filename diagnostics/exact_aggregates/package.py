"""Package only regenerated exact aggregates, never reverse engineer old bands."""
from pathlib import Path
import csv,json,hashlib,shutil,os
ROOT=Path(__file__).resolve().parents[2]
OUT=Path(os.environ.get('CLIF_EXACT_OUTPUT',ROOT/'output/run_exact_aggregates_20261030')).resolve()
SITES=['UCMC','NU','RUSH']
forbidden={'patient_id','hospitalization_id','t0','event_time','recorded_dttm','admission_dttm','discharge_dttm','sequence','marker','threshold'}
h=lambda p:hashlib.md5(p.read_bytes()).hexdigest()
consensus_ids=set()
for site in SITES:
 dest=OUT/'site_packets'/site/'share'
 if dest.exists():raise ValueError('Choose a fresh exact packet')
 for source,target,manifest in [('descriptives/share/descriptives','federated/descriptives','manifest.json'),('survival/share','federated/survival_tdauc','manifest.json'),('competing_events/share','federated/competing_events','curve_manifest.json'),('counts','exact_counts','manifest.json')]:
  src=OUT/'sites'/site/source;m=json.loads((src/manifest).read_text());assert m['site']==site and m['status']=='passed'
  cid=m.get('consensus_id',m.get('primary_consensus_id'));consensus_ids.add(cid)
  assert len(consensus_ids)==1 and cid
  for k,v in m['files'].items():assert h(src/k)==v
  for p in src.rglob('*.csv'):
   with p.open() as f:
    reader=csv.DictReader(f);assert not forbidden.intersection(reader.fieldnames)
    for row in reader:
     for base in ['n','events','patients','observed','missing','at_risk','classified','case','control']:
      lo,hi=row.get(base+'_lower',''),row.get(base+'_upper','')
      if lo and hi:assert float(lo)==float(hi),(site,p.name,base,lo,hi)
  shutil.copytree(src,dest/target)
 files=[p for p in dest.rglob('*') if p.is_file()]
 (dest/'site_return_manifest.json').write_text(json.dumps({'schema':'exact-aggregate-return-v1','status':'passed','site':site,'aggregate_policy':'exact aggregate counts; no disclosure suppression','primary_consensus_id':next(iter(consensus_ids)),'patient_level_identifiers_returned':False,'clinical_model_refits':False,'historical_bands_not_inferred':True,'files':{str(p.relative_to(dest)):h(p) for p in sorted(files)}},indent=2))
# Exact evaluation N for the absolute-AUROC display; performance estimates stay frozen.
evals={}
for site in SITES:
 with (OUT/'sites'/site/'counts/models/sofa_evaluation_denominators.csv').open() as f:
  for row in csv.DictReader(f):
   cohort,endpoint=row['cohort_endpoint'].split(' ',1);evals[(site,cohort,endpoint)]=row
src=Path(os.environ.get('CLIF_REVISION_OUTPUT',ROOT/'output/run_revision_gaps_20261026')).resolve()/'reports/tables/score_model_metadata.csv'
with src.open() as f:rows=list(csv.DictReader(f))
for row in rows:
 x=evals[(row['site'],row['cohort'],row['endpoint'])]
 for prefix,key in [('n','n'),('events','events')]:
  row[prefix]=x[key];row[prefix+'_lower']=row[prefix+'_upper']=x[key]
with (OUT/'score_model_metadata_exact.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print('Three exact aggregate packets built; no historical interval midpoints used')
