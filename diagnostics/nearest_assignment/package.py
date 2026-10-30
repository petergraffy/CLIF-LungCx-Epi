"""Complete refreshed site packets, using exact cached source diagnostics."""
from pathlib import Path
import csv,json,hashlib,shutil
ROOT=Path(__file__).resolve().parents[2];RUN=ROOT/'output/run_nearest_assignment_20261030'
h=lambda p:hashlib.md5(p.read_bytes()).hexdigest()
consensus=json.loads((RUN/'primary/coordinator/full/frozen_manifest.json').read_text())['consensus_id']
for site in ['UCMC','NU','RUSH']:
 source=RUN/'exact/site_packets'/site/'share';dest=RUN/'site_packets'/site/'share'
 if dest.exists():raise ValueError('Choose fresh complete packet')
 shutil.copytree(source,dest)
 def copy(src,relative):
  p=dest/relative;p.parent.mkdir(parents=True,exist_ok=True)
  if src.is_dir():shutil.copytree(src,p,dirs_exist_ok=True)
  else:shutil.copy2(src,p)
 copy(RUN/'primary/sites'/site/'export','federated/primary')
 copy(RUN/'revision/sites'/site/'export/scores','federated/scores')
 review=RUN/'revision/sites'/site/'export'
 for p in review.iterdir():
  if p.is_file():copy(p,Path('federated/review')/p.name)
 for v in ['h24','h48','onset','support_only','low_flow','carry0','carry12','raw_support']:
  copy(RUN/'revision/sensitivities'/v/'sites'/site/'export',Path('federated/sensitivities')/v)
 copy(RUN/'revision/sensitivities/absorbing_v2/sites'/site/'export','federated/absorbing')
 copy(RUN/'apache_adapted/sites'/site/'export','federated/apache_adapted')
 copy(RUN/'exact/sites'/site/'counts/supplement','supplement')
 for name in ['model_PH_and_denominators.csv','model_PH_terms.csv','exploratory_piecewise_associations.csv','exploratory_piecewise_metadata.csv']:
  copy(RUN/'exact/sites'/site/'model_diagnostics/share'/name,Path('federated')/name)
 copy(RUN/'exact/sites'/site/'competing_events/share/cumulative_incidence_daily.csv','federated/primary/outcomes/full/cumulative_incidence_daily.csv')
 # Inventory covers the actual refreshed exports, not the parent packet.
 inventory=dest/'supplement/return_file_inventory.csv'
 with inventory.open('w') as f:
  w=csv.DictWriter(f,fieldnames=['site','file','status']);w.writeheader()
  for p in sorted(dest.rglob('*')):
   if p.is_file() and p.name!='site_return_manifest.json':w.writerow({'site':site,'file':str(p.relative_to(dest)),'status':'captured'})
 files=[p for p in dest.rglob('*') if p.is_file() and p.name!='site_return_manifest.json']
 for p in files:
  if p.suffix!='.csv':continue
  with p.open() as f:
   r=csv.DictReader(f);assert not {'patient_id','hospitalization_id','t0','event_time','recorded_dttm','admission_dttm','discharge_dttm','sequence'}.intersection(r.fieldnames)
   for row in r:
    for base in ['n','events','observed','missing','patients','at_risk','available']:
     lo,hi=row.get(base+'_lower',''),row.get(base+'_upper','')
     if lo and hi:assert float(lo)==float(hi),(p,base,lo,hi)
    assert row.get('cell_status') not in ['suppressed_partition','suppressed_cluster']
 m={'status':'complete_local_capture_pending_site_metadata','site':site,'primary_consensus_id':consensus,'assignment_rule':'nearest template for every support-sufficient stay; distance/margin flags are diagnostic only','aggregate_policy':'exact counts; no disclosure suppression',
    'apache_adapted_extension_complete':True,'publication_descriptives_complete':True,'daily_competing_event_curves_complete':True,'survival_tdauc_complete':True,
    'frozen_templates_unchanged':True,'classification_dependent_models_refitted':True,'source_diagnostics_reused_from_frozen_same_cohort':True,'author_confirmation_required':True,
    'files':{str(p.relative_to(dest)):h(p) for p in sorted(files)}}
 (dest/'site_return_manifest.json').write_text(json.dumps(m,indent=2))
print('Complete exact site packets refreshed with nearest-assignment models and diagnostics')
