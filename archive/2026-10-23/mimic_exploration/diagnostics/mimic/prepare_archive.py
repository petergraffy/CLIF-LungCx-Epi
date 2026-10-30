"""Extract supplied MIMIC CLIF tables into an ignored site-local workspace."""
from pathlib import Path
import argparse,json,zipfile,shutil
p=argparse.ArgumentParser();p.add_argument('archive',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
a.output.mkdir(parents=True,exist_ok=True)
required={'patient','hospitalization','adt','hospital_diagnosis','respiratory_support','vitals','labs','medication_admin_continuous','patient_assessments'}
allowed=required|{'code_status','patient_assessments_raw_gcs'}
manifest=[]
with zipfile.ZipFile(a.archive) as z:
 for f in z.infolist():
  name=Path(f.filename).name
  key=name.removeprefix('clif_').removesuffix('.parquet')
  if name.endswith('.parquet') and key in allowed:
   target=a.output/name
   if target.exists():
    if target.stat().st_size!=f.file_size:raise RuntimeError('Partial or mismatched existing table')
   else:
    temp=target.with_suffix('.partial')
    with z.open(f) as src,temp.open('wb') as dst:shutil.copyfileobj(src,dst,1024*1024)
    temp.rename(target)
   manifest.append({'table':key,'archive_member':f.filename,'bytes':f.file_size,'crc32':f.CRC})
  if name.lower()=='readme.md':
   (a.output/'source_README.txt').write_bytes(z.read(f))
found={m['table'] for m in manifest}
if not required<=found:raise RuntimeError(f'Missing required tables: {sorted(required-found)}')
(a.output/'extraction_manifest.json').write_text(json.dumps({'archive':str(a.archive),'source_tables_unchanged':True,'tables':manifest},indent=2))
print('Extracted',len(manifest),'CLIF tables; patient data remain in ignored local workspace.',flush=True)
