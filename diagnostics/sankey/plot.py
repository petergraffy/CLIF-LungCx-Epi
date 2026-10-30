"""Initial ICU-hour support -> cluster Sankey from existing aggregate descriptors."""
import argparse,csv,hashlib,json,os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/clif_manuscript_mpl')
SUPPORT=['ROOM_AIR','LOW_O2','HFNC','CPAP','NIPPV','IMV','OTHER','UNKNOWN']
LABELS=dict(zip(SUPPORT,['Room air','Low-flow oxygen','HFNC','CPAP','NIPPV','IMV','Other','Unknown']))
COLORS=dict(zip(SUPPORT,['#B8D4E8','#347DB1','#4FA343','#A45BB5','#713B94','#E65D0A','#D4B15C','#999999']))
PHENOTYPES={'P1':'Room-air predominant','P2':'Low-flow oxygen','P3':'HFNC predominant','P4':'IMV predominant','I':'Insufficient information'}
PC={'P1':'#B8D4E8','P2':'#347DB1','P3':'#4FA343','P4':'#E65D0A','I':'#999999'}

def read(p):
 with p.open() as f:return list(csv.DictReader(f))

def integer(x):
 v=float(x)
 if v<0 or not v.is_integer():raise ValueError('Expected an exact nonnegative hospitalization count')
 return int(v)

def extract_counts(categorical,denominators,availability,grouping):
 """Require complete, exact partitions; never infer unavailable cells as zero."""
 cohort='primary_phenotypes' if grouping=='phenotype' else 'provisional_local_clusters'
 den=[r for r in denominators if r['cohort']==cohort and r['group']!='Overall']
 if not den or len({r['group'] for r in den})!=len(den):raise ValueError('Missing or duplicate cluster denominators')
 expected={r['group']:integer(r['n']) for r in den}
 rows=[r for r in categorical if r['cohort']==cohort and r['group']!='Overall' and r['variable']=='first_hour_support']
 available={r['group']:r['status'] for r in availability if r['cohort']==cohort and r['variable']=='first_hour_support'}
 flows={};seen=set()
 for r in rows:
  group=r['group'];support=r['category']
  if support in ['MISSING','',None]:support='UNKNOWN'
  if support not in SUPPORT:raise ValueError('Unexpected support category: '+str(support))
  if group not in expected or r['cell_status']!='observed' or str(r['partition_complete']).lower()!='true':raise ValueError('Incomplete/suppressed first-hour support partition')
  if integer(r['denominator_n'])!=expected[group]:raise ValueError('Conflicting group denominators')
  key=(support,group)
  if (r['category'],group) in seen:raise ValueError('Duplicate source category/group')
  seen.add((r['category'],group));flows[key]=flows.get(key,0)+integer(r['n'])
 for group,n in expected.items():
  if available.get(group) not in ['calculated','structural_zero_group']:raise ValueError('First-hour support unavailable for '+group)
  if sum(v for (s,g),v in flows.items() if g==group)!=n:raise ValueError('Flows do not reconcile to cluster denominator')
 if sum(expected.values())==0:raise ValueError('No hospitalizations to plot')
 return flows,expected

def load_share(path):
 candidates=[path/'federated/descriptives',path/'share/descriptives',path/'descriptives',path]
 base=next((p for p in candidates if (p/'manifest.json').exists() and (p/'categorical.csv').exists()),None)
 if base is None:raise ValueError('Cannot locate aggregate descriptors in '+str(path))
 m=json.loads((base/'manifest.json').read_text())
 if m['status']!='passed' or m['grouping'] not in ['phenotype','local_cluster']:raise ValueError('Unsupported descriptor return')
 for f,h in m['files'].items():
  p=Path(f)
  if p.is_absolute() or '..' in p.parts:raise ValueError('Unsafe descriptor path')
  if hashlib.md5((base/p).read_bytes()).hexdigest()!=h:raise ValueError('Descriptor hash mismatch')
 tables=[read(base/f) for f in ['categorical.csv','denominators.csv','availability.csv']]
 for table in tables:
  if any(r.get('site')!=m['site'] for r in table):raise ValueError('Site mismatch')
 flows,den=extract_counts(*tables,m['grouping'])
 return m,flows,den,base

def draw(flows,groups,title,subtitle,path,shared):
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from matplotlib.path import Path as MPath
 from matplotlib.patches import PathPatch,Rectangle
 import numpy as np
 sources=[s for s in SUPPORT if sum(v for (a,b),v in flows.items() if a==s)>0]
 targets=[g for g in groups if groups[g]>0]
 if shared:targets=sorted(targets,key=lambda g:(g=='I',int(g[1:]) if g.startswith('P') else 0))
 else:targets=sorted(targets,key=lambda g:int(g))
 N=sum(groups.values());gap=.022;scale=(.82-gap*(max(len(sources),len(targets))-1))/N
 def nodes(keys,totals):
  used=scale*sum(totals.values())+gap*(len(keys)-1);top=.5+used/2;out={}
  for k in keys:height=totals[k]*scale;out[k]=[top-height,top];top-=height+gap
  return out
 source_total={s:sum(v for (a,b),v in flows.items() if a==s) for s in sources}
 left=nodes(sources,source_total);right=nodes(targets,{g:groups[g] for g in targets})
 fig,ax=plt.subplots(figsize=(13,8));fig.subplots_adjust(left=.025,right=.975,top=.82,bottom=.14)
 ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
 lc={s:left[s][1] for s in sources};rc={g:right[g][1] for g in targets}
 for s in sources:
  for g in targets:
   count=flows.get((s,g),0)
   if not count:continue
   h=count*scale;y0=lc[s];y1=rc[g];lc[s]-=h;rc[g]-=h
   vertices=[(.27,y0),(.43,y0),(.57,y1),(.73,y1),(.73,y1-h),(.57,y1-h),(.43,y0-h),(.27,y0-h),(.27,y0)]
   codes=[MPath.MOVETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4,MPath.LINETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4,MPath.CLOSEPOLY]
   ax.add_patch(PathPatch(MPath(vertices,codes),facecolor=COLORS[s],edgecolor='none',alpha=.48))
 cmap=plt.get_cmap('tab10')
 for s,(lo,hi) in left.items():
  ax.add_patch(Rectangle((.25,lo),.02,hi-lo,facecolor=COLORS[s],edgecolor='white',lw=.4))
  ax.text(.239,(lo+hi)/2,f'{LABELS[s]}  {source_total[s]:,}',ha='right',va='center',fontsize=11)
 for i,g in enumerate(targets):
  lo,hi=right[g];color=PC.get(g,cmap(i%10))
  ax.add_patch(Rectangle((.73,lo),.02,hi-lo,facecolor=color,edgecolor='white',lw=.4))
  label=(g+': '+PHENOTYPES.get(g,g)) if shared else 'Local cluster '+g
  ax.text(.761,(lo+hi)/2,f'{label}\nN={groups[g]:,} ({100*groups[g]/N:.1f}%)',ha='left',va='center',fontsize=10.5)
 ax.text(.25,.98,'Support in first ICU hour',ha='center',va='bottom',fontsize=12,fontweight='bold')
 ax.text(.75,.98,'Shared 72-hour phenotype' if shared else 'Site-local 72-hour cluster',ha='center',va='bottom',fontsize=12,fontweight='bold')
 fig.suptitle(title,x=.035,ha='left',fontsize=21,y=.97)
 fig.text(.035,.90,subtitle,fontsize=12)
 fig.text(.035,.04,'Ribbon width is the exact hospitalization count; color indicates initial support. Unknown support remains explicit.\n'
          'Initial support is the analysis-panel state at hour 0, after the configured support processing; it is not the first later documented device.\n'
          'Only the primary 72-hour landmark cohort is shown. This is a classification flow, not a count of support transitions or post-landmark outcomes.',fontsize=9,color='#444444')
 for ext in ['png','pdf']:fig.savefig(path.with_suffix('.'+ext),dpi=180,facecolor='white')
 plt.close(fig)

def main():
 parser=argparse.ArgumentParser();parser.add_argument('output',type=Path);parser.add_argument('shares',type=Path,nargs='+');args=parser.parse_args()
 loaded=[load_share(p) for p in args.shares];metas=[x[0] for x in loaded]
 if len({m['site'] for m in metas})!=len(metas):raise ValueError('Duplicate sites')
 if len({m['grouping'] for m in metas})!=1:raise ValueError('Cannot mix local and shared cluster returns')
 shared=metas[0]['grouping']=='phenotype'
 if shared and (len({m['consensus_id'] for m in metas})!=1 or not metas[0]['consensus_id']):raise ValueError('Shared consensus mismatch')
 args.output.mkdir(parents=True,exist_ok=False)
 combined={};den={};export=[]
 for m,flows,groups,base in loaded:
  for (s,g),n in flows.items():
   export.append(dict(site=m['site'],initial_support=s,cluster=g,n=n,grouping=m['grouping'],consensus_id=m.get('consensus_id'),cohort='72h landmark'))
   combined[(s,g)]=combined.get((s,g),0)+n
  for g,n in groups.items():den[g]=den.get(g,0)+n
  draw(flows,groups,'Initial respiratory support to '+('shared phenotype' if shared else 'local cluster'),f'{m["site"]} | {sum(groups.values()):,} landmark hospitalizations',args.output/f'initial_support_to_cluster_{m["site"]}',shared)
 if shared:draw(combined,den,'Initial respiratory support to shared phenotype',f'Pooled {" + ".join(m["site"] for m in metas)} | {sum(den.values()):,} landmark hospitalizations',args.output/'pooled_initial_support_to_phenotype',True)
 with (args.output/'initial_support_to_cluster_counts.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(export[0]));w.writeheader();w.writerows(export)
 if shared:
  rows=[dict(initial_support=s,phenotype=g,n=n) for (s,g),n in combined.items()]
  with (args.output/'pooled_initial_support_to_phenotype_counts.csv').open('w') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 manifest=dict(status='passed',sites=[m['site'] for m in metas],grouping=metas[0]['grouping'],consensus_id=metas[0].get('consensus_id') if shared else None,n=sum(sum(x[2].values()) for x in loaded),pooled_figure_created=shared,patient_level_inputs_read=False,source='categorical.csv: first_hour_support by primary phenotype or provisional local cluster; exclude Overall to avoid double counting',denominators_reconciled=True,initial_definition='analysis-panel h=0 after configured support processing; later documentation does not replace an unknown h=0 state',early_terminal_events_included=False,source_manifest_md5={str(base/'manifest.json'):hashlib.md5((base/'manifest.json').read_bytes()).hexdigest() for m,f,d,base in loaded})
 (args.output/'manifest.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest,indent=2))
if __name__=='__main__':main()
