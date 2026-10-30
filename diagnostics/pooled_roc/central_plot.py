"""Pooled held-out ROC curves from site aggregate score-bin counts only."""
import os,csv,json,hashlib,argparse
os.environ.setdefault('MPLCONFIGDIR','/tmp/clif_manuscript_mpl')
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('output');p.add_argument('shares',nargs='+');a=p.parse_args()
out=Path(a.output);out.mkdir(parents=True,exist_ok=False)
models=['baseline','baseline_sofa','baseline_adapted','phenotype']
labels=['Clinical covariates','Clinical covariates + modified SOFA','Clinical covariates + adapted APACHE','Clinical covariates + phenotype']
colors=['#646464','#32916B','#9D5A98','#357BB5'];B=200
meta=[];data={};site_audit=[];den=[]
for share in map(Path,a.shares):
 m=json.loads((share/'manifest.json').read_text());meta.append(m)
 assert m['status']=='passed' and m['bootstrap']==B and m['bins']==100000 and not m['held_out_refitting']
 assert m['site'] not in m['training_sites'] and set(m['models'])==set(models)
 for name,h in m['files'].items():assert hashlib.md5((share/name).read_bytes()).hexdigest()==h
 with (share/'score_histograms.csv').open() as f:
  for r in csv.DictReader(f):
   assert r['site']==m['site']
   key=(r['endpoint'],r['model'],int(r['day']))
   q=data.setdefault(key,{})
   counts=np.array([float(r['cases']),float(r['controls'])]+[float(r[f'{v}_b{b:03d}']) for b in range(1,B+1) for v in ['cases','controls']]).reshape(B+1,2)
   assert np.all(counts>=0) and np.all(counts==np.floor(counts))
   ibin=int(r['bin']);assert 0<=ibin<100000
   q[ibin]=q.get(ibin,np.zeros_like(counts))+counts
 with (share/'local_validation.csv').open() as f:site_audit+=list(csv.DictReader(f))
 with (share/'cohort_counts.csv').open() as f:den+=list(csv.DictReader(f))
assert len({m['site'] for m in meta})==len(meta) and len({m['consensus_id'] for m in meta})==1
assert len(data)==24
for r in site_audit:
 if r['absolute_difference']:assert float(r['absolute_difference'])<=float(r['auc_quantization_bound'])+1e-12

def curve(c,n):
 if not c.sum() or not n.sum():return None,None,np.nan
 tp=np.r_[0,np.cumsum(c[::-1])/c.sum()];fp=np.r_[0,np.cumsum(n[::-1])/n.sum()]
 auc=(c*(np.cumsum(n)-n+.5*n)).sum()/(c.sum()*n.sum())
 assert abs((np.diff(fp)*(tp[1:]+tp[:-1])/2).sum()-auc)<1e-10
 return fp,tp,auc

def band_value(x,y,grid):
 u=np.unique(x);top=np.array([y[x==v].max() for v in u])
 return np.interp(grid,u,top)
metrics=[];rocs=[];cache={};fprgrid=np.linspace(0,1,101)
for key,bins in data.items():
 endpoint,model,day=key;z=np.stack([bins[k] for k in sorted(bins)])
 fp,tp,auc=curve(z[:,0,0],z[:,0,1]);assert fp is not None
 bs=[];bt=[]
 for b in range(1,B+1):
  x,y,v=curve(z[:,b,0],z[:,b,1]);bs.append(v)
  if x is not None:bt.append(band_value(x,y,fprgrid))
 bs=np.array(bs);bt=np.array(bt)
 ci=np.nanquantile(bs,[.025,.975]);band=np.quantile(bt,[.025,.975],axis=0)
 cases=int(z[:,0,0].sum());controls=int(z[:,0,1].sum())
 bound=(z[:,0,0]*z[:,0,1]).sum()/(2*cases*controls)
 metrics.append(dict(endpoint=endpoint,model=model,day=day,n=cases+controls,cases=cases,controls=controls,auc=auc,conf_low=ci[0],conf_high=ci[1],bootstrap_replicates=int(np.isfinite(bs).sum()),quantization_bound=bound,unrounded_auc_lower=max(0,auc-bound),unrounded_auc_upper=min(1,auc+bound)))
 for x,y in zip(fp,tp):rocs.append(dict(endpoint=endpoint,model=model,day=day,fpr=x,sensitivity=y))
 cache[key]=(fp,tp,auc,ci,band,bs)
for endpoint in ['death','death_or_hospice']:
 for day in [7,14,28]:
  rows=[r for r in metrics if r['endpoint']==endpoint and r['day']==day]
  assert len({(r['n'],r['cases'],r['controls']) for r in rows})==1
pairs=[]
for endpoint in ['death','death_or_hospice']:
 for day in [7,14,28]:
  for comparator in models[:3]:
   f=cache[(endpoint,'phenotype',day)];c=cache[(endpoint,comparator,day)];v=f[5]-c[5];lo,hi=np.nanquantile(v,[.025,.975])
   pairs.append(dict(endpoint=endpoint,day=day,comparison='phenotype_minus_'+comparator,auc_difference=f[2]-c[2],conf_low=lo,conf_high=hi))
def write(name,rows):
 with (out/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
for name,rows in [('pooled_auc.csv',metrics),('pooled_roc.csv',rocs),('paired_auc_differences.csv',pairs),('cohort_counts.csv',den),('site_validation.csv',site_audit)]:write(name,rows)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
N=next(r['n'] for r in metrics)
site_ns={r['site']:int(r['n']) for r in den if r['endpoint']=='death_or_hospice'}
for endpoint,title in [('death_or_hospice','Hospital death or hospice'),('death','Hospital death')]:
 fig,axes=plt.subplots(1,3,figsize=(15,7));fig.subplots_adjust(left=.065,right=.98,bottom=.29,top=.73,wspace=.22)
 fig.suptitle('Pooled held-out ROC curves at 7, 14 and 28 days',x=.035,ha='left',fontsize=21,y=.97)
 fig.text(.035,.90,f'{title} after the 72-hour landmark | identical evaluation hospitalizations for all four models',fontsize=12)
 fig.text(.035,.85,f'Complete modified SOFA + adapted APACHE: N={N:,} | '+', '.join(f'{s} {n}' for s,n in site_ns.items()),fontsize=11)
 for ax,day in zip(axes,[7,14,28]):
  counts=next(r for r in metrics if r['endpoint']==endpoint and r['day']==day)
  for model,label,color in zip(models,labels,colors):
   fp,tp,auc,ci,band,_=cache[(endpoint,model,day)]
   ax.plot(fp,tp,color=color,lw=2,label=f'{label}\nAUC {auc:.3f} ({ci[0]:.3f}–{ci[1]:.3f})')
  ax.plot([0,1],[0,1],ls='--',color='#BBBBBB',lw=1)
  ax.set(xlim=(0,1),ylim=(0,1),xlabel='1 − specificity',title=f'Day {day} | {counts["cases"]} cases, {counts["controls"]} controls');ax.set_aspect('equal');ax.grid(alpha=.12)
  ax.legend(frameon=False,loc='lower right',fontsize=7.4)
 axes[0].set_ylabel('Sensitivity')
 note=('Each site is evaluated using predictors frozen from the other two sites; no held-out-site refitting. Counts at common probability bins are summed across sites.\n'
       'Controls include competing hospital exits and those still free of the endpoint. Scores are frozen 28-day risks evaluated at each horizon.\n'
       'AUCs use scores grouped into fixed 0.00001-wide bins; parentheses: pointwise 95% intervals from 200 paired within-site patient-cluster bootstrap resamples.\n'
       'Unadjusted combined-cohort ROC; larger sites contribute more observations. Modified SOFA and adapted APACHE are CLIF-derived proxies.')
 fig.text(.035,.045,note,fontsize=9,color='#444444',va='bottom')
 for ext in ['png','pdf']:fig.savefig(out/f'pooled_roc_{endpoint}.{ext}',dpi=180,facecolor='white')
 plt.close(fig)
manifest=dict(status='passed',sites=[m['site'] for m in meta],consensus_id=meta[0]['consensus_id'],n=N,site_n=site_ns,patient_level_inputs_read=False,method='ROC of pooled exact counts in common 0.00001 score bins; not average site ROCs',cohort='Same matched complete-SOFA/adapted-APACHE patients across all four frozen held-out models',bootstrap=200,bootstrap_scope='Independent within-site patient-cluster resamples; same replicate across models/days/endpoints; conditional on frozen learners',unrounded_auc_ambiguity_bound_max=max(r['quantization_bound'] for r in metrics),local_exact_vs_binned_auc_max=max(float(r['absolute_difference']) for r in site_audit if r['absolute_difference']),limitations=['Complete-case selection; restricted sample','Mixture of site-held-out predictors, not one common fitted predictor','Combined-site case mix and calibration affect pooled ROC; site curves retained','Modified SOFA and adapted APACHE are not validated standard scores'],inputs={str(Path(s)/'manifest.json'):hashlib.md5((Path(s)/'manifest.json').read_bytes()).hexdigest() for s in a.shares})
(out/'manifest.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest,indent=2))
