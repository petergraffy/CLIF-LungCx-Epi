"""Pool locally cross-validated score histograms. Reads no patient records."""
import argparse, csv, hashlib, json, os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', '/tmp/clif_stage2_mpl')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('output',type=Path);p.add_argument('shares',nargs='+',type=Path);a=p.parse_args()
out=a.output;out.mkdir(parents=True,exist_ok=False)
meta=[];data={};den=[]
for share in a.shares:
    m=json.loads((share/'manifest.json').read_text());meta.append(m)
    assert m['status']=='passed' and m['validation']=='patient_grouped_local_cross_validation'
    assert not m['external_site_outcome_validation'] and m['bootstrap']==200 and m['bins']==100000
    for f,h in m['files'].items():assert hashlib.md5((share/f).read_bytes()).hexdigest()==h
    with (share/'cohort_counts.csv').open() as f:den+=list(csv.DictReader(f))
    with (share/'score_histograms.csv').open() as f:
        for r in csv.DictReader(f):
            assert r['site']==m['site']
            key=(r['cohort'],r['endpoint'],r['model']);key2=(int(r['bin']),int(r['case_day']))
            assert 0<=key2[0]<100000 and 1<=key2[1]<=29
            v=np.array([float(r['n'])]+[float(r[f'n_b{b:03d}']) for b in range(1,201)])
            assert np.all(v>=0) and np.all(v==np.floor(v))
            h=data.setdefault(key,{});h[key2]=h.get(key2,np.zeros(201))+v
assert len({m['site'] for m in meta})==len(meta) and len({m['consensus_id'] for m in meta})==1

def auc(c,n):
    den=c.sum(axis=0)*n.sum(axis=0)
    num=(c*(np.cumsum(n,axis=0)-.5*n)).sum(axis=0)
    return np.divide(num,den,out=np.full_like(num,np.nan),where=den>0)
metrics=[];curves=[];cache={}
for key,h in sorted(data.items()):
    cohort,endpoint,model=key;bins=sorted({k[0] for k in h});index={v:i for i,v in enumerate(bins)}
    total=np.zeros((len(bins),201))
    for (ib,day),v in h.items():total[index[ib]]+=v
    for day in range(1,29):
        cases=np.zeros_like(total)
        for (ib,cd),v in h.items():
            if cd<=day:cases[index[ib]]+=v
        controls=total-cases;av=auc(cases,controls);bs=av[1:];good=np.isfinite(bs)
        lo,hi=np.quantile(bs[good],[.025,.975]) if good.any() else (np.nan,np.nan)
        nc=cases[:,0].sum();nn=controls[:,0].sum()
        bound=float((cases[:,0]*controls[:,0]).sum()/(2*nc*nn)) if nc>0 and nn>0 else np.nan
        metrics.append(dict(cohort=cohort,endpoint=endpoint,model=model,day=day,n=int(nc+nn),cases=int(nc),controls=int(nn),auc=av[0],conf_low=lo,conf_high=hi,bootstrap_replicates=int(good.sum()),quantization_bound=bound,status='ok' if np.isfinite(av[0]) else 'nonestimable'))
        cache[key+(day,)]=(av,lo,hi)
        if day in (7,14,28) and nc>0 and nn>0:
            tpr=np.r_[0,np.cumsum(cases[::-1,0])/nc];fpr=np.r_[0,np.cumsum(controls[::-1,0])/nn]
            for x,y in zip(fpr,tpr):curves.append(dict(cohort=cohort,endpoint=endpoint,model=model,day=day,fpr=x,sensitivity=y))
# Verify the four models use exactly the same evaluation cohort at each horizon.
models=['baseline','baseline_sofa','baseline_adapted','phenotype']
for endpoint in ('death','death_or_hospice'):
    for day in (7,14,28):
        rows=[r for r in metrics if r['cohort']=='matched_complete_SOFA' and r['endpoint']==endpoint and r['day']==day and r['model'] in models]
        assert len(rows)==4 and len({(r['n'],r['cases'],r['controls']) for r in rows})==1
pairs=[]
for endpoint in ('death','death_or_hospice'):
    for day in (7,14,28):
        basekey=('matched_complete_SOFA',endpoint)
        for model in models[:-1]:
            v=cache[basekey+('phenotype',day)][0]-cache[basekey+(model,day)][0];bs=v[1:];good=np.isfinite(bs)
            lo,hi=np.quantile(bs[good],[.025,.975]) if good.any() else (np.nan,np.nan)
            pairs.append(dict(cohort=basekey[0],endpoint=endpoint,day=day,comparison='phenotype_minus_'+model,auc_difference=v[0],conf_low=lo,conf_high=hi))
def write(name,rows):
    if not rows:return
    with (out/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
write('pooled_tdauc.csv',metrics);write('pooled_roc.csv',curves);write('paired_auc_differences.csv',pairs);write('site_cohort_counts.csv',den)
labels=['Baseline','Baseline + SOFA','Baseline + adapted APACHE','Baseline + phenotype'];colors=['#646464','#32916B','#9D5A98','#357BB5']
for endpoint in ('death','death_or_hospice'):
    fig,axes=plt.subplots(1,3,figsize=(14,4))
    for ax,day in zip(axes,(7,14,28)):
        for model,label,color in zip(models,labels,colors):
            z=[r for r in curves if r['cohort']=='matched_complete_SOFA' and r['endpoint']==endpoint and r['model']==model and r['day']==day]
            if z:
                val=cache[('matched_complete_SOFA',endpoint,model,day)][0][0]
                ax.plot([r['fpr'] for r in z],[r['sensitivity'] for r in z],color=color,label=f'{label} ({val:.3f})')
        ax.plot([0,1],[0,1],ls='--',color='#bbbbbb');ax.set(xlim=(0,1),ylim=(0,1),xlabel='1 − specificity',ylabel='Sensitivity');ax.text(.03,.97,f'Day {day}',transform=ax.transAxes,va='top');ax.legend(frameon=False,fontsize=7)
    fig.tight_layout();fig.savefig(out/f'pooled_roc_{endpoint}.png',dpi=180);fig.savefig(out/f'pooled_roc_{endpoint}.svg');plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4))
    for model,label,color in zip(models,labels,colors):
        z=[r for r in metrics if r['cohort']=='matched_complete_SOFA' and r['endpoint']==endpoint and r['model']==model]
        ax.plot([r['day'] for r in z],[r['auc'] for r in z],label=label,color=color)
        ax.fill_between([r['day'] for r in z],[r['conf_low'] for r in z],[r['conf_high'] for r in z],color=color,alpha=.08)
    ax.set(xlim=(1,28),ylim=(0,1),xlabel='Days after 72-hour landmark',ylabel='Time-dependent AUC');ax.legend(frameon=False);fig.tight_layout();fig.savefig(out/f'pooled_tdauc_{endpoint}.png',dpi=180);fig.savefig(out/f'pooled_tdauc_{endpoint}.svg');plt.close(fig)
(out/'pooling_manifest.json').write_text(json.dumps(dict(status='passed',sites=[m['site'] for m in meta],consensus_id=meta[0]['consensus_id'],validation='patient_grouped_local_cross_validation',external_site_outcome_validation=False,bootstrap='200 within-site patient-cluster resamples of fixed OOF predictions; conditional on fitted folds and frozen phenotypes',pooling='Sum counts by common probability bin; cumulative/dynamic AUC2; competing exits are controls',ordinary_censoring='valid-outcome landmark hospitalizations only; no censoring before horizon',score='Out-of-fold 28-day risk marker evaluated over days 1-28',site_model_heterogeneity='Each site fits its own CV models; pooled discrimination is not external validation of a single shared predictor'),indent=2))
