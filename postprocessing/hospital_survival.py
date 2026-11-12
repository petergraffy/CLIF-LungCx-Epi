"""Plot exact pooled cause-specific KM from aggregate event/risk counts."""
import argparse,csv,hashlib,json,os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/clif_current_survival_mpl')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
p=argparse.ArgumentParser();p.add_argument('reports',type=Path);a=p.parse_args();out=a.reports
m=json.loads((out/'pooling_manifest.json').read_text());assert m['status']=='passed'
with (out/'pooled_survival.csv').open() as f:rows=list(csv.DictReader(f))
def select(curve):return sorted([r for r in rows if r['curve']==curve and r['endpoint']=='death'],key=lambda r:float(r['time']))
def vals(z,k):return np.array([float(r[k]) if r[k] else np.nan for r in z])
phs=['P1','P2','P3','P4'];names=['Shift toward room air','Sustained low-flow oxygen','High-flow predominant','Sustained invasive ventilation'];colors=['#7B3FA0','#C93882','#B5C83C','#35251D']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
fig=plt.figure(figsize=(14,8));gs=fig.add_gridspec(2,2,height_ratios=[3,1.05],left=.075,right=.98,top=.92,bottom=.07,hspace=.42,wspace=.22)
summary=[]
for j,title in enumerate(['A   All eligible hospitalizations','B   72-hour landmark phenotypes']):
 ax=fig.add_subplot(gs[0,j]);risk=[];labels=[]
 curves=[('icu_admission_Overall','Overall','#333333','All eligible')] if j==0 else [(f'landmark72h_{ph}',ph,c,n) for ph,c,n in zip(phs,colors,names)]
 for curve,ph,color,name in curves:
  z=select(curve);assert z;n=int(z[0]['n'])
  if not n:continue
  x=vals(z,'time');y=vals(z,'survival');assert y[0]==1 and np.all(np.diff(y[np.isfinite(y)])<=1e-12)
  ax.step(x,y,where='post',lw=2,color=color,ls='--' if ph=='I' else '-',label=f'{ph}: {name} (N={n:,})' if j else f'All eligible (N={n:,})')
  ax.fill_between(x,vals(z,'conf_low'),vals(z,'conf_high'),step='post',color=color,alpha=.09)
  labels.append(ph);risk.append([f'{int(next(r["n_risk"] for r in z if float(r["time"])==t)):,}' for t in [0,7,14,21,28]])
  for t in [7,14,28]:
   r=next(r for r in z if float(r['time'])==t);summary.append({k:r[k] for k in ['curve','time','n','n_risk','survival','conf_low','conf_high']})
 ax.set(xlim=(0,28),ylim=(0,1),xticks=[0,7,14,21,28],xlabel='Days after first ICU admission' if j==0 else 'Days after 72-hour landmark',ylabel='Hospital-death Kaplan–Meier estimate')
 ax.set_title(title,fontsize=12,pad=14);ax.yaxis.set_major_formatter(PercentFormatter(1));ax.grid(alpha=.13);ax.legend(frameon=False,fontsize=8,loc='lower left')
 ra=fig.add_subplot(gs[1,j]);ra.axis('off');ra.set_title('Number at risk',fontsize=10,pad=10)
 tb=ra.table(cellText=risk,rowLabels=labels,colLabels=['0','7','14','21','28'],loc='center',cellLoc='center');tb.auto_set_font_size(False);tb.set_fontsize(9);tb.scale(1,1.3)
 for cell in tb.get_celld().values():cell.set_edgecolor('#eeeeee')
for ext in ['png','svg']:fig.savefig(out/f'pooled_hospital_survival.{ext}',dpi=200,facecolor='white',bbox_inches='tight')
plt.close(fig)
with (out/'hospital_survival_days7_14_28.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
(out/'hospital_survival_figure_manifest.json').write_text(json.dumps(dict(status='passed',sites=m['sites'],consensus_id=m['consensus_id'],method='Cause-specific Kaplan-Meier for hospital death; hospice and live discharge censored at exit; unknown dispositions censored at observed exit',post_discharge_mortality_measured=False,clock='All eligible at ICU t0; phenotype groups enter at 72h landmark',uncertainty='200 within-site patient-cluster bootstrap resamples; pointwise percentile 95% intervals',source_csv_md5=hashlib.md5((out/'pooled_survival.csv').read_bytes()).hexdigest(),captions_and_subtitles=False),indent=2))
