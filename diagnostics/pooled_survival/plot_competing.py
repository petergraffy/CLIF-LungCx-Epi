"""Phenotype-stratified Aalen-Johansen curves from exact pooled aggregate counts."""
import argparse,csv,hashlib,json,os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/clif_manuscript_mpl')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
p=argparse.ArgumentParser();p.add_argument('reports',type=Path);a=p.parse_args();out=a.reports.resolve()
source=out/'pooled_competing_incidence.csv';manifest=json.loads((out/'pooling_manifest.json').read_text())
assert manifest['status']=='passed' and manifest['probability_partition_checked']
with source.open() as f:rows=list(csv.DictReader(f))
phs=['P1','P2','P3','P4','I'];colors=['#648BA7','#347DB1','#4FA343','#E65D0A','#777777']
names=['Room-air predominant','Low-flow oxygen','HFNC','IMV','Insufficient information']
def select(ph,endpoint):return sorted([r for r in rows if r['curve']=='landmark72h_'+ph and r['endpoint']==endpoint],key=lambda r:float(r['time']))
def values(z,key):return np.array([float(r[key]) if r[key] else np.nan for r in z])
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
fig=plt.figure(figsize=(14,8));gs=fig.add_gridspec(2,2,height_ratios=[3,1.05],left=.08,right=.98,top=.92,bottom=.06,hspace=.44,wspace=.22)
summary=[]
for j,(endpoint,title) in enumerate([('death_or_hospice','A   Hospital death or discharge to hospice'),('live_discharge','B   Discharge alive without hospice')]):
 ax=fig.add_subplot(gs[0,j]);risks=[]
 for ph,color,name in zip(phs,colors,names):
  z=select(ph,endpoint);assert z
  x=values(z,'time');y=values(z,'estimate');assert np.all(np.diff(y[np.isfinite(y)])>=-1e-12)
  ax.step(x,y,where='post',color=color,label=f'{ph}: {name} (N={int(z[0]["n"]):,})',lw=2,ls='--' if ph=='I' else '-')
  ax.fill_between(x,values(z,'conf_low'),values(z,'conf_high'),step='post',color=color,alpha=.09)
  risks.append([f'{int(next(r["n_risk"] for r in z if float(r["time"])==t)):,}' for t in [0,7,14,21,28]])
  for day in [7,14,28]:
   r=next(r for r in z if float(r['time'])==day)
   summary.append({k:r[k] for k in ['curve','endpoint','time','n','n_risk','estimate','conf_low','conf_high']})
 ax.set(xlim=(0,28),ylim=(0,1),xticks=[0,7,14,21,28],xlabel='Days after 72-hour landmark',ylabel='Cumulative incidence')
 ax.set_title(title,fontsize=12,pad=14);ax.yaxis.set_major_formatter(PercentFormatter(1));ax.grid(alpha=.13);ax.legend(frameon=False,fontsize=8,loc='upper left')
 ra=fig.add_subplot(gs[1,j]);ra.axis('off');ra.set_title('Number at risk',fontsize=10,pad=10)
 table=ra.table(cellText=risks,rowLabels=phs,colLabels=['0','7','14','21','28'],loc='center',cellLoc='center');table.auto_set_font_size(False);table.set_fontsize(9);table.scale(1,1.3)
 for cell in table.get_celld().values():cell.set_edgecolor('#eeeeee')
for ext in ['png','svg']:fig.savefig(out/f'pooled_phenotype_competing_incidence.{ext}',dpi=200,facecolor='white',bbox_inches='tight')
plt.close(fig)
with (out/'competing_incidence_days7_14_28.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
audit=dict(status='passed',consensus_id=manifest['consensus_id'],clock='72-hour ICU landmark',groups=phs,primary_n=sum(int(select(ph,'death_or_hospice')[0]['n']) for ph in phs),method='Aalen-Johansen cumulative incidence; death/hospice and live non-hospice discharge compete',uncertainty='200 within-site patient-cluster bootstrap resamples; pointwise percentile 95% intervals',unknown_dispositions='Censored at observed hospital exit',patient_level_inputs_read=False,captions_and_subtitles=False,source_md5=hashlib.md5(source.read_bytes()).hexdigest())
(out/'competing_figure_manifest.json').write_text(json.dumps(audit,indent=2));print(json.dumps(audit,indent=2))
