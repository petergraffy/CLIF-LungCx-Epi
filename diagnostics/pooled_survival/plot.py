"""Pooled KM examples from aggregate coordinator tables only."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/clif_manuscript_mpl')
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[2];OUT=Path(os.environ.get('CLIF_POOLED_SURVIVAL_OUTPUT',ROOT/'output/run_pooled_survival_20261030_v2')).resolve()/'reports'
def read(name):
 with (OUT/name).open() as f:return list(csv.DictReader(f))
def n(x):return float(x) if x else np.nan
def select(data,**kwargs):return [r for r in data if all(r.get(k)==str(v) for k,v in kwargs.items())]
pooled=read('pooled_survival.csv');site=read('site_survival.csv');discharge=read('pooled_discharge_incidence.csv');site_discharge=read('site_discharge_incidence.csv');m=json.loads((OUT/'pooling_manifest.json').read_text());assert m['status']=='passed'
PH=['P1','P2','P3','P4','I'];PC=dict(zip(PH,['#648BA7','#347DB1','#4FA343','#E65D0A','#777777']))
PN=['Room-air predominant','Low-flow oxygen','HFNC','IMV','Insufficient information']
SITES=['UCMC','NU','RUSH'];SC=dict(zip(SITES,['#247BA0','#D76345','#33966A']))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','pdf.fonttype':42})
def line(ax,z,color,label,ci=True,ls='-',lw=2,value='survival'):
 z=sorted(z,key=lambda r:n(r['time']));x=[n(r['time']) for r in z];y=[n(r[value]) for r in z]
 ax.step(x,y,where='post',color=color,label=label,lw=lw,ls=ls)
 if ci:ax.fill_between(x,[n(r['conf_low']) for r in z],[n(r['conf_high']) for r in z],color=color,alpha=.075,step='post')
 ax.set(xlim=(0,28),ylim=(0,1.02),xticks=[0,7,14,21,28]);ax.yaxis.set_major_formatter(PercentFormatter(1));ax.grid(alpha=.13)
def risk_table(ax,entries):
 ax.axis('off');values=[]
 for label,z in entries:
  values.append([f'{int(n(next(r["n_risk"] for r in z if n(r["time"])==t))):,}' for t in [0,7,14,21,28]])
 tab=ax.table(cellText=values,rowLabels=[label for label,_ in entries],colLabels=['0','7','14','21','28'],cellLoc='center',loc='center')
 tab.auto_set_font_size(False);tab.set_fontsize(9);tab.scale(1,1.35)
 for cell in tab.get_celld().values():cell.set_edgecolor('#EEEEEE')
 ax.set_title('Exact number at risk',fontsize=10,pad=10)
def save(fig,name,caption):
 fig.text(.035,.025,caption,fontsize=9,color='#444444',va='bottom')
 for ext in ['png','pdf']:fig.savefig(OUT/f'{name}.{ext}',dpi=180,facecolor='white')
 plt.close(fig)
fig=plt.figure(figsize=(14,9));gs=fig.add_gridspec(2,2,height_ratios=[3,1.1],left=.09,right=.98,top=.82,bottom=.19,hspace=.42,wspace=.20)
fig.suptitle('Pooled survival and discharge by 72-hour phenotype',x=.035,ha='left',fontsize=21,y=.97)
fig.text(.035,.91,'UCMC + NU + Rush | 3,310 landmark hospitalizations | exact pooled risk sets and event counts',fontsize=12)
for j,(endpoint,title) in enumerate([('death','A   Hospital-death Kaplan–Meier function'),('live_discharge','B   Discharge alive without hospice')]):
 ax=fig.add_subplot(gs[0,j])
 for ph,label in zip(PH,PN):
  z=select(discharge if endpoint=='live_discharge' else pooled,curve='landmark72h_'+ph,endpoint=endpoint);N=int(n(z[0]['n']))
  line(ax,z,PC[ph],f'{ph}: {label} (N={N:,})',ls='--' if ph=='I' else '-',value='estimate' if endpoint=='live_discharge' else 'survival')
 ax.set_title(title,fontsize=12,pad=16);ax.set_xlabel('Days after 72-hour landmark');ax.set_ylabel('Cumulative incidence' if endpoint=='live_discharge' else 'Kaplan–Meier estimate')
 ax.legend(frameon=False,fontsize=8,loc='upper left' if endpoint=='live_discharge' else 'lower left')
 risk_table(fig.add_subplot(gs[1,j]),[(ph,select(pooled,curve='landmark72h_'+ph,endpoint='death')) for ph in PH])
save(fig,'pooled_phenotype_survival','Unadjusted combined cohort: A, cause-specific KM; B, Aalen–Johansen probability of live non-hospice discharge, with death/hospice competing.\nShading: pointwise 95% intervals from 200 patient-cluster bootstrap resamples within each site; uncertainty is conditional on these three sites.\nA censors other hospital exits; B retains death/hospice as competing outcomes. Neither panel measures survival after discharge.\nI is an information category. Larger sites contribute more observations; site-specific curves and adjusted meta-analysis remain companions.')

fig=plt.figure(figsize=(14,9));gs=fig.add_gridspec(2,2,height_ratios=[3,1.1],left=.09,right=.98,top=.82,bottom=.19,hspace=.42,wspace=.20)
fig.suptitle('Pooled survival and discharge from ICU admission',x=.035,ha='left',fontsize=21,y=.97)
fig.text(.035,.91,'All 4,533 eligible hospitalizations enter at time zero; no 72-hour survivor requirement',fontsize=12)
for j,(endpoint,title) in enumerate([('death','A   Hospital-death Kaplan–Meier function'),('live_discharge','B   Discharge alive without hospice')]):
 ax=fig.add_subplot(gs[0,j])
 for s in SITES:line(ax,select(site_discharge if endpoint=='live_discharge' else site,curve='icu_admission_Overall',endpoint=endpoint,site=s),SC[s],s,ci=False,ls=':',lw=1.4,value='estimate' if endpoint=='live_discharge' else 'survival')
 line(ax,select(discharge if endpoint=='live_discharge' else pooled,curve='icu_admission_Overall',endpoint=endpoint),'#222222','Pooled (N=4,533)',lw=2.4,value='estimate' if endpoint=='live_discharge' else 'survival')
 ax.set_title(title,fontsize=12,pad=16);ax.set_xlabel('Days after ICU admission');ax.set_ylabel('Cumulative incidence' if endpoint=='live_discharge' else 'Kaplan–Meier estimate');ax.legend(frameon=False,fontsize=9,loc='upper left' if endpoint=='live_discharge' else 'lower left')
 entries=[('Pooled',select(pooled,curve='icu_admission_Overall',endpoint='death'))]+[(s,select(site,curve='icu_admission_Overall',endpoint='death',site=s)) for s in SITES]
 risk_table(fig.add_subplot(gs[1,j]),entries)
save(fig,'pooled_whole_cohort_survival','Solid black: exact pooled estimate; dotted lines: individual sites. A: cause-specific KM. B: Aalen–Johansen discharge incidence. Shading: pointwise 95% within-site patient-cluster bootstrap intervals.\nB: live non-hospice discharge; death/hospice compete. Unknown terminal dispositions are censored at exit; missing-outcome bias remains possible.\nThe unit is hospitalization. No identifiers or calendar timestamps are pooled, and patients cannot be linked across institutions.\nWhole-cohort admission-time survival is separate from the landmark phenotype analysis. These are unadjusted descriptive curves.')
print(OUT)
