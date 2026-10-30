"""Refresh the three-site annual missing-support figure with corrected UCMC."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/clif_manuscript_mpl')
from pathlib import Path
import csv,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'output/run_three_site_missingness_20261030';OUT.mkdir(exist_ok=True)
old=ROOT/'output/run_multisite_comparison_20261022/aggregate/unknown_by_year.csv'
new=ROOT/'output/run_UCMC_waterfall_20261026/aggregate/unknown_by_year.csv'
def read(p):
 with p.open() as f:return list(csv.DictReader(f))
rows=[r for r in read(old) if r['site']!='UCMC' and r['analysis_id']=='h72_primary']
rows += [dict(r,site='UCMC') for r in read(new) if r['analysis_id']=='h72_primary']
for r in rows:
 for k in ['admission_year','hours','original_unknown_hours','waterfall_unknown_hours']:r[k]=int(r[k])
 r['original_unknown_pct']=100*r['original_unknown_hours']/r['hours'];r['waterfall_unknown_pct']=100*r['waterfall_unknown_hours']/r['hours']
 r['hospitalizations']=r['hours']//72;assert r['hours']%72==0
source_files=[old,new];totals={}
for site in ['NU','RUSH','UCMC']:
 z=[r for r in rows if r['site']==site];hours=sum(r['hours'] for r in z);unknown=sum(r['waterfall_unknown_hours'] for r in z)
 exact=ROOT/f'output/run_exact_aggregates_20261030/sites/{site}/counts/profiles/hourly_support.csv';source_files.append(exact)
 p=read(exact);assert sum(int(r['n']) for r in p)==hours
 assert sum(int(r['n']) for r in p if r['state']=='UNKNOWN')==unknown
 totals[site]={'n':hours//72,'hours':hours,'unknown_hours':unknown,'waterfall_unknown_pct':100*unknown/hours}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig,axs=plt.subplots(1,3,figsize=(15,6.7),sharey=True)
fig.subplots_adjust(left=.075,right=.98,top=.68,bottom=.24,wspace=.12)
fig.suptitle('Unknown respiratory support across sites — updated UCMC',x=.04,y=.965,ha='left',fontsize=22)
fig.text(.04,.89,'First 72 hours after ICU arrival | same landmark cohorts | corrected UCMC source; NU and Rush unchanged',fontsize=12,color='#444444')
for ax,site in zip(axs,['NU','RUSH','UCMC']):
 z=sorted([r for r in rows if r['site']==site],key=lambda r:r['admission_year']);x=[r['admission_year'] for r in z]
 ax.plot(x,[r['original_unknown_pct'] for r in z],'-o',color='#A0A0A0',lw=2,ms=5)
 ax.plot(x,[r['waterfall_unknown_pct'] for r in z],'-o',color='#2469A6',lw=2.4,ms=5)
 t=totals[site];ax.set_title(f'{site} · N={t["n"]:,}\nAfter waterfall overall: {t["waterfall_unknown_pct"]:.2f}%',fontsize=12,pad=15)
 ax.set_xticks(x);ax.set_xticklabels(x,rotation=45);ax.set_ylim(-.8,30);ax.set_yticks([0,5,10,15,20,25,30]);ax.yaxis.set_major_formatter(PercentFormatter(100,decimals=0));ax.grid(alpha=.17)
 ax.set_xlabel('Admission year in supplied data')
axs[0].set_ylabel('Unknown support (% of hours)')
fig.legend(handles=[Line2D([],[],color='#2469A6',lw=2.4,marker='o',label='CLIF waterfall'),Line2D([],[],color='#A0A0A0',lw=2,marker='o',label='Original six-hour carry')],loc='upper center',bbox_to_anchor=(.54,.83),ncol=2,frameon=False,fontsize=11)
fig.text(.04,.055,'Year-specific denominators are 72 × the number of landmark hospitalizations admitted that year; exact counts are in the accompanying CSV.\nWaterfall uses bounded filling/inference; improved coverage does not equal additional directly observed documentation. Unknown is not assumed to be room air.\nUses the analyzed three-site extracts, including the October 26 UCMC refresh; no newer NU/Rush release was checked for this figure.',fontsize=9,color='#555555',va='bottom')
for ext in ['png','pdf']:fig.savefig(OUT/f'unknown_support_by_site_year_updated.{ext}',dpi=200,facecolor='white')
plt.close(fig)
with (OUT/'unknown_support_by_site_year_updated.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(OUT/'figure_audit.json').write_text(json.dumps({'status':'passed','counts_match_current_frozen_hourly_profiles':True,'totals':totals,'source_md5':{str(p.relative_to(ROOT)):hashlib.md5(p.read_bytes()).hexdigest() for p in source_files},'aggregate_only':True,'clinical_models_refitted':False},indent=2))
print(OUT/'unknown_support_by_site_year_updated.png')
