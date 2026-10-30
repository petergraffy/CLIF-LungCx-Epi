"""Shareable Rush-only respiratory coverage diagnostic; aggregate data only."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/clif_manuscript_mpl')
from pathlib import Path
import csv,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
ROOT=Path(__file__).resolve().parents[2];RUN=ROOT/'output/run_RUSH_respiratory_gap_20261030';OUT=RUN/'share'
def read(name):
 with (OUT/name).open() as f:return list(csv.DictReader(f))
s={r['metric']:float(r['n']) for r in read('respiratory_gap_summary.csv')}
icu={r['metric']:float(r['value']) for r in read('icu_exposure_summary.csv')}
year=read('gap_by_admission_year.csv');m=json.loads((RUN/'audit_manifest.json').read_text());assert m['status']=='passed'
for k,v in m['files'].items():assert hashlib.md5((OUT/k).read_bytes()).hexdigest()==v
N=int(s['landmark_hospitalizations']);U=int(s['all72h_unknown_support']);no=int(s['no_resp_rows_first72h_among_all_unknown'])
assert N==969 and U==126 and no==125
blue='#367C9B';orange='#D46B35';light='#E6EEF2';gray='#64747B'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','pdf.fonttype':42})
fig=plt.figure(figsize=(15,10.5));gs=fig.add_gridspec(2,2,height_ratios=[1,1.48],left=.085,right=.97,top=.78,bottom=.24,wspace=.27,hspace=.67)
fig.suptitle('Rush CLIF: respiratory-support coverage gap',x=.045,y=.97,ha='left',fontsize=24)
fig.text(.045,.913,f'{U:,} of {N:,} landmark hospitalizations ({100*U/N:.1f}%) have no classifiable support state throughout the first 72 hours.',fontsize=14,color='#293B46')
fig.text(.045,.869,'125 have no respiratory-table rows anywhere in the hospitalization; other clinical data are present.',fontsize=13,color=orange,weight='bold')
ax=fig.add_subplot(gs[0,0]);ax.set_title('A   Support-state coverage after waterfall',loc='left',pad=18,fontsize=13)
known=N-U
ax.barh([0],[known],color=blue,height=.65);ax.barh([0],[U],left=[known],color=orange,height=.65)
ax.text(known/2,0,f'≥1 known-support hour\n{known:,} / {N:,} ({100*known/N:.1f}%)',ha='center',va='center',color='white',fontsize=12)
ax.text(known+U/2,0,f'{U}\n({100*U/N:.1f}%)',ha='center',va='center',color='white',fontsize=11,weight='bold')
ax.annotate('All 72 hours unknown',xy=(known+U/2,-.34),xytext=(known+U/2,-.77),ha='center',fontsize=10,color=orange,arrowprops={'arrowstyle':'-','color':orange})
ax.set(xlim=(0,N),ylim=(-1,1),yticks=[],xticks=[0,250,500,750,969]);ax.set_xlabel('Hospitalizations (N=969)');ax.spines['left'].set_visible(False)
ax=fig.add_subplot(gs[0,1]);ax.set_title('B   Source-table audit in the 126 unknown stays',loc='left',pad=18,fontsize=13)
ax.barh([0],[no],color=orange,height=.65);ax.barh([0],[U-no],left=[no],color=gray,height=.65)
ax.text(no/2,0,f'No respiratory rows\n{no} / {U} (99.2%)',ha='center',va='center',color='white',fontsize=12,weight='bold')
ax.annotate('1 stay has 94 rows;\ndevice category missing on all 94',xy=(125.5,.34),xytext=(115,1.03),ha='right',va='bottom',fontsize=10,color=gray,arrowprops={'arrowstyle':'-','color':gray})
ax.set(xlim=(0,U),ylim=(-1,1.5),yticks=[],xticks=[0,30,60,90,126]);ax.set_xlabel('Wholly unknown hospitalizations (N=126)');ax.spines['left'].set_visible(False)
ax.text(0,-.82,'The same 125 stays have no respiratory rows during the entire hospitalization.',fontsize=9,color='#555555',ha='left',va='top')
ax=fig.add_subplot(gs[1,0]);ax.set_title('C   Gap occurs across the supplied admission years',loc='left',pad=18,fontsize=13)
x=np.arange(len(year));pct=np.array([float(r['all72h_unknown_pct']) for r in year])
ax.bar(x,pct,color=orange,width=.62)
for i,(p,r) in enumerate(zip(pct,year)):
 ax.text(i,p+.55,f'{p:.1f}%',ha='center',fontsize=10)
ax.axhline(100*U/N,color=gray,ls='--',lw=1)
ax.text(.98,.96,'Dashed line: overall 13.0%',transform=ax.transAxes,fontsize=9,color=gray,ha='right',va='top')
ax.set_xticks(x);ax.set_xticklabels([r['admission_year']+'\n'+r['all72h_unknown_n']+'/'+r['n'] for r in year],fontsize=10)
ax.set_ylim(0,24);ax.yaxis.set_major_formatter(PercentFormatter(100,decimals=0));ax.set_ylabel('Stays with all 72 hours unknown')
ax.set_xlabel('Admission year in supplied data\nUnknown stays / all landmark stays');ax.grid(axis='y',alpha=.12);ax.set_axisbelow(True)
ax=fig.add_subplot(gs[1,1]);ax.set_title('D   Other clinical records in the 126 unknown stays',loc='left',pad=18,fontsize=13)
vals=[int(s['any_spo2_first72h_among_all_unknown']),int(s['any_numeric_lab_first72h_among_all_unknown']),int(icu['at_least24_closed_ICU_hours_n'])]
labels=['SpO₂ records','Lab results','≥24 ICU hours']
ax.barh(range(3),[U]*3,color=light,height=.55);ax.barh(range(3),vals,color=blue,height=.55)
for j,v in enumerate(vals):ax.text(v-2,j,f'{v}/{U}',ha='right',va='center',color='white',weight='bold',fontsize=11)
ax.set_yticks(range(3));ax.set_yticklabels(labels,fontsize=10);ax.invert_yaxis();ax.set_xlim(0,U);ax.set_xticks([0,30,60,90,126]);ax.set_xlabel('Hospitalizations with evidence in the first 72 hours')
ax.text(0,-.32,f'Median documented ICU exposure: {icu["closed_ICU_interval_hours_median"]:.1f} hours\n(IQR {icu["closed_ICU_interval_hours_q1"]:.1f}–{icu["closed_ICU_interval_hours_q3"]:.1f}; 47 stays cover all 72 ICU hours).',transform=ax.transAxes,fontsize=10,color='#444444',va='top')
fig.text(.045,.117,'For Rush source / ETL review',fontsize=12,weight='bold',color='#293B46')
fig.text(.045,.083,'Check respiratory flowsheet extraction, admission-ID/time joins, and explicit room-air/off-oxygen device records.\nCompare the affected stays against the source EHR before attributing the gap to ETL versus absent source documentation.',fontsize=11,color='#444444',va='center')
fig.text(.045,.026,'Rush-only aggregate diagnostic; no patient identifiers. Window: first 72 hours after first ICU arrival; patients remain hospitalized at the landmark.\n“Unknown” does not mean room air. ICU exposure uses closed ADT intervals. Analyzed Rush CLIF 2.1 extract: October 22 cache; figure generated October 30, 2026.\nNo newer Rush release or source EHR extract was checked in this diagnostic.',fontsize=9,color='#666666',va='bottom')
for ext in ['png','pdf']:fig.savefig(OUT/f'RUSH_respiratory_support_coverage_gap.{ext}',dpi=200,facecolor='white')
plt.close(fig)
caption=f'''Rush CLIF respiratory-support completeness review

In the analyzed Rush CLIF 2.1 extract, {U} of {N} lung-cancer ICU landmark hospitalizations ({100*U/N:.1f}%) have unknown support for the entire first 72 hours despite waterfall preprocessing. Of these, {no} have no respiratory-table records anywhere in the hospitalization; one has 94 records, all with missing device category. All 126 have SpO2 and numeric laboratory records in the 72-hour window. Median documented ICU exposure is {icu['closed_ICU_interval_hours_median']:.1f} hours (IQR {icu['closed_ICU_interval_hours_q1']:.1f}–{icu['closed_ICU_interval_hours_q3']:.1f}); 118 have at least 24 ICU hours and 47 cover the full 72 hours.

The pattern occurs across admission years 2018–2024 in the supplied data. Please check source respiratory flowsheets, device/room-air inclusion, hospitalization-ID linkage and timestamp filtering. These aggregates demonstrate a respiratory-table coverage gap, but do not establish whether ETL omissions or absent source documentation caused it. Unknown support must not be assumed to be room air.

Source: frozen Rush CLIF 2.1 October 22 analysis cache used in the current three-site study; diagnostic generated October 30, 2026. A newer release and source EHR were not examined. This bundle contains only Rush aggregate counts.
'''
(OUT/'caption_and_team_review.txt').write_text(caption)
print(OUT/'RUSH_respiratory_support_coverage_gap.png')
