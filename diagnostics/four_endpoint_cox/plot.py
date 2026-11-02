"""Four-endpoint Table2/forest from aggregate coordinator output only."""
import argparse,csv,json,os,hashlib
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/clif_manuscript_mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('reports',type=Path);a=p.parse_args();root=a.reports
m=json.loads((root/'manifest.json').read_text());assert m['status']=='passed'
def read(p):
 with p.open() as f:return list(csv.DictReader(f))
rows=read(root/'table2_four_endpoint_core_cox.csv');meta=read(root/'primary/model_metadata_by_site.csv');part=read(root/'primary/event_partitions_by_site.csv')
endpoints=[('death','Hospital death'),('hospice','Hospice discharge'),('death_or_hospice','Hospital death or hospice'),('live_discharge','Live non-hospice discharge')]
colors={'P2':'#347DB1','P3':'#4FA343','P4':'#E65D0A'};names={'P2':'P2: Low-flow oxygen','P3':'P3: HFNC','P4':'P4: IMV'}
N=sum(int(r['n_lower']) for r in meta if r['model']=='core' and r['endpoint']=='death')
counts={ph:sum(int(r['n']) for r in part if r['phenotype']==ph) for ph in ['P1','P2','P3','P4','I']}
table=[]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig,axes=plt.subplots(2,2,figsize=(14,9));fig.subplots_adjust(left=.18,right=.97,top=.80,bottom=.23,hspace=.65,wspace=.55)
fig.suptitle('Four post-landmark cause-specific outcomes',x=.035,ha='left',fontsize=22,y=.97)
fig.text(.035,.91,f'Core-adjusted pooled Cox hazard ratios | P1 reference | {N:,} landmark hospitalizations | 28-day follow-up',fontsize=12)
for ax,(endpoint,title) in zip(axes.flat,endpoints):
 z=[r for r in rows if r['endpoint']==endpoint]
 event_n=sum(int(r['events_lower']) for r in meta if r['model']=='core' and r['endpoint']==endpoint)
 for ph,y in zip(['P2','P3','P4'],[2,1,0]):
  r=next(r for r in z if r['phenotype']==ph)
  if r['status']=='ok':
   hr,lo,hi=[float(r[k]) for k in ['hazard_ratio','hr_low','hr_high']]
   ax.plot([lo,hi],[y,y],color=colors[ph],lw=2);ax.plot(hr,y,'o',color=colors[ph],ms=7)
   label=f'{hr:.2f} ({lo:.2f}–{hi:.2f})'
  else:label='Not estimable'
  ax.text(.98,y+.20,label,transform=ax.get_yaxis_transform(),fontsize=9.5,va='center',ha='right',bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=1))
  table.append(dict(endpoint=endpoint,endpoint_label=title,phenotype=ph,reference='P1',cohort_n=N,phenotype_n=counts[ph],reference_n=counts['P1'],cohort_events=event_n,contributing_sites=r['k'],hazard_ratio=r.get('hazard_ratio',''),conf_low=r.get('hr_low',''),conf_high=r.get('hr_high',''),i2_percent=r.get('i2_percent',''),tau2=r.get('tau2',''),status=r['status']))
 ax.axvline(1,color='#888888',ls='--',lw=1);ax.set_xscale('log');ax.set(xlim=(.1,50),ylim=(-.5,2.5),yticks=[2,1,0],yticklabels=[names[p] for p in ['P2','P3','P4']],title=f'{title}\n{event_n:,} events')
 ax.set_xticks([.1,.25,.5,1,2,5,10,25,50]);ax.set_xticklabels(['0.1','0.25','0.5','1','2','5','10','25','50'],fontsize=9);ax.grid(axis='x',alpha=.13)
 ax.set_xlabel('Cause-specific hazard ratio')
fig.text(.035,.045,'Pooled site estimates: REML with modified Hartung–Knapp 95% intervals; each endpoint uses the same core covariates and cohort.\n'
         'Competing exits are censored at exit; follow-up ends at day28 after the ICU+72h landmark. I is retained in a separate baseline-hazard stratum.\n'
         'A discharge HR below1 indicates a lower live non-hospice discharge rate among patients still hospitalized. These are not absolute probabilities.\n'
         'Discharge P4 PH tests indicate time variation at all three sites; retain the exploratory days0–7 and7–28 estimates when interpreting these HRs.',fontsize=9,color='#444444')
for ext in ['png','pdf']:fig.savefig(root/f'four_endpoint_core_cox.{ext}',dpi=180,facecolor='white')
plt.close(fig)
with (root/'table2_presentation.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(table[0]));w.writeheader();w.writerows(table)
(root/'plot_manifest.json').write_text(json.dumps({'status':'passed','aggregate_inputs_only':True,'cohort_n':N,'phenotype_n':counts,'input_md5':{str(p):hashlib.md5(p.read_bytes()).hexdigest() for p in [root/'manifest.json',root/'table2_four_endpoint_core_cox.csv',root/'primary/model_metadata_by_site.csv',root/'primary/event_partitions_by_site.csv']}},indent=2))
print(root)
