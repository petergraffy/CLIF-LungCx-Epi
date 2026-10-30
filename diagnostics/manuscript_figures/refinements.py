"""Smoothed template displays, competing-event curves and absolute AUROC plots."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/clif_manuscript_mpl')
from pathlib import Path
import csv,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(os.environ.get('CLIF_FIGURE_OUTPUT',ROOT/'output/run_manuscript_figures_20261030_v2')).resolve();OUT.mkdir(exist_ok=True)
FED=Path(os.environ.get('CLIF_FEDERATED_OUTPUT',ROOT/'output/run_federated_rehearsal_20261026')).resolve()
REV=Path(os.environ.get('CLIF_REVISION_OUTPUT',ROOT/'output/run_revision_gaps_20261026')).resolve()/'reports/tables'
AP=Path(os.environ.get('CLIF_APACHE_OUTPUT',ROOT/'output/run_APACHE_adapted_20261030_v2')).resolve()/'reports/tables'
DESC=Path(os.environ.get('CLIF_DESCRIPTOR_REPORTS',ROOT/'output/run_descriptives_20261030_v2/reports')).resolve()
inputs=set()
def read(p):
 inputs.add(p)
 with p.open() as f:return list(csv.DictReader(f))
def n(v):return float(v) if v else np.nan
SITES=['UCMC','NU','RUSH'];SC=dict(zip(SITES,['#247BA0','#D76345','#33966A']))
PH=['P1','P2','P3','P4'];PN=['Room-air predominant','Low-flow oxygen predominant','HFNC predominant','IMV predominant']
PC=['#B8D4E8','#347DB1','#4FA343','#E65D0A','#A5A5A5']
# Curves need a darker P1 than the area fill for readable lines.
LC=dict(zip(PH+['I'],['#648BA7','#347DB1','#4FA343','#E65D0A','#777777']))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                     'axes.spines.right':False,'axes.titleweight':'bold','pdf.fonttype':42})
def save(fig,name,note):
 fig.text(.035,.024,note,fontsize=9,va='bottom',color='#444444')
 for ext in ['png','pdf']:fig.savefig(OUT/f'{name}.{ext}',dpi=180,facecolor='white')
 plt.close(fig)
def smooth_simplex(values,centers):
 """Common convex cubic interpolation: no negatives, overshoot or lost mass."""
 x=np.linspace(0,72,721);loc=np.clip(x,centers[0],centers[-1]);i=np.clip(np.searchsorted(centers,loc,side='right')-1,0,len(centers)-2)
 t=(loc-centers[i])/(centers[i+1]-centers[i]);w=3*t*t-2*t*t*t
 z=values[:,i]*(1-w)+values[:,i+1]*w
 assert np.all(z>=-1e-12) and np.allclose(z.sum(axis=0),1)
 assert np.allclose(smooth_at_centers(values,centers),values)
 return x,z
def smooth_at_centers(values,centers):
 # The interpolation reproduces every original interval-center value exactly.
 i=np.clip(np.searchsorted(centers,centers,side='right')-1,0,len(centers)-2)
 t=(centers-centers[i])/(centers[i+1]-centers[i]);w=3*t*t-2*t*t*t
 return values[:,i]*(1-w)+values[:,i+1]*w

den=read(DESC/'denominators_by_site.csv')
count=lambda s:int(next(r['n'] for r in den if r['site']==s and r['cohort']=='landmark_cohort' and r['group']=='Overall'))
pro=read(FED/'coordinator/full/prototypes.csv')
support=['ROOM_AIR','LOW_O2','HFNC','CPAP','NIPPV','IMV','OTHER','UNKNOWN']
support_colors=['#B8D4E8','#347DB1','#4FA343','#AA59AE','#7142A2','#E65D0A','#D6B565','#A5A5A5']
fig=plt.figure(figsize=(16,11));gs=fig.add_gridspec(4,4,height_ratios=[1.5,1,1,1.15],left=.085,right=.98,top=.85,bottom=.19,hspace=.53,wspace=.22)
fig.suptitle('Shared respiratory support–ARF phenotypes',x=.035,ha='left',fontsize=21,y=.975)
fig.text(.035,.93,'A   Smoothed display of the frozen 72-hour consensus templates',fontsize=12)
for row,(domain,states,colors,label) in enumerate([
 ('support',support,support_colors,'Respiratory\nsupport'),
 ('hypox',['H0','H1','HU'],['#B8D9B0','#C85142','#BBBBBB'],'Hypoxemia\nevidence'),
 ('hyper',['C0','C1','CU'],['#B8D9B0','#C85142','#BBBBBB'],'Hypercapnia/\nacidemia evidence')]):
 for j,ph in enumerate(PH):
  ax=fig.add_subplot(gs[row,j]);values=[]
  for st in states:
   z=sorted([r for r in pro if r['phenotype']==ph and r['domain']==domain and r['state']==st],key=lambda r:int(r['bin']))
   assert len(z)==12;values.append([n(r['probability']) for r in z])
  values=np.array(values);assert np.allclose(values.sum(axis=0),1)
  x,curves=smooth_simplex(values,np.arange(12)*6+3)
  ax.stackplot(x,curves,colors=colors,linewidth=0)
  ax.set(xlim=(0,72),ylim=(0,1),xticks=[0,24,48,72],yticks=[0,.5,1]);ax.yaxis.set_major_formatter(PercentFormatter(1));ax.tick_params(labelsize=9)
  if j==0:ax.set_ylabel(label,labelpad=18)
  else:ax.set_yticklabels([])
  if row==0:ax.set_title(ph+'\n'+PN[j],fontsize=11,pad=14)
  if row==2:ax.set_xlabel('Hours after ICU arrival')
  else:ax.set_xticklabels([])
ax=fig.add_subplot(gs[3,:]);left=np.zeros(3)
for ph,c in zip(PH+['I'],PC):
 ns=np.array([int(next(r['n'] for r in den if r['site']==s and r['cohort']=='primary_phenotypes' and r['group']==ph)) for s in SITES]);v=ns/np.array([count(s) for s in SITES])
 ax.barh(SITES,v,left=left,color=c,label=ph if ph!='I' else 'I: insufficient information',height=.62)
 for i,(nn,w) in enumerate(zip(ns,v)):
  if w>.055:ax.text(left[i]+w/2,i,f'{nn:,}',ha='center',va='center',fontsize=10,color='white' if ph in ['P2','P4'] else '#222222')
 left+=v
assert np.allclose(left,1)
ax.invert_yaxis();ax.set_xlim(0,1);ax.xaxis.set_major_formatter(PercentFormatter(1));ax.set_title('B   Shared assignments by site (hospitalizations)',loc='left',fontsize=11)
ax.legend(loc='upper center',bbox_to_anchor=(.5,-.26),ncol=5,frameon=False,fontsize=9)
fig.legend(handles=[Patch(facecolor=c,label=l) for c,l in zip(support_colors,['Room air','Low-flow O₂','HFNC','CPAP','NIPPV','IMV','Other','Unknown'])],loc='lower center',bbox_to_anchor=(.52,.105),ncol=8,frameon=False,fontsize=9)
fig.legend(handles=[Patch(facecolor=c,label=l) for c,l in zip(['#B8D9B0','#C85142','#BBBBBB'],['Evidence negative','Evidence positive','Unknown evidence'])],loc='lower center',bbox_to_anchor=(.52,.073),ncol=3,frameon=False,fontsize=9)
save(fig,'figure_2_smoothed_joint_phenotypes','Display-only interpolation between six-hour template centers; equal site weighting. Original template values and assignments are unchanged.\nCommon convex smoothing preserves nonnegative proportions summing to100%; it does not add hourly observations. I is an information category.')

# Read the newly released site-level curves, validating their manifests first.
cif=[]
for s in SITES:
 folder=OUT/'sites'/s/'share'
 if not folder.exists():folder=OUT/'sites'/s/'competing_events/share'
 mp=folder/'curve_manifest.json';inputs.add(mp);m=json.loads(mp.read_text())
 assert m['status']=='passed' and m['site']==s and m['consensus_id']==json.loads((FED/'coordinator/full/frozen_manifest.json').read_text())['consensus_id']
 for k,v in m['files'].items():assert hashlib.md5((folder/k).read_bytes()).hexdigest()==v
 cif.extend(read(folder/'cumulative_incidence_daily.csv'))
fig,axs=plt.subplots(3,3,figsize=(15,11),sharex=True,sharey=True)
fig.subplots_adjust(left=.09,right=.98,top=.84,bottom=.18,wspace=.18,hspace=.25)
fig.suptitle('Time to death, hospice or live discharge by phenotype',x=.035,ha='left',fontsize=20,y=.97)
fig.text(.035,.92,'Unadjusted competing-event cumulative incidence from the 72-hour ICU landmark',fontsize=12,color='#555555')
drawn=[]
for i,s in enumerate(SITES):
 for j,(event,title) in enumerate([('death','Hospital death'),('hospice','Hospice discharge'),('live_discharge','Other live discharge')]):
  ax=axs[i,j]
  for ph in PH+['I']:
   z=sorted([r for r in cif if r['site']==s and r['phenotype']==ph and r['event']==event],key=lambda r:n(r['day']))
   day=np.array([n(r['day']) for r in z]);v=np.array([n(r['estimate']) for r in z])
   assert np.all((v[np.isfinite(v)]>=0)&(v[np.isfinite(v)]<=1))
   # NaNs deliberately break the curve at suppressed days; never bridge gaps.
   ax.step(day,v,where='post',color=LC[ph],lw=2,ls='--' if ph=='I' else '-',label=ph)
   if np.any(np.isfinite(v)&(day>0)):drawn.append((s,event,ph))
  if i==0:ax.set_title(title,fontsize=13,pad=13)
  if j==0:ax.set_ylabel(s+'\nCumulative incidence',fontsize=11)
  if i==2:ax.set_xlabel('Days after 72-hour landmark')
  ax.set(xlim=(0,28),ylim=(0,1),xticks=[0,7,14,21,28],yticks=[0,.25,.5,.75,1]);ax.yaxis.set_major_formatter(PercentFormatter(1));ax.grid(alpha=.15)
  ax.text(.98,.04,f'Landmark N={count(s):,}',transform=ax.transAxes,ha='right',fontsize=9,color='#555555')
handles=[Line2D([],[],color=LC[p],lw=2,ls='--' if p=='I' else '-',label=p+' '+label) for p,label in zip(PH+['I'],PN+['Insufficient information'])]
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.52,.095),ncol=3,frameon=False,fontsize=10)
save(fig,'figure_3_competing_event_curves','Aalen–Johansen estimates at daily checkpoints; death, hospice and other live discharge are mutually exclusive competing events.\nEach site is displayed separately. Curves end at the last observed follow-up; exact aggregate counts are returned without disclosure masking.\nCurves are descriptive and unadjusted; the adjusted association forest plot remains a companion figure. No event timestamps are exported.')

# Absolute AUROC comparisons, with the same evaluation cases within each panel/site.
met=read(FED/'reports/tables/transport_metrics.csv');sofa=read(REV/'score_transport_metrics.csv');someta=read(REV/'score_model_metadata.csv')
ap=read(AP/'transport_metrics.csv');apcounts=read(AP/'cohort_counts.csv')
if (OUT/'score_model_metadata_exact.csv').exists():someta=read(OUT/'score_model_metadata_exact.csv')
panels=[('A   Full landmark cohort',met,None,
 [('baseline','Baseline'),('plus_phenotype','Baseline + phenotype'),('simple_summaries','Baseline + respiratory summaries'),('summaries_plus_phenotype','Baseline + summaries + phenotype')]),
 ('B   Complete modified-SOFA cohort',sofa,'measured_sofa',
 [('baseline','Baseline'),('phenotype','Baseline + phenotype'),('sofa','Baseline + modified SOFA'),('sofa_summaries','Baseline + SOFA + summaries'),('sofa_summaries_phenotype','Baseline + SOFA + summaries + phenotype')]),
 ('C   Complete adapted-APACHE cohort',ap,'complete_hb_proxy',
 [('baseline','Baseline'),('phenotype','Baseline + phenotype'),('baseline_adapted','Baseline + adapted APACHE'),('adapted_summaries','Baseline + APACHE + summaries'),('adapted_summaries_phenotype','Baseline + APACHE + summaries + phenotype')])]
fig,axs=plt.subplots(3,1,figsize=(14,13))
fig.subplots_adjust(left=.36,right=.97,top=.87,bottom=.13,hspace=.50)
fig.suptitle('Held-out AUROC comparisons',x=.035,ha='left',fontsize=22,y=.97)
fig.text(.035,.93,'28-day post-landmark hospital death or hospice; estimates and 95% patient-bootstrap intervals',fontsize=12,color='#555555')
for j,(title,data,cohort,models) in enumerate(panels):
 ax=axs[j]
 for i,s in enumerate(SITES):
  for k,(model,label) in enumerate(models):
   rr=[r for r in data if r['site']==s and r['model']==model and r['metric']=='auc' and (cohort is None or (r['cohort']==cohort and r['endpoint']=='death_or_hospice'))];assert len(rr)==1
   r=rr[0];v,lo,hi=map(n,[r['estimate'],r['conf_low'],r['conf_high']]);assert .4 <lo<=v<=hi<1
   ax.errorbar(v,k+(i-1)*.19,xerr=[[v-lo],[hi-v]],fmt='o',color=SC[s],capsize=2,ms=5,lw=1.5)
 ax.set(yticks=range(len(models)),yticklabels=[label for _,label in models],xlim=(.4,.95),xticks=np.arange(.4,1,.1),xlabel='Held-out AUROC')
 ax.invert_yaxis();ax.grid(axis='x',alpha=.18)
 if j==0:sz='; '.join(s+' '+str(count(s)) for s in SITES)
 elif j==1:
  sz='; '.join(s+' '+next((r.get('n') or (r['n_lower'] if r['n_lower']==r['n_upper'] else r['n_lower']+'–'+r['n_upper'])) for r in someta if r['site']==s and r['cohort']=='measured_sofa' and r['endpoint']=='death_or_hospice' and r['model']=='sofa_summaries_phenotype') for s in SITES)
 else:sz='; '.join(s+' '+next(r['n'] for r in apcounts if r['site']==s and r['cohort']=='complete_hb_proxy') for s in SITES)
 ax.set_title(title+'\nEvaluation N: '+sz,loc='left',fontsize=12,pad=16)
fig.legend(handles=[Line2D([],[],color=SC[s],marker='o',ls='-',label=s) for s in SITES],loc='upper right',bbox_to_anchor=(.97,.94),ncol=3,frameon=False)
save(fig,'figure_4_absolute_auroc','Same evaluation cases across models within each site/panel; cohorts differ between panels. Baseline includes demographic/clinical covariates.\nEach target site is excluded from phenotype and predictor learning. Bootstrap intervals condition on frozen learning; paired changes are reported separately.\nSOFA is modified; APACHE is a CLIF-adapted APS+age proxy. Absolute AUROCs cannot be compared across unmatched cohorts.')

# A direct same-patient score comparison: no mismatched complete-case populations.
models=[('phenotype','Baseline + phenotype'),('baseline_sofa','Baseline + modified SOFA'),('baseline_adapted','Baseline + adapted APACHE'),('sofa_summaries_phenotype','Baseline + SOFA + summaries + phenotype'),('adapted_summaries_phenotype','Baseline + APACHE + summaries + phenotype')]
fig,ax=plt.subplots(figsize=(13,7.5));fig.subplots_adjust(left=.39,right=.97,top=.72,bottom=.18)
fig.suptitle('Direct score comparison in the same patients',x=.035,ha='left',fontsize=21,y=.96)
fig.text(.035,.90,'Held-out AUROC for 28-day death or hospice; complete adapted APACHE and modified SOFA',fontsize=12)
sz='; '.join(s+' '+next(r['n'] for r in apcounts if r['site']==s and r['cohort']=='matched_complete_SOFA') for s in SITES)
ax.set_title('Matched evaluation N: '+sz,loc='left',pad=17,fontsize=12)
for i,s in enumerate(SITES):
 for j,(model,_) in enumerate(models):
  r=next(r for r in ap if r['site']==s and r['cohort']=='matched_complete_SOFA' and r['endpoint']=='death_or_hospice' and r['metric']=='auc' and r['model']==model)
  v,lo,hi=map(n,[r['estimate'],r['conf_low'],r['conf_high']]);assert .3<lo<=v<=hi<1
  ax.errorbar(v,j+(i-1)*.2,xerr=[[v-lo],[hi-v]],fmt='o',color=SC[s],capsize=2,ms=6,lw=1.6)
ax.set(yticks=range(len(models)),yticklabels=[label for _,label in models],xlim=(.3,1),xlabel='Held-out AUROC (95% interval)');ax.invert_yaxis();ax.grid(axis='x',alpha=.18)
fig.legend(handles=[Line2D([],[],color=SC[s],marker='o',ls='-',label=s) for s in SITES],loc='upper right',bbox_to_anchor=(.97,.85),ncol=3,frameon=False)
save(fig,'supplement_matched_score_auroc','All models are evaluated in identical patients within each site. Small complete-case cohorts produce wide uncertainty.\nThese modified score comparisons do not establish superiority over validated standard SOFA or APACHE II/IV; paired differences remain the inferential comparison.')

(OUT/'figure_manifest.json').write_text(json.dumps({'status':'passed','aggregate_inputs_only':True,
 'input_md5':{str(p.relative_to(ROOT)):hashlib.md5(p.read_bytes()).hexdigest() for p in sorted(inputs)},
 'smoothing':'Display-only convex cubic interpolation; all original centers unchanged; nonnegative sum-to-one verified',
 'time_event':'Daily Aalen-Johansen; competing events; withheld cells retained as gaps; no extrapolation beyond follow-up',
 'clinical_model_refits':False,'curve_series_with_observed_postzero_points':len(drawn),
 'auc':'All intervals fit plotting limits; model/site matches unique; same-case panels; separate matched score plot'},indent=2))
print(OUT)
