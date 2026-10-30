"""Render survival functions, AUC(t), and horizon-specific ROC curves from aggregates."""
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
ROOT=Path(__file__).resolve().parents[2];OUT=Path(os.environ.get('CLIF_SURVIVAL_OUTPUT',ROOT/'output/run_survival_tdauc_20261030_v2')).resolve()
SITES=['UCMC','NU','RUSH'];SC=dict(zip(SITES,['#247BA0','#D76345','#33966A']))
PH=['P1','P2','P3','P4','I'];PC=dict(zip(PH,['#648BA7','#347DB1','#4FA343','#E65D0A','#777777']))
inputs=set()
def read(p):
 inputs.add(p)
 with p.open() as f:return list(csv.DictReader(f))
def n(v):return float(v) if v else np.nan
surv=[];auc=[];roc=[]
POOL=Path(os.environ.get('CLIF_POOLED_SURVIVAL_OUTPUT',ROOT/'output/run_nearest_assignment_20261030/pooled_discharge_alive')).resolve()/'reports'
DISCHARGE=read(POOL/'site_discharge_incidence.csv')
inputs.add(POOL/'pooling_manifest.json');pool_manifest=json.loads((POOL/'pooling_manifest.json').read_text());assert pool_manifest['status']=='passed'
FIGOUT=Path(os.environ.get('CLIF_SURVIVAL_FIGURE_OUTPUT',OUT)).resolve();FIGOUT.mkdir(parents=True,exist_ok=True)
for site in SITES:
 folder=OUT/'sites'/site/'share'
 if not folder.exists():folder=OUT/'sites'/site/'survival/share'
 p=folder/'manifest.json';inputs.add(p);m=json.loads(p.read_text())
 assert m['status']=='passed' and m['site']==site
 assert m['primary_consensus_id']==pool_manifest['consensus_id']
 for k,v in m['files'].items():assert hashlib.md5((folder/k).read_bytes()).hexdigest()==v
 surv.extend(read(folder/'survival_curves.csv'));auc.extend(read(folder/'tdauc.csv'));roc.extend(read(folder/'tdroc.csv'))
discharge=[]
for r in DISCHARGE:
 q=dict(r);q['clock']='icu_admission' if r['curve'].startswith('icu_admission_') else 'landmark72h'
 q['group']=r['curve'].replace('icu_admission_','').replace('landmark72h_','')
 q['day']=r['time'];q['survival']=r['estimate'];discharge.append(q)
den=read(Path(os.environ.get('CLIF_DESCRIPTOR_REPORTS',ROOT/'output/run_descriptives_20261030_v2/reports')).resolve()/'denominators_by_site.csv')
N=lambda site,cohort:int(next(r['n'] for r in den if r['site']==site and r['cohort']==cohort and r['group']=='Overall'))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','pdf.fonttype':42})
def save(fig,name,note):
 fig.text(.035,.025,note,fontsize=9,color='#444444',va='bottom')
 for ext in ['png','pdf']:fig.savefig(FIGOUT/f'{name}.{ext}',dpi=180,facecolor='white')
 plt.close(fig)
def rows(data,**kwargs):return [r for r in data if all(str(r.get(k))==str(v) for k,v in kwargs.items())]
def km_series(ax,z,color,label,ci=False,ls='-'):
 z=sorted(z,key=lambda r:n(r['day']));x=np.array([n(r['day']) for r in z]);y=np.array([n(r['survival']) for r in z])
 ax.step(x,y,where='post',color=color,label=label,lw=2,ls=ls)
 ax.plot(x,y,'.',color=color,ms=2)
 if ci:ax.fill_between(x,[n(r['conf_low']) for r in z],[n(r['conf_high']) for r in z],color=color,alpha=.09,step='post')
 ax.set(xlim=(0,28),ylim=(0,1.02),xticks=[0,7,14,21,28]);ax.yaxis.set_major_formatter(PercentFormatter(1));ax.grid(alpha=.15)
def risks(ax,groups,clock,site=None):
 ax.axis('off');table=[]
 for g in groups:
  z=rows(surv,clock=clock,endpoint='death',group=g,**({'site':site} if site else {}))
  if site is None:z=rows(surv,clock=clock,endpoint='death',group='Overall',site=g)
  v=[]
  for t in [0,7,14,21,28]:
   q=[r for r in z if n(r['day'])==t]
   zero_group=site is not None and any(r['site']==site and r['cohort']=='primary_phenotypes' and r['group']==g and n(r['n'])==0 for r in den)
   v.append((q[0].get('n_at_risk') or (q[0]['at_risk_lower'] if q[0]['at_risk_lower']==q[0]['at_risk_upper'] else q[0]['at_risk_lower']+'–'+q[0]['at_risk_upper'])) if q else ('0' if zero_group or (z and t>max(n(r['day']) for r in z)) else '—'))
  table.append(v)
 tab=ax.table(cellText=table,rowLabels=groups,colLabels=['0','7','14','21','28'],cellLoc='center',loc='center')
 tab.auto_set_font_size(False);tab.set_fontsize(8);tab.scale(1,1.35)
 for cell in tab.get_celld().values():cell.set_edgecolor('#EEEEEE')
 ax.set_title('Number at risk',fontsize=9,pad=7)

# All eligible hospitalizations enter at ICU t0, irrespective of 72h survival.
fig=plt.figure(figsize=(14,9));gs=fig.add_gridspec(2,2,height_ratios=[3,1.1],left=.09,right=.97,bottom=.19,top=.82,hspace=.45,wspace=.20)
fig.suptitle('Whole-cohort survival from ICU admission',x=.035,ha='left',fontsize=22,y=.96)
fig.text(.035,.90,'All 4,533 eligible hospitalizations enter at time zero; no requirement to survive 72 hours',fontsize=12)
for j,(endpoint,title) in enumerate([('death','A   Hospital-death Kaplan–Meier function'),('live_discharge','B   Discharge alive without hospice')]):
 ax=fig.add_subplot(gs[0,j])
 for s in SITES:km_series(ax,rows(discharge if endpoint=='live_discharge' else surv,site=s,clock='icu_admission',group='Overall',endpoint=endpoint),SC[s],f'{s} (N={N(s,"eligible_cohort"):,})',endpoint=='death')
 ax.set_title(title,fontsize=12,pad=14);ax.set_xlabel('Days after ICU admission');ax.set_ylabel('Cumulative incidence' if endpoint=='live_discharge' else 'Kaplan–Meier estimate')
 ax.legend(frameon=False,fontsize=9,loc='lower left')
 risks(fig.add_subplot(gs[1,j]),SITES,'icu_admission')
save(fig,'whole_cohort_survival_t0','A: hospital-death KM, other exits censored. B: Aalen–Johansen live non-hospice discharge probability; death/hospice compete.\nUnknown terminal dispositions are censored at observed discharge, an explicit missing-outcome assumption. No post-discharge survival is inferred.\nDeath-KM bands are pointwise patient-clustered 95% intervals; discharge panels show point estimates. All site denominators are retained; this is separate from the 72h phenotype analysis.')

fig=plt.figure(figsize=(15,12));gs=fig.add_gridspec(3,3,height_ratios=[1,1,.60],left=.09,right=.98,top=.83,bottom=.16,hspace=.4,wspace=.2)
fig.suptitle('Survival by the frozen 72-hour phenotype',x=.035,ha='left',fontsize=21,y=.97)
fig.text(.035,.92,'Top: hospital-death KM, starts at 100%. Bottom: live-discharge incidence, starts at 0%. Clock: 72h landmark.',fontsize=12)
for i,endpoint in enumerate(['death','live_discharge']):
 for j,s in enumerate(SITES):
  ax=fig.add_subplot(gs[i,j])
  for ph in PH:km_series(ax,rows(discharge if endpoint=='live_discharge' else surv,site=s,clock='landmark72h',group=ph,endpoint=endpoint),PC[ph],ph,ls='--' if ph=='I' else '-')
  if i==0:ax.set_title(f'{s} (landmark N={N(s,"landmark_cohort"):,})',fontsize=12,pad=14)
  if j==0:ax.set_ylabel('Hospital-death KM' if endpoint=='death' else 'Live non-hospice discharge\nCumulative incidence')
  if i==1:ax.set_xlabel('Days after 72h landmark')
for j,s in enumerate(SITES):risks(fig.add_subplot(gs[2,j]),PH,'landmark72h',s)
fig.legend(handles=[Line2D([],[],color=PC[p],lw=2,ls='--' if p=='I' else '-',label=l) for p,l in zip(PH,['P1 Room-air predominant','P2 Low-flow oxygen','P3 HFNC','P4 IMV','I Insufficient information'])],loc='lower center',bbox_to_anchor=(.54,.08),ncol=5,frameon=False,fontsize=9)
save(fig,'phenotype_survival_landmark','Top: cause-specific death KM, other exits censored. Bottom: Aalen–Johansen live non-hospice discharge; death/hospice compete.\nFuture 72h phenotype labels are not assigned at ICU admission. No postdischarge survival is inferred.\nDeath KM uses daily displays; discharge uses exact event durations. Curves end when observed follow-up ends. P1–P4 are shared phenotypes; I is an information category.')

COLORS=['#666666','#357BB5','#32916B','#D96C32','#9D5A98']
def models_for(cohort,endpoint):
 if cohort=='full_landmark':return ([('baseline','Baseline'),('plus_phenotype','+ phenotype'),('simple_summaries','+ summaries'),('summaries_plus_phenotype','+ summaries + phenotype')] if endpoint=='death_or_hospice' else [('baseline','Baseline'),('phenotype','+ phenotype'),('summaries','+ summaries'),('summaries_phenotype','+ summaries + phenotype')])
 if cohort=='measured_sofa':return [('baseline','Baseline'),('phenotype','+ phenotype'),('sofa','+ SOFA'),('sofa_summaries','+ SOFA + summaries'),('sofa_summaries_phenotype','+ SOFA + summaries + phenotype')]
 if cohort=='complete_hb_proxy':return [('baseline','Baseline'),('phenotype','+ phenotype'),('baseline_adapted','+ adapted APACHE'),('adapted_summaries','+ APACHE + summaries'),('adapted_summaries_phenotype','+ APACHE + summaries + phenotype')]
 return [('phenotype','Baseline + phenotype'),('baseline_sofa','Baseline + SOFA'),('baseline_adapted','Baseline + adapted APACHE'),('sofa_summaries_phenotype','SOFA + summaries + phenotype'),('adapted_summaries_phenotype','APACHE + summaries + phenotype')]
def auc_panel(ax,cohort,endpoint,site,control,ci=True):
 for (model,label),color in zip(models_for(cohort,endpoint),COLORS):
  z=sorted(rows(auc,site=site,cohort=cohort,endpoint=endpoint,control_definition=control,model=model),key=lambda r:n(r['day']))
  assert len(z)==28
  x=[n(r['day']) for r in z];y=[n(r['estimate']) for r in z]
  ax.plot(x,y,color=color,label=label,lw=1.9)
  if ci:ax.fill_between(x,[n(r['conf_low']) for r in z],[n(r['conf_high']) for r in z],color=color,alpha=.065)
 ax.axhline(.5,color='#BBBBBB',ls='--',lw=1);ax.set(xlim=(1,28),ylim=(0,1),xticks=[1,7,14,21,28],yticks=[0,.25,.5,.75,1]);ax.grid(alpha=.13)

for control in [2,1]:
 fig,axs=plt.subplots(2,3,figsize=(15,10),sharex=True,sharey=True);fig.subplots_adjust(left=.075,right=.98,top=.82,bottom=.20,hspace=.28,wspace=.18)
 fig.suptitle('Time-dependent discrimination: AUC(t)',x=.035,ha='left',fontsize=22,y=.97)
 desc='Controls include other hospital exits and those still hospitalized' if control==2 else 'Controls are patients still free of any terminal hospital event'
 fig.text(.035,.92,'Full landmark cohort | '+desc,fontsize=12)
 for i,endpoint in enumerate(['death','death_or_hospice']):
  for j,s in enumerate(SITES):
   ax=axs[i,j];auc_panel(ax,'full_landmark',endpoint,s,control)
   if i==0:ax.set_title(f'{s} (N={N(s,"landmark_cohort"):,})',fontsize=12,pad=15)
   if j==0:ax.set_ylabel(('Hospital death' if endpoint=='death' else 'Death or hospice')+'\nAUC(t)')
   if i==1:ax.set_xlabel('Days after 72h landmark')
 fig.legend(handles=[Line2D([],[],color=c,lw=2,label=l) for c,l in zip(COLORS,['Baseline','Baseline + phenotype','Baseline + respiratory summaries','Baseline + summaries + phenotype'])],loc='lower center',bbox_to_anchor=(.54,.105),ncol=2,frameon=False)
 save(fig,'tdauc_full_cohort' if control==2 else 'supplement_tdauc_event_free_controls','Cumulative/dynamic timeROC estimates; shaded bands are pointwise 95% patient-cluster bootstrap intervals (200 resamples).\nFixed frozen 28-day risk markers are evaluated at successive horizons; models are not refitted at each day or at the held-out site.\nAUC is undefined when there are no cases or controls. AUC(t) at 28d with control definition2 reproduces the prior primary AUROC.')

fig,axs=plt.subplots(2,3,figsize=(15,11),sharex=True,sharey=True);fig.subplots_adjust(left=.08,right=.98,top=.82,bottom=.22,hspace=.43,wspace=.18)
fig.suptitle('Time-dependent AUC: severity-score comparisons',x=.035,ha='left',fontsize=21,y=.97)
fig.text(.035,.92,'28-day post-landmark hospital death or hospice | competing-event-inclusive controls',fontsize=12)
for i,cohort in enumerate(['measured_sofa','complete_hb_proxy']):
 for j,s in enumerate(SITES):
  ax=axs[i,j];auc_panel(ax,cohort,'death_or_hospice',s,2)
  ax.set_title(s+' · '+('Complete modified SOFA' if i==0 else 'Complete adapted APACHE'),fontsize=11,pad=13)
  if j==0:ax.set_ylabel('AUC(t)')
  if i==1:ax.set_xlabel('Days after 72h landmark')
 for k,(model,label) in enumerate(models_for(cohort,'death_or_hospice')):
  fig.text(.11+k*.175,.535 if i==0 else .15,label,fontsize=8,color=COLORS[k],ha='left')
save(fig,'tdauc_severity_comparisons','Same evaluation cases across models within each panel. SOFA and APACHE panels use different complete-case cohorts.\nScores are modified CLIF proxies; static frozen predictors are assessed over time. Bands are pointwise patient-bootstrap 95% intervals.\nAll detailed estimates, both endpoints and both competing-control definitions are returned as aggregates.')

fig,axs=plt.subplots(1,3,figsize=(15,6.5),sharex=True,sharey=True);fig.subplots_adjust(left=.075,right=.98,top=.70,bottom=.28,wspace=.18)
fig.suptitle('Time-dependent AUC: scores evaluated in the same patients',x=.035,ha='left',fontsize=20,y=.96)
fig.text(.035,.88,'Death or hospice | complete modified SOFA and adapted APACHE | UCMC 96, NU 198, Rush 47',fontsize=12)
for j,s in enumerate(SITES):auc_panel(axs[j],'matched_complete_SOFA','death_or_hospice',s,2);axs[j].set_title(s);axs[j].set_xlabel('Days after 72h landmark')
axs[0].set_ylabel('AUC(t)')
fig.legend(handles=[Line2D([],[],color=c,lw=2,label=l) for c,(_,l) in zip(COLORS,models_for('matched_complete_SOFA','death_or_hospice'))],loc='lower center',bbox_to_anchor=(.54,.13),ncol=3,frameon=False,fontsize=9)
save(fig,'tdauc_matched_scores','Models share identical evaluation patients within each site; these small cohorts produce wide uncertainty at early horizons.\nCumulative/dynamic AUC(t), control definition2; fixed frozen predictors, pointwise patient-bootstrap 95% intervals. No prognostic-superiority claim.')

fig,axs=plt.subplots(3,3,figsize=(14,13),sharex=True,sharey=True);fig.subplots_adjust(left=.085,right=.98,top=.86,bottom=.14,hspace=.28,wspace=.20)
fig.suptitle('Time-dependent ROC curves at 7, 14 and 28 days',x=.035,ha='left',fontsize=21,y=.97)
fig.text(.035,.92,'Held-out prediction of post-landmark hospital death or hospice; full landmark cohort',fontsize=12)
for i,t in enumerate([7,14,28]):
 for j,s in enumerate(SITES):
  ax=axs[i,j]
  for (model,label),color in zip(models_for('full_landmark','death_or_hospice'),COLORS):
   z=sorted(rows(roc,site=s,cohort='full_landmark',endpoint='death_or_hospice',day=t,control_definition=2,model=model),key=lambda r:n(r['fpr']))
   ax.plot([0]+[n(r['fpr']) for r in z],[0]+[n(r['sensitivity']) for r in z],color=color,lw=1.7,label=label+f' ({n(z[0]["auc"]):.3f})')
  ax.plot([0,1],[0,1],ls='--',color='#BBBBBB',lw=1);ax.set(xlim=(0,1),ylim=(0,1),xticks=[0,.5,1],yticks=[0,.5,1]);ax.grid(alpha=.12);ax.set_aspect('equal')
  if i==0:ax.set_title(s,fontsize=13,pad=12)
  if j==0:ax.set_ylabel(f'Day {t}\nSensitivity')
  if i==2:ax.set_xlabel('1 − specificity')
  ax.legend(frameon=False,loc='lower right',fontsize=7.4,title='Model (time-dependent AUC)',title_fontsize=8)
save(fig,'tdroc_7_14_28','Cumulative/dynamic ROC; controls include competing hospital exits. Sensitivity is sampled on a common fixed false-positive-rate grid.\nLegend AUCs are the exact timeROC estimates, not numerical areas under the sampled display. No individual predictions or thresholds leave a site.\nThese are horizon-specific ROC curves; the separate AUC(t) figure plots discrimination continuously across follow-up days.')

for name,data in [('survival_curves',surv),('tdauc',auc),('tdroc',roc)]:
 with (FIGOUT/(name+'_by_site.csv')).open('w') as f:
  writer=csv.DictWriter(f,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
(FIGOUT/'plot_manifest.json').write_text(json.dumps({'status':'passed','aggregate_inputs_only':True,'input_md5':{str(p.relative_to(ROOT)):hashlib.md5(p.read_bytes()).hexdigest() for p in sorted(inputs)},'whole_cohort_n':sum(N(s,'eligible_cohort') for s in SITES),'landmark_n':sum(N(s,'landmark_cohort') for s in SITES),'clinical_model_refits':False,'clocks_kept_separate':True,'survival_is_cause_specific_KM_not_postdischarge_allcause_survival':True,'auc_package':'timeROC','tdauc_controls_main':2,'tdauc_controls_supplement':1},indent=2))
print(OUT)
