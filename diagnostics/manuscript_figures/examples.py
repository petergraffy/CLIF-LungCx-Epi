"""Main-paper layout examples, drawn exclusively from versioned aggregate results."""
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/clif_manuscript_mpl')
import csv, json, hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'output/run_manuscript_figure_examples_20261030'
OUT.mkdir(exist_ok=True)
FED = ROOT / 'output/run_federated_rehearsal_20261026'
REV = ROOT / 'output/run_revision_gaps_20261026/reports/tables'
AP = ROOT / 'output/run_APACHE_adapted_20261030_v2/reports/tables'
DESC = ROOT / 'output/run_descriptives_20261030_v2/reports'
inputs = set()
def read(p):
    inputs.add(p)
    with p.open() as f: return list(csv.DictReader(f))
def num(x): return float(x) if x else np.nan
SITES = ['UCMC', 'NU', 'RUSH']
SC = dict(zip(SITES, ['#247BA0', '#D76345', '#33966A']))
PH = ['P1', 'P2', 'P3', 'P4']
PN = ['Room-air predominant', 'Low-flow oxygen predominant', 'HFNC predominant', 'IMV predominant']
PC = ['#B8D4E8', '#347DB1', '#4FA343', '#E65D0A', '#A5A5A5']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                     'axes.spines.right':False,'axes.titleweight':'bold','pdf.fonttype':42})
def finish(fig, name, note):
    fig.text(.035, .025, note, fontsize=9, color='#444444', va='bottom')
    for ext in ['png','pdf']: fig.savefig(OUT / f'{name}.{ext}', dpi=180, facecolor='white')
    plt.close(fig)
def box(ax, x,y,w,h,text,color='#EDF3F7',fs=11):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.012',
                              facecolor=color,edgecolor='#8CA0AD',linewidth=1))
    ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=fs)
def arrow(ax,x1,y1,x2,y2):
    ax.annotate('',xy=(x2,y2),xytext=(x1,y1),arrowprops={'arrowstyle':'->','color':'#556677','lw':1.5})

# Figure 1: exact released hospitalization denominators; residual not mislabeled as deaths.
den = read(DESC/'denominators_by_site.csv')
counts = {s:{c:int(next(r['n'] for r in den if r['site']==s and r['cohort']==c and r['group']=='Overall'))
             for c in ['eligible_cohort','landmark_cohort']} for s in SITES}
assert sum(x['landmark_cohort'] for x in counts.values())==3310
fig, axes = plt.subplots(1,2,figsize=(15,8.2),gridspec_kw={'width_ratios':[1.12,1]})
fig.subplots_adjust(left=.035,right=.97,top=.84,bottom=.13,wspace=.16)
fig.suptitle('Figure 1  |  Study cohort and federated analysis',x=.035,ha='left',fontsize=20,y=.96)
fig.text(.035,.90,'Layout example using the current UCMC, NU and Rush results',fontsize=12,color='#555555')
for ax in axes: ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
ax=axes[0];ax.set_title('A   Cohort construction',loc='left',pad=15)
box(ax,.02,.77,.94,.16,'Eligible lung-cancer ICU cohort\n4,533 hospitalizations',fs=14)
for i,s in enumerate(SITES):
    x=.025+i*.325
    arrow(ax,x+.145,.76,x+.145,.68)
    box(ax,x,.48,.29,.19,f'{s}\nEligible: {counts[s]["eligible_cohort"]:,}\nLandmark: {counts[s]["landmark_cohort"]:,}',fs=11)
    arrow(ax,x+.145,.47,.49,.36)
box(ax,.02,.17,.94,.18,'Primary analysis: ICU +72-hour landmark\n3,310 hospitalizations; unknown states retained',color='#E8F2E9',fs=13)
ax.text(.49,.06,'1,223 outside the primary landmark cohort\nRetain in cohort-flow, early-event and outcome-quality audits',ha='center',va='center',fontsize=10)
ax=axes[1];ax.set_title('B   Outcome-free learning, then outcome analysis',loc='left',pad=15)
texts=[('LOCAL','Joint support–ARF sequences, hours 0–72\nLocal clusters and aggregate fingerprints'),
       ('CENTRAL','Learn four shared phenotype templates\nFreeze definitions before outcome analyses'),
       ('LOCAL','Apply frozen templates; fit outcome models\nExport associations and prediction metrics'),
       ('CENTRAL','Pool site associations\nEvaluate each site held out of learning')]
for i,(tag,text) in enumerate(texts):
    y=.76-i*.225
    box(ax,.02,y,.94,.16,tag+'\n'+text,fs=10.5)
    if i<3:arrow(ax,.49,y-.005,.49,y-.06)
finish(fig,'figure_1_cohort_federation',
       'Hospitalizations are the analysis unit. Site-distinct patients may have repeated admissions.\nOnly aggregates move between sites; this three-site run rehearses federation locally. Detailed exclusion reasons belong in the flow supplement.')

# Figure 2: the same exported prototype marginals as the existing support/ARF plots.
pro = read(FED/'coordinator/full/prototypes.csv')
states = ['ROOM_AIR','LOW_O2','HFNC','CPAP','NIPPV','IMV','OTHER','UNKNOWN']
support_colors=['#B8D4E8','#347DB1','#4FA343','#AA59AE','#7142A2','#E65D0A','#D6B565','#A5A5A5']
fig=plt.figure(figsize=(16,11))
gs=fig.add_gridspec(4,4,height_ratios=[1.5,1,1,1.15],left=.085,right=.98,top=.85,bottom=.19,hspace=.53,wspace=.22)
fig.suptitle('Figure 2  |  Shared respiratory support–ARF phenotypes',x=.035,ha='left',fontsize=20,y=.975)
fig.text(.035,.93,'A   Consensus templates: aligned support and physiological evidence over the first 72 ICU hours',fontsize=12)
for row,(domain,order,colors,label) in enumerate([
    ('support',states,support_colors,'Respiratory\nsupport'),
    ('hypox',['H0','H1','HU'],['#B8D9B0','#C85142','#BBBBBB'],'Hypoxemia\nevidence'),
    ('hyper',['C0','C1','CU'],['#B8D9B0','#C85142','#BBBBBB'],'Hypercapnia/\nacidemia evidence')]):
    for j,ph in enumerate(PH):
        ax=fig.add_subplot(gs[row,j]);bottom=np.zeros(12)
        for state,color in zip(order,colors):
            z=sorted([r for r in pro if r['phenotype']==ph and r['domain']==domain and r['state']==state],key=lambda r:int(r['bin']))
            assert len(z)==12
            v=np.array([num(r['probability']) for r in z]);ax.bar(np.arange(12)*6,v,width=6,align='edge',bottom=bottom,color=color,linewidth=0);bottom+=v
        assert np.allclose(bottom,1)
        ax.set(xlim=(0,72),ylim=(0,1),xticks=[0,24,48,72],yticks=[0,.5,1])
        ax.yaxis.set_major_formatter(PercentFormatter(1));ax.tick_params(labelsize=9)
        if j==0:ax.set_ylabel(label,labelpad=18,fontsize=10)
        else:ax.set_yticklabels([])
        if row==0:ax.set_title(ph+'\n'+PN[j],fontsize=11,pad=14)
        if row==2:ax.set_xlabel('Hours after ICU arrival',fontsize=10)
        else:ax.set_xticklabels([])
ax=fig.add_subplot(gs[3,:]);left=np.zeros(3)
for ph,c in zip(PH+['I'],PC):
    ns=np.array([int(next(r['n'] for r in den if r['site']==s and r['cohort']=='primary_phenotypes' and r['group']==ph)) for s in SITES])
    totals=np.array([counts[s]['landmark_cohort'] for s in SITES]);v=ns/totals
    ax.barh(SITES,v,left=left,color=c,label=ph if ph!='I' else 'I: insufficient information',height=.62)
    for i,(n,w) in enumerate(zip(ns,v)):
        if w>.055:ax.text(left[i]+w/2,i,f'{n:,}',ha='center',va='center',fontsize=10,color='white' if ph in ['P2','P4'] else '#222222')
    left+=v
assert np.allclose(left,1)
ax.invert_yaxis();ax.set_xlim(0,1);ax.xaxis.set_major_formatter(PercentFormatter(1));ax.set_title('B   Shared assignments by site (numbers are hospitalizations)',loc='left',fontsize=11)
ax.legend(loc='upper center',bbox_to_anchor=(.5,-.26),ncol=5,frameon=False,fontsize=9)
support_handles=[Patch(facecolor=c,label=l) for c,l in zip(support_colors,['Room air','Low-flow O₂','HFNC','CPAP','NIPPV','IMV','Other','Unknown'])]
fig.legend(handles=support_handles,loc='lower center',bbox_to_anchor=(.52,.105),ncol=8,frameon=False,fontsize=9)
fig.legend(handles=[Patch(facecolor=c,label=l) for c,l in zip(['#B8D9B0','#C85142','#BBBBBB'],['Evidence negative','Evidence positive','Unknown evidence'])],loc='lower center',bbox_to_anchor=(.52,.073),ncol=3,frameon=False,fontsize=9)
finish(fig,'figure_2_joint_phenotypes',
       'Panel A: equal-site-weighted aggregate templates; six-hour summaries reconstructed from interval fingerprints, not individual trajectories.\nPanel B: observed frozen assignments. I is an information category, not a fifth clinical phenotype. Unknown evidence does not imply normal physiology.')

# Figure 3: actual released risk estimates and model contrasts; suppression stays missing.
risks=read(FED/'reports/tables/risks_28day.csv')
assoc=read(FED/'reports/tables/site_associations.csv')
pool=read(FED/'reports/tables/pooled_associations.csv')
fig,axs=plt.subplots(1,2,figsize=(15,9),gridspec_kw={'width_ratios':[1,1.1]})
fig.subplots_adjust(left=.19,right=.97,top=.81,bottom=.17,wspace=.42)
fig.suptitle('Figure 3  |  Outcomes after the 72-hour landmark',x=.035,ha='left',fontsize=20,y=.96)
fig.text(.035,.90,'Hospital death or hospice within 28 days after landmark',fontsize=12,color='#555555')
ax=axs[0]
for i,s in enumerate(SITES):
    for j,ph in enumerate(PH+['I']):
        r=next(r for r in risks if r['site']==s and r['phenotype']==ph and r['endpoint']=='death_or_hospice')
        if r['status']!='ok':continue
        v,lo,hi=map(num,[r['risk'],r['conf_low'],r['conf_high']])
        ax.errorbar(v,j+(i-1)*.18,xerr=[[v-lo],[hi-v]],fmt='o',color=SC[s],capsize=2,ms=6)
ax.set(yticks=range(5),yticklabels=[p+'  '+n for p,n in zip(PH,PN)]+['I  Insufficient information'],xlim=(0,.85),xlabel='Observed risk (95% interval)')
ax.invert_yaxis();ax.xaxis.set_major_formatter(PercentFormatter(1));ax.grid(axis='x',alpha=.18);ax.set_title('A   Absolute risks by site',loc='left',pad=17)
ax.legend(handles=[Line2D([],[],marker='o',ls='',color=SC[s],label=s) for s in SITES],frameon=False,loc='lower right')
ax=axs[1];yt=[];yl=[]
for j,ph in enumerate(PH[1:]):
    for i,s in enumerate(SITES):
        y=j*5+i;rs=[r for r in assoc if r['site']==s and r['phenotype']==ph and r['model']=='core' and r['endpoint']=='death_or_hospice' and r['status']=='ok']
        if rs:
            r=rs[0];est,se=num(r['estimate']),num(r['std_error']);v,lo,hi=np.exp([est,est-1.96*se,est+1.96*se]);ax.errorbar(v,y,xerr=[[v-lo],[hi-v]],fmt='o',color=SC[s],ms=5,capsize=2)
        yt.append(y);yl.append(ph+' · '+s)
    r=next(r for r in pool if r['phenotype']==ph and r['model']=='core' and r['endpoint']=='death_or_hospice')
    v,lo,hi=map(num,[r['hazard_ratio'],r['hr_low'],r['hr_high']]);y=j*5+3
    ax.errorbar(v,y,xerr=[[v-lo],[hi-v]],fmt='D',color='#202E3C',ms=6,capsize=3,lw=2)
    yt.append(y);yl.append(f'{ph} · Pooled (I² {num(r["i2_percent"]):.0f}%)')
ax.set_xscale('log');ax.set(yticks=yt,yticklabels=yl,xlim=(.5,10),xlabel='Adjusted hazard ratio (95% interval); reference P1',xticks=[.5,1,2,4,8])
ax.set_xticklabels(['0.5','1','2','4','8']);ax.invert_yaxis();ax.axvline(1,color='#888888',ls='--',lw=1);ax.grid(axis='x',alpha=.18);ax.set_title('B   Site and pooled associations',loc='left',pad=17)
finish(fig,'figure_3_outcomes',
       'Core adjustment: age, sex, race, admission year, noncancer Charlson, advanced-cancer proxy and admission code status.\nPatient-clustered site uncertainty; pooled REML / modified Hartung–Knapp intervals. Withheld risk partitions are omitted, not zero.\nI remains in the analysis with a separate baseline-hazard stratum; no I-versus-P1 hazard ratio is inferred. Associations are not causal.')

# Figure 4: paired increments, never comparing AUCs across unmatched score cohorts.
met=read(FED/'reports/tables/transport_metrics.csv')
sofa=read(REV/'score_paired_increments.csv');someta=read(REV/'score_model_metadata.csv')
ap=read(AP/'paired_increments.csv');apcounts=read(AP/'cohort_counts.csv')
fig,axs=plt.subplots(1,3,figsize=(15,7.8))
fig.subplots_adjust(left=.08,right=.98,top=.76,bottom=.27,wspace=.38)
fig.suptitle('Figure 4  |  Does phenotype add prognostic information?',x=.035,ha='left',fontsize=20,y=.96)
fig.text(.035,.90,'Each site held out of phenotype learning and predictor fitting; 28-day death or hospice',fontsize=12,color='#555555')
for j,ax in enumerate(axs):
    labels=[]
    for i,s in enumerate(SITES):
        if j==0:
            r=next(r for r in met if r['site']==s and r['model']=='phenotype_beyond_summaries' and r['metric']=='auc');size=f'n={counts[s]["landmark_cohort"]:,}'
        elif j==1:
            r=next(r for r in sofa if r['site']==s and r['cohort']=='measured_sofa' and r['endpoint']=='death_or_hospice' and r['comparison']=='phenotype_beyond_sofa' and r['metric']=='auc')
            md=next(r for r in someta if r['site']==s and r['cohort']=='measured_sofa' and r['endpoint']=='death_or_hospice' and r['model']=='sofa_summaries_phenotype');size=f'n={md["n_lower"]}–{md["n_upper"]}'
        else:
            r=next(r for r in ap if r['site']==s and r['cohort']=='complete_hb_proxy' and r['endpoint']=='death_or_hospice' and r['comparison']=='phenotype_beyond_adapted_and_summaries' and r['metric']=='auc')
            md=next(r for r in apcounts if r['site']==s and r['cohort']=='complete_hb_proxy');size=f'n={md["n"]}'
        v,lo,hi=map(num,[r['estimate'],r['conf_low'],r['conf_high']]);assert -.09 < lo <= v <= hi < .065
        ax.errorbar(v,i,xerr=[[v-lo],[hi-v]],fmt='o',color=SC[s],capsize=4,ms=8,lw=2)
        labels.append(s+'\n'+size)
    ax.set(yticks=range(3),yticklabels=labels,ylim=(2.55,-.55),xlabel='Change in held-out AUC (95% interval)')
    ax.axvline(0,color='#777777',ls='--',lw=1);ax.grid(axis='x',alpha=.18);ax.set_xlim(-.09,.065);ax.set_xticks([-.075,-.05,-.025,0,.025,.05]);ax.set_xticklabels(['−0.075','−0.050','−0.025','0','+0.025','+0.050'],fontsize=8)
    ax.set_title(['A   Beyond simple summaries\nFull landmark cohort','B   Beyond modified SOFA\n+ simple summaries','C   Beyond adapted APACHE proxy\n+ simple summaries'][j],fontsize=11,pad=22)
finish(fig,'figure_4_incremental_prediction',
       'Positive changes favor adding phenotype. Baseline covariates are included in each comparator. Paired bootstrap intervals condition on frozen learning.\nPanels use different evaluation cohorts; do not compare AUC levels between panels. SOFA N is reported in released count bands.\nModified SOFA and the CLIF-adapted APACHE APS+age proxy are not validated standard scores. Same-patient score comparisons and\nall-case missing-indicator sensitivities belong in the supplement. These results do not show consistent incremental discrimination.')

ledger={str(p.relative_to(ROOT)):hashlib.md5(p.read_bytes()).hexdigest() for p in sorted(inputs)}
(OUT/'figure_manifest.json').write_text(json.dumps({'status':'example_layouts','aggregate_inputs_only':True,'primary_consensus_id':'6245d53b0469ab7d1a5b','input_md5':ledger,'notes':['No clinical model refits','Site HR intervals reconstructed from released log-HR and clustered SE','Do not interpret templates as observed individual trajectories','Small partitions and count bands preserved']},indent=2))
print(OUT)
