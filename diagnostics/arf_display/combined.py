"""Four mutually exclusive documented-evidence categories; unknowns remain in source states."""
import argparse,csv,json,hashlib,os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/clif_manuscript_mpl')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
ROOT=Path(__file__).resolve().parents[2]
LABELS=['Hypoxemia evidence only','Hypercapnic-acidemic evidence only','Mixed positive evidence','No positive physiology documented']
COLORS=['#D85843','#5B73B4','#8C4D9A','#D3D8DA']
def category(h,c):
 if h not in ['H0','H1','HU'] or c not in ['C0','C1','CU']:raise ValueError('Unknown physiology state')
 if h=='H1' and c=='C1':return 2
 if h=='H1':return 0
 if c=='C1':return 1
 return 3

def main():
 p=argparse.ArgumentParser();p.add_argument('output',type=Path);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
 active=json.loads((ROOT/'config/active_analysis.json').read_text())
 source=ROOT/active['exact_root']/'reports/aggregate_collection/exact_counts/profiles/hourly_joint.csv'
 with source.open() as f:rows=list(csv.DictReader(f))
 groups=['Overall','P1','P2','P3','P4','I'];ncat={g:np.zeros((4,72)) for g in groups};coverage={g:np.zeros((2,72)) for g in groups}
 bysite={}
 for r in rows:
  ph=r['phenotype'];assert ph in groups[1:];hour=int(r['h']);assert 0<=hour<72;n=int(r['n'])
  h,c=r['state'].split('|')[1].split('_');k=category(h,c)
  bysite[(r['site'],ph,hour)]=bysite.get((r['site'],ph,hour),0)+n
  for g in [ph,'Overall']:
   ncat[g][k,hour]+=n
   coverage[g][0,hour]+=n*(h!='HU' and c!='CU')
   coverage[g][1,hour]+=n*(h!='HU' or c!='CU')
 den={g:ncat[g].sum(axis=0) for g in groups}
 assert np.all(den['Overall']==3310) and all(np.all(den[g]==den[g][0]) for g in groups)
 fractions={g:ncat[g]/den[g] for g in groups};assert all(np.allclose(z.sum(axis=0),1) for z in fractions.values())
 # Denominator includes every hour: unknown is grouped only for this display,
 # never overwritten in the underlying hourly hypox/hyper states.
 export=[]
 for g in groups:
  for hour in range(72):
   for k,label in enumerate(LABELS):export.append(dict(group=g,hour=hour,category=label,n=int(ncat[g][k,hour]),denominator=int(den[g][hour]),proportion=fractions[g][k,hour],both_channels_measured_n=int(coverage[g][0,hour]),at_least_one_channel_measured_n=int(coverage[g][1,hour])))
 with (a.output/'hourly_documented_evidence.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(export[0]));w.writeheader();w.writerows(export)
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
 titles={'Overall':'Overall','P1':'P1: Room-air predominant','P2':'P2: Low-flow oxygen','P3':'P3: HFNC','P4':'P4: IMV','I':'I: Insufficient support data'}
 def plot(gs,stem):
  fig=plt.figure(figsize=(max(12.5,3.3*len(gs)),8));grid=fig.add_gridspec(2,len(gs),height_ratios=[3,1],left=.065,right=.98,top=.76,bottom=.32,hspace=.24,wspace=.20)
  fig.suptitle('Documented ARF physiology during the first 72 ICU hours',x=.035,ha='left',fontsize=19,y=.97)
  fig.text(.035,.90,'Four evidence categories sum to 100% of all landmark hospitalizations at every hour',fontsize=12)
  handles=None
  for j,g in enumerate(gs):
   ax=fig.add_subplot(grid[0,j]);handles=ax.stackplot(np.arange(73),np.c_[fractions[g],fractions[g][:,-1]],colors=COLORS,step='post',linewidth=0)
   ax.set(xlim=(0,72),ylim=(0,1),xticks=[0,24,48,72],yticks=[0,.5,1]);ax.yaxis.set_major_formatter(PercentFormatter(1));ax.set_title(titles[g]+f'\nN={int(den[g][0]):,}',fontsize=10.5,pad=12);ax.tick_params(labelbottom=False)
   if j==0:ax.set_ylabel('Hospitalizations (%)')
   else:ax.set_yticklabels([])
   bx=fig.add_subplot(grid[1,j]);bx.plot(np.arange(73),np.r_[coverage[g][1]/den[g],coverage[g][1,-1]/den[g][-1]],color='#999999',lw=1.7,label='At least one channel measured')
   bx.plot(np.arange(73),np.r_[coverage[g][0]/den[g],coverage[g][0,-1]/den[g][-1]],color='#222222',lw=1.7,label='Both channels measured')
   bx.set(xlim=(0,72),ylim=(0,1),xticks=[0,24,48,72],yticks=[0,.5,1],xlabel='Hours after ICU arrival');bx.yaxis.set_major_formatter(PercentFormatter(1));bx.grid(alpha=.12)
   if j==0:bx.set_ylabel('Evidence\ncoverage')
   else:bx.set_yticklabels([])
  fig.legend(handles, LABELS,loc='lower center',bbox_to_anchor=(.52,.18),ncol=2,frameon=False,fontsize=10)
  fig.text(.035,.045,'Gray includes measured criteria not met AND unknown/unmeasured physiology; it does not mean absence or resolution of respiratory failure.\n'
           '"Only" means only that channel has positive documented evidence; the other may be negative or unmeasured. Mixed requires both positive.\n'
           'Coverage below: light line = at least one channel evaluable; dark line = both evaluable. Each includes the existing six-hour evidence lookback.\n'
           'Exact pooled hourly counts; no smoothing, no patient exclusions, and no change to ARF definitions, clustering or model assignments.',fontsize=9,color='#444444',va='bottom')
  for ext in ['png','pdf']:fig.savefig(a.output/f'{stem}.{ext}',dpi=180,facecolor='white')
  plt.close(fig)
 plot(['P1','P2','P3','P4','I'],'combined_arf_evidence_by_phenotype')
 plot(['Overall'],'combined_arf_evidence_overall')
 manifest=dict(status='passed',source=str(source.relative_to(ROOT)),source_md5=hashlib.md5(source.read_bytes()).hexdigest(),primary_consensus_id=active['primary_consensus_id'],definition_changed=False,patient_level_inputs_read=False,display_grouping_only=True,hospitalizations=3310,all_hours_sum_to_one=True,gray_definition='No positive measured physiology documented: includes both negative and unknown channels; not no clinical ARF',unknown_channels_in_original_data_unchanged=True)
 (a.output/'manifest.json').write_text(json.dumps(manifest,indent=2));print(a.output)
if __name__=='__main__':main()
