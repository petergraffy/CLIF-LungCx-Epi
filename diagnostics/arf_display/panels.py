"""Matched standalone A/B trajectory panels; Sankey remains a separate graphic."""
import argparse,csv,hashlib,json,os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/clif_manuscript_mpl')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter
from main_figure import smooth_simplex,PH,PN,SUPPORT,SC,SL
from combined import LABELS,COLORS
ROOT=Path(__file__).resolve().parents[2]
p=argparse.ArgumentParser();p.add_argument('output',type=Path);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
active=json.loads((ROOT/'config/active_analysis.json').read_text());source=ROOT/active['main_phenotype_figure'];source=source.with_name(source.stem+'_raw_counts.csv')
with source.open() as f:rows=list(csv.DictReader(f))
data={};N={}
for channel,states in [('support',SUPPORT),('arf_evidence',LABELS)]:
 for ph in PH:
  z=[r for r in rows if r['channel']==channel and r['phenotype']==ph];values=np.zeros((len(states),72));seen=set()
  for r in z:
   key=(r['state'],int(r['hour']));assert key not in seen;seen.add(key)
   values[states.index(r['state']),int(r['hour'])]=float(r['proportion']);N[ph]=int(r['denominator'])
  assert len(seen)==len(states)*72 and np.allclose(values.sum(axis=0),1);data[(channel,ph)]=smooth_simplex(values)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
def draw(fig,axes,channel,top):
 colors,labels=(SC,SL) if channel=='support' else (COLORS,LABELS)
 for j,(ax,ph) in enumerate(zip(axes,PH)):
  x,z=data[(channel,ph)];ax.stackplot(x,z,colors=colors,linewidth=0)
  ax.set(xlim=(0,72),ylim=(0,1),xticks=[0,24,48,72],yticks=[0,.5,1],xlabel='Hours after ICU arrival');ax.yaxis.set_major_formatter(PercentFormatter(1))
  ax.set_title(ph+' '+PN[j]+f'\nN={N[ph]:,}',fontsize=11,pad=12)
  if j==0:ax.set_ylabel('Hospitalizations (%)')
  else:ax.set_yticklabels([])
 return [Patch(facecolor=c,label=l) for c,l in zip(colors,labels)]
for channel,letter,title,name in [('support','A','Respiratory support trajectories','panel_A_respiratory_support'),('arf_evidence','B','Documented ARF evidence trajectories','panel_B_arf_evidence')]:
 fig,axes=plt.subplots(1,4,figsize=(14,5.7));fig.subplots_adjust(left=.075,right=.98,top=.72,bottom=.32,wspace=.2)
 fig.suptitle(letter+'   '+title,x=.035,ha='left',fontsize=21,y=.97)
 fig.text(.035,.88,f'Pooled UCMC + NU + Rush | {sum(N.values()):,} hospitalizations in P1-P4 | five-hour display smoothing',fontsize=12)
 handles=draw(fig,axes,channel,.72)
 fig.legend(handles=handles,loc='center',bbox_to_anchor=(.52,.19),ncol=8 if channel=='support' else 2,frameon=False,fontsize=10)
 note='Smoothed observed hourly proportions; five-hour centered average with convex interpolation. No changes to hourly states, assignments or models.\nThe 130 insufficient-information stays are omitted from the figure only and remain in the analysis.'
 if channel=='arf_evidence':note+='\nGray includes measured criteria not met and unmeasured physiology; it does not establish absence of respiratory failure. The six-hour evidence rule is unchanged.'
 fig.text(.035,.035,note,fontsize=9,color='#444444',va='bottom')
 for ext in ['png','svg']:fig.savefig(a.output/f'{name}.{ext}',dpi=200,facecolor='white')
 plt.close(fig)
# A publication-layout option combining the same two panels, without the Sankey.
fig=plt.figure(figsize=(14,10));fig.text(.035,.974,'A   Respiratory support trajectories',fontsize=19)
axes=[fig.add_axes([.075+j*.23,.635,.195,.25]) for j in range(4)];handles=draw(fig,axes,'support',.885)
fig.legend(handles=handles,loc='center',bbox_to_anchor=(.53,.55),ncol=8,frameon=False,fontsize=10)
fig.text(.035,.482,'B   Documented ARF evidence trajectories',fontsize=19)
axes=[fig.add_axes([.075+j*.23,.16,.195,.25]) for j in range(4)];handles=draw(fig,axes,'arf_evidence',.41)
fig.legend(handles=handles,loc='center',bbox_to_anchor=(.53,.075),ncol=2,frameon=False,fontsize=9.5)
fig.text(.035,.012,'P1-P4: 3,180 hospitalizations. Display smoothing only; six-hour physiology rules and models are unchanged. I (N=130) is omitted from this figure only.\nARF gray includes negative and unmeasured physiology; it does not mean no clinical respiratory failure.',fontsize=9,color='#444444')
for ext in ['png','svg']:fig.savefig(a.output/f'figure_1_AB.{ext}',dpi=200,facecolor='white')
plt.close(fig)
(a.output/'manifest.json').write_text(json.dumps(dict(status='passed',primary_consensus_id=active['primary_consensus_id'],source=str(source.relative_to(ROOT)),source_md5=hashlib.md5(source.read_bytes()).hexdigest(),hospitalizations_shown=sum(N.values()),phenotype_order=PH,panels={'A':'respiratory support','B':'documented ARF evidence'},sankey_included=False,definitions_or_models_changed=False),indent=2))
print(a.output)
