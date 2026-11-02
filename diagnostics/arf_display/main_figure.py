"""Main phenotype display: smoothed observed support and documented ARF evidence."""
import argparse,csv,hashlib,json,os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/clif_manuscript_mpl')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch,PathPatch,Rectangle
from matplotlib.path import Path as MPath
from matplotlib.ticker import PercentFormatter
from combined import category,LABELS,COLORS
ROOT=Path(__file__).resolve().parents[2]
PH=['P1','P2','P3','P4'];PN=['Room-air predominant','Low-flow oxygen','HFNC predominant','IMV predominant']
SUPPORT=['ROOM_AIR','LOW_O2','HFNC','CPAP','NIPPV','IMV','OTHER','UNKNOWN']
SC=['#B8D4E8','#347DB1','#4FA343','#A45BB5','#713B94','#E65D0A','#D4B15C','#999999']
SL=['Room air','Low-flow O2','HFNC','CPAP','NIPPV','IMV','Other','Unknown']

def smooth_simplex(values):
 """Common five-hour nonnegative moving-average kernel, then convex interpolation."""
 values=np.asarray(values,dtype=float);assert values.ndim==2 and values.shape[1]==72
 assert np.all(values>=0) and np.allclose(values.sum(axis=0),1)
 kernel=np.ones(5);normalizer=np.convolve(np.ones(72),kernel,mode='same')
 smoothed=np.array([np.convolve(v,kernel,mode='same')/normalizer for v in values])
 centers=np.arange(72)+.5;x=np.linspace(0,72,721);loc=np.clip(x,centers[0],centers[-1])
 i=np.clip(np.searchsorted(centers,loc,side='right')-1,0,70);t=(loc-centers[i])/(centers[i+1]-centers[i]);w=3*t*t-2*t*t*t
 z=smoothed[:,i]*(1-w)+smoothed[:,i+1]*w
 assert np.all(z>=0) and np.allclose(z.sum(axis=0),1)
 return x,z

def read(p):
 with p.open() as f:return list(csv.DictReader(f))

def create(output,active):
 output=Path(output).resolve();output.parent.mkdir(parents=True,exist_ok=True)
 base=ROOT/active['exact_root']/'reports/aggregate_collection/exact_counts/profiles'
 sources=[base/'hourly_support.csv',base/'hourly_joint.csv',ROOT/active['sankey_root']/'pooled_initial_support_to_phenotype_counts.csv']
 before={str(p.relative_to(ROOT)):hashlib.md5(p.read_bytes()).hexdigest() for p in sources}
 support={g:np.zeros((8,72)) for g in PH};arf={g:np.zeros((4,72)) for g in PH};all_den={}
 for r in read(sources[0]):
  ph=r['phenotype'];assert ph in PH+['I'],'Unexpected phenotype: update the main-figure layout';h=int(r['h']);n=int(r['n']);all_den[(ph,h)]=all_den.get((ph,h),0)+n
  if ph in PH:support[ph][SUPPORT.index(r['state']),h]+=n
 for r in read(sources[1]):
  ph=r['phenotype']
  if ph in PH:
   h,c=r['state'].split('|')[1].split('_');arf[ph][category(h,c),int(r['h'])]+=int(r['n'])
 Ns={g:int(support[g].sum(axis=0)[0]) for g in PH};shown=sum(Ns.values());omitted=all_den.get(('I',0),0)
 assert shown>0 and omitted>=0
 for g in PH:
  assert np.all(support[g].sum(axis=0)==Ns[g]) and np.all(arf[g].sum(axis=0)==Ns[g])
 flows={(r['initial_support'],r['phenotype']):int(r['n']) for r in read(sources[2]) if r['phenotype'] in PH};assert sum(flows.values())==shown
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
 fig=plt.figure(figsize=(12.5,7))
 export=[]
 for row,(matrices,colors,label,bottom,height) in enumerate([(support,SC,'Respiratory support',.62,.21),(arf,COLORS,'Documented ARF evidence',.25,.23)]):
  for j,g in enumerate(PH):
   ax=fig.add_axes([.085+j*.225,bottom,.195,height]);raw=matrices[g]/Ns[g];x,z=smooth_simplex(raw);ax.stackplot(x,z,colors=colors,linewidth=0)
   ax.set(xlim=(0,72),ylim=(0,1),xticks=[0,24,48,72],yticks=[0,.5,1]);ax.yaxis.set_major_formatter(PercentFormatter(1));ax.tick_params(labelsize=9)
   if j==0:ax.set_ylabel(label+' (%)',fontsize=10)
   else:ax.set_yticklabels([])
   if row==0:ax.set_title(g+' '+PN[j]+f'\nN={Ns[g]:,}',fontsize=11,pad=13);ax.set_xticklabels([])
   else:ax.set_xlabel('Hours after ICU arrival',fontsize=9)
   states=SUPPORT if row==0 else LABELS
   for k,st in enumerate(states):
    for h in range(72):export.append(dict(phenotype=g,channel='support' if row==0 else 'arf_evidence',state=st,hour=h,n=int(matrices[g][k,h]),denominator=Ns[g],proportion=raw[k,h]))
 fig.legend([Patch(facecolor=c) for c in SC],SL,loc='center',bbox_to_anchor=(.53,.555),ncol=8,frameon=False,fontsize=8.6)
 fig.legend([Patch(facecolor=c) for c in COLORS],LABELS,loc='center',bbox_to_anchor=(.53,.145),ncol=2,frameon=False,fontsize=9)
 fig.savefig(output,dpi=200,facecolor='white',bbox_inches='tight');fig.savefig(output.with_suffix('.svg'),facecolor='white',bbox_inches='tight');plt.close(fig)
 # Separate Sankey export, using the same displayed clinical-phenotype cohort.
 sankey_output=output.with_name('standalone_initial_support_sankey.png')
 sf=plt.figure(figsize=(12.5,6))
 ax=sf.add_axes([.035,.16,.93,.64]);ax.axis('off');ax.set(xlim=(0,1),ylim=(0,1))
 src=[s for s in SUPPORT if sum(n for (a,b),n in flows.items() if a==s)>0];tot={s:sum(n for (a,b),n in flows.items() if a==s) for s in src};scale=.65/shown;gap=.04
 def nodes(keys,counts):
  used=scale*sum(counts.values())+gap*(len(keys)-1);top=.5+used/2;out={}
  for k in keys:h=counts[k]*scale;out[k]=(top-h,top);top-=h+gap
  return out
 left=nodes(src,tot);right=nodes(PH,Ns);lc={s:left[s][1] for s in src};rc={g:right[g][1] for g in PH};color=dict(zip(SUPPORT,SC));label=dict(zip(SUPPORT,SL))
 for s in src:
  for g in PH:
   h=flows.get((s,g),0)*scale
   if not h:continue
   y0,y1=lc[s],rc[g];lc[s]-=h;rc[g]-=h
   vs=[(.26,y0),(.42,y0),(.57,y1),(.73,y1),(.73,y1-h),(.57,y1-h),(.42,y0-h),(.26,y0-h),(.26,y0)]
   codes=[MPath.MOVETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4,MPath.LINETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4,MPath.CLOSEPOLY]
   ax.add_patch(PathPatch(MPath(vs,codes),facecolor=color[s],edgecolor='none',alpha=.5))
 for s,(lo,hi) in left.items():ax.add_patch(Rectangle((.245,lo),.015,hi-lo,facecolor=color[s]));ax.text(.234,(lo+hi)/2,f'{label[s]} ({tot[s]:,})',ha='right',va='center',fontsize=8.6)
 for g,(lo,hi) in right.items():ax.add_patch(Rectangle((.73,lo),.015,hi-lo,facecolor={'P1':SC[0],'P2':SC[1],'P3':SC[2],'P4':SC[5]}[g]));ax.text(.758,(lo+hi)/2,f'{g} ({Ns[g]:,})',ha='left',va='center',fontsize=9.5)
 sf.savefig(sankey_output,dpi=200,facecolor='white',bbox_inches='tight');sf.savefig(sankey_output.with_suffix('.svg'),facecolor='white',bbox_inches='tight');plt.close(sf)
 with output.with_name(output.stem+'_raw_counts.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(export[0]));w.writeheader();w.writerows(export)
 after={str(p.relative_to(ROOT)):hashlib.md5(p.read_bytes()).hexdigest() for p in sources};assert before==after
 audit=dict(status='passed',captions_and_subtitles=False,sankey_in_trajectory_figure=False,standalone_sankey=str(sankey_output.relative_to(ROOT)),primary_consensus_id=active['primary_consensus_id'],shown_hospitalizations=shown,omitted_from_figure_only=omitted,shown_groups=PH,phenotype_n=Ns,display_source='Pooled exact observed hourly proportions for both support and physiology; not equal-site consensus templates',smoothing='Five-hour centered moving average with edge-normalized common weights; convex interpolation between hourly centers',all_smoothed_stacks_sum_to_one=True,arf_definition_changed=False,clustering_or_models_changed=False,underlying_unknown_states_unchanged=True,source_md5=before)
 output.with_name(output.stem+'_manifest.json').write_text(json.dumps(audit,indent=2));return audit
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('output',type=Path);a=p.parse_args();active=json.loads((ROOT/'config/active_analysis.json').read_text());print(json.dumps(create(a.output,active),indent=2))
