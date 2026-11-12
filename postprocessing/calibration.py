"""Plot calibration and decision curves from returned aggregate bins."""
import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("pool", type=Path)
parser.add_argument("returns", nargs="+", type=Path)
args = parser.parse_args()
out = args.pool / "figures"
out.mkdir(exist_ok=True)

def read(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))

def packet(relative_path):
    return [row for path in args.returns for row in read(path / relative_path)]

matched_n = sum(
    int(row["n"]) for row in read(args.pool / "prediction/site_cohort_counts.csv")
    if row["cohort"] == "matched_complete_SOFA" and row["endpoint"] == "death_or_hospice"
)

cohort='matched_complete_SOFA'
models={'baseline':('Clinical baseline','#424B54','--'),'baseline_sofa':('Baseline + SOFA','#217D99','-'),
        'baseline_adapted':('Baseline + APACHE II','#A35917','-'),'phenotype':('Baseline + phenotype','#71388F','-')}
ep={'death':'In-hospital death','death_or_hospice':'Death or hospice discharge'}
cal=defaultdict(lambda:[0,0.,0])
for r in packet('prediction/calibration_bins.csv'):
    if r['cohort']==cohort and r['model'] in models:
        key=(r['endpoint'],r['model'],int(r['bin']))
        cal[key][0]+=int(r['n']); cal[key][1]+=float(r['predicted_sum']); cal[key][2]+=int(r['events'])
dc=defaultdict(lambda:[0,0,0])
for r in packet('prediction/decision_curves.csv'):
    if r['cohort']==cohort and r['model'] in models:
        key=(r['endpoint'],r['model'],float(r['threshold']))
        dc[key][0]+=int(r['n']); dc[key][1]+=int(r['true_positive']); dc[key][2]+=int(r['false_positive'])
fig,grid=plt.subplots(2,2,figsize=(9.6,5.55))
fig.subplots_adjust(left=.10,right=.97,top=.90,bottom=.15,hspace=.55,wspace=.27)
for j,(endpoint,label) in enumerate(ep.items()):
    ax=grid[0,j]
    ax.plot([0,1],[0,1],color='#777777',linestyle=':',lw=1)
    for model,(name,col,ls) in models.items():
        data=sorted((key[2],*v) for key,v in cal.items() if key[:2]==(endpoint,model))
        assert sum(v[1] for v in data)==matched_n
        ax.plot([v[2]/v[1] for v in data],[v[3]/v[1] for v in data],marker='o',
                markersize=3,ls=ls,lw=1.35,color=col,label=name)
    ax.set(xlim=(0,1),ylim=(0,1),xlabel='Mean predicted risk',ylabel='Observed event proportion')
    ax.set_title(label,fontsize=10,fontweight='bold')
    ax=grid[1,j]
    ax.axhline(0,color='#777777',linestyle=':',lw=1)
    for model,(name,col,ls) in models.items():
        data=sorted((key[2],*v) for key,v in dc.items() if key[:2]==(endpoint,model))
        assert all(v[1]==matched_n for v in data)
        ax.plot([v[0] for v in data],[(v[2]-v[3]*v[0]/(1-v[0]))/v[1] for v in data],
                ls=ls,lw=1.4,color=col,label=name)
    ax.set(xlabel='Risk threshold',ylabel='Net benefit')
handles,labels=grid[0,0].get_legend_handles_labels()
fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,.024),ncol=4,frameon=False,fontsize=8.2)
fig.savefig(out/'calibration_and_net_benefit.png',dpi=220,facecolor='white')
