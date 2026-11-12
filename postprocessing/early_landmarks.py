"""Plot the frozen 24-hour and 48-hour reference profiles."""
import argparse
from pathlib import Path
from collections import defaultdict
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("bundle", type=Path)
parser.add_argument("output", type=Path)
args = parser.parse_args()
bundle = args.bundle
out = args.output
out.mkdir(parents=True, exist_ok=True)
support = {
    'ROOM_AIR': ('Room air', '#B7DAEF'),
    'LOW_O2': ('Low-flow O₂', '#2F82B5'),
    'HFNC': ('HFNC', '#009C88'),
    'CPAP': ('CPAP', '#83CCB5'),
    'NIPPV': ('NIPPV', '#146C50'),
    'IMV': ('IMV', '#193759'),
    'OTHER': ('Other', '#6E8291'),
    'UNKNOWN': ('Unknown', '#D4DBE1'),
}
arf = {
    'hypox': ('Hypoxemia only', '#DF5028'),
    'hyper': ('Hypercapnia with acidemia only', '#F9A644'),
    'mixed': ('Mixed evidence', '#8F252E'),
    'none': ('No positive evidence documented', '#F6ECD7'),
}

def arf_category(state):
    hypox, hyper = state.split('|')[1].split('_')
    if hypox == 'H1' and hyper == 'C1':
        return 'mixed'
    if hypox == 'H1':
        return 'hypox'
    if hyper == 'C1':
        return 'hyper'
    return 'none'

def handles(palette):
    return [Patch(facecolor=color, label=label) for label, color in palette.values()]

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,
                     'axes.spines.right':False,'axes.linewidth':.8,'svg.fonttype':'none'})
for hours in (24,48):
    path = bundle / f'h{hours}' / 'coordinator/full'
    manifest = json.loads((path / 'frozen_manifest.json').read_text())
    step = manifest['specification']['bin_hours']
    bins = hours // step
    phenotypes = [f'P{i}' for i in range(1, manifest['k'] + 1)]
    prob = defaultdict(float)
    with (path / 'prototypes.csv').open() as stream:
        for row in csv.DictReader(stream):
            p, b, d, s = [row[k] for k in ('phenotype','bin','domain','state')]
            if d == 'support':
                prob[p,d,s,int(b)] += float(row['probability'])
            elif d == 'joint':
                prob[p,'arf',arf_category(s),int(b)] += float(row['probability'])

    fig, axes = plt.subplots(2,len(phenotypes),figsize=(12,5.25),sharey='row',sharex=True)
    fig.subplots_adjust(left=.092,right=.987,top=.89,bottom=.15,wspace=.15,hspace=.50)
    edges = np.arange(0,hours+1,step)
    for j, p in enumerate(phenotypes):
        for i, (domain, palette) in enumerate((('support',support),('arf',arf))):
            ax = axes[i,j]
            layers = np.array([[prob[p,domain,state,b] for b in range(bins)] for state in palette])
            assert np.allclose(layers.sum(axis=0),1) and (layers>=0).all(), (hours,p,domain)
            ax.stackplot(edges,np.c_[layers,layers[:,-1]],step='post',
                         colors=[c for _,c in palette.values()],linewidth=0)
            ax.set(xlim=(0,hours),ylim=(0,1),xticks=np.arange(0,hours+1,hours//3),yticks=[0,.5,1])
            ax.yaxis.set_major_formatter(PercentFormatter(1,decimals=0))
            ax.tick_params(axis='both', labelsize=8, length=3)
            if j:
                ax.tick_params(axis='y',left=False)
                ax.spines['left'].set_visible(False)
            if i==0:
                ax.set_title(p,fontweight='bold',fontsize=11)
    axes[0,0].set_ylabel('Respiratory support\nprofile probability')
    axes[1,0].set_ylabel('ARF evidence\nprofile probability')
    fig.text(.54,.087,'Hours after first ICU entry',ha='center')
    fig.legend(handles=handles(support),loc='center',bbox_to_anchor=(.54,.505),
               ncol=8,frameon=False,columnspacing=1.2,fontsize=8.2,handlelength=1.45)
    fig.legend(handles=handles(arf),loc='lower center',bbox_to_anchor=(.54,.028),
               ncol=4,frameon=False,columnspacing=1.35,fontsize=8.2,handlelength=1.45)
    for ext in ('png','svg'):
        fig.savefig(out/f'Figure S{6 if hours==24 else 7} - {hours} Hour Frozen Profiles.{ext}',
                    dpi=220,facecolor='white')
    plt.close(fig)
