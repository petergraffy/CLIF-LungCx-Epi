import argparse
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as RibbonPath
from matplotlib.patches import Patch, PathPatch, Rectangle
from matplotlib.ticker import PercentFormatter

parser = argparse.ArgumentParser(description="Render manuscript trajectories figure.")
parser.add_argument("pool", type=Path)
args = parser.parse_args()
run = args.pool
out = run / "figures"
out.mkdir(exist_ok=True)
support = {
    'ROOM_AIR': ('Room air', '#B7DAEF'), 'LOW_O2': ('Low-flow O₂', '#2F82B5'),
    'HFNC': ('HFNC', '#009C88'), 'CPAP': ('CPAP', '#83CCB5'),
    'NIPPV': ('NIPPV', '#146C50'), 'IMV': ('IMV', '#193759'),
    'OTHER': ('Other', '#6E8291'), 'UNKNOWN': ('Unknown', '#D4DBE1'),
}
arf = {
    'Hypoxemia evidence only': ('Hypoxemia only', '#DF5028'),
    'Hypercapnic-acidemic evidence only': ('Hypercapnia with acidemia only', '#F9A644'),
    'Mixed positive evidence': ('Mixed evidence', '#8F252E'),
    'No positive physiology documented': ('No positive evidence documented', '#F6ECD7'),
}
phenotypes = {
    'P1': ('Shift toward room air', '#7B3FA0'),
    'P2': ('Sustained low-flow oxygen', '#C93882'),
    'P3': ('High-flow predominant', '#B5C83C'),
    'P4': ('Sustained invasive ventilation', '#35251D'),
}

def load(file):
    with (run / file).open() as handle:
        return list(csv.DictReader(handle))

counts = load('figures/phenotype_trajectories_raw_counts.csv')
values = {(r['phenotype'], r['channel'], r['state'], int(r['hour'])): float(r['proportion']) for r in counts}
flows = {(r['initial_support'], r['phenotype']): int(r['n'])
         for r in load('reports/aggregate_collection/exact_counts/profiles/pooled_initial_support_to_phenotype_counts.csv')}
source_totals = {s: sum(flows.get((s, p), 0) for p in phenotypes) for s in support}
phenotype_totals = {p: sum(flows.get((s, p), 0) for s in support) for p in phenotypes}
for p, total in phenotype_totals.items():
    assert {int(r['denominator']) for r in counts if r['phenotype'] == p} == {total}

plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':9,
                     'axes.spines.top':False, 'axes.spines.right':False, 'svg.fonttype':'none'})
fig = plt.figure(figsize=(14.2, 8.5), facecolor='white')
fig.text(.021, .962, 'A', fontsize=17, fontweight='bold')
fig.text(.535, .962, 'B', fontsize=17, fontweight='bold')
fig.text(.171, .960, 'Respiratory support', ha='center', fontsize=12, fontweight='bold')
fig.text(.397, .960, 'Documented ARF evidence', ha='center', fontsize=12, fontweight='bold')
gs = fig.add_gridspec(4, 2, left=.066, right=.507, bottom=.184, top=.896, wspace=.13, hspace=.43)
hours = np.arange(72)
display_hours = np.r_[0, hours + .5, 72]
kernel = np.ones(5)
edge = np.convolve(np.ones(72), kernel, mode='same')
for i, (p, (name, _)) in enumerate(phenotypes.items()):
    for j, (channel, palette) in enumerate((('support', support), ('arf_evidence', arf))):
        ax = fig.add_subplot(gs[i,j])
        raw = np.array([[values[p, channel, state, hour] for hour in hours] for state in palette])
        smooth = np.stack([np.convolve(s, kernel, mode='same') / edge for s in raw])
        assert np.allclose(smooth.sum(axis=0), 1)
        stack = np.c_[smooth[:,0], smooth, smooth[:,-1]]
        ax.stackplot(display_hours, stack, colors=[c for _, c in palette.values()], linewidth=0)
        ax.set(xlim=(0,72), ylim=(0,1), xticks=[0,24,48,72], yticks=[0,.5,1])
        ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        ax.tick_params(axis='both', labelsize=8, length=3, pad=2,
                       labelbottom=i==3, labelleft=j==0, left=j==0)
        if j:
            ax.spines['left'].set_visible(False)
        if i == 3:
            ax.set_xlabel('Hours after ICU arrival', fontsize=9)
        if not j:
            bounds = ax.get_position()
            fig.text(bounds.x0, bounds.y1 + .012, f'{p}  {name}  (N = {phenotype_totals[p]:,})',
                     ha='left', fontsize=10, fontweight='bold')
for j, palette in enumerate((support, arf)):
    legend_ax = fig.add_axes((.065 if j==0 else .285, .030, .215, .092))
    legend_ax.set_axis_off()
    legend_ax.legend(handles=[Patch(facecolor=color,label=label) for label,color in palette.values()],
                     ncol=2 if j==0 else 1, loc='upper left', fontsize=8.3, frameon=False,
                     borderaxespad=0, columnspacing=1.4, labelspacing=.45, handlelength=1.5)

ax = fig.add_axes((.535, .113, .451, .827))
ax.set_axis_off()
ax.set(xlim=(0,1), ylim=(0,1))
left, right, node_width = .242, .650, .022
ax.text(left, .967, 'Initial respiratory support', ha='center', fontweight='bold', fontsize=10.5)
ax.text(right, .967, 'Assigned phenotype', ha='center', fontweight='bold', fontsize=10.5)
height = .755
scale = height / sum(source_totals.values())
def node_positions(totals, gap):
    pos = {}
    top = .907
    for key, total in totals.items():
        lower = top - total * scale
        pos[key] = (lower, top)
        top = lower - gap
    return pos
left_nodes = node_positions(source_totals, .016)
right_nodes = node_positions(phenotype_totals, .112/3)
left_cursors = {key: top for key, (_,top) in left_nodes.items()}
right_cursors = {key: top for key, (_,top) in right_nodes.items()}
a,b = left+node_width, right
bend = .53*(b-a)
for s, (_, color) in support.items():
    for p in phenotypes:
        n = flows.get((s,p),0)
        if not n:
            continue
        y0, y1 = left_cursors[s], right_cursors[p]
        h = n * scale
        points = [(a,y0),(a+bend,y0),(b-bend,y1),(b,y1),
                  (b,y1-h),(b-bend,y1-h),(a+bend,y0-h),(a,y0-h),(a,y0)]
        codes = [RibbonPath.MOVETO] + [RibbonPath.CURVE4]*3 + [RibbonPath.LINETO] + [RibbonPath.CURVE4]*3 + [RibbonPath.CLOSEPOLY]
        ax.add_patch(PathPatch(RibbonPath(points,codes), facecolor=color, alpha=.48, edgecolor='none'))
        left_cursors[s] -= h
        right_cursors[p] -= h
for s, (label,color) in support.items():
    lo, hi = left_nodes[s]
    assert abs(left_cursors[s] - lo) < 1e-9
    ax.add_patch(Rectangle((left,lo),node_width,hi-lo,facecolor=color,edgecolor='none'))
    ax.text(left-.015,(lo+hi)/2,f'{label} ({source_totals[s]:,})',ha='right',va='center',fontsize=8.8)
for p, (label,color) in phenotypes.items():
    lo, hi = right_nodes[p]
    assert abs(right_cursors[p] - lo) < 1e-9
    ax.add_patch(Rectangle((right,lo),node_width,hi-lo,facecolor=color,edgecolor='none'))
    ax.text(right+node_width+.013,(lo+hi)/2,f'{p} {label}\nN = {phenotype_totals[p]:,}',
            ha='left',va='center',fontsize=8.6, linespacing=1.3)
for ext in ('png','svg'):
    fig.savefig(out / f'Figure 1 - Hourly Trajectories and Sankey.{ext}', dpi=210, facecolor='white')
