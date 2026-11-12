import argparse
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, PercentFormatter

parser = argparse.ArgumentParser(description="Render manuscript discrimination figure.")
parser.add_argument("pool", type=Path)
args = parser.parse_args()
root = args.pool
output = root / "figures"
output.mkdir(exist_ok=True)
models = {
    'baseline': ('Clinical baseline', '#424B54', '--'),
    'baseline_sofa': ('Baseline + SOFA', '#217D99', '-'),
    'baseline_adapted': ('Baseline + APACHE II', '#A35917', '-'),
    'phenotype': ('Baseline + phenotype', '#71388F', '-'),
}
plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.spines.top': False,
    'axes.spines.right': False, 'axes.linewidth': .85, 'svg.fonttype': 'none',
    'xtick.major.size': 3.5, 'ytick.major.size': 3.5,
})

def read(path):
    with (root / path).open() as stream:
        return list(csv.DictReader(stream))

def title(fig, ax, letter, text):
    pos = ax.get_position()
    ax.set_title(text, fontsize=11, fontweight='bold', pad=9)
    fig.text(pos.x0 - .058, pos.y1 + .022, letter, fontsize=14, fontweight='bold', ha='left')
    ax.grid(False)

def save(fig, name):
    for extension in ('svg', 'png'):
        fig.savefig(output / f'{name}.{extension}', dpi=240, facecolor='white')
    plt.close(fig)

def shared_legend(fig, handles, labels, columns=2):
    return fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(.5, .022), ncol=columns,
                      frameon=False, fontsize=9, columnspacing=3, handlelength=2.4)

# Figure 3 — matched cohort discrimination for hospital death or hospice.
def prediction_rows(name):
    return [r for r in read(f'prediction/{name}') if r['cohort'] == 'matched_complete_SOFA'
            and r['endpoint'] == 'death_or_hospice' and r['model'] in models]
roc = prediction_rows('pooled_roc.csv')
auc = prediction_rows('pooled_tdauc.csv')
cohort_n = sum(int(r['n']) for r in read('prediction/site_cohort_counts.csv')
               if r['cohort'] == 'matched_complete_SOFA' and r['endpoint'] == 'death_or_hospice')
assert cohort_n > 0 and {int(r['n']) for r in auc} == {cohort_n}
assert len(auc) == 4 * 28 and all(r['status']=='ok' for r in auc)
aucs = {(r['model'], int(r['day'])): float(r['auc']) for r in auc}
fig = plt.figure(figsize=(9.4, 7.5), facecolor='white')
panels = [fig.add_axes((x, y, .395, .332)) for y in (.602, .162) for x in (.095, .570)]
checks = []
for i, (ax, day) in enumerate(zip(panels, (7, 14, 28, None))):
    if day is not None:
        ax.plot([0, 1], [0, 1], color='#888888', ls=':', lw=1.05, zorder=0)
        handles = []
        for m, (_, color, ls) in models.items():
            points = [r for r in roc if r['model'] == m and int(r['day']) == day]
            assert points
            line, = ax.plot([float(r['fpr']) for r in points], [float(r['sensitivity']) for r in points],
                            color=color, ls=ls, lw=1.75)
            handles.append(line)
        ax.set(xlim=(0, 1), ylim=(0, 1), xlabel='False-positive rate', ylabel='Sensitivity')
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_locator(MultipleLocator(.25))
            axis.set_major_formatter(FormatStrFormatter('%g'))
        title(fig, ax, 'ABCD'[i], f'Day {day}')
        checks.append((ax, handles, day))
    else:
        ax.axhline(.5, color='#888888', ls=':', lw=1.05, zorder=0)
        for m, (label, color, ls) in models.items():
            ax.plot(range(1, 29), [aucs[m, d] for d in range(1, 29)], label=label,
                    color=color, ls=ls, lw=1.75)
        ax.set(xlim=(1, 28), ylim=(min(.47, min(aucs.values())-.02), max(.72, max(aucs.values())+.02)),
               xticks=[1, 7, 14, 21, 28],
               xlabel='Days after 72-hour landmark', ylabel='AUC')
        ax.yaxis.set_major_locator(MultipleLocator(.05))
        title(fig, ax, 'D', 'Time-dependent AUC')
shared_legend(fig, *panels[-1].get_legend_handles_labels(), columns=4)
for ax, handles, day in checks:
    labels = [f'{aucs[m, day]:.3f}' for m in models]
    for location in ('lower right', 'upper left', 'lower left', 'upper right'):
        box = ax.legend(handles, labels, title='AUC', loc=location, alignment='left',
                        fontsize=8.2, title_fontsize=8.2, handlelength=1.7, borderpad=.55,
                        labelspacing=.33, framealpha=1, facecolor='white', edgecolor='#C9CDD0')
        fig.canvas.draw()
        bbox = box.get_window_extent(fig.canvas.get_renderer()).expanded(1.04, 1.04)
        if not any(line.get_path().transformed(line.get_transform()).intersects_bbox(bbox, filled=False)
                   for line in ax.lines):
            break
    else:
        box.remove()
        # The numerical comparison remains visible when curves use all four corners.
        ax.legend(handles, labels, title='AUC', loc='center', bbox_to_anchor=(.5, -.26),
                  ncol=4, fontsize=7, title_fontsize=7, handlelength=1.3,
                  columnspacing=.8, frameon=False)
save(fig, 'Figure 3 - Pooled Discrimination')
