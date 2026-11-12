import argparse
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, PercentFormatter

parser = argparse.ArgumentParser(description="Render manuscript competing risks figure.")
parser.add_argument("pool", type=Path)
args = parser.parse_args()
root = args.pool
output = root / "figures"
output.mkdir(exist_ok=True)
phenotypes = {
    'P1': ('Shift toward room air', '#7B3FA0'),
    'P2': ('Sustained low-flow oxygen', '#C93882'),
    'P3': ('High-flow predominant', '#B5C83C'),
    'P4': ('Sustained invasive ventilation', '#35251D'),
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

# Figure 2 — pooled cumulative incidence and risk sets.
reports = 'survival/reports'
incidence = read(f'{reports}/pooled_competing_incidence.csv')
fig = plt.figure(figsize=(10.4, 5.6), facecolor='white')
panels = [fig.add_axes((x, .38, .40, .51)) for x in (.092, .573)]
risk_by_endpoint = {}
risk_days = (0, 7, 14, 21, 28)
for i, (ax, (endpoint, label)) in enumerate(zip(panels, (
    ('death_or_hospice', 'Death or hospice discharge'),
    ('live_discharge', 'Live discharge without hospice'),
))):
    risk_cells = []
    for p, (name, color) in phenotypes.items():
        data = sorted((r for r in incidence if r['curve'] == 'landmark72h_' + p and r['endpoint'] == endpoint),
                      key=lambda r: float(r['time']))
        times = [float(r['time']) for r in data]
        vals = np.asarray([float(r['estimate']) for r in data])
        assert vals[0] == 0 and np.all(np.diff(vals) >= -1e-12)
        lo = [float(r['conf_low']) if r['conf_low'] else np.nan for r in data]
        hi = [float(r['conf_high']) if r['conf_high'] else np.nan for r in data]
        ax.fill_between(times, lo, hi, step='post', color=color, alpha=.13, linewidth=0)
        ax.step(times, vals, where='post', lw=1.85, color=color, label=f'{p} {name}')
        risk = {float(r['time']): int(r['n_risk']) for r in data}
        risk_cells.append([f'{risk[d]:,}' for d in risk_days])
    ax.set(xlim=(0, 28), ylim=(0, 1), xticks=risk_days, yticks=np.linspace(0, 1, 6),
           xlabel='Days after 72-hour landmark', ylabel='Cumulative incidence')
    ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    title(fig, ax, 'AB'[i], label)
    risk_by_endpoint[endpoint] = risk_cells

assert risk_by_endpoint['death_or_hospice'] == risk_by_endpoint['live_discharge']
risk_cells = risk_by_endpoint['death_or_hospice']
risk_ax = fig.add_axes((.265, .128, .47, .123))
risk_ax.set_axis_off()
risk_ax.text(.5, 1.10, 'Number at risk', fontsize=8.2, fontweight='bold', ha='center', va='bottom',
             transform=risk_ax.transAxes)
table = risk_ax.table(cellText=risk_cells, colLabels=[f'Day {d}' for d in risk_days],
                      rowLabels=list(phenotypes), cellLoc='center', rowLoc='right',
                      bbox=[.07, -.01, .93, 1.05])
table.auto_set_font_size(False)
table.set_fontsize(7.8)
for (row, col), cell in table.get_celld().items():
    cell.visible_edges = 'B' if row == 0 else ''
    cell.set_edgecolor('#A3A9AD')
    cell.set_linewidth(.45)
    if row == 0 or col == -1:
        cell.set_text_props(fontweight='bold')
shared_legend(fig, *panels[0].get_legend_handles_labels())
save(fig, 'Figure 2 - Aalen Johansen Competing Risks')
