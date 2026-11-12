"""Pool exact hourly support, joint states, and initial-support flows."""
from pathlib import Path
from collections import Counter, defaultdict
import argparse
import csv
import hashlib
import json
import sys

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("output", type=Path)
parser.add_argument("bundle", type=Path)
parser.add_argument("returns", nargs="+", type=Path)
args = parser.parse_args()
base = args.output.resolve()
repo = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo / 'diagnostics/arf_display'))
from main_figure import create

paths = args.returns
profiles = base / 'reports/aggregate_collection/exact_counts/profiles'
profiles.mkdir(parents=True, exist_ok=True)
for name, keys in [
    ('hourly_support', ('phenotype', 'h', 'state')),
    ('hourly_joint', ('phenotype', 'h', 'state')),
    ('initial_support_to_phenotype', ('initial_support', 'phenotype')),
]:
    counts = Counter()
    for path in paths:
        manifest = json.loads((path / 'return_manifest.json').read_text())
        rel = f'profiles/{name}.csv'
        source = path / rel
        assert hashlib.md5(source.read_bytes()).hexdigest() == manifest['files'][rel]
        with source.open() as handle:
            rows = list(csv.DictReader(handle))
        assert rows and all(row['site'] == manifest['site'] for row in rows)
        if name.startswith('hourly'):
            sums = defaultdict(int)
            for row in rows:
                sums[row['phenotype'], row['h']] += int(row['n'])
            assert all(sums[row['phenotype'], row['h']] == int(row['denominator']) for row in rows)
        for row in rows:
            counts[tuple(row[key] for key in keys)] += int(row['n'])
    filename = f'{name}.csv' if name.startswith('hourly') else 'pooled_initial_support_to_phenotype_counts.csv'
    with (profiles / filename).open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow([*keys, 'n'])
        writer.writerows((*key, n) for key, n in sorted(counts.items()))
active = {
    'exact_root': str(base.relative_to(repo)),
    'sankey_root': str(profiles.relative_to(repo)),
     'primary_consensus_id': json.loads(
        (args.bundle / "bundle_manifest.json").read_text()
    )["primary_consensus_id"],
}
print(json.dumps(create(base / 'figures/phenotype_trajectories.png', active), indent=2))
