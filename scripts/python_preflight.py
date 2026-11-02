"""Check the pinned waterfall runtime without reading clinical data."""
import hashlib, importlib.metadata, inspect, json, sys
from pathlib import Path
from clifpy.utils.waterfall import process_resp_support_waterfall
if sys.version_info[:2] != (3, 12):
    raise SystemExit('Use Python 3.12 for the waterfall environment.')
for line in Path('config/waterfall_requirements.txt').read_text().splitlines():
    if not line or line.startswith('#'):
        continue
    name, expected = line.split('==')
    actual = importlib.metadata.version(name)
    if actual != expected:
        raise SystemExit(f'{name}: expected {expected}, found {actual}; restore the requirements file.')
print(json.dumps({'python':sys.version.split()[0], 'clifpy':importlib.metadata.version('clifpy'),
                  'waterfall_sha256':hashlib.sha256(Path(inspect.getfile(process_resp_support_waterfall)).read_bytes()).hexdigest()}))
