"""Verify the shared tree contains runnable source, not results or retired data."""
from pathlib import Path
import re, subprocess
root=Path(__file__).resolve().parents[1]
tracked=subprocess.check_output(['git','ls-files'],cwd=root,text=True).splitlines()
tracked_set=set(tracked)
banned=('output/','sites/','exposome/','resources/','outlier-thresholds/','tools/','archive/','tmp/','renv/','docx_render_')
for name in tracked:
    assert not name.startswith(banned), f'Generated or retired material tracked: {name}'
    assert Path(name).suffix.lower() not in {'.rds','.parquet','.fst','.rdata','.zip','.shp','.dbf','.png','.tif'}, name
    if Path(name).suffix in {'.R','.py','.sh','.json'}:
        text=(root/name).read_text()
        assert not re.search(r'/(?:Users|home)/[^/\s]+/',text), f'Private machine path in {name}'
        for dependency in re.findall(r'''source\(['"]((?:code|scripts|diagnostics|tests)/[^'"]+\.R)['"]''',text):
            assert dependency in tracked_set, f'{name} imports untracked {dependency}'
for name in ['README.md','code/README.md','config/README.md']:
    text=(root/name).read_text()
    for target in re.findall(r'\]\(([^)]+)\)',text):
        if not target.startswith(('http:','https:','#')):
            resolved=(root/name).parent/target
            rel=resolved.resolve().relative_to(root.resolve()).as_posix()
            assert rel in tracked_set, f'{name} links to untracked {target}'
print(f'PASS: {len(tracked)} tracked files; no output/site/exposome data, retired assets, private machine paths or missing literal R imports')
