"""Project unused code_status hospitalization ID out of a local CLIF table view.

Use only for the current loader's selected-column/filter mismatch. Source
files and every code status row/value are preserved; other tables symlink to
source. No analytical code or clinical rule is changed.
"""
import json
from pathlib import Path
import sys
import pyarrow.parquet as pq

config_path = Path(sys.argv[1]).resolve()
config = json.loads(config_path.read_text())
root = config_path.parent
manifest = root / "schema_adapter.json"
if manifest.exists():
    raise SystemExit("Schema view already exists; inspect schema_adapter.json")
source = Path(config["tables_path"])
view = root / "clif_schema_view"
view.mkdir(exist_ok=False)
(root / "source_config.json").write_text(json.dumps(config, indent=2))
for path in source.glob("*.parquet"):
    dest = view / path.name
    if path.name == "clif_code_status.parquet":
        table = pq.read_table(path)
        assert "hospitalization_id" in table.column_names
        pq.write_table(table.drop(["hospitalization_id"]), dest)
    else:
        dest.symlink_to(path)
config["tables_path"] = str(view)
config_path.write_text(json.dumps(config, indent=2))
manifest.write_text(json.dumps({
    "source_tables_path": str(source), "view_tables_path": str(view),
    "change": "Drop unused hospitalization_id from code_status projection; preserve patient_id, timestamps and all code status values. All other tables are symlinks to unchanged source.",
    "reason": "Core loader selects only patient-level code status fields, but tries filtering a discarded hospitalization_id when present.",
    "code_status_source_rows": len(table), "code_status_output_rows": len(table)
}, indent=2))
