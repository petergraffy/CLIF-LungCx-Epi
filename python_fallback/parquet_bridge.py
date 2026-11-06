"""Local-only prototype: Parquet <-> CSV strings + binary doubles.

No patient data go to stdout. Requests and intermediate files stay local.
This is deliberately outside the released Stage 2 source contract.
"""
import csv
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq


def field_info(field):
    t = field.type
    if pa.types.is_dictionary(t):
        t = t.value_type
    if pa.types.is_timestamp(t):
        kind = "timestamp"
    elif pa.types.is_date(t):
        kind = "date"
    elif pa.types.is_boolean(t):
        kind = "boolean"
    elif pa.types.is_integer(t):
        kind = "integer"
    elif pa.types.is_floating(t) or pa.types.is_decimal(t):
        kind = "double"
    elif pa.types.is_string(t) or pa.types.is_large_string(t) or pa.types.is_null(t):
        kind = "string"
    else:
        kind = "unsupported"
    return dict(name=field.name, type=str(field.type), kind=kind,
                timezone=getattr(t, "tz", None), unit=getattr(t, "unit", None),
                bits=getattr(t, "bit_width", None) if kind == "integer" else None)


def read(req):
    dataset = ds.dataset(req["path"], format="parquet")
    columns = req.get("columns") or dataset.schema.names
    fields = [field_info(dataset.schema.field(n)) for n in columns]
    predicate = None
    if req.get("id_column"):
        # Compare identifiers as strings so uint64/int64 identifiers stay exact.
        with open(req["ids_file"], encoding="utf-8") as handle:
            ids = json.load(handle)
        predicate = ds.field(req["id_column"]).cast(pa.string()).isin(pa.array(ids, type=pa.string()))
    with open(req["metadata"], "w", encoding="utf-8") as handle:
        json.dump(fields, handle)
    if req.get("schema_only"):
        return
    if any(f["kind"] == "unsupported" for f in fields):
        raise TypeError("Unsupported projected Parquet field type")
    for i, f in enumerate(fields):
        if f["kind"] in ("timestamp", "double"):
            Path(req["numeric_dir"], f"{i}.bin").write_bytes(b"")
            Path(req["numeric_dir"], f"{i}.null").write_bytes(b"")
    with open(req["csv"], "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        if req.get("distinct"):
            batches = dataset.to_table(columns=columns, filter=predicate).group_by(columns, use_threads=False).aggregate([]).to_batches(max_chunksize=65536)
        else:
            batches = dataset.scanner(columns=columns, filter=predicate,
                                      batch_size=65536, use_threads=False).to_batches()
        for batch in batches:
            values = []
            for i, (col, f) in enumerate(zip(batch.columns, fields)):
                if f["kind"] == "timestamp":
                    scale = {"s": 1, "ms": 1000, "us": 1000000, "ns": 1000000000}[f["unit"]]
                    numbers = pc.cast(col, pa.int64()).to_numpy(zero_copy_only=False).astype("<f8") / scale
                elif f["kind"] == "double":
                    numbers = pc.cast(col, pa.float64()).to_numpy(zero_copy_only=False).astype("<f8")
                if f["kind"] in ("timestamp", "double"):
                    with open(Path(req["numeric_dir"], f"{i}.bin"), "ab") as binary:
                        binary.write(numbers.tobytes())
                    with open(Path(req["numeric_dir"], f"{i}.null"), "ab") as mask:
                        mask.write(col.is_null().to_numpy(zero_copy_only=False).astype("uint8").tobytes())
                    vals = [""] * len(col)
                elif f["kind"] == "date":
                    vals = [None if v is None else v.isoformat() for v in col.to_pylist()]
                else:
                    vals = col.to_pylist()
                values.append(vals)
            for row in zip(*values):
                if any(isinstance(v, str) and v == req["null"] for v in row):
                    raise ValueError("Null marker collision")
                writer.writerow([req["null"] if v is None else
                                 format(v, ".17g") if isinstance(v, float) else v for v in row])


def write(req):
    fields = json.loads(Path(req["metadata"]).read_text())
    data = {f["name"]: [] for f in fields}
    with open(req["csv"], newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            for f in fields:
                v = row[f["name"]]
                if v == req["null"]:
                    v = None
                elif f["kind"] in ("timestamp", "double"):
                    v = float(v)
                elif f["kind"] == "integer":
                    v = int(v)
                elif f["kind"] == "boolean":
                    v = v.upper() == "TRUE"
                data[f["name"]].append(v)
    for i, f in enumerate(fields):
        if f["kind"] in ("timestamp", "double"):
            numbers = np.fromfile(Path(req["numeric_dir"], f"{i}.bin"), dtype="<f8")
            nulls = np.fromfile(Path(req["numeric_dir"], f"{i}.null"), dtype="uint8")
            if len(numbers) != len(data[f["name"]]) or len(nulls) != len(numbers):
                raise ValueError("Numeric transport row count mismatch")
            data[f["name"]] = [None if missing else v for v, missing in zip(numbers.tolist(), nulls)]
    arrays, names = [], []
    for f in fields:
        vals = data[f["name"]]
        kind = f["kind"]
        if kind == "timestamp":
            vals = [None if v is None else round(v * 1000000) for v in vals]
            a = pa.array(vals, type=pa.int64()).cast(pa.timestamp("us", tz=f["timezone"]))
        elif kind == "date":
            a = pc.cast(pa.array(vals, type=pa.string()), pa.date32())
        else:
            t = dict(integer=pa.int64(), double=pa.float64(), boolean=pa.bool_(), string=pa.string())[kind]
            a = pa.array(vals, type=t)
        arrays.append(a)
        names.append(f["name"])
    pq.write_table(pa.Table.from_arrays(arrays, names=names), req["path"])


if __name__ == "__main__":
    if sys.argv[1:] == ["--probe"]:
        with tempfile.TemporaryDirectory(prefix="clif_parquet_probe_") as directory:
            path = Path(directory, "probe.parquet")
            pq.write_table(pa.table({"probe": [1, None]}), path)
            assert ds.dataset(path, format="parquet").to_table()["probe"].to_pylist() == [1, None]
        print(json.dumps(dict(pyarrow=pa.__version__, parquet=True, dataset=True)))
        sys.exit(0)
    request = json.loads(Path(sys.argv[1]).read_text())
    if request["operation"] == "read":
        read(request)
    elif request["operation"] == "write":
        write(request)
    else:
        raise ValueError("Unknown operation")
