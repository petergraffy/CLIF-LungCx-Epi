"""Official CLIFpy waterfall with encounter-specific input/output boundaries."""
import os
from pathlib import Path
import json
import inspect
import hashlib
import importlib.metadata
import warnings
import pandas as pd
from clifpy.utils.waterfall import process_resp_support_waterfall

TEXT = ["device_category", "device_name", "mode_category", "mode_name"]
NUMERIC = ["tracheostomy", "fio2_set", "lpm_set", "peep_set", "tidal_volume_set",
           "resp_rate_set", "resp_rate_obs", "pressure_support_set", "peak_inspiratory_pressure_set"]

def bounded_waterfall(raw, bounds):
    bounds = bounds.copy()
    for col in ("start", "end"):
        bounds[col] = pd.to_datetime(bounds[col], utc=True)
    assert not bounds.hospitalization_id.duplicated().any()
    raw = raw.copy()
    raw["recorded_dttm"] = pd.to_datetime(raw.recorded_dttm, utc=True)
    scoped = raw.merge(bounds[["hospitalization_id", "start", "end"]], on="hospitalization_id", how="inner")
    scoped = scoped.loc[(scoped.recorded_dttm >= scoped.start) & (scoped.recorded_dttm < scoped.end)]
    scoped = scoped.drop(columns=["start", "end"])
    for col in TEXT:
        if col not in scoped:
            scoped[col] = None
    for col in NUMERIC:
        if col not in scoped:
            scoped[col] = float("nan")
    original = scoped.copy()
    # Empty boundary markers extend the scaffold, not the source observations.
    # Never duplicate a real source timestamp: CLIFpy's duplicate filter would
    # otherwise discard a device-missing row containing real numeric values.
    markers = pd.concat([
        bounds[["hospitalization_id", "start"]].rename(columns={"start": "recorded_dttm"}),
        bounds.assign(last=bounds.end - pd.Timedelta(microseconds=1))[["hospitalization_id", "last"]]
        .rename(columns={"last": "recorded_dttm"})], ignore_index=True)
    keys = pd.MultiIndex.from_frame(scoped[["hospitalization_id", "recorded_dttm"]])
    markers = markers.loc[~pd.MultiIndex.from_frame(markers).isin(keys)]
    data = pd.concat([scoped, markers], ignore_index=True)
    for col in TEXT:
        # The reference function's change-ID helper expects object strings:
        # nullable StringDtype propagates NA through its first comparison.
        data[col] = data[col].astype(object).where(data[col].notna(), None)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = process_resp_support_waterfall(data, bfill=False, verbose=False)
    # Scaffold rows may extend to the end of a clock hour; trim again to the
    # exact patient-specific half-open interval after processing.
    result = result.merge(bounds[["hospitalization_id", "start", "end"]], on="hospitalization_id", how="left")
    result = result.loc[(result.recorded_dttm >= result.start) & (result.recorded_dttm < result.end)].copy()
    assert result.start.notna().all() and result.end.notna().all()
    assert not result.duplicated(["hospitalization_id", "recorded_dttm"]).any()
    result = result.drop(columns=["start", "end"])
    observed = original.dropna(subset=["device_category"]).copy()
    observed["device_category"] = observed.device_category.str.lower().str.strip()
    categories = observed.groupby(["hospitalization_id", "recorded_dttm"])["device_category"].agg(set).to_dict()
    result["raw_device_observed_at_timestamp"] = [
        pd.notna(cat) and str(cat).lower().strip() in categories.get((hid, time), set())
        for hid, time, cat in zip(result.hospitalization_id, result.recorded_dttm, result.device_category)]
    fio = original.dropna(subset=["fio2_set"]).groupby(["hospitalization_id", "recorded_dttm"])["fio2_set"].agg(set).to_dict()
    result["raw_fio2_observed_at_timestamp"] = [
        pd.notna(value) and value in fio.get((hid, time), set())
        for hid, time, value in zip(result.hospitalization_id, result.recorded_dttm, result.fio2_set)]
    return result, sorted(set(str(w.message) for w in caught)), len(original)

def main():
    out = Path(os.environ.get("CLIF_WATERFALL_OUTPUT", "output/run_UCMC_waterfall_20261022"))
    raw = pd.read_parquet(out / "local/raw_respiratory.parquet")
    bounds = pd.read_parquet(out / "local/window_bounds.parquet")
    provenance = {"clifpy_version": importlib.metadata.version("clifpy"), "numeric_bfill": False,
                  "categorical_within_window_lookahead": True,
                  "waterfall_source_sha256": hashlib.sha256(Path(inspect.getfile(process_resp_support_waterfall)).read_bytes()).hexdigest(),
                  "pandas_version": pd.__version__, "jobs": []}
    for job, window in bounds.groupby("job", sort=False):
        print("Processing", job, "hospitalizations", len(window), flush=True)
        result, notices, n_raw = bounded_waterfall(raw, window)
        result.to_parquet(out / "local" / f"{job}_waterfall.parquet", index=False)
        provenance["jobs"].append({"job": job, "hospitalizations": len(window), "raw_rows_in_window": n_raw,
                                   "processed_rows": len(result), "warnings": notices})
        print("Completed", job, "processed rows", len(result), flush=True)
    (out / "local/waterfall_manifest.json").write_text(json.dumps(provenance, indent=2))

if __name__ == "__main__":
    main()
