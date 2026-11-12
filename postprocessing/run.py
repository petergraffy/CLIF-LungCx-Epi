"""Pool existing Stage 2 aggregate returns and build the coordinator displays."""
import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CODE = REPO / "postprocessing"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="Fresh run directory under output/")
    parser.add_argument("returns", nargs="+", type=Path, help="Completed Stage 2 return directories")
    parser.add_argument("--bundle", type=Path, default=REPO / "frozen/stage2_v1")
    args = parser.parse_args()
    output = args.output.resolve()
    bundle = args.bundle.resolve()
    returns = [path.resolve() for path in args.returns]
    if output == REPO / "output" or not output.is_relative_to(REPO / "output"):
        parser.error("Choose a fresh run directory under the ignored output/ directory.")
    if output.exists():
        parser.error(f"Output already exists: {output}")
    for path in returns:
        if not (path / "return_manifest.json").is_file():
            parser.error(f"No return manifest in {path}")
    output.mkdir(parents=True)
    env = os.environ.copy()
    env.setdefault("MPLCONFIGDIR", str(output / ".matplotlib"))

    def step(name, interpreter, script, *values):
        log_path = output / f"{name}.log"
        command = [*interpreter, str(REPO / script), *map(str, values)]
        print(f"{name} → {log_path}", flush=True)
        with log_path.open("w") as log:
            status = subprocess.run(command, cwd=REPO, env=env, stdout=log, stderr=subprocess.STDOUT)
        if status.returncode:
            raise SystemExit(f"{name} failed. See {log_path}")

    r = ["Rscript", "--vanilla"]
    py = [sys.executable]
    figures = output / "figures"
    survival = output / "survival"
    prediction = output / "prediction"
    reconstruction = survival / "reconstruction"
    reports = survival / "reports"

    step("associations", r, CODE / "analyze.R", bundle, output, *returns)
    with (output / "packet_audit.csv").open(newline="") as stream:
        audit = list(csv.DictReader(stream))
    sites = [row["site"] for row in audit]
    step("descriptives", r, CODE / "collect_descriptives.R", output / "publication_tables", *returns)
    step("profiles", py, CODE / "pool_profiles.py", output, bundle, *returns)
    step("survival_reconstruction", r, "diagnostics/stage2/central_survival.R",
         reconstruction, *(path / "survival" for path in returns))
    step("survival_pooling", r, "diagnostics/pooled_survival/central_pool.R",
         reports, reconstruction / "grid", *(reconstruction / "sites" / site for site in sites))
    manifest_path = reports / "pooling_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["bootstrap"] = (
        "200 independent patient-cluster resamples within each site; all admissions carried together; "
        "pointwise percentile intervals; conditional on the participating sites"
    )
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    step("prediction", py, "diagnostics/stage2/central_prediction.py",
         prediction, *(path / "prediction" for path in returns))
    for script in ("figure_1_trajectories", "figure_2_competing_risks",
                   "figure_3_discrimination", "publication_tables"):
        step(script, py, CODE / f"{script}.py", output)
    step("cohort_flow", r, CODE / "cohort_flow.R", output, figures, *returns)
    step("calibration", py, CODE / "calibration.py", output, *returns)
    step("hospital_survival", py, CODE / "hospital_survival.py", reports)
    step("early_landmarks", py, CODE / "early_landmarks.py", bundle, figures)

    result = {
        "status": "complete",
        "sites": sites,
        "patient_level_inputs_read": False,
        "operational_contract_exceptions": [
            row["site"] for row in audit if row["issued_contract"] != "TRUE"
        ],
        "primary_consensus_id": manifest["consensus_id"],
    }
    (output / "postprocessing_manifest.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"Postprocessing complete → {output}")


if __name__ == "__main__":
    main()
