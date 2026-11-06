"""Build a verified, isolated copy of the released Stage 2 analysis.

Explicit overrides change Parquet I/O and the environment contract only.
The original checkout, bundle, native workspaces and CLIF files stay untouched.
"""
import hashlib
import json
import shutil
import sys
from pathlib import Path


def md5(path):
    return hashlib.md5(Path(path).read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def bundle_manifest(parent, release):
    manifest = json.loads(json.dumps(parent))
    manifest["source_hash"] = release["effective_source_hash"]
    manifest["workspace_id"] = "stage2_python_v1_independent"
    manifest["python_fallback"] = dict(parent_bundle_md5=release["parent_bundle_md5"],
                                       parent_source_hash=release["parent_source_hash"],
                                       r_arrow_required=False, parquet_backend="pyarrow")
    return json.dumps(manifest, indent=2) + "\n"


def build(repo, config, run_root):
    repo, config, run_root = map(lambda p: Path(p).resolve(), (repo, config, run_root))
    release = json.loads((repo / "python_fallback/release.json").read_text())
    for rel, expected in release["distribution_files_md5"].items():
        if md5(repo / rel) != expected:
            raise RuntimeError(f"Fallback release file differs: {rel}; use the issued checkout.")
    for rel, expected in release["parent_source_files_md5"].items():
        if md5(repo / rel) != expected:
            raise RuntimeError(f"Parent source differs: {rel}; restore the coordinator-issued file.")
    for rel, expected in release["parent_aux_files_md5"].items():
        if md5(repo / rel) != expected:
            raise RuntimeError(f"Parent configuration template differs: {rel}")
    parent = repo / "frozen/stage2_v1"
    if hashlib.md5((parent / "bundle_manifest.json").read_bytes()).hexdigest() != release["parent_bundle_md5"]:
        raise RuntimeError("The parent frozen bundle differs from the fallback release.")
    cfg = json.loads(config.read_text())
    # Preserve executable symlinks: resolving the executable itself loses venv identity.
    for key in ("repo", "tables_path", "python"):
        p = Path(cfg[key]).expanduser()
        cfg[key] = str(p.absolute()) if key == "python" else str(p.resolve())
    cfg["analysis"]["run_id"] = "stage2_python_v1"
    cfg.pop("stage1_root", None)
    request = dict(config=json.loads(json.dumps(cfg)), release_md5=md5(repo / "python_fallback/release.json"))
    state = run_root / "local/build_state.json"
    workspace = run_root / "local/pipeline"
    if state.exists():
        if json.loads(state.read_text()) != request:
            raise RuntimeError("Fallback config or release changed; preserve this workspace and use a fresh CLIF_FALLBACK_OUTPUT directory.")
        for rel, expected in release["effective_source_files_md5"].items():
            if md5(workspace / rel) != expected:
                raise RuntimeError(f"Staged fallback source differs: {rel}")
        for rel, expected in release["parent_aux_files_md5"].items():
            if md5(workspace / rel) != expected:
                raise RuntimeError(f"Staged configuration template differs: {rel}")
        expected_config = json.loads(json.dumps(cfg))
        expected_config["repo"] = str(workspace)
        if json.loads((workspace / "config/config.json").read_text()) != expected_config:
            raise RuntimeError("Staged site config changed; edit the original config only and use a fresh fallback workspace.")
        if json.loads((workspace / "fallback/source_inventory.json").read_text()) != list(release["effective_source_files_md5"]):
            raise RuntimeError("Staged source inventory changed.")
        if hashlib.md5((workspace / "frozen/stage2_python_v1/bundle_manifest.json").read_bytes()).hexdigest() != release["fallback_bundle_md5"]:
            raise RuntimeError("The staged fallback bundle manifest differs.")
        return workspace
    if workspace.exists():
        # A failed build is never used as analysis input. Do not erase anything.
        raise RuntimeError("Incomplete pipeline build; preserve it and choose a fresh CLIF_FALLBACK_OUTPUT directory.")
    workspace.mkdir(parents=True)
    for rel in release["parent_source_files_md5"]:
        destination = workspace / rel
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo / rel, destination)
    (workspace / "config").mkdir(exist_ok=True)
    for rel in release["parent_aux_files_md5"]:
        (workspace / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo / rel, workspace / rel)
    for source in (repo / "python_fallback/overrides").rglob("*"):
        if source.is_file():
            destination = workspace / source.relative_to(repo / "python_fallback/overrides")
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
    (workspace / "fallback").mkdir()
    for name in ("parquet.R", "parquet_bridge.py"):
        shutil.copy2(repo / "python_fallback" / name, workspace / "fallback" / name)
    inventory = list(release["effective_source_files_md5"])
    (workspace / "fallback/source_inventory.json").write_text(json.dumps(inventory))
    for rel, expected in release["effective_source_files_md5"].items():
        if md5(workspace / rel) != expected:
            raise RuntimeError(f"Built source failed verification: {rel}")
    bundle = workspace / "frozen/stage2_python_v1"
    shutil.copytree(parent, bundle)
    manifest = bundle_manifest(json.loads((bundle / "bundle_manifest.json").read_text()), release)
    if hashlib.md5(manifest.encode()).hexdigest() != release["fallback_bundle_md5"]:
        raise RuntimeError("Fallback bundle contract differs from its release.")
    (bundle / "bundle_manifest.json").write_text(manifest)
    cfg["repo"] = str(workspace)
    (workspace / "config/config.json").write_text(json.dumps(cfg, indent=2))
    state.write_text(json.dumps(request, indent=2))
    return workspace


if __name__ == "__main__":
    if sys.argv[1:] == ["--manifest-md5"]:
        release = json.loads(Path("python_fallback/release.json").read_text())
        parent = json.loads(Path("frozen/stage2_v1/bundle_manifest.json").read_text())
        print(hashlib.md5(bundle_manifest(parent, release).encode()).hexdigest())
    else:
        print(build(*sys.argv[1:]))
