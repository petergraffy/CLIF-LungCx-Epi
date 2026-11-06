"""Exercise build/resume and tamper checks without opening clinical tables."""
import importlib.util
import json
import sys
import tempfile
from pathlib import Path

repo = Path.cwd()
spec = importlib.util.spec_from_file_location("fallback_build", repo / "python_fallback/build_pipeline.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
(repo / "output").mkdir(exist_ok=True)
with tempfile.TemporaryDirectory(prefix="run_fallback_build_test_", dir=repo / "output") as directory:
    config = Path(directory) / "site.json"
    cfg = json.loads((repo / "config/config_template.json").read_text())
    cfg.update(site_name="hopkins", repo=str(repo), tables_path=str(Path(directory) / "unused_source_tables"), python=sys.executable)
    config.write_text(json.dumps(cfg))
    workspace = module.build(repo, config, Path(directory) / "run with spaces")
    assert module.build(repo, config, Path(directory) / "run with spaces") == workspace
    cases = [workspace / "config/config.json", workspace / "config/config_template.json",
             workspace / "fallback/source_inventory.json", workspace / "code/lib/common.R",
             workspace / "frozen/stage2_python_v1/bundle_manifest.json"]
    for path in cases:
        original = path.read_bytes()
        if path.name == "source_inventory.json":
            path.write_text("[]")
        elif path.name.endswith(".json"):
            value = json.loads(original)
            value["unexpected_test_change"] = True
            path.write_text(json.dumps(value))
        else:
            path.write_bytes(original + b"\n# Unexpected source change\n")
        try:
            module.build(repo, config, Path(directory) / "run with spaces")
        except RuntimeError:
            pass
        else:
            raise AssertionError(f"Tampered input was accepted: {path.name}")
        finally:
            path.write_bytes(original)
print("PASS: isolated build, unchanged resume, paths with spaces, and rejection of changed config/template/inventory/source/bundle")
