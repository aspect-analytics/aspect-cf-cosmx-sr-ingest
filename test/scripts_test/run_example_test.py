"""
Test Any scripts in the `src/scripts/` directory.
"""
import subprocess
from pathlib import Path

import pytest


@pytest.fixture
def run_script() -> Path:
    project_root = Path(__file__).resolve().parent.parent.parent
    script = project_root / "src" / "scripts" / "run_example.py"
    assert script.exists() and script.is_file()
    return script


def test_script_execution(run_script: Path):
    process = subprocess.Popen(
        ["python", "-O", str(run_script), "--molecule", "H2O"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert process.stdout
    stdout = process.stdout.read().decode("utf-8").strip()
    assert stdout == "The mass number of H2O is 18"
    assert process.stderr
    stderr = process.stderr.read().decode("utf-8").strip()
    assert stderr == ""
    process.stdout.close()
    process.stderr.close()
    process.terminate()
    process.wait(timeout=1)
    process.kill()
