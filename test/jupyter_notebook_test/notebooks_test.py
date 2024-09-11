"""
Test the notebooks using [nbconvert](https://nbconvert.readthedocs.io/en/latest/).
"""
import logging
import tempfile
from pathlib import Path
from typing import Iterator

import nbformat
import pytest
from nbconvert.preprocessors import CellExecutionError, ExecutePreprocessor
from nbformat import NotebookNode

_log = logging.Logger(__name__)


TWO_MINUTES = 120
JUPYTER_VERSION = 4


def run_notebook(notebook_path: Path) -> NotebookNode:
    """
    Run a notebook and return the executed notebook object.
    """
    _log.info(f"Runnin notebook tests: '{notebook_path!s}'.")
    with notebook_path.open("r") as f:
        nb = nbformat.read(f, as_version=JUPYTER_VERSION)
    ep = ExecutePreprocessor(timeout=TWO_MINUTES)
    try:
        nb, _ = ep.preprocess(nb)
        return nb
    except CellExecutionError as exc:
        _log.error(f"Error executing the notebook: {notebook_path}")
        raise exc


@pytest.fixture
def failure_notebook() -> Iterator[Path]:
    with tempfile.TemporaryDirectory() as tmp_dir:
        nb_path = Path(tmp_dir) / "failure.ipynb"
        nb = nbformat.v4.new_notebook()
        nb["cells"] = [
            nbformat.v4.new_code_cell(
                "raise RuntimeError('This is an expected failure.')"
            )
        ]
        with nb_path.open("w") as f:
            nbformat.write(nb, f)
        yield nb_path


@pytest.mark.notebook_test
def test_failure_detected(failure_notebook: Path):
    """Test that detecting failures works as expected."""
    with pytest.raises(CellExecutionError):
        try:
            run_notebook(failure_notebook)
        except CellExecutionError as exc:
            assert exc.ename == "RuntimeError"
            assert exc.evalue == "This is an expected failure."
            _log.info(f"Expected failure: {exc.ename} - {exc.evalue}")
            raise exc


@pytest.mark.notebook_test
def test_notebooks():
    """
    Iterates over the notebooks in the project's notebook folder and verifies they run with no errors
    """
    notebooks_dir = Path(__file__).parent.parent.parent / "notebooks"
    for notebook in notebooks_dir.rglob("*.ipynb"):
        _log.info(f"Testing notebook: {notebook}")
        nb = run_notebook(notebook)
        nb_cells = nb["cells"]
        for cell in nb_cells:
            if "outputs" in cell:
                assert "error" not in cell["outputs"], _log.error(
                    f"Notebook failed: {notebook} with output {cell['outputs']}"
                )
