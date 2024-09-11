# Tests

Tests are run with [PyTest](https://docs.pytest.org/).

To run all basic tests run the following from the project repo root:
```
./test/run_tests.sh
```


## test directory

* [aspect_py_project_template_test](aspect_py_project_template_test) - Contains all basic tests for the project.
  * Run with [run_tests.sh](run_tests.sh)
* [jupyter_notebook_test](jupyter_notebook_test) - Tests for Jupyter notebooks.
  * Run with [run_notebook_tests.sh](run_notebook_tests.sh)
* [pytest.ini](pytest.ini) - [PyTest configuration](https://docs.pytest.org/en/stable/reference/customize.html#pytest-ini) file.


### Structure
The `./test/` directory follows the `./src/` directory but with every file and folder name having the suffix `_test` appended. This is because:
* It enables [PyTest test discovery](https://docs.pytest.org/en/latest/explanation/goodpractices.html#test-discovery)
* It does replicate the project directory without resulting in naming conflicts.
* It keeps the `test` directory separated from the `src` directory so tests don't have to be shipped in the resulting Python package.
