# Python Environment

This directory contains files to replicate the Python environment for this project.


## Local Development Conda Environment

The provided conda environment can be created locally for development if you have [Mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html), which is a [conda](https://docs.conda.io/en/latest/) replacement, installed by running the following from the project root:

```
./env/setup_local_conda_env.sh --name example_env
conda activate example_env
```


And can be cleaned up afterwards with:
```
conda deactivate && conda remove --name example_env --all
```

### Local development install
The local development environment provided at [env/](env/) will install this Python package via `pip install --editable .`, which means it completely references to the root directory of this repository and is editable (re-importing might be required after an update). This helps with imports and makes sure local package imports work irregardless of the directory you're in.


## Package versions
The package versions are defined in the following files:
* [`env_base_conda.yml`](env_base_conda.yml) - Base conda environment. Defines the Python and Pip version used, as well as dependencies that cannot be installed via pip.
* [env_pip_requirements.txt](env_pip_requirements.txt) - Necessary dependencies for the environment that can be installed via pip/PyPI packages. Try to put all dependencies in this file, and only use conda for packages that are not available on PyPI.
* [env_testing.txt](env_testing.txt) - Development tools that are not necessary for the environment, but are useful for development. E.g. Testing, formatting & checking.
* [env_notebook.txt](env_notebook.txt) - Jupyter notebook and visualization dependencies used for debugging and analysis, not core to the environment.

We would recommend bounding the Python package versions for [reproducibility](https://pip.pypa.io/en/stable/topics/repeatable-installs/). We use the [Requirements File Format](https://pip.pypa.io/en/stable/reference/requirements-file-format/#requirements-file-format) to define dependencies and their versions.




## FAQ

### Solving environment: Found conflicts! Looking for incompatible packages.
Reduce the version bounds on your package versions in the environment files. E.g. if you have `numpy==1.19.5` and `pandas==1.2.4` in your environment file, try changing it to `numpy>=1.19.5,<2` and `pandas>=1.2.4,<2`.

### Updating packages
To list outdated packages you can use `pip list --outdated` and/or `conda update --all --dry-run --channel conda-forge` in the activated environment.

### Mamba installation issues (environment cannot be activated)
Installing Mamba (or conda) will modify your `.bashrc` or `.zshrc` file. In case of any issues:
Please verify that:
1. Mamba is installed correctly by running `mamba --version`.
2. Remove anything that mamba/conda added to your `.bashrc` or `.zshrc` file and replace it with the following (assuming you installed MambaForge in your home directory):
```
# Mamba/Conda
#######################################
source ~/mambaforge/etc/profile.d/conda.sh
source ~/mambaforge/etc/profile.d/mamba.sh
```

### Why conda (instead of x/y/z)?
For complex environments we have a preference towards using conda/mamba because:
* conda has less issues with dependency conflicts
* conda can handle non-python dependencies for Python libraries
If you are new to conda you can learn about the differences between conda and venv at ["Conda: Myths and Misconceptions"](https://jakevdp.github.io/blog/2016/08/25/conda-myths-and-misconceptions/).

When using conda make sure to avoid the default `anaconda` channel. Preferably use the `conda-forge` channel to avoid any licencing issues. E.g. from the [Anaconda.com ToS](https://conda-forge.org/blog/posts/2020-11-20-anaconda-tos/):
> What makes it non-surprising is that, at the moment, any third party channel like conda-forge is free. The TOS change does not apply to conda-forge, nor to other channels hosted on anaconda.org; the TOS change in question applies only to the “defaults” channel and other software hosted on repo.anaconda.com.
