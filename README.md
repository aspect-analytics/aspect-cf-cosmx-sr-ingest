# Python Project Template
Template for a new Python project.

To use this template:
* Read through this README.
* Start a new repository from this template by clicking "Use this template" on the [repo's GitHub page](https://github.com/aspect-analytics/python_project_template).
  * !["Use this template" button in green](https://github.com/aspect-analytics/python_project_template/assets/951093/e8d132b4-004c-4826-a9ea-f76b6403aa1b)
* Name your repository source, this will be the name to import (e.g. `import aspect_py_project_template`).
  * Replace `./src/aspect_py_project_template` with your desired package name.
  * Replace all text instances of `aspect_py_project_template` with your desired package name.
* Replace all text instances of `aspect-py-project-template` with your desired package name.
  * This is the name of the package that will be published to the internal PyPI repository. e.g. installable via `pip install aspect-py-project-template`.
* Replace all text instances of `example_env` with your desired environment name for local development.
* Stage the changes to be committed to git, and commit them.
* Add an initial git tag to the commit, this will be the first version of the package.
  * `git tag v0.0.1`
* Push the changes to the remote repository.
  * `git push && git push --tags`

If you want to create a new cloudfunction, look at the [`cloudfunctions_py_example_app`](https://github.com/aspect-analytics/cloudfunctions_py_example_app) project instead as a starting point!



## Structure of the project:
* [`./env`](env/): Python environment requirements and scripts to setup a local development environment.
* [`./infra`](infra/): Docker containers and other project specific infrastructure code.
* [`./notebooks`](notebooks/): Jupyter notebooks (if needed)
* [`./src`](src/): Source code directory
    * [`./src/aspect_py_project_template`](src/aspect_py_project_template/): Package of this repo, change `aspect_py_project_template` to your desired package name.
    * [`./src/scripts`](src/scripts/): Any scripts that should not be part of the installable package itself.
* [`./test/`](test/): Project specific tests
* [`./tools/`](tools/): Tools for code checking, formatting, and testing.

See each subdirectory for its specific documentation.

Please remove any directories and files from this template that are not needed for your project.



## Python Package
This repository defines an example installable Python package with example name `aspect-py-project-template`.

To automatically publish new builds of this package to our internal PyPI repository, Set `pythonPublish: true` in the [`Jenkinsfile`](Jenkinsfile).

The package config is defined in the [pyproject.toml](pyproject.toml) file, which follows the [pyproject.toml specification](https://pip.pypa.io/en/stable/reference/build-system/pyproject-toml/).
The package itself is build in our cloud platform by our CI/CD build system as configured in the [`Jenkinsfile`](Jenkinsfile). The CI/CD pipeline uses the container provided at [infra/Dockerfile](infra/Dockerfile) to build the package using the the [build](https://pypa-build.readthedocs.io/en/stable/#) tool.


### Aspect Internal PyPI repository
If the tests and build passes, and the [`Jenkinsfile`](Jenkinsfile) has `pythonPublish: true` a new version of the package is automatically published to our [internal PyPI repository](https://nexus.aspect.aspect-analytics.ninja/#browse/browse:pypi-aspect-analytics-private) each time when successfully merged to the Git `main` branch. It gets assigned a version according to an environment variable injected by the CI system which should correspond to the [Git tag](https://git-scm.com/book/en/v2/Git-Basics-Tagging).

You can find more information on how to use the internal PyPI repositories in the [Aspect-Analytics - Nexus3 docs](https://docs.aspect.aspect-analytics.ninja/tooling/nexus3/), you might need VPN access to view this.


### Install from Aspect Internal PyPI repository
You can install packages build by the CI/CD system from pip. To do this you will have to setup credentials to the [Aspect Analytics internal PyPI repository](https://nexus.aspect.aspect-analytics.ninja/#browse/browse:pypi-aspect-analytics-private).

You can find the credentials under your Aspect Analytics 1Password, in the "Shared" vault under [`Aspect-Analytics - Nexus3 engineer`](https://start.1password.com/open/i?a=4FVQAO23LNCLVEMFCCQUWRLI3U&v=pzytlb6yoa3e3qjm4kxapuallq&i=kvhhqfi6fbwrf2cv3pezjfsi6a&h=aspect-analytics.1password.com).

You can use these credentials to pip install a package using [`--extra-index-url`](https://pip.pypa.io/en/stable/cli/pip_install/#cmdoption-extra-index-url):
```
pip install --extra-index-url <url_with_credentials> <package_name>
```
Where `<url_with_credentials>` is the URL to the internal PyPI repository with login credentials of [`Aspect-Analytics - Nexus3 engineer`](https://start.1password.com/open/i?a=4FVQAO23LNCLVEMFCCQUWRLI3U&v=pzytlb6yoa3e3qjm4kxapuallq&i=kvhhqfi6fbwrf2cv3pezjfsi6a&h=aspect-analytics.1password.com), and is of the form `https://<account_name>:<password>@nexus.aspect.aspect-analytics.ninja/repository/pypi-aspect-analytics-private/simple`.


#### Configure credentials
Instead of always having to specify the credentials you can configure them by setting up your [~/.config/pip/pip.conf](https://pip.pypa.io/en/stable/topics/configuration/) as detailed below.

If the `~/.config/pip/pip.conf` does not exist yet for you you can create it by running the following commands in your terminal:
```
mkdir -p ~/.config/pip/ && touch ~/.config/pip/pip.conf
```

Then open the file in your favorite text editor. Once in text editor, add the following credentials:
```ini
[global]
extra-index-url = <url_with_credentials>
```

##### Test that credentials have been correctly configured
To install this example package build by the CI/CD system run:
```
pip install aspect-py-project-template
```
You can get all available versions of the package by running: `pip index versions aspect-py-project-template` with the latest version of pip.



## Contribute

### Environment

A Python environment to locally develop on this repo is provided at [`./env`](env/), for more info see the [env Readme](env/README.md). This environment can be setup locally by running:
```
./env/setup_local_conda_env.sh --name example_env
```
If you've chosen a different environment name, make sure to do a search and replace `example_env` to the name of your choosing in `.pre-commit-config.yaml`.


### Tools

Tooling to standardize, check and test the code are provided at [`./tools`](tools/). To run all checks:
```
./tools/run_all.sh
```


### Git pre-commit hooks

[Pre-commit](https://pre-commit.com/) hooks are setup when setting up the local environment from [./env](env/). Once setup, they will run automatically when committing changes to the repository. To run them manually run:
```
pre-commit run --all-files
```

The pre-commit hooks will mainly take care of standardizing the code by running formatters. By default they won't run tests. To ensure that your code works run tests locally before committing. Additionally, automated tests are part of the container build in [./infra](infra/) and will run automatically on the CI/CD system.



## CI/CD
This repository's container is automatically built as is defined in the [`Jenkinsfile`](Jenkinsfile). The build includes running all checks and tests and should fail if any of these fail. The build results can be found in [Jenkins](https://jenkins.aspect.aspect-analytics.ninja/job/github/job/python_project_template/) (this requires VPN access).


### Docker Container
An example Docker image including tests is provided at [`./infra/Dockerfile`](infra/Dockerfile). More info on how to build and run this at [`./infra/README.md`](infra/README.md).


### Versioning

Use [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/) git messages to version the repository automatically.
