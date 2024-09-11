# Docker container

A Docker container image is provided for CI/CD and package deployment.


## Build container
Build the container from the local package with:
```
docker buildx build \
    --progress=plain \
    --secret id=PIPCONF,src=${HOME}/.config/pip/pip.conf \
    --tag python-project-template-app \
    --file ./infra/Dockerfile \
    --build-arg="VERSION=$(git describe --tags --always --dirty=-dev)" \
    .
```

This will build the container with the example run script, and will trigger all tests to run.
Please note: the build will fail if the tests won't pass!

### Access internal repository
If you need to access the Aspect Analytics private repository, you can pass the credentials as a secret file in the build command using `--secret id=PIPCONF,src=${HOME}/.config/pip/pip.conf`, each docker build `RUN` command that needs to accees the repository needs to start with `RUN --mount=type=secret,id=PIPCONF,dst=/etc/secrets/pip.conf export PIP_CONFIG_FILE="/etc/secrets/pip.conf"`. This will mount the secret file in the container and set the [`PIP_CONFIG_FILE` environment variable](https://pip.pypa.io/en/stable/topics/configuration/#pip-config-file) to point to the secret file.

Make sure you have locally setup your `~/.config/pip/pip.conf` file by looking at [these instructions](https://github.com/aspect-analytics/python_project_template#install).


## Conda environment

The conda environment is replicated in the container using [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html).


### Activate environment

We don't activate the environment explicitly, but to run a comment in the container you can use the [`conda run`](https://docs.conda.io/projects/conda/en/latest/commands/run.html) equivalent: `micromamba run --prefix $PYTHON_ENV_DIR`. For example see [`./entrypoint.sh`](entrypoint.sh).

Activating a conda environment typically sets a few variables, this can be more than just updating the `$PATH` variable see the [conda documentation for more info](https://docs.conda.io/projects/conda/en/latest/dev-guide/deep-dives/activation.html).

The downside of using `conda run` is that it will still need conda/mamba/micromamba to be present. Ideally our final container is minimal (for efficiency and security reasons). In the future we could think of using something like [conda-pack](https://conda.github.io/conda-pack/) to create a self-contained environment. Alternatively, we could just manually source `$PYTHON_ENV_DIR/etc/conda/activate.d` ([similar to conda-pack](https://github.com/conda/conda-pack/blob/49f6ca9ac8138b8c5b5daec791ccde5700acbd20/conda_pack/scripts/posix/activate#L60)) to initialize the environment.
