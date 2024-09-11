# syntax=docker/dockerfile:1.4
#
# Dockerfile containing the conda enviroment + internal library installs + pipeline scripts.
# Build from root directory of this repository (which should be one level up from this file.)


# Directory for the Python environment
ARG PYTHON_ENV_DIR=/python_env/
ARG BASE_IMG=ubuntu:22.04
# Package version, this is used to tag the package
ARG VERSION


# Container with micromamba, updated ubuntu
##############################################################################
FROM ${BASE_IMG} as mamba-base

ARG DEBIAN_FRONTEND=noninteractive
ARG PYTHONDONTWRITEBYTECODE=1

ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

RUN cd /tmp/ \
    # Update ubuntu
    && apt-get update \
    && apt-get --yes upgrade \
    # Install packages and build dependencies
    && build_deps="\
        bzip2 \
        curl" \
    && apt-get update \
    && apt-get install --quiet --yes --no-install-recommends \
        # General requirements (tar cannot be uninstalled so it's not in the build_deps)
        ca-certificates \
        tar \
        # Build dependencies (temporarily installed)
        ${build_deps} \
    # Install mamba (micromamba)
    # https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html
    && curl -sL https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj --directory /usr/ bin/micromamba \
    # Link micromamba to mamba
    && ln -s /usr/bin/micromamba /usr/bin/mamba \
    # Cleanup
    && mamba clean --all --force-pkgs-dirs --yes \
    # apt-get purge also removes config files (which remove doesn't)
    && apt-get purge --yes ${build_deps} \
    && apt-get autoremove --purge --yes \
    && apt-get clean \
    && rm --recursive --force ~/.cache/pip/ /tmp/* /var/tmp/* /var/lib/apt/lists/*


# Container with base Python environment (not including the package)
##############################################################################
FROM mamba-base as python-base

# Temporary build ARGs, not permanent ENVs
ARG DEBIAN_FRONTEND=noninteractive
ARG PYTHONDONTWRITEBYTECODE=1

ARG PYTHON_ENV_DIR
ENV PYTHON_ENV_DIR=${PYTHON_ENV_DIR}

# Install Python environment
COPY ./env/env_base_conda.yml /tmp/conda_export/env_base_conda.yml
# Mount pip.conf to be able to install from internal repo
# `PIP_CONFIG_FILE` defines where pip looks for the config file
RUN --mount=type=secret,id=PIPCONF,dst=/etc/secrets/pip.conf export PIP_CONFIG_FILE="/etc/secrets/pip.conf" \
    && cd /tmp/ \
    # Setup Python environment
    && mamba env create \
        --file /tmp/conda_export/env_base_conda.yml \
        --prefix ${PYTHON_ENV_DIR} \
    # Cleanup
    && mamba clean --all --force-pkgs-dirs --yes \
    && rm --recursive --force ~/.cache/pip/ /tmp/* /var/tmp/*


# Container with package installed
##############################################################################
FROM python-base as python-package-env

# Temporary build ARGs, not permanent ENVs
ARG DEBIAN_FRONTEND=noninteractive
ARG PYTHONDONTWRITEBYTECODE=1

ARG VERSION

# Install repo package in Python environment
COPY ./ /tmp/repo_clone
# Mount pip.conf to be able to install from internal repo
# `PIP_CONFIG_FILE` defines where pip looks for the config file
RUN --mount=type=secret,id=PIPCONF,dst=/etc/secrets/pip.conf export PIP_CONFIG_FILE="/etc/secrets/pip.conf" \
    && echo "Pip config file found at: \"$(mamba run --prefix ${PYTHON_ENV_DIR} pip config get ':env:.config-file')\"" \
    # Ensure aspect-analytics repo is found
    && if [ $(mamba run --prefix ${PYTHON_ENV_DIR} pip config list | grep -ic "aspect-analytics") -ge 1 ]; \
        then echo "Aspect Analytics Python repository is found in pip config."; \
        else echo "Aspect Analytics Python repository is not found in pip config!"; \
        mamba run --prefix ${PYTHON_ENV_DIR} pip config debug; \
        exit 1; fi \
    && cd /tmp/ \
    # Setup Python environment by installing the local package
    && mamba run --prefix ${PYTHON_ENV_DIR} python -m pip install \
        /tmp/repo_clone \
    # Env cleanup to save space
    && find ${PYTHON_ENV_DIR}/ -follow -type f -name '*.a' -delete \
    && find ${PYTHON_ENV_DIR}/ -follow -type d -name '__pycache__' -exec rm -rf {} + -depth \
    # Cleanup
    && rm --recursive --force ~/.cache/pip/ /tmp/* /var/tmp/*


# Container with tests
##############################################################################
FROM python-package-env as testing

# Temporary build ARGs, not permanent ENVs
ARG DEBIAN_FRONTEND=noninteractive
ARG PYTHONDONTWRITEBYTECODE=1

COPY ./ /tmp/repo_clone
# Mount pip.conf to be able to install from internal repo
# `PIP_CONFIG_FILE` defines where pip looks for the config file
RUN --mount=type=secret,id=PIPCONF,dst=/etc/secrets/pip.conf export PIP_CONFIG_FILE="/etc/secrets/pip.conf" \
    && cd /tmp/repo_clone \
    # Setup Python testing environment
    && mamba run --prefix ${PYTHON_ENV_DIR} python -m pip install \
        --requirement ./env/env_testing.txt \
    # Run check and tests
    && mamba run --prefix ${PYTHON_ENV_DIR} ./tools/run_linters.sh \
    && mamba run --prefix ${PYTHON_ENV_DIR} ./test/run_tests.sh \
    # Write file iff tests are ok.
    # This file is used in a next build stage to assert tests pass.
    && touch /test.ok
    # Cleanup not needed (Nothing builds on this layer)


# Container to build package, this should be the last build stage
##############################################################################
FROM python-base as pkg_builder

# Temporary build ARGs, not permanent ENVs
ARG DEBIAN_FRONTEND=noninteractive
ARG PYTHONDONTWRITEBYTECODE=1
ARG VERSION

# Ensure tests pass
# Trigger testing continer by ensuring this output (empty file written when tests are ok)
COPY --from=testing --chmod=555 /test.ok /test.ok

# Make sure `python` refers to our installed Python environment
ENV PATH=${PYTHON_ENV_DIR}/bin:$PATH

# Copy repository to container and set as working directory
ENV REPO_PATH=/app/repository
COPY ./ ${REPO_PATH}
WORKDIR ${REPO_PATH}

# Build package (but don't push to package repository, pushing is done by the CI system)
RUN --mount=type=secret,id=PIPCONF,dst=/etc/secrets/pip.conf export PIP_CONFIG_FILE="/etc/secrets/pip.conf" \
    && mamba run --prefix ${PYTHON_ENV_DIR} \
        python -m pip install --no-compile --no-cache-dir --upgrade \
            build \
            twine \
    && mamba run --prefix ${PYTHON_ENV_DIR} python -m build \
    && [ -d "${REPO_PATH}/dist" ] && echo "Directory ${REPO_PATH}/dist exists." \
    && ls -al ${REPO_PATH}/dist
