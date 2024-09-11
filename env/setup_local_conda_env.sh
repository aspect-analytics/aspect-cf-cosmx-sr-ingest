#!/usr/bin/env bash

set -eo pipefail

# Colors for output formatting
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

# Parse Script Arguments ###########################################
# Help parse the repository root path
# https://stackoverflow.com/a/18443300
realpath() (
    OURPWD="$PWD"
    cd "$(dirname "$1")"
    LINK=$(readlink "$(basename "$1")")
    while [ "$LINK" ]; do
        cd "$(dirname "$LINK")"
        LINK=$(readlink "$(basename "$1")")
    done
    REALPATH="$PWD/$(basename "$1")"
    cd "$OURPWD"
    echo -e "$REALPATH"
)
# Parse repository root path
REPO_ROOT_PATH=$(dirname "$(dirname "$(realpath "$0")")")


# Parse command line arguments
PARAMS=""
ENV_NAME=""
while (( "$#" )); do
  case "$1" in
    -n|--name)
      if [ -n "$2" ] && [ "${2:0:1}" != "-" ]; then
        ENV_NAME=$2
        shift 2
      else
        echo -e "${RED}Error${NC}: Argument for $1 is missing" >&2
        exit 1
      fi
      ;;
    --*|-*) # unsupported flags
      echo -e "${RED}Error${NC}: Unsupported flag $1" >&2
      exit 1
      ;;
    *) # preserve positional arguments
      PARAMS="$PARAMS $1"
      shift
      ;;
  esac
done
# set positional arguments in their proper place
eval set -- "$PARAMS"

# Check environment name
if [ -z "${ENV_NAME}" ]
then
    echo -e "${RED}Error${NC}: A Conda environment name must be given."
    exit 1
fi
echo -e "Trying to create an environment with name '${ENV_NAME}'"


# Setup and check Conda/Mamba ######################################
# Check if mamba is installed, and if not fallback to conda
if command -v mamba > /dev/null; then
    echo -e "Mamba installed, using Mamba to create environment."
    conda(){ command mamba "$@"; }  # Redefine conda to use mamba
elif command -v micromamba > /dev/null; then
    echo -e "MicroMamba installed, using Mamba to create environment."
    conda(){ command micromamba "$@"; } # Redefine conda to use micromamba
elif command -v conda > /dev/null; then
    echo -e "Mamba not installed, falling back to Conda to create environment."
else
    echo -e "${RED}Error${NC}: Mamba nor Conda are installed!. Preferably install Mamba to be able to create the Python environment."
    exit 1
fi


# Cleanup old
if { conda env list | grep "\b${ENV_NAME}\b"; } >/dev/null 2>&1; then
    echo -e ""
    echo -e "${RED}Error${NC}: Environment with name '${ENV_NAME}' already exists!"
    echo -e "Exiting environment setup."
    echo -e "You can remove the old environment using: \`conda deactivate && conda remove --name ${ENV_NAME} --all\`"
    exit 1
fi


# Check if Aspect repo is available ################################
# Check pip.conf
if command -v pip > /dev/null; then
  if [ "$(pip config list | grep -ic "nexus.aspect.aspect-analytics.ninja")" -ge 1 ]
  then
      echo -e "Aspect Analytics Python repository is found in pip config."
  else
      echo -e "${RED}Error${NC}: Aspect Analytics Python repository is not found in pip config!"
      echo -e "pip config:"
      pip config list
      echo -e "Please look at https://github.com/aspect-analytics/python_project_template#configure-credentials to see how to setup the Aspect Python repository."
      exit 1
  fi
else
  echo -e "pip is not installed by itself, skipping pip config check."
  echo -e "Make sure the Aspect Analytics Python repository is setup properly or the following install might fail."
  echo -e "Please look at https://github.com/aspect-analytics/python_project_template#configure-credentials to see how to setup the Aspect Python repository."
fi


# Setup new Python environment #####################################
# Create new env
echo -e "Creating a Mamba/Conda environment named '${ENV_NAME}':"
conda env create \
    --file "${REPO_ROOT_PATH}/env/env_base_conda.yml" \
    --name "${ENV_NAME}"
# Setup dependencies & dev tools & jupyter
echo -e "Running pip install in '${ENV_NAME}':"
conda run --name "${ENV_NAME}" python -m pip install \
    --requirement "${REPO_ROOT_PATH}/env/env_pip_requirements.txt" \
    --requirement "${REPO_ROOT_PATH}/env/env_testing.txt" \
    --requirement "${REPO_ROOT_PATH}/env/env_notebook.txt"

# Install local repo as package (helps with imports)
VERSION="$(git describe --tags --always --dirty=-local)"
export VERSION
echo -e "Installing local repository as package with local version '${VERSION}'"
conda run --name "${ENV_NAME}" python -m \
  pip install --editable "${REPO_ROOT_PATH}"

echo -e "Environment '${ENV_NAME}' created at $(conda run --name "${ENV_NAME}" which python)"

# Setup Git pre-commit
echo -e "Installing Git pre-commit hooks"
conda run --name "${ENV_NAME}" python -m pip install pre-commit
echo -e "Setting up Git pre-commit hooks"
cd "${REPO_ROOT_PATH}"
conda run --name "${ENV_NAME}" pre-commit install

echo -e "${GREEN}Environment '${ENV_NAME}' created successfully.${NC}"
echo -e "To activate the environment run: \`conda activate ${ENV_NAME}\`"
