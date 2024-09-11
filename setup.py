import os
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).parent / "src"))

from aspect_py_dev_tools.version_parser import correct_version_safe
from setuptools import setup

VERSION_ENV_VAR = "VERSION"


# Set version number for package build based on environment variable.
# https://stackoverflow.com/a/77001804
version = os.environ.get(VERSION_ENV_VAR)
if version is None:
    print(f"{VERSION_ENV_VAR!r} environment variable is not set!")
    version = "version-not-set"
# Correct version string to be PEP440 compliant.
# Set default version if version is not set.
version = correct_version_safe(version)
print(f"Using version: {version!r}.")

setup(version=version)
