#!/usr/bin/env bash
# Run with `./tools/run_checks.sh`

ROOT_PATH=$(dirname "$(dirname "$(readlink -f "$0")")")  # Root directory path
cd "${ROOT_PATH}"

set -eu

echo -e "\nRun ruff code formatting checks"
ruff check \
    "${ROOT_PATH}/src" \
    "${ROOT_PATH}/test"

echo -e "\nRun type hinting checks"
echo "Run mypy checks for ./src/ and ./test/"
mypy \
    --install-types \
    --non-interactive \
    "${ROOT_PATH}/src/" \
    "${ROOT_PATH}/test/"
