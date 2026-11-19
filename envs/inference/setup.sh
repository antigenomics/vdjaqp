#!/bin/sh
# Isolated pyABC + existing source dependencies; leave baseline runtime unchanged.
set -eu
repo=$(CDPATH= cd -- "$(dirname "$0")/../.." && pwd)
uv venv --python "$repo/.venv/bin/python" "$repo/.tool-envs/inference"
uv pip install --python "$repo/.tool-envs/inference/bin/python" -r "$repo/envs/inference/requirements.lock" -e "$repo"
uv pip freeze --python "$repo/.tool-envs/inference/bin/python" > "$repo/.tool-envs/inference/installed.txt"
