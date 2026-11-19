#!/bin/sh
# Published-tool runtime, isolated from the baseline scoring environment.
set -eu
[ -x .tool-envs/generation/bin/python ] || uv venv --python 3.13.14 .tool-envs/generation
uv pip sync --require-hashes --python .tool-envs/generation/bin/python envs/generation/requirements.txt
