#!/bin/sh
# Fresh publication analysis; Make never reuses computed analysis intermediates.
set -eu
exec make -C "$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)" analyze "$@"
