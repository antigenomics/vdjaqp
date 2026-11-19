#!/bin/sh
# Build published SAMtools/HTSlib 1.25 locally; baseline .venv and compiler required.
set -eu
cd "$(dirname "$0")/../.."
threads=$(.venv/bin/python -c 'import sys; from vdjaqp import resolve_cpu_budget; print(resolve_cpu_budget(int(sys.argv[1]) if sys.argv[1] else None))' "${1:-}")
target="$PWD/.tool-envs/samtools"
if [ -e "$target" ]; then echo 'SAMtools target exists; refusing overwrite' >&2; exit 1; fi
mkdir -p "$target/bin"
curl -fL --max-time 60 https://github.com/samtools/samtools/releases/download/1.25/samtools-1.25.tar.bz2 -o "$target/source.tar.bz2"
.venv/bin/python - "$target/source.tar.bz2" <<'PY'
import hashlib,sys
from pathlib import Path
assert hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest()=='8bbd2b3cb447953149a74ab0b00e3ec4e88a789a27b8373a31fa6b5b2dee02b1','release archive hash mismatch'
PY
tar -xjf "$target/source.tar.bz2" -C "$target"
(cd "$target/samtools-1.25" && ./configure --without-curses --disable-libcurl --disable-s3 --disable-gcs) > "$target/configure.log" 2>&1
.venv/bin/python - "$target" "$threads" <<'PY'
import os,signal,subprocess,sys
from pathlib import Path
target=Path(sys.argv[1])
with (target/'build.log').open('w') as log:
 process=subprocess.Popen(['make','-j'+sys.argv[2],'samtools'],cwd=target/'samtools-1.25',stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
 try: result=process.wait(timeout=300)
 except subprocess.TimeoutExpired:
  os.killpg(process.pid,signal.SIGKILL);process.wait();raise
 if result:raise subprocess.CalledProcessError(result,process.args)
PY
cp "$target/samtools-1.25/samtools" "$target/bin/samtools"
"$target/bin/samtools" --version
make -C "$target/samtools-1.25" clean > "$target/clean.log" 2>&1
