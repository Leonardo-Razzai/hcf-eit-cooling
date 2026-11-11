"""Run one line of a scan list (used by run_array.sbatch).

A scan list has one run per line: a study script in simulation/studies and its arguments, e.g.
    zeeman2d.py run --scheme te01_bz --case A --s 1.0 --d1 0.4 --d2 0.4 --gamma 0.0 --N 6
Blank lines and lines starting with '#' are ignored; --line counts the remaining lines from 1.
Extra arguments after the list are appended to the command (e.g. --tmax 300 for a test).

    python simulation/cluster/run_list.py simulation/cluster/scans/zeeman2d_test.txt --line 1
"""
import argparse
import shlex
import subprocess
import sys
from pathlib import Path

STUDIES = Path(__file__).resolve().parents[1] / 'studies'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('list')
    ap.add_argument('--line', type=int, required=True)
    ap.add_argument('--dry', action='store_true', help='print the command only')
    a, extra = ap.parse_known_args()
    lines = [l.strip() for l in open(a.list) if l.strip() and not l.lstrip().startswith('#')]
    if not 1 <= a.line <= len(lines):
        sys.exit('line %d out of range (1-%d)' % (a.line, len(lines)))
    cmd = [sys.executable] + shlex.split(lines[a.line - 1]) + extra
    print('run:', ' '.join(cmd), flush=True)
    if a.dry:
        return
    sys.exit(subprocess.call(cmd, cwd=STUDIES))


if __name__ == '__main__':
    main()
