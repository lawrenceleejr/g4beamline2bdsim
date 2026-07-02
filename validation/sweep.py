#!/usr/bin/env python3
"""Sweep a set of G4beamline files: convert each, check that BDSIM accepts the
output, optionally run G4beamline, and classify the level of closure.

For every input it reports:
  * convert      - did the converter run (and how many elements / warnings)
  * bdsim_run    - does BDSIM accept and run the generated GMAD
  * g4bl_run     - does G4beamline run the original
  * trackable    - does the original have a beam + a detector
  * class        - closure classification (see CLASS_* below)

Quantitative plane-by-plane closure is produced separately by compare.py; this
sweep is the breadth pass over many examples.

Usage:
    python3 validation/sweep.py <dir-or-files...> [--g4bl] [--md out.md]
"""

from __future__ import annotations

import argparse
import glob
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from typing import List, Optional

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

from g4beamline2bdsim.converter import Converter            # noqa: E402
from g4beamline2bdsim.gmad_writer import GmadWriter         # noqa: E402
from g4beamline2bdsim.parser import parse_g4bl              # noqa: E402

G4BL_IMAGE = os.environ.get("G4BL_IMAGE", "g4beamline:3.06")
BDSIM_IMAGE = os.environ.get("BDSIM_IMAGE", "bdsim/ubuntu24-g4.11.3-bdsim:develop")
BDSIM_ENV = ("source /usr/local/bin/geant4.sh 2>/dev/null; "
             "source /bdsim/root/bin/thisroot.sh 2>/dev/null")

# Features G4beamline supports that the converter intentionally does not map to
# the 1-D BDSIM beamline (so closure is not expected if present).
# Features still not converted to the 1-D BDSIM beamline (fieldmap/fieldexpr are
# now converted to BDSIM field maps, so they are no longer listed here).
UNSUPPORTED = ["spacecharge", "helicaldipole", "fieldlines"]


@dataclass
class Result:
    name: str
    n_elements: int = 0
    n_warnings: int = 0
    convert_ok: bool = False
    convert_err: str = ""
    bdsim_ok: Optional[bool] = None
    bdsim_err: str = ""
    g4bl_ok: Optional[bool] = None
    has_beam: bool = False
    has_detector: bool = False
    unsupported: List[str] = field(default_factory=list)
    klass: str = ""


def classify(r: Result) -> str:
    if not r.convert_ok:
        return "CONVERT-ERROR"
    if not r.has_beam and not r.has_detector:
        return "NON-TRACKING"          # field-line / visualisation only
    if r.unsupported:
        return "PARTIAL (unsupported features)"
    if r.bdsim_ok is False:
        return "BDSIM-ERROR"
    if r.has_beam and r.has_detector and r.bdsim_ok:
        return "TRACKABLE"
    return "CONVERTS"


def analyse(path: str, run_bdsim: bool, run_g4bl: bool, workroot: str) -> Result:
    name = os.path.basename(path)
    r = Result(name=name)
    text = open(path, "r", encoding="utf-8", errors="replace").read()

    low = text.lower()
    r.has_beam = any(l.strip().startswith("beam ")
                     for l in text.splitlines())
    r.has_detector = any(l.strip().split()[0] in ("virtualdetector", "detector")
                         for l in text.splitlines() if l.strip())
    r.unsupported = [u for u in UNSUPPORTED if u in low]

    # Convert.
    try:
        commands, _ = parse_g4bl(text, base_dir=os.path.dirname(os.path.abspath(path)))
        model = Converter(commands, source_name=name,
                          base_dir=os.path.dirname(os.path.abspath(path))).convert()
        r.convert_ok = True
        r.n_elements = len(model.elements)
        r.n_warnings = len(model.warnings)
    except Exception as exc:  # noqa: BLE001
        r.convert_err = str(exc)
        r.klass = classify(r)
        return r

    work = os.path.join(workroot, os.path.splitext(name)[0])
    os.makedirs(work, exist_ok=True)
    gmad = os.path.splitext(name)[0] + ".gmad"
    GmadWriter(model).write(os.path.join(work, gmad))
    # Write auxiliary files (GDML geometry, field maps) next to the GMAD.
    model.write_aux_files(work)

    if run_bdsim and r.n_elements > 0:
        script = (f"{BDSIM_ENV}; bdsim --file={gmad} --outfile=o --batch "
                  f"--ngenerate=10 > b.log 2>&1; echo RC=$?")
        try:
            res = subprocess.run(
                ["docker", "run", "--rm", "--entrypoint", "bash",
                 "-v", f"{work}:/work", "-w", "/work", BDSIM_IMAGE, "-lc", script],
                capture_output=True, text=True, timeout=600)
            log = ""
            logp = os.path.join(work, "b.log")
            if os.path.exists(logp):
                log = open(logp).read()
            r.bdsim_ok = ("RC=0" in res.stdout) and ("End of Run" in log)
            if not r.bdsim_ok:
                for line in log.splitlines():
                    low = line.lower()
                    # Skip the harmless ROOT autoload/cling warnings.
                    if any(s in line for s in
                           ("cling::", "FileEntry", "autoload")):
                        continue
                    if any(k in low for k in
                           ("error", "malformed", "unexpected", "fatal",
                            "required to be set", "must be", "unknown")):
                        r.bdsim_err = line.strip()[:160]
                        break
        except subprocess.TimeoutExpired:
            r.bdsim_ok = False
            r.bdsim_err = "timeout"

    if run_g4bl and r.has_beam:
        # Cap events for speed by overriding via a wrapper file is complex;
        # just run as-is with a timeout.
        g4 = name
        shutil.copy(path, os.path.join(work, g4))
        script = (f"source $G4BL_DIR/bin/g4bl-setup.sh 2>/dev/null; "
                  f"g4bl {g4} first=1 last=5 > g.log 2>&1; echo RC=$?")
        try:
            res = subprocess.run(
                ["docker", "run", "--rm", "-v", f"{work}:/work", "-w", "/work",
                 G4BL_IMAGE, "bash", "-lc", script],
                capture_output=True, text=True, timeout=600)
            r.g4bl_ok = "RC=0" in res.stdout
        except subprocess.TimeoutExpired:
            r.g4bl_ok = False

    r.klass = classify(r)
    return r


def collect(paths: List[str]) -> List[str]:
    files: List[str] = []
    for p in paths:
        if os.path.isdir(p):
            files += sorted(glob.glob(os.path.join(p, "*.g4bl")))
        else:
            files.append(p)
    return files


def to_markdown(results: List[Result]) -> str:
    out = ["| example | elems | warns | beam | det | bdsim | g4bl | unsupported | class |",
           "|---|---|---|---|---|---|---|---|---|"]
    def b(x):
        return "-" if x is None else ("ok" if x else "FAIL")
    for r in results:
        out.append(
            f"| {r.name} | {r.n_elements} | {r.n_warnings} | "
            f"{'Y' if r.has_beam else '.'} | {'Y' if r.has_detector else '.'} | "
            f"{b(r.bdsim_ok)} | {b(r.g4bl_ok)} | "
            f"{','.join(r.unsupported) or '-'} | {r.klass} |")
    return "\n".join(out)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("paths", nargs="+")
    ap.add_argument("--bdsim", action="store_true", help="run BDSIM on each")
    ap.add_argument("--g4bl", action="store_true", help="run G4beamline on each")
    ap.add_argument("--workroot", default="/tmp/g4bl2bdsim_sweep")
    ap.add_argument("--md", help="write a markdown table to this path")
    args = ap.parse_args(argv)

    os.makedirs(args.workroot, exist_ok=True)
    results = []
    for path in collect(args.paths):
        r = analyse(path, args.bdsim, args.g4bl, args.workroot)
        print(f"{r.klass:28s} {r.name:28s} elems={r.n_elements} "
              f"warns={r.n_warnings} bdsim={r.bdsim_ok} g4bl={r.g4bl_ok}"
              + (f"  [{r.bdsim_err}]" if r.bdsim_err else ""))
        results.append(r)

    table = to_markdown(results)
    print("\n" + table)
    if args.md:
        open(args.md, "w").write(table + "\n")

    # summary counts
    from collections import Counter
    counts = Counter(r.klass for r in results)
    print("\n=== class counts ===")
    for k, v in counts.most_common():
        print(f"  {v:3d}  {k}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
