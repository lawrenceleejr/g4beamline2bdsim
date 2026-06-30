#!/usr/bin/env python3
"""Validate dipole conversion by comparing the G4beamline lab-frame deflection
to the analytic bend angle (which is what BDSIM realises).

BDSIM samplers live in the bent local frame, so an on-design particle has
x=xp=0 after a bend -- a direct sampler-to-sampler comparison is degenerate.
Instead we run only G4beamline (lab frame), measure the outgoing angle of the
reference particle, and compare it with theta = B * L / Brho, the angle BDSIM's
rbend/sbend produces for the same field, length and rigidity.

Usage:
    python3 validation/check_bend.py [--B 0.1 --L 0.5 --p 1000]
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import g4bl_ascii  # noqa: E402

G4BL_IMAGE = os.environ.get("G4BL_IMAGE", "g4beamline:3.06")

TEMPLATE = """\
physics QGSP_BERT
param Pz={p}
reference referenceMomentum=$Pz particle=proton beamZ=0
beam gaussian particle=proton nEvents={n} meanMomentum=$Pz beamZ=0 \\
     sigmaX=0.05 sigmaY=0.05 sigmaXp=0.00003 sigmaYp=0.00003 sigmaP=0
genericbend D fieldWidth=800 fieldHeight=600 fieldLength={Lmm} ironLength=0 fringe=0 kill=0
virtualdetector Det radius=3000 length=1 format=ascii
place D rename=B1 By={B} z=1000
place Det rename=DetEnd z={detz}
"""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=float, default=0.1, help="field [T]")
    ap.add_argument("--L", type=float, default=0.5, help="length [m]")
    ap.add_argument("--p", type=float, default=1000.0, help="momentum [MeV/c]")
    ap.add_argument("--n", type=int, default=400)
    ap.add_argument("--rtol", type=float, default=0.03)
    args = ap.parse_args(argv)

    brho = (args.p * 1e-3) / 0.299792458
    theta = args.B * args.L / brho

    work = tempfile.mkdtemp(prefix="bendcheck_")
    g4bl = os.path.join(work, "bend.g4bl")
    with open(g4bl, "w") as fh:
        fh.write(TEMPLATE.format(
            p=args.p, n=args.n, Lmm=args.L * 1000.0, B=args.B,
            detz=1000 + args.L * 1000.0 / 2 + 700))

    script = "source $G4BL_DIR/bin/g4bl-setup.sh 2>/dev/null; g4bl bend.g4bl >/dev/null 2>&1; echo done"
    subprocess.run(
        ["docker", "run", "--rm", "-v", f"{work}:/work", "-w", "/work",
         G4BL_IMAGE, "bash", "-lc", script],
        check=True, capture_output=True, text=True, timeout=600)

    tracks = g4bl_ascii.primaries(
        g4bl_ascii.parse_bltrackfile(os.path.join(work, "DetEnd.txt")), 2212)
    s = g4bl_ascii.stats(tracks)
    g4bl_theta = abs(s["mean_xp"])
    rel = abs(g4bl_theta - theta) / theta if theta else 0.0

    print(f"B={args.B} T  L={args.L} m  p={args.p} MeV/c  Brho={brho:.5f} T*m")
    print(f"  analytic theta = B*L/Brho = {theta:.6f} rad  (BDSIM rbend angle)")
    print(f"  G4beamline lab deflection = {g4bl_theta:.6f} rad  (n={s['n']})")
    print(f"  relative difference = {rel*100:.2f}%")
    ok = rel <= args.rtol
    print("  ->", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
