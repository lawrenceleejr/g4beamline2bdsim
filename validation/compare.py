#!/usr/bin/env python3
"""End-to-end validation: run G4beamline and BDSIM on the same lattice and
compare the beam at equivalent detector / sampler planes.

Steps
-----
1. Run the original ``.g4bl`` file in the G4beamline container -> ASCII
   ``<detector>.txt`` NTuples.
2. Convert the ``.g4bl`` to ``.gmad`` with the converter under test.
3. Run the ``.gmad`` in the BDSIM container -> ROOT, then dump each sampler to
   ``<sampler>.dat`` (via ``bdsim_dump.py`` inside the container).
4. For every detector present in both, compare the primary-particle statistics
   (mean/sigma of x, y, xp, yp) and report pass/fail against tolerances.

This module is meant to be run on the host (it shells out to ``docker``).
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from typing import Dict, List, Optional

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

from g4beamline2bdsim.converter import Converter            # noqa: E402
from g4beamline2bdsim.gmad_writer import GmadWriter         # noqa: E402
from g4beamline2bdsim.parser import parse_g4bl              # noqa: E402
import g4bl_ascii                                           # noqa: E402

G4BL_IMAGE = os.environ.get("G4BL_IMAGE", "g4beamline:3.06")
BDSIM_IMAGE = os.environ.get("BDSIM_IMAGE", "bdsim/ubuntu24-g4.11.3-bdsim:develop")

# Environment setup commands inside each container.
G4BL_ENV = "source $G4BL_DIR/bin/g4bl-setup.sh 2>/dev/null"
BDSIM_ENV = (
    "source /usr/local/bin/geant4.sh 2>/dev/null; "
    "source /bdsim/root/bin/thisroot.sh 2>/dev/null"
)


def run(cmd: List[str], **kw) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def _force_ascii_detectors(text: str) -> str:
    """Ensure every virtualdetector/detector definition writes ASCII output."""
    out = []
    for line in text.splitlines():
        stripped = line.lstrip()
        first = stripped.split()[0] if stripped.split() else ""
        if first in ("virtualdetector", "detector") and "format=" not in line:
            line = line.rstrip() + " format=ascii"
        out.append(line)
    return "\n".join(out) + "\n"


def docker_run(image: str, workdir_mount: str, script: str,
               entrypoint: Optional[str] = "bash", extra: Optional[List[str]] = None,
               timeout: int = 1200) -> subprocess.CompletedProcess:
    cmd = ["docker", "run", "--rm"]
    if entrypoint:
        cmd += ["--entrypoint", entrypoint]
    cmd += ["-v", f"{workdir_mount}:/work", "-w", "/work"]
    if extra:
        cmd += extra
    cmd += [image, "-lc", script]
    return run(cmd, timeout=timeout)


def run_g4beamline(work: str, g4bl_name: str, n_events: int) -> None:
    script = f"{G4BL_ENV}; g4bl {g4bl_name} 2>&1 | tail -3"
    res = docker_run(G4BL_IMAGE, work, script)
    if res.returncode != 0:
        print(res.stdout)
        print(res.stderr, file=sys.stderr)
        raise RuntimeError("G4beamline run failed")


def run_bdsim(work: str, gmad_name: str, n_events: int, dump_script: str) -> None:
    script = (
        f"{BDSIM_ENV}; "
        f"bdsim --file={gmad_name} --outfile=bdsim_out --batch "
        f"--ngenerate={n_events} > bdsim.log 2>&1; "
        f"python3 {dump_script} bdsim_out.root . > dump.log 2>&1; "
        f"echo BDSIM_DONE"
    )
    res = docker_run(BDSIM_IMAGE, work, script)
    if "BDSIM_DONE" not in res.stdout:
        print(res.stdout)
        print(res.stderr, file=sys.stderr)
        # surface the in-container logs
        for log in ("bdsim.log", "dump.log"):
            p = os.path.join(work, log)
            if os.path.exists(p):
                print(f"--- {log} ---")
                print(open(p).read()[-2000:])
        raise RuntimeError("BDSIM run failed")


def bdsim_stats(path: str, pdgid: int) -> Dict[str, float]:
    """Stats from a BDSIM sampler .dat (x,y in m -> mm; xp,yp rad)."""
    xs, ys, xps, yps, ps = [], [], [], [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            x, y, z, xp, yp, energy = (float(v) for v in p[:6])
            partid, parentid = int(p[6]), int(p[7])
            if parentid != 0:
                continue
            if pdgid is not None and partid != pdgid:
                continue
            xs.append(x * 1e3)   # m -> mm
            ys.append(y * 1e3)
            xps.append(xp)
            yps.append(yp)
            ps.append(energy)
    n = len(xs)
    if n == 0:
        return {"n": 0}

    def mean(v):
        return sum(v) / len(v)

    def sig(v, m):
        return (sum((a - m) ** 2 for a in v) / (len(v) - 1)) ** 0.5 if len(v) > 1 else 0.0

    mx, my, mxp, myp = mean(xs), mean(ys), mean(xps), mean(yps)
    return {
        "n": n,
        "mean_x_mm": mx, "mean_y_mm": my, "mean_xp": mxp, "mean_yp": myp,
        "sigma_x_mm": sig(xs, mx), "sigma_y_mm": sig(ys, my),
        "sigma_xp": sig(xps, mxp), "sigma_yp": sig(yps, myp),
    }


def compare_case(g4bl_path: str, n_events: int, workroot: str,
                 atol_mm: float, atol_ang: float, rtol: float,
                 smoke: bool = False) -> bool:
    name = os.path.splitext(os.path.basename(g4bl_path))[0]
    work = os.path.join(workroot, name)
    if os.path.exists(work):
        shutil.rmtree(work)
    os.makedirs(work)

    # Stage inputs.  Force ASCII detector output so we can parse it on the host
    # (g4beamline defaults to ROOT otherwise).
    g4bl_name = os.path.basename(g4bl_path)
    raw = open(g4bl_path).read()
    with open(os.path.join(work, g4bl_name), "w") as fh:
        fh.write(_force_ascii_detectors(raw))
    shutil.copy(os.path.join(HERE, "bdsim_dump.py"),
                os.path.join(work, "bdsim_dump.py"))

    # Convert.
    text = open(g4bl_path).read()
    commands, _ = parse_g4bl(text)
    converter = Converter(commands, source_name=g4bl_name,
                          base_dir=os.path.dirname(os.path.abspath(g4bl_path)))
    model = converter.convert()
    gmad_name = name + ".gmad"
    GmadWriter(model).write(os.path.join(work, gmad_name))
    # Write auxiliary files (GDML geometry, field maps) next to the GMAD.
    for fname, content in model.aux_files.items():
        with open(os.path.join(work, fname), "w") as fh:
            fh.write(content)

    beam_particle = model.beam.get("particle", "proton")
    pdgid = g4bl_ascii.PDGID.get(beam_particle, 2212)
    samplers = [s for s in model.samplers if s != "all"]

    print(f"\n=== case: {name} (particle={beam_particle}, n={n_events}) ===")
    print(f"    detectors/samplers: {samplers}")

    # Run both codes.
    run_g4beamline(work, g4bl_name, n_events)
    run_bdsim(work, gmad_name, n_events, "bdsim_dump.py")

    ok = True
    for sampler in samplers:
        g4_txt = os.path.join(work, sampler + ".txt")
        bd_dat = os.path.join(work, sampler + ".dat")
        if not os.path.exists(g4_txt):
            print(f"  [{sampler}] MISSING g4beamline output ({sampler}.txt)")
            ok = False
            continue
        if not os.path.exists(bd_dat):
            print(f"  [{sampler}] MISSING bdsim output ({sampler}.dat)")
            ok = False
            continue
        g_tracks = g4bl_ascii.primaries(
            g4bl_ascii.parse_bltrackfile(g4_txt), pdgid)
        g_stats = g4bl_ascii.stats(g_tracks)
        b_stats = bdsim_stats(bd_dat, pdgid)
        ok &= _report(sampler, g_stats, b_stats, atol_mm, atol_ang, rtol, smoke)
    return ok


def _report(sampler, g, b, atol_mm, atol_ang, rtol, smoke=False) -> bool:
    print(f"  [{sampler}] g4bl n={g.get('n')}  bdsim n={b.get('n')}")
    if g.get("n", 0) == 0 or b.get("n", 0) == 0:
        # In smoke mode, both codes ran; a plane with no primaries (e.g. fully
        # stopped in iron/material) is reported but not a hard failure.
        print("    -> no primaries at this plane"
              + (" (smoke: informational)" if smoke else ""))
        return smoke
    if smoke:
        for key in ("mean_x_mm", "mean_y_mm", "mean_xp", "mean_yp",
                    "sigma_x_mm", "sigma_y_mm"):
            print(f"    -- {key:12s} g4bl={g.get(key,0):+.5g} "
                  f"bdsim={b.get(key,0):+.5g}")
        return True

    import math as _m
    ng, nb = g.get("n", 1), b.get("n", 1)
    # (field, abs-tol, unit, sigma-key for statistical error, stat factor)
    fields = [
        ("mean_x_mm", atol_mm, "mm", "sigma_x_mm", 1.0),
        ("mean_y_mm", atol_mm, "mm", "sigma_y_mm", 1.0),
        ("mean_xp", atol_ang, "rad", "sigma_xp", 1.0),
        ("mean_yp", atol_ang, "rad", "sigma_yp", 1.0),
        # sigma estimator error ~ sigma/sqrt(2N).
        ("sigma_x_mm", atol_mm, "mm", "sigma_x_mm", 0.7071),
        ("sigma_y_mm", atol_mm, "mm", "sigma_y_mm", 0.7071),
    ]
    ok = True
    for key, atol, unit, sigkey, fac in fields:
        gv = g.get(key, 0.0)
        bv = b.get(key, 0.0)
        diff = abs(gv - bv)
        # 3-sigma combined statistical uncertainty from both samples.
        sg = g.get(sigkey, 0.0) * fac
        sb = b.get(sigkey, 0.0) * fac
        stat = 3.0 * _m.sqrt(sg * sg / max(ng, 1) + sb * sb / max(nb, 1))
        tol = atol + rtol * max(abs(gv), abs(bv)) + stat
        status = "ok " if diff <= tol else "BAD"
        if diff > tol:
            ok = False
        print(f"    {status} {key:12s} g4bl={gv:+.5g} bdsim={bv:+.5g} "
              f"|d|={diff:.3g} {unit} (tol {tol:.3g}, stat {stat:.3g})")
    return ok


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("cases", nargs="+", help="G4beamline case files (.g4bl)")
    ap.add_argument("-n", "--n-events", type=int, default=2000)
    ap.add_argument("--workroot", default="/tmp/g4bl2bdsim_validation")
    ap.add_argument("--atol-mm", type=float, default=0.5,
                    help="absolute position tolerance [mm]")
    ap.add_argument("--atol-ang", type=float, default=5e-4,
                    help="absolute angle tolerance [rad]")
    ap.add_argument("--rtol", type=float, default=0.05,
                    help="relative tolerance")
    ap.add_argument("--smoke", action="store_true",
                    help="smoke mode: just confirm both codes run and produce "
                         "output; report stats without enforcing tolerances")
    args = ap.parse_args(argv)

    os.makedirs(args.workroot, exist_ok=True)
    all_ok = True
    results = []
    for case in args.cases:
        try:
            ok = compare_case(case, args.n_events, args.workroot,
                              args.atol_mm, args.atol_ang, args.rtol,
                              smoke=args.smoke)
        except Exception as exc:  # noqa: BLE001
            print(f"  ERROR: {exc}")
            ok = False
        results.append((case, ok))
        all_ok &= ok

    print("\n=== summary ===")
    for case, ok in results:
        print(f"  {'PASS' if ok else 'FAIL'}  {os.path.basename(case)}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
