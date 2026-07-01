"""Dump BDSIM sampler data to a simple ASCII table.

Runs *inside* the BDSIM container (needs ROOT + libbdsimRootEvent).  For each
requested sampler it writes ``<outdir>/<sampler>.dat`` with one row per particle:

    x y z xp yp energy partID parentID eventID

Units follow BDSIM: x,y,z in metres; xp,yp in radians; energy in GeV.

Usage (inside container):
    python3 bdsim_dump.py <rootfile> <outdir> <sampler1> [sampler2 ...]
If no samplers are given, every sampler branch in the Event tree is dumped.
"""

import os
import sys

import ROOT

# The shared library carrying the BDSIM ROOT-event class dictionaries.
for cand in (
    "/usr/local/lib/libbdsimRootEvent.so",
    "libbdsimRootEvent.so",
):
    if ROOT.gSystem.Load(cand) >= 0:
        break

# Sampler branches are these BDSIM classes; everything else in Event is skipped.
_NON_SAMPLER = {
    "Summary.", "Primary.", "PrimaryGlobal.", "Eloss.", "PrimaryFirstHit.",
    "PrimaryLastHit.", "ApertureImpacts.", "Histos.", "Info.", "Trajectory.",
    "ElossVacuum.", "ElossTunnel.", "ElossWorld.", "ElossWorldContents.",
    "ElossWorldExit.",
}


def sampler_branches(tree):
    names = []
    for b in tree.GetListOfBranches():
        nm = b.GetName()
        if nm in _NON_SAMPLER:
            continue
        names.append(nm)
    return names


def dump(rootfile, outdir, samplers):
    f = ROOT.TFile(rootfile)
    tree = f.Get("Event")
    if not tree:
        raise SystemExit(f"no Event tree in {rootfile}")

    if not samplers:
        samplers = sampler_branches(tree)

    os.makedirs(outdir, exist_ok=True)
    nentries = tree.GetEntries()

    for sampler in samplers:
        branchname = sampler if sampler.endswith(".") else sampler + "."
        attr = sampler.rstrip(".")
        outpath = os.path.join(outdir, attr + ".dat")
        rows = 0
        with open(outpath, "w") as out:
            out.write("# x y z xp yp energy partID parentID eventID\n")
            out.write("# m m m rad rad GeV - - -\n")
            for i in range(nentries):
                tree.GetEntry(i)
                s = getattr(tree, attr, None)
                if s is None:
                    continue
                n = s.n
                # In a BDSIM sampler x,y,xp,yp,energy,partID,parentID are
                # per-particle vectors, but z (the sampler plane position) is a
                # scalar.  The event index is the tree entry number i.
                z = float(s.z)
                for j in range(n):
                    out.write(
                        "%.9g %.9g %.9g %.9g %.9g %.9g %d %d %d\n"
                        % (
                            s.x[j], s.y[j], z,
                            s.xp[j], s.yp[j], s.energy[j],
                            s.partID[j], s.parentID[j], i,
                        )
                    )
                    rows += 1
        print(f"wrote {outpath} ({rows} rows)")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise SystemExit("usage: bdsim_dump.py <rootfile> <outdir> [samplers...]")
    dump(sys.argv[1], sys.argv[2], sys.argv[3:])
