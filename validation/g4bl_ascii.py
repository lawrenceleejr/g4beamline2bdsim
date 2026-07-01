"""Parse G4beamline ASCII (BLTrackFile) NTuple output.

Format (written by ``virtualdetector ... format=ascii``)::

    #BLTrackFile VirtualDetector/DetEnd
    #x y z Px Py Pz t PDGid EventID TrackID ParentID Weight
    #mm mm mm MeV/c MeV/c MeV/c ns - - - - -
    -2.849 -0.44 999.5 -1.838 -0.43 999.998 4.57 2212 1 1 0 1
    ...

Coordinates are centreline x,y,z [mm] and momentum components [MeV/c].
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional


@dataclass
class Track:
    x: float          # mm
    y: float          # mm
    z: float          # mm
    px: float         # MeV/c
    py: float         # MeV/c
    pz: float         # MeV/c
    t: float          # ns
    pdgid: int
    event_id: int
    track_id: int
    parent_id: int
    weight: float

    @property
    def xp(self) -> float:
        return self.px / self.pz if self.pz else 0.0

    @property
    def yp(self) -> float:
        return self.py / self.pz if self.pz else 0.0

    @property
    def p(self) -> float:
        return (self.px ** 2 + self.py ** 2 + self.pz ** 2) ** 0.5


def parse_bltrackfile(path: str) -> List[Track]:
    tracks: List[Track] = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 12:
                continue
            try:
                tracks.append(
                    Track(
                        x=float(parts[0]),
                        y=float(parts[1]),
                        z=float(parts[2]),
                        px=float(parts[3]),
                        py=float(parts[4]),
                        pz=float(parts[5]),
                        t=float(parts[6]),
                        pdgid=int(float(parts[7])),
                        event_id=int(float(parts[8])),
                        track_id=int(float(parts[9])),
                        parent_id=int(float(parts[10])),
                        weight=float(parts[11]),
                    )
                )
            except ValueError:
                continue
    return tracks


def primaries(tracks: List[Track], pdgid: Optional[int] = None) -> List[Track]:
    out = [t for t in tracks if t.parent_id == 0]
    if pdgid is not None:
        out = [t for t in out if t.pdgid == pdgid]
    return out


def stats(tracks: List[Track]) -> Dict[str, float]:
    """Mean and sigma of x,y [mm] and xp,yp [rad] over a track list."""
    n = len(tracks)
    if n == 0:
        return {"n": 0}

    def mean(vals):
        return sum(vals) / len(vals)

    def sigma(vals, m):
        if len(vals) < 2:
            return 0.0
        return (sum((v - m) ** 2 for v in vals) / (len(vals) - 1)) ** 0.5

    xs = [t.x for t in tracks]
    ys = [t.y for t in tracks]
    xps = [t.xp for t in tracks]
    yps = [t.yp for t in tracks]
    ps = [t.p for t in tracks]
    mx, my, mxp, myp, mp = mean(xs), mean(ys), mean(xps), mean(yps), mean(ps)
    return {
        "n": n,
        "mean_x_mm": mx,
        "mean_y_mm": my,
        "mean_xp": mxp,
        "mean_yp": myp,
        "mean_p_MeV": mp,
        "sigma_x_mm": sigma(xs, mx),
        "sigma_y_mm": sigma(ys, my),
        "sigma_xp": sigma(xps, mxp),
        "sigma_yp": sigma(yps, myp),
    }


# PDG ids for the common beam particles.
PDGID = {
    "proton": 2212,
    "e-": 11,
    "e+": -11,
    "mu-": 13,
    "mu+": -13,
    "pi+": 211,
    "pi-": -211,
    "gamma": 22,
    "neutron": 2112,
}
