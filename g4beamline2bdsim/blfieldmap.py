"""Parse G4beamline BLFieldMap files and convert them to BDSIM field maps.

BLFieldMap ``grid`` format (from G4beamline ``BLFieldMap.cc``)::

    param normB=1.0 ...            # optional global scaling
    grid nX=.. nY=.. nZ=.. dX=.. dY=.. dZ=.. X0=.. Y0=.. Z0=.. \
         extendX=0 extendY=0 extendZ=0
    data
    X Y Z Bx By Bz [Ex Ey Ez]      # positions mm, B Tesla, E MV/m
    ...

Both codes store positions in mm and B in Tesla, so the conversion is mostly a
re-grid + re-header into the BDSIM 3D format.  Only the magnetic field is
carried across (BDSIM would need an ebmap for E).  ``extend*`` mirror-symmetry
flags are honoured by reflecting the stored data onto the full grid.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


@dataclass
class BLFieldGrid:
    xs: List[float] = field(default_factory=list)      # sorted unique X (mm)
    ys: List[float] = field(default_factory=list)
    zs: List[float] = field(default_factory=list)
    # (ix,iy,iz) -> (Bx,By,Bz) in Tesla
    b: Dict[Tuple[int, int, int], Tuple[float, float, float]] = field(default_factory=dict)
    has_efield: bool = False
    norm_b: float = 1.0

    def field_fn(self):
        xi = {round(v, 6): i for i, v in enumerate(self.xs)}
        yi = {round(v, 6): i for i, v in enumerate(self.ys)}
        zi = {round(v, 6): i for i, v in enumerate(self.zs)}

        def fn(x, y, z):
            key = (xi.get(round(x, 6)), yi.get(round(y, 6)), zi.get(round(z, 6)))
            bx, by, bz = self.b.get(key, (0.0, 0.0, 0.0))
            return bx * self.norm_b, by * self.norm_b, bz * self.norm_b

        return fn


def _tokenize(line: str) -> List[str]:
    return line.replace(",", " ").split()


def _named(tokens: List[str]) -> Dict[str, str]:
    out = {}
    for t in tokens:
        if "=" in t:
            k, _, v = t.partition("=")
            out[k] = v
    return out


def parse(path: str) -> Optional[BLFieldGrid]:
    """Parse a BLFieldMap ``grid`` file.  Returns None if not a grid map."""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()

    grid = BLFieldGrid()
    extend = {"extendX": 0, "extendY": 0, "extendZ": 0}
    section = None
    is_grid = False
    raw_points: List[Tuple[float, float, float, float, float, float]] = []

    for line in lines:
        s = line.strip()
        if not s or s[0] in "#*":
            continue
        tokens = _tokenize(s)
        head = tokens[0]
        if head == "param":
            nm = _named(tokens[1:])
            if "normB" in nm:
                try:
                    grid.norm_b = float(nm["normB"])
                except ValueError:
                    pass
            continue
        if head == "grid":
            is_grid = True
            nm = _named(tokens[1:])
            for k in extend:
                if k in nm:
                    try:
                        extend[k] = int(float(nm[k]))
                    except ValueError:
                        pass
            section = None
            continue
        if head == "cylinder":
            return None            # cylindrical maps not supported here
        if head == "data":
            section = "data"
            continue
        if head and head[0].isalpha():
            # Any other keyword line ends the data section.
            section = None
            continue
        if section == "data":
            vals = [float(v) for v in tokens]
            if len(vals) < 6:
                continue
            x, y, z, bx, by, bz = vals[:6]
            if len(vals) > 6:
                grid.has_efield = True
            raw_points.append((x, y, z, bx, by, bz))

    if not is_grid or not raw_points:
        return None

    # Apply mirror symmetry (extend*) to build the full set of points.
    points = list(raw_points)
    for axis, flagkey in ((0, "extendX"), (1, "extendY"), (2, "extendZ")):
        if extend[flagkey]:
            mirrored = []
            for p in points:
                if p[axis] != 0.0:
                    q = list(p)
                    q[axis] = -q[axis]
                    # B parity under reflection: the component along the mirror
                    # axis is even, the transverse components flip sign.
                    for bcomp in (3, 4, 5):
                        if (bcomp - 3) != axis:
                            q[bcomp] = -q[bcomp]
                    mirrored.append(tuple(q))
            points.extend(mirrored)

    xs = sorted({round(p[0], 6) for p in points})
    ys = sorted({round(p[1], 6) for p in points})
    zs = sorted({round(p[2], 6) for p in points})
    grid.xs, grid.ys, grid.zs = xs, ys, zs
    xi = {v: i for i, v in enumerate(xs)}
    yi = {v: i for i, v in enumerate(ys)}
    zi = {v: i for i, v in enumerate(zs)}
    for x, y, z, bx, by, bz in points:
        grid.b[(xi[round(x, 6)], yi[round(y, 6)], zi[round(z, 6)])] = (bx, by, bz)
    return grid
