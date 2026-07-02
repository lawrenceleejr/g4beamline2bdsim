"""Write BDSIM-format field-map files and sample fields onto grids.

BDSIM field-map ASCII format — verified against the BDSIM loader source
(``src/BDSFieldLoaderBDSIM.cc``):

* header lines ``xmin> <v>`` / ``xmax> <v>`` / ``nx> <n>`` (and ``y``/``z``/``t``
  for higher dimensions) — **positions in centimetres** (the loader multiplies
  header values and coordinate columns by ``CLHEP::cm``);
* an optional ``loopOrder> xyzt`` header (default ``xyzt``);
* a column header comment ``! X Y Z Fx Fy Fz``;
* one whitespace/tab-separated row per grid point — field values are raw and
  interpreted as **Tesla** for magnetic maps;
* **loop order ``xyzt``: x varies fastest (innermost), then y, then z, then t**.
  The coordinate columns are read but the grid indexing comes purely from the
  header, so the row ORDER is what matters.

A field object is declared in GMAD as::

    fm: field, type="bmap3d", magneticFile="bdsim3d:fm.dat",
        magneticInterpolator="linear";

and attached to a beamline element with ``fieldAll="fm"`` (name quoted).

This module's public API takes positions in **millimetres** (the G4beamline
convention used throughout the converter) and converts to cm on write.
"""

from __future__ import annotations

from typing import Callable, List, Sequence, Tuple

# A field function maps (x, y, z) in mm -> (Bx, By, Bz) in Tesla.
FieldFn = Callable[[float, float, float], Tuple[float, float, float]]

_MM_TO_CM = 0.1


def linspace(a: float, b: float, n: int) -> List[float]:
    if n <= 1:
        return [0.5 * (a + b)]
    step = (b - a) / (n - 1)
    return [a + step * i for i in range(n)]


def _fmt(v: float) -> str:
    return f"{v:.8E}"


def build_3d(xs: Sequence[float], ys: Sequence[float], zs: Sequence[float],
             field_fn: FieldFn) -> str:
    """Return 3D BDSIM magnetic field-map text.

    ``xs``/``ys``/``zs`` are the grid node positions in **mm**; ``field_fn``
    is called with (x, y, z) in mm and must return (Bx, By, Bz) in Tesla.
    Positions are written in cm and rows in ``xyzt`` order (x fastest).
    """
    c = _MM_TO_CM
    lines = [
        f"xmin> {xs[0]*c:.6g}", f"xmax> {xs[-1]*c:.6g}", f"nx> {len(xs)}",
        f"ymin> {ys[0]*c:.6g}", f"ymax> {ys[-1]*c:.6g}", f"ny> {len(ys)}",
        f"zmin> {zs[0]*c:.6g}", f"zmax> {zs[-1]*c:.6g}", f"nz> {len(zs)}",
        "loopOrder> xyzt",
        "! X\tY\tZ\tFx\tFy\tFz",
    ]
    out = "\n".join(lines) + "\n"
    rows = []
    for z in zs:                     # z outermost
        for y in ys:                 # y middle
            for x in xs:             # x innermost (fastest)
                bx, by, bz = field_fn(x, y, z)
                rows.append(f"{_fmt(x*c)}\t{_fmt(y*c)}\t{_fmt(z*c)}\t"
                            f"{_fmt(bx)}\t{_fmt(by)}\t{_fmt(bz)}")
    return out + "\n".join(rows) + "\n"


def build_2d(xs: Sequence[float], ys: Sequence[float], field_fn: FieldFn,
             z: float = 0.0) -> str:
    """Return 2D BDSIM magnetic field-map text over (x, y) at fixed z (mm)."""
    c = _MM_TO_CM
    lines = [
        f"xmin> {xs[0]*c:.6g}", f"xmax> {xs[-1]*c:.6g}", f"nx> {len(xs)}",
        f"ymin> {ys[0]*c:.6g}", f"ymax> {ys[-1]*c:.6g}", f"ny> {len(ys)}",
        "loopOrder> xyzt",
        "! X\tY\tFx\tFy\tFz",
    ]
    out = "\n".join(lines) + "\n"
    rows = []
    for y in ys:                     # y outer
        for x in xs:                 # x fastest
            bx, by, bz = field_fn(x, y, z)
            rows.append(f"{_fmt(x*c)}\t{_fmt(y*c)}\t"
                        f"{_fmt(bx)}\t{_fmt(by)}\t{_fmt(bz)}")
    return out + "\n".join(rows) + "\n"


def gmad_field_object(name: str, filename: str, dims: int,
                      interpolator: str = "linear") -> str:
    """The GMAD ``field`` object declaration for a written map file."""
    return (f'{name}: field, type="bmap{dims}d", '
            f'magneticFile="bdsim{dims}d:{filename}", '
            f'magneticInterpolator="{interpolator}";')
