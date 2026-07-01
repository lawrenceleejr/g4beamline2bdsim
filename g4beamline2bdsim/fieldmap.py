"""Write BDSIM-format field-map files and sample fields onto grids.

BDSIM field-map ASCII format (verified empirically against BDSIM 1.7.7):

* header lines ``xmin> <v>`` / ``xmax> <v>`` / ``nx> <n>`` (and ``y``/``z``/``t``
  for higher dimensions), **positions in millimetres**;
* a column header comment ``! X Y Z Fx Fy Fz``;
* one whitespace/tab-separated row per grid point, **field in Tesla**;
* **loop order: first coordinate outermost, last coordinate innermost**
  (3D: X outer, Y middle, Z fastest; 2D: X outer, Y fastest).

A field object is declared in GMAD as::

    fm: field, type="bmap3d", magneticFile="bdsim3d:fm.dat",
        magneticInterpolator="linear";

and attached to a beamline element with ``fieldAll="fm"`` (name quoted).
"""

from __future__ import annotations

from typing import Callable, List, Sequence, Tuple

# A field function maps (x, y, z) in mm -> (Bx, By, Bz) in Tesla.
FieldFn = Callable[[float, float, float], Tuple[float, float, float]]


def linspace(a: float, b: float, n: int) -> List[float]:
    if n <= 1:
        return [0.5 * (a + b)]
    step = (b - a) / (n - 1)
    return [a + step * i for i in range(n)]


def _fmt(v: float) -> str:
    return f"{v:.8E}"


def write_3d(path: str, xs: Sequence[float], ys: Sequence[float],
             zs: Sequence[float], field_fn: FieldFn) -> None:
    """Write a 3D BDSIM magnetic field map (positions mm, field Tesla)."""
    lines = [
        f"xmin> {xs[0]:.6g}", f"xmax> {xs[-1]:.6g}", f"nx> {len(xs)}",
        f"ymin> {ys[0]:.6g}", f"ymax> {ys[-1]:.6g}", f"ny> {len(ys)}",
        f"zmin> {zs[0]:.6g}", f"zmax> {zs[-1]:.6g}", f"nz> {len(zs)}",
        "! X\tY\tZ\tFx\tFy\tFz",
    ]
    out = "\n".join(lines) + "\n"
    rows = []
    for x in xs:                     # X outermost
        for y in ys:                 # Y middle
            for z in zs:             # Z innermost (fastest)
                bx, by, bz = field_fn(x, y, z)
                rows.append(f"{_fmt(x)}\t{_fmt(y)}\t{_fmt(z)}\t"
                            f"{_fmt(bx)}\t{_fmt(by)}\t{_fmt(bz)}")
    return out + "\n".join(rows) + "\n"


def write_2d(path: str, xs: Sequence[float], ys: Sequence[float],
             field_fn: FieldFn, z: float = 0.0) -> str:
    """Write a 2D BDSIM magnetic field map over (x, y) at fixed z."""
    lines = [
        f"xmin> {xs[0]:.6g}", f"xmax> {xs[-1]:.6g}", f"nx> {len(xs)}",
        f"ymin> {ys[0]:.6g}", f"ymax> {ys[-1]:.6g}", f"ny> {len(ys)}",
        "! X\tY\tFx\tFy\tFz",
    ]
    out = "\n".join(lines) + "\n"
    rows = []
    for x in xs:                     # X outermost
        for y in ys:                 # Y fastest
            bx, by, bz = field_fn(x, y, z)
            rows.append(f"{_fmt(x)}\t{_fmt(y)}\t{_fmt(bx)}\t{_fmt(by)}\t{_fmt(bz)}")
    return out + "\n".join(rows) + "\n"


def build_3d(xs: Sequence[float], ys: Sequence[float], zs: Sequence[float],
             field_fn: FieldFn) -> str:
    """Return the 3D map text (positions mm, field Tesla)."""
    return write_3d("", xs, ys, zs, field_fn)


def gmad_field_object(name: str, filename: str, dims: int,
                      interpolator: str = "linear") -> str:
    """The GMAD ``field`` object declaration for a written map file."""
    return (f'{name}: field, type="bmap{dims}d", '
            f'magneticFile="bdsim{dims}d:{filename}", '
            f'magneticInterpolator="{interpolator}";')
