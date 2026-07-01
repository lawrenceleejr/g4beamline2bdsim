"""Realistic solenoid field from a coil + current density (G4beamline-style).

G4beamline models a solenoid as its physical ``coil`` (inner radius, outer
radius, length) carrying a current density, and computes the field by summing
current sheets.  This module reproduces that: it sums circular current loops
(radial * axial subdivisions) using the exact single-loop field expressed with
complete elliptic integrals, giving the true on- and off-axis field including
the end/fringe fall-off -- far better than BDSIM's ideal hard-edge ``ks``.

The result is sampled onto a grid and written as a BDSIM field map
(:mod:`~g4beamline2bdsim.fieldmap`), so BDSIM tracks through the real field.

Units: lengths in mm internally converted to m; current density in A/mm^2;
field returned in Tesla.
"""

from __future__ import annotations

import math
from typing import Callable, List, Tuple

MU0 = 4.0e-7 * math.pi          # T*m / A


def _agm_ke(m: float) -> Tuple[float, float]:
    """Complete elliptic integrals K(m), E(m) (parameter m = k^2) via the AGM.

    Valid for 0 <= m < 1.
    """
    if m < 0:
        m = 0.0
    if m >= 1.0:
        m = 1.0 - 1e-12
    # E(m)/K(m) = 1 - sum_{n=0}^inf 2^{n-1} c_n^2, with c_0^2 = m.
    a = 1.0
    b = math.sqrt(1.0 - m)
    s = 0.5 * m            # n=0 term: 2^{-1} c_0^2 = m/2
    p = 1.0                # 2^{n-1} for n=1 is 2^0 = 1
    for _ in range(60):
        an = 0.5 * (a + b)
        bn = math.sqrt(a * b)
        cn = 0.5 * (a - b)
        s += p * cn * cn
        p *= 2.0
        a, b = an, bn
        if abs(cn) < 1e-15:
            break
    K = math.pi / (2.0 * a)
    E = K * (1.0 - s)
    return K, E


def _loop_field(a: float, r: float, z: float, current: float) -> Tuple[float, float]:
    """Field (Br, Bz) [T] of a single circular loop of radius ``a`` [m] carrying
    ``current`` [A], at cylindrical (r, z) [m].  Axis along z, loop at z=0.
    """
    if a <= 0.0:
        return 0.0, 0.0
    # On-axis: closed form, avoids r=0 singularity.
    if r < 1e-9:
        denom = (a * a + z * z) ** 1.5
        bz = MU0 * current * a * a / (2.0 * denom) if denom > 0 else 0.0
        return 0.0, bz
    q = (a + r) ** 2 + z * z
    sq = math.sqrt(q)
    m = 4.0 * a * r / q                     # = k^2
    K, E = _agm_ke(m)
    d = (a - r) ** 2 + z * z
    if d < 1e-18:
        d = 1e-18
    c = MU0 * current / (2.0 * math.pi)
    bz = c / sq * (K + (a * a - r * r - z * z) / d * E)
    br = c * z / (r * sq) * (-K + (a * a + r * r + z * z) / d * E)
    return br, bz


def make_coil_field(inner_r_mm: float, outer_r_mm: float, length_mm: float,
                    current_A_per_mm2: float, n_radial: int = 8,
                    n_axial: int = 0) -> Callable[[float, float], Tuple[float, float]]:
    """Return ``B(r_mm, z_mm) -> (Br, Bz)`` [T] for a thick finite solenoid.

    The coil (a1..a2, length L, centred at z=0) is divided into
    ``n_radial * n_axial`` current loops carrying the current density.
    """
    a1 = inner_r_mm * 1e-3
    a2 = outer_r_mm * 1e-3
    L = length_mm * 1e-3
    J = current_A_per_mm2 * 1e6              # A/mm^2 -> A/m^2
    if n_axial <= 0:
        # Resolve the length reasonably; more sheets over a longer coil.
        n_axial = max(20, min(120, int(L / max(a2 - a1, 1e-3)) * 4 + 20))
    da = (a2 - a1) / n_radial
    dz = L / n_axial
    # Loop radii (radial layer centres) and axial positions.
    radii = [a1 + da * (i + 0.5) for i in range(n_radial)]
    zpos = [-0.5 * L + dz * (j + 0.5) for j in range(n_axial)]
    # Current per loop = J * (cross-sectional area element).
    dI = J * da * dz

    def field(r_mm: float, z_mm: float) -> Tuple[float, float]:
        r = r_mm * 1e-3
        z = z_mm * 1e-3
        br = 0.0
        bz = 0.0
        for a in radii:
            for zc in zpos:
                lbr, lbz = _loop_field(a, r, z - zc, dI)
                br += lbr
                bz += lbz
        return br, bz

    return field


def central_field(inner_r_mm: float, outer_r_mm: float, length_mm: float,
                  current_A_per_mm2: float) -> float:
    """Central on-axis field [T] (thick-solenoid closed form)."""
    a1 = inner_r_mm
    a2 = outer_r_mm
    if a1 <= 0 or a2 <= a1 or length_mm <= 0:
        return 0.0
    J = current_A_per_mm2 * 1e6
    alpha = a2 / a1
    beta = (length_mm / 2.0) / a1
    f = beta * math.log((alpha + math.sqrt(alpha ** 2 + beta ** 2))
                        / (1.0 + math.sqrt(1.0 + beta ** 2)))
    return MU0 * J * (a1 * 1e-3) * f


def sampled_3d_map(field_rz, half_x_mm: float, half_y_mm: float,
                   z_min_mm: float, z_max_mm: float,
                   nx: int, ny: int, nz: int):
    """Build a Cartesian field function from an axisymmetric ``B(r,z)``.

    Returns ``(field_fn, xs, ys, zs)`` ready for :func:`fieldmap.write_3d`.
    """
    from .fieldmap import linspace

    xs = linspace(-half_x_mm, half_x_mm, nx)
    ys = linspace(-half_y_mm, half_y_mm, ny)
    zs = linspace(z_min_mm, z_max_mm, nz)

    # Precompute a (r, z) table and bilinearly interpolate for speed.
    r_max = math.hypot(half_x_mm, half_y_mm)
    nr = max(nx, ny) + 4
    rs = linspace(0.0, r_max, nr)
    table: List[List[Tuple[float, float]]] = [
        [field_rz(r, z) for z in zs] for r in rs
    ]
    dr = rs[1] - rs[0] if nr > 1 else 1.0

    def lookup(r: float, z_index: int) -> Tuple[float, float]:
        if r <= 0:
            return table[0][z_index]
        fr = r / dr
        i = int(fr)
        if i >= nr - 1:
            return table[nr - 1][z_index]
        t = fr - i
        br = table[i][z_index][0] * (1 - t) + table[i + 1][z_index][0] * t
        bz = table[i][z_index][1] * (1 - t) + table[i + 1][z_index][1] * t
        return br, bz

    z_index_of = {round(z, 6): k for k, z in enumerate(zs)}

    def field_fn(x: float, y: float, z: float) -> Tuple[float, float, float]:
        r = math.hypot(x, y)
        k = z_index_of[round(z, 6)]
        br, bz = lookup(r, k)
        if r > 1e-9:
            return br * x / r, br * y / r, bz
        return 0.0, 0.0, bz

    return field_fn, xs, ys, zs
