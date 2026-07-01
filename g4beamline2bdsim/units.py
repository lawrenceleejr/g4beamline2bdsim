"""Unit handling and conversion between G4beamline and BDSIM conventions.

G4beamline / Geant4 base units
------------------------------
* length         : millimetre (mm)
* energy/momentum: mega-electron-volt (MeV)
* magnetic field : tesla (T)
* angle          : degree for the user-facing rotation arguments of ``place``
                   and the half-angle args of ``genericbend`` are radians in
                   some places; the converter is explicit at each call site.

BDSIM / GMAD base units
-----------------------
* length : metre (m)
* energy : giga-electron-volt (GeV)
* angle  : radian (rad)
* field  : tesla (T)

This module only provides simple scalar helpers; the physics-aware conversions
(field gradient -> normalised k1, etc.) live in :mod:`converter` because they
depend on the beam rigidity.
"""

from __future__ import annotations

# Length: G4beamline works in mm, BDSIM in m.
MM_TO_M = 1.0e-3
M_TO_MM = 1.0e3

# Energy: G4beamline works in MeV, BDSIM in GeV.
MEV_TO_GEV = 1.0e-3
GEV_TO_MEV = 1.0e3

# Speed of light (m/s), matching GMAD's clight constant.
C_LIGHT = 2.99792458e8


def mm_to_m(value: float) -> float:
    return value * MM_TO_M


def mev_to_gev(value: float) -> float:
    return value * MEV_TO_GEV


def brho_from_momentum_mev(momentum_mev: float, charge: float = 1.0) -> float:
    """Magnetic rigidity B*rho [T*m] for a particle of given momentum.

    For a singly charged particle ``Brho [T m] = p [GeV/c] / 0.299792458``.
    More generally ``Brho = p / (0.299792458 * |q|)`` with ``p`` in GeV/c and
    ``q`` in units of the elementary charge.
    """
    p_gev = momentum_mev * MEV_TO_GEV
    if charge == 0:
        return 0.0
    return p_gev / (0.299792458 * abs(charge))
