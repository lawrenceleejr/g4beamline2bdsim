"""Generate GDML geometry for G4beamline passive volumes (``box``, ``tubs``).

BDSIM cannot express an arbitrary material block in GMAD, but it *can* place
external GDML geometry as an ``element``.  Converting a G4beamline target /
absorber to GDML (instead of a drift) lets the material actually interact with
the beam in BDSIM -- restoring the energy loss, scattering and secondary
production that a drift would drop.

The GDML is written without an ``xs:noNamespaceSchemaLocation`` so it does not
depend on a (possibly offline) schema URL; Geant4 emits a non-fatal validation
note and reads the geometry regardless.
"""

from __future__ import annotations

from typing import Optional, Tuple

# G4beamline / common material names -> Geant4 NIST material names.
# z padding [mm] each side between the target solid and the GDML world box;
# the converter inflates the BDSIM element length by the same amount so the
# loaded geometry always fits inside the element.
WORLD_Z_PAD_MM = 0.1

_MATERIAL_MAP = {
    "vacuum": "G4_Galactic",
    "air": "G4_AIR",
    "cu": "G4_Cu", "copper": "G4_Cu",
    "fe": "G4_Fe", "iron": "G4_Fe",
    "al": "G4_Al", "aluminum": "G4_Al", "aluminium": "G4_Al",
    "w": "G4_W", "tungsten": "G4_W",
    "pb": "G4_Pb", "lead": "G4_Pb",
    "c": "G4_C", "carbon": "G4_C", "graphite": "G4_GRAPHITE",
    "be": "G4_Be", "beryllium": "G4_Be",
    "ti": "G4_Ti", "au": "G4_Au", "ag": "G4_Ag",
    "h2o": "G4_WATER", "water": "G4_WATER",
    "ss": "G4_STAINLESS-STEEL", "stainlesssteel": "G4_STAINLESS-STEEL",
    "kapton": "G4_KAPTON", "mylar": "G4_MYLAR",
    "scintillator": "G4_PLASTIC_SC_VINYLTOLUENE",
    "poly": "G4_POLYETHYLENE", "polyethylene": "G4_POLYETHYLENE",
    "concrete": "G4_CONCRETE", "si": "G4_Si", "silicon": "G4_Si",
}


def map_material(name: Optional[str]) -> Tuple[str, bool]:
    """Map a G4beamline material to a Geant4 name.

    Returns (geant4_name, recognised).  Unknown names are passed through with a
    ``G4_`` prefix (if not already present) and ``recognised=False`` so the
    caller can warn.
    """
    if not name:
        return "G4_Galactic", True
    raw = name.strip()
    key = raw.lower().lstrip("g4_")
    if raw.startswith("G4_"):
        return raw, True
    if raw.lower() in _MATERIAL_MAP:
        return _MATERIAL_MAP[raw.lower()], True
    if key in _MATERIAL_MAP:
        return _MATERIAL_MAP[key], True
    return "G4_" + raw, False


def is_vacuum(material: Optional[str]) -> bool:
    if not material:
        return True
    g4, _ = map_material(material)
    return g4 in ("G4_Galactic", "G4_vacuum")


def material_block(name: str, user_defs: dict) -> Tuple[str, str, list]:
    """Build a GDML ``<materials>`` block for a user-defined material.

    ``user_defs`` maps material name -> dict with optional keys ``Z``, ``A``,
    ``density`` [g/cm^3], ``state`` ('g'/'l'/'s'), and ``components`` (list of
    ``(ref_name, mass_fraction)``).  Components referencing other user-defined
    materials are included recursively; anything else is mapped to a Geant4
    NIST name (resolvable by Geant4's GDML reader).

    Returns ``(materials_xml, ref_name, warnings)`` -- ``materials_xml`` is
    empty when *name* is not user-defined (use the NIST mapping instead).
    """
    warnings: list = []
    if name not in user_defs:
        return "", name, warnings

    emitted: list = []
    seen = set()
    state_map = {"g": "gas", "l": "liquid", "s": "solid"}

    def emit(n: str) -> str:
        if n in seen:
            return n
        if n not in user_defs:
            g4, known = map_material(n)
            if not known:
                warnings.append(
                    f"material component '{n}' mapped to '{g4}'; verify it is "
                    f"a valid Geant4/NIST name")
            return g4
        seen.add(n)
        d = user_defs[n]
        density = d.get("density", 1.0)
        attrs = ""
        if d.get("state") in state_map:
            attrs += f' state="{state_map[d["state"]]}"'
        comps = d.get("components") or []
        if comps:
            # Mixture by mass fraction; emit components first.
            refs = [(emit(cn), cf) for cn, cf in comps]
            body = "".join(
                f'\n      <fraction n="{cf:.6g}" ref="{cn}"/>' for cn, cf in refs)
            emitted.append(
                f'    <material name="{n}"{attrs}>\n'
                f'      <D value="{density:.6g}"/>{body}\n    </material>')
        else:
            z = d.get("Z", 1.0)
            a = d.get("A", 1.0)
            emitted.append(
                f'    <material name="{n}" Z="{z:.6g}"{attrs}>\n'
                f'      <D value="{density:.6g}"/>\n'
                f'      <atom value="{a:.6g}"/>\n    </material>')
        return n

    ref = emit(name)
    xml = "  <materials>\n" + "\n".join(emitted) + "\n  </materials>\n"
    return xml, ref, warnings


def _gdml_document(solids: str, target_solid: str, material: str,
                   world_x: float, world_y: float, world_z: float,
                   materials_xml: str = "") -> str:
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<gdml>
{materials_xml}  <solids>
{solids}
    <box name="world_solid" x="{world_x:.6g}" y="{world_y:.6g}" z="{world_z:.6g}" lunit="mm"/>
  </solids>
  <structure>
    <volume name="target_vol">
      <materialref ref="{material}"/>
      <solidref ref="{target_solid}"/>
    </volume>
    <volume name="world_vol">
      <materialref ref="G4_Galactic"/>
      <solidref ref="world_solid"/>
      <physvol>
        <volumeref ref="target_vol"/>
      </physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="world_vol"/>
  </setup>
</gdml>
"""


def box_gdml(width_mm: float, height_mm: float, length_mm: float,
             material: str, materials_xml: str = "") -> str:
    """GDML for a box (full dimensions in mm), centred on the origin."""
    solid = (f'    <box name="target_solid" x="{width_mm:.6g}" '
             f'y="{height_mm:.6g}" z="{length_mm:.6g}" lunit="mm"/>')
    wx = max(width_mm, 1.0) * 1.2 + 20.0
    wy = max(height_mm, 1.0) * 1.2 + 20.0
    wz = max(length_mm, 1.0) + 2.0 * WORLD_Z_PAD_MM
    return _gdml_document(solid, "target_solid", material, wx, wy, wz,
                          materials_xml)


def sphere_gdml(inner_r_mm: float, outer_r_mm: float, material: str,
                start_phi_deg: float = 0.0, delta_phi_deg: float = 360.0,
                start_theta_deg: float = 0.0, delta_theta_deg: float = 180.0,
                materials_xml: str = "") -> str:
    """GDML for a sphere/spherical shell (mm, degrees), centred on the origin."""
    solid = (f'    <sphere name="target_solid" rmin="{inner_r_mm:.6g}" '
             f'rmax="{outer_r_mm:.6g}" startphi="{start_phi_deg:.6g}" '
             f'deltaphi="{delta_phi_deg:.6g}" starttheta="{start_theta_deg:.6g}" '
             f'deltatheta="{delta_theta_deg:.6g}" aunit="deg" lunit="mm"/>')
    w = max(2.0 * outer_r_mm, 1.0) * 1.2 + 20.0
    wz = max(2.0 * outer_r_mm, 1.0) + 2.0 * WORLD_Z_PAD_MM
    return _gdml_document(solid, "target_solid", material, w, w, wz,
                          materials_xml)


def polycone_gdml(zs_mm, rin_mm, rout_mm, material: str,
                  start_phi_deg: float = 0.0, delta_phi_deg: float = 360.0,
                  materials_xml: str = "") -> str:
    """GDML for a polycone.  z positions are recentred about the origin."""
    zc = 0.5 * (min(zs_mm) + max(zs_mm))
    planes = "\n".join(
        f'      <zplane z="{z - zc:.6g}" rmin="{ri:.6g}" rmax="{ro:.6g}"/>'
        for z, ri, ro in zip(zs_mm, rin_mm, rout_mm))
    solid = (f'    <polycone name="target_solid" startphi="{start_phi_deg:.6g}" '
             f'deltaphi="{delta_phi_deg:.6g}" aunit="deg" lunit="mm">\n'
             f'{planes}\n    </polycone>')
    rmax = max(rout_mm) if rout_mm else 1.0
    zext = max(zs_mm) - min(zs_mm) if len(zs_mm) > 1 else 1.0
    w = max(2.0 * rmax, 1.0) * 1.2 + 20.0
    wz = max(zext, 1.0) + 2.0 * WORLD_Z_PAD_MM
    return _gdml_document(solid, "target_solid", material, w, w, wz,
                          materials_xml)


def tubs_gdml(inner_r_mm: float, outer_r_mm: float, length_mm: float,
              material: str, start_phi_deg: float = 0.0,
              delta_phi_deg: float = 360.0, materials_xml: str = "") -> str:
    """GDML for a tube/cylinder (mm, degrees), centred on the origin."""
    solid = (f'    <tube name="target_solid" rmin="{inner_r_mm:.6g}" '
             f'rmax="{outer_r_mm:.6g}" z="{length_mm:.6g}" '
             f'startphi="{start_phi_deg:.6g}" deltaphi="{delta_phi_deg:.6g}" '
             f'aunit="deg" lunit="mm"/>')
    w = max(2.0 * outer_r_mm, 1.0) * 1.2 + 20.0
    wz = max(length_mm, 1.0) + 2.0 * WORLD_Z_PAD_MM
    return _gdml_document(solid, "target_solid", material, w, w, wz,
                          materials_xml)
