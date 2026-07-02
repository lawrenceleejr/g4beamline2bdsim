"""Convert parsed G4beamline commands into a BDSIM (GMAD) model.

Conversion strategy
-------------------
G4beamline describes geometry by *defining* elements and then *placing* them at
absolute ``z`` positions (or sequentially) in centreline coordinates.  BDSIM
instead describes a *sequence* of elements whose lengths tile the beamline, with
``drift`` sections filling the gaps.

The converter therefore:

1. Records every element definition and every ``place`` command.
2. Computes each placement's longitudinal extent ``[z_entry, z_exit]`` from the
   element's physical length and placement ``z`` (handling ``front=1`` and the
   sector-bend "front-face" convention).
3. Sorts placements by ``z_entry`` and inserts ``drift`` elements to fill gaps.
4. Emits GMAD element definitions, the ``line``, ``beam``, ``option`` and
   ``sample`` commands.

Units
-----
Lengths/positions/apertures are kept in **millimetres** and emitted with the
``*mm`` suffix so the exact G4beamline numbers are preserved.  Magnet field
strengths are converted to BDSIM's normalised ``k`` values using the beam
rigidity derived from the reference/beam momentum (see :mod:`units`).
"""

from __future__ import annotations

import math
import os
import re
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from . import blfieldmap as _blfieldmap
from . import fieldmap as _fieldmap
from . import gdml as _gdml
from . import solenoid as _solenoid
from .expr import evaluate as _eval_expr
from .model import BdsimModel, Element
from .parser import G4BLCommand
from .units import C_LIGHT

# ---------------------------------------------------------------------------
# Particle data: charge (in e) for rigidity, and name mapping G4bl -> BDSIM.
# ---------------------------------------------------------------------------
_PARTICLE_CHARGE: Dict[str, float] = {
    "proton": 1.0,
    "anti_proton": -1.0,
    "e-": -1.0,
    "e+": 1.0,
    "mu-": -1.0,
    "mu+": 1.0,
    "pi+": 1.0,
    "pi-": -1.0,
    "kaon+": 1.0,
    "kaon-": -1.0,
    "gamma": 0.0,
    "neutron": 0.0,
}

# G4beamline particle name -> BDSIM particle name (mostly identical).
_PARTICLE_NAME: Dict[str, str] = {
    "proton": "proton",
    "anti_proton": "anti_proton",
    "e-": "e-",
    "e+": "e+",
    "mu-": "mu-",
    "mu+": "mu+",
    "pi+": "pi+",
    "pi-": "pi-",
    "kaon+": "kaon+",
    "kaon-": "kaon-",
    "gamma": "gamma",
    "neutron": "neutron",
}

# A tiny tolerance (mm) below which gaps/overlaps are ignored.
_Z_TOL = 1e-6

# GMAD reserved words we must not collide with when naming elements.
_GMAD_RESERVED = {
    "drift", "sbend", "rbend", "quadrupole", "sextupole", "octupole",
    "decapole", "multipole", "solenoid", "rcol", "ecol", "rfcavity", "rf",
    "marker", "line", "use", "beam", "option", "sample", "field",
}


def _to_float(value: Optional[str], default: float = 0.0) -> float:
    if value is None:
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        # Allow numeric expressions that survived parameter expansion, evaluated
        # safely (no arbitrary code execution).
        result = _eval_expr(value)
        return result if result is not None else default


def _parse_float_list(value: Optional[str]) -> List[float]:
    """Parse a comma-separated numeric list (G4beamline list argument)."""
    if not value:
        return []
    out = []
    for tok in value.replace(";", ",").split(","):
        tok = tok.strip()
        if tok:
            out.append(_to_float(tok))
    return out


def _sanitize(name: str) -> str:
    """Make a valid, non-reserved GMAD identifier from a G4beamline name."""
    cleaned = re.sub(r"[^A-Za-z0-9_]", "_", name)
    if not cleaned or not cleaned[0].isalpha():
        cleaned = "e_" + cleaned
    if cleaned.lower() in _GMAD_RESERVED:
        cleaned = cleaned + "_"
    return cleaned


@dataclass
class _Placement:
    """A single placed instance along the beamline."""

    name: str            # final (unique) name
    element: Element     # the BDSIM element to emit
    length_mm: float     # physical length (mm)
    z_entry_mm: float    # longitudinal entry position (mm)
    is_sampler: bool = False


class Converter:
    """Convert a list of :class:`G4BLCommand` into a :class:`BdsimModel`."""

    def __init__(self, commands: List[G4BLCommand], source_name: str = "",
                 emit_gdml: bool = True, emit_field_maps: bool = True,
                 solenoid_field_map: bool = True, base_dir: str = "") -> None:
        self.commands = commands
        self.source_name = source_name
        # Directory for resolving external files referenced by the input (e.g.
        # a BLFieldMap file named by a ``fieldmap`` command).
        self.base_dir = base_dir
        # When True, passive material volumes (box/tubs) become BDSIM elements
        # with external GDML geometry (material interacts); otherwise drifts.
        self.emit_gdml = emit_gdml
        # When True, fieldexpr formulas are converted to BDSIM field maps.
        self.emit_field_maps = emit_field_maps
        # When True (default, since we always convert from G4beamline), a
        # solenoid becomes a full 3D field map of the G4beamline coil field --
        # the faithful "g4bl-style" field including the end fringe.  Set False
        # for a native BDSIM solenoid with ks from the peak field (hard-edge,
        # more robust for very strong solenoids).
        self.solenoid_field_map = solenoid_field_map
        self.model = BdsimModel()

        # element-definition name -> (G4BLCommand, type)
        self._defs: Dict[str, G4BLCommand] = {}
        self._coils: Dict[str, G4BLCommand] = {}
        # User-defined materials (g4bl 'material' command) for GDML export.
        self._user_materials: Dict[str, dict] = {}
        self._placements: List[_Placement] = []
        self._used_names: Dict[str, int] = {}
        self._place_counts: Dict[str, int] = {}
        # Largest transverse aperture/detector radius seen [mm]; used to size the
        # BDSIM beampipe so a wide (vacuum) G4beamline beam is not clipped.
        self._max_aperture_mm: float = 0.0

        # Beam / rigidity state.
        self._momentum_mev: Optional[float] = None
        self._particle: str = "mu+"           # G4beamline default
        self._charge: float = 1.0
        self._n_events: Optional[int] = None
        self._beam_cmd: Optional[G4BLCommand] = None
        self._physics: Optional[str] = None

        # Running z for sequential (no-z) placement.
        self._z_cursor_mm: float = 0.0

    # -- public API ---------------------------------------------------------
    def convert(self) -> BdsimModel:
        if self.source_name:
            self.model.header_comments.append(f"Source: {self.source_name}")

        # First pass: collect definitions, beam/physics, and placements.
        for cmd in self.commands:
            self._dispatch(cmd)

        # Resolve rigidity now that beam/reference are known.
        self._finalize_beam()

        # Convert genericquad/multipole field strengths to k-values, which need
        # the rigidity -- done lazily here by revisiting placements.
        self._apply_rigidity()

        # Build the ordered beamline with drifts.
        self._build_beamline()

        # Beam, options, samplers.
        self._build_beam_block()
        self._build_options_block()

        return self.model

    # -- dispatch -----------------------------------------------------------
    def _dispatch(self, cmd: G4BLCommand) -> None:
        name = cmd.name
        if name in _DEFINITION_TYPES:
            self._record_definition(cmd)
        elif name == "coil":
            if cmd.args:
                self._coils[cmd.args[0]] = cmd
        elif name == "material":
            self._record_material(cmd)
        elif name == "place":
            self._handle_place(cmd)
        elif name == "beam":
            self._beam_cmd = cmd
            self._handle_beam(cmd)
        elif name == "reference":
            self._handle_reference(cmd)
        elif name == "physics":
            if cmd.args:
                self._physics = cmd.args[0]
        elif name in ("corner", "cornerarc", "start", "group", "endgroup"):
            # Geometry/centreline helpers: the bend magnets themselves carry the
            # deflection in BDSIM, so these are informational only.
            if name in ("corner", "cornerarc"):
                self.model.warn(
                    f"'{name}' (line {cmd.line_no}) ignored: bending is taken "
                    f"from the magnet field/angle in BDSIM."
                )
        elif name in _UNSUPPORTED_COMMANDS:
            self.model.warn(
                f"'{name}' (line {cmd.line_no}) not converted: "
                f"{_UNSUPPORTED_COMMANDS[name]}"
            )
        # Everything else (g4ui, trace, param, tune, output, ...) is ignored.

    def _record_material(self, cmd: G4BLCommand) -> None:
        """Record a g4bl 'material' definition (element or mass-fraction mix)."""
        if not cmd.args:
            return
        name = cmd.args[0]
        d: dict = {}
        z = cmd.get("Z") or cmd.get("z")
        a = cmd.get("A") or cmd.get("a") or cmd.get("aaa")
        if z is not None:
            d["Z"] = _to_float(z)
        if a is not None:
            d["A"] = _to_float(a)
        if cmd.get("density") is not None:
            d["density"] = _to_float(cmd.get("density"))
        if cmd.get("state"):
            d["state"] = cmd.get("state")
        comps = []
        for tok in cmd.args[1:]:
            if "," in tok:
                cname, _, frac = tok.partition(",")
                comps.append((cname.strip(), _to_float(frac)))
        if comps:
            d["components"] = comps
        self._user_materials[name] = d

    def _material_for_gdml(self, material: Optional[str]):
        """Resolve a material for GDML export.

        Returns (ref_name, materials_xml, known): user-defined materials get a
        <materials> block; everything else maps to a Geant4/NIST name.
        """
        if material and material in self._user_materials:
            xml, ref, warns = _gdml.material_block(material, self._user_materials)
            for w in warns:
                self.model.warn(w)
            return ref, xml, True
        g4mat, known = _gdml.map_material(material)
        return g4mat, "", known

    def _record_definition(self, cmd: G4BLCommand) -> None:
        if not cmd.args:
            return
        self._defs[cmd.args[0]] = cmd

    # -- beam / reference / rigidity ---------------------------------------
    def _handle_beam(self, cmd: G4BLCommand) -> None:
        particle = cmd.get("particle")
        if particle:
            self._particle = particle
        p = (
            cmd.get("meanMomentum")
            or cmd.get("meanP")
            or cmd.get("P")
        )
        if p is not None:
            self._momentum_mev = _to_float(p)
        n = cmd.get("nEvents")
        if n is not None:
            self._n_events = int(_to_float(n))

    def _handle_reference(self, cmd: G4BLCommand) -> None:
        particle = cmd.get("particle")
        if particle:
            self._particle = particle
        p = (
            cmd.get("referenceMomentum")
            or cmd.get("meanMomentum")
            or cmd.get("P")
        )
        if p is not None and self._momentum_mev is None:
            self._momentum_mev = _to_float(p)

    def _finalize_beam(self) -> None:
        self._charge = _PARTICLE_CHARGE.get(self._particle, 1.0)
        if self._momentum_mev is None:
            self.model.warn(
                "No reference/beam momentum found; magnet k-values cannot be "
                "computed. Field strengths left as comments. Set the beam "
                "momentum and re-run, or fill in k-values by hand."
            )

    def _brho(self) -> Optional[float]:
        """Rigidity magnitude B*rho [T*m], or None if unknown/uncharged."""
        if self._momentum_mev is None:
            return None
        q = self._charge if self._charge != 0 else 1.0
        p_gev = self._momentum_mev * 1.0e-3
        return p_gev / (0.299792458 * abs(q))

    # -- placement bookkeeping ---------------------------------------------
    def _handle_place(self, cmd: G4BLCommand) -> None:
        if not cmd.args:
            return
        def_name = cmd.args[0]
        definition = self._defs.get(def_name)
        if definition is None:
            self.model.warn(
                f"place of undefined element '{def_name}' (line {cmd.line_no}) "
                f"skipped."
            )
            return

        if cmd.get("parent") is not None:
            self.model.warn(
                f"'{def_name}' (line {cmd.line_no}) placed inside a parent "
                f"volume; nested geometry is not converted to the beamline."
            )
            return

        # Determine final name.  G4beamline replaces '#' in a rename with the
        # (1-based) placement number of that element definition, and a leading
        # '+' prepends the parent name (no parent here, so just drop it).
        rename = cmd.get("rename")
        self._place_counts[def_name] = self._place_counts.get(def_name, 0) + 1
        if rename:
            base = rename.lstrip("+")
            if "#" in base:
                base = base.replace("#", str(self._place_counts[def_name]))
        else:
            base = def_name
        final_name = self._unique_name(_sanitize(base))

        length_mm = self._element_length_mm(definition)
        self._record_aperture(definition)

        # Convert the element (field strengths carried as raw values, resolved
        # later once rigidity is known).
        element = self._build_element(final_name, definition, cmd, length_mm)
        if element is None:
            return
        self._apply_transverse_placement(element, cmd)

        # Longitudinal position.
        z_entry = self._placement_z_entry(cmd, definition, length_mm)
        self._z_cursor_mm = z_entry + length_mm

        is_sampler = definition.name in ("virtualdetector", "detector")
        self._placements.append(
            _Placement(
                name=final_name,
                element=element,
                length_mm=length_mm,
                z_entry_mm=z_entry,
                is_sampler=is_sampler,
            )
        )

    # Matches one rotation token like Y90 or Z-30.5 (g4bl rotation syntax).
    _ROT_TOKEN = re.compile(r"^([XYZ])(-?[0-9.eE+]+)$")

    def _apply_transverse_placement(self, element: Element,
                                    cmd: G4BLCommand) -> None:
        """Map ``place x=/y=/rotation=`` onto the BDSIM element.

        * quadrupole offsets become feed-down steering (an offset quad is an
          on-axis quad plus a dipole kick); BDSIM's ``offsetX`` does NOT move
          the tracking field, so thin kickers are injected instead (resolved
          in :meth:`_apply_rigidity` once k1 is known).
        * GDML ``element`` offsets displace the geometry -- exactly right for
          a passive target -- so ``offsetX``/``offsetY`` are used directly.
        * a rotation about the beam axis (``Z<deg>``) maps to ``tilt`` [rad];
          X/Y rotations have no 1-D equivalent and are warned (for bends the
          ``Y(angle/2)`` convention accompanies a ``corner`` and is expected).
        """
        if element.type == "marker":
            return                      # samplers cannot carry offsets
        x = _to_float(cmd.get("x"))
        y = _to_float(cmd.get("y"))
        if x != 0.0 or y != 0.0:
            if element.type == "quadrupole":
                element.params["__offset_dx_mm"] = x
                element.params["__offset_dy_mm"] = y
            elif element.type == "element":
                if x != 0.0:
                    element.set("offsetX", (x, "mm"))
                if y != 0.0:
                    element.set("offsetY", (y, "mm"))
            else:
                self.model.warn(
                    f"'{element.name}': transverse placement offset "
                    f"(x={x:g}, y={y:g} mm) on a '{element.type}' is not "
                    f"mapped (only quadrupole feed-down and GDML geometry "
                    f"offsets are supported)."
                )
        rotation = cmd.get("rotation")
        if not rotation:
            return
        tilt = 0.0
        unmapped = []
        for token in rotation.split(","):
            token = token.strip()
            if not token:
                continue
            m = self._ROT_TOKEN.match(token)
            if m and m.group(1) == "Z":
                tilt += math.radians(_to_float(m.group(2)))
            else:
                unmapped.append(token)
        if tilt != 0.0:
            element.set("tilt", tilt)
        if unmapped:
            self.model.warn(
                f"'{element.name}': rotation component(s) "
                f"{','.join(unmapped)} about X/Y cannot be mapped to a 1-D "
                f"beamline element; ignored."
            )

    def _placement_z_entry(
        self, cmd: G4BLCommand, definition: G4BLCommand, length_mm: float
    ) -> float:
        z = cmd.get("z")
        front = _to_float(cmd.get("front"), 0.0) != 0.0
        sector = definition.name in ("idealsectorbend", "genericsectorbend")
        if z is None:
            # Sequential placement: just downstream of the previous element.
            return self._z_cursor_mm
        z_val = _to_float(z)
        if front or sector:
            # 'z' marks the front face.
            return z_val
        # 'z' marks the geometric centre.
        return z_val - length_mm / 2.0

    def _unique_name(self, base: str) -> str:
        count = self._used_names.get(base, 0)
        self._used_names[base] = count + 1
        if count == 0:
            return base
        return f"{base}_{count}"

    def _record_aperture(self, definition: G4BLCommand) -> None:
        """Track the largest transverse *aperture* radius referenced.

        Only genuine transverse apertures are used -- NOT yoke/coil outer radii
        (``ironRadius``, ``outerRadius``) or sector-bend arc radii
        (``fieldInnerRadius``/``fieldCenterRadius``/``fieldOuterRadius``), which
        describe curvature, not the beam aperture.
        """
        g = definition.get
        for key in ("apertureRadius", "radius"):
            v = g(key)
            if v is not None:
                self._max_aperture_mm = max(self._max_aperture_mm, _to_float(v))
        # Half-extents of rectangular field/box apertures.
        for key in ("fieldHeight", "fieldWidth", "height", "width"):
            v = g(key)
            if v is not None:
                self._max_aperture_mm = max(self._max_aperture_mm,
                                            _to_float(v) / 2.0)

    # -- per-element length -------------------------------------------------
    def _element_length_mm(self, definition: G4BLCommand) -> float:
        t = definition.name
        g = definition.get
        if t == "genericbend":
            return _to_float(g("fieldLength"))
        if t == "genericquad":
            return _to_float(g("fieldLength"))
        if t == "multipole":
            return _to_float(g("fieldLength"))
        if t in ("idealsectorbend", "genericsectorbend"):
            angle_deg = _to_float(g("angle"))
            radius = _to_float(g("fieldCenterRadius"))
            return abs(math.radians(angle_deg)) * radius
        if t == "solenoid":
            coil_name = g("coilName") or g("coil")
            coil = self._coils.get(coil_name) if coil_name else None
            if coil is not None:
                clen = _to_float(coil.get("length"))
                if self.solenoid_field_map:
                    # Extend the field region to capture the end fringe, which
                    # carries the solenoid's focusing (radial-field impulse).
                    return clen + 2.0 * self._solenoid_margin(coil)
                return clen
            return 0.0
        if t in ("pillbox", "rfdevice"):
            return _to_float(g("innerLength"))
        if t in ("tubs", "cylinder", "box"):
            return self._gdml_padded(definition, _to_float(g("length")))
        if t == "fieldexpr":
            return _to_float(g("length"))
        if t == "sphere":
            return self._gdml_padded(definition, 2.0 * _to_float(g("outerRadius")))
        if t == "polycone":
            zs = _parse_float_list(g("z"))
            span = (max(zs) - min(zs)) if len(zs) > 1 else 0.0
            return self._gdml_padded(definition, span)
        if t == "fieldmap":
            grid = self._parse_blfieldmap(definition)
            if grid and len(grid.zs) > 1:
                return grid.zs[-1] - grid.zs[0]
            return 0.0
        if t in ("virtualdetector", "detector"):
            return _to_float(g("length"), 1.0)
        return 0.0

    def _gdml_padded(self, definition: G4BLCommand, solid_len: float) -> float:
        """Element length for a passive solid: padded when it becomes GDML.

        The GDML world box is the solid plus ``WORLD_Z_PAD_MM`` each side, and
        BDSIM requires the loaded geometry to fit inside the element length, so
        the layout length must carry the same padding (and a 1 mm floor for
        very thin solids, matching the world-box floor).
        """
        material = definition.get("material")
        if self.emit_gdml and not _gdml.is_vacuum(material):
            return max(solid_len, 1.0) + 2.0 * _gdml.WORLD_Z_PAD_MM
        return solid_len

    # -- element construction ----------------------------------------------
    def _build_element(
        self,
        name: str,
        definition: G4BLCommand,
        place: G4BLCommand,
        length_mm: float,
    ) -> Optional[Element]:
        t = definition.name
        builder = getattr(self, f"_conv_{t}", None)
        if builder is None:
            self.model.warn(
                f"element type '{t}' (placed as '{name}') is not supported; "
                f"represented as a drift of the same length."
            )
            return self._drift(name, length_mm)
        return builder(name, definition, place, length_mm)

    def _drift(self, name: str, length_mm: float) -> Element:
        el = Element(name=name, type="drift")
        el.set("l", (length_mm, "mm"))
        return el

    # field strengths are stashed on the element under private keys and
    # resolved in _apply_rigidity once Brho is known.
    def _conv_genericquad(self, name, definition, place, length_mm) -> Element:
        el = Element(name=name, type="quadrupole")
        el.set("l", (length_mm, "mm"))
        ap = _to_float(definition.get("apertureRadius"))
        if ap > 0:
            el.set("apertureType", "circular")
            el.set("aper1", (ap, "mm"))
        gradient = _to_float(self._tuned(definition, place, "gradient"))
        el.params["__gradient_T_per_m"] = gradient
        return el

    def _conv_genericbend(self, name, definition, place, length_mm) -> Element:
        # Rectangular box field -> rbend, field given directly in Tesla.  For
        # strong bends the implied angle B*L/Brho can exceed what BDSIM's rbend
        # geometry supports (pole faces overlap), so the final rbend-vs-sbend
        # choice is made in _resolve_strengths once the rigidity is known.
        el = Element(name=name, type="rbend")
        el.set("l", (length_mm, "mm"))
        by = _to_float(self._tuned(definition, place, "By"))
        el.set("B", (by, "T"))
        el.params["__bend_check"] = True
        el.params["__bend_fw"] = _to_float(definition.get("fieldWidth"))
        el.params["__bend_fh"] = _to_float(definition.get("fieldHeight"))
        return el

    def _conv_idealsectorbend(self, name, definition, place, length_mm) -> Element:
        el = Element(name=name, type="sbend")
        el.set("l", (length_mm, "mm"))
        angle_deg = _to_float(self._tuned(definition, place, "angle"))
        if angle_deg != 0.0:
            el.set("angle", math.radians(angle_deg))
        else:
            by = _to_float(self._tuned(definition, place, "By"))
            el.set("B", (by, "T"))
        return el

    def _conv_genericsectorbend(self, name, definition, place, length_mm) -> Element:
        el = Element(name=name, type="sbend")
        el.set("l", (length_mm, "mm"))
        angle_deg = _to_float(self._tuned(definition, place, "angle"))
        if angle_deg != 0.0:
            el.set("angle", math.radians(angle_deg))
        dip = _to_float(definition.get("DipoleField"))
        if angle_deg == 0.0 and dip != 0.0:
            el.set("B", (dip, "T"))
        quad = _to_float(definition.get("QuadrupoleField"))
        if quad != 0.0:
            el.params["__gradient_T_per_m"] = quad
        return el

    def _conv_multipole(self, name, definition, place, length_mm) -> Element:
        g = definition.get
        dipole = _to_float(g("dipole"))
        quad = _to_float(g("quadrupole"))
        sext = _to_float(g("sextupole"))
        octo = _to_float(g("octopole") or g("octupole"))
        deca = _to_float(g("decapole"))
        # Choose the dominant single-order representation; combined functions
        # are represented as a 'multipole' with knl handled in rigidity pass.
        nonzero = [v for v in (quad, sext, octo, deca) if v != 0.0]
        if dipole != 0.0 and not nonzero:
            el = Element(name=name, type="rbend")
            el.set("l", (length_mm, "mm"))
            el.set("B", (dipole, "T"))
            el.params["__bend_check"] = True
            return el
        if len(nonzero) == 1 and quad != 0.0:
            el = Element(name=name, type="quadrupole")
            el.set("l", (length_mm, "mm"))
            el.params["__gradient_T_per_m"] = quad
        elif len(nonzero) == 1 and sext != 0.0:
            el = Element(name=name, type="sextupole")
            el.set("l", (length_mm, "mm"))
            el.params["__sext_T_per_m2"] = sext
        elif len(nonzero) == 1 and octo != 0.0:
            el = Element(name=name, type="octupole")
            el.set("l", (length_mm, "mm"))
            el.params["__oct_T_per_m3"] = octo
        else:
            # General multipole (possibly combined-function).
            el = Element(name=name, type="multipole")
            el.set("l", (length_mm, "mm"))
            el.params["__mp_quad"] = quad
            el.params["__mp_sext"] = sext
            el.params["__mp_oct"] = octo
            el.params["__mp_deca"] = deca
        return el

    def _conv_solenoid(self, name, definition, place, length_mm) -> Element:
        coil_name = definition.get("coilName") or definition.get("coil")
        coil = self._coils.get(coil_name) if coil_name else None
        current = _to_float(self._tuned(definition, place, "current"))

        # Experimental: full 3D field map of the ported coil field.
        if self.solenoid_field_map and coil is not None and length_mm > 0:
            el = self._solenoid_field_map(name, coil, current, length_mm)
            if el is not None:
                return el

        # Default: native BDSIM solenoid, with ks derived from the peak on-axis
        # field of the ported (G4beamline-matched) coil field.  ks needs the
        # rigidity, so stash the field and resolve ks in _apply_rigidity.
        el = Element(name=name, type="solenoid")
        el.set("l", (length_mm, "mm"))
        b_peak = self._solenoid_peak_field(coil, current) if coil else None
        if b_peak:
            el.params["__solenoid_bpeak"] = b_peak
        else:
            self.model.warn(
                f"solenoid '{name}': could not model field (missing coil); "
                f"set 'ks' manually."
            )
        return el

    @staticmethod
    def _solenoid_peak_field(coil, current) -> float:
        """Peak on-axis field [T] of the coil (validated against G4beamline)."""
        return _solenoid.central_field(
            _to_float(coil.get("innerRadius")),
            _to_float(coil.get("outerRadius")),
            _to_float(coil.get("length")),
            current)

    def _solenoid_margin(self, coil: G4BLCommand) -> float:
        """Fringe margin [mm] each side so the field map decays to ~0 at its
        edges (an abrupt truncation would inject a spurious radial kick).

        Returns the distance beyond the coil end where the on-axis field falls
        below 0.1% of the peak (the truncated tail integral then contributes a
        negligible Larmor-phase error), capped to keep the element sane.
        """
        a1 = _to_float(coil.get("innerRadius"))
        a2 = _to_float(coil.get("outerRadius"))
        clen = _to_float(coil.get("length"))
        current = _to_float(coil.get("current")) or 1.0
        if a1 <= 0 or a2 <= a1 or clen <= 0:
            return max(a2, 1.0)
        field_rz = _solenoid.make_coil_field(a1, a2, clen, current)
        b0 = abs(field_rz(0.0, 0.0)[1]) or 1.0
        end = clen / 2.0
        cap = max(6.0 * a2, 2.0 * clen)     # sane upper bound
        # Bisect for the distance d where |Bz(0, end+d)| = 0.001*b0.
        lo, hi = 0.0, cap
        for _ in range(28):
            mid = 0.5 * (lo + hi)
            if abs(field_rz(0.0, end + mid)[1]) > 0.001 * b0:
                lo = mid
            else:
                hi = mid
        return min(hi, cap)

    def _solenoid_field_map(self, name, coil, current, length_mm):
        """Build a BDSIM field-map element from the coil's computed field."""
        a1 = _to_float(coil.get("innerRadius"))
        a2 = _to_float(coil.get("outerRadius"))
        clen = _to_float(coil.get("length"))
        if a1 <= 0 or a2 <= a1 or clen <= 0:
            return None
        field_rz = _solenoid.make_coil_field(a1, a2, clen, current)
        # Transverse extent = the bore; particles stay within it.
        half = max(a1 * 0.98, 1.0)
        nx = ny = 33
        # length_mm already includes the fringe margins (see _element_length_mm);
        # resolve z finely enough over the whole region.
        nz = int(max(101, min(261, length_mm / max(a1, 1.0) * 12 + 101)))
        field_fn, xs, ys, zs = _solenoid.sampled_3d_map(
            field_rz, half, half, -length_mm / 2.0, length_mm / 2.0, nx, ny, nz)
        fmap = _fieldmap.build_3d(xs, ys, zs, field_fn)
        fname = f"{name}.dat"
        self.model.aux_files[fname] = fmap
        fobj = f"{name}_field"
        # Cubic interpolation is smoother than linear for the axisymmetric
        # coil field, reducing divergence errors in the transverse focusing.
        self.model.field_objects.append(
            _fieldmap.gmad_field_object(fobj, fname, 3, "cubic"))
        el = Element(name=name, type="drift")
        el.set("l", (length_mm, "mm"))
        el.set("fieldAll", fobj)
        bc = field_rz(0.0, 0.0)[1]
        el.comment = f"g4bl coil field (map {fname}), B0~{bc:.3g} T"
        self.model.warn(
            f"solenoid '{name}': full G4beamline coil field (incl. end fringe) "
            f"exported as BDSIM field map {fname}, B0~{bc:.3g} T. Tracking "
            f"validated against G4beamline to ~1% Larmor phase; use "
            f"--no-solenoid-field-map for a hard-edge native solenoid instead."
        )
        return el

    def _estimate_solenoid_B(
        self, coil: Optional[G4BLCommand], current_A_per_mm2: float
    ) -> Optional[float]:
        """Central on-axis field of a thick finite solenoid (Tesla).

        Uses B0 = mu0 * J * a1 * beta * ln[(alpha + sqrt(alpha^2+beta^2)) /
        (1 + sqrt(1+beta^2))], with alpha=a2/a1, beta=L/(2 a1).
        """
        if coil is None:
            return None
        a1 = _to_float(coil.get("innerRadius"))   # mm
        a2 = _to_float(coil.get("outerRadius"))   # mm
        length = _to_float(coil.get("length"))    # mm
        if a1 <= 0 or a2 <= a1 or length <= 0:
            return None
        mu0 = 4.0e-7 * math.pi                     # T*m/A
        J = current_A_per_mm2 * 1.0e6              # A/m^2
        a1_m = a1 * 1.0e-3
        alpha = a2 / a1
        beta = (length / 2.0) / a1
        f = beta * math.log(
            (alpha + math.sqrt(alpha * alpha + beta * beta))
            / (1.0 + math.sqrt(1.0 + beta * beta))
        )
        return mu0 * J * a1_m * f

    def _conv_pillbox(self, name, definition, place, length_mm) -> Element:
        return self._rf_element(name, definition, place, length_mm)

    def _conv_rfdevice(self, name, definition, place, length_mm) -> Element:
        return self._rf_element(name, definition, place, length_mm)

    def _rf_element(self, name, definition, place, length_mm) -> Element:
        el = Element(name=name, type="rfcavity")
        el.set("l", (length_mm, "mm"))
        max_grad = _to_float(self._tuned(definition, place, "maxGradient"))  # MV/m
        freq_ghz = _to_float(definition.get("frequency"))                    # GHz
        phase_deg = _to_float(definition.get("phaseAcc"))                    # deg
        length_m = length_mm * 1.0e-3
        if max_grad != 0.0 and length_m > 0.0:
            # Peak voltage [MV] = gradient [MV/m] * length [m].
            el.set("E", (max_grad * length_m, "MV"))
        if freq_ghz != 0.0:
            el.set("frequency", (freq_ghz, "GHz"))
        if phase_deg != 0.0:
            el.set("phase", math.radians(phase_deg))
        self.model.warn(
            f"rf cavity '{name}': phase convention differs between G4beamline "
            f"(phaseAcc) and BDSIM; verify the phase."
        )
        return el

    def _conv_tubs(self, name, definition, place, length_mm) -> Element:
        return self._passive_geometry(name, definition, length_mm, "tubs")

    def _conv_cylinder(self, name, definition, place, length_mm) -> Element:
        return self._passive_geometry(name, definition, length_mm, "tubs")

    def _conv_box(self, name, definition, place, length_mm) -> Element:
        return self._passive_geometry(name, definition, length_mm, "box")

    def _passive_geometry(self, name, definition, length_mm, shape) -> Element:
        material = definition.get("material")
        # Vacuum (or GDML disabled): a drift preserves the beamline length.
        if not self.emit_gdml or _gdml.is_vacuum(material) or length_mm <= 0:
            if not _gdml.is_vacuum(material):
                self.model.warn(
                    f"passive volume '{name}' ({definition.name}, "
                    f"material={material}) converted to a drift; enable GDML "
                    f"export to make the material interact with the beam."
                )
            return self._drift(name, length_mm)

        # Material volume -> external GDML geometry so it interacts in BDSIM.
        g4mat, mats_xml, known = self._material_for_gdml(material)
        if not known:
            self.model.warn(
                f"material '{material}' on '{name}' mapped to '{g4mat}'; verify "
                f"it is a valid Geant4/NIST material name."
            )
        # length_mm is the padded ELEMENT length; the solid keeps its true size.
        g = definition.get
        solid_len = max(_to_float(g("length")), 1.0)
        if shape == "box":
            width = _to_float(g("width"))
            height = _to_float(g("height"))
            gdml_text = _gdml.box_gdml(width, height, solid_len, g4mat,
                                       materials_xml=mats_xml)
            extent = max(width, height)
        else:
            inner = _to_float(g("innerRadius"))
            outer = _to_float(g("outerRadius") or g("radius"))
            phi0 = _to_float(g("initialPhi"))
            phi1 = g("finalPhi")
            dphi = (_to_float(phi1) - phi0) if phi1 is not None else 360.0
            gdml_text = _gdml.tubs_gdml(inner, outer, solid_len, g4mat,
                                        phi0, dphi if dphi else 360.0,
                                        materials_xml=mats_xml)
            extent = 2.0 * outer
        return self._gdml_element(name, gdml_text, length_mm, extent,
                                  definition.name, g4mat, True)

    def _conv_fieldexpr(self, name, definition, place, length_mm) -> Element:
        """Sample a G4beamline analytic field expression onto a BDSIM map.

        Box region uses coords {x,y,z}; cylinder region uses {r,z} (Br/Bphi/Bz).
        Field values are Tesla; coordinates mm.  E-field expressions are not
        exported (BDSIM would need an ebmap) -- a warning is emitted.
        """
        from .expr import compile_expr, ExprError

        g = definition.get
        factor_b = _to_float(g("factorB"), 1.0)
        if any(g(k) for k in ("Ex", "Ey", "Ez", "Er")):
            self.model.warn(
                f"fieldexpr '{name}': electric-field expressions are not "
                f"exported (only magnetic). Add a BDSIM ebmap by hand if needed."
            )
        cylinder = g("radius") is not None
        # Compile the provided expressions (missing -> zero).
        def comp(key):
            e = g(key)
            if e is None:
                return None
            try:
                return compile_expr(e)
            except ExprError:
                self.model.warn(f"fieldexpr '{name}': cannot parse {key}={e}")
                return None

        if cylinder:
            br_f, bphi_f, bz_f = comp("Br"), comp("Bphi"), comp("Bz")
            radius = _to_float(g("radius"))
            half_x = half_y = radius

            def field_fn(x, y, z):
                import math as _m
                r = _m.hypot(x, y)
                br = factor_b * br_f(r=r, z=z, x=x, y=y) if br_f else 0.0
                bphi = factor_b * bphi_f(r=r, z=z, x=x, y=y) if bphi_f else 0.0
                bz = factor_b * bz_f(r=r, z=z, x=x, y=y) if bz_f else 0.0
                if r > 1e-9:
                    return (br * x / r - bphi * y / r,
                            br * y / r + bphi * x / r, bz)
                return 0.0, 0.0, bz
        else:
            bx_f, by_f, bz_f = comp("Bx"), comp("By"), comp("Bz")
            half_x = _to_float(g("width")) / 2.0 or 1.0
            half_y = _to_float(g("height")) / 2.0 or 1.0

            def field_fn(x, y, z):
                bx = factor_b * bx_f(x=x, y=y, z=z) if bx_f else 0.0
                by = factor_b * by_f(x=x, y=y, z=z) if by_f else 0.0
                bz = factor_b * bz_f(x=x, y=y, z=z) if bz_f else 0.0
                return bx, by, bz

        nx = ny = 21
        nz = int(max(41, min(121, length_mm / 20.0 + 41)))
        xs = _fieldmap.linspace(-half_x, half_x, nx)
        ys = _fieldmap.linspace(-half_y, half_y, ny)
        zs = _fieldmap.linspace(-length_mm / 2.0, length_mm / 2.0, nz)
        fmap = _fieldmap.build_3d(xs, ys, zs, field_fn)
        fname = f"{name}.dat"
        self.model.aux_files[fname] = fmap
        fobj = f"{name}_field"
        self.model.field_objects.append(
            _fieldmap.gmad_field_object(fobj, fname, 3, "linear"))
        el = Element(name=name, type="drift")
        el.set("l", (length_mm, "mm"))
        el.set("fieldAll", fobj)
        el.comment = f"fieldexpr sampled to map {fname}"
        self.model.warn(
            f"fieldexpr '{name}': analytic field sampled onto a BDSIM field map "
            f"({fname}, {nx}x{ny}x{nz})."
        )
        return el

    def _conv_sphere(self, name, definition, place, length_mm) -> Element:
        material = definition.get("material")
        if not self.emit_gdml or _gdml.is_vacuum(material):
            return self._drift(name, length_mm)
        g4mat, mats_xml, known = self._material_for_gdml(material)
        g = definition.get
        gdml_text = _gdml.sphere_gdml(
            _to_float(g("innerRadius")), _to_float(g("outerRadius")), g4mat,
            _to_float(g("initialPhi")),
            (_to_float(g("finalPhi")) - _to_float(g("initialPhi"))) if g("finalPhi") else 360.0,
            _to_float(g("initialTheta")),
            (_to_float(g("finalTheta")) - _to_float(g("initialTheta"))) if g("finalTheta") else 180.0,
            materials_xml=mats_xml)
        extent = 2.0 * _to_float(g("outerRadius"))
        return self._gdml_element(name, gdml_text, length_mm, extent, "sphere", g4mat, known)

    def _conv_polycone(self, name, definition, place, length_mm) -> Element:
        material = definition.get("material")
        g = definition.get
        zs = _parse_float_list(g("z"))
        routs = _parse_float_list(g("outerRadius"))
        rins = _parse_float_list(g("innerRadius")) or [0.0] * len(zs)
        if not self.emit_gdml or _gdml.is_vacuum(material) or len(zs) < 2:
            return self._drift(name, length_mm)
        if len(rins) < len(zs):
            rins = rins + [0.0] * (len(zs) - len(rins))
        g4mat, mats_xml, known = self._material_for_gdml(material)
        gdml_text = _gdml.polycone_gdml(
            zs, rins, routs, g4mat, _to_float(g("initialPhi")),
            (_to_float(g("finalPhi")) - _to_float(g("initialPhi"))) if g("finalPhi") else 360.0,
            materials_xml=mats_xml)
        extent = 2.0 * (max(routs) if routs else 1.0)
        return self._gdml_element(name, gdml_text, length_mm, extent, "polycone", g4mat, known)

    def _gdml_element(self, name, gdml_text, length_mm, extent, kind, g4mat, known) -> Element:
        """Create a BDSIM element backed by generated GDML geometry."""
        if not known:
            self.model.warn(
                f"material on '{name}' mapped to '{g4mat}'; verify it is a "
                f"valid Geant4/NIST material name.")
        fname = f"{name}.gdml"
        self.model.aux_files[fname] = gdml_text
        el = Element(name=name, type="element")
        el.set("geometryFile", f"gdml:{fname}")
        el.set("l", (max(length_mm, 1e-3), "mm"))
        # 1 mm clearance beyond the GDML world box each side.
        el.set("horizontalWidth", (round(max(extent, 1.0) * 1.2 + 22.0, 3), "mm"))
        el.comment = f"{kind} target, material {g4mat}"
        self.model.warn(
            f"passive volume '{name}' ({kind}, material={g4mat}) exported as "
            f"GDML geometry ({fname}); placed on-axis.")
        return el

    def _parse_blfieldmap(self, definition: G4BLCommand):
        """Parse (and cache) the BLFieldMap file named by a fieldmap command."""
        fname = definition.get("filename") or definition.get("file")
        if not fname:
            return None
        if not hasattr(self, "_blmaps"):
            self._blmaps = {}
        if fname in self._blmaps:
            return self._blmaps[fname]
        path = fname if os.path.isabs(fname) else os.path.join(self.base_dir, fname)
        grid = None
        try:
            if os.path.exists(path):
                grid = _blfieldmap.parse(path)
        except Exception:  # noqa: BLE001
            grid = None
        self._blmaps[fname] = grid
        return grid

    def _conv_fieldmap(self, name, definition, place, length_mm) -> Element:
        grid = self._parse_blfieldmap(definition)
        fname_in = definition.get("filename") or definition.get("file") or "?"
        if grid is None:
            self.model.warn(
                f"fieldmap '{name}': BLFieldMap file '{fname_in}' not found or "
                f"not a supported (grid) map; represented as a drift.")
            return self._drift(name, length_mm)
        # 'current'/'gradient' on the fieldmap scale B/E (Tunable).
        scale = _to_float(self._tuned(definition, place, "current"), 1.0) or 1.0
        base_fn = grid.field_fn()
        field_fn = (lambda x, y, z: tuple(scale * c for c in base_fn(x, y, z))) \
            if scale != 1.0 else base_fn
        fmap = _fieldmap.build_3d(grid.xs, grid.ys, grid.zs, field_fn)
        fname = f"{name}.dat"
        self.model.aux_files[fname] = fmap
        fobj = f"{name}_field"
        self.model.field_objects.append(
            _fieldmap.gmad_field_object(fobj, fname, 3, "linear"))
        el = Element(name=name, type="drift")
        el.set("l", (length_mm, "mm"))
        el.set("fieldAll", fobj)
        el.comment = f"BLFieldMap {fname_in} -> BDSIM map {fname}"
        note = (" (E-field dropped)" if grid.has_efield else "")
        self.model.warn(
            f"fieldmap '{name}': BLFieldMap '{fname_in}' converted to BDSIM "
            f"field map ({fname}, {len(grid.xs)}x{len(grid.ys)}x{len(grid.zs)})"
            f"{note}.")
        return el

    def _conv_virtualdetector(self, name, definition, place, length_mm) -> Element:
        el = Element(name=name, type="marker")
        el.comment = "virtualdetector -> sampler"
        return el

    def _conv_detector(self, name, definition, place, length_mm) -> Element:
        el = Element(name=name, type="marker")
        el.comment = "detector -> sampler"
        return el

    # -- tunable parameter resolution --------------------------------------
    @staticmethod
    def _tuned(definition: G4BLCommand, place: G4BLCommand, key: str) -> Optional[str]:
        """A tunable param may be set on the definition or overridden on place."""
        if key in place.params:
            return place.params[key]
        return definition.params.get(key)

    # -- rigidity application ----------------------------------------------
    def _apply_rigidity(self) -> None:
        brho = self._brho()
        for pl in list(self._placements):
            el = pl.element
            self._resolve_strengths(el, brho)
            self._inject_feeddown_kickers(pl)

    def _inject_feeddown_kickers(self, pl: _Placement) -> None:
        """Replace a quadrupole's transverse offset with feed-down steering.

        An offset quad is an on-axis quad plus a dipole component
        ``B = -g*offset``; represented as two thin kickers, half the kick at
        the entrance and half at the exit (symmetric thin-lens split).
        Sign convention validated against G4beamline: for gradient g>0 and
        offset dx>0 a positive particle receives hkick=+k1*L*dx.
        """
        el = pl.element
        dx_mm = el.params.pop("__offset_dx_mm", 0.0)
        dy_mm = el.params.pop("__offset_dy_mm", 0.0)
        if dx_mm == 0.0 and dy_mm == 0.0:
            return
        k1 = el.params.get("k1")
        if k1 is None:
            self.model.warn(
                f"'{el.name}': offset quadrupole feed-down needs k1 (beam "
                f"momentum); offset ignored.")
            return
        length_m = self._length_m(el)
        hkick = k1 * length_m * (dx_mm * 1e-3)
        vkick = -k1 * length_m * (dy_mm * 1e-3)
        idx = self._placements.index(pl)
        halves = []
        for tag, z in (("in", pl.z_entry_mm), ("out", pl.z_entry_mm + pl.length_mm)):
            kick = Element(name=f"{el.name}_okick_{tag}", type="kicker")
            kick.set("l", 0)
            if hkick:
                kick.set("hkick", 0.5 * hkick)
            if vkick:
                kick.set("vkick", 0.5 * vkick)
            kick.comment = f"feed-down of {el.name} offset ({dx_mm:g},{dy_mm:g}) mm"
            halves.append(_Placement(name=kick.name, element=kick,
                                     length_mm=0.0, z_entry_mm=z))
        # Entrance kicker just before the quad, exit kicker just after.
        self._placements.insert(idx, halves[0])
        self._placements.insert(idx + 2, halves[1])
        self.model.warn(
            f"'{el.name}': transverse offset ({dx_mm:g},{dy_mm:g}) mm converted "
            f"to feed-down steering kicks (hkick={hkick:.4g}, "
            f"vkick={vkick:.4g} split entrance/exit)."
        )

    # Above this implied bend angle [rad] an rbend's pole faces (each at
    # angle/2) would overlap for typical lengths; use an sbend instead.
    _RBEND_MAX_ANGLE = 0.3

    def _resolve_strengths(self, el: Element, brho: Optional[float]) -> None:
        params = el.params
        if params.pop("__bend_check", None):
            # genericbend: switch a strong rbend to an sbend so BDSIM can
            # build the geometry (rbend faces are at angle/2 to the chord).
            fw = params.pop("__bend_fw", 0.0)
            fh = params.pop("__bend_fh", 0.0)
            b_val = params.get("B")
            b_tesla = b_val[0] if isinstance(b_val, tuple) else _to_float(b_val)
            length_m = self._length_m(el)
            if brho and length_m > 0 and b_tesla:
                angle = abs(b_tesla) * length_m / brho
                if angle > self._RBEND_MAX_ANGLE:
                    el.type = "sbend"
                    # A tightly-curved sbend cannot be wider than a fraction of
                    # its bend radius (BDSIM's face-overlap check trips at
                    # ~0.65*rho empirically), so fit the magnet body and
                    # aperture to the bend radius and g4bl field aperture.
                    rho_mm = brho / abs(b_tesla) * 1000.0
                    hw = 0.5 * rho_mm
                    if fw > 0:
                        hw = min(hw, max(fw, fh) * 1.2 + 40.0)
                    aper1 = min(fw / 2.0 if fw > 0 else hw / 2.0 - 10.0,
                                hw / 2.0 - 10.0)
                    aper2 = min(fh / 2.0 if fh > 0 else aper1,
                                hw / 2.0 - 10.0)
                    el.set("apertureType", "rectangular")
                    el.set("aper1", (round(aper1, 3), "mm"))
                    el.set("aper2", (round(aper2, 3), "mm"))
                    el.set("horizontalWidth", (round(hw, 3), "mm"))
                    self.model.warn(
                        f"bend '{el.name}': implied angle {angle:.3g} rad is too "
                        f"large for BDSIM rbend pole-face geometry; emitted as "
                        f"sbend (arc) with aperture/width fitted to the bend "
                        f"radius (may clip more than the g4bl aperture)."
                    )
        if "__solenoid_bpeak" in params:
            b_peak = params.pop("__solenoid_bpeak")
            if brho:
                ks = b_peak / brho          # ks = B / (B*rho)  [m^-1]
                el.set("ks", ks)
                el.comment = (f"native solenoid, ks from G4beamline coil field "
                              f"B0={b_peak:.4g} T (hard-edge approx of the "
                              f"fringe field)")
            else:
                el.comment = f"solenoid B0={b_peak:.4g} T (set ks once Brho known)"
        if "__gradient_T_per_m" in params:
            grad = params.pop("__gradient_T_per_m")
            if brho:
                el.set("k1", grad / brho)
            else:
                el.comment = f"gradient={grad} T/m (set k1 once Brho known)"
        if "__sext_T_per_m2" in params:
            s = params.pop("__sext_T_per_m2")
            if brho:
                el.set("k2", 2.0 * s / brho)
            else:
                el.comment = f"sextupole={s} T/m^2"
        if "__oct_T_per_m3" in params:
            o = params.pop("__oct_T_per_m3")
            if brho:
                el.set("k3", 6.0 * o / brho)
            else:
                el.comment = f"octopole={o} T/m^3"
        if "__mp_quad" in params:
            quad = params.pop("__mp_quad")
            sext = params.pop("__mp_sext", 0.0)
            octo = params.pop("__mp_oct", 0.0)
            deca = params.pop("__mp_deca", 0.0)
            length_m = self._length_m(el)
            self._set_multipole_knl(el, brho, length_m, quad, sext, octo, deca)

    @staticmethod
    def _length_m(el: Element) -> float:
        l = el.params.get("l")
        if isinstance(l, tuple):
            value, unit = l
            scale = {"mm": 1e-3, "cm": 1e-2, "m": 1.0}.get(unit, 1.0)
            return value * scale
        return _to_float(str(l)) if l is not None else 0.0

    def _set_multipole_knl(self, el, brho, length_m, quad, sext, octo, deca) -> None:
        # BDSIM 'multipole' uses integrated normalised strengths knl = kn * L.
        if not brho or length_m <= 0:
            el.comment = (
                f"combined-function multipole: quad={quad} T/m, sext={sext} "
                f"T/m^2, oct={octo} T/m^3 (needs Brho & length to normalise)"
            )
            return
        k1 = quad / brho
        k2 = 2.0 * sext / brho
        k3 = 6.0 * octo / brho
        k4 = 24.0 * deca / brho
        knl = [k1 * length_m, k2 * length_m, k3 * length_m, k4 * length_m]
        # Trim trailing zeros.
        while knl and knl[-1] == 0.0:
            knl.pop()
        if knl:
            el.set("knl", "{" + ", ".join(f"{v:.10g}" for v in knl) + "}")

    # -- beamline assembly --------------------------------------------------
    def _build_beamline(self) -> None:
        placements = sorted(self._placements, key=lambda p: p.z_entry_mm)
        line: List[str] = []
        drift_index = 0
        prev_exit: Optional[float] = None

        # Leading drift so the BDSIM beamline starts at z=0, matching the
        # G4beamline centreline origin.
        if placements and placements[0].z_entry_mm > _Z_TOL:
            drift_index += 1
            dname = self._unique_name(f"drift_{drift_index}")
            drift = self._drift(dname, placements[0].z_entry_mm)
            self.model.add_element(drift)
            line.append(dname)
            prev_exit = placements[0].z_entry_mm

        for pl in placements:
            if prev_exit is not None:
                gap = pl.z_entry_mm - prev_exit
                if gap > _Z_TOL:
                    drift_index += 1
                    dname = self._unique_name(f"drift_{drift_index}")
                    drift = self._drift(dname, gap)
                    self.model.add_element(drift)
                    line.append(dname)
                elif gap < -_Z_TOL:
                    self.model.warn(
                        f"element '{pl.name}' overlaps the previous element by "
                        f"{-gap:.4g} mm; no drift inserted (check geometry)."
                    )

            self.model.add_element(pl.element)
            line.append(pl.name)
            if pl.is_sampler:
                self.model.samplers.append(pl.name)

            exit_z = pl.z_entry_mm + pl.length_mm
            prev_exit = exit_z if prev_exit is None else max(prev_exit, exit_z)

        self.model.line = line

    # -- beam / options blocks ---------------------------------------------
    def _build_beam_block(self) -> None:
        beam = self.model.beam
        beam["particle"] = _PARTICLE_NAME.get(self._particle, self._particle)
        if self._momentum_mev is not None:
            beam["momentum"] = (self._momentum_mev, "MeV")

        cmd = self._beam_cmd
        if cmd is None:
            return

        # Beam centroid offsets: G4beamline beamX/beamY [mm] and meanXp/meanYp
        # [slope] -> BDSIM X0/Y0 [m] and Xp0/Yp0.
        bx = cmd.get("beamX") or cmd.get("x")
        by = cmd.get("beamY") or cmd.get("y")
        mxp = cmd.get("meanXp") or cmd.get("beamXp")
        myp = cmd.get("meanYp") or cmd.get("beamYp")
        if bx is not None and _to_float(bx) != 0.0:
            beam["X0"] = (_to_float(bx), "mm")
        if by is not None and _to_float(by) != 0.0:
            beam["Y0"] = (_to_float(by), "mm")
        if mxp is not None and _to_float(mxp) != 0.0:
            beam["Xp0"] = _to_float(mxp)
        if myp is not None and _to_float(myp) != 0.0:
            beam["Yp0"] = _to_float(myp)

        btype = cmd.args[0] if cmd.args else "gaussian"
        if btype == "gaussian":
            beam["distrType"] = "gauss"
            self._copy_sigma(cmd, beam)
        elif btype == "rectangular":
            beam["distrType"] = "square"
            w = cmd.get("beamWidth")
            h = cmd.get("beamHeight")
            if w is not None:
                beam["envelopeX"] = (_to_float(w) / 2.0, "mm")
            if h is not None:
                beam["envelopeY"] = (_to_float(h) / 2.0, "mm")
        elif btype in ("ascii", "root"):
            beam["distrType"] = "userfile"
            self.model.warn(
                "file-based beam ('%s') needs a BDSIM userfile column spec; "
                "set distrFile/columns manually." % btype
            )
        else:
            beam["distrType"] = "reference"

    def _copy_sigma(self, cmd: G4BLCommand, beam) -> None:
        # In G4beamline a NEGATIVE sigma denotes a flat/uniform distribution
        # with |sigma| as the half-width.  BDSIM's gauss needs a positive RMS,
        # so we use |sigma| (approximating flat as gaussian) and warn once.
        sx = cmd.get("sigmaX")
        sy = cmd.get("sigmaY")
        sxp = cmd.get("sigmaXp")
        syp = cmd.get("sigmaYp")
        st = cmd.get("sigmaT")
        sp = cmd.get("sigmaP")
        flat = any(s is not None and _to_float(s) < 0
                   for s in (sx, sy, sxp, syp, st, sp))
        if flat:
            self.model.warn(
                "beam has negative sigma(s) (G4beamline flat distribution); "
                "converted to a gaussian with |sigma| as RMS."
            )
        if sx is not None:
            beam["sigmaX"] = (abs(_to_float(sx)), "mm")
        if sy is not None:
            beam["sigmaY"] = (abs(_to_float(sy)), "mm")
        if sxp is not None:
            beam["sigmaXp"] = abs(_to_float(sxp))
        if syp is not None:
            beam["sigmaYp"] = abs(_to_float(syp))
        if st is not None:
            beam["sigmaT"] = (abs(_to_float(st)), "ns")
        if sp is not None and self._momentum_mev:
            # Relative momentum spread ~ relative energy spread (ultra-rel.).
            beam["sigmaE"] = abs(_to_float(sp)) / self._momentum_mev

    def _build_options_block(self) -> None:
        opt = self.model.options
        if self._physics:
            opt["physicsList"] = self._map_physics(self._physics)
        if self._n_events is not None:
            opt["ngenerate"] = self._n_events
        # BDSIM always builds a beam pipe; G4beamline does not.  Size the pipe to
        # the largest aperture/detector so a wide (vacuum) beam is not clipped or
        # scattered by the pipe wall, matching G4beamline's free-space tracking.
        if self._max_aperture_mm > 0:
            r = round(self._max_aperture_mm, 3)
            opt["beampipeRadius"] = (r, "mm")
            # BDSIM requires the magnet outer width (horizontalWidth) to exceed
            # 2*(aper1 + beampipe thickness); enlarge it to match the big pipe.
            opt["horizontalWidth"] = (round(2 * r + 200.0, 3), "mm")

    def _map_physics(self, name: str) -> str:
        if name.lower() == "default":
            return "g4FTFP_BERT"
        if name.startswith("g4"):
            return name
        # Strip a trailing Geant4 EM-option suffix (_EMV/_EMX/.../_SS) to check
        # the base reference-list name.
        base = re.sub(r"_(EMV|EMX|EMY|EMZ|LIV|PEN|GS|SS|WVI|LE)$", "", name)
        if base in _VALID_G4_LISTS:
            # BDSIM wraps Geant4 reference lists with a 'g4' prefix.
            return "g4" + name
        self.model.warn(
            f"physics list '{name}' is not a valid Geant4 reference list in "
            f"BDSIM; using g4FTFP_BERT. Set option,physicsList=... manually if "
            f"you need a specific list."
        )
        return "g4FTFP_BERT"


# Valid Geant4 reference physics-list base names accepted by BDSIM (with the
# 'g4' prefix).  Bare 'QGSP' and similar removed-from-Geant4 names are excluded.
_VALID_G4_LISTS = {
    "FTFP_BERT", "FTFP_BERT_HP", "FTFP_BERT_TRV", "FTFP_INCLXX", "FTF_BIC",
    "LBE", "QBBC", "QGSP_BERT", "QGSP_BERT_HP", "QGSP_BIC", "QGSP_BIC_HP",
    "QGSP_BIC_AllHP", "QGSP_FTFP_BERT", "QGSP_INCLXX", "QGSP_INCLXX_HP",
    "QGS_BIC", "Shielding", "ShieldingLEND", "NuBeam", "FTFQGSP_BERT",
}


# G4beamline commands with no BDSIM/GMAD equivalent -- warn when encountered so
# the conversion is transparent (see validation/LIMITATIONS.md for workarounds).
_UNSUPPORTED_COMMANDS = {
    "spacecharge": "BDSIM is a single-particle tracker; space charge / "
                   "collective effects are not modelled.",
    "helicaldipole": "no standard BDSIM helical-dipole element; approximate "
                     "with a field map or a sequence of dipoles.",
    "tune": "BDSIM does not auto-tune; set magnet strengths explicitly or "
            "match externally with pybdsim.",
    "fieldlines": "field-line visualisation is a G4beamline-only feature.",
    "printfield": "field printing is a G4beamline-only feature.",
    "particlefilter": "per-region particle filtering has no direct equivalent; "
                      "use collimators or minimumKineticEnergy.",
}


# Element-definition command names recognised by the converter.
_DEFINITION_TYPES = {
    "genericbend", "genericquad", "idealsectorbend", "genericsectorbend",
    "multipole", "solenoid", "pillbox", "rfdevice", "tubs", "cylinder",
    "box", "virtualdetector", "detector", "fieldexpr", "sphere", "polycone",
    "fieldmap",
}


def convert_commands(commands: List[G4BLCommand], source_name: str = "") -> BdsimModel:
    return Converter(commands, source_name=source_name).convert()
