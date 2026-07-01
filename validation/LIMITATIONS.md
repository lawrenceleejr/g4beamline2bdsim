# Feature gaps: G4beamline → BDSIM

This is the map of what does and does not carry across, learned from converting
and running all 44 example lattices (see `RESULTS.md`).  It separates **genuine
BDSIM capability gaps** from **different-mechanism** cases and **converter
limitations**, and gives a workaround for each.

Legend for "status":
- 🟢 converted automatically
- 🟡 converted with a documented approximation / needs a check
- 🔵 possible in BDSIM but through a different mechanism (may need manual work)
- 🔴 no BDSIM equivalent

## 1. Genuinely absent in BDSIM (🔴)

| G4beamline feature | Why it doesn't map | Workaround |
|---|---|---|
| **Space charge / collective effects** (`spacecharge`) | BDSIM is a single-particle Geant4 tracker; particles don't see each other. | None in BDSIM. Use a space-charge code (e.g. IMPACT, GPT) for that regime. |
| **Analytic field expressions** (`fieldexpr`, e.g. `By=1.45*cos(z/L)`) | GMAD cannot express a field as a formula. | Sample the expression onto a grid and load it as a BDSIM field map; or approximate with magnet elements. |
| **Coil + current-density solenoid fields** (`coil`+`solenoid current=`) | BDSIM's solenoid is an ideal hard-edge `ks`/`B`; it has no coil field solver. | The converter estimates the central field (§ below). For the true end/fringe field, compute it in G4beamline (`printfield`) and load as a BDSIM field map. |
| **Reference-particle / strength auto-tuning** (`tune`, `tuneMomentum`, `tuneZ`) | BDSIM only forward-simulates. | Match externally with `pybdsim`/`pymadx`, or set strengths explicitly. |
| **In-language scripting** (`do` loops, `if/else`, `define` macros, `param` arithmetic beyond constants) | GMAD has variables/expressions and `include`, but no loops/conditionals/macros. | Generate the lattice programmatically with `pybdsim` (the BDSIM-recommended route). |
| **Field-line / field visualisation** (`fieldlines`, `printfield`) | Diagnostic, G4beamline-only. | Use the BDSIM/Geant4 visualiser; `printfield` has no equivalent. |
| **Helical dipole** (`helicaldipole`) | No standard BDSIM element. | Build from a field map or a sequence of rotated dipoles. |
| **Per-region particle filters** (`particlefilter`, `trackcuts keep=`) | No fine-grained keep/kill by species mid-lattice. | Use collimators, `minimumKineticEnergy`, or the kill flag on elements. |

The converter now emits an explicit **warning** for each of these when it sees
them, so nothing is silently dropped.

## 2. Present in BDSIM, different mechanism (🔵)

| G4beamline feature | BDSIM equivalent | Converter status |
|---|---|---|
| **Material targets / absorbers** (`box`, `tubs`, `cylinder` of a real material) | External **GDML** geometry placed as an `element` | 🟢 **Now converted** — see §4. |
| **Tabulated field maps** (`fieldmap`, TOSCA/BLFieldMap) | BDSIM `field` object (`bdsim1d/2d/3d/4d`, `poisson`, …) attached via `fieldAll` | 🔴 not yet — the BLFieldMap file must be reformatted; the converter warns and points here. |
| **Arbitrary 3-D geometry** (polycone, extrusion, sphere, `group`/`parent` trees) | External GDML / `placement` | 🟡 only `box`/`tubs`/`cylinder` are exported to GDML today. |
| **Flexible centreline steering** (`corner`, `cornerarc`, rings via `start`) | BDSIM derives geometry from the beamline sequence; rings via `nturns` | 🟢 the bend elements carry the deflection; `corner`/`cornerarc` are ignored (with a warning). |
| **File-driven beams** (`beam ascii/root file=`) | BDSIM `distrType="userfile"` + column spec | 🟡 flagged with a warning; the column spec must be set by hand. |

### Field-map conversion recipe (manual, for `fieldmap`)

1. Convert the BLFieldMap file to a BDSIM field-map format (BDSIM ships
   converters in `pybdsim.Field`; e.g. write a `bdsim2d`/`bdsim3d` file).
2. Define a field object and attach it to a drift spanning the region:
   ```gmad
   fm: field, type="bmap3d", magneticFile="bdsim3d:map.dat",
       magneticInterpolator="cubic";
   reg: drift, l=<region length>*mm, fieldAll=fm;
   ```
3. Insert `reg` in the line at the mapped z-range.

## 3. Modelling differences (🟡 — closure caveats)

| Area | G4beamline | BDSIM | Effect |
|---|---|---|---|
| Magnet fringe/edge | Enge fringe on by default | hard-edge by default (optional edges) | small focusing differences; the optics tests disable fringe to match. |
| RF phase | `phaseAcc` (0° = rising zero-crossing), rich timing | `rfcavity` `phase` | phase zero differs — verify after conversion. |
| Dipole geometry | box field (`genericbend`) or sector | `rbend`/`sbend` | strong short bends can exceed BDSIM's `rbend` face geometry (see MICE example). |
| Beam pipe / world | only what you place | always a pipe + yokes | converter sizes the pipe to the apertures to avoid clipping. |
| Geant4 version | 10.5 (image) | 11.3 (image) | different scattering/energy-loss models at the % level. |

## 4. What the converter does about material (GDML export)

By default, a `box`/`tubs`/`cylinder` with a real material is exported as a
**GDML solid** and placed as a BDSIM `element`, so the material actually
interacts with the beam (energy loss, scattering, secondaries) instead of being
dropped to a drift.  Verified against G4beamline for a 100 mm tungsten target
(2 GeV/c protons):

| | protons surviving | mean E | σx |
|---|---|---|---|
| G4beamline | 62 / 200 | 1.978 GeV | 17.48 mm |
| BDSIM (GDML target) | 67 / 200 | 1.980 GeV | 17.49 mm |

Energy loss and scattering now agree to < 0.1 %; before, BDSIM saw the target as
a drift and showed ~0 scattering.

Notes:
- Material names are mapped to Geant4/NIST names (`W`→`G4_W`, `Cu`→`G4_Cu`, …);
  unknown names are passed through as `G4_<name>` with a warning.
- Transverse offsets/rotations of the volume are placed on-axis.
- Vacuum volumes stay drifts. Use `--no-gdml` to force the old drift behaviour.
- Only `box`, `tubs`, `cylinder` are exported; other solids are not yet.
