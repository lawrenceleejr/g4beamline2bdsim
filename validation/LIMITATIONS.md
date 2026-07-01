# Feature gaps: G4beamline → BDSIM

This is the map of what does and does not carry across, learned from converting
and running all the example lattices (see `RESULTS.md`).  It separates **genuine
BDSIM capability gaps** from **different-mechanism** cases and **modelling
differences**, and gives a workaround for each.

Legend for "status":
- 🟢 converted automatically
- 🟡 converted with a documented approximation / needs a check
- 🔵 possible in BDSIM through a different mechanism
- 🔴 no BDSIM equivalent

## 1. Genuinely absent in BDSIM (🔴)

| G4beamline feature | Why it doesn't map | Workaround |
|---|---|---|
| **Space charge / collective effects** (`spacecharge`) | BDSIM is a single-particle Geant4 tracker; particles don't see each other. | None in BDSIM. Use a space-charge code (IMPACT, GPT). |
| **Reference-particle / strength auto-tuning** (`tune`, `tuneMomentum`, `tuneZ`) | BDSIM only forward-simulates. | Match externally with `pybdsim`/`pymadx`, or set strengths explicitly. |
| **In-language scripting** (`do` loops, `if/else`, `define` macros) | GMAD has variables/expressions and `include`, but no loops/conditionals/macros. | Generate the lattice with `pybdsim`. (Numeric `param` **arithmetic is** evaluated by the converter, incl. `sqrt`, `^`, `if`, trig.) |
| **Field-line / field visualisation** (`fieldlines`, `printfield`) | Diagnostic, G4beamline-only. | Use the BDSIM/Geant4 visualiser. |
| **Helical dipole** (`helicaldipole`) | No standard BDSIM element. | Build from a field map or rotated dipoles. |
| **Per-region particle filters** (`particlefilter`, `trackcuts keep=`) | No fine-grained keep/kill by species mid-lattice. | Use collimators, `minimumKineticEnergy`, or element kill flags. |

The converter emits an explicit **warning** for each of these, so nothing is
silently dropped.

## 2. Now handled by the converter (previously gaps)

These G4beamline features are no longer gaps — the converter translates them,
mostly by generating BDSIM **field maps** or **GDML geometry**:

| Feature | How it is converted | Validation |
|---|---|---|
| **Material targets/absorbers** (`box`, `tubs`, `cylinder`, `sphere`, `polycone`) | Exported as **GDML** geometry placed as a BDSIM `element` (material interacts) | 100 mm W target: energy loss & scattering match G4beamline to < 0.1 % |
| **Analytic field expressions** (`fieldexpr`) | Formula(s) `Bx/By/Bz` (box `{x,y,z}`) or `Br/Bphi/Bz` (cylinder `{r,z}`) **auto-sampled onto a BDSIM 3-D field map** on a drift | uniform `By`: deflection matches to 5 sig figs |
| **BLFieldMap files** (`fieldmap filename=...`) | Parsed (grid format; `normB` scaling, `extend*` mirror symmetry) and **re-written as a BDSIM field map** | round-trip uniform `By` reproduces `B·L/Bρ` exactly |
| **Coil + current-density solenoids** (`coil` + `solenoid current=`) | Field computed from the coil by **current-loop superposition** (elliptic integrals) → native BDSIM `solenoid` with `ks` from the exact peak field; `--solenoid-field-map` writes the full 3-D map | on-axis `Bz` matches G4beamline `printfield` to 4 sig figs |

### The field-map bridge

`fieldexpr`, `fieldmap`, and (optionally) `solenoid` all funnel through one
mechanism: sample/parse a field onto a grid and write a BDSIM field map
(`g4beamline2bdsim/fieldmap.py`), then attach it to a drift with
`fieldAll="…"`.  The BDSIM ASCII format was reverse-engineered and verified
empirically (positions **mm**, field **Tesla**, loop order X-outer→Z-inner).
Only **magnetic** fields are exported; `fieldexpr` electric-field terms are
warned (BDSIM would need an `ebmap`).

## 3. Different mechanism / still manual (🔵🟡)

| G4beamline feature | BDSIM equivalent | Status |
|---|---|---|
| Arbitrary 3-D geometry (`extrusion`, `trap`, `torus`, `group`/`parent` trees) | External GDML / `placement` | 🟡 `box`/`tubs`/`cylinder`/`sphere`/`polycone` are exported; other solids → drift + warning |
| Flexible centreline steering (`corner`, `cornerarc`, rings via `start`) | BDSIM derives geometry from the sequence; rings via `nturns` | 🟢 bends carry the deflection; corners ignored (warned) |
| File-driven beams (`beam ascii/root file=`) | `distrType="userfile"` + column spec | 🟡 warned; column spec set by hand |
| Cylindrical BLFieldMap (`cylinder` grid) | 3-D field map | 🟡 not parsed yet (grid maps only); warned |

## 4. Solenoid: what "porting the field" means, and its caveat

The G4beamline solenoid field **is** ported exactly — `solenoid.py` reproduces
the coil field (`printfield` on-axis `Bz` agrees to 4 sig figs, including the
end fringe).  Two delivery modes:

* **Default** — a native BDSIM `solenoid` with `ks = B_peak/Bρ`, `B_peak` from
  the ported field.  Robust; tracks correctly.  Because it is hard-edged it
  differs from the real fringe field by the usual hard-edge amount (~10–20 % on
  the transverse centroid for a strong solenoid).
* **`--solenoid-field-map`** (experimental) — the full 3-D coil field as a BDSIM
  field map.  This carries the true fringe, but tracking a **strong**
  axisymmetric field through a Cartesian map is delicate: linear interpolation
  of an axisymmetric field is not exactly divergence-free, which can perturb the
  focusing.  Use with a fine grid and validate against G4beamline.

## 5. Modelling differences (🟡 — closure caveats)

| Area | G4beamline | BDSIM | Effect |
|---|---|---|---|
| Magnet fringe/edge | Enge fringe on by default | hard-edge by default | small focusing differences |
| RF phase | `phaseAcc` (0° = rising zero-crossing) | `rfcavity` `phase` | phase zero differs — verify |
| Dipole geometry | box field or sector | `rbend`/`sbend` | strong short bends can exceed BDSIM's `rbend` face geometry |
| Beam pipe / world | only what you place | always a pipe + yokes | converter sizes the pipe to the apertures |
| Geant4 version | 10.5 (image) | 11.3 (image) | %-level scattering/energy-loss differences |
