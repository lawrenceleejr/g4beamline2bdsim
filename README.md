# g4beamline2bdsim

Convert [G4beamline](http://g4beamline.muonsinc.com/) input files into
[BDSIM](https://bdsim-collaboration.github.io/bdsim/sphinx/index.html) (GMAD)
configurations.

G4beamline and BDSIM are both Geant4-based accelerator simulation codes, but
they describe a machine very differently:

| | G4beamline | BDSIM (GMAD) |
|---|---|---|
| Lattice | elements *placed* at absolute `z` | *sequence* of elements that tile the line |
| Lengths | millimetres | metres |
| Energy/momentum | MeV, MeV/c | GeV |
| Magnets | physical fields (T, T/m, …) | normalised strengths (`k1` [m⁻²], …) |
| Bends | `By` [T] + geometry | `angle` [rad] or `B` [T] |

This tool parses a `.g4bl` file, builds an ordered beamline (inserting drifts to
fill the gaps between placed elements), converts each magnet's physical field to
the normalised strength BDSIM expects using the beam rigidity, and writes a
valid `.gmad` model.

## Installation

```bash
pip install -e .
```

Pure Python, no runtime dependencies (Python ≥ 3.8).

## Usage

```bash
g4beamline2bdsim input.g4bl              # writes input.gmad
g4beamline2bdsim input.g4bl -o out.gmad
g4beamline2bdsim input.g4bl --sample-all # 'sample, all;' instead of per-detector
cat input.g4bl | g4beamline2bdsim -      # stdin -> stdout
```

Conversion warnings (unsupported elements, field estimates, convention
mismatches) are printed to stderr and embedded as comments at the top of the
generated GMAD.

As a library:

```python
from g4beamline2bdsim import parse_g4bl, Converter, GmadWriter

commands, _ = parse_g4bl(open("input.g4bl").read())
model = Converter(commands).convert()
print(GmadWriter(model).to_string())
```

## What is converted

| G4beamline | BDSIM | Notes |
|---|---|---|
| `genericquad` | `quadrupole` | `k1 = gradient[T/m] / Brho` |
| `genericbend` | `rbend` | uniform box field → `B` [T] |
| `idealsectorbend` | `sbend` | geometric `angle` (deg→rad); arc length `R·θ` |
| `genericsectorbend` | `sbend` (+`k1`) | dipole `angle`/`B`, combined `k1` |
| `multipole` | `quadrupole`/`sextupole`/`octupole` or combined `multipole` (`knl`) | `kₙ = n!·strengthₙ / Brho` |
| `solenoid` (+`coil`) | field map on a `drift` (default) or native `solenoid` | coil field computed by current-loop superposition (matches G4beamline `printfield` to 4 sig figs); **default** = full 3-D field map (g4bl-style); `--no-solenoid-field-map` = native `solenoid` with `ks` from the peak field |
| `pillbox`, `rfdevice` | `rfcavity` | `E = maxGradient·L`, `frequency` [GHz], `phase` — phase convention differs |
| `box`, `tubs`, `cylinder`, `sphere`, `polycone` | `element` + **GDML** (or `drift`) | material volumes exported as GDML geometry so they interact with the beam; vacuum → drift; `--no-gdml` forces drifts |
| `fieldexpr` | field map on a `drift` | analytic `Bx/By/Bz` (or `Br/Bphi/Bz`) sampled onto a BDSIM 3-D field map |
| `fieldmap` (BLFieldMap) | field map on a `drift` | BLFieldMap grid file parsed and re-written as a BDSIM field map |
| `virtualdetector`, `detector` | `marker` + `sample` | sampler plane |
| `beam`, `reference` | `beam` | particle, momentum, gaussian σ's, centroid offsets |
| `physics` | `option, physicsList=...` | reference lists prefixed `g4` |

The magnetic rigidity `Brho = p[GeV/c] / (0.299792458·|q|)` is taken from the
`reference`/`beam` momentum.  If no momentum is found, `k`-values cannot be
computed and the field is left in a comment with a warning.

### Units

Lengths, positions and apertures are kept in **millimetres** and emitted with
`*mm`, preserving the exact G4beamline numbers.  Momentum is emitted as `*MeV`,
RF frequency as `*GHz`, times as `*ns`; angles are radians.

## Limitations / things to check by hand

See [`validation/LIMITATIONS.md`](validation/LIMITATIONS.md) for the full
feature-gap matrix (what maps, what needs a workaround, what has no BDSIM
equivalent).  Highlights:

* **Solenoids**: the G4beamline coil field is ported exactly (matches
  `printfield` to 4 sig figs). **By default** the full 3-D coil field (with end
  fringe) is written as a BDSIM field map — the faithful g4bl-style field, with
  tracking validated against G4beamline (≤ 0.1 mm centroid for a 2 T solenoid,
  ~1 % Larmor phase for a 6 T one). `--no-solenoid-field-map` gives a hard-edge
  native BDSIM `solenoid` with `ks` from the peak field instead.
* **RF phase**: G4beamline `phaseAcc` (0° = rising zero-crossing) and BDSIM
  `phase` do not share a zero — check the phase.
* **Bends and 3-D geometry**: BDSIM follows the reference orbit automatically, so
  `corner`/`cornerarc` are ignored. Transverse `x`/`y` placement offsets and
  rotations are not mapped to the 1-D beamline.
* **Targets/collimators/solids** (`box`, `tubs`, `cylinder`, `sphere`,
  `polycone`): exported as GDML geometry so the material interacts (validated);
  `--no-gdml` forces plain drifts.
* **Fields**: `fieldexpr` formulas and `fieldmap` (BLFieldMap) files are
  converted to BDSIM field maps automatically.
* **Not converted** (a warning is emitted): `spacecharge`, `helicaldipole`,
  `tune`, in-language `do`/`if`/`define`, `parent=` nested placements.

## Validation

The converter is cross-checked by running **both** codes in Docker and comparing
the beam at equivalent planes — see [`validation/`](validation/) and
[`validation/RESULTS.md`](validation/RESULTS.md).  Summary:

* **Ideal optics** (drift, single quad, FD doublet, sextupole): mean trajectory
  agrees to **~10⁻⁷ rad / <0.05 mm**; the quadrupole result confirms
  `k1 = gradient/Brho` and the sextupole confirms the `k2 = 2·strength/Brho`
  factorial factor.
* **Dipole**: G4beamline's deflection matches the analytic `B·L/Brho` (the angle
  BDSIM realises) to **~0.5 %**.
* **A 21-element, 8-cell FODO channel**: beam centroid and envelope agree across
  **five detectors** spanning the whole line.
* **Full example lattices** with iron, RF and solenoids run end-to-end in both
  codes.

```bash
./validation/build_images.sh   # one-off: build G4beamline, pull BDSIM (slow)
./validation/run_all.sh        # run the whole comparison suite
```

## Development

```bash
pip install pytest
pytest                         # 26 unit tests for parser/converter/writer
```

## Repository layout

```
g4beamline2bdsim/   parser, model, converter, gmad_writer, cli
examples/           sample .g4bl inputs and their converted .gmad
tests/              pytest unit tests
validation/         Docker-based cross-check against G4beamline + BDSIM
```
