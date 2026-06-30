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
| `solenoid` (+`coil`) | `solenoid` | central `B` **estimated** from current density (thick-solenoid formula) — verify |
| `pillbox`, `rfdevice` | `rfcavity` | `E = maxGradient·L`, `frequency` [GHz], `phase` — phase convention differs |
| `tubs`, `cylinder`, `box` | `drift` | passive geometry; length kept, material **not** transferred |
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

* **Solenoids**: BDSIM wants `ks`/`B`; G4beamline defines a coil + current and
  computes the field map.  We estimate the central on-axis field analytically —
  good for a starting point, but verify against G4beamline's field.
* **RF phase**: G4beamline `phaseAcc` (0° = rising zero-crossing) and BDSIM
  `phase` do not share a zero — check the phase.
* **Bends and 3-D geometry**: BDSIM follows the reference orbit automatically, so
  `corner`/`cornerarc` are ignored (the magnet carries the bend).  Transverse
  `x`/`y` placement offsets and rotations are not mapped to the 1-D beamline.
* **Targets/collimators** (`box`, `tubs`): represented as drifts to keep the
  geometry length; replace with a BDSIM collimator or custom geometry if the
  block should interact with the beam.
* `parent=` (nested) placements are not flattened into the beamline.

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
