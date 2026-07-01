# Validation harness

This directory cross-checks the converter by running **both** simulation codes in
Docker and comparing the beam at equivalent planes:

* **G4beamline** on the original `.g4bl` input
  (image built from [`JeffersonLab/docker-g4beamline`](https://github.com/JeffersonLab/docker-g4beamline),
  G4beamline 3.06 / Geant4 10.5).
* **BDSIM** on the converted `.gmad`
  (official image `bdsim/ubuntu24-g4.11.3-bdsim`, Geant4 11.3).

## What is compared

For optics validation we track a **reference particle** with a small transverse
offset and compare its position `x, y` and angle `xp = px/pz, yp = py/pz` at a
detector / sampler plane.  A correct `k1` / bend conversion reproduces the same
transfer through the element.

To make the comparison meaningful the optics test inputs:

* run in vacuum with iron and fringe fields disabled (so both codes see an ideal
  hard-edged field), and
* use a single reference particle (no stochastic processes).

## Scripts

| Script | Purpose |
|---|---|
| `build_images.sh` | Build the G4beamline image and pull the BDSIM image. |
| `g4bl_ascii.py`   | Parse a G4beamline ASCII `virtualdetector` NTuple. |
| `bdsim_dump.py`   | (runs *inside* the BDSIM container) dump a sampler to ASCII via ROOT. |
| `compare.py`      | End-to-end: run both codes for a case and compare. |

## Usage

```bash
./validation/build_images.sh                 # one-off, slow
python3 validation/compare.py validation/cases/optics_quad.g4bl
```
