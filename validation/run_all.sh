#!/usr/bin/env bash
# Run the full validation suite: quantitative optics comparison, the dipole
# bend-angle check, and smoke tests of the full example lattices.
#
# Prerequisites: ./validation/build_images.sh (builds/pulls the two images).
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(dirname "$HERE")"
cd "$REPO"

N="${N:-4000}"
rc=0

echo "######################################################################"
echo "# 1. Quantitative optics comparison (frame-aligned, ideal fields)"
echo "######################################################################"
python3 validation/compare.py \
    validation/cases/optics_drift.g4bl \
    validation/cases/optics_quad.g4bl \
    validation/cases/optics_doublet.g4bl \
    validation/cases/optics_sextupole.g4bl \
    validation/cases/fodo_channel.g4bl \
    -n "$N" --atol-mm 0.6 || rc=1

echo
echo "######################################################################"
echo "# 2. Dipole bend-angle check (G4beamline deflection vs analytic angle)"
echo "######################################################################"
python3 validation/check_bend.py || rc=1
python3 validation/check_bend.py --B 0.2 --L 1.0 --p 2000 || rc=1

echo
echo "######################################################################"
echo "# 3. Smoke tests of full example lattices (iron / RF / solenoid)"
echo "######################################################################"
python3 validation/compare.py \
    examples/fodo.g4bl \
    examples/bend_line.g4bl \
    examples/rf_solenoid.g4bl \
    -n 500 --smoke || rc=1

echo
if [ "$rc" -eq 0 ]; then
    echo "ALL VALIDATION PASSED"
else
    echo "SOME VALIDATION FAILED (rc=$rc)"
fi
exit "$rc"
