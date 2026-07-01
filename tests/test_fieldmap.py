import math

from g4beamline2bdsim import blfieldmap, fieldmap, solenoid
from g4beamline2bdsim.converter import Converter
from g4beamline2bdsim.parser import parse_g4bl


def test_fieldmap_writer_3d_format():
    xs, ys, zs = [-1.0, 0.0, 1.0], [-1.0, 1.0], [-2.0, 0.0, 2.0]
    text = fieldmap.build_3d(xs, ys, zs, lambda x, y, z: (0.0, 0.5, 0.0))
    assert text.startswith("xmin> -1")
    assert "nx> 3" in text and "ny> 2" in text and "nz> 3" in text
    assert "! X\tY\tZ\tFx\tFy\tFz" in text
    # loop order: X outer, Y mid, Z inner -> first 3 data rows share x,y, vary z
    rows = [l for l in text.splitlines() if l and (l[0].isdigit() or l[0] == "-")]
    assert len(rows) == 3 * 2 * 3


def test_gmad_field_object():
    s = fieldmap.gmad_field_object("f1", "f1.dat", 3, "cubic")
    assert s == ('f1: field, type="bmap3d", magneticFile="bdsim3d:f1.dat", '
                 'magneticInterpolator="cubic";')


def test_blfieldmap_parse_grid(tmp_path):
    p = tmp_path / "map.txt"
    lines = ["# test", "grid nX=3 nY=3 nZ=3 dX=10 dY=10 dZ=10 X0=-10 Y0=-10 Z0=-10",
             "data"]
    for ix in range(3):
        for iy in range(3):
            for iz in range(3):
                x, y, z = -10 + 10 * ix, -10 + 10 * iy, -10 + 10 * iz
                lines.append(f"{x} {y} {z} 0 0.3 0")
    p.write_text("\n".join(lines) + "\n")
    g = blfieldmap.parse(str(p))
    assert g is not None
    assert g.xs == [-10.0, 0.0, 10.0]
    assert g.field_fn()(0.0, 0.0, 0.0) == (0.0, 0.3, 0.0)


def test_blfieldmap_normB_scaling(tmp_path):
    p = tmp_path / "map.txt"
    p.write_text("param normB=2.0\ngrid nX=2 nY=2 nZ=2 dX=10 dY=10 dZ=10 "
                 "X0=0 Y0=0 Z0=0\ndata\n"
                 "0 0 0 0 0.1 0\n10 0 0 0 0.1 0\n0 10 0 0 0.1 0\n10 10 0 0 0.1 0\n"
                 "0 0 10 0 0.1 0\n10 0 10 0 0.1 0\n0 10 10 0 0.1 0\n10 10 10 0 0.1 0\n")
    g = blfieldmap.parse(str(p))
    assert g.field_fn()(0.0, 0.0, 0.0)[1] == 0.2   # 0.1 * normB(2.0)


def test_blfieldmap_cylinder_unsupported(tmp_path):
    p = tmp_path / "cyl.txt"
    p.write_text("cylinder nR=2 nZ=2 dR=10 dZ=10 Z0=0\ndata\n")
    assert blfieldmap.parse(str(p)) is None


def _convert(text, **kw):
    cmds, _ = parse_g4bl(text)
    return Converter(cmds, **kw).convert()


def test_fieldmap_command_to_map(tmp_path):
    mp = tmp_path / "m.txt"
    lines = ["grid nX=2 nY=2 nZ=3 dX=10 dY=10 dZ=10 X0=-10 Y0=-10 Z0=-10", "data"]
    for ix in range(2):
        for iy in range(2):
            for iz in range(3):
                lines.append(f"{-10+20*ix} {-10+20*iy} {-10+10*iz} 0 0.2 0")
    mp.write_text("\n".join(lines) + "\n")
    text = """
    reference referenceMomentum=1000 particle=proton
    fieldmap FM filename=m.txt
    place FM rename=FM1 z=500
    """
    m = _convert(text, base_dir=str(tmp_path))
    fm1 = m.get_element("FM1")
    assert fm1.type == "drift"
    assert fm1.params["fieldAll"] == "FM1_field"
    assert "FM1.dat" in m.aux_files


def test_solenoid_field_matches_g4beamline_profile():
    # The ported coil field must reproduce the on-axis Bz profile.
    f = solenoid.make_coil_field(150, 250, 600, 20)
    # values cross-checked against G4beamline printfield
    assert abs(f(0.0, 0.0)[1] - 2.0906) < 0.01
    assert abs(f(0.0, 300.0)[1] - 1.1913) < 0.02
    assert abs(f(0.0, 600.0)[1] - 0.1808) < 0.02
