import math

import pytest

from g4beamline2bdsim.converter import Converter, convert_commands
from g4beamline2bdsim.parser import parse_g4bl


def convert(text):
    cmds, _ = parse_g4bl(text)
    return Converter(cmds).convert()


def get(model, name):
    return model.get_element(name)


def test_quad_k1_from_gradient():
    text = """
    reference referenceMomentum=200 particle=mu+
    genericquad Q fieldLength=300 apertureRadius=100
    place Q rename=QF gradient=4.5 z=300
    """
    m = convert(text)
    qf = get(m, "QF")
    assert qf.type == "quadrupole"
    brho = 0.2 / 0.299792458
    assert qf.params["k1"] == pytest.approx(4.5 / brho)
    # aperture preserved
    assert qf.params["aper1"] == (100.0, "mm")
    assert qf.params["apertureType"] == "circular"


def test_quad_length_in_mm():
    text = """
    reference referenceMomentum=1000 particle=proton
    genericquad Q fieldLength=250
    place Q rename=Q1 gradient=1 z=200
    """
    m = convert(text)
    assert get(m, "Q1").params["l"] == (250.0, "mm")


def test_drifts_fill_gaps():
    text = """
    reference referenceMomentum=1000 particle=proton
    genericquad Q fieldLength=200
    place Q rename=Q1 gradient=1 z=300
    place Q rename=Q2 gradient=-1 z=1000
    """
    m = convert(text)
    # Q1 center 300, len 200 -> entry 200, exit 400.  Leading drift 0..200.
    # Q2 center 1000 -> entry 900.  Gap drift 400..900 = 500 mm.
    names = m.line
    assert names[0].startswith("drift")          # leading drift
    assert "Q1" in names and "Q2" in names
    drifts = [get(m, n) for n in names if n.startswith("drift")]
    lengths = sorted(d.params["l"][0] for d in drifts)
    assert 200.0 in lengths   # leading drift
    assert 500.0 in lengths   # gap between Q1 and Q2


def test_genericbend_to_rbend():
    text = """
    reference referenceMomentum=1000 particle=proton
    genericbend D fieldLength=500
    place D rename=B1 By=0.8 z=1000
    """
    m = convert(text)
    b1 = get(m, "B1")
    assert b1.type == "rbend"
    assert b1.params["B"] == (0.8, "T")
    assert b1.params["l"] == (500.0, "mm")


def test_idealsectorbend_angle():
    text = """
    reference referenceMomentum=1000 particle=proton
    idealsectorbend D angle=15 By=0.5 fieldCenterRadius=2000
    place D z=4000
    """
    m = convert(text)
    d = get(m, "D")
    assert d.type == "sbend"
    assert d.params["angle"] == pytest.approx(math.radians(15))
    # arc length = R * angle_rad
    assert d.params["l"][0] == pytest.approx(2000 * math.radians(15))


def test_sectorbend_front_face_positioning():
    # idealsectorbend z marks the FRONT face, not the centre.
    text = """
    reference referenceMomentum=1000 particle=proton
    idealsectorbend D angle=10 By=0.5 fieldCenterRadius=1000
    place D z=2000
    virtualdetector V radius=10 length=1
    place V rename=V1 z=5000
    """
    m = convert(text)
    arc = 1000 * math.radians(10)
    # The drift before D should be exactly 2000 mm (D starts at z=2000).
    lengths = [get(m, n).params["l"][0] for n in m.line if n.startswith("drift")]
    assert any(abs(l - 2000.0) < 1e-6 for l in lengths)


def test_pillbox_to_rfcavity():
    text = """
    reference referenceMomentum=200 particle=mu+
    pillbox RF maxGradient=16 frequency=0.201 innerLength=430 phaseAcc=40
    place RF rename=RF1 z=500
    """
    m = convert(text)
    rf = get(m, "RF1")
    assert rf.type == "rfcavity"
    assert rf.params["frequency"] == (0.201, "GHz")
    # E = gradient[MV/m] * length[m]
    assert rf.params["E"][0] == pytest.approx(16 * 0.430)
    assert rf.params["phase"] == pytest.approx(math.radians(40))


def test_multipole_combined_function():
    text = """
    reference referenceMomentum=200 particle=mu+
    multipole MP fieldLength=200 quadrupole=2.0 sextupole=50.0
    place MP rename=M1 z=300
    """
    m = convert(text)
    m1 = get(m, "M1")
    assert m1.type == "multipole"
    brho = 0.2 / 0.299792458
    k1l = (2.0 / brho) * 0.2
    k2l = (2.0 * 50.0 / brho) * 0.2
    assert f"{k1l:.10g}" in m1.params["knl"]
    assert f"{k2l:.10g}" in m1.params["knl"]


def test_single_sextupole():
    text = """
    reference referenceMomentum=200 particle=mu+
    multipole MP fieldLength=200 sextupole=50.0
    place MP rename=S1 z=300
    """
    m = convert(text)
    s1 = get(m, "S1")
    assert s1.type == "sextupole"
    brho = 0.2 / 0.299792458
    assert s1.params["k2"] == pytest.approx(2.0 * 50.0 / brho)


def test_solenoid_field_estimate():
    text = """
    reference referenceMomentum=200 particle=mu+
    coil C1 innerRadius=300 outerRadius=400 length=500 material=Cu
    solenoid S coilName=C1 current=80
    place S rename=S1 z=500
    """
    m = convert(text)
    s1 = get(m, "S1")
    assert s1.type == "solenoid"
    assert s1.params["B"][0] == pytest.approx(5.857, abs=0.01)
    assert any("solenoid" in w for w in m.warnings)


def test_beam_and_options():
    text = """
    physics QGSP_BERT
    reference referenceMomentum=200 particle=mu+
    beam gaussian particle=mu+ nEvents=1000 meanMomentum=200 sigmaX=5 sigmaXp=0.002 sigmaP=2
    genericquad Q fieldLength=300
    place Q rename=Q1 gradient=1 z=300
    """
    m = convert(text)
    assert m.beam["particle"] == "mu+"
    assert m.beam["momentum"] == (200.0, "MeV")
    assert m.beam["distrType"] == "gauss"
    assert m.beam["sigmaX"] == (5.0, "mm")
    assert m.beam["sigmaE"] == pytest.approx(2.0 / 200.0)
    assert m.options["physicsList"] == "g4QGSP_BERT"
    assert m.options["ngenerate"] == 1000


def test_virtualdetector_becomes_sampler():
    text = """
    reference referenceMomentum=200 particle=mu+
    virtualdetector V radius=100 length=1
    place V rename=Det z=500
    """
    m = convert(text)
    assert "Det" in m.samplers
    assert get(m, "Det").type == "marker"


def test_missing_momentum_warns():
    text = """
    genericquad Q fieldLength=300
    place Q rename=Q1 gradient=4.5 z=300
    """
    m = convert(text)
    assert any("momentum" in w for w in m.warnings)
    # k1 not set, gradient preserved in a comment
    assert get(m, "Q1").comment is not None


def test_rename_hash_numbering():
    # rename=Det# -> Det1, Det2, ... (1-based placement number).
    text = """
    reference referenceMomentum=200 particle=mu+
    virtualdetector Det radius=100 length=1
    place Det z=1000 rename=Det#
    place Det z=2000 rename=Det#
    place Det z=3000 rename=Det#
    """
    m = convert(text)
    assert "Det1" in m.samplers
    assert "Det2" in m.samplers
    assert "Det3" in m.samplers


def test_beampipe_sized_from_aperture():
    text = """
    reference referenceMomentum=1000 particle=proton
    genericquad Q fieldLength=300 apertureRadius=120 ironRadius=400
    place Q rename=Q1 gradient=1 z=300
    virtualdetector D radius=200 length=1
    place D rename=Det z=600
    """
    m = convert(text)
    # beampipe radius takes the largest genuine aperture (detector radius=200),
    # NOT the iron outer radius (400).
    assert m.options["beampipeRadius"] == (200.0, "mm")
    # horizontalWidth must exceed 2*beampipeRadius.
    hw = m.options["horizontalWidth"][0]
    assert hw > 2 * 200.0


def test_sectorbend_arc_radius_not_used_as_aperture():
    # fieldOuterRadius is the bend arc radius, not a transverse aperture.
    text = """
    reference referenceMomentum=1000 particle=proton
    idealsectorbend B angle=10 By=0.5 fieldCenterRadius=2000 \
        fieldOuterRadius=2200 fieldHeight=200
    place B z=2000
    """
    m = convert(text)
    # Only fieldHeight/2 (=100) should drive the beampipe, not 2200.
    assert m.options["beampipeRadius"][0] <= 100.0


def test_invalid_physics_list_fallback():
    text = """
    physics QGSP
    reference referenceMomentum=200 particle=mu+
    genericquad Q fieldLength=300
    place Q rename=Q1 gradient=1 z=300
    """
    m = convert(text)
    assert m.options["physicsList"] == "g4FTFP_BERT"
    assert any("physics list" in w for w in m.warnings)


def test_valid_physics_list_prefixed():
    text = """
    physics QGSP_BERT
    reference referenceMomentum=200 particle=mu+
    genericquad Q fieldLength=300
    place Q rename=Q1 gradient=1 z=300
    """
    m = convert(text)
    assert m.options["physicsList"] == "g4QGSP_BERT"


def test_param_expression_drives_k1():
    # gradient given as an expression referencing a param.
    text = """
    param G=3.0
    reference referenceMomentum=1000 particle=proton
    genericquad Q fieldLength=300
    place Q rename=Q1 gradient=$G*2 z=300
    """
    m = convert(text)
    brho = 1.0 / 0.299792458
    assert get(m, "Q1").params["k1"] == pytest.approx(6.0 / brho)


def test_negative_particle_charge_rigidity():
    # electrons: rigidity uses |q|, so k1 stays positive for positive gradient.
    text = """
    reference referenceMomentum=1000 particle=e-
    genericquad Q fieldLength=300
    place Q rename=Q1 gradient=3 z=300
    """
    m = convert(text)
    brho = 1.0 / 0.299792458
    assert get(m, "Q1").params["k1"] == pytest.approx(3.0 / brho)
