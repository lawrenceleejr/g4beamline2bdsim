import pytest

from g4beamline2bdsim import gdml
from g4beamline2bdsim.converter import Converter
from g4beamline2bdsim.parser import parse_g4bl


def convert(text, **kw):
    cmds, _ = parse_g4bl(text)
    return Converter(cmds, **kw).convert()


def test_material_mapping():
    assert gdml.map_material("W") == ("G4_W", True)
    assert gdml.map_material("Cu") == ("G4_Cu", True)
    assert gdml.map_material("G4_STAINLESS-STEEL") == ("G4_STAINLESS-STEEL", True)
    name, known = gdml.map_material("Unobtanium")
    assert name == "G4_Unobtanium" and known is False


def test_is_vacuum():
    assert gdml.is_vacuum("Vacuum")
    assert gdml.is_vacuum(None)
    assert not gdml.is_vacuum("W")


def test_box_gdml_contains_solid_and_material():
    g = gdml.box_gdml(60, 60, 100, "G4_W")
    assert '<box name="target_solid"' in g
    assert 'x="60"' in g and 'z="100"' in g
    assert '<materialref ref="G4_W"/>' in g
    assert "<world ref=" in g


def test_tubs_gdml_tube():
    g = gdml.tubs_gdml(0, 5, 200, "G4_W")
    assert '<tube name="target_solid"' in g
    assert 'rmax="5"' in g and 'z="200"' in g


def test_box_material_becomes_gdml_element():
    text = """
    reference referenceMomentum=2000 particle=proton
    box T width=60 height=60 length=100 material=W
    place T rename=T1 z=500
    """
    m = convert(text)
    el = m.get_element("T1")
    assert el.type == "element"
    assert el.params["geometryFile"] == "gdml:T1.gdml"
    assert el.params["l"] == (100.0, "mm")
    assert "T1.gdml" in m.aux_files
    assert "G4_W" in m.aux_files["T1.gdml"]


def test_cylinder_material_becomes_gdml_element():
    text = """
    reference referenceMomentum=2000 particle=proton
    cylinder T outerRadius=5 length=200 material=W
    place T rename=T1 z=500
    """
    m = convert(text)
    el = m.get_element("T1")
    assert el.type == "element"
    assert "<tube" in m.aux_files["T1.gdml"]


def test_vacuum_box_stays_drift():
    text = """
    reference referenceMomentum=2000 particle=proton
    box T width=100 height=100 length=10 material=Vacuum
    place T rename=T1 z=500
    """
    m = convert(text)
    assert m.get_element("T1").type == "drift"
    assert not m.aux_files


def test_no_gdml_flag_falls_back_to_drift():
    text = """
    reference referenceMomentum=2000 particle=proton
    box T width=60 height=60 length=100 material=W
    place T rename=T1 z=500
    """
    m = convert(text, emit_gdml=False)
    assert m.get_element("T1").type == "drift"
    assert not m.aux_files
    assert any("drift" in w for w in m.warnings)


# --- fieldexpr auto-map tests ---
def test_fieldexpr_box_to_map():
    text = """
    reference referenceMomentum=1000 particle=proton
    fieldexpr F width=600 height=400 length=500 By=0.1
    place F rename=F1 z=1000
    """
    m = convert(text)
    el = m.get_element("F1")
    assert el.type == "drift"
    assert el.params["fieldAll"] == "F1_field"
    assert "F1.dat" in m.aux_files
    body = m.aux_files["F1.dat"]
    assert body.startswith("xmin>")
    # By ~ 0.1 T everywhere; check a data row has 0.1 in the Fy column.
    assert "1.00000000E-01" in body


def test_fieldexpr_varying_expression():
    text = """
    reference referenceMomentum=1000 particle=proton
    fieldexpr F width=200 height=200 length=400 By=0.2*cos(z/400)
    place F rename=F1 z=1000
    """
    m = convert(text)
    assert "F1.dat" in m.aux_files
    # field object declared before elements
    assert any("F1_field: field" in fo for fo in m.field_objects)


def test_fieldexpr_efield_warns():
    text = """
    reference referenceMomentum=1000 particle=proton
    fieldexpr F width=200 height=200 length=400 Ez=1.0
    place F rename=F1 z=1000
    """
    m = convert(text)
    assert any("electric-field" in w for w in m.warnings)
