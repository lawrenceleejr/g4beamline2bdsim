from g4beamline2bdsim.converter import convert_commands
from g4beamline2bdsim.gmad_writer import GmadWriter, _format_value
from g4beamline2bdsim.parser import parse_g4bl


def render(text):
    cmds, _ = parse_g4bl(text)
    model = convert_commands(cmds)
    return GmadWriter(model).to_string()


def test_format_value_units():
    assert _format_value((5.0, "mm")) == "5*mm"
    assert _format_value((0.201, "GHz")) == "0.201*GHz"
    assert _format_value("copper") == '"copper"'
    assert _format_value(1000) == "1000"


def test_full_render_is_valid_gmad_shape():
    text = """
    physics QGSP_BERT
    reference referenceMomentum=200 particle=mu+
    beam gaussian particle=mu+ nEvents=10 meanMomentum=200 sigmaX=5
    genericquad Q fieldLength=300 apertureRadius=100
    place Q rename=Q1 gradient=4.5 z=300
    virtualdetector V radius=100 length=1
    place V rename=Det z=1000
    """
    out = render(text)
    # Every statement ends with a semicolon (ignoring inline ! comments).
    for line in out.splitlines():
        s = line.strip()
        if not s or s.startswith("!"):
            continue
        code = s.split("!", 1)[0].strip()
        if not code:
            continue
        assert code.endswith(";") or code.endswith(",") or code.endswith("("), line
    assert "use, mainline;" in out
    assert "beam," in out
    assert "option," in out
    assert "sample, range=Det;" in out
    assert "Q1: quadrupole" in out


def test_long_line_wraps():
    # Build a long beamline and ensure the line= wraps without breaking tokens.
    parts = ["reference referenceMomentum=200 particle=mu+",
             "genericquad Q fieldLength=100"]
    for i in range(20):
        parts.append(f"place Q rename=Q{i} gradient=1 z={200 + i*300}")
    out = render("\n".join(parts))
    assert "line=(" in out
    # No line should exceed a sane width by a huge margin.
    assert all(len(l) < 200 for l in out.splitlines())
