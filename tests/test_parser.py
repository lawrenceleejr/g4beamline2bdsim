import math

from g4beamline2bdsim.parser import parse_g4bl


def test_basic_command():
    cmds, params = parse_g4bl("genericquad Q fieldLength=300 gradient=4.5")
    assert len(cmds) == 1
    c = cmds[0]
    assert c.name == "genericquad"
    assert c.args == ["Q"]
    assert c.params["fieldLength"] == "300"
    assert c.params["gradient"] == "4.5"


def test_comments_and_blank_lines():
    text = "# a comment\n\ngenericquad Q  # trailing comment\nplace Q z=10\n"
    cmds, _ = parse_g4bl(text)
    assert [c.name for c in cmds] == ["genericquad", "place"]


def test_param_substitution():
    text = "param P=200\nbeam gaussian meanMomentum=$P\n"
    cmds, params = parse_g4bl(text)
    assert params["P"] == "200"
    beam = [c for c in cmds if c.name == "beam"][0]
    assert beam.params["meanMomentum"] == "200"


def test_param_unset_does_not_override():
    text = "param P=200\nparam -unset P=999\n"
    _, params = parse_g4bl(text)
    assert params["P"] == "200"


def test_line_continuation():
    text = "beam gaussian particle=mu+ \\\n     meanMomentum=200 sigmaX=5\n"
    cmds, _ = parse_g4bl(text)
    assert len(cmds) == 1
    assert cmds[0].params["sigmaX"] == "5"
    assert cmds[0].params["meanMomentum"] == "200"


def test_quoted_value():
    cmds, _ = parse_g4bl('box B material="G4_Cu" color="1 0 0"')
    assert cmds[0].params["material"] == "G4_Cu"
    assert cmds[0].params["color"] == "1 0 0"


def test_nested_param_reference():
    text = "param A=2\nparam B=$A\nbox X width=$B\n"
    cmds, params = parse_g4bl(text)
    assert params["B"] == "2"
    box = [c for c in cmds if c.name == "box"][0]
    assert box.params["width"] == "2"
