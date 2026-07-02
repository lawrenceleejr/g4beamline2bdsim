"""Tests for G4beamline control flow: do/enddo, if/elseif/else/endif,
define macros, and include."""

import math

import pytest

from g4beamline2bdsim.parser import G4BLParser, parse_g4bl


def names(cmds):
    return [c.name for c in cmds]


def test_do_loop_expands():
    text = """
    genericquad Q fieldLength=100
    do i 0 2
      place Q rename=Q$i gradient=1 z=$i*300+200
    enddo
    """
    cmds, params = parse_g4bl(text)
    places = [c for c in cmds if c.name == "place"]
    assert len(places) == 3
    assert [p.params["rename"] for p in places] == ["Q0", "Q1", "Q2"]
    # z carries the substituted loop variable, evaluated downstream
    assert places[1].params["z"] in ("1*300+200", "500", "500.0")


def test_do_loop_with_expression_bounds():
    text = """
    param N=4
    do i 0 $N-1
      place X z=$i
    enddo
    """
    cmds, _ = parse_g4bl(text)
    assert len([c for c in cmds if c.name == "place"]) == 4


def test_do_loop_increment_and_negative():
    text = "do i 10 0 -5\nplace X z=$i\nenddo\n"
    cmds, _ = parse_g4bl(text)
    zs = [c.params["z"] for c in cmds if c.name == "place"]
    assert zs == ["10", "5", "0"]


def test_nested_do():
    text = """
    do i 1 2
    do j 1 2
    place X rename=P$i$j
    enddo
    enddo
    """
    cmds, _ = parse_g4bl(text)
    renames = [c.params["rename"] for c in cmds if c.name == "place"]
    assert renames == ["P11", "P12", "P21", "P22"]


def test_if_block_true_false():
    text = """
    param A=1
    if $A==1
      place X rename=yes
    else
      place X rename=no
    endif
    """
    cmds, _ = parse_g4bl(text)
    assert [c.params["rename"] for c in cmds if c.name == "place"] == ["yes"]


def test_if_elseif_else():
    text = """
    param A=2
    if $A==1
      place X rename=one
    elseif $A==2
      place X rename=two
    else
      place X rename=other
    endif
    """
    cmds, _ = parse_g4bl(text)
    assert [c.params["rename"] for c in cmds if c.name == "place"] == ["two"]


def test_if_else_fallthrough():
    text = """
    param A=9
    if $A==1
      place X rename=one
    elseif $A==2
      place X rename=two
    else
      place X rename=other
    endif
    """
    cmds, _ = parse_g4bl(text)
    assert [c.params["rename"] for c in cmds if c.name == "place"] == ["other"]


def test_inline_if():
    text = 'param A=1\nif $A==1 "place X rename=inline"\n'
    cmds, _ = parse_g4bl(text)
    assert [c.params.get("rename") for c in cmds if c.name == "place"] == ["inline"]


def test_if_inside_do():
    # The official multipole.g4bl pattern: if/else inside a do loop.
    text = """
    param N=3
    do i 0 $N-1
      if $i==$N-1
        param exit=1
      else
        param exit=0
      endif
      place X rename=P$i tag=$exit
    enddo
    """
    cmds, params = parse_g4bl(text)
    tags = [c.params["tag"] for c in cmds if c.name == "place"]
    assert tags == ["0", "0", "1"]


def test_define_macro_positional_args():
    text = """
    define Cell "place QF rename=$1_F z=$2" "place QD rename=$1_D z=$2+500"
    Cell CellA 1000
    Cell CellB 3000
    """
    cmds, _ = parse_g4bl(text)
    places = [c for c in cmds if c.name == "place"]
    assert len(places) == 4
    assert places[0].params["rename"] == "CellA_F"
    assert places[0].params["z"] == "1000"
    assert places[3].params["rename"] == "CellB_D"
    assert places[3].params["z"] in ("3000+500", "3500")


def test_define_macro_expansion_counter():
    text = """
    define M "place X rename=U$#"
    M
    M
    """
    cmds, _ = parse_g4bl(text)
    renames = [c.params["rename"] for c in cmds if c.name == "place"]
    assert renames == ["U1", "U2"]


def test_define_dollar_dollar_deferred():
    # $$P is expanded at invocation time, $P at define time.
    text = """
    param P=early
    define M "place X a=$P b=$$P"
    param P=late
    M
    """
    cmds, _ = parse_g4bl(text)
    place = [c for c in cmds if c.name == "place"][0]
    assert place.params["a"] == "early"
    assert place.params["b"] == "late"


def test_include(tmp_path):
    inc = tmp_path / "common.g4bl"
    inc.write_text("genericquad Q fieldLength=100\n")
    main = "include common.g4bl\nplace Q z=100\n"
    parser = G4BLParser(base_dir=str(tmp_path))
    cmds = parser.parse(main)
    assert names(cmds) == ["genericquad", "place"]


def test_include_missing_warns():
    parser = G4BLParser(base_dir="/nonexistent")
    parser.parse("include nothere.g4bl\n")
    assert any("not found" in w for w in parser.warnings)


def test_star_and_ui_lines_skipped():
    text = "* banner comment\n/vis/viewer refresh\n!ls\nplace X z=1\n"
    cmds, _ = parse_g4bl(text)
    assert names(cmds) == ["place"]


def test_runaway_loop_guard():
    # A pathological giant loop must raise, not hang.
    text = "do i 0 999999\nplace X z=$i\nenddo\n"
    from g4beamline2bdsim.parser import G4BLParseError
    with pytest.raises(G4BLParseError):
        parse_g4bl(text)
