import math

import pytest

from g4beamline2bdsim.expr import evaluate, has_operator
from g4beamline2bdsim.parser import parse_g4bl


def test_basic_arithmetic():
    assert evaluate("2+3*4") == 14
    assert evaluate("(2+3)*4") == 20
    assert evaluate("10/4") == 2.5


def test_power_caret():
    assert evaluate("2^10") == 1024
    assert evaluate("3^2 + 1") == 10


def test_functions():
    assert evaluate("sqrt(16)") == 4
    assert evaluate("pow(2,8)") == 256
    assert evaluate("abs(-5)") == 5
    assert evaluate("if(1>0, 7, 9)") == 7
    assert evaluate("if(0, 7, 9)") == 9


def test_constants():
    assert evaluate("pi") == pytest.approx(math.pi)
    assert evaluate("180*degree") == pytest.approx(math.pi)


def test_non_numeric_returns_none():
    assert evaluate("copper") is None
    assert evaluate("$P") is None          # unresolved reference
    assert evaluate("") is None
    assert evaluate("1 0 0") is None       # color triple, not an expression


def test_no_code_execution():
    # Must not evaluate arbitrary Python.
    assert evaluate("__import__('os')") is None
    assert evaluate("open('x')") is None


def test_has_operator():
    assert has_operator("2+3")
    assert has_operator("sqrt(2)")
    assert not has_operator("938.272")
    assert not has_operator("proton")


def test_param_expression_evaluated_in_parser():
    text = "param KE=200\nparam M=938.272\n" \
           "param P=sqrt(($KE+$M)*($KE+$M)-$M*$M)\n" \
           "beam gaussian particle=proton meanMomentum=$P\n"
    cmds, params = parse_g4bl(text)
    expected = math.sqrt((200 + 938.272) ** 2 - 938.272 ** 2)
    assert float(params["P"]) == pytest.approx(expected)
    beam = [c for c in cmds if c.name == "beam"][0]
    assert float(beam.params["meanMomentum"]) == pytest.approx(expected)


def test_param_arithmetic_in_placement():
    text = "param G=4\ngenericquad Q fieldLength=300\n" \
           "place Q rename=Q1 gradient=$G*1.5 z=$G*100\n"
    cmds, _ = parse_g4bl(text)
    place = [c for c in cmds if c.name == "place"][0]
    # gradient is substituted to '4*1.5' (an arg expr, evaluated by converter)
    assert place.params["gradient"] in ("4*1.5", "6", "6.0")
