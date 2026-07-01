"""Safe evaluator for G4beamline numeric expressions.

G4beamline evaluates arithmetic in argument values, e.g.::

    param P=sqrt(($KE+$M)*($KE+$M)-$M*$M)
    place Q gradient=$G*1.5 z=$L/2

Supported (per the G4beamline User's Guide, sec. 5.1):

* operators ``+ - * / ^`` where ``^`` is power, and parentheses;
* comparisons ``< <= > >= == !=`` -> 1.0/0.0; logical ``&& ||``;
* functions ``abs min max sqrt pow sin cos tan asin acos atan atan2 sinh
  cosh tanh exp log log10 floor ceil if``;
* constants ``pi e gamma radian rad degree deg``.

Evaluation is done with a restricted ``ast`` walk -- no names, attributes or
calls outside the whitelist are allowed, so arbitrary code cannot run.
"""

from __future__ import annotations

import ast
import math
import operator
import re
from typing import Optional

_CONSTANTS = {
    "pi": math.pi,
    "e": math.e,
    "gamma": 0.5772156649015329,   # Euler-Mascheroni, as in G4beamline
    "radian": 1.0,
    "rad": 1.0,
    "degree": math.pi / 180.0,
    "deg": math.pi / 180.0,
}


def _g4_if(cond, a, b):
    return a if cond else b


_FUNCTIONS = {
    "abs": abs,
    "min": min,
    "max": max,
    "sqrt": math.sqrt,
    "pow": math.pow,
    "sin": math.sin,
    "cos": math.cos,
    "tan": math.tan,
    "asin": math.asin,
    "acos": math.acos,
    "atan": math.atan,
    "atan2": math.atan2,
    "sinh": math.sinh,
    "cosh": math.cosh,
    "tanh": math.tanh,
    "exp": math.exp,
    "log": math.log,
    "log10": math.log10,
    "floor": math.floor,
    "ceil": math.ceil,
    "if": _g4_if,
    "_g4if": _g4_if,    # 'if' is a Python keyword; rewritten before parsing
}

_BINOPS = {
    ast.Add: operator.add,
    ast.Sub: operator.sub,
    ast.Mult: operator.mul,
    ast.Div: operator.truediv,
    ast.Pow: operator.pow,        # '^' is rewritten to '**' before parsing
    ast.Mod: operator.mod,
    ast.BitXor: operator.pow,     # in case a raw '^' slips through
}

_UNARYOPS = {
    ast.UAdd: operator.pos,
    ast.USub: operator.neg,
    ast.Not: lambda x: 1.0 if not x else 0.0,
}

_CMPOPS = {
    ast.Lt: operator.lt,
    ast.LtE: operator.le,
    ast.Gt: operator.gt,
    ast.GtE: operator.ge,
    ast.Eq: operator.eq,
    ast.NotEq: operator.ne,
}


class ExprError(Exception):
    pass


def _eval(node, variables):
    if isinstance(node, ast.Expression):
        return _eval(node.body, variables)
    if isinstance(node, ast.Constant):
        if isinstance(node.value, bool):
            return 1.0 if node.value else 0.0
        if isinstance(node.value, (int, float)):
            return float(node.value)
        raise ExprError(f"non-numeric constant: {node.value!r}")
    if isinstance(node, ast.Name):
        if variables and node.id in variables:
            return float(variables[node.id])
        if node.id in _CONSTANTS:
            return _CONSTANTS[node.id]
        raise ExprError(f"unknown name: {node.id}")
    if isinstance(node, ast.BinOp):
        op = _BINOPS.get(type(node.op))
        if op is None:
            raise ExprError(f"unsupported operator: {type(node.op).__name__}")
        return op(_eval(node.left, variables), _eval(node.right, variables))
    if isinstance(node, ast.UnaryOp):
        op = _UNARYOPS.get(type(node.op))
        if op is None:
            raise ExprError("unsupported unary operator")
        return op(_eval(node.operand, variables))
    if isinstance(node, ast.BoolOp):
        vals = [_eval(v, variables) for v in node.values]
        if isinstance(node.op, ast.And):
            return 1.0 if all(vals) else 0.0
        return 1.0 if any(vals) else 0.0
    if isinstance(node, ast.Compare):
        if len(node.ops) != 1:
            raise ExprError("chained comparisons not supported")
        op = _CMPOPS.get(type(node.ops[0]))
        if op is None:
            raise ExprError("unsupported comparison")
        return (1.0 if op(_eval(node.left, variables),
                          _eval(node.comparators[0], variables)) else 0.0)
    if isinstance(node, ast.Call):
        if not isinstance(node.func, ast.Name):
            raise ExprError("only named functions allowed")
        fn = _FUNCTIONS.get(node.func.id)
        if fn is None:
            raise ExprError(f"unknown function: {node.func.id}")
        return float(fn(*[_eval(a, variables) for a in node.args]))
    raise ExprError(f"unsupported syntax: {type(node).__name__}")


def _prepare(text: str) -> str:
    s = text.strip()
    s = s.replace("^", "**")
    # 'if' is a Python keyword -> rewrite the G4beamline if(...) call.
    return re.sub(r"\bif\s*\(", "_g4if(", s)


def evaluate(text: str, variables: Optional[dict] = None) -> Optional[float]:
    """Evaluate *text* as a G4beamline numeric expression.

    ``variables`` optionally supplies values for names (e.g. ``x``/``y``/``z``).
    Returns the float value, or ``None`` if it is not a valid numeric
    expression (a string literal, a material name, an unresolved ``$ref``).
    """
    if text is None:
        return None
    s = text.strip()
    if not s:
        return None
    if "$" in s:               # unresolved parameter reference
        return None
    try:
        tree = ast.parse(_prepare(s), mode="eval")
        return float(_eval(tree, variables))
    except (ExprError, SyntaxError, ValueError, TypeError, ZeroDivisionError,
            OverflowError):
        return None


def compile_expr(text: str):
    """Compile a G4beamline expression to a callable ``f(**variables)``.

    Raises :class:`ExprError` if the text is not a valid expression.  The
    returned callable evaluates fast (AST is parsed once).
    """
    if text is None:
        raise ExprError("empty expression")
    tree = ast.parse(_prepare(text), mode="eval")

    def fn(**variables):
        return float(_eval(tree, variables))

    return fn


def has_operator(text: str) -> bool:
    """True if *text* looks like an expression (contains an operator/func)."""
    if text is None:
        return False
    return any(c in text for c in "+-*/^()") or any(
        fn + "(" in text for fn in _FUNCTIONS)
