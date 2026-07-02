"""Parser/interpreter for G4beamline input files.

G4beamline input is a line-oriented command language.  Each (logical) line is a
command consisting of a command name followed by positional arguments and
``name=value`` keyword arguments::

    command arg1 arg2 name1=value1 name2=value2

Key syntactic features handled here:

* ``#`` introduces a comment at line start or after whitespace (a ``#``
  embedded in a token, e.g. ``rename=Det#``, is not a comment).  ``*`` at line
  start is a comment echoed by g4bl; ``/`` and ``!`` lines are Geant4-UI /
  shell escapes (ignored here).
* A backslash ``\\`` at the end of a line continues the command.
* ``param name=value ...`` defines parameters referenced with ``$name`` (or
  ``${name}``); a value containing an operator is evaluated numerically.
  ``param -unset`` only sets a value if not already defined.

Control flow (semantics verified against the G4beamline 3.06 source —
``BLCMDdo.cc``, ``BLCMDif.cc``, ``BLCMDdefine.cc``, ``BLCMDinclude.cc``):

* ``do i first last [incr]`` … ``enddo`` — loop, setting param ``i``.
* ``if <expr>`` … ``elseif <expr>`` … ``else`` … ``endif`` — block form; and
  the inline form ``if <expr> <cmd> [<cmd> …]`` executing each extra argument
  as a command.
* ``define name <line1> <line2> …`` — macro; ``$0``–``$9`` and ``$#`` are
  substituted at invocation, ``$param`` at define time, ``$$param`` at
  invocation.
* ``include <file>`` — reads another input file (resolved against
  ``base_dir``).

The parser is intentionally permissive: anything it does not understand is kept
verbatim so the converter can decide how to handle (or warn about) it.
"""

from __future__ import annotations

import os
import re
import shlex
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .expr import evaluate as _eval_expr, has_operator as _has_operator


@dataclass
class G4BLCommand:
    """A single parsed G4beamline command."""

    name: str
    args: List[str] = field(default_factory=list)
    params: Dict[str, str] = field(default_factory=dict)
    # Preserve order of keyword params for nicer round-tripping / messages.
    param_order: List[str] = field(default_factory=list)
    line_no: int = 0
    raw: str = ""

    def get(self, key: str, default: Optional[str] = None) -> Optional[str]:
        return self.params.get(key, default)

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        parts = [self.name, *self.args]
        parts += [f"{k}={self.params[k]}" for k in self.param_order]
        return "G4BLCommand(" + " ".join(parts) + ")"


_PARAM_REF = re.compile(r"\$\{?([A-Za-z_][A-Za-z0-9_]*)\}?")
_MACRO_REF = re.compile(r"\$([0-9#])")
# Sentinel protecting '$$name' (deferred expansion) during '$name' expansion.
_DOLLAR_SENTINEL = "\x00DOLLAR\x00"

# Hard cap on executed logical lines - guards against runaway do-loops.
_MAX_EXECUTED_LINES = 500_000


def _format_eval(value: float) -> str:
    """Render an evaluated number compactly (integers without a trailing .0)."""
    if value == int(value) and abs(value) < 1e15:
        return str(int(value))
    return repr(value)


class G4BLParseError(Exception):
    pass


class G4BLParser:
    """Interpret G4beamline text into a list of :class:`G4BLCommand`."""

    def __init__(self, base_dir: str = "") -> None:
        self.base_dir = base_dir
        self.params: Dict[str, str] = {}
        self.commands: List[G4BLCommand] = []
        self.macros: Dict[str, List[str]] = {}
        self.warnings: List[str] = []
        self._n_expansions = 0
        self._executed = 0

    # -- public API ---------------------------------------------------------
    def parse(self, text: str) -> List[G4BLCommand]:
        lines = self._logical_lines(text)
        self._run(lines)
        return self.commands

    # -- lexing -------------------------------------------------------------
    @staticmethod
    def _strip_comment(line: str) -> str:
        """Remove a ``#`` comment, respecting quotes.

        ``#`` starts a comment only at the start of a line (after optional
        whitespace) or when preceded by whitespace - a ``#`` embedded in a
        token (e.g. ``rename=Det#``) is *not* a comment.
        """
        out = []
        quote = None
        prev = ""
        for i, ch in enumerate(line):
            if quote:
                out.append(ch)
                if ch == quote:
                    quote = None
            else:
                if ch in ("'", '"'):
                    quote = ch
                    out.append(ch)
                elif ch == "#" and (i == 0 or prev in (" ", "\t")):
                    break
                else:
                    out.append(ch)
            prev = ch
        return "".join(out)

    def _logical_lines(self, text: str) -> List[Tuple[int, str]]:
        """Return (line_number, logical_line), merging backslash continuations."""
        result: List[Tuple[int, str]] = []
        buffer = ""
        start_no = 0
        for idx, raw in enumerate(text.splitlines(), start=1):
            stripped = raw.lstrip()
            # '*' comment lines, Geant4-UI '/' lines and shell '!' lines.
            if buffer == "" and stripped[:1] in ("*", "/", "!"):
                continue
            line = self._strip_comment(raw)
            if buffer == "":
                start_no = idx
            if line.rstrip().endswith("\\"):
                buffer += line.rstrip()[:-1] + " "
                continue
            buffer += line
            result.append((start_no, buffer))
            buffer = ""
        if buffer:
            result.append((start_no, buffer))
        return result

    @staticmethod
    def _tokenize(line: str) -> List[str]:
        try:
            return shlex.split(line, posix=True)
        except ValueError:
            # Unbalanced quotes etc. - fall back to a naive split.
            return line.split()

    # -- execution ----------------------------------------------------------
    def _run(self, lines: List[Tuple[int, str]]) -> None:
        i = 0
        n = len(lines)
        while i < n:
            self._executed += 1
            if self._executed > _MAX_EXECUTED_LINES:
                raise G4BLParseError(
                    "input executes too many lines (runaway do loop?)")
            line_no, text = lines[i]
            stripped = text.strip()
            if not stripped:
                i += 1
                continue
            head = stripped.split(None, 1)[0]

            if head == "do":
                i = self._exec_do(lines, i)
            elif head == "if":
                i = self._exec_if(lines, i)
            elif head in ("enddo", "endif", "else") or head.startswith("elseif"):
                # Stray terminator (should be consumed by do/if handlers).
                self.warnings.append(
                    f"unexpected '{head}' at line {line_no}; ignored")
                i += 1
            elif head == "define":
                self._exec_define(stripped, line_no)
                i += 1
            elif head == "include":
                self._exec_include(stripped, line_no)
                i += 1
            elif head in self.macros:
                self._exec_macro(stripped, line_no)
                i += 1
            else:
                cmd = self._parse_line(stripped, line_no)
                if cmd is not None:
                    if cmd.name == "param":
                        self._handle_param(cmd)
                    self.commands.append(cmd)
                i += 1

    # -- do ... enddo --------------------------------------------------------
    @staticmethod
    def _is_block_open(stripped: str, kind: str,
                       macros: Dict[str, List[str]]) -> bool:
        head = stripped.split(None, 1)[0] if stripped.split() else ""
        if kind == "do":
            return head == "do"
        if kind == "if":
            # Only the BLOCK form of if (exactly one argument) opens a region.
            if head != "if":
                return False
            toks = G4BLParser._tokenize(stripped)
            return len(toks) == 2
        return False

    def _find_matching(self, lines, start: int, kind: str) -> int:
        """Index of the matching ``enddo``/``endif`` for the block at start."""
        closer = "enddo" if kind == "do" else "endif"
        level = 1
        j = start + 1
        while j < len(lines):
            s = lines[j][1].strip()
            if not s:
                j += 1
                continue
            head = s.split(None, 1)[0]
            if self._is_block_open(s, kind, self.macros):
                level += 1
            elif head == closer:
                level -= 1
                if level == 0:
                    return j
            j += 1
        raise G4BLParseError(
            f"unterminated '{kind}' block starting at line {lines[start][0]}")

    def _exec_do(self, lines, i: int) -> int:
        line_no, text = lines[i]
        toks = self._tokenize(text.strip())
        # do var first last [incr]
        if len(toks) < 4 or len(toks) > 5:
            self.warnings.append(f"invalid do command at line {line_no}")
            return i + 1
        var = toks[1]
        spec = [self._eval_int(self.expand(t), line_no) for t in toks[2:]]
        first, last = spec[0], spec[1]
        incr = spec[2] if len(spec) > 2 else 1
        if incr == 0:
            self.warnings.append(f"do increment 0 at line {line_no}; skipped")
            incr = 1
            first, last = 1, 0
        end = self._find_matching(lines, i, "do")
        body = lines[i + 1:end]
        value = first
        while (value <= last) if incr > 0 else (value >= last):
            self.params[var] = _format_eval(float(value))
            self._run(body)
            value += incr
        return end + 1

    def _eval_int(self, token: str, line_no: int) -> int:
        result = _eval_expr(token)
        if result is None:
            try:
                return int(float(token))
            except (TypeError, ValueError):
                self.warnings.append(
                    f"cannot evaluate '{token}' at line {line_no}; using 0")
                return 0
        return int(result)

    # -- if / elseif / else / endif -------------------------------------------
    def _truthy(self, expr_text: str, line_no: int) -> bool:
        value = _eval_expr(self.expand(expr_text))
        if value is None:
            self.warnings.append(
                f"cannot evaluate condition '{expr_text}' at line {line_no}; "
                f"treated as false")
            return False
        return abs(value) > 1.0e-12

    def _exec_if(self, lines, i: int) -> int:
        line_no, text = lines[i]
        toks = self._tokenize(text.strip())
        if len(toks) < 2:
            self.warnings.append(f"invalid if command at line {line_no}")
            return i + 1
        if len(toks) > 2:
            # Inline form: if <expr> <cmd> [<cmd> ...]
            if self._truthy(toks[1], line_no):
                for cmd_text in toks[2:]:
                    self._run([(line_no, cmd_text)])
            return i + 1

        # Block form: scan branches at level 1 up to the matching endif.
        end = self._find_matching(lines, i, "if")
        branches: List[Tuple[Optional[str], int, int]] = []  # (cond, start, stop)
        cond: Optional[str] = toks[1]
        seg_start = i + 1
        level = 0
        j = i + 1
        while j < end:
            s = lines[j][1].strip()
            head = s.split(None, 1)[0] if s.split() else ""
            if self._is_block_open(s, "if", self.macros) or head == "do":
                level += 1
            elif head in ("endif", "enddo") and level > 0:
                level -= 1
            elif level == 0 and head == "else":
                branches.append((cond, seg_start, j))
                cond = None                      # else branch
                seg_start = j + 1
            elif level == 0 and head == "elseif":
                branches.append((cond, seg_start, j))
                rest = s.split(None, 1)
                cond = self._tokenize(rest[1])[0] if len(rest) > 1 else "0"
                seg_start = j + 1
            j += 1
        branches.append((cond, seg_start, end))

        for branch_cond, start, stop in branches:
            taken = (branch_cond is None or
                     self._truthy(branch_cond, line_no))
            if taken:
                self._run(lines[start:stop])
                break
        return end + 1

    # -- define / macro expansion ---------------------------------------------
    def _exec_define(self, stripped: str, line_no: int) -> None:
        toks = self._tokenize(stripped)
        if len(toks) < 2:
            self.warnings.append(f"invalid define at line {line_no}")
            return
        name = toks[1]
        # $param expanded now; $$param survives (as $param) until invocation.
        body = [self.expand(t) for t in toks[2:] if t]
        self.macros[name] = body

    def _exec_macro(self, stripped: str, line_no: int) -> None:
        toks = self._tokenize(stripped)
        name = toks[0]
        args = toks[1:]
        self._n_expansions += 1

        def repl(match: "re.Match[str]") -> str:
            key = match.group(1)
            if key == "#":
                return str(self._n_expansions)
            idx = int(key)
            if idx == 0:
                return name
            return args[idx - 1] if idx - 1 < len(args) else ""

        expanded: List[Tuple[int, str]] = []
        for body_line in self.macros[name]:
            expanded.append((line_no, _MACRO_REF.sub(repl, body_line)))
        self._run(expanded)

    # -- include ---------------------------------------------------------------
    def _exec_include(self, stripped: str, line_no: int) -> None:
        toks = self._tokenize(stripped)
        if len(toks) < 2:
            self.warnings.append(f"invalid include at line {line_no}")
            return
        fname = self.expand(toks[1])
        path = fname if os.path.isabs(fname) else os.path.join(self.base_dir, fname)
        if not os.path.exists(path):
            self.warnings.append(
                f"include file '{fname}' (line {line_no}) not found; skipped")
            return
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            self._run(self._logical_lines(fh.read()))

    # -- ordinary commands ------------------------------------------------------
    def _parse_line(self, line: str, line_no: int) -> Optional[G4BLCommand]:
        tokens = self._tokenize(line)
        if not tokens:
            return None

        name = tokens[0]
        args: List[str] = []
        params: Dict[str, str] = {}
        order: List[str] = []
        for tok in tokens[1:]:
            if "=" in tok and not tok.startswith("="):
                key, _, value = tok.partition("=")
                params[key] = self.expand(value)
                if key not in order:
                    order.append(key)
            else:
                args.append(self.expand(tok))
        return G4BLCommand(
            name=name,
            args=args,
            params=params,
            param_order=order,
            line_no=line_no,
            raw=line,
        )

    def expand(self, value: str) -> str:
        """Substitute ``$name``/``${name}`` params; ``$$name`` -> ``$name``."""
        protected = value.replace("$$", _DOLLAR_SENTINEL)

        def repl(match: "re.Match[str]") -> str:
            key = match.group(1)
            return self.params.get(key, match.group(0))

        prev = None
        cur = protected
        for _ in range(10):
            if cur == prev:
                break
            prev = cur
            cur = _PARAM_REF.sub(repl, cur)
        return cur.replace(_DOLLAR_SENTINEL, "$")

    def _handle_param(self, cmd: G4BLCommand) -> None:
        # ``param -unset name=value`` / ``param -default name=value`` only set
        # the value if not already present.  Plain ``param`` always sets it.
        only_if_unset = any(a in ("-unset", "-default") for a in cmd.args)
        for key in cmd.param_order:
            if only_if_unset and key in self.params:
                continue
            value = cmd.params[key]
            # G4beamline evaluates a param value when it is a numeric
            # expression with an operator (guide sec. 5.1).
            if _has_operator(value):
                result = _eval_expr(value)
                if result is not None:
                    value = _format_eval(result)
            self.params[key] = value


def parse_g4bl(text: str, base_dir: str = "") -> Tuple[List[G4BLCommand], Dict[str, str]]:
    """Parse G4beamline *text*; return (commands, resolved-params)."""
    parser = G4BLParser(base_dir=base_dir)
    commands = parser.parse(text)
    return commands, parser.params


def parse_g4bl_file(path: str) -> Tuple[List[G4BLCommand], Dict[str, str]]:
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        return parse_g4bl(fh.read(), base_dir=os.path.dirname(os.path.abspath(path)))
