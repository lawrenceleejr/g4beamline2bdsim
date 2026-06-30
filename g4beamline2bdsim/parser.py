"""Parser for G4beamline input files.

G4beamline input is a line-oriented command language.  Each (logical) line is a
command consisting of a command name followed by positional arguments and
``name=value`` keyword arguments::

    command arg1 arg2 name1=value1 name2=value2

Key syntactic features handled here:

* ``#`` introduces a comment that runs to the end of the line.
* A backslash ``\\`` at the end of a line continues the command on the next
  line.
* ``param name=value ...`` defines parameters that are later referenced with
  ``$name`` (or ``${name}``).  ``param -unset`` / ``param -default`` only set a
  value if it is not already defined.
* Whitespace separates tokens; values containing spaces may be quoted with
  single or double quotes.

The parser is intentionally permissive: anything it does not understand is kept
verbatim so the converter can decide how to handle (or warn about) it.
"""

from __future__ import annotations

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


def _format_eval(value: float) -> str:
    """Render an evaluated number compactly (integers without a trailing .0)."""
    if value == int(value) and abs(value) < 1e15:
        return str(int(value))
    return repr(value)


class G4BLParser:
    """Parse G4beamline text into a list of :class:`G4BLCommand` objects."""

    def __init__(self) -> None:
        self.params: Dict[str, str] = {}
        self.commands: List[G4BLCommand] = []

    # -- public API ---------------------------------------------------------
    def parse(self, text: str) -> List[G4BLCommand]:
        for line_no, logical_line in self._logical_lines(text):
            stripped = logical_line.strip()
            if not stripped:
                continue
            cmd = self._parse_line(stripped, line_no)
            if cmd is None:
                continue
            if cmd.name == "param":
                self._handle_param(cmd)
            self.commands.append(cmd)
        return self.commands

    # -- helpers ------------------------------------------------------------
    @staticmethod
    def _strip_comment(line: str) -> str:
        """Remove a ``#`` comment, respecting quotes.

        Per the G4beamline User's Guide, ``#`` starts a comment only at the start
        of a line (after optional whitespace) or when preceded by whitespace -
        a ``#`` embedded in a token (e.g. ``rename=Det#``) is *not* a comment.
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
        """Yield (line_number, logical_line) merging backslash continuations."""
        result: List[Tuple[int, str]] = []
        buffer = ""
        start_no = 0
        for idx, raw in enumerate(text.splitlines(), start=1):
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

    def _parse_line(self, line: str, line_no: int) -> Optional[G4BLCommand]:
        try:
            tokens = shlex.split(line, posix=True)
        except ValueError:
            # Unbalanced quotes etc. - fall back to a naive split.
            tokens = line.split()
        if not tokens:
            return None

        name = tokens[0]
        args: List[str] = []
        params: Dict[str, str] = {}
        order: List[str] = []
        for tok in tokens[1:]:
            if "=" in tok and not tok.startswith("="):
                key, _, value = tok.partition("=")
                # Expand parameter references in the value.
                value = self.expand(value)
                params[key] = value
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
        """Substitute ``$name`` / ``${name}`` parameter references."""

        def repl(match: "re.Match[str]") -> str:
            key = match.group(1)
            return self.params.get(key, match.group(0))

        # Repeat until stable (params may reference other params).
        prev = None
        cur = value
        for _ in range(10):
            if cur == prev:
                break
            prev = cur
            cur = _PARAM_REF.sub(repl, cur)
        return cur

    def _handle_param(self, cmd: G4BLCommand) -> None:
        # ``param -unset name=value`` / ``param -default name=value`` only set
        # the value if not already present.  Plain ``param name=value`` always
        # sets it.
        only_if_unset = any(a in ("-unset", "-default") for a in cmd.args)
        for key in cmd.param_order:
            if only_if_unset and key in self.params:
                continue
            value = cmd.params[key]
            # G4beamline evaluates a param value when it is a numeric expression
            # with an operator (sec. 5.1).  Store the evaluated number so later
            # $references resolve to a value, not an expression string.
            if _has_operator(value):
                result = _eval_expr(value)
                if result is not None:
                    value = _format_eval(result)
            self.params[key] = value


def parse_g4bl(text: str) -> Tuple[List[G4BLCommand], Dict[str, str]]:
    """Parse G4beamline *text*; return (commands, resolved-params)."""
    parser = G4BLParser()
    commands = parser.parse(text)
    return commands, parser.params


def parse_g4bl_file(path: str) -> Tuple[List[G4BLCommand], Dict[str, str]]:
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        return parse_g4bl(fh.read())
