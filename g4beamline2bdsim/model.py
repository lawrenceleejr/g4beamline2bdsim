"""Intermediate representation of a BDSIM (GMAD) model.

The converter populates these structures from parsed G4beamline commands and the
:class:`~g4beamline2bdsim.gmad_writer.GmadWriter` serialises them to GMAD text.

Parameter values stored on an :class:`Element`/``beam``/``option`` may be:

* a number (int/float)  -> written verbatim (GMAD default units),
* a ``(value, unit)`` tuple -> written as ``value*unit`` (e.g. ``5*cm``),
* a string -> written double-quoted (e.g. material names, distrType).

Use the helper constructors in :mod:`converter` rather than building these by
hand.
"""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

# A parameter value is a number, a (number, unit-string) tuple, or a string.
ParamValue = Any


@dataclass
class Element:
    """A single GMAD component definition (``name: type, params...;``)."""

    name: str
    type: str
    params: "OrderedDict[str, ParamValue]" = field(default_factory=OrderedDict)
    comment: Optional[str] = None

    def set(self, key: str, value: ParamValue) -> "Element":
        self.params[key] = value
        return self


@dataclass
class BdsimModel:
    """Container for a full BDSIM model."""

    elements: List[Element] = field(default_factory=list)
    # Ordered list of element names making up the main beamline.
    line: List[str] = field(default_factory=list)
    line_name: str = "mainline"
    beam: "OrderedDict[str, ParamValue]" = field(default_factory=OrderedDict)
    options: "OrderedDict[str, ParamValue]" = field(default_factory=OrderedDict)
    samplers: List[str] = field(default_factory=list)  # element names, or ["all"]
    header_comments: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    # name -> Element index for de-duplication / lookup.
    _by_name: Dict[str, int] = field(default_factory=dict, repr=False)

    def add_element(self, element: Element) -> Element:
        if element.name in self._by_name:
            # Element already defined; keep the first definition.
            return self.elements[self._by_name[element.name]]
        self._by_name[element.name] = len(self.elements)
        self.elements.append(element)
        return element

    def get_element(self, name: str) -> Optional[Element]:
        idx = self._by_name.get(name)
        return self.elements[idx] if idx is not None else None

    def warn(self, message: str) -> None:
        self.warnings.append(message)
