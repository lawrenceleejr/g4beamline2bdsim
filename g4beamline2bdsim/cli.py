"""Command-line interface for g4beamline2bdsim."""

from __future__ import annotations

import argparse
import os
import sys
from typing import List, Optional

from . import __version__
from .converter import Converter
from .gmad_writer import GmadWriter
from .parser import parse_g4bl


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="g4beamline2bdsim",
        description="Convert a G4beamline input file into a BDSIM (GMAD) model.",
    )
    p.add_argument("input", help="G4beamline input file (use '-' for stdin)")
    p.add_argument(
        "-o",
        "--output",
        help="Output GMAD file (default: <input>.gmad, or stdout for '-')",
    )
    p.add_argument(
        "--line-name",
        default="mainline",
        help="Name of the generated beamline (default: mainline)",
    )
    p.add_argument(
        "--sample-all",
        action="store_true",
        help="Emit 'sample, all;' instead of per-detector samplers",
    )
    p.add_argument(
        "--no-gdml",
        action="store_true",
        help="Convert box/tubs material volumes to drifts instead of exporting "
             "GDML geometry (material will not interact with the beam)",
    )
    p.add_argument(
        "--no-field-maps",
        action="store_true",
        help="Do not generate BDSIM field maps for fieldexpr formulas",
    )
    p.add_argument(
        "--no-solenoid-field-map",
        action="store_true",
        help="Use a native BDSIM solenoid with ks (hard-edge) instead of the "
             "default full 3D G4beamline coil field map (more robust for very "
             "strong solenoids)",
    )
    p.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress conversion warnings on stderr",
    )
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


def run(argv: Optional[List[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)

    # Read input.
    if args.input == "-":
        text = sys.stdin.read()
        source_name = "<stdin>"
    else:
        if not os.path.exists(args.input):
            print(f"error: input file not found: {args.input}", file=sys.stderr)
            return 2
        with open(args.input, "r", encoding="utf-8", errors="replace") as fh:
            text = fh.read()
        source_name = os.path.basename(args.input)

    base_dir = "" if args.input == "-" else os.path.dirname(os.path.abspath(args.input))
    commands, _params = parse_g4bl(text, base_dir=base_dir)
    converter = Converter(commands, source_name=source_name,
                          emit_gdml=not args.no_gdml,
                          emit_field_maps=not args.no_field_maps,
                          solenoid_field_map=not args.no_solenoid_field_map,
                          base_dir=base_dir)
    model = converter.convert()
    model.line_name = args.line_name
    if args.sample_all:
        model.samplers = ["all"]

    writer = GmadWriter(model)
    gmad_text = writer.to_string()

    # Determine output path.
    if args.output == "-" or (args.output is None and args.input == "-"):
        sys.stdout.write(gmad_text)
        if model.aux_files and not args.quiet:
            print(
                f"note: {len(model.aux_files)} auxiliary file(s) (GDML/field "
                f"maps) not written to stdout; use -o to emit them: "
                f"{', '.join(model.aux_files)}",
                file=sys.stderr,
            )
    else:
        out_path = args.output
        if out_path is None:
            base, _ = os.path.splitext(args.input)
            out_path = base + ".gmad"
        with open(out_path, "w", encoding="utf-8") as fh:
            fh.write(gmad_text)
        # Write auxiliary files (GDML geometry, field maps) next to the GMAD.
        model.write_aux_files(os.path.dirname(os.path.abspath(out_path)))
        if not args.quiet:
            extra = (f" (+ {len(model.aux_files)} auxiliary file(s))"
                     if model.aux_files else "")
            print(f"Wrote {out_path}{extra}", file=sys.stderr)

    # Report warnings.
    if model.warnings and not args.quiet:
        print(
            f"\n{len(model.warnings)} conversion warning(s):", file=sys.stderr
        )
        for w in model.warnings:
            print(f"  - {w}", file=sys.stderr)

    return 0


def main() -> None:  # pragma: no cover
    sys.exit(run())


if __name__ == "__main__":  # pragma: no cover
    main()
