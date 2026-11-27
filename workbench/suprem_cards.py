"""Read the command tree used by the local SUPREM4 runtime."""

from dataclasses import dataclass, field
from io import BytesIO
import math
import os
from pathlib import Path
import re
import struct

from runner import ROOT, FACTOR, NUMBER

REAL, INTEGER, STRING, BOOLEAN, COMMAND, CHOICE = (0x40, 0x80, 0x100, 0x200, 0x400, 0x1000)
PARAMETER_MASK = 0x083F
VALIDATED_KEY = ROOT / "suprem4/build/validated-native64/data/suprem.uk"
SOURCE_KEY = ROOT / "suprem4/data/suprem.uk"


@dataclass
class Parameter:
    name: str
    kind: int
    description: str
    children: list["Parameter"] = field(default_factory=list)

    def __getitem__(self, name):
        return next(child for child in self.children if child.name == name)


def key_file():
    if os.environ.get("TCAD_SUPREM4"):
        path = os.environ.get("KEYFILE") or os.environ.get("SUP4KEYFILE")
        if path:
            return Path(path)
    return VALIDATED_KEY if VALIDATED_KEY.is_file() else SOURCE_KEY


def catalog(path):
    # Layout follows the native reader in src/misc/reader.c: depth, 12-byte
    # name, type, default, two strings, and a recursive validity expression.
    stream = BytesIO(path.read_bytes())

    def string():
        chars = bytearray()
        while (char := stream.read(1)) != b"\0":
            if not char:
                raise ValueError("Unterminated string in the SUPREM4 command catalog.")
            chars.extend(char)
        return chars.decode("ascii", errors="replace")

    def integer():
        return struct.unpack("=i", stream.read(4))[0]

    def skip_expression():
        kind = integer()
        if kind == 1:  # Operator with one operand for !, two for the others
            op = integer()
            skip_expression()
            if op != ord("!"):
                skip_expression()
        elif kind == REAL:
            integer()  # A four-byte float
        elif kind == 3:  # Parameter reference
            string()
        elif kind != 0:
            raise ValueError("Unknown validity expression in the SUPREM4 command catalog.")

    cards = []
    stack = [cards]
    try:
        while depth_bytes := stream.read(4):
            depth = struct.unpack("=i", depth_bytes)[0]
            if not 0 <= depth < len(stack):
                raise ValueError("Invalid nesting in the SUPREM4 command catalog.")
            name = stream.read(12).partition(b"\0")[0].decode("ascii").lower()
            kind = integer() & ~PARAMETER_MASK
            if kind == STRING:
                string()
            elif kind in {REAL, INTEGER, BOOLEAN, COMMAND, CHOICE}:
                integer()
            else:
                raise ValueError("Unknown parameter type in the SUPREM4 command catalog.")
            unit = string()
            message = string()
            skip_expression()
            node = Parameter(name, kind, unit or message)
            stack[depth].append(node)
            stack[depth + 1:] = [node.children]
    except (struct.error, RecursionError, UnicodeError) as error:
        raise ValueError(f"Could not read SUPREM4 command catalog {path}") from error
    return {card.name: card for card in cards}


def parameter_token(parameter, value):
    value = value.strip()
    if not value:
        return ""
    name = parameter.name
    if parameter.kind == BOOLEAN:
        return name if value == "On" else f"{name}=false"
    if parameter.kind == INTEGER and not re.fullmatch(r"-?\d+", value):
        raise ValueError(f"{name} requires an integer.")
    if "\n" in value or "\r" in value or "\x00" in value:
        raise ValueError(f"{name} requires a single-line value.")
    if parameter.kind == STRING:
        for quote in ('"', "'"):
            if quote not in value:
                return f"{name}={quote}{value}{quote}"
        raise ValueError(f"{name} contains both quote styles. Edit that value in the deck.")
    if parameter.kind == REAL and len(value.split()) > 1:
        value = f"({value})"
    return f"{name}={value}"


def polygon_cards(material, points, spacing):
    rows = [re.split(r"[\s,]+", line.strip()) for line in points.splitlines() if line.strip()]
    if len(rows) < 3:
        raise ValueError("Enter at least three polygon vertices.")
    if any(len(row) != 2 or any(numeric_literal(value) is None for value in row) for row in rows):
        raise ValueError("Enter one numeric X, Y pair per line.")
    if spacing and (numeric_literal(spacing) is None or numeric_literal(spacing) <= 0):
        raise ValueError("Polygon spacing must be a positive number.")
    return "\n".join(
        " ".join(filter(None, (
            "etch", material if index == 0 else "",
            "start" if index == 0 else "done" if index == len(rows) - 1 else "continue",
            f"x={x}", f"y={y}", f"spacing={spacing}" if spacing else "",
        )))
        for index, (x, y) in enumerate(rows)
    )


def numeric_literal(value):
    factor = FACTOR.fullmatch(value)
    literal = factor[2] if factor else value
    if re.fullmatch(NUMBER, literal):
        number = float(literal)
        if not math.isfinite(number):
            raise ValueError("Numeric values must be finite.")
        return number
    return None
