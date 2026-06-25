"""Parser for CSPToolbox numeric Z-matrix files."""

from __future__ import annotations

from collections import Counter
import math
from pathlib import Path
from typing import Iterable

from .model import ZMatrixAtom, ZMatrixDocument


ELEMENT_SYMBOLS = frozenset(
    """
    H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni
    Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe
    Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg
    Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg
    Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og
    """.split()
)


def load_zmatrix(path: str | Path) -> ZMatrixDocument:
    """Load a CSPToolbox numeric Z-matrix file."""

    source_path = Path(path)
    return parse_zmatrix_text(
        source_path.read_text(encoding="utf-8"),
        source_name=str(source_path),
    )


def parse_zmatrix_text(text: str, *, source_name: str | None = None) -> ZMatrixDocument:
    """Parse a CSPToolbox ``# ZMAT v1`` style numeric Z-matrix.

    Current CSPToolbox files store element symbols and 1-based references. This
    parser also accepts a backward-compatible labelled row form:
    ``label element ...``.
    """

    title = Path(source_name).stem if source_name else "Z-matrix molecule"
    labels_from_comment: list[str] | None = None
    explicit_bonds: set[tuple[int, int]] = set()
    row_specs: list[tuple[int, list[str]]] = []

    for line_number, raw_line in enumerate(text.splitlines(), start=1):
        stripped = raw_line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            comment = stripped[1:].strip()
            lower = comment.lower()
            if lower.startswith("title:"):
                parsed_title = comment.split(":", 1)[1].strip()
                if parsed_title:
                    title = parsed_title
            elif lower.startswith(("labels:", "atom-labels:")):
                labels_from_comment = comment.split(":", 1)[1].split()
            elif lower.startswith("bonds:"):
                explicit_bonds.update(_parse_bonds_comment(comment, line_number))
            continue

        row_text = stripped.split("#", 1)[0].strip()
        if row_text:
            row_specs.append((line_number, row_text.split()))

    if not row_specs:
        raise ValueError("No Z-matrix atom rows found.")

    parsed_rows = [_parse_row(tokens, line_number) for line_number, tokens in row_specs]
    labels = _resolve_labels(parsed_rows, labels_from_comment)
    atoms = tuple(
        _build_atom(row_index, line_number, element, label, values)
        for row_index, ((line_number, element, _row_label, values), label) in enumerate(
            zip(parsed_rows, labels),
            start=1,
        )
    )

    for atom in atoms:
        _validate_atom(atom)

    return ZMatrixDocument(
        title=title,
        atoms=atoms,
        explicit_bonds=frozenset(_normal_pair(*pair) for pair in explicit_bonds),
        source_name=source_name,
    )


def _parse_row(
    tokens: list[str],
    line_number: int,
) -> tuple[int, str, str | None, list[str]]:
    if not tokens:
        raise ValueError(f"Line {line_number}: empty Z-matrix row.")

    row_label: str | None = None
    element: str
    values: list[str]

    if _is_element(tokens[0]) and (len(tokens) - 1) % 2 == 0:
        element = _normal_element(tokens[0])
        values = tokens[1:]
    elif len(tokens) >= 2 and _is_element(tokens[1]) and (len(tokens) - 2) % 2 == 0:
        row_label = tokens[0]
        element = _normal_element(tokens[1])
        values = tokens[2:]
    elif (
        len(tokens) >= 2
        and _is_element(tokens[0])
        and not _looks_int(tokens[1])
        and (len(tokens) - 2) % 2 == 0
    ):
        element = _normal_element(tokens[0])
        row_label = tokens[1]
        values = tokens[2:]
    else:
        raise ValueError(
            f"Line {line_number}: expected 'element ...' or 'label element ...'."
        )

    if len(values) not in (0, 2, 4, 6):
        raise ValueError(f"Line {line_number}: internal-coordinate fields must be pairs.")

    return line_number, element, row_label, values


def _resolve_labels(
    parsed_rows: Iterable[tuple[int, str, str | None, list[str]]],
    labels_from_comment: list[str] | None,
) -> list[str]:
    rows = list(parsed_rows)
    row_labels = [row_label for _line, _element, row_label, _values in rows]

    if labels_from_comment is not None and len(labels_from_comment) != len(rows):
        raise ValueError(
            "# labels count does not match the number of Z-matrix atom rows "
            f"({len(labels_from_comment)} labels for {len(rows)} rows)."
        )

    element_counts: Counter[str] = Counter()
    labels: list[str] = []
    for index, (_line, element, row_label, _values) in enumerate(rows):
        element_counts[element] += 1
        generated = f"{element}{element_counts[element]}"
        label = row_label or (labels_from_comment[index] if labels_from_comment else generated)
        labels.append(label)

    if len(labels) != len(set(labels)):
        raise ValueError("Z-matrix atom labels are not unique.")
    if any(label is None for label in row_labels) and labels_from_comment is None:
        return labels
    return labels


def _build_atom(
    row_index: int,
    line_number: int,
    element: str,
    label: str,
    values: list[str],
) -> ZMatrixAtom:
    pairs = [_parse_pair(values[index : index + 2], line_number) for index in range(0, len(values), 2)]
    bond_to, bond_length = pairs[0] if len(pairs) >= 1 else (None, None)
    angle_to, angle_degrees = pairs[1] if len(pairs) >= 2 else (None, None)
    dihedral_to, dihedral_degrees = pairs[2] if len(pairs) >= 3 else (None, None)
    return ZMatrixAtom(
        row_index=row_index,
        line_number=line_number,
        label=label,
        element=element,
        bond_to=bond_to,
        bond_length=bond_length,
        angle_to=angle_to,
        angle_degrees=angle_degrees,
        dihedral_to=dihedral_to,
        dihedral_degrees=dihedral_degrees,
    )


def _parse_pair(tokens: list[str], line_number: int) -> tuple[int, float]:
    if len(tokens) != 2:
        raise ValueError(f"Line {line_number}: expected reference/value pair.")
    try:
        reference = int(tokens[0])
    except ValueError as error:
        raise ValueError(f"Line {line_number}: reference {tokens[0]!r} is not an integer.") from error
    try:
        value = float(tokens[1])
    except ValueError as error:
        raise ValueError(f"Line {line_number}: value {tokens[1]!r} is not numeric.") from error
    if not math.isfinite(value):
        raise ValueError(f"Line {line_number}: value {tokens[1]!r} is not finite.")
    return reference, value


def _validate_atom(atom: ZMatrixAtom) -> None:
    row = atom.row_index
    expected_pairs = 0 if row == 1 else 1 if row == 2 else 2 if row == 3 else 3
    observed_pairs = sum(
        value is not None for value in (atom.bond_to, atom.angle_to, atom.dihedral_to)
    )
    if observed_pairs != expected_pairs:
        raise ValueError(
            f"Line {atom.line_number}: row {row} has {observed_pairs} "
            f"reference pairs; expected {expected_pairs}."
        )

    for name, reference in (
        ("bond", atom.bond_to),
        ("angle", atom.angle_to),
        ("dihedral", atom.dihedral_to),
    ):
        if reference is None:
            continue
        if not 1 <= reference < row:
            raise ValueError(
                f"Line {atom.line_number}: {name} reference {reference} "
                f"does not point to an earlier row."
            )

    if atom.bond_length is not None and atom.bond_length <= 0.0:
        raise ValueError(f"Line {atom.line_number}: bond length must be positive.")


def _parse_bonds_comment(comment: str, line_number: int) -> set[tuple[int, int]]:
    payload = comment.split(":", 1)[1].strip()
    bonds: set[tuple[int, int]] = set()
    if not payload:
        return bonds
    for token in payload.split():
        if "-" not in token:
            raise ValueError(f"Line {line_number}: invalid bond token {token!r}.")
        left, right = token.split("-", 1)
        try:
            bonds.add(_normal_pair(int(left), int(right)))
        except ValueError as error:
            raise ValueError(f"Line {line_number}: invalid bond token {token!r}.") from error
    return bonds


def _normal_pair(left: int, right: int) -> tuple[int, int]:
    if left == right:
        raise ValueError("A bond cannot connect an atom to itself.")
    return (left, right) if left < right else (right, left)


def _is_element(value: str) -> bool:
    return _normal_element(value) in ELEMENT_SYMBOLS


def _normal_element(value: str) -> str:
    if not value:
        return value
    return value[0].upper() + value[1:].lower()


def _looks_int(value: str) -> bool:
    try:
        int(value)
    except ValueError:
        return False
    return True

