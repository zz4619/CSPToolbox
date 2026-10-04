"""Space-group symmetry: spglib detection, symmetry-operation parsing and
formatting, SHELX LATT/SYMM records and their expansion to an explicit cell."""

from __future__ import annotations

from fractions import Fraction
import math
import re
from typing import TYPE_CHECKING, Iterable

import numpy as np
from ase.cell import Cell
from ase.data import atomic_numbers
import spglib

from .lattice import canonicalize_fractional
from .records import AtomRecord

if TYPE_CHECKING:
    from .structure import CrystalStructure


_CSORM_SUPPORTED_TRANSLATIONS = frozenset(
    {
        Fraction(0, 1),
        Fraction(1, 4),
        Fraction(1, 3),
        Fraction(1, 2),
        Fraction(2, 3),
        Fraction(3, 4),
    }
)
_CSORM_DECIMAL_TRANSLATION_MAP = {
    "2": Fraction(1, 4),
    "3": Fraction(1, 3),
    "5": Fraction(1, 2),
    "6": Fraction(2, 3),
    "7": Fraction(3, 4),
}
_RES_TRANSLATION_ROUNDING_TOLERANCE = 0.01
_RES_TRANSLATION_TARGETS = (
    0.0,
    1.0 / 4.0,
    1.0 / 3.0,
    1.0 / 2.0,
    2.0 / 3.0,
    3.0 / 4.0,
    1.0,
)


def apply_symmetry_operation(
    operation: str,
    frac_coords: tuple[float, float, float],
) -> tuple[float, float, float]:
    """Apply an 'x,y,z'-style operation to fractional coordinates."""

    x, y, z = frac_coords
    components = [part.strip().replace(" ", "").lower() for part in operation.split(",")]
    values = [_evaluate_symmetry_component(component, {"x": x, "y": y, "z": z}) for component in components]
    return tuple(float(value) for value in values)  # type: ignore[return-value]


def expand_shelx_atoms(
    *,
    atoms: list[AtomRecord],
    cell: np.ndarray,
    latt_value: int,
    symmetry_operations: Iterable[str],
) -> list[AtomRecord]:
    """Generate every symmetry copy of the atoms from SHELX LATT/SYMM records."""

    generated_atoms: list[AtomRecord] = []
    centrosymmetric = latt_value > 0
    centering_translations = _shelx_centering_translations(latt_value)
    normalized_operations = [str(operation).strip() for operation in symmetry_operations]

    for atom in atoms:
        frac = tuple(float(value) for value in Cell(cell).scaled_positions(np.array([atom.coordinates], dtype=float))[0])
        generated_fractionals = [frac]
        generated_fractionals.extend(
            apply_symmetry_operation(operation, frac)
            for operation in normalized_operations
        )
        if centrosymmetric:
            generated_fractionals.extend(
                tuple(-value for value in coords)
                for coords in list(generated_fractionals)
            )

        unique_positions: list[tuple[float, float, float]] = []
        seen_positions: set[tuple[float, float, float]] = set()
        for coords in generated_fractionals:
            for translation in centering_translations:
                canonical = canonicalize_fractional(
                    tuple(float(left + right) for left, right in zip(coords, translation))
                )
                if canonical in seen_positions:
                    continue
                seen_positions.add(canonical)
                unique_positions.append(canonical)

        use_suffix = len(unique_positions) > 1
        for index, canonical in enumerate(unique_positions, start=1):
            cart = np.dot(np.array(canonical, dtype=float), cell)
            generated_atoms.append(
                AtomRecord(
                    label=f"{atom.label}_{index}" if use_suffix else atom.label,
                    element=atom.element,
                    coordinates=tuple(float(value) for value in cart),
                )
            )

    return generated_atoms


def normalize_space_group_symbol(symbol: str) -> str:
    """Normalise spacing in a Hermann-Mauguin symbol ('P21/c' -> 'P 21/c')."""

    normalized = " ".join(str(symbol).strip().split())
    match = re.match(r"^([A-Z])([^\s].*)$", normalized)
    if match:
        return f"{match.group(1)} {match.group(2)}"
    return normalized


def shelx_symmetry_records(
    space_group_symbol: str,
    *,
    hall_number: int | None = None,
    rounding: bool = False,
) -> tuple[int, list[str]]:
    """LATT value and SYMM operations for a space group (Hall number or symbol)."""

    if hall_number is None:
        hall_number = _space_group_hall_number(space_group_symbol)
    symmetry = spglib.get_symmetry_from_database(hall_number)
    return shelx_symmetry_records_from_symmetry(
        space_group_symbol=space_group_symbol,
        rotations=np.asarray(symmetry["rotations"], dtype=int),
        translations=np.asarray(symmetry["translations"], dtype=float),
        rounding=rounding,
    )


def shelx_symmetry_records_from_symmetry(
    *,
    space_group_symbol: str,
    rotations: np.ndarray,
    translations: np.ndarray,
    rounding: bool = False,
) -> tuple[int, list[str]]:
    """LATT value and SYMM operations from explicit rotations and translations."""

    rotations = np.asarray(rotations, dtype=int)
    translations = np.asarray(translations, dtype=float)

    inversion_index = None
    for index, rotation in enumerate(rotations):
        if np.array_equal(rotation, -np.eye(3, dtype=int)):
            inversion_index = index
            break

    centrosymmetric = inversion_index is not None
    latt_value = _shelx_latt_value(
        space_group_symbol=space_group_symbol,
        centrosymmetric=centrosymmetric,
    )
    centering_translations = _shelx_centering_translations(latt_value)

    operations: list[tuple[np.ndarray, np.ndarray]] = []
    seen_keys: set[tuple[tuple[int, ...], tuple[int, int, int]]] = set()
    inversion_translation = (
        translations[inversion_index] if inversion_index is not None else None
    )

    for rotation, translation in zip(rotations, translations):
        if np.array_equal(rotation, np.eye(3, dtype=int)):
            continue
        if centrosymmetric and np.array_equal(rotation, -np.eye(3, dtype=int)):
            continue

        key = _centering_reduced_symmetry_operation_key(
            rotation,
            translation,
            centering_translations=centering_translations,
        )
        if key in seen_keys:
            continue
        seen_keys.add(key)

        if centrosymmetric and inversion_translation is not None:
            inversion_mate_key = _centering_reduced_symmetry_operation_key(
                -rotation,
                translation + inversion_translation,
                centering_translations=centering_translations,
            )
            seen_keys.add(inversion_mate_key)

        operations.append((rotation, translation))

    operation_strings = [
        _format_shelx_symmetry_operation(rotation, translation, rounding=rounding)
        for rotation, translation in operations
    ]
    return latt_value, operation_strings


def deduplicate_shelx_symmetry_operations(
    symmetry_operations: Iterable[str],
    *,
    latt_value: int,
    rounding: bool = False,
) -> list[str]:
    """Drop SYMM operations that repeat one another up to centring."""

    centering_translations = _shelx_centering_translations(latt_value)
    deduplicated: list[str] = []
    seen_keys: set[tuple[tuple[int, ...], tuple[int, int, int]]] = set()

    for operation in symmetry_operations:
        normalized = _normalize_shelx_symmetry_operation(operation, rounding=rounding)
        rotation_rows: list[tuple[int, int, int]] = []
        translation_terms: list[float] = []
        for component in normalized.split(","):
            rotation_row, translation = _parse_general_symmetry_component(component)
            rotation_rows.append(rotation_row)
            translation_terms.append(translation)
        key = _centering_reduced_symmetry_operation_key(
            np.asarray(rotation_rows, dtype=int),
            np.asarray(translation_terms, dtype=float),
            centering_translations=centering_translations,
        )
        if key in seen_keys:
            continue
        seen_keys.add(key)
        deduplicated.append(normalized if rounding else str(operation).strip())

    return deduplicated


def normalize_csorm_symmetry_operation(operation: str) -> str:
    """Rewrite one operation as CSO-RM parses it; ValueError if CSO-RM cannot."""

    components = [component.strip() for component in operation.split(",")]
    if len(components) != 3:
        raise ValueError(f"Expected 3 symmetry components, got {len(components)} in {operation!r}")

    rotations: list[tuple[int, int, int]] = []
    translations: list[float] = []
    for component in components:
        rotation_row, translation = _parse_csorm_symmetry_component(component)
        rotations.append(rotation_row)
        translations.append(float(translation))

    return _format_shelx_symmetry_operation(
        np.asarray(rotations, dtype=int),
        np.asarray(translations, dtype=float),
        rounding=True,
    )


def spglib_cell(structure: CrystalStructure) -> tuple[np.ndarray, np.ndarray, list[int]]:
    """Return the (lattice, fractional positions, atomic numbers) cell spglib expects."""

    atoms = structure.to_ase_atoms()
    return (
        np.asarray(atoms.cell, dtype=float),
        np.asarray(atoms.get_scaled_positions(wrap=True), dtype=float),
        [int(atomic_numbers[atom.element]) for atom in structure.atoms],
    )


def symmetry_dataset(
    structure: CrystalStructure,
    *,
    symprec: float,
    angle_tolerance: float,
) -> spglib.SpglibDataset:
    """Run spglib symmetry detection; ValueError if no dataset is found."""

    dataset = spglib.get_symmetry_dataset(
        spglib_cell(structure),
        symprec=symprec,
        angle_tolerance=angle_tolerance,
    )
    if dataset is None:
        raise ValueError("Could not determine space-group symmetry for this structure.")
    return dataset


def _normalize_space_group_lookup_key(symbol: str) -> str:
    return normalize_space_group_symbol(symbol).replace(" ", "").replace("_", "").upper()


def _space_group_hall_number(space_group_symbol: str) -> int:
    target = _normalize_space_group_lookup_key(space_group_symbol)
    for hall_number in range(1, 531):
        space_group_type = spglib.get_spacegroup_type(hall_number)
        candidates: set[str] = set()
        for raw in (
            space_group_type.international_short,
            space_group_type.international_full,
            space_group_type.international,
        ):
            for piece in str(raw).split("="):
                candidates.add(_normalize_space_group_lookup_key(piece))
        if target in candidates:
            return hall_number
    raise ValueError(f"Unsupported or unknown space group for SHELX output: {space_group_symbol!r}")


def _shelx_centering_translations(latt_value: int) -> tuple[tuple[float, float, float], ...]:
    centering = abs(latt_value)
    centering_map: dict[int, tuple[tuple[float, float, float], ...]] = {
        1: ((0.0, 0.0, 0.0),),
        2: ((0.0, 0.0, 0.0), (0.5, 0.5, 0.5)),
        3: ((0.0, 0.0, 0.0), (2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0), (1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0)),
        4: ((0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)),
        5: ((0.0, 0.0, 0.0), (0.0, 0.5, 0.5)),
        6: ((0.0, 0.0, 0.0), (0.5, 0.0, 0.5)),
        7: ((0.0, 0.0, 0.0), (0.5, 0.5, 0.0)),
    }
    if centering not in centering_map:
        raise ValueError(f"Unsupported SHELX LATT value: {latt_value}")
    return centering_map[centering]


def _shelx_latt_value(space_group_symbol: str, *, centrosymmetric: bool) -> int:
    leading = normalize_space_group_symbol(space_group_symbol).split()[0].upper()
    centering_map = {
        "P": 1,
        "I": 2,
        "R": 3,
        "F": 4,
        "A": 5,
        "B": 6,
        "C": 7,
    }
    if leading not in centering_map:
        raise ValueError(f"Unsupported lattice centering for SHELX output: {space_group_symbol!r}")
    value = centering_map[leading]
    return value if centrosymmetric else -value


def _symmetry_operation_key(
    rotation: np.ndarray,
    translation: np.ndarray | Iterable[float],
) -> tuple[tuple[int, ...], tuple[int, int, int]]:
    normalized_translation = tuple(
        int(round(value * 24)) % 24
        for value in np.mod(np.asarray(list(translation), dtype=float), 1.0)
    )
    return tuple(rotation.reshape(-1).tolist()), normalized_translation


def _centering_reduced_symmetry_operation_key(
    rotation: np.ndarray,
    translation: np.ndarray | Iterable[float],
    *,
    centering_translations: Iterable[tuple[float, float, float]],
) -> tuple[tuple[int, ...], tuple[int, int, int]]:
    translation_array = np.asarray(list(translation), dtype=float)
    candidate_keys = [
        _symmetry_operation_key(rotation, translation_array - np.asarray(centering, dtype=float))
        for centering in centering_translations
    ]
    return min(candidate_keys)


def _normalize_shelx_symmetry_operation(operation: str, *, rounding: bool = False) -> str:
    components = [component.strip() for component in operation.split(",")]
    if len(components) != 3:
        raise ValueError(f"Expected 3 symmetry components, got {len(components)} in {operation!r}")

    rotations: list[tuple[int, int, int]] = []
    translations: list[float] = []
    for component in components:
        rotation_row, translation = _parse_general_symmetry_component(component)
        rotations.append(rotation_row)
        translations.append(float(translation))

    return _format_shelx_symmetry_operation(
        np.asarray(rotations, dtype=int),
        np.asarray(translations, dtype=float),
        rounding=rounding,
    )


def _parse_general_symmetry_component(component: str) -> tuple[tuple[int, int, int], float]:
    expr = component.replace(" ", "").upper()
    if not expr:
        raise ValueError("Empty symmetry component.")

    coeffs = {"X": 0, "Y": 0, "Z": 0}
    translation = 0.0
    expr = expr.replace("-", "+-")
    if expr.startswith("+-"):
        expr = "-" + expr[2:]
    terms = [term for term in expr.split("+") if term]
    if not terms:
        raise ValueError(f"Could not parse symmetry component {component!r}")

    for term in terms:
        if term in ("X", "Y", "Z"):
            coeffs[term] += 1
            continue
        if term in ("-X", "-Y", "-Z"):
            coeffs[term[1:]] -= 1
            continue
        translation += _parse_general_translation_token(term)

    rotation = (coeffs["X"], coeffs["Y"], coeffs["Z"])
    if any(abs(value) > 1 for value in rotation):
        raise ValueError(f"Unsupported symmetry rotation term in {component!r}")

    return rotation, translation


def _parse_general_translation_token(token: str) -> float:
    cleaned = token.strip().upper()
    sign = 1.0
    if cleaned.startswith("-"):
        sign = -1.0
        cleaned = cleaned[1:]
    elif cleaned.startswith("+"):
        cleaned = cleaned[1:]

    if not cleaned:
        raise ValueError("Empty translation token.")

    if "/" in cleaned:
        return sign * float(Fraction(cleaned))
    return sign * float(cleaned)


def _parse_csorm_symmetry_component(component: str) -> tuple[tuple[int, int, int], Fraction]:
    expr = component.replace(" ", "").upper()
    if not expr:
        raise ValueError("Empty symmetry component.")

    coeffs = {"X": 0, "Y": 0, "Z": 0}
    translation = Fraction(0, 1)
    expr = expr.replace("-", "+-")
    if expr.startswith("+-"):
        expr = "-" + expr[2:]
    terms = [term for term in expr.split("+") if term]
    if not terms:
        raise ValueError(f"Could not parse symmetry component {component!r}")

    for term in terms:
        if term in ("X", "Y", "Z"):
            coeffs[term] += 1
            continue
        if term in ("-X", "-Y", "-Z"):
            coeffs[term[1:]] -= 1
            continue
        translation += _parse_csorm_translation_token(term)

    normalized_translation = translation % 1
    if normalized_translation not in _CSORM_SUPPORTED_TRANSLATIONS:
        raise ValueError(
            f"Unsupported CSORM translation {component!r}: final offset {normalized_translation}"
        )

    rotation = (coeffs["X"], coeffs["Y"], coeffs["Z"])
    if any(abs(value) > 1 for value in rotation):
        raise ValueError(f"Unsupported CSORM rotation term in {component!r}")

    return rotation, normalized_translation


def _parse_csorm_translation_token(token: str) -> Fraction:
    cleaned = token.strip().upper()
    sign = 1
    if cleaned.startswith("-"):
        sign = -1
        cleaned = cleaned[1:]
    elif cleaned.startswith("+"):
        cleaned = cleaned[1:]

    if not cleaned:
        raise ValueError("Empty translation token.")

    if "/" in cleaned:
        fraction = Fraction(cleaned)
        normalized = fraction % 1
        if normalized not in _CSORM_SUPPORTED_TRANSLATIONS:
            raise ValueError(f"Unsupported CSORM rational translation {token!r}")
        return sign * fraction

    if "." in cleaned:
        dot_index = cleaned.find(".")
        if dot_index == len(cleaned) - 1:
            raise ValueError(f"Invalid CSORM decimal translation {token!r}")
        decimal_digit = cleaned[dot_index + 1]
        if decimal_digit not in _CSORM_DECIMAL_TRANSLATION_MAP:
            raise ValueError(f"Unsupported CSORM decimal translation {token!r}")
        integer_part = cleaned[:dot_index]
        whole = int(integer_part) if integer_part else 0
        return sign * (Fraction(whole, 1) + _CSORM_DECIMAL_TRANSLATION_MAP[decimal_digit])

    return sign * Fraction(cleaned)


def _format_shelx_symmetry_operation(
    rotation: np.ndarray,
    translation: np.ndarray,
    *,
    rounding: bool = False,
) -> str:
    axes = ("X", "Y", "Z")
    components: list[str] = []
    for row_index in range(3):
        terms: list[str] = []
        offset = _format_fractional_offset(float(translation[row_index]), rounding=rounding)
        if offset:
            terms.append(offset)

        for coeff, axis in zip(rotation[row_index], axes):
            if coeff == 0:
                continue
            if coeff == 1:
                if terms:
                    terms.append(f"+{axis}")
                else:
                    terms.append(axis)
            elif coeff == -1:
                if terms:
                    terms.append(f"-{axis}")
                else:
                    terms.append(f"-{axis}")
            else:
                raise ValueError("Only unit rotation coefficients are supported for SHELX symmetry output.")

        component = "".join(terms) if terms else "0"
        components.append(component)
    return ",".join(components)


def _format_fractional_offset(value: float, *, rounding: bool = False) -> str:
    normalized = float(np.mod(value, 1.0))
    if rounding:
        normalized = _round_res_fractional_translation(normalized)
    if abs(normalized) < 1e-8 or abs(normalized - 1.0) < 1e-8:
        return ""
    fraction = Fraction(normalized).limit_denominator(24)
    if abs(float(fraction) - normalized) > 1e-6:
        return f"{normalized:.6f}".rstrip("0").rstrip(".")
    if fraction.denominator == 1:
        return str(fraction.numerator)
    return f"{fraction.numerator}/{fraction.denominator}"


def _round_res_fractional_translation(value: float) -> float:
    normalized = float(np.mod(value, 1.0))
    for target in _RES_TRANSLATION_TARGETS:
        if abs(normalized - target) <= _RES_TRANSLATION_ROUNDING_TOLERANCE:
            return 0.0 if math.isclose(target, 1.0, abs_tol=1e-12) else target
    return normalized


def _evaluate_symmetry_component(component: str, variables: dict[str, float]) -> float:
    """Evaluate one operation component such as ``-x+1/2`` or ``x-y``.

    Terms are summed left to right with the same floating-point operations
    Python arithmetic would use, so results are bit-identical to evaluating the
    text as an expression. Only x/y/z, integers, decimals and fractions are
    accepted; anything else raises ValueError (operations come from files and
    must never be executed as code).
    """

    if not component:
        raise ValueError("Empty symmetry-operation component.")
    total: float | None = None
    index = 0
    while index < len(component):
        sign = 1.0
        while index < len(component) and component[index] in "+-":
            if component[index] == "-":
                sign = -sign
            index += 1
        end = index
        while end < len(component) and component[end] not in "+-":
            end += 1
        operand = component[index:end]
        if not operand:
            raise ValueError(f"Missing operand in symmetry component {component!r}.")
        if operand in variables:
            value = variables[operand]
        elif "/" in operand:
            numerator, _, denominator = operand.partition("/")
            try:
                value = float(numerator) / float(denominator)
            except ValueError as error:
                raise ValueError(f"Invalid term {operand!r} in symmetry component {component!r}.") from error
        else:
            try:
                value = float(operand)
            except ValueError as error:
                raise ValueError(f"Invalid term {operand!r} in symmetry component {component!r}.") from error
        term = -value if sign < 0 else value
        total = term if total is None else total + term
        index = end
    return total
