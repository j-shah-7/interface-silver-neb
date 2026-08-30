#!/usr/bin/env python3
"""Add z-only VASP constraints for a metal / defective h-BN interface.

Designed for a structure containing one metal species, B, N, a layered metal
slab, and one isolated migrating metal atom. By default, the two metal layers
nearest h-BN are fixed in Cartesian z, the other two metal layers and migrating
metal are fully free, B/N atoms at an in-plane distance strictly below 3.5 A
from the migrating atom are fully free, and the remaining B/N atoms are fixed
in z. All x/y coordinates remain free.

This version uses only the Python standard library; NumPy is not required.
"""

from __future__ import annotations

import argparse
import math
import statistics
from pathlib import Path


def row_times_matrix(vector, matrix):
    return [
        sum(vector[row] * matrix[row][column] for row in range(3))
        for column in range(3)
    ]


def inverse_3x3(matrix):
    a, b, c = matrix[0]
    d, e, f = matrix[1]
    g, h, i = matrix[2]
    determinant = (
        a * (e * i - f * h)
        - b * (d * i - f * g)
        + c * (d * h - e * g)
    )
    if abs(determinant) < 1.0e-14:
        raise ValueError("The POSCAR lattice matrix is singular.")
    return [
        [(e * i - f * h) / determinant, (c * h - b * i) / determinant, (b * f - c * e) / determinant],
        [(f * g - d * i) / determinant, (a * i - c * g) / determinant, (c * d - a * f) / determinant],
        [(d * h - e * g) / determinant, (b * g - a * h) / determinant, (a * e - b * d) / determinant],
    ]


def fractional_to_cartesian(fractional, lattice):
    return [row_times_matrix(position, lattice) for position in fractional]


def read_poscar(path: Path):
    lines = path.read_text().splitlines()
    if len(lines) < 9:
        raise ValueError("The input is too short to be a valid POSCAR.")
    scale = float(lines[1].split()[0])
    if scale <= 0:
        raise ValueError("This script currently requires a positive POSCAR scale.")
    lattice = [
        [float(value) * scale for value in lines[index].split()[:3]]
        for index in range(2, 5)
    ]
    symbols = lines[5].split()
    counts = [int(value) for value in lines[6].split()]
    if len(symbols) != len(counts):
        raise ValueError("A VASP-5 POSCAR with element names is required.")

    line = 7
    has_selective_dynamics = lines[line].strip().lower().startswith("s")
    if has_selective_dynamics:
        line += 1
    coordinate_mode = lines[line].strip().lower()
    direct = coordinate_mode.startswith("d")
    cartesian = coordinate_mode.startswith(("c", "k"))
    if not direct and not cartesian:
        raise ValueError("Coordinate mode must be Direct or Cartesian.")
    start = line + 1
    number_of_atoms = sum(counts)
    coordinates = [
        [float(value) for value in lines[start + index].split()[:3]]
        for index in range(number_of_atoms)
    ]
    existing_flags = None
    if has_selective_dynamics:
        existing_flags = []
        for index in range(number_of_atoms):
            fields = lines[start + index].split()
            if len(fields) < 6:
                raise ValueError(
                    f"Selective-dynamics flags are missing for atom #{index + 1}."
                )
            atom_flags = [value.upper() for value in fields[3:6]]
            if any(value not in ("T", "F") for value in atom_flags):
                raise ValueError(
                    f"Invalid selective-dynamics flags for atom #{index + 1}."
                )
            existing_flags.append(atom_flags)
    if direct:
        fractional = [position[:] for position in coordinates]
    else:
        inverse_lattice = inverse_3x3(lattice)
        cartesian_coordinates = [
            [value * scale for value in position] for position in coordinates
        ]
        fractional = [
            row_times_matrix(position, inverse_lattice)
            for position in cartesian_coordinates
        ]

    elements = []
    for symbol, count in zip(symbols, counts):
        elements.extend([symbol] * count)
    return lines[0], lattice, symbols, counts, elements, fractional, existing_flags


def migrating_metal_index(elements, cartesian, metal_symbol, bn_plane, layer_tolerance):
    metal_ids = [index for index, element in enumerate(elements) if element == metal_symbol]
    z_coordinates = [cartesian[index][2] for index in metal_ids]
    neighbor_counts = [
        sum(abs(other_z - atom_z) <= layer_tolerance for other_z in z_coordinates)
        for atom_z in z_coordinates
    ]
    isolated = [atom for atom, neighbors in zip(metal_ids, neighbor_counts) if neighbors <= 3]
    if isolated:
        return min(isolated, key=lambda atom: abs(cartesian[atom][2] - bn_plane))
    raise ValueError(
        "Could not identify one isolated migrating metal atom. "
        "Try a smaller --layer-tolerance."
    )


def xy_distance(fractional_a, fractional_b, lattice):
    difference = [fractional_a[index] - fractional_b[index] for index in range(3)]
    difference[0] -= round(difference[0])
    difference[1] -= round(difference[1])
    difference[2] = 0.0
    cartesian_difference = row_times_matrix(difference, lattice)
    return math.sqrt(sum(value * value for value in cartesian_difference))


def split_framework_layers(cartesian, framework, number_of_layers):
    ordered = sorted(framework, key=lambda atom: cartesian[atom][2])
    gaps = [
        cartesian[ordered[index + 1]][2] - cartesian[ordered[index]][2]
        for index in range(len(ordered) - 1)
    ]
    if len(gaps) < number_of_layers - 1:
        raise ValueError("Not enough framework atoms to identify metal layers.")
    largest_gap_indices = sorted(
        range(len(gaps)), key=lambda index: gaps[index], reverse=True
    )[: number_of_layers - 1]
    split_positions = sorted(index + 1 for index in largest_gap_indices)
    layers = []
    beginning = 0
    for ending in split_positions + [len(ordered)]:
        layers.append(ordered[beginning:ending])
        beginning = ending
    layers.sort(key=lambda atoms: statistics.median(cartesian[atom][2] for atom in atoms))
    if len(layers) != number_of_layers or min(map(len, layers)) < 4:
        raise ValueError("Metal-layer detection failed; inspect the input structure.")
    return layers


def analyze_planes(elements, cartesian, metal_symbol, migrant, number_of_layers):
    framework = [
        index for index, element in enumerate(elements)
        if element == metal_symbol and index != migrant
    ]
    bn_ids = [index for index, element in enumerate(elements) if element in ("B", "N")]
    bn_plane = float(statistics.median(cartesian[index][2] for index in bn_ids))
    layers = split_framework_layers(cartesian, framework, number_of_layers)
    layer_planes = [
        float(statistics.median(cartesian[index][2] for index in layer))
        for layer in layers
    ]
    order_near_bn = sorted(
        range(len(layers)), key=lambda index: abs(layer_planes[index] - bn_plane)
    )
    surface_plane = layer_planes[order_near_bn[0]]
    return layers, layer_planes, order_near_bn, bn_ids, surface_plane, bn_plane


def write_poscar(path, comment, lattice, symbols, counts, fractional, flags):
    with path.open("w") as handle:
        handle.write(comment.rstrip() + " | z constraints\n")
        handle.write("1.0\n")
        for vector in lattice:
            handle.write("  " + "  ".join(f"{value:20.14f}" for value in vector) + "\n")
        handle.write("  " + "  ".join(symbols) + "\n")
        handle.write("  " + "  ".join(str(value) for value in counts) + "\n")
        handle.write("Selective dynamics\n")
        handle.write("Direct\n")
        for position, atom_flags in zip(fractional, flags):
            handle.write(
                "  " + "  ".join(f"{value:20.16f}" for value in position)
                + "   " + " ".join(atom_flags) + "\n"
            )


def main():
    parser = argparse.ArgumentParser(
        description="Apply z-only constraints to a metal/defective-h-BN POSCAR."
    )
    parser.add_argument("input", type=Path, help="Input POSCAR or CONTCAR")
    parser.add_argument("output", type=Path, help="Output constrained POSCAR")
    parser.add_argument("--target-gap", type=float, default=None,
                        help="Optional far-field surface-to-hBN spacing in Angstrom")
    parser.add_argument("--local-radius", type=float, default=3.5,
                        help="Strict in-plane cutoff (< radius) for fully free local B/N atoms")
    parser.add_argument("--layer-tolerance", type=float, default=0.35,
                        help="z tolerance used to identify an isolated migrating metal atom")
    parser.add_argument(
        "--migrant-atom",
        type=int,
        default=None,
        help=(
            "One-based overall POSCAR atom number of the migrating metal. "
            "Use this when it is still part of a dense surface layer."
        ),
    )
    parser.add_argument(
        "--flags-reference",
        type=Path,
        default=None,
        help=(
            "Constrained POSCAR whose atom-by-atom T/F flags are copied exactly "
            "to the output; coordinates always come from the input file."
        ),
    )
    parser.add_argument("--metal-layers", type=int, default=4,
                        help="Number of dense framework metal layers to detect (default: 4)")
    parser.add_argument("--fixed-near-layers", type=int, default=2,
                        help="Number of metal layers nearest h-BN to fix in z (default: 2)")
    args = parser.parse_args()

    if args.input.resolve() == args.output.resolve():
        parser.error(
            "input and output must be different files; use, for example, "
            "'POSCAR POSCAR_zfix'"
        )

    comment, lattice, symbols, counts, elements, fractional, _ = read_poscar(args.input)
    if "B" not in symbols or "N" not in symbols:
        raise ValueError("B and N must be present in the element line.")
    metal_symbols = [symbol for symbol in symbols if symbol not in ("B", "N")]
    if len(metal_symbols) != 1:
        raise ValueError("Expected exactly one metal species plus B and N.")
    if not 0 < args.fixed_near_layers <= args.metal_layers:
        raise ValueError("--fixed-near-layers must be between 1 and --metal-layers.")
    if args.local_radius <= 0:
        raise ValueError("--local-radius must be positive.")
    metal = metal_symbols[0]

    cartesian = fractional_to_cartesian(fractional, lattice)
    bn_ids = [index for index, element in enumerate(elements) if element in ("B", "N")]
    bn_plane = float(statistics.median(cartesian[index][2] for index in bn_ids))
    if args.migrant_atom is not None:
        if not 1 <= args.migrant_atom <= len(elements):
            raise ValueError(
                f"--migrant-atom must be between 1 and {len(elements)}."
            )
        migrant = args.migrant_atom - 1
        if elements[migrant] != metal:
            raise ValueError(
                f"Atom #{args.migrant_atom} is {elements[migrant]}, not {metal}."
            )
        migrant_method = "explicit --migrant-atom override"
    else:
        migrant = migrating_metal_index(
            elements, cartesian, metal, bn_plane, args.layer_tolerance
        )
        migrant_method = "automatic isolated-atom detection"
    layers, layer_planes, near_order, bn_ids, surface_plane, bn_plane = analyze_planes(
        elements, cartesian, metal, migrant, args.metal_layers
    )
    initial_gap = abs(bn_plane - surface_plane)

    if args.target_gap is not None:
        if args.target_gap <= 0:
            raise ValueError("--target-gap must be positive.")
        direction = 1.0 if bn_plane >= surface_plane else -1.0
        dz = surface_plane + direction * args.target_gap - bn_plane
        fractional_shift = row_times_matrix([0.0, 0.0, dz], inverse_3x3(lattice))
        for atom in set(bn_ids + [migrant]):
            fractional[atom] = [
                fractional[atom][axis] + fractional_shift[axis] for axis in range(3)
            ]
            fractional[atom][0] %= 1.0
            fractional[atom][1] %= 1.0
        cartesian = fractional_to_cartesian(fractional, lattice)
        layers, layer_planes, near_order, bn_ids, surface_plane, bn_plane = analyze_planes(
            elements, cartesian, metal, migrant, args.metal_layers
        )

    flags = [["T", "T", "T"] for _ in elements]
    fixed_layer_numbers = set(near_order[:args.fixed_near_layers])
    fixed_metal = [
        atom for layer_number in sorted(fixed_layer_numbers) for atom in layers[layer_number]
    ]
    for atom in fixed_metal:
        flags[atom][2] = "F"
    local_bn = [
        atom for atom in bn_ids
        if xy_distance(fractional[atom], fractional[migrant], lattice) < args.local_radius
    ]
    local_bn_set = set(local_bn)
    far_bn = [atom for atom in bn_ids if atom not in local_bn_set]
    for atom in far_bn:
        flags[atom][2] = "F"

    constraints_source = "calculated from this input structure"
    if args.flags_reference is not None:
        (
            _,
            _,
            reference_symbols,
            reference_counts,
            reference_elements,
            _,
            reference_flags,
        ) = read_poscar(args.flags_reference)
        if reference_symbols != symbols or reference_counts != counts:
            raise ValueError(
                "--flags-reference must have the same species order and atom counts "
                "as the input structure."
            )
        if reference_elements != elements:
            raise ValueError(
                "--flags-reference atom ordering does not match the input structure."
            )
        if reference_flags is None:
            raise ValueError(
                "--flags-reference does not contain a Selective dynamics section."
            )
        flags = [atom_flags[:] for atom_flags in reference_flags]
        constraints_source = f"copied exactly from {args.flags_reference}"

        if flags[migrant] != ["T", "T", "T"]:
            raise ValueError(
                "The migrating atom is not T T T in --flags-reference."
            )
        fixed_layer_numbers = set()
        for layer_number, layer in enumerate(layers):
            layer_z_flags = {flags[atom][2] for atom in layer}
            if len(layer_z_flags) != 1:
                raise ValueError(
                    "A detected metal layer has mixed z flags in --flags-reference."
                )
            if layer_z_flags == {"F"}:
                fixed_layer_numbers.add(layer_number)
        fixed_metal = [
            atom
            for atom, element in enumerate(elements)
            if element == metal and atom != migrant and flags[atom][2] == "F"
        ]
        local_bn = [atom for atom in bn_ids if flags[atom] == ["T", "T", "T"]]
        far_bn = [atom for atom in bn_ids if atom not in set(local_bn)]

    write_poscar(args.output, comment, lattice, symbols, counts, fractional, flags)
    final_gap = abs(bn_plane - surface_plane)
    metal_number = sum(element == metal for element in elements[:migrant + 1])
    print(f"Metal species              : {metal}")
    print(f"Migrating atom             : overall #{migrant + 1}, {metal} #{metal_number}")
    print(f"Migrant selection          : {migrant_method}")
    print(f"Constraints source         : {constraints_source}")
    print(f"Initial surface-hBN gap    : {initial_gap:.3f} Angstrom")
    print(f"Output surface-hBN gap     : {final_gap:.3f} Angstrom")
    print("Detected metal layers      :")
    for number, (layer, plane) in enumerate(zip(layers, layer_planes), start=1):
        state = "T T F (z fixed)" if number - 1 in fixed_layer_numbers else "T T T (fully free)"
        print(
            f"  layer {number}: n={len(layer):3d}, z={plane:8.3f} A, "
            f"distance to h-BN={abs(plane - bn_plane):6.3f} A, {state}"
        )
    print(f"Fixed metal-layer atoms    : {len(fixed_metal)}")
    print(f"Migrating {metal} flags          : T T T")
    if args.flags_reference is None:
        print(f"Local B/N (<{args.local_radius:g} A) free  : {len(local_bn)}")
        print(f"Outside-radius B/N z-fixed : {len(far_bn)}")
    else:
        print(f"Reference B/N copied T T T : {len(local_bn)}")
        print(f"Reference B/N copied T T F : {len(far_bn)}")
    print(f"Output                     : {args.output}")


if __name__ == "__main__":
    main()
