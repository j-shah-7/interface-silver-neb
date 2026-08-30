# Selective z Constraints for Metal/h-BN Interfaces

This repository contains a Python utility for adding VASP `Selective dynamics`
constraints to layered metal/h-BN structures. It is intended for endpoint
relaxations and NEB preparation in which the metal–h-BN interface distance
should remain approximately fixed while atoms near a vacancy are allowed to
relax locally.

The included `POSCAR_zfix_5p5` file is an example output for a four-layer
Ag/h-BN structure with a B vacancy, one migrating Ag atom, and a 5.5 Å
far-field interface distance.

## Repository contents

| File | Purpose |
| --- | --- |
| `fix_interface_z.py` | Detects the metal layers, interface, migrating metal, and local vacancy region; then writes a constrained VASP structure. |
| `POSCAR_zfix_5p5` | Example constrained POSCAR generated for the Ag/h-BN structure at a 5.5 Å interface distance. |

## Constraint strategy

For the default four-layer metal slab, the program assigns:

| Structural region | VASP flags | Allowed motion |
| --- | --- | --- |
| Two metal layers farthest from h-BN | `T T T` | Fully free |
| Two metal layers closest to h-BN | `T T F` | Free in x/y; fixed in z |
| Migrating metal atom | `T T T` | Fully free |
| B/N atoms with in-plane distance `< 3.5 Å` from the vacancy center | `T T T` | Fully free |
| Remaining B/N atoms | `T T F` | Free in x/y; fixed in z |

The migrating metal atom's x/y projection is used as the vacancy center. The
in-plane distance calculation includes periodic boundary conditions.

This arrangement anchors the two metal planes nearest the interface and the
far-field h-BN plane, but still permits the migrating atom and nearby B/N atoms
to move vertically, pucker, and form local bonds.

## Supported structures

The default settings are designed for structures with:

- one metal species, such as Ag, Au, or Cu;
- B and N as the remaining species;
- a planar metal/h-BN interface normal to Cartesian z;
- four dense metal layers;
- one metal atom removed from a surface layer and positioned as an isolated
  migrating atom near h-BN;
- a vacancy whose in-plane position coincides with the migrating metal atom;
- an h-BN sheet on either side of the metal slab; and
- interface distances such as 5.2–8.2 Å.

The script does not select layers using fixed atom numbers or fixed z
coordinates. It detects the four dense metal planes from the input structure,
measures their distances from the median h-BN plane, and fixes the two nearest
planes. It therefore continues to work when the interface spacing changes.

The input can use `Direct` or `Cartesian` coordinates. A VASP 5-style element
line is required, for example:

```text
Ag B N
192 63 64
```

## Requirements

- Python 3

The script uses only the Python standard library. NumPy, ASE, pymatgen, and
other third-party packages are not required. This makes it suitable for HPC
login and compute nodes that provide a minimal Python installation.

## Basic usage

```bash
python3 fix_interface_z.py INPUT_STRUCTURE OUTPUT_STRUCTURE
```

Both positional filenames are required. Running only
`python3 fix_interface_z.py` displays a missing-arguments error because the
program needs to know which structure to read and where to write the result.

For a file named `POSCAR`:

```bash
python3 fix_interface_z.py POSCAR POSCAR_zfix
```

This command detects and reports the existing interface distance but does not
change it. It adds the appropriate `Selective dynamics` flags and writes the
coordinates in `Direct` format.

After checking the report, the generated file can be used as the next VASP
input:

```bash
cp POSCAR_zfix POSCAR
```

Keep a backup of the unconstrained structure before replacing it.

Do not use the same filename for both arguments. The script intentionally
rejects commands such as `python3 fix_interface_z.py POSCAR POSCAR` to protect
the original structure from being overwritten.

## Initial structures with no displaced metal atom

If the migrating metal is still located inside a complete or nearly complete
surface plane, it is not an isolated atom and cannot be identified reliably
from z coordinates alone. Supply its one-based overall POSCAR atom number with
`--migrant-atom`:

```bash
python3 fix_interface_z.py POSCAR POSCAR_zfix --migrant-atom 170
```

For the included Ag/B-vacancy structure, the migrating atom is overall atom
**#170**, which is also Ag #170 because Ag is the first species. Its x/y
position is aligned with the B vacancy and its three nearest N neighbors.

Use the same explicit atom number for the initial endpoint, final endpoint, and
every NEB image so that the constraints remain atom-for-atom consistent:

```bash
python3 fix_interface_z.py POSCAR POSCAR_zfix --migrant-atom 170
```

If the migrating atom is already displaced from the dense Ag layer, the script
can normally identify it automatically. The explicit option is still preferred
for a related NEB series because it removes ambiguity.

## Structures already prepared from 5.2 to 8.2 Å

If each input POSCAR already has the desired interface spacing, run the script
without `--target-gap`:

```bash
python3 fix_interface_z.py POSCAR POSCAR_zfix
```

Typical output includes the detected interface distance, atom number of the
migrating metal, population and z position of each metal plane, constraints
assigned to each plane, and numbers of local and far-field B/N atoms.

The program was tested on the included structural family at both 5.2 and
8.2 Å. At both limits it selected the same two interface-side Ag layers while
correctly detecting the changed separation.

## Setting an exact interface distance

Use `--target-gap` only when the input structure does not already have the
desired spacing:

```bash
python3 fix_interface_z.py POSCAR POSCAR_zfix_5p2 --target-gap 5.2
python3 fix_interface_z.py POSCAR POSCAR_zfix_6p2 --target-gap 6.2
python3 fix_interface_z.py POSCAR POSCAR_zfix_8p2 --target-gap 8.2
```

The target distance is defined between:

1. the median z coordinate of the nearest dense metal framework layer; and
2. the median z coordinate of all B/N atoms.

To reach the target, the program translates h-BN and the migrating metal
together as one rigid block along Cartesian z. The metal framework remains in
place, and the relative geometry between h-BN and the migrating atom is
preserved.

## Command-line options

Display all available options:

```bash
python3 fix_interface_z.py --help
```

### Change the local vacancy radius

The default local radius is a strict `< 3.5 Å` in-plane cutoff:

```bash
python3 fix_interface_z.py POSCAR POSCAR_zfix --local-radius 3.0
```

Atoms exactly equal to the specified cutoff are classified as outside the
local region and receive `T T F`.

### Use a different number of metal layers

For a six-layer metal slab while still fixing the two layers nearest h-BN:

```bash
python3 fix_interface_z.py POSCAR POSCAR_zfix \
    --metal-layers 6 \
    --fixed-near-layers 2
```

To fix three of the six layers nearest h-BN:

```bash
python3 fix_interface_z.py POSCAR POSCAR_zfix \
    --metal-layers 6 \
    --fixed-near-layers 3
```

### Adjust metal-layer detection

The default z tolerance used to distinguish an isolated migrating atom from a
dense metal layer is 0.35 Å:

```bash
python3 fix_interface_z.py POSCAR POSCAR_zfix \
    --layer-tolerance 0.25
```

A smaller tolerance can help when the migrating atom lies relatively close to
the interface metal plane. A slightly larger tolerance can help when atoms
within one metal plane have substantial z roughness. Always inspect the
reported layer populations after changing this value.

## Included example

`POSCAR_zfix_5p5` was generated with:

```bash
python3 fix_interface_z.py POSCAR POSCAR_zfix_5p5 --target-gap 5.5
```

For this Ag/h-BN example, the program detects:

| Category | Count | Flags |
| --- | ---: | --- |
| Two Ag layers farthest from h-BN | 96 | `T T T` |
| Two Ag layers nearest h-BN | 95 | `T T F` |
| Migrating Ag atom | 1 | `T T T` |
| B/N atoms within `< 3.5 Å` | 12 | `T T T` |
| B/N atoms outside the local radius | 115 | `T T F` |

The nearest Ag layer contains 47 rather than 48 framework atoms because one Ag
atom from that layer is treated as the migrating atom.

## Using the constraints in NEB calculations

For NEB calculations:

- preserve identical atom ordering in the initial state, final state, and all
  intermediate images;
- ensure the same atoms receive the same `Selective dynamics` flags in every
  image;
- keep the migrating metal atom `T T T`, otherwise it cannot cross the
  vacancy;
- use the same `--migrant-atom` value for every image when the migrating atom
  is not isolated in one or more images;
- verify that every processed image reports the same migrating atom number and
  the same layer populations; and
- preferably prepare constrained endpoints first and then regenerate the NEB
  images from those endpoints.

The fixed far-field atoms maintain the intended interface geometry, while the
local B/N atoms remain able to follow the migrating atom through the vacancy.

### Recommended master-reference workflow

Do not calculate the local-radius selection independently for both endpoints.
Atoms can move across the numerical cutoff between the initial and final
structures, producing different constraints even when the same radius is used.

Instead, create the `00` endpoint first and use it as the atom-by-atom flag
master for `07`.

In `00`:

```bash
python3 fix_interface_z.py POSCAR-init POSCAR \
    --migrant-atom 170 \
    --local-radius 3.5
```

In `07`:

```bash
python3 fix_interface_z.py POSCAR-init POSCAR \
    --migrant-atom 170 \
    --flags-reference ../00/POSCAR
```

`--flags-reference` copies every `T/F` flag from the reference POSCAR while
retaining all coordinates and lattice vectors from the current input. It also
checks that the element order and atom counts match and that the migrating atom
is `T T T` in the reference.

For the supplied endpoints, this workflow gives 13 fully free local B/N atoms
in both `00` and `07`, 114 z-fixed B/N atoms, and exactly matching flags for all
319 atoms.

After making the endpoint constraints identical, generate the intermediate NEB
images. If constraints must be reapplied to an existing intermediate image,
use the same master:

```bash
python3 fix_interface_z.py POSCAR-init POSCAR \
    --migrant-atom 170 \
    --flags-reference ../00/POSCAR
```

The reference must already contain a `Selective dynamics` section.

## How the detection works

The program performs the following operations:

1. Reads the lattice, element names, atom counts, and coordinates.
2. Calculates the median h-BN z plane.
3. Finds the isolated metal atom closest to the h-BN plane and identifies it as
   the migrating atom.
4. Removes that atom temporarily from the metal framework.
5. Divides the remaining metal atoms into the requested number of layers using
   the largest gaps in their sorted z coordinates.
6. Ranks the layers by distance from h-BN.
7. Fixes z for the requested number of nearest metal layers.
8. Calculates periodic in-plane B/N distances from the migrating atom.
9. Leaves atoms inside the local cutoff fully free and fixes z outside it.
10. Writes a clean VASP structure with `Selective dynamics` flags.

## Important limitations

The current version is not intended for:

- alloy electrodes containing multiple metal species;
- two separate electrodes or metal/h-BN/metal sandwich structures in the same
  POSCAR;
- multiple isolated migrating metal atoms;
- interfaces whose normal is not Cartesian z;
- strongly curved or highly tilted h-BN sheets;
- structures where the migrating atom is not isolated from the dense metal
  planes; or
- vacancies that are not aligned in x/y with the migrating metal atom.

The output contains the structural coordinates and constraints only. If the
input is a `CONTCAR` containing appended velocity blocks, those velocity lines
are not copied.

Existing selective-dynamics flags are replaced by the flags calculated by this
program.

## Troubleshooting

### `ModuleNotFoundError: No module named 'numpy'`

Download the current version of `fix_interface_z.py`. The present version has
no NumPy dependency and runs with the standard `python3` command:

```bash
python3 fix_interface_z.py POSCAR POSCAR_zfix
```

If the error continues, an older copy of the script is still being executed.
Confirm that the active file does not import NumPy:

```bash
grep -n numpy fix_interface_z.py
```

The current file should produce no matching import line.

### Missing `input` and `output` arguments

The command must contain both filenames:

```bash
python3 fix_interface_z.py POSCAR POSCAR_zfix
```

Use the input first and output second. The script never overwrites the input
unless the same filename is deliberately supplied for both arguments; using a
different output name is strongly recommended.

### `Could not identify one isolated migrating metal atom`

This occurs when the migrating metal is still part of a dense surface layer.
Provide its one-based overall POSCAR atom number explicitly:

```bash
python3 fix_interface_z.py POSCAR POSCAR_zfix --migrant-atom 170
```

For this Ag/h-BN structural series, atom #170 has been confirmed as the metal
atom aligned with the B vacancy.

## Validation before running VASP

Before submitting a calculation, confirm that:

- four sensible metal layers were detected;
- the two layers closest to h-BN are labeled `T T F`;
- the two farther layers are labeled `T T T`;
- the correct migrating metal atom is labeled `T T T`;
- the reported interface distance matches the intended structure;
- the local B/N count is consistent across related structures; and
- no atoms were unintentionally moved when `--target-gap` was omitted.

For a standard relaxation, use an INCAR with ionic relaxation enabled, such as
an appropriate positive `NSW` and suitable `IBRION`, while keeping the lattice
fixed unless cell relaxation is specifically intended.
