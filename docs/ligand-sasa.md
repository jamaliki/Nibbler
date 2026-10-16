# Ligand solvent-accessible surface area

`nibbler.mmcif.ligand_sasa` calculates heavy-atom SASA for ligands in their
protein context. It returns `{source_atom_row: area}` with zero-based
`atom_site` rows and areas in square angstroms. Sum the values for the total.

Reuse an existing semantic model to avoid decoding coordinates twice:

```python
import nibbler

model = nibbler.mmcif.read("complex.cif")
areas = nibbler.mmcif.ligand_sasa(model)
total_a2 = sum(areas.values())
```

Paths, bytes, streams and existing `CifDocument` objects also work. Documents
need no second parse or CCD lookup, but SASA must decode their coordinate
fields. This costs more than using a model whose coordinates are already
decoded. Do not construct a semantic model solely to calculate SASA; use the
document path when the rest of the application does not need one.

The operation is explicit. Importing Nibbler or reading a structure does not
calculate SASA, build its grid or initialize its sphere points. Native work
releases the GIL and leaves the input unchanged.

## Selection and chemical convention

- Automatic selection uses `_entity.type` and `label_entity_id`. It includes
  non-polymer and branched entities, including ions and glycans. Modified
  polymer residues remain polymer atoms.
- `atom_rows=[...]` chooses particular source rows, including covalently bound
  ligand atoms. It does not remove the surrounding protein or other ligands
  from the occluding context. Results keep source-row order.
- `model_number=1` selects the deposited model number. Other models never
  occlude it. A missing model raises `ValueError`.
- Water, H, D, T and zero-occupancy atoms neither receive areas nor occlude
  surfaces. Nonzero occupancy counts as a present atom; areas are not
  occupancy-weighted. This heavy-atom convention does not infer missing H.
- Competing alternate labels within a residue raise `ValueError`. Select a
  conformer first; the chosen labels can remain in the selected document:

  ```python
  document = nibbler.cif.read("complex.cif")
  selected = nibbler.mmcif.select_conformer(document).document
  areas = nibbler.mmcif.ligand_sasa(selected)
  ```

- Symmetry copies affect the result only when explicit coordinates contain
  them. Call `mmcif.assembly` first for a biological assembly.
- Missing or invalid coordinates, unknown radii and invalid options raise
  `ValueError`; Nibbler does not substitute a radius for an unknown element.

The default elemental radii follow the Bondi table used by
[Biopython's SASA implementation](https://github.com/biopython/biopython/blob/master/Bio/PDB/SASA.py):
C 1.70, N 1.55, O 1.52, F 1.47, P/S 1.80, Cl 1.75, Br 1.85, I 1.98,
Se 1.90, He 1.40, Na 2.27, Mg 1.73, K 2.75, Ca 2.31, Ni 1.63,
Cu 1.40, Zn 1.39, Cd 1.58 and Hg 1.55 A. Pass `radii={"FE": 2.0}` to set
an element's radius. Python keys ignore case. This convention differs from
ProtOr, which adjusts heavy-atom radii for bonded hydrogens.

## Method and accuracy

The native kernel uses Shrake-Rupley sampling with `probe_radius=1.4` A and
`point_count=256` Fibonacci sphere points. Each selected atom gets a sphere
of radius `vdw_radius + probe_radius`. Nearby spheres hide sample points;
the exposed fraction times `4*pi*r*r` gives its area. Nonselected atoms
receive no samples. A linked-cell grid finds possible occluders without an
all-pairs distance matrix. Compact cell ranges use direct indexing; widely
separated coordinates use a sparse grid to avoid large allocations.

The estimate depends on point count and orientation. Raise `point_count` for
more accurate areas, or lower it when coarse exposure labels suffice. An
atom's sampling increment is `4*pi*r*r/point_count`; this is not an error
bound. Coincident spheres occlude one another, including equal-radius
duplicates. Deduplicate overlapping assembly copies before this calculation.

The tests compare the kernel with analytic sphere intersections, a brute-force
sampler and [Biotite](https://www.biotite-python.org/latest/apidoc/biotite.structure.sasa.html)
using the same radii and sphere points. Matching point sets tests the kernel;
it does not establish convergence at 256 points for every ligand.

For already decoded coordinate arrays, use the same kernel directly:

```python
areas = nibbler.geometry.sasa(
    coordinates, vdw_radii, ligand_atom_indices,
    probe_radius=1.4, point_count=256,
)
```

This array API accepts every supplied atom as an occluder, including H or
water. The caller must filter context, models and conformers. It returns a
list in `ligand_atom_indices` order. Python sequences require native-array
conversion, so reuse `MmcifModel` when one already exists.

## Performance

Run matched, in-memory reads with and without SASA:

```console
python -m benchmarks.ligand_sasa --file complex.cif.gz --warmups 5 --samples 51
```

The benchmark selects one conformer before timing, records source and timed
CIF hashes, alternates measurement order, and reports p50 and p90. It excludes
disk I/O, decompression and conformer selection. It keeps object destruction
outside both read measurements. Compare both total read cost and the separate
SASA measurement; a single timing difference can contain measurement noise.

SASA has a measured cost when requested. A fast syntax-only parse can show a
large relative increase even when SASA takes less than a millisecond. Reusing
decoded model coordinates reduces that added cost. Ligand count, local atom
density, point count and structure size still affect it.

On 2026-10-17, a release build on macOS 26.7 ARM64 with CPython 3.12.13 gave
the following results at 256 points, with seven warmups and 101 samples.
The [JSON report](../benchmarks/ligand-sasa-results/macos-arm64.json) contains
input hashes and all read and SASA timing summaries.

| Entry | Ligand atoms | Model SASA p50 (ms) | Model SASA p90 (ms) | Added model read cost | Added document read cost |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1a1e | 0 | 0.038 | 0.052 | 2.1% | 30.0% |
| 1dpn | 4 | 0.025 | 0.027 | 2.1% | 19.3% |
| 1xan | 74 | 0.181 | 0.215 | 4.9% | 46.4% |
| 1cul | 124 | 0.343 | 0.398 | 6.5% | 49.9% |
| 1a3g | 45 | 0.220 | 0.265 | 3.6% | 43.2% |
| 1qfe | 22 | 0.115 | 0.154 | 3.7% | 46.6% |
| 6w13 | 21 | 0.092 | 0.130 | 3.1% | 35.6% |

These entries test small ligand counts and up to 124 ligand atoms. They do not
establish a bound for ligand-rich assemblies or dense overlapping structures.
