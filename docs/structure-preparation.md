# Explicit structure preparation

Nibbler owns source identities, conformer selection, assembly copies, and linked
chemical products. Reduce3 owns hydrogen placement and optimization. Ordinary CIF
and semantic reads keep the deposited state; they do not run these repairs.

The rules and regression set draw on AtomWorks at commit
`1b731b247f750e2f7117806a96b5144f94c898b2`, especially its
[migration guide](https://github.com/RosettaCommons/atomworks/blob/1b731b247f750e2f7117806a96b5144f94c898b2/docs/migration/migration-long.md),
[linked-product rules](https://github.com/RosettaCommons/atomworks/blob/1b731b247f750e2f7117806a96b5144f94c898b2/src/atomworks/io/utils/link_chemistry.py),
and [product regressions](https://github.com/RosettaCommons/atomworks/blob/1b731b247f750e2f7117806a96b5144f94c898b2/tests/io/utils/test_link_product_regressions.py).
Nibbler implements the operations in Rust without a runtime RDKit dependency.
These tests check selected rules; they do not establish full AtomWorks equivalence.
Public test CIF and CCD files download from RCSB into an ignored cache, using the
[hash-pinned test manifest](../tests/preparation_corpus.json). Pytest fetches only
requested coordinates and builds the CCD subset when a product test needs it.
The [test-data instructions](preparation-test-data.md) describe
offline runs and cache checks. Ordinary library calls never fetch this data.

## Linked products

```python
import nibbler

registry = nibbler.components.read("components.cif.gz")
source = nibbler.cif.read("entry.cif.gz")
selected = nibbler.mmcif.select_conformer(source, seed=0)
assert selected.model_number is not None
model = nibbler.mmcif.read(selected.document)
product = nibbler.mmcif.prepare_chemistry(
    model, registry=registry, model_number=selected.model_number
)
print(product.schema_version, product.atom_count, product.bond_count)
```

`linked-product.v1` contains the complete expected component graph, including
unobserved template sites with no coordinates. Removed leaving atoms and bonds
remain in the graph with `active=False`, so indices stay stable. Atom records keep
source rows, template atom names, deposited charges, and initial template charges.
Bond records distinguish component, `struct_conn`, sequence, and hydrogen-transfer
evidence. They also keep authored connection IDs and declared orders. Edits record
each changed neighborhood and its rule; diagnostics record excluded links.

The resolver handles leaving groups, lost bond-order restoration, additions
across double and triple bonds, donor hydrogen transfer, terminal resonance,
unobserved amide substitution, coordinated phosphoryl addition, and supported
charge changes. It excludes coordination from covalent valence. It rejects unknown
linked templates, ambiguous reactions, unsupported orders, nonidentity covalent
symmetry links, and unsupported radical-free product valences. The checked valence
table comes from RDKit 2026.03.6; `tools/gen_product_valence.py` reproduces it.

Sequence links require a declared polymer and consecutive label sequence IDs.
They do not replace an existing attachment to the next residue. Observed endpoints
more than 2.4 A apart do not receive an inferred sequence link. Unobserved endpoints
can connect through declared sequence evidence. Authored links keep their stated
orders; glycosylation roles and carbohydrate links have a recorded 2.4 A distance
exception. The exception does not apply to all covalent links.

Hydrogen counts describe template/product chemistry before a pH assignment.
They do not assert a measured protonation state and do not place missing atoms.
The graph validates changed neighborhoods; it does not sanitize arbitrary unchanged
CCD chemistry. `to_dict()` copies and serializes the graph only on request. A bulk
pipeline should not use that Python export in its per-atom processing loop.

The current kaveh-data ligand-feature builder does not yet derive donor/acceptor,
hybridization, and other RDKit features from this graph. This API provides a tested
chemical state, but it does not by itself fix those downstream features. Wiring that
state into Reduce3 and feature export needs a separate contract and cost check.

## Conformer selection

`mmcif.select_conformer(source, model_number=None, seed=0)` selects the lowest
deposited model unless a number is supplied. `atom_rows` indexes the source table;
`document` constructs the filtered view only when requested. Shared atoms stay.
The view preserves source IDs, filters anisotropic rows, and removes authored links
whose observed endpoints belong only to discarded alternatives or models. Links
to truly unobserved sites remain. Unsupported mixed-category coordinate loops fail.

`clash-aware-altloc.v1` prefers summed occupancy within each residue, with a seeded
tie order. Sparse spatial cells identify heavy-atom clashes below 1.5 A between
residues. Hydrogen occupancies still count, but H, D, and T do not enter contact tests.
Authored covalent, disulfide, and coordination bonds within one crystal image and
adjacent peptide C-N bonds
do not count as clashes. Nearby cross-chain alternatives prefer matching letters
on occupancy ties. A bounded backtracking search finds a consistent assignment;
it does not claim the global maximum-occupancy solution.

No competing alternatives means no contact search. The solver stops after 100,000
choices or a correlated component of more than 512 residues. A conflict raises
`PDBX_CONFORMER_CONFLICT`; the work limit raises `PDBX_CONFORMER_LIMIT`. The selector
does not replace a modified residue with its parent. `6QHP` therefore rejects in
the pinned regression set instead of changing chemical identity.

`mmcif.select_atom_conformers` accepts already parsed arrays and bond exclusions.
This avoids another CIF parse in kaveh-data. The caller must group alternate
component identities into one residue choice and check selected identity coherence.
Its required `clash_mask` excludes hydrogens from contact constraints without
changing their occupancy scores.

## Hydrogen removal

`mmcif.without_hydrogens` copies a one-block document without its deposited
hydrogen sites. A site is a hydrogen when its `type_symbol` is H, D or T in any
case; a site without a symbol stays. `atom_site_anisotrop` keeps the rows whose
`id` names a retained site, so tensors of removed or absent sites go too. A table
that keeps rows keeps its loop or single-row form, an emptied table is dropped,
and `struct_conn` and every other category are copied as authored. A document
without such rows is returned itself. Mixed-category coordinate loops and missing
`atom_site` or `type_symbol` raise `PDBX_HYDROGEN_REMOVAL`.

## Assembly overlap

Assembly expansion removes whole copies of the same free non-polymer unit when
at least half its occupied atoms lie within 1 A of an earlier retained copy. Atom
permutations in a symmetric ligand can satisfy that test. Identity transforms get
priority. Distinct source units never merge. Covalently attached units and
multi-atom coordinated units remain separate because geometric overlap does not
prove that their atom-specific attachments agree.

The rule supplements Nibbler's existing 0.2 A same-source-atom coincidence rule.
It keeps a complete representative ligand and remaps authored connections to the
retained copy. `7MUB` and `1XAN` cover symmetric whole-ligand overlap.

```python
expanded, merges = nibbler.mmcif.assembly_with_mapping(model, "1")
for source_row, removed_label_chain, retained_label_chain in merges:
    print(source_row, removed_label_chain, retained_label_chain)
```

The mapping indexes the input atom_site table and reports whole-unit overlap
merges. It excludes the older partial coincidence rule. It comes from the same
expansion, without a second parse or geometry pass. The array API
`mmcif.deduplicate_assembly_copies` returns retained rows and same-source-atom aliases
for callers that already own expanded coordinates.

kaveh-data uses the native array check after assembly expansion, stores aliases in
shard metadata, and remaps covalent, hydrogen-bond, and coordination edges by
operator path. Its AP Novo preparation audits keep the native assembly map before
Reduce3 runs. Existing shards remain unchanged and need a rebuild to use these rules.

## Cost and qualification

These operations run in preparation, without added training-time I/O or shard-cache
limits. Components resolve once per component ID; authored endpoints use indices.
Neighbor queries borrow adjacency lists, and final hydrogen counts do not allocate
one temporary list per atom. Conformer and overlap contacts use cell grids rather
than dense distance matrices. Python output construction stays explicit.

`benchmarks/preparation_cost.py` measures document reads, semantic reads, selection
including its document view, assembly expansion, and product preparation separately.
It records input/output hashes, all samples, p50/p90, warmups, Python/platform data,
and peak RSS in bytes. Inputs and any CCD registry load before measured operations.
Full-product JSON export runs outside the timed region and is not free.

Use a release extension and the same inputs/build settings for baseline and current
workers. Run workers in alternating order without other CPU tests. Keep ordinary
parse/read comparisons separate from the cost of new optional stages. The local
qualification results belong beside the saved JSON measurements.
The [local qualification report](preparation-performance.md) records latency,
memory, tests, and deployment limits.
