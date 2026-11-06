# Nine audit rows that went from ok to excluded in the 2026-11-06 PDB rebuild

The rebuild of 2026-11-06 moved 291 PDB assemblies from excluded to ok and 9 audit rows
from ok to excluded (v7 to v9, `status-regress`). This note says what the nine are, which
change caused each, and whether the exclusion is right.

**Result.** Four entries are correct exclusions (the new selection is coherent and the
deposited chemistry is a pentacoordinate phosphorus that Nibbler does not support). Two are
regressions: `7V1G` (a merge defect in the component definition, exposed by the new
selection) and `2BWX` (a defect of the complete-conformer rule added in PR #18). Both are
fixed in this branch, each with a test that fails on `main` and passes here. None of the
nine is caused by PR #17 (peptide links) or by the Reduce3 pin.

## The nine rows

`status_regress.py` keys a row by (entry, assembly, model, first altloc seed that gave
the state). States are de-duplicated by their atom rows, so an excluded state can appear
as a new row instead of replacing an old one; a key counts as ok to excluded when the same
key was ok in v7. All nine reproduce locally (see Method).

| Entry | Rows (assembly, seed) | v7 | v9 | Cause | Verdict | Fix |
| --- | --- | --- | --- | --- | --- | --- |
| `5YUX` | a2 s1 | ok | excluded: product valence 6 at `K/900/TTP/PA` | PR #18, `cf77362` (cross-residue site rule) | correct | none |
| `5YUZ` | a2 s1 | ok | excluded: valence 6 at `M/901/TTP/PA` | same | correct | none |
| `5YV0` | a2 s1 | ok | excluded: valence 6 at `N/901/TTP/PA` | same | correct | none |
| `6IG1` | a2 s1 (a2 s3 is a new key, also excluded) | ok | excluded: valence 6 at `N/901/TTP/PA` | same | correct | none |
| `7V1G` | a1 s0, a1 s3 | ok | excluded: valence 10 at `E/304/5GP/P` | `cf77362` selects the Arg23 NH1 A / 5GP A state together; it hits an older defect in `chemistry/input.rs::supplement` (embedded and CCD atom sets are unioned by name, aliases ignored) | **regression** | alias-aware `supplement`; test `a_legacy_embedded_definition_and_the_ccd_aliases_describe_one_atom_set` |
| `2BWX` | a1 s0, s1, s3 | ok | excluded: conflicting component identities (all four seeds) | PR #18, `268b2ef` (prefer the complete conformer) | **regression** | label that re-names shared atoms is not selectable; test `a_microheterogeneous_residue_keeps_one_component_on_every_seed` |

Counting by key: 4 (TTP entries) + 2 (`7V1G`) + 3 (`2BWX`) = 9 rows in 6 assemblies of 6
entries.

## Method

Each entry's mmCIF (the local `~/data/pdb` copy) went through kaveh-data's own PDB index
step, `kaveh_data.novo.indexing.index_pdb_entry`, with the production settings of
`configs/novo-build.yaml` (`exclude_audited`, altloc seeds 0 to 3, `molecular=True`), once
per Nibbler build, and every audit row was compared by the `status_regress` key. Builds:

| Label | Nibbler | Adds |
| --- | --- | --- |
| old | `70fe3fa` (v7 pin) | |
| pl | `6da996d` | PR #17 peptide links |
| c1 | `cf77362` | PR #18: cross-residue site rule |
| c2 | `753e5ec` | PR #18: assembly-copy overlap rule |
| c3 | `268b2ef` | PR #18: prefer the complete conformer |
| ds | `1ae550d` | PR #18 as merged |
| new | `9f86414` (v9 pin) | version bump |
| fix | this branch | the two fixes below |

The index step never runs Reduce3, so the Reduce3 pin cannot change an audit row. kaveh-data
was the rebuild branch (`1d9c1aa`). Local data stands in for the production archives
(`chem_data` 2026.9, `ccd-catalog-v3`, the `~/data/pdb` mmCIF copies) and was not varied, so
a row-for-row match with the Lustre audit is not claimed beyond the nine rows counted above,
which reproduce exactly.

| Entry | old / pl | c1 / c2 | c3 / ds / new | fix |
| --- | --- | --- | --- | --- |
| `5YUX` (and `5YUZ`, `5YV0`) | 4 ok | 7 ok, 1 excluded | same | same |
| `6IG1` | 4 ok | 6 ok, 2 excluded | same | same |
| `7V1G` | 4 ok | 2 ok, 2 excluded | same | 4 ok |
| `2BWX` | 3 ok, 1 excluded | same as old | 4 excluded | 4 ok |
| `5JIW`, `7ZOH`, `2J58` | excluded, excluded, ok | unchanged | unchanged | unchanged |

## The four correct exclusions (`5YUX`, `5YUZ`, `5YV0`, `6IG1`)

`5YUX` was examined atom by atom. The other three show the same message at `PA` of a TTP
and the same status pattern across seeds; I did not inspect their conformers one by one.

These are primer-template complexes that model the incoming TTP and the product of its
incorporation in the same place. The TTP (letter A, occupancy 0.9 for the copy `K` of
`5YUX`) is overlapped by product conformers written with letter B: the extended primer
(`DC 18` B and `DT 19` B, 0.1 each) and the pyrophosphate (`DPO` N, B, 0.1), with contacts
below 1.6 A to the TTP atoms. `struct_conn` covale2 bonds `O3'` of `DC 18` (letter A) to
`PA` of the TTP at 1.58 A, with the TTP's own `PB`-`O3A` bond still intact.

- **v7.** The selection resolved residues one at a time. For seed 1 it kept `DC 18` and
  `DT 19` as letter B and the TTP as letter A. The link names letter A on its primer end,
  so it was not satisfied, was dropped, and the graph built: an incoherent state that holds
  substrate and product within one angstrom of each other, with no link between them.
- **v9.** The site rule (`cf77362`) treats the TTP, the primer residues and the `DPO` as
  one site and keeps one letter. For seed 1 of assembly 2 it keeps `DC 18` A and the TTP
  `K`, so the deposited link is live (O3'-PA 1.58 A). `PA` then carries `O1A` (double), `O2A`, `O3A`,
  `O5'` and the primer `O3'`: valence 6, a phosphorane. Neutral phosphorus accepts 3 or 5
  in Nibbler's valence table, so the exclusion is `PDBX_LINK_PRODUCT`, as designed.
- **Verdict.** The v9 state is the coherent reading of the deposition and the chemistry
  it declares (a covalent bond to a phosphorus that keeps all four of its own substituents)
  is a trapped transition state, not a supported molecule. Excluding it is right. The v7
  `ok` row was never the molecule the entry describes. The assemblies keep their other
  states: the other seeds, in which the primer takes letter B, are still ok (`5YUX`: 7 ok
  and 1 excluded of 8; `6IG1`: 6 and 2 of 8). The cost is one adduct state per entry, not
  the entry.
- Not changed: the chemistry rule (`valence::accepted`) and the site rule. If a transition
  state should be representable, that is a new rule for phosphoranes, not a repair of this
  exclusion.

## `7V1G`: regression in the component merge

`7V1G` has two conformers of Arg23: A (0.41), bonded through `NH1` to the phosphorus of the
ligand `5GP` (chain `E`, letter A, 0.41), and B (0.59) with no ligand.

- **v7.** The residue search chose Arg23 B and the single-letter `5GP` stayed because a
  residue that has one letter was never resolved, so the link (letter A on both ends) was
  not satisfied and was dropped. Again an incoherent state: ligand A with Arg B.
- **v9.** The site rule makes the choice coherent, so seeds 0 and 3 keep Arg23 A with
  `5GP` and the link is live, and seeds 1 and 2 keep Arg23 B without the ligand. The first
  two fail with `Unsupported product valence 10 ... at E/304/5GP/P`.
- **Cause of the valence.** The file embeds `5GP` with the legacy names `O1P`, `O2P`,
  `O3P`. The current CCD calls the same atoms `OP1`, `OP3`, `OP2` and lists the legacy names
  as `alt_atom_id`, with `OP3` the leaving atom. `supplement` merged the two definitions by
  atom name and appended every CCD atom whose name was not in the embedded set, so the
  phosphorus carried both sets of oxygens. A phosphoramidate of `5GP` is a valid valence-5
  phosphorus (two oxygens kept, the ester oxygen, the nitrogen). Renaming the three atoms
  in the file to the current names builds the entry, which isolates the cause.
- **Verdict.** A false rejection, and an old defect that the new selection exposes.
- **Fix.** A CCD atom whose `alt_atom_id` is the name of an embedded atom (and whose
  element agrees) is that atom: it is not appended, the embedded atom takes its leaving
  flag and charge, and its bonds map onto the embedded names. A name that exists in both
  definitions still matches as before. With the fix all four seeds of `7V1G` are ok, and
  two of them now carry the Arg23-`5GP` product.

## `2BWX`: regression in the complete-conformer rule

Residue 249 of chain `A` is a microheterogeneous site: the backbone `N`, `CA`, `C`, `O` is
unlabelled and typed `CYS` (0.7), `CB` and `SG` are `CYS` A (0.4) and B (0.3), and `CSO` C
(0.3) has its own `N`, `CA`, `CB`, `SG`, `C`, `O` at the same coordinates as the shared
ones plus `OD`.

- **v7.** Seeds 0, 1 and 3 chose A or B and built. Seed 2 chose C, kept the `CYS`
  backbone beside it and failed with "conflicting component identities". That already
  excluded one seed in 4.
- **v9.** `268b2ef` ranks conformers by completeness and counts only the lettered heavy
  atoms. A and B name two atoms; C names seven, so only C is complete and is chosen on
  every seed, with the same failure. All four seeds, and so the entry, are excluded.
- **Verdict.** A false rejection of a whole entry, introduced by PR #18.
- **Fix.** In `conformer.rs::groups`, a label that names a heavy atom the shared rows of
  its residue already hold is another identity written over them. It is blocked while
  another label of the residue is free, and it does not define what a complete conformer
  is. All four seeds then build, and seeds 0, 1 and 3 reproduce the v7 states exactly
  (same state hashes). What is lost is the `CSO` state, which neither v7 nor v9 could
  build; keeping it would mean dropping the shared rows when C is chosen, which is a
  larger change (follow-up below).

## The saturated-backbone-link rows

The brief asked whether v9 changes the exclusions with "Sequence does not justify a
saturated backbone link" (163 rows in the audit).

- `5JIW` (`A/370/SNN/C`, 4 rows) and `7ZOH` (the same residue in chains A to D, 16 rows)
  are excluded with that message in every build I ran, old through new and fix.
  `2J58`, which the earlier finding contrasts them with, is ok in every build: it has no
  `SNN` residue (its modified residues are `MSE`), so it is not a counterexample.
- In `5JIW` and `7ZOH` the succinimide `SNN` is linked to its predecessor through `C`-`N`
  and to its successor through `N1`-`CA` (1.43 A), not through its own `C`. That is what
  the message is about, but I did not trace the rule in the code, so the reading that the
  sequence implies a `C`-`N` link the saturated `C` cannot take is an inference. Whether
  the entry should be excluded or modelled as an imide is a separate decision.
- v9 therefore changes none of the three. I could not count the 163 rows: the per-entry
  rows are in the Lustre audit only (below). The rule is not touched by PR #17 or #18.

## Gate

`make check` on this branch (the dev environment's Rust 1.97.1; the Rust suite without
features also passes on 1.88, the MSRV): lint, typecheck, Rust tests for the three feature sets, `pytest` (536 passed, 74 skipped), robustness, docs, corpus,
fixture and release validation, and benchmarks all complete with no failure. The two new
tests fail on `6d10620` with the production messages (`Selected shared and alternate atoms
have conflicting component identities`; `Unsupported product valence`) and pass here.
`.cache/chem-data-v1` was copied from the corpus-rebuild worktree: a fresh fetch failed
`verify_content` on `data_SF4.cif` (17,360 bytes received, 17,228 expected), so the upstream
chemical-component data has moved since the pin; that is unrelated to this change.

## What was not done

- **Lustre files.** The audit rows exist only as `pdb-assembly-audit.jsonl` in each PDB
  build root on Tokyo, which this session could not reach. To check the nine rows against
  the real audit, and to count the 163 saturated-link rows by entry, the two files needed
  are `pdb-assembly-audit.jsonl` of the v7 and the v9 PDB build roots (or the
  `status-regress.json` that `status_regress.py` writes from them).
- **Blast radius of the fixes.** Measured on random samples of the local mmCIF copy
  (253,091 entries, seed 20261106), not on the 4,998 audited entries:
  - Selection, 4,000 entries, five configurations each (`occupancy`; `random_clash_aware`,
    seeds 0 to 3): no selection changes between `9f86414` and this branch. `2BWX`, run as a
    positive control, goes from an error in all five configurations (the occupancy policy
    included) to a selection in all five. The class is rare: it did not occur in the sample.
  - Chemistry, the index step on 240 of those entries (944 audit rows): every state id and
    every status is identical except one row. `9T1J` assembly 1 (`6MZ` at `C/2030`, a modified
    nucleotide in an RNA chain, written with `O1P` in the file) is excluded with `Unsupported
    product valence 9 ... at C/2030/6MZ/P` under v7 and v9 and is `filtered` under the fix.
    That the fix removes it points to the same cause, so the merge defect is older than
    PR #18 and wider than the nine rows: an entry whose embedded component uses legacy names
    and has a link at its phosphorus is at risk.
  - Not measured: the 4,998-entry audit set, and the silent effect of the alias union on
    residues with legacy names that are not linked (extra unobserved atoms in their
    template). Run the comparison of `docs/duplicate-sites-comparison.md` on the audit set,
    with and without this branch, before the corpus is rebuilt on it.
- **The attribution of the 82 link-product excluded-to-ok moves to PR #17** was not tested.
  What this replay shows is narrower: PR #17 (`6da996d`) changes none of the nine rows,
  and the nine entries change at the PR #18 commits. Check the 82 against both PRs
  before repeating the attribution.

## Follow-ups

- Run the audit-set comparison above and keep its result next to the duplicate-sites one.
- Select the `CSO` state of `2BWX` when its occupancy wins, by dropping the shared rows
  of the same names.
- Measure the silent effect of the alias union on residues with legacy names (a scan of
  embedded `_chem_comp_atom` blocks against the CCD's `alt_atom_id`).
- A rule for phosphorane adducts (`5YUX` family), if those states are wanted.
- Decide the `SNN` imide entries (`5JIW`, `7ZOH`, and the other of the 163 rows).
