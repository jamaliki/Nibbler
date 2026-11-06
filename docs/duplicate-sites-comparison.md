# One site stored twice: before and after

This note records what the duplicate-site fixes change on the 5,000 entries of the second
data audit, so that the corpus rebuild can be checked against it. Three rules changed:

- **site** (conformer selection): single-letter residues that store the competing
  conformers of one site, as separate ligands, residues or chains, keep one letter
  (`docs/structure-preparation.md`, Conformer selection).
- **copy** (assembly expansion): an assembly copy that is wholly partially occupied and
  lies on an earlier retained copy of the same unit is not emitted
  (`docs/structure-preparation.md`, Assembly overlap).
- **whole** (conformer selection): a residue does not keep a conformer that lacks heavy
  atoms another conformer has while a complete one is free (the second audit's class 5).

## Method

Before is `f1939a5`, after is the commit that carries this note. Each build selected
every entry five ways (`occupancy`; `random_clash_aware` with seeds 0 to 3) and expanded
every assembly of the seed-0 random view with the deposited hydrogens removed. A
second pass expanded 104 entries with their hydrogens kept. Each selection is recorded as
source rows, each assembly as its chains with the atom count and a hash of every
coordinate. Two entries, `9l7q` and `9opf`, do not read in either build.

## Result

- **4,913 of 4,998 readable entries are unchanged**: all five selections keep the same
  rows and every assembly has the same chains, atom counts and coordinate hashes.
- **85 entries change**: 55 through the `site` and `copy` rules, and 30 more through
  `whole`.
- Rule by rule, in the build that adds each rule to the previous one:

  | Rule | Entries | Occupancy policy: rows removed / added | Random seed 0: rows removed / added |
  | --- | ---: | --- | --- |
  | site | 49 | 1,642 / 0 | 1,900 / 19 |
  | copy | 7 | 0 / 0 | 0 / 0 |
  | whole | 35 | 14 / 18 | 115 / 218 |

  `copy` changes no selection: it removes 40 assembly copies, 10,252 heavy atoms in nine
  assemblies of seven entries (`2B2E` alone 10,098: the half-occupied RNA stem loop whose
  capsid symmetry copies overlap). Seed 0 ends at 2,015 rows removed and 237 added
  over the 68 entries it changes; the occupancy policy at 1,656 removed and 18 added over
  50.
- The expanded assemblies of the 85 entries lose 12,214 heavy atoms and gain 222.

## Checks that nothing else moves

- Every removed row of the `site` rule belongs to a residue that carries one letter and
  has a heavy atom closer than 1.5 A to a heavy atom of another residue with another
  letter (1,642 rows in the occupancy policy, 1,881 at seed 0, all five configurations
  alike). Every other removed row is a conformer swap inside a residue with several
  letters, and every added row sits in such a residue. No other row is removed or
  added: 0 unexplained in the five configurations.
- A swap arises two ways. A residue next to a removed loser stops being blocked by it
  and takes the conformer the loser had blocked. Or `whole` replaces an incomplete
  conformer with a complete one. The added rows are source rows with their deposited
  coordinates.
- Across the 129 assemblies of the 85 changed entries, 817,438 atoms are in both
  builds. No atom appears in both with other coordinates: comparing the atoms removed
  with the atoms added, no key (chain, residue, atom, letter, occupancy) is on both
  sides.
- The `copy` rule removes only copies whose unit is wholly partial and has an atom within
  1 A of a partially occupied atom of the retained copy, which keeps all its atoms
  (`tests/assembly_overlap.rs`, `tests/assembly_duplicate_sites.rs`).

## Audited hydrogens

The close-hydrogen audit lists 9,204 hydrogens, which the data pipeline adds to heavy
atoms. Counting a hydrogen as cleared when its parent atom or the atom it touches is no
longer in the view, the fixes clear 1,517 of them with the occupancy policy (1,506 with
seed 0):

| Cleared by | Hydrogens |
| --- | ---: |
| a removed assembly copy (`copy`) | 278 |
| a loser chain, or atoms dropped from the selection (`site`) | 1,239 |

The rest stay for reasons the fixes do not touch. Of the audit's 844 hydrogens that touch
another copy of their own chain, 70 are cleared; 514 belong to `2IZ9`. Its protein copies
interpenetrate at full occupancy: the closest atoms belong to different residues, only 1.7%
of a copy's atoms lie within 1 A of an earlier copy, and no copy is a second instance of
the same residues. They are deposited clashes between full-occupancy molecules, so
the rules leave `2IZ9` unchanged and its hydrogens are not duplicates of one site. Others
are not one site by the 1.5 A rule, such as a side-chain hydrogen 1.3 A from the oxygen of
a water of another letter whose parent atom lies 1.7 A away.

## Remaining incomplete conformers

Counting residues whose kept letter lacks a heavy atom that another letter of the
residue has: 4 residues in 4 entries with the occupancy policy before, 13 in 10 entries at
seed 0 (8 to 18 across seeds 1 to 3); after, 2 residues in 2 entries in every
configuration. They are `3LF3` NRQ 66 and `8W6X` CSO/CYS 546. In both the complete conformer has heavy
atoms within 1.5 A of another residue's (a water at 0.4 A in `8W6X`), which would make the
search fall back; whether an authored bond exempts them was not traced.

## Changed entries

`rows` are source `atom_site` rows removed / added. Chains are those missing from the
expanded assembly, naming copies as `chain-n`. Atoms are the net heavy atoms removed from
the expanded assemblies of the entry.

| Entry | Rules | Occupancy rows | Seed 0 rows | Chains removed | Atoms |
| --- | --- | --- | --- | --- | ---: |
| 13rt | site | 15 / 0 | 53 / 0 | C | 37 |
| 1cjx | whole | 0 / 0 | 3 / 4 | - | -1 |
| 1d4x | copy | 0 / 0 | 0 / 0 | F-2 | 3 |
| 1ecf | copy | 0 / 0 | 0 / 0 | H-2 | 9 |
| 1kf4 | site | 1 / 0 | 1 / 0 | - | 1 |
| 1oko | site, whole | 36 / 0 | 36 / 0 | L, O, S | 36 |
| 1sfs | site | 1 / 0 | 1 / 0 | - | 1 |
| 1u77 | whole | 0 / 0 | 0 / 0 | - | 0 |
| 1vlq | whole | 0 / 0 | 0 / 0 | - | 0 |
| 1z1p | whole | 0 / 0 | 0 / 0 | - | 0 |
| 1zo4 | site | 2 / 0 | 2 / 0 | - | 2 |
| 2aot | whole | 0 / 0 | 6 / 7 | - | -1 |
| 2b2e | copy | 0 / 0 | 0 / 0 | B-16, B-23, B-24, B-26, B-28, B-29, B-3, B-33, B-34, B-36, B-38, B-39, B-42, B-43, B-44, B-45, B-46, B-47, B-48, B-49, B-5, B-50, B-52, B-53, B-54, B-55, B-56, B-57, B-58, B-59, B-6, B-60 | 10098 |
| 2bwc | site | 6 / 0 | 6 / 0 | L | 6 |
| 2ckr | site | 1 / 0 | 1 / 0 | - | 1 |
| 2ntx | whole | 0 / 0 | 7 / 8 | - | -3 |
| 2vs2 | whole | 0 / 0 | 6 / 17 | - | -11 |
| 2y65 | site, copy | 7 / 0 | 7 / 0 | G-2 | 95 |
| 3c58 | whole | 0 / 0 | 0 / 0 | - | 0 |
| 3g7z | site | 64 / 0 | 77 / 0 | - | 77 |
| 3gyc | whole | 0 / 0 | 6 / 8 | - | -2 |
| 3gza | site | 1 / 0 | 1 / 0 | - | 1 |
| 3lf3 | whole | 0 / 0 | 11 / 12 | - | -1 |
| 3mu0 | site | 11 / 0 | 20 / 9 | - | 11 |
| 3pe9 | whole | 0 / 0 | 12 / 18 | - | -6 |
| 3qvt | site | 6 / 0 | 17 / 0 | F, F-2, F-3, F-4 | 68 |
| 3s6e | site | 1 / 0 | 1 / 0 | - | 1 |
| 3t0k | site | 2 / 0 | 5 / 3 | E, E-2, F, F-2 | 4 |
| 3uvm | whole | 0 / 0 | 0 / 0 | - | 0 |
| 4c6i | site, whole | 16 / 0 | 18 / 0 | I, I-2 | 24 |
| 4ecu | site | 31 / 0 | 31 / 0 | E, J | 31 |
| 4in6 | whole | 0 / 0 | 0 / 0 | - | 0 |
| 4pnc | whole | 0 / 0 | 6 / 9 | - | -3 |
| 4qm8 | site | 7 / 0 | 7 / 0 | F | 7 |
| 4r34 | site | 52 / 0 | 52 / 0 | E, E-2, O, O-2 | 104 |
| 5a8j | copy | 0 / 0 | 0 / 0 | H-2 | 13 |
| 5avt | site | 2 / 0 | 2 / 0 | G, J | 2 |
| 5dw3 | site | 5 / 0 | 5 / 0 | N | 5 |
| 5kv2 | whole | 0 / 0 | 3 / 4 | - | -1 |
| 5l6u | whole | 0 / 0 | 0 / 0 | - | 0 |
| 5ltr | site, whole | 1 / 0 | 1 / 0 | - | 1 |
| 5mny | whole | 0 / 0 | 0 / 0 | - | 0 |
| 5moo | whole | 0 / 0 | 1 / 4 | - | -3 |
| 5mop | whole | 0 / 0 | 6 / 31 | - | -22 |
| 5pou | site | 1 / 0 | 10 / 6 | F | 4 |
| 5pq7 | site | 1 / 0 | 4 / 0 | G | 4 |
| 5s2e | site | 1 / 0 | 1 / 0 | - | 1 |
| 5s3g | site | 1 / 0 | 1 / 0 | - | 1 |
| 5sou | site | 74 / 0 | 74 / 0 | C, D | 46 |
| 5sp0 | site | 17 / 0 | 42 / 0 | C | 32 |
| 5w2u | site | 20 / 0 | 20 / 0 | H, H-2, H-3, H-4 | 80 |
| 6fji | whole | 0 / 0 | 2 / 11 | - | -6 |
| 6fxq | site | 393 / 0 | 397 / 0 | H, L, O, S, V | 239 |
| 6gew | site | 26 / 0 | 27 / 0 | C, C-2 | 54 |
| 6h85 | copy | 0 / 0 | 0 / 0 | D-2 | 5 |
| 6hc1 | whole | 0 / 0 | 4 / 7 | - | -3 |
| 6i82 | site | 1 / 0 | 1 / 0 | - | 1 |
| 6mb0 | whole | 0 / 0 | 0 / 0 | - | 0 |
| 6s48 | site | 449 / 0 | 450 / 0 | F | 450 |
| 6zrx | whole | 0 / 0 | 12 / 18 | - | -6 |
| 7aza | whole | 0 / 0 | 0 / 0 | - | 0 |
| 7hc5 | site | 17 / 0 | 45 / 0 | C | 31 |
| 7hfw | site | 12 / 0 | 43 / 0 | C | 25 |
| 7ij8 | site | 11 / 0 | 42 / 0 | C | 27 |
| 7jor | whole | 0 / 0 | 12 / 15 | - | -3 |
| 7l0a | whole | 0 / 0 | 0 / 0 | - | 0 |
| 7m45 | site | 29 / 0 | 29 / 0 | F | 29 |
| 7m89 | site | 1 / 0 | 30 / 0 | H | 30 |
| 7m8a | site | 31 / 0 | 31 / 0 | D, J | 31 |
| 7m8c | site | 31 / 0 | 31 / 0 | E, K | 31 |
| 7t5e | whole | 0 / 0 | 1 / 11 | - | -4 |
| 7zcp | site, whole | 6 / 0 | 6 / 0 | - | 16 |
| 8am6 | site, whole | 13 / 11 | 23 / 1 | DA, G, H | 22 |
| 8brn | site | 27 / 0 | 27 / 0 | MD | 27 |
| 8bt1 | whole | 0 / 0 | 4 / 8 | - | -4 |
| 8e8e | site | 31 / 0 | 31 / 0 | G | 31 |
| 8e8f | site | 31 / 0 | 31 / 0 | G | 31 |
| 8e8g | site | 31 / 0 | 31 / 0 | G | 31 |
| 8he6 | whole | 0 / 0 | 1 / 8 | - | -7 |
| 8p73 | site | 30 / 0 | 30 / 0 | D | 30 |
| 8sci | site | 62 / 0 | 62 / 0 | G, L | 62 |
| 9cnt | copy | 0 / 0 | 0 / 0 | E-5 | 36 |
| 9fvy | site | 8 / 0 | 10 / 0 | C | 10 |
| 9kvl | site | 27 / 0 | 27 / 0 | B, B-2, B-3 | 30 |
| 9rej | whole | 6 / 7 | 12 / 18 | - | -6 |
