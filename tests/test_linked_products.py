"""Pinned real-entry regressions for native linked-product preparation."""

import pyarrow as pa
import pytest

import nibbler
from benchmarks.linked_products import CASES


@pytest.fixture(scope="module")
def registry(preparation_data):
    return nibbler.components.read(preparation_data.components())


@pytest.fixture(scope="module")
def product(preparation_data, registry):
    def prepare(case, alternate_id=None):
        return nibbler.mmcif.prepare_chemistry(
            preparation_data.entry(case), registry=registry, alternate_id=alternate_id
        ).to_dict()

    return prepare


def atom(graph, component, name, chain=None, sequence=None):
    matches = [
        (i, a)
        for i, a in enumerate(graph["atoms"])
        if a["component_id"] == component
        and a["atom_id"] == name
        and (chain is None or a["asym_id"] == chain)
        and (sequence is None or a["auth_sequence_id"] == sequence)
    ]
    assert matches, (component, name, chain)
    return matches[0]


def order(graph, a, b):
    return next(
        edge["order"]
        for edge in graph["bonds"]
        if edge["active"] and set(edge["atom_indices"]) == {a, b}
    )


@pytest.mark.parametrize(
    "case", [case for case in CASES if case not in {"1vq7", "7p1d"}]
)
def test_real_entries_preserve_observed_heavy_atoms(case, registry, preparation_data):
    file = preparation_data.entry(case)
    document = nibbler.chomp(file)
    before = document.to_canonical()
    graph = nibbler.mmcif.prepare_chemistry(document, registry=registry).to_dict()
    assert document.to_canonical() == before
    sites = (
        pa.RecordBatchReader.from_stream(nibbler.chomp(file, category="atom_site"))
        .read_all()
        .to_pylist()
    )
    for a in graph["atoms"]:
        if a["source_row"] is not None and a["element"] not in {"H", "D", "T"}:
            row = sites[a["source_row"]]
            assert a["coordinates"] == [
                float(row[key]) for key in ("Cartn_x", "Cartn_y", "Cartn_z")
            ]
    assert all(
        graph["atoms"][i]["active"]
        for edge in graph["bonds"]
        if edge["active"]
        for i in edge["atom_indices"]
    )


@pytest.mark.parametrize(
    "case,detail", [("1vq7", "valence 7"), ("7p1d", "saturated backbone")]
)
def test_unsupported_products_fail(case, detail, product):
    with pytest.raises(nibbler.ChemistryError, match=detail):
        product(case)


@pytest.mark.parametrize(
    "case,ligand,center,partner,expected_h",
    [
        ("8qia", "FMN", "C4A", "N5", 1),
        ("1u9x", "IHJ", "C12", "N13", 1),
        ("2zz6", "6AZ", "C6", "C5", 2),
    ],
)
def test_pi_additions_transfer_hydrogens(
    case, ligand, center, partner, expected_h, product
):
    graph = product(case)
    c, _ = atom(graph, ligand, center)
    n, a = atom(graph, ligand, partner)
    assert order(graph, c, n) == (2 if case == "1u9x" else 1)
    assert a["hydrogen_count"] == expected_h
    assert any(edit["reason"] == "hydrogen_transfer" for edit in graph["edits"])


def test_photodimer_consumes_both_pi_bonds(product):
    graph = product("1n4e")
    altered = [
        edge
        for edge in graph["bonds"]
        if edge["active"]
        and edge["template_order"] == 2
        and edge["order"] == 1
        and all(graph["atoms"][i]["component_id"] == "DT" for i in edge["atom_indices"])
    ]
    by_chain = {}
    for edge in altered:
        names = {graph["atoms"][i]["atom_id"] for i in edge["atom_indices"]}
        assert names == {"C5", "C6"}
        chain = graph["atoms"][edge["atom_indices"][0]]["asym_id"]
        by_chain[chain] = by_chain.get(chain, 0) + 1
    assert by_chain and set(by_chain.values()) == {2}


def test_remediated_ver_template_has_one_set_of_hydrogens(product):
    graph = product("1twr")
    for chain in ("C", "E"):
        for name, count in (("C1", 3), ("CMD", 3), ("CBD", 2)):
            _, carbon = atom(graph, "VER", name, chain)
            assert carbon["hydrogen_count"] == count


def test_missing_asn_nitrogen_leaves_but_observed_oxygen_stays(product):
    graph = product("6n0a")
    edits = [
        edit
        for edit in graph["edits"]
        if edit["reason"] == "unresolved_amide_substitution"
    ]
    assert edits
    for edit in edits:
        nitrogen = [
            graph["atoms"][i]
            for i in edit["atom_indices"]
            if graph["atoms"][i]["atom_id"] == "ND2"
        ]
        assert nitrogen and all(
            not a["active"] and a["source_row"] is None for a in nitrogen
        )


@pytest.mark.parametrize(
    "case,component,name,charge",
    [
        ("8r2z", "XR9", "B2", -1),
        ("5lcl", "8AF", "P", 1),
        ("1q29", "G", "O5'", 1),
        ("8q2p", "IRY", "SE1", 1),
    ],
)
def test_product_charges(case, component, name, charge, product):
    graph = product(case)
    values = [
        a["formal_charge"]
        for a in graph["atoms"]
        if a["component_id"] == component and a["atom_id"] == name
    ]
    assert charge in values


def test_coordinated_phosphoryl_product(product):
    graph = product("1cul")
    p, _ = atom(graph, "3PO", "PA")
    o, oxygen = atom(graph, "3PO", "O1A")
    assert order(graph, p, o) == 1
    assert oxygen["formal_charge"] == -1
    assert any(
        edge["active"] and edge["order"] == 0 and o in edge["atom_indices"]
        for edge in graph["bonds"]
    )


def test_glycosylation_distance_exception_is_recorded(product):
    graph = product("2odp")
    assert any("Kept long authored link" in text for text in graph["diagnostics"])


def test_terminal_carboxylate_moves_its_double_bond(product):
    graph = product("1ppv")
    sites = {
        a["atom_id"]: (i, a)
        for i, a in enumerate(graph["atoms"])
        if a["asym_id"] == "A" and a["auth_sequence_id"] == "116"
    }
    carbon, _ = sites["CD"]
    linked, left = sites["OE1"]
    sibling, right = sites["OE2"]
    assert (order(graph, carbon, linked), order(graph, carbon, sibling)) == (1, 2)
    assert (left["formal_charge"], right["formal_charge"]) == (0, 0)


def test_nucleotide_link_removes_missing_op3_and_keeps_observed_op2(product):
    graph = product("6w13")
    _, missing = atom(graph, "8OG", "OP3", "B")
    _, observed = atom(graph, "8OG", "OP2", "B")
    assert missing["source_row"] is None and not missing["active"]
    assert observed["source_row"] is not None and observed["active"]


def test_alternate_links_keep_covalent_and_coordination_states_separate(product):
    graphs = {alternate: product("6dp5", alternate) for alternate in ("B", "C")}
    for alternate, graph in graphs.items():
        p, _ = atom(graph, "C", "P", "C")
        terminal, _ = atom(graph, "U", "O3'", "B")
        _, op3 = atom(graph, "C", "OP3", "C")
        links = [
            edge
            for edge in graph["bonds"]
            if edge["active"] and set(edge["atom_indices"]) == {p, terminal}
        ]
        assert bool(links) == (alternate == "C")
        assert op3["active"] == (alternate == "B")
        assert op3["alternate_id"] == alternate
        if alternate == "B":
            op3_index, _ = atom(graph, "C", "OP3", "C")
            assert any(
                edge["active"]
                and edge["order"] == 0
                and op3_index in edge["atom_indices"]
                for edge in graph["bonds"]
            )


def test_schiff_base_retains_double_link_and_neutral_carbon(product):
    graph = product("1a3g")
    nitrogen, _ = atom(graph, "LYS", "NZ", "A", "159")
    carbon, ligand = atom(graph, "PLP", "C4A", "D")
    assert order(graph, nitrogen, carbon) == 2
    assert ligand["formal_charge"] == 0


@pytest.mark.parametrize("case", ["1ivo", "4js1"])
def test_linked_amide_nitrogen_is_neutral_with_one_hydrogen(case, product):
    graph = product(case)
    linked = set(
        i
        for edge in graph["bonds"]
        if edge["active"] and edge["source"] == "struct_conn"
        for i in edge["atom_indices"]
    )
    amides = [
        graph["atoms"][i]
        for i in linked
        if graph["atoms"][i]["component_id"] == "ASN"
        and graph["atoms"][i]["atom_id"] == "ND2"
    ]
    assert amides and all(
        a["formal_charge"] == 0 and a["hydrogen_count"] == 1 for a in amides
    )


def test_lysine_amine_link_displaces_hydrogen_without_charge_two(product):
    graph = product("1qfe")
    _, nitrogen = atom(graph, "LYS", "NZ", "A", "170")
    assert nitrogen["formal_charge"] == 1 and nitrogen["hydrogen_count"] == 2
