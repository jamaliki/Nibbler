"""Capture processed CCTBX restraints and Probe2 results for native parity."""

import argparse
import contextlib
import hashlib
import importlib.metadata
import io
import json
import os
import sys
import time
from pathlib import Path


def digest(file):
    return hashlib.sha256(Path(file).read_bytes()).hexdigest()


def chemistry(root):
    import libtbx.load_env
    from libtbx.env_config import get_installed_env
    from libtbx.path import absolute_path
    from mmtbx.monomer_library import server

    root = root.resolve(strict=True)
    for env in (libtbx.env, get_installed_env()):
        if env is not None:
            env.repository_paths = [absolute_path(str(root.parent))] + [
                value
                for value in env.repository_paths
                if abs(value) != str(root.parent)
            ]
    for name in ("MMTBX_CCP4_MONOMER_LIB", "CLIBD_MON"):
        os.environ.pop(name, None)
    result = {}
    for parts in [
        ("list", "mon_lib_list.cif"),
        ("list", "geostd_list.cif"),
        ("ener_lib.cif",),
    ]:
        file = server.find_mon_lib_file(relative_path_components=list(parts))
        if file is None or not Path(file).resolve().is_relative_to(root):
            raise ValueError(f"Chemistry lookup escaped the declared tree for {parts}")
        result["/".join(parts)] = dict(file=str(file), sha256=digest(file))
    legacy = root / "geostd/list/mon_lib_list.cif"
    result["reduce_legacy_list"] = dict(file=str(legacy), sha256=digest(legacy))
    return result


def atom_rows(atoms):
    rows = []
    for i, atom in enumerate(atoms):
        labels = atom.fetch_labels()
        rows.append(
            dict(
                index=i,
                i_seq=int(atom.i_seq),
                id=atom.id_str(),
                name=atom.name.strip(),
                element=atom.element.strip(),
                xyz=list(atom.xyz),
                occupancy=atom.occ,
                b_iso=atom.b,
                hetero=atom.hetero,
                charge=atom.charge.strip(),
                model=labels.model_id,
                chain=labels.chain_id,
                resseq=labels.resseq,
                icode=labels.icode,
                component=labels.resname,
                altloc=labels.altloc,
            )
        )
    return rows


def proxy_rows(proxies, attributes):
    result = []
    for proxy in proxies:
        row = {}
        for name in attributes:
            if hasattr(proxy, name):
                value = getattr(proxy, name)
                if callable(value):
                    value = value()
                if not isinstance(value, (bool, float, int, str, type(None))):
                    try:
                        value = list(value)
                    except TypeError:
                        value = str(value)
                row[name] = value
        result.append(row)
    return result


def restraint_graph(model):
    from types import SimpleNamespace

    from mmtbx.validation.restraints import chirality

    geometry = model.get_restraints_manager().geometry
    pairs = geometry.pair_proxies(sites_cart=model.get_sites_cart())
    common = ("i_seqs", "weight", "origin_id")
    result = dict(
        bonds=proxy_rows(
            pairs.bond_proxies.simple,
            (*common, "distance_ideal", "slack", "limit", "top_out"),
        ),
        bonds_asu=proxy_rows(
            pairs.bond_proxies.asu,
            (
                "i_seq",
                "j_seq",
                "j_sym",
                "distance_ideal",
                "weight",
                "slack",
                "origin_id",
            ),
        ),
        angles=proxy_rows(geometry.angle_proxies, (*common, "angle_ideal")),
        dihedrals=proxy_rows(
            geometry.dihedral_proxies,
            (
                *common,
                "angle_ideal",
                "periodicity",
                "alt_angle_ideals",
                "limit",
                "top_out",
            ),
        ),
        chiralities=proxy_rows(
            geometry.chirality_proxies, (*common, "volume_ideal", "both_signs")
        ),
        planarities=proxy_rows(
            geometry.planarity_proxies, ("i_seqs", "weights", "origin_id")
        ),
        symmetry=str(model.get_xray_structure().crystal_symmetry()),
        scope="Final processed geometry manager includes applied link/modification and "
        "CDL target/weight changes. Simple proxy weights and per-atom plane weights "
        "are captured; no unexposed covariance matrix is assumed.",
    )
    atoms = model.get_hierarchy().atoms()
    for proxy in result["chiralities"]:
        atom = atoms[proxy["i_seqs"][0]]
        proxy["pseudochiral"] = chirality(
            atoms_info=[
                SimpleNamespace(resname=atom.fetch_labels().resname, name=atom.name)
            ],
            delta=0.0,
            sigma=1.0,
        ).is_pseudochiral()
    return result


def probe_capture(file, nuclear):
    from iotbx.cli_parser import CCTBXParser
    from mmtbx.programs import probe2
    from mmtbx.programs.clashscore2 import Program
    from mmtbx.validation.clashscore2 import clashscore2

    log = io.StringIO()
    parser = CCTBXParser(program_class=Program, logger=log)
    args = [
        str(file),
        "keep_hydrogens=True",
        "do_flips=False",
        f"nuclear={nuclear}",
        "ignore_missing_restraints=False",
    ]
    parser.parse_args(args)
    before = atom_rows(parser.data_manager.get_model().get_hierarchy().atoms())
    params = parser.working_phil.extract()
    initial, runs = [], []
    original_info = probe2.Helpers.getExtraAtomInfo
    original_run = probe2.Program.run
    fields = dict(
        vdw_radius="vdwRadius",
        is_acceptor="isAcceptor",
        is_donor="isDonor",
        is_dummy_hydrogen="isDummyHydrogen",
        is_ion="isIon",
        charge="charge",
        alt="altLoc",
    )

    def capture_info(*args, **kwargs):
        result = original_info(*args, **kwargs)
        atoms = kwargs["model"].get_hierarchy().atoms()
        indices = {atom: index for index, atom in enumerate(atoms)}
        initial.append(
            dict(
                atoms=atom_rows(atoms),
                extra_atom_info=[
                    dict(
                        index=index,
                        **{
                            name: getattr(
                                result.extraAtomInfo.getMappingFor(atom), attr
                            )
                            for name, attr in fields.items()
                        },
                    )
                    for index, atom in enumerate(atoms)
                ],
                bonds=sorted(
                    {
                        tuple(sorted((indices[atom], indices[neighbor])))
                        for atom, neighbors in kwargs["bondedNeighborLists"].items()
                        for neighbor in neighbors
                        if atom in indices and neighbor in indices
                    }
                ),
            )
        )
        return result

    def capture_run(program):
        value = original_run(program)
        selected = sorted(
            set(program._source_atoms_sorted) | set(program._target_atoms_sorted),
            key=lambda atom: atom.i_seq,
        )
        runs.append(
            dict(
                effective_phil=program.master_phil.format(
                    python_object=program.params
                ).as_str(),
                raw_output=value[1],
                final_selected_atoms=atom_rows(selected),
                final_selected_extra_atom_info=[
                    dict(
                        index=index,
                        **{
                            name: getattr(
                                program._extraAtomInfo.getMappingFor(atom), attr
                            )
                            for name, attr in fields.items()
                        },
                    )
                    for index, atom in enumerate(selected)
                ],
            )
        )
        return value

    probe2.Helpers.getExtraAtomInfo, probe2.Program.run = capture_info, capture_run
    try:
        value = clashscore2(
            params.probe,
            data_manager=parser.data_manager,
            keep_hydrogens=True,
            nuclear=nuclear,
            do_flips=False,
            ignore_missing_restraints=False,
            out=log,
        )
    finally:
        probe2.Helpers.getExtraAtomInfo, probe2.Program.run = (
            original_info,
            original_run,
        )
    return dict(
        clashscore=value.get_clashscore(),
        pairs=[json.loads(row.as_JSON()) for row in value.results],
        input_atoms=before,
        output_atoms=atom_rows(value.hydrogenated_model.get_hierarchy().atoms()),
        initial_atom_info_calls=initial,
        program_runs=runs,
        command_parameters=args,
        effective_phil=parser.working_phil.as_str(),
        log_tail=log.getvalue()[-8192:],
    )


def capture(file, chem_data, nuclear=False, use_cdl=True):
    started, log = time.perf_counter(), io.StringIO()
    result = dict(
        schema="nibbler.native-validation.reference.v1",
        input=dict(file=str(file.resolve()), sha256=digest(file)),
        nuclear=nuclear,
        use_cdl=use_cdl,
        python=sys.version,
        executable=str(Path(sys.executable).resolve()),
        cctbx_version=importlib.metadata.version("cctbx-base"),
        complete=False,
    )
    try:
        with contextlib.redirect_stdout(log):
            result["chemistry_sources"] = chemistry(chem_data)
            import iotbx.pdb
            import mmtbx.model
            from mmtbx.validation import restraints

            source = iotbx.pdb.input(file_name=str(file))
            original = source.construct_hierarchy()
            original.atoms().reset_i_seq()
            result["source_atoms"] = atom_rows(original.atoms())
            model = mmtbx.model.manager(
                model_input=source,
                expand_with_mtrix=False,
                process_biomt=False,
                stop_for_unknowns=True,
                log=log,
            )
            params = model.get_default_pdb_interpretation_params()
            params.pdb_interpretation.use_neutron_distances = nuclear
            params.pdb_interpretation.restraints_library.cdl = use_cdl
            result["interpretation_phil"] = (
                model.get_default_pdb_interpretation_scope()
                .format(python_object=params)
                .as_str()
            )
            model.process(
                pdb_interpretation_params=params,
                make_restraints=True,
                run_clash_guard=False,
            )
            result["atoms"] = atom_rows(model.get_hierarchy().atoms())
            source_sites = {row["id"]: row["xyz"] for row in result["source_atoms"]}
            result["processing_coordinate_changes"] = [
                dict(atom=row["id"], before=source_sites[row["id"]], after=row["xyz"])
                for row in result["atoms"]
                if row["id"] in source_sites and row["xyz"] != source_sites[row["id"]]
            ]
            result["graph"] = restraint_graph(model)
            evaluated = restraints.combined(
                pdb_hierarchy=model.get_hierarchy(),
                xray_structure=model.get_xray_structure(),
                geometry_restraints_manager=model.get_restraints_manager().geometry,
                ignore_hd=False,
                sigma_cutoff=0.0,
            )
            result["restraint_results"] = {
                name: dict(
                    n_total=getattr(evaluated, name).n_total,
                    rows=[
                        json.loads(row.as_JSON())
                        for row in getattr(evaluated, name).results
                    ],
                )
                for name in ("bonds", "angles", "chiralities", "planarities")
            }
            result["probe2"] = probe_capture(file, nuclear)
            result["complete"] = True
    except Exception as error:
        result["error"] = f"{type(error).__name__}: {error}"
    result.update(
        elapsed_seconds=time.perf_counter() - started, log_tail=log.getvalue()[-8192:]
    )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--chem-data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--nuclear", action="store_true")
    parser.add_argument("--no-cdl", action="store_true")
    args = parser.parse_args()
    report = capture(args.input, args.chem_data, args.nuclear, not args.no_cdl)
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(json.dumps({key: report[key] for key in ("complete", "elapsed_seconds")}))
    raise SystemExit(0 if report["complete"] else 1)
