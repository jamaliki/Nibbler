from collections.abc import Mapping, Sequence
from typing import Literal, TypeAlias, TypedDict

PredicateSpec: TypeAlias = tuple[str, str, list[str]]
ErrorFields: TypeAlias = tuple[
    str,
    str,
    str | None,
    int | None,
    int | None,
    int | None,
    int | None,
]
DiagnosticFields: TypeAlias = tuple[str, str, str, list[str]]
BondRestraint: TypeAlias = tuple[int, int, float, float, float, bool]
AngleRestraint: TypeAlias = tuple[int, int, int, float, float]
ChiralRestraint: TypeAlias = tuple[int, int, int, int, float, float, bool, bool]
PlaneRestraint: TypeAlias = tuple[Sequence[int], Sequence[float]]

class GeometryCheck(TypedDict):
    status: Literal["checked"]
    count: int

class GeometryEvaluation(TypedDict):
    near_cutoff: bool
    checks: dict[str, GeometryCheck]
    failures: list[dict[str, object]]

class RestraintGraph:
    def __init__(
        self,
        atom_ids: Sequence[str],
        hydrogen_flags: Sequence[bool],
        bonds: Sequence[BondRestraint],
        angles: Sequence[AngleRestraint],
        chiralities: Sequence[ChiralRestraint],
        planes: Sequence[PlaneRestraint],
    ) -> None: ...
    def evaluate(
        self, coordinates: Sequence[Sequence[float]], sigma_cutoff: float = 5.0
    ) -> GeometryEvaluation: ...

ValidationFields: TypeAlias = tuple[
    str,
    str,
    list[str],
    list[DiagnosticFields],
]

class CifDocument:
    @property
    def block_count(self) -> int: ...
    @property
    def schema(self) -> str | None: ...
    def project(
        self,
        category: str,
        *,
        columns: Sequence[str] | None = None,
        schema: str | None = None,
    ) -> CifTable: ...
    def to_canonical(self) -> str: ...
    def to_preserving(self) -> str: ...
    def to_binary(self) -> bytes: ...

class CifTable:
    @property
    def category(self) -> str: ...
    @property
    def columns(self) -> list[str]: ...
    def __len__(self) -> int: ...
    def with_missing(
        self, missing: Literal["collapse", "columns", "extension"]
    ) -> CifTable: ...
    def __arrow_c_stream__(self, requested_schema: object = None) -> object: ...

class PreparedChemistry:
    @property
    def schema_version(self) -> str: ...
    @property
    def atom_count(self) -> int: ...
    @property
    def bond_count(self) -> int: ...
    def to_dict(self) -> dict[str, object]: ...

class ConformerSelection:
    @property
    def atom_rows(self) -> list[int]: ...
    @property
    def search_nodes(self) -> int: ...
    @property
    def model_number(self) -> int | None: ...
    @property
    def schema_version(self) -> str: ...
    @property
    def document(self) -> CifDocument: ...

class MmcifModel:
    @property
    def profile(self) -> Literal["pdbx", "modelcif"]: ...
    @property
    def entry_id(self) -> str: ...
    @property
    def entity_count(self) -> int: ...
    @property
    def asym_unit_count(self) -> int: ...
    @property
    def component_count(self) -> int: ...
    @property
    def atom_site_count(self) -> int: ...
    @property
    def connection_count(self) -> int: ...
    @property
    def entity_kinds(self) -> list[str]: ...
    @property
    def component_ids(self) -> list[str]: ...
    @property
    def prediction_model_count(self) -> int: ...
    @property
    def target_entity_count(self) -> int: ...
    @property
    def template_count(self) -> int: ...
    @property
    def qa_metric_count(self) -> int: ...
    @property
    def qa_value_count(self) -> int: ...
    @property
    def software_names(self) -> list[str]: ...
    @property
    def qa_metric_names(self) -> list[str]: ...
    @property
    def qa_metric_modes(self) -> list[str]: ...
    def source_document(self) -> CifDocument: ...

class Registry:
    def __len__(self) -> int: ...

class _NativeScan:
    def __init__(
        self,
        workers: int,
        category: str | None = None,
        columns: Sequence[str] | None = None,
        predicates: Sequence[PredicateSpec] = (),
        schema: str | None = None,
    ) -> None: ...
    def submit_file(self, file: str) -> None: ...
    def submit_bytes(self, source_name: str, bytes: bytes) -> None: ...
    def next(
        self,
    ) -> tuple[int, CifDocument | CifTable | None, ErrorFields | None] | None: ...
    def close(self) -> None: ...
    def cancel(self) -> None: ...

def native_version() -> str: ...
def contract_version() -> int: ...
def build_profile() -> str: ...
def read_file(
    file: str,
    category: str | None = None,
    columns: Sequence[str] | None = None,
    predicates: Sequence[PredicateSpec] = (),
    schema: str | None = None,
) -> CifDocument | CifTable: ...
def read_bytes(
    source_name: str,
    bytes: bytes,
    category: str | None = None,
    columns: Sequence[str] | None = None,
    predicates: Sequence[PredicateSpec] = (),
    schema: str | None = None,
) -> CifDocument | CifTable: ...
def validate_document(
    document: CifDocument,
    schema: str,
) -> ValidationFields: ...
def build_mmcif_model(
    document: CifDocument,
    profile: Literal["pdbx", "modelcif"],
    registry: Registry | None = None,
) -> MmcifModel: ...
def prepare_mmcif_chemistry(
    model: MmcifModel,
    registry: Registry | None = None,
    model_number: int = 1,
    alternate_id: str | None = None,
) -> PreparedChemistry: ...
def select_mmcif_conformer(
    model: MmcifModel, model_number: int | None = None, seed: int = 0
) -> ConformerSelection: ...
def select_document_conformer(
    document: CifDocument, model_number: int | None = None, seed: int = 0
) -> ConformerSelection: ...
def select_conformer_arrays(
    coordinates: list[list[float]],
    occupancy: list[float],
    residue_ids: list[int],
    chain_ids: list[str],
    sequence_ids: list[int | None],
    atom_names: list[str],
    alternate_ids: list[str | None],
    bonds: list[list[int]],
    clash_mask: list[bool],
    seed: int = 0,
) -> ConformerSelection: ...
def deduplicate_assembly_arrays(
    positions: list[list[float]],
    source_atoms: list[int],
    units: list[int],
    copies: list[int],
    eligible: list[bool],
    occupied: list[bool],
) -> tuple[list[int], list[list[int]]]: ...
def _model_document(
    model: MmcifModel, mirror_local_qa_metric: int | None = None
) -> CifDocument: ...
def build_assembly_document(
    model: MmcifModel, assembly_id: str | None = None
) -> CifDocument: ...
def build_assembly_with_mapping(
    model: MmcifModel, assembly_id: str | None = None
) -> tuple[CifDocument, list[tuple[int, str, str]]]: ...
def build_component_registry(document: CifDocument) -> Registry: ...
def validate_mmcif_document(
    document: CifDocument, profile: Literal["pdbx", "modelcif"]
) -> ValidationFields: ...
def validate_mmcif_model(model: MmcifModel) -> ValidationFields: ...
def reduce_document(
    document: CifDocument, options: Mapping[str, object]
) -> tuple[CifDocument, str, str]: ...
def reduce_required_parent_sources(document: CifDocument, expected_keys: str, *,
                                   model: int = 1, chem_data: str | None = None) -> list[int]: ...
def validate_geometry_document(
    document: CifDocument,
    chem_data: str,
    nuclear: bool = False,
    use_cdl: bool = True,
    sigma_cutoff: float = 5.0,
    include_graph: bool = False,
    restraint_dictionaries: list[CifDocument] | None = None,
    restraint_links: list[tuple[str, int, int]] | None = None,
) -> str: ...
class PreparedObservationSource:
    """Frozen selected-model snapshot, unchanged by later document edits."""
    def ligand_components(self) -> list[str]: ...
    def missing_hydrogen_dictionaries(self, *, chem_data: str) -> list[str]: ...

class ObservationPacker:
    def __init__(self, schema: str, elements: list[str], excluded: list[str]) -> None: ...
    def add_component(self, component: str, features: dict[str, tuple[int, int, int, int]]) -> None: ...
    def add_hydrogen_parents(self, component: str, parents: dict[str, str]) -> None: ...
    def afdb_eligibility(self, payload: object, *, min_chain_residues: int = 50,
                         crop_size: int = 384) -> dict[str, tuple[str, list[int], bytearray]]: ...
    @staticmethod
    def crop_permutation(before: bytes, after: bytes, chains: bytes, ligand: bytes) -> list[int]: ...
    @staticmethod
    def crop_plans(positions: bytes, chains: bytes, ligand: bytes, usable: bytes, *, crop_size: int = 384,
                   centre_eligible: bytes | None = None
                   ) -> tuple[dict[str, tuple[str, list[int], bytearray]], float]: ...
    def ligand_components(self, document: CifDocument, *, model: int = 1) -> list[str]: ...
    def missing_hydrogen_dictionaries(self, document: CifDocument, *, chem_data: str,
                                     model: int = 1) -> list[str]: ...
    def prepare_source(self, document: CifDocument, *, model: int = 1) -> PreparedObservationSource:
        """Copy the selected source once for its queries and target packing."""
        ...
    def pack_prepared(self, source: PreparedObservationSource, *, complete: bool = False,
                      topology_policy: str = "observed_connectivity", protein_only: bool = False,
                      complete_tokens: bool = False, expected_keys: str | None = None,
                      validate: bool = False, with_keys: bool = True, chem_data: str | None = None,
                      target_uncertainties: str | None = None
                      ) -> tuple[dict[str, tuple[str, list[int], bytearray]], list[list[str]]]:
        """Pack a frozen source with the packer's current feature tables."""
        ...
    def protein_document(self, positions: bytes, masks: bytes, identities: bytes,
                         identifier: str) -> tuple[CifDocument, int, str]: ...
    def pack(self, document: CifDocument, *, model: int = 1, complete: bool = False,
             topology_policy: str = "observed_connectivity", protein_only: bool = False,
             complete_tokens: bool = False, expected_keys: str | None = None,
             validate: bool = False, with_keys: bool = True, chem_data: str | None = None,
             target_uncertainties: str | None = None
             ) -> tuple[dict[str, tuple[str, list[int], bytearray]], list[list[str]]]: ...
