//! Read selected atom rows and the covalent links used by training targets.

use std::borrow::Cow;
use std::collections::{BTreeSet, HashMap};

use crate::cif::{
    BlockCategory, CifBlock, CifDocument, CifLoop, CifValueRef, FastMap, TextColumn, numeric,
};
use crate::pdbx::{AtomSite, ShortText};

pub(crate) struct Atom {
    pub source_row: u32,
    pub name: ShortText,
    pub element: ShortText,
    pub chain: ShortText,
    pub component: ShortText,
    pub entity: ShortText,
    pub auth_sequence: ShortText,
    pub auth_chain: ShortText,
    pub auth_component: ShortText,
    pub insertion: ShortText,
    pub alternate: ShortText,
    pub sequence: i64,
    pub position: [f32; 3],
    pub occupancy: f32,
    pub protein: bool,
    pub group: usize,
}

#[derive(Default)]
pub(crate) struct Source {
    pub atoms: Vec<Atom>,
    pub groups: Vec<Vec<usize>>,
    pub bonds: Vec<(usize, usize)>,
    pub hydrogen_bonds: Vec<(usize, Option<usize>)>,
    pub connections: Vec<Connection>,
    pub has_protein_annotations: bool,
    pub lengths: HashMap<String, i64>,
    pub entity_types: HashMap<String, String>,
}

pub(crate) struct Connection {
    pub row: usize,
    pub atoms: [Option<usize>; 2],
    pub kind: ShortText,
    pub order: ShortText,
    pub symmetries: [ShortText; 2],
}
type Connections = (
    Vec<(usize, usize)>,
    Vec<(usize, Option<usize>)>,
    Vec<Connection>,
);
type Occurrence<'a> = (
    Cow<'a, str>,
    i64,
    Cow<'a, str>,
    Cow<'a, str>,
    Cow<'a, str>,
    Cow<'a, str>,
);
type Residue<'a> = (Cow<'a, str>, i64, Cow<'a, str>, Cow<'a, str>);
type Endpoint<'a> = (&'a str, i64, &'a str, &'a str);
/// The atom names that connections name in each residue, each with its candidate slot.
type Wanted<'a> = FastMap<(&'a str, i64, &'a str), Vec<(&'a str, usize)>>;

fn value<'a>(
    table: &BlockCategory<'a>,
    row: usize,
    column: Option<usize>,
) -> Option<CifValueRef<'a>> {
    column.and_then(|c| table.value(row, c))
}

/// Return `?` for unknown, `.` for not applicable or absent values, and otherwise the
/// text or a number's source spelling. Text borrows from the document.
fn text_of(value: Option<CifValueRef<'_>>) -> Cow<'_, str> {
    match value {
        Some(CifValueRef::Unknown) => Cow::Borrowed("?"),
        Some(CifValueRef::NotApplicable) | None => Cow::Borrowed("."),
        Some(CifValueRef::Text(text)) => Cow::Borrowed(text.as_str()),
        Some(CifValueRef::Integer(number, original)) => original.map_or_else(
            || Cow::Owned(number.to_string()),
            |text| Cow::Borrowed(text.as_str()),
        ),
        Some(CifValueRef::Float(number, _, original)) => original.map_or_else(
            || Cow::Owned(number.to_string()),
            |text| Cow::Borrowed(text.as_str()),
        ),
    }
}

fn text(table: &BlockCategory<'_>, row: usize, column: Option<usize>) -> String {
    text_of(value(table, row, column)).into_owned()
}

fn normal(value: &str) -> &str {
    // Text that starts and ends with a printable ASCII byte has nothing to trim.
    let printable = |byte: Option<&u8>| byte.is_some_and(|byte| (0x21..0x80).contains(byte));
    let bytes = value.as_bytes();
    let value = if printable(bytes.first()) && printable(bytes.last()) {
        value
    } else {
        value.trim()
    };
    if matches!(value, "" | "." | "?") {
        ""
    } else {
        value
    }
}

/// The bytes of `text` when [`normal`] keeps all of them. Reading the bytes skips the
/// UTF-8 check that `as_str` makes.
fn plain(text: &ShortText) -> Option<&[u8]> {
    let bytes = text.as_bytes();
    let printable = |byte: Option<&u8>| byte.is_some_and(|byte| (0x21..0x80).contains(byte));
    (printable(bytes.first()) && printable(bytes.last()) && bytes != b"." && bytes != b"?")
        .then_some(bytes)
}

/// Apply [`normal`] without copying borrowed text.
fn normal_cow<'a>(value: &Cow<'a, str>) -> Cow<'a, str> {
    match value {
        Cow::Borrowed(text) => Cow::Borrowed(normal(text)),
        Cow::Owned(text) => Cow::Owned(normal(text).to_owned()),
    }
}

/// One `atom_site` column, read directly when a columnar loop stores it as text.
#[derive(Clone, Copy)]
enum Reader<'a> {
    Absent,
    Text(TextColumn<'a>),
    Table(usize),
}

impl<'a> Reader<'a> {
    fn new(table: &BlockCategory<'a>, column: Option<usize>) -> Self {
        match column {
            None => Self::Absent,
            Some(column) => table
                .text_column(column)
                .map_or(Self::Table(column), Self::Text),
        }
    }

    #[inline]
    fn get(self, table: &BlockCategory<'a>, row: usize) -> Option<CifValueRef<'a>> {
        match self {
            Self::Absent => None,
            Self::Text(column) => column.value(row),
            Self::Table(column) => table.value(row, column),
        }
    }
}

fn number_of(value: Option<CifValueRef<'_>>) -> Option<f64> {
    match value? {
        CifValueRef::Integer(value, _) => Some(value as f64),
        CifValueRef::Float(value, _, _) => value.is_finite().then_some(value),
        CifValueRef::Text(value) => numeric::parse_float(value.as_str()),
        CifValueRef::Unknown | CifValueRef::NotApplicable => None,
    }
}

fn integer_of(value: Option<CifValueRef<'_>>) -> Option<i64> {
    match value? {
        CifValueRef::Integer(value, _) => Some(value),
        CifValueRef::Text(value) => numeric::parse_integer(value.as_str()),
        _ => None,
    }
}

fn integer(table: &BlockCategory<'_>, row: usize, column: Option<usize>) -> Option<i64> {
    integer_of(value(table, row, column))
}

fn checked_number(
    value: Option<CifValueRef<'_>>,
    row: usize,
    name: &str,
) -> Result<Option<f64>, String> {
    let number = number_of(value);
    let present = value
        .is_some_and(|value| !matches!(value, CifValueRef::Unknown | CifValueRef::NotApplicable));
    if present && number.is_none() {
        return Err(format!(
            "Invalid {name} number in atom_site row {}",
            row + 1
        ));
    }
    Ok(number)
}

fn entity_types(entries: &[crate::cif::CifEntry]) -> HashMap<String, String> {
    let mut result = HashMap::new();
    if let Some(table) = BlockCategory::find(entries, "entity_poly") {
        let entity = table.column("entity_id");
        let kind = table.column("type");
        for row in 0..table.row_count() {
            let entity = text(&table, row, entity);
            let kind = text(&table, row, kind);
            if !normal(&entity).is_empty() && !normal(&kind).is_empty() {
                result.insert(entity, kind);
            }
        }
    }
    result
}

fn entity_lengths(entries: &[crate::cif::CifEntry]) -> HashMap<String, i64> {
    let mut result = HashMap::new();
    if let Some(table) = BlockCategory::find(entries, "entity_poly_seq") {
        let entity = table.column("entity_id");
        let num = table.column("num");
        for row in 0..table.row_count() {
            if let Some(number) = integer(&table, row, num) {
                let entity = text_of(value(&table, row, entity));
                if let Some(length) = result.get_mut(entity.as_ref()) {
                    *length = number.max(*length);
                } else {
                    result.insert(entity.into_owned(), number.max(0));
                }
            }
        }
    }
    result
}

/// Columns that identify one struct_conn partner.
struct PartnerColumns {
    chain: Option<usize>,
    component: Option<usize>,
    name: Option<usize>,
    sequence: Option<usize>,
    alternate: Option<usize>,
    author: Option<usize>,
    author_chain: Option<usize>,
    author_component: Option<usize>,
    insertion: Option<usize>,
}

impl PartnerColumns {
    fn find(table: &BlockCategory<'_>, partner: usize) -> Self {
        let column = |name: &str| table.column(&format!("ptnr{partner}_{name}"));
        Self {
            chain: column("label_asym_id"),
            component: column("label_comp_id"),
            name: column("label_atom_id"),
            sequence: column("label_seq_id"),
            alternate: table.column(&format!("pdbx_ptnr{partner}_label_alt_id")),
            author: column("auth_seq_id"),
            author_chain: column("auth_asym_id"),
            author_component: column("auth_comp_id"),
            insertion: table.column(&format!("pdbx_ptnr{partner}_PDB_ins_code")),
        }
    }
}

/// One partner's normalized identity. `key` is `None` when its chain, component or
/// atom name is missing.
struct Partner<'a> {
    key: Option<Residue<'a>>,
    alternate: Cow<'a, str>,
    author: Cow<'a, str>,
    author_chain: Cow<'a, str>,
    author_component: Cow<'a, str>,
    insertion: Cow<'a, str>,
}

impl<'a> Partner<'a> {
    fn read(table: &BlockCategory<'a>, row: usize, columns: &PartnerColumns) -> Self {
        let field = |column| normal_cow(&text_of(value(table, row, column)));
        let chain = field(columns.chain);
        let component = field(columns.component);
        let name = field(columns.name);
        let key = (!chain.is_empty() && !component.is_empty() && !name.is_empty()).then(|| {
            let sequence = numeric::parse_integer(&field(columns.sequence)).unwrap_or(-1);
            (chain, sequence, component, name)
        });
        Self {
            key,
            alternate: field(columns.alternate),
            author: field(columns.author),
            author_chain: field(columns.author_chain),
            author_component: field(columns.author_component),
            insertion: field(columns.insertion),
        }
    }

    fn resolve(
        &self,
        row: usize,
        partner: usize,
        atoms: &[Atom],
        candidates: &[usize],
    ) -> Result<Option<usize>, String> {
        let matches = |field: &str, atom: &str| field.is_empty() || normal(atom) == field;
        let mut selected = None;
        for &i in candidates {
            let atom = &atoms[i];
            if matches(&self.alternate, &atom.alternate)
                && matches(&self.author, &atom.auth_sequence)
                && matches(&self.author_chain, &atom.auth_chain)
                && matches(&self.author_component, &atom.auth_component)
                && matches(&self.insertion, &atom.insertion)
            {
                if selected.is_some() {
                    return Err(format!(
                        "Connection row {} partner {partner} has ambiguous atom identity",
                        row + 1
                    ));
                }
                selected = Some(i);
            }
        }
        Ok(selected)
    }
}

fn connections(
    entries: &[crate::cif::CifEntry],
    atoms: &[Atom],
    groups: &[Vec<usize>],
) -> Result<Connections, String> {
    let mut pairs = BTreeSet::new();
    let mut hydrogen_pairs = BTreeSet::new();
    let mut records = Vec::new();
    let Some(table) = BlockCategory::find(entries, "struct_conn") else {
        return Ok((Vec::new(), Vec::new(), Vec::new()));
    };
    let columns = [
        PartnerColumns::find(&table, 1),
        PartnerColumns::find(&table, 2),
    ];
    let partners: Vec<[Partner<'_>; 2]> = (0..table.row_count())
        .map(|row| {
            [
                Partner::read(&table, row, &columns[0]),
                Partner::read(&table, row, &columns[1]),
            ]
        })
        .collect();
    // Give each atom identity that some connection names a slot, and list the names
    // that each named residue needs.
    let mut slots: FastMap<Endpoint<'_>, usize> = FastMap::default();
    let mut residues = Wanted::default();
    let partner_slots: Vec<[Option<usize>; 2]> = partners
        .iter()
        .map(|pair| {
            pair.each_ref().map(|partner| {
                let key = partner.key.as_ref()?;
                let next = slots.len();
                let slot = *slots
                    .entry((&*key.0, key.1, &*key.2, &*key.3))
                    .or_insert(next);
                if slot == next {
                    residues
                        .entry((&*key.0, key.1, &*key.2))
                        .or_default()
                        .push((&*key.3, slot));
                }
                Some(slot)
            })
        })
        .collect();
    let mut candidates = vec![Vec::new(); slots.len()];
    if !residues.is_empty() {
        // Atoms of one group share chain, sequence and component, so one lookup per group
        // skips the groups that no connection names. A partner takes the one atom that
        // matches it, so the order of the candidates does not matter.
        for rows in groups {
            let Some(atom) = rows.first().map(|&first| &atoms[first]) else {
                continue;
            };
            let residue = (
                normal(&atom.chain),
                atom.sequence.max(-1),
                normal(&atom.component),
            );
            let Some(names) = residues.get(&residue) else {
                continue;
            };
            for &index in rows {
                let name = &atoms[index].name;
                let name = plain(name).unwrap_or_else(|| normal(name).as_bytes());
                if let Some(&(_, slot)) = names.iter().find(|(wanted, _)| wanted.as_bytes() == name)
                {
                    candidates[slot].push(index);
                }
            }
        }
    }
    let kind_column = table.column("conn_type_id");
    let order_column = table.column("pdbx_value_order");
    let symmetry_columns = [
        table.column("ptnr1_symmetry"),
        table.column("ptnr2_symmetry"),
    ];
    let resolve = |row: usize, partner: usize| -> Result<Option<usize>, String> {
        match partner_slots[row][partner - 1] {
            Some(slot) => {
                partners[row][partner - 1].resolve(row, partner, atoms, &candidates[slot])
            }
            None => Ok(None),
        }
    };
    for row in 0..partners.len() {
        // Short text holds these values inline, so a record allocates nothing.
        let kind = text_of(value(&table, row, kind_column));
        let kind = match normal(&kind) {
            text if text.bytes().any(|byte| byte.is_ascii_uppercase()) => {
                ShortText::new(&text.to_ascii_lowercase())
            }
            text => ShortText::new(text),
        };
        let symmetries = symmetry_columns.map(|column| text_of(value(&table, row, column)));
        let covalent = kind.starts_with("covale") || kind.starts_with("disulf");
        let identity = symmetries
            .iter()
            .all(|value| matches!(normal(value), "" | "1_555"));
        let endpoints = [resolve(row, 1), resolve(row, 2)];
        let [first, second] = if covalent && identity {
            let [first, second] = endpoints;
            [first?, second?]
        } else {
            endpoints.map(|value| value.ok().flatten())
        };
        records.push(Connection {
            row,
            atoms: [first, second],
            kind,
            order: ShortText::new(&text_of(value(&table, row, order_column))),
            symmetries: symmetries.each_ref().map(|value| ShortText::new(value)),
        });
        if !covalent || !identity {
            continue;
        }
        let hydrogen = |row: usize| matches!(atoms[row].element.as_str(), "H" | "D" | "T");
        if let (Some(a), Some(b)) = (first, second) {
            if hydrogen(a) != hydrogen(b) {
                hydrogen_pairs.insert(if hydrogen(a) {
                    (a, Some(b))
                } else {
                    (b, Some(a))
                });
            }
            if atoms[a].group != atoms[b].group {
                pairs.insert((a.min(b), a.max(b)));
            }
        } else if let Some(child) = first.or(second).filter(|&row| hydrogen(row)) {
            hydrogen_pairs.insert((child, None));
        }
    }
    Ok((
        pairs.into_iter().collect(),
        hydrogen_pairs.into_iter().collect(),
        records,
    ))
}

/// The `atom_site` values of one row of the chosen model that [`read`] keeps.
struct SiteRow<'a> {
    row: usize,
    /// The columns of [`FIELDS`], as [`text_of`] gives them.
    fields: [Cow<'a, str>; 11],
    sequence: i64,
    numbers: Numbers<'a>,
}

const FIELDS: [&str; 11] = [
    "label_atom_id",
    "type_symbol",
    "label_asym_id",
    "label_comp_id",
    "label_entity_id",
    "auth_seq_id",
    "pdbx_PDB_ins_code",
    "label_alt_id",
    "group_PDB",
    "auth_asym_id",
    "auth_comp_id",
];

/// The cells of a row's coordinates and occupancy, still to check.
struct Numbers<'a>([Option<CifValueRef<'a>>; 4]);

impl Numbers<'_> {
    fn position(&self, row: usize) -> Result<[f32; 3], String> {
        let mut position = [f32::NAN; 3];
        for (axis, slot) in position.iter_mut().enumerate() {
            *slot = checked_number(self.0[axis], row, "coordinate")?.map_or(f32::NAN, |n| n as f32);
        }
        Ok(position)
    }

    fn occupancy(&self, row: usize) -> Result<Option<f64>, String> {
        checked_number(self.0[3], row, "occupancy")
    }
}

/// Return a site's optional text, or for an unset value the `?` or `.` of its cell.
fn site_text<'a>(
    text: Option<&'a ShortText>,
    column: Reader<'a>,
    table: &BlockCategory<'a>,
    row: usize,
) -> Cow<'a, str> {
    match text {
        Some(text) => Cow::Borrowed(text.as_str()),
        None => text_of(column.get(table, row)),
    }
}

/// The atoms, their occupancies and their residue groups.
type CollectedAtoms = (Vec<Atom>, Vec<f32>, Vec<Vec<usize>>);

/// Build the atoms of `rows` and group them by residue occurrence; also return each
/// atom's occupancy, 0 when absent.
fn collect_atoms<'a>(
    rows: impl Iterator<Item = SiteRow<'a>>,
    capacity: usize,
    types: &HashMap<String, String>,
    protein_types: &FastMap<&str, bool>,
) -> Result<CollectedAtoms, String> {
    let mut atoms = Vec::with_capacity(capacity);
    let mut occupancy = Vec::with_capacity(capacity);
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut occurrences = FastMap::<Occurrence<'a>, usize>::default();
    let mut alternatives = FastMap::<Residue<'a>, Cow<'a, str>>::default();
    let mut previous: Option<(Occurrence<'a>, usize)> = None;
    for SiteRow {
        row,
        fields,
        sequence,
        numbers,
    } in rows
    {
        let chain = normal_cow(&fields[2]);
        let insertion = normal_cow(&fields[6]);
        let author = normal_cow(&fields[5]);
        let alternate = normal_cow(&fields[7]);
        if !alternate.is_empty() {
            let residue = (
                chain.clone(),
                sequence.max(-1),
                insertion.clone(),
                author.clone(),
            );
            match alternatives.get(&residue) {
                Some(choice) if *choice != alternate => {
                    return Err(
                        "Resolve alternate conformers before native target packing".to_owned()
                    );
                }
                Some(_) => {}
                None => {
                    alternatives.insert(residue, alternate);
                }
            }
        }
        let occurrence = (
            chain,
            sequence.max(-1),
            normal_cow(&fields[3]),
            normal_cow(&fields[4]),
            insertion,
            author,
        );
        // Atoms of one residue are usually adjacent, so most rows skip the map.
        let group = match &previous {
            Some((key, group)) if *key == occurrence => *group,
            _ => {
                let next = groups.len();
                let group = *occurrences.entry(occurrence.clone()).or_insert(next);
                if group == next {
                    groups.push(Vec::new());
                }
                previous = Some((occurrence, group));
                group
            }
        };
        let protein = if types.is_empty() {
            normal(&fields[8]).eq_ignore_ascii_case("ATOM")
        } else {
            protein_types
                .get(normal(&fields[4]))
                .copied()
                .unwrap_or(false)
        };
        let position = numbers.position(row)?;
        let source_row =
            u32::try_from(row).map_err(|_| "atom_site row exceeds the source index limit")?;
        let occupied = numbers.occupancy(row)?;
        groups[group].push(atoms.len());
        let [
            name,
            element,
            chain,
            component,
            entity,
            auth_sequence,
            insertion,
            alternate,
            _,
            auth_chain,
            auth_component,
        ] = fields;
        atoms.push(Atom {
            source_row,
            name: ShortText::new(&name),
            element: ShortText::new(&element),
            chain: ShortText::new(&chain),
            component: ShortText::new(&component),
            entity: ShortText::new(&entity),
            auth_sequence: ShortText::new(&auth_sequence),
            auth_chain: ShortText::new(&auth_chain),
            auth_component: ShortText::new(&auth_component),
            insertion: ShortText::new(&insertion),
            alternate: ShortText::new(&alternate),
            sequence,
            position,
            occupancy: occupied.map_or(1.0, |n| n as f32),
            protein,
            group,
        });
        occupancy.push(occupied.map_or(0.0, |n| n as f32));
    }
    Ok((atoms, occupancy, groups))
}

/// Whether [`normal`] maps `text` to the empty string. Most values are one of the three
/// spellings tested first, which compare as fixed-size arrays.
fn blank(text: &ShortText, spellings: &[ShortText; 3]) -> bool {
    spellings.contains(text) || normal(text).is_empty()
}

/// Whether [`plain`] accepts the names of the atoms of `rows` and no two share one, so the
/// choice below keeps every row. Larger groups take the general path.
fn distinct_names(atoms: &[Atom], rows: &[usize]) -> bool {
    const LIMIT: usize = 64;
    if rows.len() <= 1 {
        return true;
    }
    if rows.len() > LIMIT {
        return false;
    }
    // A name of up to seven bytes packs into one integer that starts with its length, so
    // two such names are equal exactly when their integers are.
    let mut keys = [0_u64; LIMIT];
    let mut short = true;
    for (key, &row) in keys.iter_mut().zip(rows) {
        let Some(name) = plain(&atoms[row].name) else {
            return false;
        };
        short &= name.len() <= 7;
        *key = name
            .iter()
            .fold(name.len() as u64, |key, &byte| key << 8 | u64::from(byte));
    }
    if !short {
        return rows.iter().enumerate().all(|(index, &row)| {
            rows[..index]
                .iter()
                .all(|&other| atoms[other].name != atoms[row].name)
        });
    }
    let keys = &mut keys[..rows.len()];
    keys.sort_unstable();
    keys.windows(2).all(|pair| pair[0] != pair[1])
}

/// Group the atoms of `atoms` by their `group` fields, in atom order.
fn group_rows(atoms: &[Atom], count: usize) -> Vec<Vec<usize>> {
    let mut sizes = vec![0_usize; count];
    for atom in atoms {
        sizes[atom.group] += 1;
    }
    let mut groups: Vec<Vec<usize>> = sizes.into_iter().map(Vec::with_capacity).collect();
    for (index, atom) in atoms.iter().enumerate() {
        groups[atom.group].push(index);
    }
    groups
}

/// Build the atoms of the `model` sites as [`collect_atoms`] does for the same rows, but
/// copy each text from its site and skip the residue lookup when a row repeats the
/// identity of the row before it.
///
/// Returns `None` when a site of any model holds a coordinate or occupancy that is not
/// finite. The caller then reads the cells, as it did before it checked the sites.
fn collect_sites<'a>(
    sites: &'a [AtomSite],
    model: i64,
    table: &BlockCategory<'a>,
    types: &HashMap<String, String>,
    protein_types: &FastMap<&str, bool>,
) -> Result<Option<CollectedAtoms>, String> {
    let finite = |site: &AtomSite| {
        [site.x, site.y, site.z, site.occupancy]
            .iter()
            .all(|n| n.is_finite())
    };
    // An error stands only when the sites were usable; otherwise the cells decide.
    let fail = |row: usize, message: String| {
        if sites[row..].iter().all(finite) {
            Err(message)
        } else {
            Ok(None)
        }
    };
    let insertion_column = Reader::new(table, table.column("pdbx_PDB_ins_code"));
    let alternate_column = Reader::new(table, table.column("label_alt_id"));
    let spellings = [ShortText::new(""), ShortText::new("."), ShortText::new("?")];
    let short = |text: Option<&ShortText>, column: Reader<'a>, row: usize| match text {
        Some(text) => text.clone(),
        None => ShortText::new(&text_of(column.get(table, row))),
    };
    let mut atoms: Vec<Atom> = Vec::with_capacity(sites.len());
    let mut occupancy = Vec::with_capacity(sites.len());
    // Whether each group is a protein residue, when the entities decide it.
    let mut proteins: Vec<bool> = Vec::new();
    let mut occurrences = FastMap::<Occurrence<'a>, usize>::default();
    let mut alternatives = FastMap::<Residue<'a>, Cow<'a, str>>::default();
    let mut previous: Option<(Occurrence<'a>, usize)> = None;
    for (row, site) in sites.iter().enumerate() {
        if !finite(site) {
            return Ok(None);
        }
        if site.model_number != model {
            continue;
        }
        let sequence = site.label_seq_id.unwrap_or(-1);
        let insertion = short(site.insertion_code.as_ref(), insertion_column, row);
        let alternate = short(site.label_alt_id.as_ref(), alternate_column, row);
        let text =
            |value: Option<&'a ShortText>, column: Reader<'a>| site_text(value, column, table, row);
        if !blank(&alternate, &spellings) {
            let alternate = normal_cow(&text(site.label_alt_id.as_ref(), alternate_column));
            let residue = (
                Cow::Borrowed(normal(&site.label_asym_id)),
                sequence.max(-1),
                normal_cow(&text(site.insertion_code.as_ref(), insertion_column)),
                Cow::Borrowed(normal(&site.auth_seq_id)),
            );
            match alternatives.get(&residue) {
                Some(choice) if *choice != alternate => {
                    return fail(
                        row,
                        "Resolve alternate conformers before native target packing".to_owned(),
                    );
                }
                Some(_) => {}
                None => {
                    alternatives.insert(residue, alternate);
                }
            }
        }
        // Equal texts normalize equally, so a row that repeats the identity of the row
        // before it joins that row's group.
        let last = atoms.last().filter(|last| {
            last.chain == site.label_asym_id
                && last.sequence == sequence
                && last.component == site.label_comp_id
                && last.entity == site.label_entity_id
                && last.insertion == insertion
                && last.auth_sequence == site.auth_seq_id
        });
        let (group, protein) = match last {
            Some(last) => (last.group, last.protein),
            None => {
                let occurrence = (
                    Cow::Borrowed(normal(&site.label_asym_id)),
                    sequence.max(-1),
                    Cow::Borrowed(normal(&site.label_comp_id)),
                    Cow::Borrowed(normal(&site.label_entity_id)),
                    normal_cow(&text(site.insertion_code.as_ref(), insertion_column)),
                    Cow::Borrowed(normal(&site.auth_seq_id)),
                );
                let group = match &previous {
                    Some((key, group)) if *key == occurrence => *group,
                    _ => {
                        let next = proteins.len();
                        let group = *occurrences.entry(occurrence.clone()).or_insert(next);
                        if group == next {
                            proteins.push(
                                protein_types
                                    .get(occurrence.3.as_ref())
                                    .copied()
                                    .unwrap_or(false),
                            );
                        }
                        previous = Some((occurrence, group));
                        group
                    }
                };
                (group, proteins[group])
            }
        };
        let protein = if types.is_empty() {
            normal(&site.group_pdb).eq_ignore_ascii_case("ATOM")
        } else {
            protein
        };
        let Ok(source_row) = u32::try_from(row) else {
            return fail(
                row,
                "atom_site row exceeds the source index limit".to_owned(),
            );
        };
        atoms.push(Atom {
            source_row,
            name: site.label_atom_id.clone(),
            element: site.element.clone(),
            chain: site.label_asym_id.clone(),
            component: site.label_comp_id.clone(),
            entity: site.label_entity_id.clone(),
            auth_sequence: site.auth_seq_id.clone(),
            auth_chain: site.auth_asym_id.clone(),
            auth_component: site.auth_comp_id.clone(),
            insertion,
            alternate,
            sequence,
            position: [site.x as f32, site.y as f32, site.z as f32],
            occupancy: site.occupancy as f32,
            protein,
            group,
        });
        occupancy.push(site.occupancy as f32);
    }
    let groups = group_rows(&atoms, proteins.len());
    Ok(Some((atoms, occupancy, groups)))
}

pub(crate) fn read(document: &CifDocument, model: i64) -> Result<Source, String> {
    read_with_sites(document, model, None)
}

/// Read as [`read`] does. `sites` are the loop and atom sites that
/// [`crate::pdbx::PdbxModel::loop_sites`] gives for `document`; when they describe the
/// same cells, the rows are read from them instead of being decoded again.
pub(crate) fn read_with_sites(
    document: &CifDocument,
    model: i64,
    sites: Option<(&CifLoop, &[AtomSite])>,
) -> Result<Source, String> {
    let block = document
        .blocks()
        .iter()
        .find(|block| BlockCategory::find(block.entries(), "atom_site").is_some())
        .ok_or_else(|| "Source has no atom_site table".to_owned())?;
    if BlockCategory::find(block.entries(), "pdbx_struct_assembly_gen")
        .is_some_and(|table| table.column("assembly_id").is_some() && table.row_count() > 0)
    {
        return Err("Expand PDB assemblies before tokenization and protonation".to_owned());
    }
    let table = BlockCategory::find(block.entries(), "atom_site")
        .ok_or_else(|| "Source has no atom_site table".to_owned())?;
    let types = entity_types(block.entries());
    let protein_types: FastMap<&str, bool> = types
        .iter()
        .map(|(entity, kind)| {
            (
                entity.as_str(),
                kind.to_ascii_lowercase().starts_with("polypeptide"),
            )
        })
        .collect();
    let row_count = table.row_count();
    let sequence_index = table.column("label_seq_id");
    let model_index = table.column("pdbx_PDB_model_num");
    // A site gives the text `text_of` gives, except for `?` and `.`, which only the cells
    // tell apart. A model reads a typed float as an integer, and `integer_of` does not.
    let sites = sites.and_then(|(only, sites)| {
        let BlockCategory::Loop(found) = &table else {
            return None;
        };
        let fits = std::ptr::eq(*found, only)
            && [sequence_index, model_index]
                .into_iter()
                .flatten()
                .all(|column| !only.holds_floats(column));
        fits.then_some(sites)
    });
    let collected = match sites {
        Some(sites) => collect_sites(sites, model, &table, &types, &protein_types)?,
        None => None,
    };
    let (atoms, occupancy, groups) = match collected {
        Some(collected) => collected,
        None => {
            let columns = FIELDS.map(|name| Reader::new(&table, table.column(name)));
            let model_column = Reader::new(&table, model_index);
            let sequence_column = Reader::new(&table, sequence_index);
            let coordinates = ["Cartn_x", "Cartn_y", "Cartn_z"]
                .map(|name| Reader::new(&table, table.column(name)));
            let occupancy_column = Reader::new(&table, table.column("occupancy"));
            let rows = (0..row_count)
                .filter(|&row| integer_of(model_column.get(&table, row)).unwrap_or(-1) == model)
                .map(|row| SiteRow {
                    row,
                    fields: columns.map(|column| text_of(column.get(&table, row))),
                    sequence: integer_of(sequence_column.get(&table, row)).unwrap_or(-1),
                    numbers: Numbers([
                        coordinates[0].get(&table, row),
                        coordinates[1].get(&table, row),
                        coordinates[2].get(&table, row),
                        occupancy_column.get(&table, row),
                    ]),
                });
            collect_atoms(rows, row_count, &types, &protein_types)?
        }
    };
    finish(block, model, atoms, &occupancy, groups, types)
}

/// The `_atom_site` items of a decoded site, in the order of the masks of [`Copied`].
#[cfg(feature = "reduce3")]
const ITEMS: [&str; 21] = [
    "group_pdb",
    "id",
    "type_symbol",
    "label_atom_id",
    "label_alt_id",
    "label_comp_id",
    "label_asym_id",
    "label_entity_id",
    "label_seq_id",
    "pdbx_pdb_ins_code",
    "cartn_x",
    "cartn_y",
    "cartn_z",
    "occupancy",
    "b_iso_or_equiv",
    "pdbx_formal_charge",
    "auth_seq_id",
    "auth_comp_id",
    "auth_asym_id",
    "auth_atom_id",
    "pdbx_pdb_model_num",
];

/// The items, by their bits in [`ITEMS`] masks, that every row must copy from its source
/// row: the residue fields and the model number.
#[cfg(feature = "reduce3")]
const RESIDUE: u32 =
    1 | 1 << 5 | 1 << 6 | 1 << 7 | 1 << 8 | 1 << 9 | 1 << 16 | 1 << 17 | 1 << 18 | 1 << 20;

/// The items the linked product reads from a heavy site, other than its coordinates.
#[cfg(feature = "reduce3")]
const IDENTITY: u32 = RESIDUE & !1 | 1 << 2 | 1 << 3 | 1 << 4 | 1 << 13 | 1 << 15;

/// The coordinates and the occupancy.
#[cfg(feature = "reduce3")]
const NUMBERS: u32 = 1 << 10 | 1 << 11 | 1 << 12 | 1 << 13;

/// The source that [`read_with_sites`] reads from the model of a document that Reduce3
/// wrote from the document of another model, built from the source of that document.
#[cfg(feature = "reduce3")]
pub(crate) struct Copied<'a> {
    sites: &'a [AtomSite],
    rows: &'a [u32],
    /// The rows that wrote a coordinate or the occupancy, with the values a model decodes.
    numbers: Vec<(u32, [f64; 4])>,
    /// Whether each row copies its atom name.
    names: Vec<bool>,
    heavy: Vec<usize>,
    pub(crate) source: Source,
}

#[cfg(feature = "reduce3")]
impl<'a> Copied<'a> {
    /// Read `document`, whose row `i` copies the site `rows[i]` of `model_of` in the columns
    /// `made[i]` leaves clear (see [`crate::pdbx::build_model_with_copies`]). `base` is the
    /// source of model `model` that [`read_with_sites`] read from `model_of`.
    ///
    /// Returns `None` unless the model of `document` would decode and give the source this
    /// returns. So every row must copy its residue fields and model number, `model`, and
    /// name no alternate; every heavy row must copy the fields the linked product reads,
    /// apart from coordinates; every site of `model_of` must be copied; and the values a row
    /// wrote itself must decode as the model decodes them. Outside `_atom_site` and
    /// `_atom_type`, which the model does not read, the documents must hold the same
    /// entries; and `base` must hold an atom for each site, in order.
    pub(crate) fn new(
        document: &'a CifDocument,
        model_of: &'a crate::pdbx::PdbxModel,
        rows: &'a [u32],
        made: &'a [u64],
        model: i64,
        base: &Source,
    ) -> Option<Self> {
        use crate::pdbx::fields::category_tables;
        let ([block], [original]) = (document.blocks(), model_of.source_document().blocks()) else {
            return None;
        };
        if block.kind() != crate::cif::BlockKind::Data {
            return None;
        }
        let (own, sites) = model_of.loop_sites()?;
        let only = crate::pdbx::only_atom_site_loop(block)?;
        let table = BlockCategory::find(block.entries(), "atom_site")?;
        let BlockCategory::Loop(found) = &table else {
            return None;
        };
        let columns = [
            table.column("label_seq_id"),
            table.column("pdbx_PDB_model_num"),
        ];
        if !std::ptr::eq(*found, only)
            || only.row_count() != rows.len()
            || rows.len() != made.len()
            || only.tags() != own.tags()
            || only.holds_number_spellings()
            || columns
                .into_iter()
                .flatten()
                .any(|column| only.holds_floats(column))
            || base.atoms.len() != sites.len()
            || !base
                .atoms
                .iter()
                .enumerate()
                .all(|(index, atom)| atom.source_row as usize == index)
            || !same_entries(block, original, only, own)
        {
            return None;
        }
        let columns = ITEMS.map(|item| {
            only.tags().iter().position(|tag| {
                crate::cif::split_tag(tag).is_some_and(|(category, name)| {
                    category.eq_ignore_ascii_case("atom_site") && name.eq_ignore_ascii_case(item)
                })
            })
        });
        let bits = columns.map(|column| {
            column
                .and_then(|column| u32::try_from(column).ok())
                .filter(|bit| *bit < u64::BITS)
        });
        let held_by = |made: u64| {
            bits.iter().rev().fold(0, |held, bit| {
                held << 1 | u32::from(bit.is_some_and(|bit| made >> bit & 1 == 0))
            })
        };
        let tables = category_tables(
            crate::cif::CategoryView::new(block.entries()),
            "atom_site",
            ITEMS,
        );
        let [values] = tables.as_slice() else {
            return None;
        };
        if values.rows() != rows.len() {
            return None;
        }
        let (mut element, mut name, mut alternate) =
            (values.column(2), values.column(3), values.column(4));
        let mut numbers = [10, 11, 12, 13].map(|field| values.column(field));
        let (mut b_iso, mut charge) = (values.column(14), values.column(15));
        let [id, auth_name, alternate_column] =
            [1, 19, 4].map(|item| Reader::new(&table, columns[item]));
        // A required text item decodes unless its value is absent, `?` or `.`.
        let present = |column: Reader<'a>, row: usize| {
            !matches!(
                column.get(&table, row),
                None | Some(CifValueRef::Unknown | CifValueRef::NotApplicable)
            )
        };
        let category = "atom_site";
        let mut atoms = Vec::with_capacity(rows.len());
        let mut occupancy = Vec::with_capacity(rows.len());
        let mut copies = Self {
            sites,
            rows,
            numbers: Vec::new(),
            names: Vec::with_capacity(rows.len()),
            heavy: Vec::new(),
            source: Source::default(),
        };
        let mut copied = vec![false; sites.len()];
        // The group of each group of `base`, numbered as the rows first name them.
        let mut groups = vec![usize::MAX; base.groups.len()];
        let mut count = 0;
        let mut last = (u64::MAX, 0);
        for (row, (&source_row, &mask)) in rows.iter().zip(made).enumerate() {
            let source_row = usize::try_from(source_row).ok()?;
            let (site, atom) = (sites.get(source_row)?, base.atoms.get(source_row)?);
            *copied.get_mut(source_row)? = true;
            if last.0 != mask {
                last = (mask, held_by(mask));
            }
            let held = last.1;
            if held & RESIDUE != RESIDUE || site.model_number != model {
                return None;
            }
            let wrote = |item: u32| held >> item & 1 == 0;
            // A model rejects a document whose written values do not decode, so these
            // decodes only check them.
            if wrote(1) && !present(id, row) || wrote(19) && !present(auth_name, row) {
                return None;
            }
            if wrote(14) {
                b_iso
                    .optional_float(row, row, category, "b_iso_or_equiv")
                    .ok()?;
            }
            if wrote(15) {
                charge
                    .optional_charge(row, row, category, "pdbx_formal_charge")
                    .ok()?;
            }
            let element = if wrote(2) {
                element
                    .required_short(row, row, category, "type_symbol")
                    .ok()?
            } else {
                site.element.clone()
            };
            let name = if wrote(3) {
                name.required_short(row, row, category, "label_atom_id")
                    .ok()?
            } else {
                site.label_atom_id.clone()
            };
            // A row with an alternate needs a conformer choice that may drop rows.
            let alternate = if wrote(4) {
                if alternate.optional_short(row).is_some() {
                    return None;
                }
                ShortText::new(&text_of(alternate_column.get(&table, row)))
            } else {
                if site.label_alt_id.is_some() {
                    return None;
                }
                atom.alternate.clone()
            };
            let mut decoded = [site.x, site.y, site.z, site.occupancy];
            if held & NUMBERS != NUMBERS {
                for ((item, column), value) in (10..).zip(&mut numbers).zip(&mut decoded) {
                    if wrote(item) {
                        *value = column
                            .required_float(row, row, category, "coordinate")
                            .ok()?;
                    }
                }
                copies.numbers.push((u32::try_from(row).ok()?, decoded));
            }
            // The read of the model's sites leaves a source with a value that is not finite
            // to the cells.
            if !decoded.iter().all(|number| number.is_finite()) {
                return None;
            }
            if !crate::pdbx::hydrogen_element(&element) {
                if held & IDENTITY != IDENTITY {
                    return None;
                }
                copies.heavy.push(row);
            }
            copies.names.push(!wrote(3));
            let group = groups.get_mut(atom.group)?;
            if *group == usize::MAX {
                *group = count;
                count += 1;
            }
            // The row copies every residue field of the site that gave `atom`.
            atoms.push(Atom {
                source_row: u32::try_from(row).ok()?,
                name,
                element,
                chain: atom.chain.clone(),
                component: atom.component.clone(),
                entity: atom.entity.clone(),
                auth_sequence: atom.auth_sequence.clone(),
                auth_chain: atom.auth_chain.clone(),
                auth_component: atom.auth_component.clone(),
                insertion: atom.insertion.clone(),
                alternate,
                sequence: atom.sequence,
                position: [decoded[0] as f32, decoded[1] as f32, decoded[2] as f32],
                occupancy: decoded[3] as f32,
                protein: atom.protein,
                group: *group,
            });
            occupancy.push(decoded[3] as f32);
        }
        if !copied.iter().all(|&copied| copied) {
            return None;
        }
        let groups = group_rows(&atoms, count);
        copies.source = finish(
            block,
            model,
            atoms,
            &occupancy,
            groups,
            entity_types(block.entries()),
        )
        .ok()?;
        Some(copies)
    }

    /// The heavy rows, in order.
    pub(crate) fn heavy(&self) -> &[usize] {
        &self.heavy
    }

    /// The site of the model that heavy row `row` copies in every field the product reads,
    /// apart from coordinates.
    pub(crate) fn identity(&self, row: usize) -> &'a AtomSite {
        &self.sites[self.rows[row] as usize]
    }

    /// The coordinates of row `row`, as the model of the document decodes them.
    pub(crate) fn position(&self, row: usize) -> [f64; 3] {
        let [x, y, z, _] = match self
            .numbers
            .binary_search_by_key(&row, |&(written, _)| written as usize)
        {
            Ok(index) => self.numbers[index].1,
            Err(_) => {
                let site = self.identity(row);
                [site.x, site.y, site.z, site.occupancy]
            }
        };
        [x, y, z]
    }

    /// Whether row `row` copies the atom name of site `site`, and so holds the same text.
    pub(crate) fn copies_name(&self, row: usize, site: usize) -> bool {
        self.names.get(row) == Some(&true)
            && self
                .rows
                .get(row)
                .is_some_and(|&copied| copied as usize == site)
    }
}

/// Whether two blocks hold the same entries, apart from the `_atom_site` loops `loop_` and
/// `other_loop` and the `_atom_type` entries.
#[cfg(feature = "reduce3")]
fn same_entries(block: &CifBlock, other: &CifBlock, loop_: &CifLoop, other_loop: &CifLoop) -> bool {
    use crate::cif::{CifEntry, in_category};
    fn kept(entry: &CifEntry, skipped: &CifLoop) -> bool {
        match entry {
            CifEntry::Loop(values) => {
                !std::ptr::eq(values, skipped)
                    && !values
                        .tags()
                        .first()
                        .is_some_and(|tag| in_category(tag, "atom_type"))
            }
            CifEntry::Item(item) => !in_category(item.tag(), "atom_type"),
            CifEntry::Frame(_) => true,
        }
    }
    block
        .entries()
        .iter()
        .filter(|entry| kept(entry, loop_))
        .eq(other
            .entries()
            .iter()
            .filter(|entry| kept(entry, other_loop)))
}

/// Keep one atom per name in each residue, then read the connections and entity lengths.
fn finish(
    block: &CifBlock,
    model: i64,
    atoms: Vec<Atom>,
    occupancy: &[f32],
    mut groups: Vec<Vec<usize>>,
    types: HashMap<String, String>,
) -> Result<Source, String> {
    if atoms.is_empty() {
        return Err(format!("model_id {model} is not present in the source"));
    }
    // A selected document can retain shared and selected alternate rows for one name.
    let mut selected = vec![false; atoms.len()];
    {
        let mut by_name = FastMap::<Cow<'_, str>, usize>::default();
        for rows in &groups {
            if distinct_names(&atoms, rows) {
                for &row in rows {
                    selected[row] = true;
                }
                continue;
            }
            by_name.clear();
            for &row in rows {
                let name = normal(&atoms[row].name);
                let name = if name.is_empty() {
                    Cow::Owned(format!("__row_{row}"))
                } else {
                    Cow::Borrowed(name)
                };
                let previous = by_name.get(&name).copied();
                let score = |i: usize| {
                    (
                        !normal(&atoms[i].alternate).is_empty(),
                        if occupancy[i].is_finite() {
                            -occupancy[i]
                        } else {
                            0.0
                        },
                        i,
                    )
                };
                if previous.is_none_or(|i| score(row) < score(i)) {
                    by_name.insert(name, row);
                }
            }
            for &row in by_name.values() {
                selected[row] = true;
            }
        }
    }
    // Rows are rebuilt only when the choice above discarded an atom.
    let atoms: Vec<_> = if selected.iter().all(|&kept| kept) {
        atoms
    } else {
        let atoms: Vec<_> = atoms
            .into_iter()
            .enumerate()
            .filter_map(|(row, atom)| selected[row].then_some(atom))
            .collect();
        for rows in &mut groups {
            rows.clear();
        }
        for (row, atom) in atoms.iter().enumerate() {
            groups[atom.group].push(row);
        }
        atoms
    };
    let (bonds, hydrogen_bonds, connections) = connections(block.entries(), &atoms, &groups)?;
    Ok(Source {
        atoms,
        groups,
        bonds,
        hydrogen_bonds,
        connections,
        has_protein_annotations: !types.is_empty(),
        lengths: entity_lengths(block.entries()),
        entity_types: types,
    })
}

#[cfg(test)]
#[allow(clippy::expect_used)]
#[path = "source_tests.rs"]
mod tests;
