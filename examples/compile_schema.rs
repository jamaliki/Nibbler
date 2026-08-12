//! Compile one verified DDL2 source into a Nibbler release artifact.

use std::env;
use std::error::Error;
use std::fs;
use std::path::Path;

use _core::cif::{compile_dictionary, encode_dictionary, parse};

fn main() -> Result<(), Box<dyn Error>> {
    let arguments = env::args().skip(1).collect::<Vec<_>>();
    let [
        schema_name,
        expected_name,
        expected_version,
        source_sha256,
        input,
        output,
    ] = arguments.as_slice()
    else {
        return Err("usage: compile_schema SCHEMA DICTIONARY VERSION SHA256 INPUT OUTPUT".into());
    };

    let source = fs::read(input)?;
    let document = parse(&source)?;
    let dictionary = compile_dictionary(&document, schema_name, source_sha256)?;
    if dictionary.metadata().dictionary_name() != expected_name {
        return Err(format!(
            "dictionary title mismatch: expected {expected_name:?}, found {:?}",
            dictionary.metadata().dictionary_name()
        )
        .into());
    }
    if dictionary.metadata().version() != expected_version {
        return Err(format!(
            "dictionary version mismatch: expected {expected_version:?}, found {:?}",
            dictionary.metadata().version()
        )
        .into());
    }

    let artifact = encode_dictionary(&dictionary)?;
    let output = Path::new(output);
    if let Some(parent) = output.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::write(output, artifact)?;
    println!(
        "{} {}: {} categories, {} items",
        dictionary.metadata().schema_name(),
        dictionary.metadata().version(),
        dictionary.categories().len(),
        dictionary.items().len(),
    );
    Ok(())
}
