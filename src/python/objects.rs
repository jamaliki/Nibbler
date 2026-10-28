//! Convert `Serialize` values to Python objects without writing JSON text.
//!
//! The result equals `json.loads(serde_json::to_string(value))`. Structs and maps become
//! dicts, sequences become lists, and a non-finite float becomes `None` because
//! `serde_json` writes it as `null`. Equal short texts and field names share one Python
//! string; strings are immutable, so only an identity test can tell.

use std::fmt;

use pyo3::IntoPyObjectExt;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyDict, PyFloat, PyList, PyString};
use serde::ser::{self, Serialize};

use crate::cif::FastMap;

/// Return `value` as Python objects.
pub(super) fn to_python<'py, T: Serialize + ?Sized>(
    py: Python<'py>,
    value: &T,
) -> PyResult<Bound<'py, PyAny>> {
    let mut converter = Converter {
        py,
        names: FastMap::default(),
        texts: FastMap::default(),
        templates: FastMap::default(),
        fields: Vec::new(),
    };
    value
        .serialize(&mut converter)
        .map_err(|error| match error {
            Error::Python(error) => error,
            Error::Message(message) => PyValueError::new_err(message),
        })
}

/// A Python exception, or a message that `serde_json` would also report.
#[derive(Debug)]
pub(super) enum Error {
    Python(PyErr),
    Message(String),
}

impl fmt::Display for Error {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Python(error) => fmt::Display::fmt(error, formatter),
            Self::Message(message) => formatter.write_str(message),
        }
    }
}

impl std::error::Error for Error {}

impl ser::Error for Error {
    fn custom<T: fmt::Display>(message: T) -> Self {
        Self::Message(message.to_string())
    }
}

impl From<PyErr> for Error {
    fn from(error: PyErr) -> Self {
        Self::Python(error)
    }
}

fn key_must_be_a_string() -> Error {
    Error::Message("key must be a string".to_owned())
}

struct Converter<'py> {
    py: Python<'py>,
    /// Field and variant names. They are `'static`, so their address and length
    /// identify them.
    names: FastMap<(usize, usize), Bound<'py, PyString>>,
    /// Texts of at most seven bytes, keyed by the bytes and the length in one word.
    texts: FastMap<u64, Bound<'py, PyString>>,
    /// The first field list seen for each struct name.
    templates: FastMap<(usize, usize), Template<'py>>,
    /// Fields of the structs under construction. Serialization is depth-first, so each
    /// open struct owns a suffix of this list.
    fields: Vec<(&'static str, Bound<'py, PyAny>)>,
}

/// The field names of one struct and a dict that maps them to `None`.
///
/// A dict built by insertion grows its table twice on the way to 17 keys. A copy of the
/// template has the final table already, so filling it only replaces values.
struct Template<'py> {
    names: Vec<&'static str>,
    keys: Vec<Bound<'py, PyString>>,
    dict: Bound<'py, PyDict>,
}

impl Template<'_> {
    fn matches(&self, fields: &[(&'static str, Bound<'_, PyAny>)]) -> bool {
        self.names.len() == fields.len()
            && self.names.iter().zip(fields).all(|(&name, &(field, _))| {
                (name.as_ptr() == field.as_ptr() && name.len() == field.len()) || name == field
            })
    }
}

impl<'py> Converter<'py> {
    fn name(&mut self, name: &'static str) -> Bound<'py, PyString> {
        let py = self.py;
        self.names
            .entry((name.as_ptr() as usize, name.len()))
            .or_insert_with(|| PyString::new(py, name))
            .clone()
    }

    fn text(&mut self, text: &str) -> Bound<'py, PyString> {
        let bytes = text.as_bytes();
        if bytes.len() > 7 {
            return PyString::new(self.py, text);
        }
        let key = bytes
            .iter()
            .enumerate()
            .fold((bytes.len() as u64) << 56, |key, (index, &byte)| {
                key | u64::from(byte) << (8 * index)
            });
        let py = self.py;
        self.texts
            .entry(key)
            .or_insert_with(|| PyString::new(py, text))
            .clone()
    }

    fn none(&self) -> Bound<'py, PyAny> {
        self.py.None().into_bound(self.py)
    }

    fn float(&self, value: f64) -> Bound<'py, PyAny> {
        if value.is_finite() {
            PyFloat::new(self.py, value).into_any()
        } else {
            self.none()
        }
    }

    fn integer<T: IntoPyObject<'py>>(&self, value: T) -> Result<Bound<'py, PyAny>, Error> {
        Ok(value.into_bound_py_any(self.py)?)
    }

    /// Return a dict of the fields from `start` on, and remove them from `fields`.
    fn finish_struct(
        &mut self,
        name: &'static str,
        start: usize,
    ) -> Result<Bound<'py, PyDict>, Error> {
        let identity = (name.as_ptr() as usize, name.len());
        if let Some(template) = self.templates.get(&identity)
            && template.matches(&self.fields[start..])
        {
            let dict = template.dict.copy()?;
            for (key, (_, value)) in template.keys.iter().zip(self.fields.drain(start..)) {
                dict.set_item(key, value)?;
            }
            return Ok(dict);
        }
        let names: Vec<&'static str> = self.fields[start..].iter().map(|&(name, _)| name).collect();
        let keys: Vec<Bound<'py, PyString>> = names.iter().map(|name| self.name(name)).collect();
        let dict = PyDict::new(self.py);
        for (key, (_, value)) in keys.iter().zip(self.fields.drain(start..)) {
            // A repeated key keeps its first position and takes the last value, as in
            // `json.loads`.
            dict.set_item(key, value)?;
        }
        if !self.templates.contains_key(&identity) {
            let template = PyDict::new(self.py);
            for key in &keys {
                template.set_item(key, self.py.None())?;
            }
            self.templates.insert(
                identity,
                Template {
                    names,
                    keys,
                    dict: template,
                },
            );
        }
        Ok(dict)
    }

    /// Return `{name: value}`, the form `serde_json` gives an enum variant with data.
    fn variant(
        &mut self,
        name: &'static str,
        value: Bound<'py, PyAny>,
    ) -> Result<Bound<'py, PyAny>, Error> {
        let dict = PyDict::new(self.py);
        dict.set_item(self.name(name), value)?;
        Ok(dict.into_any())
    }
}

impl<'a, 'py> ser::Serializer for &'a mut Converter<'py> {
    type Ok = Bound<'py, PyAny>;
    type Error = Error;
    type SerializeSeq = Sequence<'a, 'py>;
    type SerializeTuple = Sequence<'a, 'py>;
    type SerializeTupleStruct = Sequence<'a, 'py>;
    type SerializeTupleVariant = Sequence<'a, 'py>;
    type SerializeMap = Map<'a, 'py>;
    type SerializeStruct = Struct<'a, 'py>;
    type SerializeStructVariant = Struct<'a, 'py>;

    fn serialize_bool(self, value: bool) -> Result<Self::Ok, Error> {
        Ok(PyBool::new(self.py, value).to_owned().into_any())
    }

    fn serialize_i8(self, value: i8) -> Result<Self::Ok, Error> {
        self.integer(value)
    }

    fn serialize_i16(self, value: i16) -> Result<Self::Ok, Error> {
        self.integer(value)
    }

    fn serialize_i32(self, value: i32) -> Result<Self::Ok, Error> {
        self.integer(value)
    }

    fn serialize_i64(self, value: i64) -> Result<Self::Ok, Error> {
        self.integer(value)
    }

    fn serialize_i128(self, value: i128) -> Result<Self::Ok, Error> {
        self.integer(value)
    }

    fn serialize_u8(self, value: u8) -> Result<Self::Ok, Error> {
        self.integer(value)
    }

    fn serialize_u16(self, value: u16) -> Result<Self::Ok, Error> {
        self.integer(value)
    }

    fn serialize_u32(self, value: u32) -> Result<Self::Ok, Error> {
        self.integer(value)
    }

    fn serialize_u64(self, value: u64) -> Result<Self::Ok, Error> {
        self.integer(value)
    }

    fn serialize_u128(self, value: u128) -> Result<Self::Ok, Error> {
        self.integer(value)
    }

    fn serialize_f32(self, value: f32) -> Result<Self::Ok, Error> {
        if !value.is_finite() {
            return Ok(self.none());
        }
        // `serde_json` writes the shortest decimal that names the `f32`, and
        // `json.loads` reads that decimal as an `f64`; `f64::from` could differ.
        let text = written(value)?;
        let parsed = text
            .parse::<f64>()
            .map_err(|error| Error::Message(error.to_string()))?;
        Ok(self.float(parsed))
    }

    fn serialize_f64(self, value: f64) -> Result<Self::Ok, Error> {
        Ok(self.float(value))
    }

    fn serialize_char(self, value: char) -> Result<Self::Ok, Error> {
        Ok(self.text(value.encode_utf8(&mut [0; 4])).into_any())
    }

    fn serialize_str(self, value: &str) -> Result<Self::Ok, Error> {
        Ok(self.text(value).into_any())
    }

    fn serialize_bytes(self, value: &[u8]) -> Result<Self::Ok, Error> {
        Ok(PyList::new(self.py, value)?.into_any())
    }

    fn serialize_none(self) -> Result<Self::Ok, Error> {
        Ok(self.none())
    }

    fn serialize_some<T: Serialize + ?Sized>(self, value: &T) -> Result<Self::Ok, Error> {
        value.serialize(self)
    }

    fn serialize_unit(self) -> Result<Self::Ok, Error> {
        Ok(self.none())
    }

    fn serialize_unit_struct(self, _name: &'static str) -> Result<Self::Ok, Error> {
        Ok(self.none())
    }

    fn serialize_unit_variant(
        self,
        _name: &'static str,
        _index: u32,
        variant: &'static str,
    ) -> Result<Self::Ok, Error> {
        Ok(self.name(variant).into_any())
    }

    fn serialize_newtype_struct<T: Serialize + ?Sized>(
        self,
        _name: &'static str,
        value: &T,
    ) -> Result<Self::Ok, Error> {
        value.serialize(self)
    }

    fn serialize_newtype_variant<T: Serialize + ?Sized>(
        self,
        _name: &'static str,
        _index: u32,
        variant: &'static str,
        value: &T,
    ) -> Result<Self::Ok, Error> {
        let value = value.serialize(&mut *self)?;
        self.variant(variant, value)
    }

    fn serialize_seq(self, length: Option<usize>) -> Result<Self::SerializeSeq, Error> {
        Ok(Sequence {
            converter: self,
            items: Vec::with_capacity(length.unwrap_or(0)),
            variant: None,
        })
    }

    fn serialize_tuple(self, length: usize) -> Result<Self::SerializeTuple, Error> {
        self.serialize_seq(Some(length))
    }

    fn serialize_tuple_struct(
        self,
        _name: &'static str,
        length: usize,
    ) -> Result<Self::SerializeTupleStruct, Error> {
        self.serialize_seq(Some(length))
    }

    fn serialize_tuple_variant(
        self,
        _name: &'static str,
        _index: u32,
        variant: &'static str,
        length: usize,
    ) -> Result<Self::SerializeTupleVariant, Error> {
        Ok(Sequence {
            converter: self,
            items: Vec::with_capacity(length),
            variant: Some(variant),
        })
    }

    fn serialize_map(self, _length: Option<usize>) -> Result<Self::SerializeMap, Error> {
        let dict = PyDict::new(self.py);
        Ok(Map {
            converter: self,
            dict,
            key: None,
        })
    }

    fn serialize_struct(
        self,
        name: &'static str,
        _length: usize,
    ) -> Result<Self::SerializeStruct, Error> {
        let start = self.fields.len();
        Ok(Struct {
            converter: self,
            name,
            start,
            variant: None,
        })
    }

    fn serialize_struct_variant(
        self,
        name: &'static str,
        _index: u32,
        variant: &'static str,
        _length: usize,
    ) -> Result<Self::SerializeStructVariant, Error> {
        let start = self.fields.len();
        Ok(Struct {
            converter: self,
            name,
            start,
            variant: Some(variant),
        })
    }
}

/// A list under construction; with `variant`, it ends as `{variant: list}`.
pub(super) struct Sequence<'a, 'py> {
    converter: &'a mut Converter<'py>,
    items: Vec<Bound<'py, PyAny>>,
    variant: Option<&'static str>,
}

impl<'py> Sequence<'_, 'py> {
    fn push<T: Serialize + ?Sized>(&mut self, value: &T) -> Result<(), Error> {
        let item = value.serialize(&mut *self.converter)?;
        self.items.push(item);
        Ok(())
    }

    fn finish(self) -> Result<Bound<'py, PyAny>, Error> {
        let list = PyList::new(self.converter.py, self.items)?.into_any();
        match self.variant {
            Some(variant) => self.converter.variant(variant, list),
            None => Ok(list),
        }
    }
}

impl<'py> ser::SerializeSeq for Sequence<'_, 'py> {
    type Ok = Bound<'py, PyAny>;
    type Error = Error;

    fn serialize_element<T: Serialize + ?Sized>(&mut self, value: &T) -> Result<(), Error> {
        self.push(value)
    }

    fn end(self) -> Result<Self::Ok, Error> {
        self.finish()
    }
}

impl<'py> ser::SerializeTuple for Sequence<'_, 'py> {
    type Ok = Bound<'py, PyAny>;
    type Error = Error;

    fn serialize_element<T: Serialize + ?Sized>(&mut self, value: &T) -> Result<(), Error> {
        self.push(value)
    }

    fn end(self) -> Result<Self::Ok, Error> {
        self.finish()
    }
}

impl<'py> ser::SerializeTupleStruct for Sequence<'_, 'py> {
    type Ok = Bound<'py, PyAny>;
    type Error = Error;

    fn serialize_field<T: Serialize + ?Sized>(&mut self, value: &T) -> Result<(), Error> {
        self.push(value)
    }

    fn end(self) -> Result<Self::Ok, Error> {
        self.finish()
    }
}

impl<'py> ser::SerializeTupleVariant for Sequence<'_, 'py> {
    type Ok = Bound<'py, PyAny>;
    type Error = Error;

    fn serialize_field<T: Serialize + ?Sized>(&mut self, value: &T) -> Result<(), Error> {
        self.push(value)
    }

    fn end(self) -> Result<Self::Ok, Error> {
        self.finish()
    }
}

/// A dict under construction for a map.
pub(super) struct Map<'a, 'py> {
    converter: &'a mut Converter<'py>,
    dict: Bound<'py, PyDict>,
    key: Option<Bound<'py, PyString>>,
}

impl<'py> ser::SerializeMap for Map<'_, 'py> {
    type Ok = Bound<'py, PyAny>;
    type Error = Error;

    fn serialize_key<T: Serialize + ?Sized>(&mut self, key: &T) -> Result<(), Error> {
        self.key = Some(key.serialize(MapKey {
            converter: &mut *self.converter,
        })?);
        Ok(())
    }

    fn serialize_value<T: Serialize + ?Sized>(&mut self, value: &T) -> Result<(), Error> {
        let key = self
            .key
            .take()
            .ok_or_else(|| Error::Message("map value without a key".to_owned()))?;
        let value = value.serialize(&mut *self.converter)?;
        // A repeated key keeps its first position and takes the last value, as in
        // `json.loads`.
        self.dict.set_item(key, value)?;
        Ok(())
    }

    fn end(self) -> Result<Self::Ok, Error> {
        Ok(self.dict.into_any())
    }
}

/// A struct under construction; with `variant`, it ends as `{variant: dict}`.
pub(super) struct Struct<'a, 'py> {
    converter: &'a mut Converter<'py>,
    name: &'static str,
    start: usize,
    variant: Option<&'static str>,
}

impl<'py> Struct<'_, 'py> {
    fn push<T: Serialize + ?Sized>(&mut self, name: &'static str, value: &T) -> Result<(), Error> {
        let value = value.serialize(&mut *self.converter)?;
        self.converter.fields.push((name, value));
        Ok(())
    }

    fn finish(self) -> Result<Bound<'py, PyAny>, Error> {
        let dict = self
            .converter
            .finish_struct(self.name, self.start)?
            .into_any();
        match self.variant {
            Some(variant) => self.converter.variant(variant, dict),
            None => Ok(dict),
        }
    }
}

impl<'py> ser::SerializeStruct for Struct<'_, 'py> {
    type Ok = Bound<'py, PyAny>;
    type Error = Error;

    fn serialize_field<T: Serialize + ?Sized>(
        &mut self,
        name: &'static str,
        value: &T,
    ) -> Result<(), Error> {
        self.push(name, value)
    }

    fn end(self) -> Result<Self::Ok, Error> {
        self.finish()
    }
}

impl<'py> ser::SerializeStructVariant for Struct<'_, 'py> {
    type Ok = Bound<'py, PyAny>;
    type Error = Error;

    fn serialize_field<T: Serialize + ?Sized>(
        &mut self,
        name: &'static str,
        value: &T,
    ) -> Result<(), Error> {
        self.push(name, value)
    }

    fn end(self) -> Result<Self::Ok, Error> {
        self.finish()
    }
}

/// Turn a map key into the string `serde_json` writes for it, or refuse it as
/// `serde_json` does.
struct MapKey<'a, 'py> {
    converter: &'a mut Converter<'py>,
}

/// Return the text `serde_json` writes for a float.
fn written<T: Serialize>(value: T) -> Result<String, Error> {
    serde_json::to_string(&value).map_err(|error| Error::Message(error.to_string()))
}

impl<'py> ser::Serializer for MapKey<'_, 'py> {
    type Ok = Bound<'py, PyString>;
    type Error = Error;
    type SerializeSeq = ser::Impossible<Self::Ok, Error>;
    type SerializeTuple = ser::Impossible<Self::Ok, Error>;
    type SerializeTupleStruct = ser::Impossible<Self::Ok, Error>;
    type SerializeTupleVariant = ser::Impossible<Self::Ok, Error>;
    type SerializeMap = ser::Impossible<Self::Ok, Error>;
    type SerializeStruct = ser::Impossible<Self::Ok, Error>;
    type SerializeStructVariant = ser::Impossible<Self::Ok, Error>;

    fn serialize_bool(self, value: bool) -> Result<Self::Ok, Error> {
        Ok(self.converter.text(if value { "true" } else { "false" }))
    }

    fn serialize_i8(self, value: i8) -> Result<Self::Ok, Error> {
        Ok(self.converter.text(&value.to_string()))
    }

    fn serialize_i16(self, value: i16) -> Result<Self::Ok, Error> {
        Ok(self.converter.text(&value.to_string()))
    }

    fn serialize_i32(self, value: i32) -> Result<Self::Ok, Error> {
        Ok(self.converter.text(&value.to_string()))
    }

    fn serialize_i64(self, value: i64) -> Result<Self::Ok, Error> {
        Ok(self.converter.text(&value.to_string()))
    }

    fn serialize_i128(self, value: i128) -> Result<Self::Ok, Error> {
        Ok(self.converter.text(&value.to_string()))
    }

    fn serialize_u8(self, value: u8) -> Result<Self::Ok, Error> {
        Ok(self.converter.text(&value.to_string()))
    }

    fn serialize_u16(self, value: u16) -> Result<Self::Ok, Error> {
        Ok(self.converter.text(&value.to_string()))
    }

    fn serialize_u32(self, value: u32) -> Result<Self::Ok, Error> {
        Ok(self.converter.text(&value.to_string()))
    }

    fn serialize_u64(self, value: u64) -> Result<Self::Ok, Error> {
        Ok(self.converter.text(&value.to_string()))
    }

    fn serialize_u128(self, value: u128) -> Result<Self::Ok, Error> {
        Ok(self.converter.text(&value.to_string()))
    }

    fn serialize_f32(self, value: f32) -> Result<Self::Ok, Error> {
        if !value.is_finite() {
            return Err(Error::Message(
                "float key must be finite (got NaN or +/-inf)".to_owned(),
            ));
        }
        let text = written(value)?;
        Ok(self.converter.text(&text))
    }

    fn serialize_f64(self, value: f64) -> Result<Self::Ok, Error> {
        if !value.is_finite() {
            return Err(Error::Message(
                "float key must be finite (got NaN or +/-inf)".to_owned(),
            ));
        }
        let text = written(value)?;
        Ok(self.converter.text(&text))
    }

    fn serialize_char(self, value: char) -> Result<Self::Ok, Error> {
        Ok(self.converter.text(value.encode_utf8(&mut [0; 4])))
    }

    fn serialize_str(self, value: &str) -> Result<Self::Ok, Error> {
        Ok(self.converter.text(value))
    }

    fn serialize_bytes(self, _value: &[u8]) -> Result<Self::Ok, Error> {
        Err(key_must_be_a_string())
    }

    fn serialize_none(self) -> Result<Self::Ok, Error> {
        Err(key_must_be_a_string())
    }

    fn serialize_some<T: Serialize + ?Sized>(self, value: &T) -> Result<Self::Ok, Error> {
        value.serialize(self)
    }

    fn serialize_unit(self) -> Result<Self::Ok, Error> {
        Err(key_must_be_a_string())
    }

    fn serialize_unit_struct(self, _name: &'static str) -> Result<Self::Ok, Error> {
        Err(key_must_be_a_string())
    }

    fn serialize_unit_variant(
        self,
        _name: &'static str,
        _index: u32,
        variant: &'static str,
    ) -> Result<Self::Ok, Error> {
        Ok(self.converter.name(variant))
    }

    fn serialize_newtype_struct<T: Serialize + ?Sized>(
        self,
        _name: &'static str,
        value: &T,
    ) -> Result<Self::Ok, Error> {
        value.serialize(self)
    }

    fn serialize_newtype_variant<T: Serialize + ?Sized>(
        self,
        _name: &'static str,
        _index: u32,
        _variant: &'static str,
        _value: &T,
    ) -> Result<Self::Ok, Error> {
        Err(key_must_be_a_string())
    }

    fn serialize_seq(self, _length: Option<usize>) -> Result<Self::SerializeSeq, Error> {
        Err(key_must_be_a_string())
    }

    fn serialize_tuple(self, _length: usize) -> Result<Self::SerializeTuple, Error> {
        Err(key_must_be_a_string())
    }

    fn serialize_tuple_struct(
        self,
        _name: &'static str,
        _length: usize,
    ) -> Result<Self::SerializeTupleStruct, Error> {
        Err(key_must_be_a_string())
    }

    fn serialize_tuple_variant(
        self,
        _name: &'static str,
        _index: u32,
        _variant: &'static str,
        _length: usize,
    ) -> Result<Self::SerializeTupleVariant, Error> {
        Err(key_must_be_a_string())
    }

    fn serialize_map(self, _length: Option<usize>) -> Result<Self::SerializeMap, Error> {
        Err(key_must_be_a_string())
    }

    fn serialize_struct(
        self,
        _name: &'static str,
        _length: usize,
    ) -> Result<Self::SerializeStruct, Error> {
        Err(key_must_be_a_string())
    }

    fn serialize_struct_variant(
        self,
        _name: &'static str,
        _index: u32,
        _variant: &'static str,
        _length: usize,
    ) -> Result<Self::SerializeStructVariant, Error> {
        Err(key_must_be_a_string())
    }
}
