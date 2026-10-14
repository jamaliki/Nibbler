//! Check complete-chain crop eligibility directly in lossless A14P bytes.

use super::{Schema, residue::distance};

const LETTERS: &[u8; 20] = b"ARNDCQEGHILKMFPSTWYV";

pub(crate) struct Eligibility {
    pub eligible: Vec<u8>,
    pub complete_residues: Vec<i64>,
    pub max_chain_residues: Vec<i64>,
    pub reconstructed_chains: Vec<i64>,
}

struct Packed<'a, F> {
    byte: &'a F,
    samples: usize,
    sequence: usize,
    coordinates: usize,
    plane: usize,
}

impl<F: Fn(usize) -> u8> Packed<'_, F> {
    fn integer(&self, offset: usize) -> usize {
        u32::from_le_bytes(std::array::from_fn(|i| (self.byte)(offset + i))) as usize
    }

    fn bits(&self, row: usize, slot: usize, axis: usize) -> u32 {
        let index = self.coordinates + row * 42 + slot * 3 + axis;
        u32::from_le_bytes(std::array::from_fn(|plane| {
            (self.byte)(index + plane * self.plane)
        }))
    }

    fn finite(&self, row: usize, slot: usize, axis: usize) -> bool {
        let index = self.coordinates + row * 42 + slot * 3 + axis;
        (self.byte)(index + 3 * self.plane) & 0x7f != 0x7f
            || (self.byte)(index + 2 * self.plane) & 0x80 == 0
    }

    fn position(&self, row: usize, slot: usize) -> [f32; 3] {
        std::array::from_fn(|axis| f32::from_bits(self.bits(row, slot, axis)))
    }

    fn observed(&self, schema: &Schema, row: usize, sample: usize) -> Result<[bool; 14], String> {
        let letter = (self.byte)(self.sequence + row);
        let aa = LETTERS
            .iter()
            .position(|&value| value == letter)
            .ok_or_else(|| {
                format!("A14P sample {sample} residue {row} has a noncanonical identity")
            })?;
        let mut observed = [false; 14];
        for (slot, value) in observed.iter_mut().enumerate() {
            let physical = !schema.residues[aa].names[slot].is_empty();
            if !physical {
                if (0..3).any(|axis| self.bits(row, slot, axis) & 0x7fffffff != 0) {
                    return Err(format!(
                        "A14P sample {sample} residue {row} slot {slot} has invalid padding"
                    ));
                }
                continue;
            }
            let finite = [0, 1, 2].map(|axis| self.finite(row, slot, axis));
            if (0..3)
                .any(|axis| !finite[axis] && self.bits(row, slot, axis) & 0x7fffffff == 0x7f800000)
            {
                return Err(format!(
                    "A14P sample {sample} residue {row} slot {slot} has infinite coordinates"
                ));
            }
            if finite.iter().any(|&value| value) != finite.iter().all(|&value| value) {
                return Err(format!(
                    "A14P sample {sample} residue {row} slot {slot} has partially missing coordinates"
                ));
            }
            *value = finite[0];
        }
        if !observed[1] {
            return Err(format!(
                "A14P sample {sample} residue {row} lacks a CA atom"
            ));
        }
        Ok(observed)
    }

    fn sample(&self, schema: &Schema, sample: usize) -> Result<(i64, i64, i64), String> {
        let start = self.integer(16 + sample * 4);
        let end = self.integer(20 + sample * 4);
        let mut total = 0;
        let mut count = 0;
        let mut maximum = 0;
        let mut chains = 0;
        let mut previous = None;
        for row in start..end {
            let observed = self.observed(schema, row, sample)?;
            let position = self.position(row, 0);
            let continuous = previous.is_some_and(|c| {
                let value = distance(c, position);
                observed[0] && value > 0.8 && value < 1.8
            });
            if !continuous {
                maximum = maximum.max(count);
                count = 0;
                chains += 1;
            }
            let aa = LETTERS
                .iter()
                .position(|&value| value == (self.byte)(self.sequence + row))
                .ok_or("Missing canonical identity")?;
            let complete =
                (0..14).all(|slot| schema.residues[aa].names[slot].is_empty() || observed[slot]);
            count += i64::from(complete);
            total += i64::from(complete);
            previous = observed[2].then(|| self.position(row, 2));
        }
        maximum = maximum.max(count);
        Ok((total, maximum, chains))
    }
}

pub(crate) fn scan(
    schema: &Schema,
    length: usize,
    byte: impl Fn(usize) -> u8,
    min_chain: usize,
    crop_size: usize,
) -> Result<Eligibility, String> {
    if min_chain == 0 || crop_size < min_chain {
        return Err("AFDB eligibility needs a positive chain minimum within the crop size".into());
    }
    if length < 16 || [0, 1, 2, 3, 4, 5, 6, 7].map(&byte) != *b"A14P\x01\x20\x01\x00" {
        return Err("Unsupported A14P eligibility payload".into());
    }
    let integer = |offset| u32::from_le_bytes(std::array::from_fn(|i| byte(offset + i))) as usize;
    let samples = integer(8);
    let residues = integer(12);
    let sequence = 16usize
        .checked_add((samples + 1).checked_mul(4).ok_or("A14P size overflow")?)
        .ok_or("A14P size overflow")?;
    let expected = sequence
        .checked_add(residues.checked_mul(169).ok_or("A14P size overflow")?)
        .ok_or("A14P size overflow")?;
    if samples == 0 || residues < samples || length != expected {
        return Err("Invalid A14P eligibility payload length".into());
    }
    let input = Packed {
        byte: &byte,
        samples,
        sequence,
        coordinates: sequence + residues,
        plane: residues * 42,
    };
    if input.integer(16) != 0 || input.integer(16 + samples * 4) != residues {
        return Err("Invalid A14P residue offsets".into());
    }
    let mut result = Eligibility {
        eligible: Vec::with_capacity(samples),
        complete_residues: Vec::with_capacity(samples),
        max_chain_residues: Vec::with_capacity(samples),
        reconstructed_chains: Vec::with_capacity(samples),
    };
    for sample in 0..input.samples {
        let start = input.integer(16 + sample * 4);
        let end = input.integer(20 + sample * 4);
        if end <= start || end > residues {
            return Err(format!("Invalid A14P residue offsets at sample {sample}"));
        }
        if end - start > crop_size {
            return Err(format!(
                "A14P sample {sample} length {} exceeds crop_size {crop_size}; spatial-crop eligibility needs a separate plan",
                end - start
            ));
        }
        let (total, maximum, chains) = input.sample(schema, sample)?;
        result
            .eligible
            .push(u8::from(maximum as usize >= min_chain));
        result.complete_residues.push(total);
        result.max_chain_residues.push(maximum);
        result.reconstructed_chains.push(chains);
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::Packed;
    use proptest::prelude::*;

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(4096))]
        #[test]
        fn exponent_bytes_match_float32_isfinite(bits in any::<u32>()) {
            let raw = bits.to_le_bytes();
            let read = |i| raw[i];
            let view = Packed { byte: &read, samples: 1, sequence: 0, coordinates: 0, plane: 1 };
            prop_assert_eq!(view.finite(0, 0, 0), f32::from_bits(bits).is_finite());
        }
    }
}
