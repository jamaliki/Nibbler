//! A fast hasher for maps keyed by short text and indices.
//!
//! The standard SipHash resists collision attacks, which these maps do not
//! need: their keys come from one parsed structure and never leave the call.

use std::collections::{HashMap, HashSet};
use std::hash::{BuildHasherDefault, Hasher};

const SEED: u64 = 0x51_7c_c1_b7_27_22_0a_95;

/// The word-at-a-time multiply-rotate hash used by rustc.
#[derive(Clone, Copy, Default)]
pub(crate) struct Fx(u64);

impl Fx {
    #[inline]
    fn add(&mut self, word: u64) {
        self.0 = (self.0.rotate_left(5) ^ word).wrapping_mul(SEED);
    }
}

impl Hasher for Fx {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        let mut chunks = bytes.chunks_exact(8);
        for chunk in &mut chunks {
            let mut word = [0; 8];
            word.copy_from_slice(chunk);
            self.add(u64::from_le_bytes(word));
        }
        let rest = chunks.remainder();
        if !rest.is_empty() {
            let mut word = [0; 8];
            word[..rest.len()].copy_from_slice(rest);
            self.add(u64::from_le_bytes(word));
        }
    }

    #[inline]
    fn write_u8(&mut self, value: u8) {
        self.add(u64::from(value));
    }

    #[inline]
    fn write_u32(&mut self, value: u32) {
        self.add(u64::from(value));
    }

    #[inline]
    fn write_u64(&mut self, value: u64) {
        self.add(value);
    }

    #[inline]
    fn write_usize(&mut self, value: usize) {
        self.add(value as u64);
    }

    #[inline]
    fn finish(&self) -> u64 {
        // The table picks buckets from the low bits, which the multiply leaves
        // dependent on the low input bytes alone.
        let mixed = (self.0 ^ (self.0 >> 32)).wrapping_mul(0x9e37_79b9_7f4a_7c15);
        mixed ^ (mixed >> 29)
    }
}

/// A `HashMap` with the fast hasher.
pub(crate) type FastMap<K, V> = HashMap<K, V, BuildHasherDefault<Fx>>;
/// A `HashSet` with the fast hasher.
pub(crate) type FastSet<K> = HashSet<K, BuildHasherDefault<Fx>>;
