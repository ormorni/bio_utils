use rand::prelude::{IteratorRandom, StdRng};
use crate::alignment::AminoAcid;

/// A class to hold sequence data. Currently just a vector of amino acids.
pub type Sequence = Vec<AminoAcid>;

/// Generates a random sequence.
pub fn random_seq(length: usize, rng: &mut StdRng) -> Sequence {
    let mut seq = Sequence::new();
    for _ in 0..length {
        seq.push(*AminoAcid::iter().choose(rng).unwrap());
    }
    seq
}
