use fxhash::FxHashMap;

use crate::alignment::amino::AminoAcid;
use crate::alignment::Sequence;


/// A class used to store multiple sequence alignments.
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct Alignment {
    pub data: FxHashMap<String, Sequence>,
}


impl Alignment {
    /// Initializes a new empty alignment.
    pub fn new() -> Alignment {
        Alignment {
            data: FxHashMap::default(),
        }
    }



    /// Gets an item in the alignment.
    pub fn get(&self, index: &str) -> Option<&Vec<AminoAcid>> {
        self.data.get(index)
    }

    /// Gets an item in the alignment as a mutable string.
    pub fn get_mut(&mut self, index: &str) -> Option<&mut Vec<AminoAcid>> {
        self.data.get_mut(index)
    }

    /// Inserts an item to the alignment.
    pub fn insert(&mut self, index: &String, value: &Vec<AminoAcid>) {
        self.data.insert(index.clone(), value.clone());
    }

    /// Trims columns from the alignment.
    pub fn trim<F>(&self, f: F) -> Alignment
    where
        F: Fn(usize) -> bool,
    {
        let mut res = Alignment::new();

        let mut good_indices: Vec<bool> = Vec::new();

        for (name, seq) in self.data.iter() {
            while good_indices.len() < seq.len() {
                good_indices.push(f(good_indices.len()));
            }

            let mut trimmed_seq = Sequence::new();
            for (good_ind, aa) in good_indices.iter().zip(seq.iter()) {
                if *good_ind {
                    trimmed_seq.push(*aa);
                }
            }

            res.insert(&name, &trimmed_seq);
        }

        res
    }

    /// Returns the number of sequences in the alignment.
    pub fn size(&self) -> usize {
        self.data.len()
    }

    /// Returns the length of the sequences in the alignment.
    /// Assumes the alignment is an MSA: all sequences are of the same length.
    pub fn len(&self) -> usize {
        self.data.values().nth(0).unwrap().len()
    }
}

#[cfg(test)]
mod tests {
    use crate::alignment::AminoAcid;

    /// Testing that the amino acid int conversions work properly.
    #[test]
    fn test_amino_acid_int() {
        for i in 0..21 {
            assert_eq!(i as usize, AminoAcid::from_index(i) as usize);
        }
    }

    #[test]
    fn test_amino_acid_char() {
        for (i, c) in "ARNDCQEGHILKMFPSTWYV-".chars().enumerate() {
            assert_eq!(AminoAcid::from_char(c), AminoAcid::from_index(i as usize));
        }
    }
}
