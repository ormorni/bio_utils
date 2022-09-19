use std::collections::hash_map::{Iter, IterMut, Keys, Values};
use fxhash::FxHashMap;

use crate::alignment::{AminoAcid, Sequence};


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
    pub fn get(&self, index: &str) -> Option<&Sequence> {
        self.data.get(index)
    }

    /// Gets an item in the alignment as a mutable string.
    pub fn get_mut(&mut self, index: &str) -> Option<&mut Sequence> {
        self.data.get_mut(index)
    }

    /// Inserts an item to the alignment.
    pub fn insert(&mut self, index: &str, value: &Sequence) {
        self.data.insert(String::from(index), value.clone());
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
        self.data.values().nth(0).map(|seq|seq.len()).unwrap_or(0)
    }

    pub fn iter(&self) -> Iter<String, Sequence> {
        self.data.iter()
    }
    pub fn iter_mut(&mut self) -> IterMut<String, Sequence> {
        self.data.iter_mut()
    }
    pub fn values(&self) -> Values<String, Sequence> {
        self.data.values()
    }
    pub fn keys(&self) -> Keys<String, Sequence> {
        self.data.keys()
    }
}

impl <'key, 'val> FromIterator<(&'key str, &'val Sequence)> for Alignment {
    fn from_iter<T: IntoIterator<Item=(&'key str, &'val Sequence)>>(iter: T) -> Self {
        let mut res = Alignment::new();
        for (name, seq) in iter {
            res.insert(name, seq);
        }
        res
    }
}
impl <'key, 'val> FromIterator<(&'key String, &'val Sequence)> for Alignment {
    fn from_iter<T: IntoIterator<Item=(&'key String, &'val Sequence)>>(iter: T) -> Self {
        let mut res = Alignment::new();
        for (name, seq) in iter {
            res.insert(name, seq);
        }
        res
    }
}