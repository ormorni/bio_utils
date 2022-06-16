use fxhash::FxHashMap;
use rand::prelude::*;
use std::fmt::{Display, Formatter};
use std::fs::File;
use std::io::{BufReader, Error, ErrorKind, Read};
use std::path::Path;

const EMPTY_EXT: String = String::new();

/// A class used to store multiple sequence alignments.
#[derive(Debug, Clone)]
pub struct Alignment {
    pub data: FxHashMap<String, Sequence>,
}

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

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum AminoAcid {
    A = 0,
    R = 1,
    N = 2,
    D = 3,
    C = 4,
    Q = 5,
    E = 6,
    G = 7,
    H = 8,
    I = 9,
    L = 10,
    K = 11,
    M = 12,
    F = 13,
    P = 14,
    S = 15,
    T = 16,
    W = 17,
    Y = 18,
    V = 19,
    GAP = 20,
}

impl AminoAcid {
    pub fn from_u8(letter: u8) -> AminoAcid {
        match letter.to_ascii_uppercase() {
            b'A' => AminoAcid::A,
            b'R' => AminoAcid::R,
            b'N' => AminoAcid::N,
            b'D' => AminoAcid::D,
            b'C' => AminoAcid::C,
            b'E' => AminoAcid::E,
            b'Q' => AminoAcid::Q,
            b'G' => AminoAcid::G,
            b'H' => AminoAcid::H,
            b'I' => AminoAcid::I,
            b'K' => AminoAcid::K,
            b'L' => AminoAcid::L,
            b'M' => AminoAcid::M,
            b'F' => AminoAcid::F,
            b'P' => AminoAcid::P,
            b'S' => AminoAcid::S,
            b'T' => AminoAcid::T,
            b'W' => AminoAcid::W,
            b'Y' => AminoAcid::Y,
            b'V' => AminoAcid::V,
            b'_' | b'*' | b'-' => AminoAcid::GAP,

            // Odd codes.
            b'B' => AminoAcid::N, // Aspargine or Aspartic Acid
            b'Z' => AminoAcid::Q, // Glutamine or Glutamic Acid
            b'U' => AminoAcid::S, // Selenocysteine
            b'X' => AminoAcid::A, // Unrecognized amino acid

            _ => panic!(
                "Unrecognized amino acid: {}",
                String::from_utf8_lossy(&[letter])
            ),
        }
    }

    pub fn from_char(c: char) -> AminoAcid {
        let mut buf = [b'X'];
        c.encode_utf8(&mut buf);
        AminoAcid::from_u8(buf[0])
    }

    pub fn string(&self) -> &str {
        match self {
            AminoAcid::A => "A",
            AminoAcid::R => "R",
            AminoAcid::N => "N",
            AminoAcid::D => "D",
            AminoAcid::C => "C",
            AminoAcid::E => "E",
            AminoAcid::Q => "Q",
            AminoAcid::G => "G",
            AminoAcid::H => "H",
            AminoAcid::I => "I",
            AminoAcid::K => "K",
            AminoAcid::L => "L",
            AminoAcid::M => "M",
            AminoAcid::F => "F",
            AminoAcid::P => "P",
            AminoAcid::S => "S",
            AminoAcid::T => "T",
            AminoAcid::W => "W",
            AminoAcid::Y => "Y",
            AminoAcid::V => "V",
            AminoAcid::GAP => "-",
        }
    }

    pub fn from_index(index: usize) -> AminoAcid {
        match index {
            0 => AminoAcid::A,
            1 => AminoAcid::R,
            2 => AminoAcid::N,
            3 => AminoAcid::D,
            4 => AminoAcid::C,
            5 => AminoAcid::Q,
            6 => AminoAcid::E,
            7 => AminoAcid::G,
            8 => AminoAcid::H,
            9 => AminoAcid::I,
            10 => AminoAcid::L,
            11 => AminoAcid::K,
            12 => AminoAcid::M,
            13 => AminoAcid::F,
            14 => AminoAcid::P,
            15 => AminoAcid::S,
            16 => AminoAcid::T,
            17 => AminoAcid::W,
            18 => AminoAcid::Y,
            19 => AminoAcid::V,
            20 => AminoAcid::GAP,
            _ => panic!("Illegal amino acid index: {}", index),
        }
    }

    /// Returns an iterator over all amino acids.
    pub fn iter() -> impl Iterator<Item = &'static AminoAcid> + Clone {
        [
            AminoAcid::A,
            AminoAcid::R,
            AminoAcid::N,
            AminoAcid::D,
            AminoAcid::C,
            AminoAcid::Q,
            AminoAcid::E,
            AminoAcid::G,
            AminoAcid::H,
            AminoAcid::I,
            AminoAcid::L,
            AminoAcid::K,
            AminoAcid::M,
            AminoAcid::F,
            AminoAcid::P,
            AminoAcid::S,
            AminoAcid::T,
            AminoAcid::W,
            AminoAcid::Y,
            AminoAcid::V,
            AminoAcid::GAP,
        ]
        .iter()
    }
}

impl Display for AminoAcid {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.string())
    }
}

impl Alignment {
    /// Initializes a new empty alignment.
    pub fn new() -> Alignment {
        Alignment {
            data: FxHashMap::default(),
        }
    }

    /// Reads an alignment from a given address.
    pub fn from_address(address: &str) -> std::io::Result<Alignment> {
        let path = Path::new(address);
        let file = File::open(path)?;
        let mut buf_reader = BufReader::new(file);
        let mut data = String::new();

        buf_reader.read_to_string(&mut data)?;

        return match path
            .extension()
            .map(|os| os.to_str())
            .unwrap_or(Some(&EMPTY_EXT))
            .unwrap_or(&EMPTY_EXT)
        {
            "fa" | "fasta" | "afa" => Ok(Alignment::from_fasta(&data)),
            _ => Err(Error::new(ErrorKind::Other, "Unrecognized extension!")),
        };
    }

    /// Initializes an MSA from a string in the FASTA format.
    pub fn from_fasta(data: &String) -> Alignment {
        let mut res = Alignment::new();

        let mut last_key: Option<String> = None;

        for raw_ln in data.lines() {
            let ln = raw_ln.trim();
            if ln.len() == 0 {
                continue;
            }
            if &ln[0..1] == ">" {
                last_key = Some(String::from(&ln[1..]));
                res.data.insert(String::from(&ln[1..]), Sequence::new());
            } else if let Some(key) = last_key.as_ref() {
                res.get_mut(key)
                    .expect("Last key not matching any key!")
                    .extend(ln.chars().map(AminoAcid::from_char));
            }
        }

        res
    }

    /// Gets an item in the alignment.
    pub fn get(&self, index: &String) -> Option<&Vec<AminoAcid>> {
        self.data.get(index)
    }

    /// Gets an item in the alignment as a mutable string.
    pub fn get_mut(&mut self, index: &String) -> Option<&mut Vec<AminoAcid>> {
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
