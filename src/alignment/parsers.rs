use std::fs::File;
use std::io::{BufReader, BufWriter, Error, ErrorKind, Read, Write};
use std::path::Path;
use itertools::Itertools;
use crate::alignment::{Alignment, AminoAcid, Sequence};

const STOCKHOLM_HEADER: &str = "# STOCKHOLM 1.0\n";
const STOCKHOLM_FOOTER: &str = "//\n";

#[derive(Eq, PartialEq, Hash, Copy, Clone, Debug)]
pub enum AlignmentFormat {
    FASTA,
    STOCKHOLM,
    OTHER
}

impl AlignmentFormat {
    /// Initializes the alignment format from an extension.
    pub fn from_ext(ext: &str) -> AlignmentFormat {
        match ext {
            "fa" | "fasta" | "afa" | "raptorx" => AlignmentFormat::FASTA,
            "sto" | "stockholm" | "stk" => AlignmentFormat::STOCKHOLM,
            _ => AlignmentFormat::OTHER,
        }
    }

    /// Initializes the alignment format from a path, using its extension.
    pub fn from_path(path: &Path) -> AlignmentFormat {
        AlignmentFormat::from_ext(path.extension()
            .map(|s|s.to_str().unwrap_or(""))
            .unwrap_or(""))
    }
    /// Gets an extension fitting the format.
    pub fn get_ext(&self) -> &str {
        match self {
            AlignmentFormat::FASTA => "fa",
            AlignmentFormat::STOCKHOLM => "sto",
            AlignmentFormat::OTHER => panic!(),
        }
    }
}

/// Implementations for MSA reading and writing from various file formats.
impl Alignment {
    /// Reads an alignment from a given address.
    pub fn from_address(address: &Path) -> std::io::Result<Alignment> {
        let file = File::open(address)?;
        let mut buf_reader = BufReader::new(file);
        let mut data = String::new();

        buf_reader.read_to_string(&mut data)?;

        match AlignmentFormat::from_path(address) {
            AlignmentFormat::FASTA => Ok(Alignment::from_fasta(&data)),
            AlignmentFormat::STOCKHOLM => Ok(Alignment::from_stockholm(&data)),
            _ => Err(Error::new(ErrorKind::Other, "Unrecognized extension!")),
        }
    }

    /// Reads an alignment from a given address.
    pub fn to_address(&self, address: &str) -> std::io::Result<usize> {
        let path = Path::new(address);
        let file = File::create(path)?;

        let mut buf_writer = BufWriter::new(file);

        match AlignmentFormat::from_path(path) {
            AlignmentFormat::FASTA => buf_writer.write(self.to_fasta().as_bytes()),
            AlignmentFormat::STOCKHOLM => buf_writer.write(self.to_stockholm().as_bytes()),
            _ => Err(Error::new(ErrorKind::Other, "Unrecognized extension!")),
        }
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
            if ln.starts_with(">") {
                let name = ln[1..].split_whitespace().next().unwrap();
                last_key = Some(String::from(name));
                res.data.insert(String::from(name), Sequence::new());
            } else if let Some(key) = last_key.as_ref() {
                res.get_mut(key)
                    .expect(&format!("Last key not matching any key: {}!", key))
                    .extend(ln.chars().map(AminoAcid::from_char));
            }
        }

        res
    }

    /// Encodes the alignment in the FASTA format.
    pub fn to_fasta(&self) -> String {
        let mut res = String::new();
        for (name, seq) in self.data.iter() {
            res.push('>');
            res.extend(name.chars());
            res.push('\n');
            res.extend(seq.iter().map(|aa|aa.string()));
            res.push('\n');
        }
        res
    }

    /// Initializes an MSA from a string in the FASTA format.
    pub fn from_stockholm(data: &String) -> Alignment {
        let mut res = Alignment::new();

        for raw_ln in data.lines() {
            let ln = raw_ln.trim();
            // Empty lines.
            if ln.len() == 0 {
                continue;
            }
            // Metadata lines.
            if ln.starts_with("#") || ln.starts_with("/") {
                continue;
            }
            let parts = ln.split_whitespace().collect_vec();
            assert_eq!(parts.len(), 2);

            let name = String::from(parts[0]);
            let seq = parts[1];
            res.data.entry(name).or_default().extend(seq.chars().map(AminoAcid::from_char));
        }

        res
    }

    /// Encodes the alignment in the Stockholm format.
    pub fn to_stockholm(&self) -> String {
        let mut res = String::new();
        res.extend(STOCKHOLM_HEADER.chars());
        for (name, seq) in self.data.iter() {
            res.extend(name.chars());
            res.push(' ');
            res.extend(seq.iter().map(|aa|aa.string()));
            res.push('\n');
        }
        res.extend(STOCKHOLM_FOOTER.chars());
        res
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;
    use rand::prelude::{IteratorRandom, Rng, SeedableRng, StdRng};
    use crate::alignment::{Alignment, AminoAcid};
    use crate::test_utils::SEED;

    const ALNUM: &str = "abcdefghijklmnopqrstuvwxysABCDEFGHIJKLMNOPQRSTUVWXYZ_|0123456789";

    fn random_alignment(rng: &mut StdRng) -> Alignment {
        let mut ali = Alignment::new();

        for _ in 0..10 {
            let name = (0..10).map(|_|ALNUM.chars().nth(rng.gen_range(0..ALNUM.len())).unwrap()).collect::<String>();
            let seq = (0..50).map(|_|AminoAcid::iter().choose(rng).unwrap().clone()).collect_vec();
            ali.insert(&name, &seq);
        }

        ali
    }

    #[test]
    fn test_fasta() {
        let mut rng = StdRng::from_seed(SEED);

        for _ in 0..50 {
            let ali = random_alignment(&mut rng);
            let encoded = ali.to_fasta();
            let decoded = Alignment::from_fasta(&encoded);
            assert_eq!(ali, decoded);
        }
    }

    #[test]
    fn test_stockholm() {
        let mut rng = StdRng::from_seed(SEED);

        for _ in 0..50 {
            let ali = random_alignment(&mut rng);
            let encoded = ali.to_stockholm();
            let decoded = Alignment::from_stockholm(&encoded);
            assert_eq!(ali, decoded);
        }
    }
}