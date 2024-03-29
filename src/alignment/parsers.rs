use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Error, ErrorKind, Lines, Read, Write};
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

enum FastaIterState {
    Name,
    Sequence,
}

pub struct FastaIter {
    reader: Option<BufReader<File>>,
    buf: Vec<u8>,
}

/// An iterator iterating over the sequences in a FASTA file, without storing the whole file in memory.
/// Useful for large files.
impl FastaIter {
    /// An internal function possibly succeeding in making an iterator.
    fn make_reader(address: &Path) -> Option<FastaIter> {
        let file = File::open(address).ok()?;
        let mut temp = Vec::new();
        let mut reader = BufReader::new(file);
        reader.read_until(b'>', &mut temp).ok()?;
        Some(FastaIter {
            reader: Some(reader),
            buf: Vec::with_capacity(65536),
        })
    }

    pub fn new(address: &Path) -> FastaIter {
        if let Some(f_iter) = FastaIter::make_reader(address) {
            f_iter
        } else {
            FastaIter {
                reader: None,
                buf: Vec::new(),
            }
        }
    }
}

impl Iterator for FastaIter {
    type Item = (String, String);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(reader) = self.reader.as_mut() {
            self.buf.truncate(0);
            let name_length = reader.read_until(b'\n', &mut self.buf).unwrap();
            self.buf.pop();
            let name = String::from(String::from_utf8_lossy(&self.buf));
            self.buf.truncate(0);
            let seq_length = reader.read_until(b'>', &mut self.buf).unwrap();
            let next_exists = self.buf[seq_length - 1] == b'>';

            if next_exists {
                self.buf.truncate(seq_length - 1);
            }
            let seq = String::from(String::from_utf8_lossy(&self.buf));

            self.buf.truncate(0);
            if !next_exists {
                self.reader = None;
            }
            return Some((name, seq))
        }
        None
    }
}

pub fn string_to_seq(data: &[u8]) -> Sequence {
    let mut seq = Sequence::with_capacity(data.len());
    for aa_utf8 in data.iter().cloned() {
        if aa_utf8 != b'\n' {
            seq.push(AminoAcid::from_u8(aa_utf8));
        }
    }
    seq
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
    pub fn to_address(&self, address: &Path) -> std::io::Result<usize> {
        let file = File::create(address)?;

        let mut buf_writer = BufWriter::new(file);

        match AlignmentFormat::from_path(address) {
            AlignmentFormat::FASTA => buf_writer.write(self.to_fasta().as_bytes()),
            AlignmentFormat::STOCKHOLM => buf_writer.write(self.to_stockholm().as_bytes()),
            _ => Err(Error::new(ErrorKind::Other, "Unrecognized extension!")),
        }
    }

    /// Initializes an MSA from a string in the FASTA format.
    pub fn from_fasta(data: &str) -> Alignment {
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
    pub fn from_stockholm(data: &str) -> Alignment {
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
    use crate::utils::tests::SEED;

    const ALPHANUMERIC: &str = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_|0123456789";

    fn random_alignment(rng: &mut StdRng) -> Alignment {
        let mut ali = Alignment::new();

        for _ in 0..10 {
            let name = (0..10).map(|_| ALPHANUMERIC.chars().nth(rng.gen_range(0..ALPHANUMERIC.len())).unwrap()).collect::<String>();
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