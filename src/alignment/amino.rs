use std::fmt::{Debug, Display, Formatter};

#[derive(Debug, Copy, Clone, Eq, PartialEq, Hash)]
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
            b'_' | b'*' | b'-' | b'.' => AminoAcid::GAP,

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