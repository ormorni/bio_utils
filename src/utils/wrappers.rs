use itertools::Itertools;
use regex::Regex;
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use crate::alignment::Alignment;
use crate::alignment::parsers::AlignmentFormat;
use crate::phylogeny::tree::Tree;

/// A trait for objects that hold an Alignment that can be represented as a file.
pub trait AliPath {
    /// Gets the path of the file that stores the alignment.
    fn get_path(&self, format: AlignmentFormat) -> PathBuf;
}

/// A path to an alignment automatically satisfies this.
/// However, it may store the alignment in the wrong format, and would have to be converted.
impl AliPath for &Path {
    fn get_path(&self, format: AlignmentFormat) -> PathBuf {
        if AlignmentFormat::from_path(self) == format {
            self.to_path_buf()
        } else {
            (&Alignment::from_address(self).unwrap()).get_path(format)
        }
    }
}

/// A path to an alignment automatically satisfies this.
/// However, it may store the alignment in the wrong format, and would have to be converted.
impl AliPath for &PathBuf {
    fn get_path(&self, format: AlignmentFormat) -> PathBuf {
        if AlignmentFormat::from_path(self) == format {
            PathBuf::from(self)
        } else {
            (&Alignment::from_address(self).unwrap()).get_path(format)
        }
    }
}

/// To represent an alignment as a file, the alignment is written to a temporary file.
impl AliPath for Alignment {
    fn get_path(&self, format: AlignmentFormat) -> PathBuf {
        let mut temp_dir = std::env::temp_dir();
        let uuid = uuid::Uuid::new_v4();
        temp_dir.push(format!("{}.{}", uuid, format.get_ext()));
        self.to_address(&temp_dir).unwrap();
        temp_dir
    }
}

impl AliPath for &Alignment {
    fn get_path(&self, format: AlignmentFormat) -> PathBuf {
        (*self).get_path(format)
    }
}

/// Runs a single linux command.
/// On Windows environments, the command uses WSL to run linux.
pub fn run_cmd(cmd: &str) -> Result<Output, std::io::Error> {
    if cfg!(windows) {
        let cmd = cmd.replace(">", "^>");
        let cmd = cmd.split_whitespace().map(|part| {
            if part.starts_with("C:") {
                part.replace("C:", "/mnt/c").replace("\\", "/")
            } else {
                part.replace("\\", "/")
            }
        });

        Command::new("wsl").args(cmd).output()
    } else {
        let mut args = cmd.split_whitespace();
        let name = args.next().expect("run_cmd called with no argument!");
        Command::new(name).args(args).output()
    }
}

/// A struct holding the results of an HmmSearch query.
#[derive(Copy, Clone, Debug)]
pub struct HmmerResult {
    /// The start position of the match in the query HMM.
    pub q_start: usize,
    /// The end position of the match in the query HMM.
    pub q_end: usize,
    /// The start position of the match in the target sequence.
    pub t_start: usize,
    /// The end position of the match in the target sequence.
    pub t_end: usize,
    /// The E-Value of the match.
    pub e_value: f64,
}

/// A wrapper around HMMER's HmmSearch method.
pub fn hmmsearch(hmm: &Path, database: &Path) -> HashMap<String, Vec<HmmerResult>> {
    assert!(hmm.exists(), "Query HMM doesn't exist: {:?}", hmm);
    assert!(database.exists(), "Database doesn't exist: {:?}", database);

    let output = run_cmd(&format!(
        "hmmsearch {} {}",
        hmm.to_str().unwrap(),
        database.to_str().unwrap()
    ))
    .expect("HmmSearch command failed!");

    let mut res = HashMap::new();

    let mut curr_name = String::new();

    let match_pattern: Regex = Regex::new(r"^\s*\d+\s+[?!].*").unwrap();

    for row in String::from_utf8(output.stdout)
        .expect("Failed to parse HmmSearch output!")
        .lines()
    {
        if row.starts_with(">>") {
            // A row starting the matches of a new protein.
            curr_name = String::from(row.split_whitespace().skip(1).next().unwrap());
            res.insert(curr_name.clone(), Vec::new());
        } else if match_pattern.is_match(row) {
            // A row containing the data on the match of the HMM to a sequence in the database.
            let parts = row.split_whitespace().collect_vec();

            let match_ = HmmerResult {
                q_start: parts[6].parse::<usize>().unwrap(),
                q_end: parts[7].parse::<usize>().unwrap(),
                t_start: parts[9].parse::<usize>().unwrap(),
                t_end: parts[10].parse::<usize>().unwrap(),
                e_value: parts[4].parse::<f64>().unwrap(),
            };

            res.entry(curr_name.clone()).or_default().push(match_);
        }
    }

    res
}

/// A wrapper around MAFFT.
/// Generates an MSA from an alignment.
pub fn mafft<AliPathRef: AliPath>(ali: AliPathRef) -> Alignment {
    let path_str = ali.get_path(AlignmentFormat::FASTA);
    let ali = Path::new(&path_str);
    assert!(ali.exists(), "Alignment doesn't exist: {:?}", ali);

    let res = run_cmd(&format!("mafft --quiet {}", ali.to_str().unwrap()))
        .expect("MAFFT command failed!");
    let out = String::from_utf8(res.stdout).unwrap();
    Alignment::from_fasta(&out)
}

/// A wrapper around HmmBuild.
/// Generates an HMM from an MSA in the stockholm format.
pub fn hmmbuild<AliPathRef: AliPath>(msa: AliPathRef, hmm_address: &Path) {
    let msa = msa.get_path(AlignmentFormat::STOCKHOLM);
    assert!(msa.exists(), "Alignment doesn't exist: {:?}", msa);

    run_cmd(&format!(
        "hmmbuild {} {}",
        hmm_address.to_str().unwrap(),
        msa.to_str().unwrap()
    ))
    .expect("HmmBuild command failed!");
}

/// Uses the CD-HIT program to cluster sequences.
pub fn cd_hit<AliPathRef: AliPath>(ali: AliPathRef, tar_address: &Path, cluster_cutoff: f64) {
    let ali = ali.get_path(AlignmentFormat::FASTA);
    assert!(
        Path::new(&ali).exists(),
        "Alignment doesn't exist: {:?}",
        ali
    );

    let word_size = match cluster_cutoff {
        _ if (0.0..0.4).contains(&cluster_cutoff) => {
            panic!("CD-HIT can't be used with cutoff={}", cluster_cutoff)
        }
        _ if (0.4..0.6667).contains(&cluster_cutoff) => 2,
        _ if (0.6667..0.75).contains(&cluster_cutoff) => 3,
        _ if (0.75..0.8).contains(&cluster_cutoff) => 4,
        _ => 5,
    };

    run_cmd(&format!(
        "cd-hit -i {} -o {} -c {:.2} -n {}",
        ali.to_str().unwrap(),
        tar_address.to_str().unwrap(),
        cluster_cutoff,
        word_size,
    ))
    .expect("CD-HIT command failed!");
}

/// Uses the FastTree program to compute a phylogenetic tree for the given multiple sequence alignment.
pub fn fasttree<AliPathRef: AliPath>(ali: AliPathRef) -> Tree<String, f64> {
    let path = ali.get_path(AlignmentFormat::FASTA);

    let out =
        run_cmd(&format!("fasttree {}", path.to_str().unwrap())).expect("Failed to run FastTree!");

    let tree_data = String::from_utf8(out.stdout).expect("Failed to parse tree!");
    if tree_data.is_empty() {
        // An error has occurred.
        eprintln!("Error: {}", String::from_utf8(out.stderr).expect("Failed to parse error message!"));
        panic!("Error in FastTree!");
    }
    Tree::<String, f64>::from_newick(&tree_data)
}

pub fn main() {
    // // Testing HmmSearch
    // let res = hmmsearch(Path::new("data/theme1ex2.hmm"), Path::new("data/pdb.fa"));
    // for (key, val) in res.iter() {
    //     println!("{}: {}", key, val.len());
    // }

    // hmmbuild(Path::new("data/1bxw1.afa"), Path::new("data/1bxw1_test.hmm"));
    // let res = hmmsearch(Path::new("data/1bxw1_test.hmm"), Path::new("data/pdb.fa"));
    // for name in res.keys() {
    //     println!("{}", name);
    // }

    // let ali = Alignment::from_address(Path::new("data/1bxw1.fa")).unwrap();
    // ali.to_address(Path::new("data/1bxw1.sto"));
    // cd_hit(
    //     Path::new("data/1bxw1.sto"),
    //     Path::new("data/1bxw1_90.fa"),
    //     0.5,
    // );
}
