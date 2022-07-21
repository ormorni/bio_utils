use super::defs::*;
use crate::alignment::Alignment;
use crate::alignment::AminoAcid;
use crate::phylogeny::tree::{Node, Tree};
use crate::transition::table::TransitionTable;
use std::collections::HashMap;

/// A struct used to sample assignments of amino acids for a position on the phylogenetic tree.
pub struct TreeSampler<'t> {
    /// The phylogenetic tree sampled.
    pub tree: &'t Tree<String, f64>,
    /// The transition table used.
    pub table: TransitionTable,
    /// The probabilities of seeing each amino acid in each node, given only the amino acids below the noden.
    pub(crate) probs_given_below: HashMap<&'t Node<String, f64>, Vec21>,
}

impl<'t> TreeSampler<'t> {
    /// Calculates the probability of the assignment of each amino acid to each node recursively.
    fn calculate_node_probs(&mut self, ali: &Alignment, index: usize, node: &'t Node<String, f64>) {
        // Initializing the probabilities of leaf nodes.
        if node.child_nodes.is_empty() {
            if let Some(seq) = ali.get(&node.data) {
                let mut res = zero_vec();
                res[seq[index] as usize] = 1.;
                self.probs_given_below.insert(node, res);
            } else {
                panic!("Leaf node has no sequence!");
            }
        } else {
            let mut probs = uniform();
            for (child_node, edge_length) in node.child_nodes.iter() {
                self.calculate_node_probs(ali, index, child_node);
                let child_probs = self.probs_given_below.get(child_node).unwrap();

                let evolved_probs = self.table.get_table(*edge_length) * child_probs;
                for i in 0..MAT_SIZE {
                    probs[i] *= evolved_probs[i];
                }
                probs /= probs.sum();
            }
            self.probs_given_below.insert(node, probs);
        }
    }

    /// Initializes a tree sampler from an index of an MSA and a phylogenetic tree.
    pub fn from_ali(
        ali: &Alignment,
        index: usize,
        tree: &'t Tree<String, f64>,
        table: &TransitionTable,
    ) -> TreeSampler<'t> {
        let mut sampler = TreeSampler {
            tree,
            table: table.clone(),
            probs_given_below: HashMap::new(),
        };
        sampler.calculate_node_probs(ali, index, &sampler.tree.root);
        sampler
    }

    /// Samples an assignment of amino acids
    pub fn sample(&self) -> HashMap<&'t Node<String, f64>, AminoAcid> {
        let mut res = HashMap::new();
        let mut parent_effect = HashMap::new();
        parent_effect.insert(&self.tree.root, uniform());

        for node in self.tree.preorder_iter() {
            // Calculating the probabilities.
            let parent_probs = parent_effect.get(node).unwrap();
            let mut probs = self.probs_given_below.get(node).unwrap().clone();

            for i in 0..MAT_SIZE {
                probs[i] *= parent_probs[i];
            }
            probs /= probs.sum();

            // Choosing the amino acid.
            let mut res_aa = AminoAcid::GAP;
            let mut chooser: f64 = rand::random();
            for i in 0..MAT_SIZE {
                chooser -= probs[i];
                if chooser <= 0. {
                    res_aa = AminoAcid::from_index(i as usize);
                    break;
                }
            }

            res.insert(node, res_aa);

            // Propagating the effect.
            for (child_node, dis) in node.child_nodes.iter() {
                let evolved_table = self.table.get_table(*dis);
                let prob_col = evolved_table.column(res_aa as usize);
                parent_effect.insert(child_node, Vec21::from(prob_col));
            }
        }

        res
    }
}
