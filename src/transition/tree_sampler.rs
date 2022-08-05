use super::defs::*;
use crate::alignment::Alignment;
use crate::alignment::AminoAcid;
use crate::phylogeny::tree::{Node, Tree};
use crate::transition::table::TransitionTable;
use std::collections::HashMap;
use itertools::{Itertools, all};
use lazy_static::lazy_static;
use rand::distributions::{Uniform, Distribution};
use rand::{Rng, RngCore};

/// A struct used to sample assignments of amino acids for a position on the phylogenetic tree.
pub struct TreeSampler<'t> {
    /// The phylogenetic tree sampled.
    pub tree: &'t Tree<String, f64>,
    /// The transition table used.
    pub table: TransitionTable,
    /// The probabilities of seeing each amino acid in each node, given only the amino acids below the noden.
    pub(crate) probs_given_below: HashMap<&'t Node<String, f64>, Vec21>,
}

fn choose_by_probs(v: &Vec21, rng: &mut impl Rng) -> usize {
    let mut uniform = rand_distr::Uniform::new(0., 1.);
    let mut chooser = uniform.sample(rng);
    for i in 0..20 {
        chooser -= v[i];
        if chooser < 0. {
            return i;
        }
    }
    20
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
                panic!("Leaf node has no sequence: {}!", node.data);
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

const EDGE_RESOLUTION: f64 = 0.01;
/// An accelerated tree sampling function.
pub fn fast_sample<'tree>(ali: &Alignment, pos: usize, tree: &'tree Tree<String, f64>, table: &TransitionTable, rng: &mut impl Rng) -> HashMap<&'tree Node<String, f64>, AminoAcid> {

    let node_count = tree.nodes().count();
    let base_id = tree.root.id;

    let mut cons = 0;
    let mut consensus_aa = vec![None; node_count];

    for node in tree.leaf_nodes() {
        consensus_aa[node.id - base_id] = Some(ali.get(&node.data).unwrap()[pos]);
    }

    for node in tree.postorder_iter() {
        if node.is_leaf() {
            continue;
        }
        let aa = consensus_aa[node.child_nodes[0].0.id - base_id];

        if all(node.child_nodes.iter(), |(child, _)|consensus_aa[child.id - base_id] == aa) {
            consensus_aa[node.id - base_id] = aa;
            cons += 1;
        }
    }

    let empty_table = table.get_table(0.);
    let mut table_sources = vec![];
    let mut tables = vec![&empty_table; node_count];

    let sorted: Vec<(f64, &Node<String, f64>)> = tree.edges()
        .filter(|edge|consensus_aa[edge.parent.id - base_id].is_none())
        .map(|edge|(*edge.data, edge.child))
        .sorted_by_key(|t| (t.0 * 100. / EDGE_RESOLUTION).round() as isize)
        .collect_vec();

    let mut start_idx = 0;
    let mut end_idx = 1;

    while start_idx < sorted.len() {
        while (end_idx < sorted.len()) && (sorted[end_idx].0 - sorted[start_idx].0 < 2. * EDGE_RESOLUTION) {
            end_idx += 1;
        }

        let mid_edge_length = (if end_idx < sorted.len() {sorted[end_idx].0} else {sorted[sorted.len() - 1].0} + sorted[start_idx].0) / 2.;
        table_sources.push(table.get_table(mid_edge_length));
        start_idx = end_idx;
        end_idx = start_idx + 1;
    }

    let mut start_idx = 0;
    let mut end_idx = 1;
    let mut table_ind = 0;

    while start_idx < sorted.len() {
        while (end_idx < sorted.len()) && (sorted[end_idx].0 - sorted[start_idx].0 < 2. * EDGE_RESOLUTION) {
            end_idx += 1;
        }
        for idx in start_idx..end_idx {
            tables[sorted[idx].1.id - base_id] = &table_sources[table_ind];
        }
        start_idx = end_idx;
        end_idx = start_idx + 1;
        table_ind += 1;
    }

    let mut probs: Vec<Vec21> = vec![Vec21::zeros(); tree.nodes().count()];

    // Initializing the probabilities.
    for node in tree.nodes() {
        if let Some(aa) = consensus_aa[node.id - base_id] {
            probs.get_mut(node.id - base_id).unwrap()[aa as usize] = 1.;
        } else {
            for i in 0..21 {
                probs.get_mut(node.id - base_id).unwrap()[i] = 1.;
            }
        }
    }

    // Updating the probabilities of parent nodes to take the subtree into account.
    for edge in tree.postorder_edge_iter() {
        // If the parent node is known, we skip the computation.
        if consensus_aa[edge.parent.id - base_id].is_some() {
            continue;
        }
        let child_probs = probs[edge.child.id - base_id];
        let evolved_probs: Vec21 = tables[edge.child.id - base_id] * child_probs;
        let parent_probs = &mut probs[edge.parent.id - base_id];
        for i in 0..21 {
            parent_probs[i] *= evolved_probs[i];
        }
        let tot = parent_probs.iter().sum::<f64>();
        for i in 0..21 {
            parent_probs[i] /= tot;
        }
    }

    // Assigning amino acids to all nodes.
    consensus_aa[0] = Some(AminoAcid::from_index(choose_by_probs(&probs[0], rng)));
    for edge in tree.preorder_edge_iter() {
        // If the child node already has an assignment, we don't have to calculate it anew.
        if consensus_aa[edge.child.id - base_id].is_some() {
            continue;
        }

        let choice = consensus_aa[edge.parent.id - base_id].unwrap();
        let table = &tables[edge.child.id - base_id];
        let col = table.column(choice as usize);

        let probs = &mut probs[edge.child.id - base_id];
        for i in 0..21 {
            probs[i] *= col[i];
        }
        let tot = probs.iter().sum::<f64>();
        for i in 0..21 {
            probs[i] /= tot;
        }

        consensus_aa[edge.child.id - base_id] = Some(AminoAcid::from_index(choose_by_probs(probs, rng)));
    }

    tree.nodes().map(|node|(node, consensus_aa[node.id - base_id].unwrap())).collect()
}
