use crate::transition::defs::Float;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{BufReader, Read};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

/// An implementation of a phylogenetic tree.
#[derive(Debug, Eq, PartialEq)]
pub struct Tree {
    pub root: Node,
}

static NODE_COUNTER: AtomicUsize = AtomicUsize::new(0);

/// A node in the phylogenetic tree.
#[derive(Debug, Clone, PartialEq)]
pub struct Node {
    pub name: String,
    pub child_nodes: Vec<(Node, Option<Float>)>,
    pub id: usize,
}

/// An edge in the phylogenetic tree.
/// Since the nodes are heavy, holds only a reference to the nodes.
#[derive(Debug, Copy, Clone)]
pub struct Edge<'t> {
    pub parent: &'t Node,
    pub child: &'t Node,
    pub length: Option<Float>,
}

impl <'t> Edge<'t> {
    pub fn new(parent: &'t Node, child: &'t Node, length: Option<Float>) -> Edge<'t> {
        Edge {parent, child, length}
    }
}

impl Tree {
    /// Initializing a new empty tree.
    pub fn new() -> Tree {
        Tree { root: Node::new() }
    }

    /// Initializing the phylogenetic tree from a path.
    pub fn from_address(address: &str) -> Tree {
        // Reading the data.
        let path = Path::new(address);
        let file = File::open(path).expect(&*format!("File not found: {}", address));
        let mut buf_reader = BufReader::new(file);
        let mut data = String::new();
        match buf_reader.read_to_string(&mut data) {
            Err(e) => panic!("Error in reading file: {:?}", e),
            _ => {}
        }

        return match path
            .extension()
            .expect(&*format!("No extension found for address {}", address))
            .to_str()
            .expect("Conversion of extension to str failed!")
        {
            "newick" | "treefile" => Tree::from_newick(&data),
            _ => panic!("Unrecognized format: {}", address),
        };
    }

    /// Parsing a slice containing a phylogenetic tree in newick format.
    fn parse_newick(data: &[u8]) -> Option<(Node, Option<Float>, usize)> {
        let mut res = Node::new();
        let mut dis: Option<Float> = None;
        let mut ind = 0;
        let mut next_ind;

        while data.get(ind)?.is_ascii_whitespace() {
            ind += 1;
        }

        // Parsing child node data.
        if *data.get(ind)? == b'(' {
            ind += 1;
            while *data.get(ind)? != b')' {
                let next_data = Tree::parse_newick(&data[ind..])?;
                res.child_nodes.push((next_data.0, next_data.1));
                ind += next_data.2;
            }
            assert_eq!(
                *data.get(ind).unwrap_or(&b'*'),
                b')',
                "No close-parenthesis at the end of the child list."
            );
            ind += 1;
        }
        // Checking for name.
        if !b":),;".contains(data.get(ind)?) {
            next_ind = ind;
            while !b":),;".contains(data.get(next_ind)?) {
                next_ind += 1;
            }
            res.name = String::from(String::from_utf8_lossy(&data[ind..next_ind]).trim());
            ind = next_ind;
        }
        // Checking for length.
        if *data.get(ind)? == b':' as u8 {
            next_ind = ind;
            while !b"),;".contains(data.get(next_ind)?) {
                next_ind += 1;
            }
            let trimmed_data =
                String::from(String::from_utf8_lossy(&data[ind + 1..next_ind]).trim());
            dis = trimmed_data.parse::<Float>().ok();
            ind = next_ind;
        }
        // Checking for termination.
        if *data.get(ind)? == b',' {
            ind += 1;
        }

        return Some((res, dis, ind));
    }

    /// Parses data in newick format to a phylogenetic tree.
    pub fn from_newick(data: &str) -> Tree {
        let (node, dis, _) = Tree::parse_newick(data.trim().as_bytes()).unwrap();
        assert!(dis.is_none(), "Distance set on root node!");
        return Tree { root: node };
    }
}

impl Node {
    /// Initializes an empty node.
    pub fn new() -> Node {
        Node {
            name: String::new(),
            child_nodes: vec![],
            id: NODE_COUNTER.fetch_add(1, Ordering::SeqCst),
        }
    }
    /// Initializes an empty node with a given name.
    pub fn from_name(name: String) -> Node {
        Node {
            name,
            child_nodes: vec![],
            id: NODE_COUNTER.fetch_add(1, Ordering::SeqCst),
        }
    }
    /// Returns if the node has no child nodes.
    pub fn is_leaf(&self) -> bool {
        self.child_nodes.is_empty()
    }
}

impl Hash for Node {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.id.hash(state);
    }
}

impl std::cmp::Eq for Node {}

#[cfg(test)]
mod tests {
    use crate::phylogeny::tree::{Node, Tree};

    /// Testing that the tree initialization from newick strings works as expected, and that
    #[test]
    fn test_newick_initialization() {
        let mut cons_tree = Tree::new();
        cons_tree.root.name = String::from("root");
        cons_tree.root.child_nodes.push((Node::new(), Some(1.)));
        cons_tree.root.child_nodes.push((Node::new(), Some(0.1)));
        cons_tree.root.child_nodes.push((Node::new(), None));
        cons_tree.root.child_nodes[0].0.name = String::from("A");
        cons_tree.root.child_nodes[1].0.name = String::from("B");
        cons_tree.root.child_nodes[2].0.name = String::from("C");

        let data = String::from("  (A : 1 , B : .1 ,  C ) root ;  \n");

        let parsed_tree = Tree::from_newick(&data);

        assert_eq!(&cons_tree.root.name, &parsed_tree.root.name);
        assert_eq!(parsed_tree.root.child_nodes.len(), 3);

        for i in 0..3 {
            assert_eq!(
                cons_tree.root.child_nodes[i].0.name,
                parsed_tree.root.child_nodes[i].0.name
            );
            assert_eq!(
                cons_tree.root.child_nodes[i].1,
                parsed_tree.root.child_nodes[i].1
            );
        }
    }
}

impl<'t> PartialEq<Self> for Edge<'t> {
    fn eq(&self, other: &Self) -> bool {
        self.parent == other.parent && self.child == other.child && self.length == other.length
    }
}

impl<'t> Eq for Edge<'t> {}

impl<'t> Hash for Edge<'t> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.parent.hash(state);
        self.child.hash(state);
    }
}
