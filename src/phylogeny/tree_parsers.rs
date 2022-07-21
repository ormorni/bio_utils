use std::fmt::{Display, Write as FmtWrite};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;
use crate::phylogeny::tree::{Node, Tree};

/// Reads the data in the path to a string.
fn read_path(address: &Path) -> String {
    assert!(address.exists(), "Attempted to read nonexistent file: {}", address.to_str().unwrap());
    let file = File::open(address).expect(&format!("File not found: {}", address.to_str().unwrap()));
    let mut buf_reader = BufReader::new(file);
    let mut data = String::new();

    match buf_reader.read_to_string(&mut data) {
        Err(e) => panic!("Error in reading file: {:?}", e),
        _ => {}
    };

    data
}

fn write_path(address: &Path, data: &str) {
    let file = File::open(address).expect(&format!("Failed to open file: {}", address.to_str().unwrap()));
    let mut buf_writer = BufWriter::new(file);
    assert_eq!(buf_writer.write(data.as_bytes()).unwrap(), data.len());
}

/// Parsing a slice containing a phylogenetic tree in newick format.
fn parse_newick(data: &[u8]) -> Option<(Node<String, Option<f64>>, Option<f64>, usize)> {
    let mut res = Node::new();
    let mut dis: Option<f64> = None;
    let mut ind = 0;
    let mut next_ind;

    while data.get(ind)?.is_ascii_whitespace() {
        ind += 1;
    }

    // Parsing child node data.
    if *data.get(ind)? == b'(' {
        ind += 1;
        while *data.get(ind)? != b')' {
            let next_data = parse_newick(&data[ind..])?;
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
        res.data = String::from(String::from_utf8_lossy(&data[ind..next_ind]).trim());
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
        dis = trimmed_data.parse::<f64>().ok();
        ind = next_ind;
    }
    // Checking for termination.
    if *data.get(ind)? == b',' {
        ind += 1;
    }

    return Some((res, dis, ind));
}

fn node_to_newick<NodeData: Display, EdgeData: Display>(node: &Node<NodeData, EdgeData>, buf: &mut String) {
    if !node.child_nodes.is_empty() {
        buf.write_str("(");
        for (child_node, edge_data) in node.child_nodes.iter() {
            let edge_repr = edge_data.to_string();
            node_to_newick(child_node, buf);
            if !edge_repr.is_empty() {
                buf.write_str(":");
                buf.write_str(&edge_data.to_string());
            }
            buf.write_str(", ");
        }
        buf.pop();
        buf.pop();
        buf.push(')');
    }

    buf.write_str(&node.data.to_string());
}

pub fn to_newick<NodeData: Display, EdgeData: Display>(tree: &Tree<NodeData, EdgeData>) -> String {
    let mut res = String::new();
    node_to_newick(&tree.root, &mut res);
    res.push(';');
    res
}

pub fn to_newick_path<NodeData: Display, EdgeData: Display>(path: &Path, tree: &Tree<NodeData, EdgeData>) {
    write_path(path, &to_newick(tree));
}

/// Parses data in newick format to a phylogenetic tree.
pub fn from_newick(data: &str) -> Tree<String, Option<f64>> {
    let (node, dis, _) = parse_newick(data.trim().as_bytes()).unwrap();
    assert!(dis.is_none(), "Distance set on root node!");
    return Tree { root: node };
}

pub fn from_newick_path(path: &Path) -> Tree<String, Option<f64>> {
    assert_eq!(path.extension().unwrap().to_str().unwrap(), "newick");
    from_newick(&read_path(path))
}

#[cfg(test)]
mod tests {
    use crate::phylogeny::tree::{Node, Tree};
    use crate::phylogeny::tree_parsers::{from_newick, to_newick};

    /// Testing that the tree initialization from newick strings works as expected, and that
    #[test]
    fn test_newick_initialization() {
        let mut cons_tree = Tree::from_root(Node::from_data("root"));
        cons_tree.root.child_nodes.push((Node::from_data("A"), Some(1.)));
        cons_tree.root.child_nodes.push((Node::from_data("B"), Some(0.1)));
        cons_tree.root.child_nodes.push((Node::from_data("C"), None));

        let data = String::from("  (A : 1 , B : .1 ,  C ) root ;  \n");

        let parsed_tree = from_newick(&data);

        assert_eq!(&cons_tree.root.data, &parsed_tree.root.data);
        assert_eq!(parsed_tree.root.child_nodes.len(), 3);

        for i in 0..3 {
            assert_eq!(
                cons_tree.root.child_nodes[i].0.data,
                parsed_tree.root.child_nodes[i].0.data
            );
            assert_eq!(
                cons_tree.root.child_nodes[i].1,
                parsed_tree.root.child_nodes[i].1
            );
        }
    }

    /// Testing postorder iteration.
    #[test]
    fn test_newick_reparse() {
        let data = String::from("(A, (B, C)D, E, F)G;");

        let parsed_tree = from_newick(&data);
        let newick = to_newick(&parsed_tree.map_edge_data(|_|""));

        assert_eq!(data, newick);
    }
}