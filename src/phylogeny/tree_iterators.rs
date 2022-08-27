use super::tree::*;

/// An iterator over the nodes of a phylogenetic tree.
/// Yields first the current node, then its child nodes.
#[derive(Clone)]
pub struct PreorderNodeIter<'t, NodeData, EdgeData> {
    /// The stack of nodes to process.
    node_stack: Vec<&'t Node<NodeData, EdgeData>>,
}

impl<'t, NodeData, EdgeData> PreorderNodeIter<'t, NodeData, EdgeData> {
    fn new(node: &'t Node<NodeData, EdgeData>) -> PreorderNodeIter<'t, NodeData, EdgeData> {
        PreorderNodeIter {
            node_stack: vec![node],
        }
    }
}

impl<'t, NodeData, EdgeData> Iterator for PreorderNodeIter<'t, NodeData, EdgeData> {
    type Item = &'t Node<NodeData, EdgeData>;

    fn next(&mut self) -> Option<Self::Item> {
        let node = self.node_stack.pop()?;
        for (child_node, _) in node.child_nodes.iter().rev() {
            self.node_stack.push(child_node);
        }
        return Some(node);
    }
}

/// An iterator over the nodes of a phylogenetic tree.
/// First iterates over the child nodes of the current node, then returns the node.
#[derive(Clone)]
pub struct PostorderNodeIter<'t, NodeData, EdgeData> {
    /// The stack of nodes to process. Since nodes are processed only after the child nodes are
    /// processed, keeps also a list of booleans specifying whether the node should be processed
    /// or expanded.
    node_stack: Vec<(&'t Node<NodeData, EdgeData>, bool)>,
}

impl<'t, NodeData, EdgeData> PostorderNodeIter<'t, NodeData, EdgeData> {
    fn new(node: &'t Node<NodeData, EdgeData>) -> PostorderNodeIter<'t, NodeData, EdgeData> {
        PostorderNodeIter {
            node_stack: vec![(node, false)],
        }
    }
}

impl<'t, NodeData, EdgeData> Iterator for PostorderNodeIter<'t, NodeData, EdgeData> {
    type Item = &'t Node<NodeData, EdgeData>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let (node, should_return) = self.node_stack.pop()?;
            if should_return {
                return Some(node);
            } else {
                self.node_stack.push((node, true));
                for (child_node, _) in node.child_nodes.iter().rev() {
                    self.node_stack.push((child_node, false));
                }
            }
        }
    }
}
#[derive(Clone)]
pub struct PreorderEdgeIterator<'t, NodeData, EdgeData> {
    node_stack: Vec<Edge<'t, NodeData, EdgeData>>,
}

/// A preorder iterator over the edges of the tree.
impl<'t, NodeData, EdgeData> Iterator for PreorderEdgeIterator<'t, NodeData, EdgeData> {
    type Item = Edge<'t, NodeData, EdgeData>;

    fn next(&mut self) -> Option<Self::Item> {
        let curr_edge = self.node_stack.pop()?;
        let child_node = curr_edge.child;
        for (grandchild_node, data) in child_node.child_nodes.iter().rev() {
            self.node_stack.push(Edge {
                parent: child_node,
                child: grandchild_node,
                data,
            });
        }

        Some(curr_edge)
    }
}

#[derive(Clone)]
pub struct PostorderEdgeIterator<'t, NodeData, EdgeData> {
    node_iter: PostorderNodeIter<'t, NodeData, EdgeData>,
    curr_node: &'t Node<NodeData, EdgeData>,
    curr_index: usize,
}

/// A Postorder iterator over the edges of the tree.
impl<'t, NodeData, EdgeData> Iterator for PostorderEdgeIterator<'t, NodeData, EdgeData> {
    type Item = Edge<'t, NodeData, EdgeData>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // We are iterating over
            while self.curr_index >= self.curr_node.child_nodes.len() {
                self.curr_node = self.node_iter.next()?;
                self.curr_index = 0;
            }
            let (child_node, data) = &self.curr_node.child_nodes[self.curr_index];
            self.curr_index += 1;
            return Some(Edge {
                parent: self.curr_node,
                child: child_node,
                data,
            });
        }
    }
}

impl<'t, NodeData, EdgeData> Tree<NodeData, EdgeData> {
    /// Returns an iterator over the nodes of the tree, in preorder iteration order.
    pub fn preorder_iter(&'t self) -> PreorderNodeIter<'t, NodeData, EdgeData> {
        self.root.preorder_iter()
    }
    /// Returns an iterator over the nodes of the tree, in postorder iteration order.
    pub fn postorder_iter(&'t self) -> PostorderNodeIter<'t, NodeData, EdgeData> {
        self.root.postorder_iter()
    }
    /// A generic iterator over the nodes, with no guaranteed order.
    pub fn nodes(&'t self) -> PreorderNodeIter<'t, NodeData, EdgeData> {
        self.root.nodes()
    }
    /// Returns an iterator over the edges of the tree, in preorder iteration order.
    pub fn preorder_edge_iter(&'t self) -> PreorderEdgeIterator<'t, NodeData, EdgeData> {
        self.root.preorder_edge_iter()
    }
    /// Returns an iterator over the edges of the tree, in posrtorder iteration order.
    pub fn postorder_edge_iter(&'t self) -> PostorderEdgeIterator<'t, NodeData, EdgeData> {
        self.root.postorder_edge_iter()
    }
    /// A generic iterator over the edges, with no guaranteed iteration order.
    pub fn edges(&'t self) -> PreorderEdgeIterator<'t, NodeData, EdgeData> {
        self.root.edges()
    }
    /// An iterator over the leaf nodes of the tree.
    pub fn leaf_nodes(&'t self) -> impl Iterator<Item = &'t Node<NodeData, EdgeData>> {
        self.root.leaf_nodes()
    }
}

impl<'t, NodeData, EdgeData> Node<NodeData, EdgeData> {
    /// Returns an iterator over the nodes of the tree, in preorder iteration order.
    pub fn preorder_iter(&'t self) -> PreorderNodeIter<'t, NodeData, EdgeData> {
        PreorderNodeIter::new(self)
    }

    /// Returns an iterator over the nodes of the tree, in postorder iteration order.
    pub fn postorder_iter(&'t self) -> PostorderNodeIter<'t, NodeData, EdgeData> {
        PostorderNodeIter::new(self)
    }

    /// A generic iterator over the nodes, with no guaranteed order.
    pub fn nodes(&'t self) -> PreorderNodeIter<'t, NodeData, EdgeData> {
        self.preorder_iter()
    }

    /// Returns an iterator over the edges of the tree, in preorder iteration order.
    pub fn preorder_edge_iter(&'t self) -> PreorderEdgeIterator<'t, NodeData, EdgeData> {
        PreorderEdgeIterator {
            node_stack: self
                .child_nodes
                .iter()
                .map(|(child, data)| Edge {
                    parent: self,
                    child,
                    data,
                })
                .rev()
                .collect(),
        }
    }

    /// Returns an iterator over the edges of the tree, in posrtorder iteration order.
    pub fn postorder_edge_iter(&'t self) -> PostorderEdgeIterator<'t, NodeData, EdgeData> {
        let mut node_iter = self.postorder_iter();
        let node = node_iter.next();
        PostorderEdgeIterator {
            node_iter,
            curr_node: node.unwrap_or(&self),
            curr_index: 0,
        }
    }

    /// A generic iterator over the edges, with no guaranteed iteration order.
    pub fn edges(&'t self) -> PreorderEdgeIterator<'t, NodeData, EdgeData> {
        self.preorder_edge_iter()
    }

    /// An iterator over the leaf nodes of the tree.
    pub fn leaf_nodes(&'t self) -> impl Iterator<Item = &'t Node<NodeData, EdgeData>> {
        self.nodes().filter(|node| node.is_leaf())
    }
}

#[cfg(test)]
mod tests {
    use crate::phylogeny::tree_parsers::{from_newick};

    /// Testing postorder iteration.
    #[test]
    fn test_postorder() {
        let data = String::from("(A, (B, C)D, E, F)G;");

        let parsed_tree = from_newick(&data);

        let expected_names = vec!["A", "B", "C", "D", "E", "F", "G"];
        let mut iter_names = Vec::new();
        for node in parsed_tree.postorder_iter() {
            iter_names.push(&node.data);
        }

        assert_eq!(expected_names, iter_names);
    }

    /// Testing preorder iteration.
    #[test]
    fn test_preorder() {
        let data = String::from("(A, (B, C)D, E, F)G;");

        let parsed_tree = from_newick(&data);

        let expected_names = vec!["G", "A", "D", "B", "C", "E", "F"];
        let mut iter_names = Vec::new();
        for node in parsed_tree.preorder_iter() {
            iter_names.push(&node.data);
        }

        assert_eq!(expected_names, iter_names);
    }
}
