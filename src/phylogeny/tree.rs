use std::fmt::Debug;
use std::hash::{Hash, Hasher};
use std::sync::atomic::{AtomicUsize, Ordering};

/// An implementation of a phylogenetic tree.
#[derive(Debug, Eq, PartialEq, Default)]
pub struct Tree<NodeData, EdgeData> {
    pub root: Node<NodeData, EdgeData>,
}

static NODE_COUNTER: AtomicUsize = AtomicUsize::new(0);

/// A node in the phylogenetic tree.
#[derive(Debug, Clone, Default)]
pub struct Node<NodeData, EdgeData> {
    pub data: NodeData,
    pub child_nodes: Vec<(Node<NodeData, EdgeData>, EdgeData)>,
    pub id: usize,
}

/// An edge in the phylogenetic tree.
/// Since the nodes are heavy, holds only a reference to the nodes.
#[derive(Debug, Copy, Clone)]
pub struct Edge<'t, NodeData, EdgeData> {
    pub parent: &'t Node<NodeData, EdgeData>,
    pub child: &'t Node<NodeData, EdgeData>,
    pub data: &'t EdgeData,
}

impl <'t, NodeData, EdgeData> Edge<'t, NodeData, EdgeData> {
    pub fn new(parent: &'t Node<NodeData, EdgeData>, child: &'t Node<NodeData, EdgeData>, data: &'t EdgeData) -> Edge<'t, NodeData, EdgeData> {
        Edge {parent, child, data }
    }
}

impl <NodeData: Default, EdgeData> Tree<NodeData, EdgeData> {
    /// Initializing a new empty tree.
    pub fn new() -> Tree<NodeData, EdgeData> {
        Tree { root: Node::new() }
    }
}

impl <NodeData, EdgeData> Tree<NodeData, EdgeData> {
    /// Initializing a new empty tree.
    pub fn from_root(root: Node<NodeData, EdgeData>) -> Tree<NodeData, EdgeData> {
        Tree { root }
    }
}

impl <NodeData, EdgeData: Clone> Tree<NodeData, EdgeData> {
    pub fn map_node_data<NewData, F : Clone + Fn(&NodeData) -> NewData>(&self, func: F) -> Tree<NewData, EdgeData> {
        Tree::from_root(self.root.map_node_data(func))
    }
}

impl <NodeData: Clone, EdgeData> Tree<NodeData, EdgeData> {
    pub fn map_edge_data<NewData, F : Clone + Fn(&EdgeData) -> NewData>(&self, func: F) -> Tree<NodeData, NewData> {
        Tree::from_root(self.root.map_edge_data(func))
    }
}

impl <NodeData: Default, EdgeData> Node<NodeData, EdgeData> {
    /// Initializes an empty node.
    pub fn new() -> Node<NodeData, EdgeData> {
        Node {
            data: NodeData::default(),
            child_nodes: vec![],
            id: NODE_COUNTER.fetch_add(1, Ordering::SeqCst),
        }
    }
}

impl <NodeData, EdgeData> Node<NodeData, EdgeData> {

    /// Initializes an empty node with a given name.
    pub fn from_data(data: NodeData) -> Node<NodeData, EdgeData> {
        Node {
            data,
            child_nodes: vec![],
            id: NODE_COUNTER.fetch_add(1, Ordering::SeqCst),
        }
    }
    /// Returns if the node has no child nodes.
    pub fn is_leaf(&self) -> bool {
        self.child_nodes.is_empty()
    }

    pub fn child_edges(&self) -> impl Iterator<Item=Edge<NodeData, EdgeData>> + Clone {
        self.child_nodes.iter()
            .map(|(node, data)|Edge::new(self, node, data))
    }
}

impl <NodeData, EdgeData: Clone> Node<NodeData, EdgeData> {
    fn map_node_data<NewData, F : Clone + Fn(&NodeData) -> NewData>(&self, func: F) -> Node<NewData, EdgeData> {
        let mut res = Node::from_data(func(&self.data));
        for (child_node, child_edge) in self.child_nodes.iter() {
            res.child_nodes.push((child_node.map_node_data(func.clone()), child_edge.clone()));
        }
        res
    }
}

impl <NodeData: Clone, EdgeData> Node<NodeData, EdgeData> {
    fn map_edge_data<NewData, F : Clone + Fn(&EdgeData) -> NewData>(&self, func: F) -> Node<NodeData, NewData> {
        let mut res = Node::from_data(self.data.clone());
        for (child_node, child_edge) in self.child_nodes.iter() {
            res.child_nodes.push((child_node.map_edge_data(func.clone()), func(child_edge)));
        }
        res
    }
}


impl <NodeData, EdgeData> Hash for Node<NodeData, EdgeData> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.id.hash(state);
    }
}

impl<NodeData, EdgeData> PartialEq<Self> for Node<NodeData, EdgeData> {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl <NodeData, EdgeData> std::cmp::Eq for Node<NodeData, EdgeData> {}


impl<'t, NodeData, EdgeData> PartialEq<Self> for Edge<'t, NodeData, EdgeData> {
    fn eq(&self, other: &Self) -> bool {
        self.parent == other.parent && self.child == other.child
    }
}

impl<'t, NodeData, EdgeData> Eq for Edge<'t, NodeData, EdgeData> {}

impl<'t, NodeData, EdgeData> Hash for Edge<'t, NodeData, EdgeData> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.parent.hash(state);
        self.child.hash(state);
    }
}
