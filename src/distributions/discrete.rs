use std::collections::HashMap;
use std::fmt::Debug;
use std::hash::Hash;
use rand::Rng;

/// A distribution over a discrete set of elements.
pub struct DiscreteDist<T> where T: Hash + Eq {
    data: HashMap<T, f64>,
}

impl <T> DiscreteDist<T> where T: Hash + Eq + Debug {
    /// Initializes a discrete distribution from an iterator over the keys and a
    /// function mapping each key to its probability density.
    pub fn from_fn(iter: impl Iterator<Item=T>, mapper: impl Fn(&T) -> f64) -> DiscreteDist<T> {
        let data: HashMap<T, f64> = iter.map(|item|{
            let prob = mapper(&item);
            assert!(prob >= 0., "DiscreteDist initialized with negative probability: {:?}: {}", item, prob);
            (item, prob)
        }).collect();

        DiscreteDist::from_map(data)
    }

    /// Initializes a discrete distribution from a map, with the probability space being the keys
    /// and the probability density the values.
    pub fn from_map(mut data: HashMap<T, f64>) -> DiscreteDist<T> {
        let total_weight = data.values().sum::<f64>();

        for i in data.values_mut() {
            *i /= total_weight;
        }

        DiscreteDist {
            data
        }
    }

    /// Samples an element according to the distribution.
    pub fn sample(&self, rng: &mut impl Rng) -> &T {
        let mut chooser: f64 = rng.gen();

        for (item, prob) in self.data.iter() {
            chooser -= prob;
            if chooser <= 0. {
                return item;
            }
        }

        // In case of numerical errors, return the first element.
        self.data.iter().next().unwrap().0
    }
}