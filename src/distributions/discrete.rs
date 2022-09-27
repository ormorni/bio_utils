use std::fmt::Debug;
use std::hash::Hash;
use fxhash::FxHashMap;
use rand::Rng;

/// A distribution over a discrete set of elements.
#[derive(Clone, PartialEq)]
pub struct DiscreteDist<T> where T: Hash + Eq + Clone {
    data: FxHashMap<T, f64>,
}

impl <T> DiscreteDist<T> where T: Hash + Eq + Debug + Clone {
    /// Initializes a discrete distribution from an iterator over the keys and a
    /// function mapping each key to its probability density.
    pub fn from_fn(iter: impl Iterator<Item=T>, mapper: impl Fn(&T) -> f64) -> DiscreteDist<T> {
        let data: FxHashMap<T, f64> = iter.map(|item|{
            let prob = mapper(&item);
            assert!(prob >= 0., "DiscreteDist initialized with negative probability: {:?}: {}", item, prob);
            (item, prob)
        }).collect();

        DiscreteDist::from_map(data)
    }

    /// Initializes a discrete distribution from a map, with the probability space being the keys
    /// and the probability density the values.
    pub fn from_map(mut data: FxHashMap<T, f64>) -> DiscreteDist<T> {
        let total_weight = data.values().sum::<f64>();

        for i in data.values_mut() {
            *i /= total_weight;
        }

        DiscreteDist {
            data
        }
    }

    /// Initializes a discrete distribution from a slice of samples.
    pub fn from_sample(sample: &[T]) -> DiscreteDist<T> {
        let mut data = FxHashMap::default();
        for t in sample {
            *data.entry(t.clone()).or_default() += 1./sample.len() as f64;
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

    /// Returns an integral over the function using the distribution as a measure.
    pub fn integrate(&self, function: impl Fn(&T) -> f64) -> f64 {
        self.data.iter().map(|(key, val)|*val * function(key)).sum::<f64>()
    }
}