use std::collections::HashSet;

use filecoin_hashers::poseidon::PoseidonHasher;
use storage_proofs_core::{crypto::feistel, drgraph::BASE_DEGREE, is_legacy_porep_id, PoRepID};

use storage_proofs_porep::stacked::{StackedBucketGraph, EXP_DEGREE};

#[test]
fn test_is_legacy() {
    fn p(v: u64) -> PoRepID {
        let mut res = [0u8; 32];
        res[..8].copy_from_slice(&v.to_le_bytes());
        res
    }

    assert!(is_legacy_porep_id(p(0)));
    assert!(is_legacy_porep_id(p(1)));
    assert!(is_legacy_porep_id(p(4)));
    assert!(!is_legacy_porep_id(p(5)));
}

// Test that 3 (or more) rounds of the Feistel cipher can be used
// as a pseudorandom permutation, that is, each input will be mapped
// to a unique output (and though not test here, since the cipher
// is symmetric, the decryption rounds also work as the inverse
// permutation), for more details see:
// https://en.wikipedia.org/wiki/Feistel_cipher#Theoretical_work.
#[test]
fn test_shuffle() {
    let n = 2_u64.pow(10);
    let d = EXP_DEGREE as u64;
    // Use a relatively small value of `n` as Feistel is expensive (but big
    // enough that `n >> d`).

    let mut shuffled: HashSet<u64> = HashSet::with_capacity((n * d) as usize);

    let feistel_keys = &[1, 2, 3, 4];
    let feistel_precomputed = feistel::precompute((n * d) as feistel::Index);

    for i in 0..n {
        for k in 0..d {
            let permuted = feistel::permute(n * d, i * d + k, feistel_keys, feistel_precomputed);

            // Since the permutation implies a one-to-one correspondence,
            // traversing the entire input space should generate the entire
            // output space (in `shuffled`) without repetitions (since a duplicate
            // output would imply there is another output that wasn't generated
            // and the permutation would be incomplete).
            assert!(shuffled.insert(permuted));
        }
    }

    // Actually implied by the previous `assert!` this is left in place as an
    // extra safety check that indeed the permutation preserved all the output
    // space (of `n * d` nodes) without repetitions (which the `HashSet` would
    // have skipped as duplicates).
    assert_eq!(shuffled.len(), (n * d) as usize);
}

#[test]
/// The initial implementation had a bug which prevented parents from ever falling in the later half of a sector.
/// In fact, it is even worse than that, in the case of 64GiB sectors.
/// This test demonstrates conclusively that non-legacy graphs do not suffer from this pathology.
/// It also suggests, inconclusively, that legacy graphds do suffer from it (which we already know).
fn test_graph_distribution_pathology() {
    let sector32_nodes: u32 = 1 << 30;
    let sector64_nodes: u32 = 1 << 31;

    let porep_id = |id: u8| {
        let mut porep_id = [0u8; 32];
        porep_id[0] = id;

        porep_id
    };

    test_pathology_aux(porep_id(3), sector32_nodes);
    test_pathology_aux(porep_id(4), sector64_nodes);

    test_pathology_aux(porep_id(8), sector32_nodes);
    test_pathology_aux(porep_id(9), sector64_nodes);
}

fn test_pathology_aux(porep_id: PoRepID, nodes: u32) {
    // In point of fact, the concrete graphs expected to be non-pathological
    // appear to demonstrate this immediately (i.e. in the first node). We
    // test more than that just to make the tentative diagnosis of pathology
    // more convincing in the cases where we expect it. In the interest of
    // keeping the tests brief, we keep this fairly small, though, since we
    // already know the previous porep_ids exhibit the problem. The main
    // reason to test those cases at all is to convince ourselves the test
    // is sound.
    let test_n = 1_000;

    let expect_pathological = is_legacy_porep_id(porep_id);

    let graph = StackedBucketGraph::<PoseidonHasher>::new_stacked(
        nodes as usize,
        BASE_DEGREE,
        EXP_DEGREE,
        porep_id,
    )
    .unwrap();

    // If a parent index is not less than half the total node count, then
    // the parent falls in the second half of the previous layer. By the
    // definition of 'pathology' used here, that means the graph producing
    // this parent is not pathological.
    let demonstrably_large_enough = |p: &u32| *p >= (nodes / 2);

    dbg!(&porep_id, &nodes, &expect_pathological);
    for i in 0..test_n {
        let mut expanded_parents = [0u32; EXP_DEGREE];
        graph.expanded_parents(i, &mut expanded_parents).unwrap();

        if expect_pathological {
            // If we ever see a large-enough parent, then this graph is not
            // pathological, so the test fails.
            assert!(
                !expanded_parents.iter().any(demonstrably_large_enough),
                "Expected pathological graph but found large-enough parent."
            );
        } else {
            if expanded_parents.iter().any(demonstrably_large_enough) {
                // If we ever see a large-enough parent, then this graph is
                // not pathological, and the test succeeds. This is the only
                // way for a test expecting a non-pathological graph to
                // succeed, so there is no risk of false negatives (i.e.
                // failure to identify pathological graphs when unexpected).
                return;
            }
        }
    }

    // If we get here, we did not observe a parent large enough to conclude
    // that the graph is not pathological. In that case, the test fails if we
    // expected a non-pathological graph and succeeds otherwise. NOTE: this
    // could lead us to conclude that an actually non-pathological graph is
    // pathological, if `test_n` is set too low. Since the primary purpose
    // of this test is to assure us that newer graphs are not pathological,
    // it suffices to set `test_n` high enough to detect that.
    assert!(expect_pathological, "Did not expect pathological graph, but did not see large-enough parent to prove otherwise.");
}

// Tests that the set of expander edges has not been truncated.
#[test]
fn test_high_parent_bits() {
    // 64GiB sectors have 2^31 nodes.
    const N_NODES: usize = 1 << 31;

    // `u32` truncation would reduce the expander edge bit-length from 34 bits to 32 bits, thus
    // the first parent truncated would be the node at index `2^32 / EXP_DEGREE = 2^29`.
    const FIRST_TRUNCATED_PARENT: u32 = 1 << 29;

    // The number of child nodes to test before failing. This value was chosen arbitrarily and
    // can be changed.
    const N_CHILDREN_SAMPLED: usize = 3;

    // Non-legacy porep-id.
    let mut porep_id = [0u8; 32];
    porep_id[..8].copy_from_slice(&5u64.to_le_bytes());

    let graph = StackedBucketGraph::<PoseidonHasher>::new_stacked(
        N_NODES,
        BASE_DEGREE,
        EXP_DEGREE,
        porep_id,
    )
    .unwrap();

    let mut exp_parents = [0u32; EXP_DEGREE];
    for v in 0..N_CHILDREN_SAMPLED {
        graph.expanded_parents(v, &mut exp_parents[..]).unwrap();
        if exp_parents.iter().any(|u| *u >= FIRST_TRUNCATED_PARENT) {
            return;
        }
    }
    assert!(false);
}

// Checks that the distribution of parent node indexes within a sector is within a set bound.
#[test]
fn test_exp_parent_histogram() {
    // 64GiB sectors have 2^31 nodes.
    const N_NODES: usize = 1 << 31;

    // The number of children used to construct the histogram. This value is chosen
    // arbitrarily and can be changed.
    const N_CHILDREN_SAMPLED: usize = 10000;

    // The number of bins used to partition the set of sector nodes. This value was chosen
    // arbitrarily and can be changed to any integer that is a multiple of `EXP_DEGREE` and
    // evenly divides `N_NODES`.
    const N_BINS: usize = 32;
    const N_NODES_PER_BIN: u32 = (N_NODES / N_BINS) as u32;
    const PARENT_COUNT_PER_BIN_UNIFORM: usize = N_CHILDREN_SAMPLED * EXP_DEGREE / N_BINS;

    // This test will pass if every bin's parent count is within the bounds:
    // `(1 +/- FAILURE_THRESHOLD) * PARENT_COUNT_PER_BIN_UNIFORM`.
    const FAILURE_THRESHOLD: f32 = 0.4;
    const MAX_PARENT_COUNT_ALLOWED: usize =
        ((1.0 + FAILURE_THRESHOLD) * PARENT_COUNT_PER_BIN_UNIFORM as f32) as usize - 1;
    const MIN_PARENT_COUNT_ALLOWED: usize =
        ((1.0 - FAILURE_THRESHOLD) * PARENT_COUNT_PER_BIN_UNIFORM as f32) as usize + 1;

    // Non-legacy porep-id.
    let mut porep_id = [0u8; 32];
    porep_id[..8].copy_from_slice(&5u64.to_le_bytes());

    let graph = StackedBucketGraph::<PoseidonHasher>::new_stacked(
        N_NODES,
        BASE_DEGREE,
        EXP_DEGREE,
        porep_id,
    )
    .unwrap();

    // Count the number of parents in each bin.
    let mut hist = [0usize; N_BINS];
    let mut exp_parents = [0u32; EXP_DEGREE];
    for sample_index in 0..N_CHILDREN_SAMPLED {
        let v = sample_index * N_NODES / N_CHILDREN_SAMPLED;
        graph.expanded_parents(v, &mut exp_parents[..]).unwrap();
        for u in exp_parents.iter() {
            let bin_index = (u / N_NODES_PER_BIN) as usize;
            hist[bin_index] += 1;
        }
    }

    let success = hist.iter().all(|&n_parents| {
        n_parents >= MIN_PARENT_COUNT_ALLOWED && n_parents <= MAX_PARENT_COUNT_ALLOWED
    });

    assert!(success);
}
