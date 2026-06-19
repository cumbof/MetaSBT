//! DeltaTree: High-performance Rust backend for MetaSBT v2.0
//! 
//! This module implements the core algorithms for the Delta-SBT architecture,
//! which drastically reduces the disk footprint of Sequence Bloom Trees by 
//! vertically stripping redundant "Core" k-mers from child nodes, storing
//! only strictly disjoint "Delta" (accessory) k-mers at the leaves.
//! 
//! It leverages:
//! - `needletail` for blazing fast FASTA parsing.
//! - `twox-hash` for rapid k-mer hashing.
//! - `roaring` (Roaring Bitmaps) for highly compressed, bitwise-operable sets.
//! - `pyo3` to expose these functions as a native Python extension.

use pyo3::prelude::*;
use pyo3::exceptions::{PyIOError, PyValueError};
use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use roaring::RoaringBitmap;
use needletail::parse_fastx_file;
use twox_hash::XxHash64;
use std::hash::Hasher;

/// A simple helper function to compute the reverse complement of a DNA sequence.
/// It maps A <-> T and C <-> G, defaulting unknown characters to N.
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev() // Reverse the sequence
        .map(|&c| match c {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        })
        .collect()
}

/// Generate a FracMinHash sketch and save it as a compressed Roaring Bitmap.
///
/// FracMinHash sub-samples the k-mer space deterministically. Instead of storing 
/// every k-mer in a genome, it only stores hashes that fall below a certain threshold.
/// This reduces the sketch size by a factor of `scaled` while preserving distance metrics.
///
/// # Arguments
/// * `filepath` - Path to the input FASTA file.
/// * `out_filepath` - Path to save the serialized RoaringBitmap.
/// * `kmer_size` - The length of k-mers to extract (e.g., 21, 31).
/// * `scaled` - The scale factor (e.g., 1000 means keep 1 in 1000 k-mers).
/// * `_threads` - Number of threads (currently unused, reserved for future Rayon integration).
#[pyfunction]
fn sketch(filepath: &str, out_filepath: &str, kmer_size: usize, scaled: u32, _threads: usize) -> PyResult<()> {
    let mut bitmap = RoaringBitmap::new();
    
    // Calculate the maximum hash value allowed for the given scaled fraction.
    // We cast u32::MAX to f64 to avoid overflow during division.
    // RoaringBitmaps store 32-bit integers, making u32 the perfect target type.
    let max_hash = (u32::MAX as f64 / scaled as f64) as u32;

    // Open the FASTA file using needletail for fast parsing
    let mut reader = parse_fastx_file(filepath)
        .map_err(|e| PyIOError::new_err(format!("Failed to open {}: {}", filepath, e)))?;

    while let Some(record) = reader.next() {
        let seqrec = record.map_err(|e| PyValueError::new_err(e.to_string()))?;
        let seq = seqrec.seq();
        
        // Skip sequences shorter than the k-mer size
        if seq.len() < kmer_size {
            continue;
        }

        // Iterate over sliding windows of size k
        for kmer in seq.windows(kmer_size) {
            // Ignore k-mers containing ambiguous bases (N)
            if kmer.iter().any(|&b| b == b'N' || b == b'n') {
                continue;
            }

            // Determine the canonical k-mer (the lexicographically smaller of the forward
            // sequence and its reverse complement). This ensures strand independence.
            let revcomp = reverse_complement(kmer);
            let canonical = if kmer < revcomp.as_slice() { kmer } else { revcomp.as_slice() };

            // Hash the canonical k-mer using XxHash64 (extremely fast, non-cryptographic)
            let mut hasher = XxHash64::with_seed(0);
            hasher.write(canonical);
            let h = hasher.finish() as u32; // Truncate down to u32 for RoaringBitmap compatibility

            // Apply the FracMinHash condition: only keep the hash if it is <= max_hash
            if h <= max_hash {
                bitmap.insert(h);
            }
        }
    }

    // Serialize the highly compressed RoaringBitmap to disk
    let mut out_file = File::create(out_filepath)
        .map_err(|e| PyIOError::new_err(format!("Failed to create {}: {}", out_filepath, e)))?;
    
    bitmap.serialize_into(&mut out_file)
        .map_err(|e| PyIOError::new_err(format!("Failed to serialize bitmap to {}: {}", out_filepath, e)))?;

    Ok(())
}

/// Compute Containment-based Average Nucleotide Identity (ANI) 
/// between a focus sketch and a list of target sketches.
///
/// Because we use FracMinHash, we can estimate ANI purely mathematically 
/// without needing full alignments.
#[pyfunction]
fn containment_ani(focus: &str, targets: Vec<&str>, kmer_size: usize) -> PyResult<HashMap<String, f64>> {
    // 1. Load the focus query sketch
    let mut focus_file = File::open(focus)
        .map_err(|e| PyIOError::new_err(format!("Failed to open focus sketch {}: {}", focus, e)))?;
    
    let focus_bm = RoaringBitmap::deserialize_from(&mut focus_file)
        .map_err(|e| PyIOError::new_err(format!("Failed to deserialize focus sketch {}: {}", focus, e)))?;
    
    let focus_len = focus_bm.len() as f64;
    let mut results = HashMap::new();

    // Edge case: If the focus sketch is completely empty, distance is maximum (1.0)
    if focus_len == 0.0 {
        for t in targets { results.insert(t.to_string(), 1.0); }
        return Ok(results);
    }

    // 2. Compare against every target sketch
    for target in targets {
        let mut t_file = File::open(target)
            .map_err(|e| PyIOError::new_err(format!("Failed to open target sketch {}: {}", target, e)))?;
        
        let target_bm = RoaringBitmap::deserialize_from(&mut t_file)
            .map_err(|e| PyIOError::new_err(format!("Failed to deserialize target sketch {}: {}", target, e)))?;
        
        // Fast bitwise intersection length computed directly on the compressed RoaringBitmap!
        let intersection = focus_bm.intersection_len(&target_bm) as f64;
        
        // Containment Index (C) = |Focus ∩ Target| / |Focus|
        let containment = intersection / focus_len;
        
        // Estimate ANI using the Mash distance formula: ANI ≈ 1 + (1/k) * ln(C)
        let mut ani = 1.0;
        if containment > 0.0 {
            ani = 1.0 + (1.0 / kmer_size as f64) * containment.ln();
        } else {
            ani = 0.0;
        }

        // MetaSBT expects a distance metric where 0.0 is identical and 1.0 is completely different
        let distance = if ani <= 0.0 { 
            1.0 
        } else if ani >= 1.0 { 
            0.0 
        } else { 
            1.0 - ani 
        };

        results.insert(target.to_string(), distance);
    }

    Ok(results)
}

/// Build the Delta-SBT tree structure from a list of children sketches.
///
/// In a Delta-SBT, an internal node (Core) is the mathematical intersection of its children.
/// The children are then modified to only contain their unique "Delta" (Accessory) genes.
#[pyfunction]
fn build_delta_tree(sketches_list: &str, out_tree: &str, is_flat: bool) -> PyResult<()> {
    let list_file = File::open(sketches_list)
        .map_err(|e| PyIOError::new_err(format!("Cannot open sketches list: {}", e)))?;
    let reader = BufReader::new(list_file);

    let mut child_paths = Vec::new();
    let mut child_bms = Vec::new();
    let mut core_bitmap = RoaringBitmap::new();
    let mut first = true;

    // PASS 1 (Bottom-Up Step): Load all children and compute their strictly shared Core (Intersection)
    for line in reader.lines() {
        let filepath = line.map_err(|e| PyIOError::new_err(e.to_string()))?;
        let filepath = filepath.trim().to_string();
        if filepath.is_empty() { continue; }

        let mut t_file = File::open(&filepath)
            .map_err(|e| PyIOError::new_err(format!("Failed to open sketch {}: {}", filepath, e)))?;
        
        let t_bm = RoaringBitmap::deserialize_from(&mut t_file)
            .map_err(|e| PyIOError::new_err(format!("Failed to deserialize {}: {}", filepath, e)))?;
        
        if first {
            core_bitmap = t_bm.clone();
            first = false;
        } else {
            // Core is the strict mathematical intersection of all children: C = C ∩ Child_i
            core_bitmap &= &t_bm;
        }

        child_paths.push(filepath);
        child_bms.push(t_bm);
    }

    // PASS 2 (Top-Down Step): Strip the Core from the children, leaving strictly disjoint Deltas
    // We skip this if the database is configured as 'flat' (no hierarchical stripping)
    if !is_flat {
        for (path, mut bm) in child_paths.into_iter().zip(child_bms.into_iter()) {
            // DELTA = CHILD \ CORE
            // This removes all shared sequences, shrinking the child sketch by ~95%
            bm -= &core_bitmap;

            // Overwrite the child sketch on disk with its new, stripped Delta version
            let mut out_child = File::create(&path)
                .map_err(|e| PyIOError::new_err(format!("Failed to overwrite delta child {}: {}", path, e)))?;
            bm.serialize_into(&mut out_child)
                .map_err(|e| PyIOError::new_err(format!("Failed to serialize delta child: {}", e)))?;
        }
    }

    // PASS 3: Save the Core intersection as the signpost filter for this parent node
    let mut out_file = File::create(out_tree)
        .map_err(|e| PyIOError::new_err(format!("Cannot create tree out file: {}", e)))?;
    
    core_bitmap.serialize_into(&mut out_file)
        .map_err(|e| PyIOError::new_err(format!("Cannot serialize tree: {}", e)))?;

    Ok(())
}

/// The Accumulator Search Algorithm for Delta-SBT traversal.
///
/// Standard SBTs search full filters. Delta-SBTs must logically piece filters back together 
/// as they traverse down. Because the sets are strictly disjoint (thanks to `build_delta_tree`),
/// we can use the distributive property: |Q ∩ (Core ∪ Delta)| = |Q ∩ Core| + |Q ∩ Delta|.
/// 
/// We keep a running tally (`accumulated_score`) and never have to decompress or rebuild
/// the full Bloom Filters in memory!
///
/// # Arguments
/// * `query_sketch` - Path to the query genome's sketch.
/// * `tree_root` - Path to the database root Core filter.
/// * `tree_topology` - A dictionary mapping a parent node path to a list of its children `(path, level)`.
/// * `kmer_size` - Used for ANI estimation.
/// * `theta` - The minimum containment threshold to continue searching down a branch.
/// * `_uncertainty` - Currently unused in core rust loop, reserved for advanced pruning.
#[pyfunction]
fn accumulator_search(
    query_sketch: &str,
    tree_root: &str,
    tree_topology: HashMap<String, Vec<(String, String)>>,
    kmer_size: usize,
    theta: f64,
    _uncertainty: f64
) -> PyResult<HashMap<String, HashMap<String, f64>>> {
    
    // Load the query sketch into memory once
    let mut q_file = File::open(query_sketch).map_err(|e| PyIOError::new_err(e.to_string()))?;
    let query_bm = RoaringBitmap::deserialize_from(&mut q_file).map_err(|e| PyIOError::new_err(e.to_string()))?;
    let query_len = query_bm.len() as f64;

    // This will hold the final profile mapping: { level_name: { node_path: distance } }
    let mut profiles: HashMap<String, HashMap<String, f64>> = HashMap::new();

    // Initialize the Breadth-First Search (BFS) queue.
    // Tuple holds: (Current Node File Path, Running Tally of Matches, Taxonomic Level)
    let mut queue = VecDeque::new();
    queue.push_back((tree_root.to_string(), 0.0, "db".to_string()));

    while let Some((node_path, mut accumulated_score, level_name)) = queue.pop_front() {
        if !Path::new(&node_path).exists() { continue; }

        let mut node_file = File::open(&node_path).map_err(|e| PyIOError::new_err(e.to_string()))?;
        let node_bm = RoaringBitmap::deserialize_from(&mut node_file).map_err(|e| PyIOError::new_err(e.to_string()))?;

        // 1. Accumulate Score: 
        // Compute how many query k-mers hit this specific Delta, and add it to the tally
        // inherited from the parent nodes.
        let node_intersection = query_bm.intersection_len(&node_bm) as f64;
        accumulated_score += node_intersection;

        // 2. Evaluate Threshold (Pruning Phase)
        let containment = accumulated_score / query_len;
        if containment < theta {
            // The query does not have enough k-mers in this branch. 
            // We stop pushing its children to the queue (Pruning).
            continue; 
        }

        // 3. Compute Distance for Profiling
        // Estimate ANI dynamically based on the accumulated score
        let mut ani = 1.0;
        if containment > 0.0 {
            ani = 1.0 + (1.0 / kmer_size as f64) * containment.ln();
        } else { 
            ani = 0.0; 
        }
        
        let distance = if ani <= 0.0 { 1.0 } else if ani >= 1.0 { 0.0 } else { 1.0 - ani };

        // Record the valid match
        profiles.entry(level_name.clone())
            .or_insert_with(HashMap::new)
            .insert(node_path.clone(), distance);

        // 4. Graph Traversal: Look up this node's children in the topology dict
        // and enqueue them for processing, passing the current accumulated score down to them.
        if let Some(children) = tree_topology.get(&node_path) {
            for (child_path, next_level) in children {
                queue.push_back((child_path.clone(), accumulated_score, next_level.clone()));
            }
        }
    }

    Ok(profiles)
}

/// Dynamic Rebalancing / Delta update for a new genome insertion.
///
/// Inserting a new genome into a disjoint Delta tree requires shifting the consensus.
/// If the new genome is missing a "Core" gene, that gene is no longer core! It must be
/// evicted from the parent node and pushed down into the Deltas of the existing siblings.
///
/// # Arguments
/// * `new_sketch` - Path to the newly sketched genome (contains ALL its k-mers).
/// * `species_node_path` - Path to the species Core filter it is being added to.
/// * `sibling_sketches` - Paths to the Delta filters of all existing strains in this species.
#[pyfunction]
fn update_delta_tree(new_sketch: &str, species_node_path: &str, sibling_sketches: Vec<&str>) -> PyResult<()> {
    // 1. Load the new genome sketch (Full complement of k-mers, N)
    let mut n_file = File::open(new_sketch).map_err(|e| PyIOError::new_err(e.to_string()))?;
    let mut new_bm = RoaringBitmap::deserialize_from(&mut n_file).map_err(|e| PyIOError::new_err(e.to_string()))?;

    // 2. Load the current species Core (Intersection of all existing siblings, C)
    let mut s_file = File::open(species_node_path).map_err(|e| PyIOError::new_err(e.to_string()))?;
    let species_bm = RoaringBitmap::deserialize_from(&mut s_file).map_err(|e| PyIOError::new_err(e.to_string()))?;

    // 3. Rebalancing Math:
    // The strictly new core (C') MUST be the intersection of the Old Core and the New Genome
    // C' = C ∩ N
    let mut new_core = species_bm.clone();
    new_core &= &new_bm;

    // The Lost Core (L) = C \ C'
    // These are k-mers that were previously shared by all strains, but the new strain lacks them!
    let mut lost_core = species_bm.clone();
    lost_core -= &new_core;

    // 4. Execute the Rebalance IF the core shifted (L is not empty)
    if !lost_core.is_empty() {
        // Push the Lost Core down into the existing Delta of every sibling
        // This ensures no genetic information is accidentally deleted from the tree
        for sibling_path in sibling_sketches {
            let mut sib_file = File::open(sibling_path).map_err(|e| PyIOError::new_err(e.to_string()))?;
            let mut sib_bm = RoaringBitmap::deserialize_from(&mut sib_file).map_err(|e| PyIOError::new_err(e.to_string()))?;
            
            sib_bm |= &lost_core; // Bitwise OR: Inject lost core into sibling delta
            
            // Overwrite sibling delta on disk
            let mut sib_out = File::create(sibling_path).map_err(|e| PyIOError::new_err(e.to_string()))?;
            sib_bm.serialize_into(&mut sib_out).map_err(|e| PyIOError::new_err(e.to_string()))?;
        }

        // Overwrite the species node with the new, mathematically strict (smaller) core
        let mut core_out = File::create(species_node_path).map_err(|e| PyIOError::new_err(e.to_string()))?;
        new_core.serialize_into(&mut core_out).map_err(|e| PyIOError::new_err(e.to_string()))?;
    }

    // 5. Finally, compute the Delta for the newly inserted genome
    // D_new = N \ C'
    new_bm -= &new_core;

    // Overwrite the original full-genome sketch with its highly compressed Leaf Delta
    let mut out_file = File::create(new_sketch).map_err(|e| PyIOError::new_err(e.to_string()))?;
    new_bm.serialize_into(&mut out_file).map_err(|e| PyIOError::new_err(e.to_string()))?;

    Ok(())
}

/// Read a serialized RoaringBitmap sketch file and return its cardinality.
///
/// Cardinality is the exact number of elements (subsampled k-mer hashes) stored
/// in the FracMinHash sketch. In the Delta-SBT architecture this replaces the
/// obsolete "density" metric (ratio of set bits to total bits) that was only
/// meaningful for fixed-size Bloom filters.
///
/// For a delta-encoded tree, comparing the cardinality of child nodes against
/// their parent's Core provides a direct measure of compression efficiency.
#[pyfunction]
fn sketch_cardinality(sketch_path: &str) -> PyResult<u64> {
    let mut file = File::open(sketch_path)
        .map_err(|e| PyIOError::new_err(format!("Failed to open sketch {}: {}", sketch_path, e)))?;
    let bitmap = RoaringBitmap::deserialize_from(&mut file)
        .map_err(|e| PyIOError::new_err(format!("Failed to deserialize sketch {}: {}", sketch_path, e)))?;
    Ok(bitmap.len())
}

/// The Python Module Definition.
/// 
/// This macro creates the entry points that `pyo3` and `maturin` will compile 
/// into the Python extension module. The name of the function `deltatree` MUST 
/// match the `lib.name` setting in the `Cargo.toml`.
#[pymodule]
fn deltatree(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sketch, m)?)?;
    m.add_function(wrap_pyfunction!(sketch_cardinality, m)?)?;
    m.add_function(wrap_pyfunction!(containment_ani, m)?)?;
    m.add_function(wrap_pyfunction!(accumulator_search, m)?)?;
    m.add_function(wrap_pyfunction!(build_delta_tree, m)?)?;
    m.add_function(wrap_pyfunction!(update_delta_tree, m)?)?;
    Ok(())
}
