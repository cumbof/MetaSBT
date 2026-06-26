//! DeltaTree: High-performance Rust backend for MetaSBT v2.0
//! 
//! This module implements the core algorithms for the Delta-SBT architecture,
//! which drastically reduces the disk footprint of Sequence Bloom Trees by 
//! vertically stripping redundant "Core" k-mers from child nodes, storing
//! only strictly disjoint "Delta" (accessory) k-mers at the leaves.
//! 
//! It leverages:
//! - `needletail` for blazing fast FASTA parsing.
//! - `nthash` for an O(1) rolling hash over DNA k-mers, plus a fused rolling
//!   polynomial hash over the 3-frame amino-acid translation.
//! - `roaring` (Roaring Bitmaps) for highly compressed, bitwise-operable sets.
//! - `pyo3` to expose these functions as a native Python extension.

use pyo3::prelude::*;
use pyo3::exceptions::{PyIOError, PyValueError};
use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader, Cursor, Read, Write};
use std::path::Path;
use std::sync::OnceLock;

use rayon::prelude::*;
use roaring::RoaringBitmap;
use needletail::parse_fastx_file;
use nthash::NtHashIterator;

/// Multiplicative base for the rolling polynomial hash over amino-acid k-mers
/// (the 64-bit FNV prime; any odd constant works as the polynomial base).
const AA_HASH_BASE: u64 = 0x0000_0100_0000_01b3;

/// SplitMix64 finalizer. A rolling polynomial hash is cheap but its low bits are
/// poorly distributed, which would bias FracMinHash sampling (we keep a hash when
/// its low 32 bits fall under `max_hash`). Running each rolling value through this
/// avalanche step restores a near-uniform distribution at O(1) cost, so the kept
/// fraction stays an unbiased ~1/scaled sample — exactly what ANI/AAI estimation
/// assumes.
#[inline]
fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

// One Rayon thread pool per process, sized on first sketch() call.
// Each mp.Pool worker process initialises its own independent copy.
static RAYON_POOL: OnceLock<rayon::ThreadPool> = OnceLock::new();

fn get_pool(nthreads: usize) -> &'static rayon::ThreadPool {
    RAYON_POOL.get_or_init(|| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(nthreads)
            .build()
            .expect("Failed to build Rayon thread pool")
    })
}

/// Translate a single DNA codon (3 bases) to an amino acid using the standard
/// genetic code. Returns b'X' for unknown codons.
fn dna_to_aa(codon: &[u8]) -> u8 {
    if codon.len() < 3 { return b'X'; }
    let a = codon[0].to_ascii_uppercase();
    let b = codon[1].to_ascii_uppercase();
    let c = codon[2].to_ascii_uppercase();
    match (a, b, c) {
        (b'T', b'T', b'T') | (b'T', b'T', b'C') => b'F',
        (b'T', b'T', b'A') | (b'T', b'T', b'G') => b'L',
        (b'T', b'C', b'T') | (b'T', b'C', b'C') | (b'T', b'C', b'A') | (b'T', b'C', b'G') => b'S',
        (b'T', b'A', b'T') | (b'T', b'A', b'C') => b'Y',
        (b'T', b'A', b'A') | (b'T', b'A', b'G') => b'*',
        (b'T', b'G', b'T') | (b'T', b'G', b'C') => b'C',
        (b'T', b'G', b'A') | (b'T', b'G', b'G') => b'W',
        (b'C', b'T', b'T') | (b'C', b'T', b'C') | (b'C', b'T', b'A') | (b'C', b'T', b'G') => b'L',
        (b'C', b'C', b'T') | (b'C', b'C', b'C') | (b'C', b'C', b'A') | (b'C', b'C', b'G') => b'P',
        (b'C', b'A', b'T') | (b'C', b'A', b'C') => b'H',
        (b'C', b'A', b'A') | (b'C', b'A', b'G') => b'Q',
        (b'C', b'G', b'T') | (b'C', b'G', b'C') | (b'C', b'G', b'A') | (b'C', b'G', b'G') => b'R',
        (b'A', b'T', b'T') | (b'A', b'T', b'C') | (b'A', b'T', b'A') => b'I',
        (b'A', b'T', b'G') => b'M',
        (b'A', b'C', b'T') | (b'A', b'C', b'C') | (b'A', b'C', b'A') | (b'A', b'C', b'G') => b'T',
        (b'A', b'A', b'T') | (b'A', b'A', b'C') => b'N',
        (b'A', b'A', b'A') | (b'A', b'A', b'G') => b'K',
        (b'A', b'G', b'T') | (b'A', b'G', b'C') => b'S',
        (b'A', b'G', b'A') | (b'A', b'G', b'G') => b'R',
        (b'G', b'T', b'T') | (b'G', b'T', b'C') | (b'G', b'T', b'A') | (b'G', b'T', b'G') => b'V',
        (b'G', b'C', b'T') | (b'G', b'C', b'C') | (b'G', b'C', b'A') | (b'G', b'C', b'G') => b'A',
        (b'G', b'A', b'T') | (b'G', b'A', b'C') => b'D',
        (b'G', b'A', b'A') | (b'G', b'A', b'G') => b'E',
        (b'G', b'G', b'T') | (b'G', b'G', b'C') | (b'G', b'G', b'A') | (b'G', b'G', b'G') => b'G',
        _ => b'X',
    }
}

/// Write a dual-payload file containing two serialized RoaringBitmaps
/// (DNA and AA), each prefixed by an 8-byte length in little-endian.
fn write_bitmap_pair(path: &str, dna: &RoaringBitmap, aa: &RoaringBitmap) -> PyResult<()> {
    let mut dna_buf = Vec::new();
    dna.serialize_into(&mut dna_buf)
        .map_err(|e| PyIOError::new_err(format!("Failed to serialize DNA bitmap: {}", e)))?;
    let mut aa_buf = Vec::new();
    aa.serialize_into(&mut aa_buf)
        .map_err(|e| PyIOError::new_err(format!("Failed to serialize AA bitmap: {}", e)))?;

    let mut file = File::create(path)
        .map_err(|e| PyIOError::new_err(format!("Failed to create {}: {}", path, e)))?;

    let dna_len = dna_buf.len() as u64;
    file.write_all(&dna_len.to_le_bytes())
        .map_err(|e| PyIOError::new_err(format!("Failed to write DNA length: {}", e)))?;
    file.write_all(&dna_buf)
        .map_err(|e| PyIOError::new_err(format!("Failed to write DNA bitmap: {}", e)))?;

    let aa_len = aa_buf.len() as u64;
    file.write_all(&aa_len.to_le_bytes())
        .map_err(|e| PyIOError::new_err(format!("Failed to write AA length: {}", e)))?;
    file.write_all(&aa_buf)
        .map_err(|e| PyIOError::new_err(format!("Failed to write AA bitmap: {}", e)))?;

    Ok(())
}

/// Read a dual-payload file and return (dna_bitmap, aa_bitmap).
fn read_bitmap_pair(path: &str) -> PyResult<(RoaringBitmap, RoaringBitmap)> {
    let mut file = File::open(path)
        .map_err(|e| PyIOError::new_err(format!("Failed to open {}: {}", path, e)))?;

    let mut dna_len_buf = [0u8; 8];
    file.read_exact(&mut dna_len_buf)
        .map_err(|e| PyIOError::new_err(format!("Failed to read DNA bitmap length: {}", e)))?;
    let dna_len = u64::from_le_bytes(dna_len_buf) as usize;

    let mut dna_buf = vec![0u8; dna_len];
    file.read_exact(&mut dna_buf)
        .map_err(|e| PyIOError::new_err(format!("Failed to read DNA bitmap: {}", e)))?;
    let dna = RoaringBitmap::deserialize_from(&mut Cursor::new(dna_buf))
        .map_err(|e| PyIOError::new_err(format!("Failed to deserialize DNA bitmap: {}", e)))?;

    let mut aa_len_buf = [0u8; 8];
    file.read_exact(&mut aa_len_buf)
        .map_err(|e| PyIOError::new_err(format!("Failed to read AA bitmap length: {}", e)))?;
    let aa_len = u64::from_le_bytes(aa_len_buf) as usize;

    let mut aa_buf = vec![0u8; aa_len];
    file.read_exact(&mut aa_buf)
        .map_err(|e| PyIOError::new_err(format!("Failed to read AA bitmap: {}", e)))?;
    let aa = RoaringBitmap::deserialize_from(&mut Cursor::new(aa_buf))
        .map_err(|e| PyIOError::new_err(format!("Failed to deserialize AA bitmap: {}", e)))?;

    Ok((dna, aa))
}

/// Select one bitmap from a dual-payload file based on mode.
fn select_bitmap(path: &str, mode: &str) -> PyResult<RoaringBitmap> {
    let (dna, aa) = read_bitmap_pair(path)?;
    match mode {
        "dna" | "DNA" => Ok(dna),
        "aa" | "AA" => Ok(aa),
        _ => Err(PyValueError::new_err(format!(
            "Invalid mode '{}': expected 'dna' or 'aa'", mode
        ))),
    }
}

/// Generate a dual-payload FracMinHash sketch and save it as two compressed
/// Roaring Bitmaps (DNA + AA translation in 3 forward frames).
///
/// FracMinHash sub-samples the k-mer space deterministically. Instead of storing
/// every k-mer in a genome, it only stores hashes that fall below a certain threshold.
/// This reduces the sketch size by a factor of `scaled` while preserving distance metrics.
///
/// DNA k-mers are hashed with ntHash, a rolling hash specialised for nucleotide k-mers:
/// it derives each k-mer's canonical (strand-independent) hash from the previous one in
/// O(1) time, which is substantially faster than re-hashing every k-mer byte by byte.
///
/// # Arguments
/// * `filepath` - Path to the input FASTA file.
/// * `out_filepath` - Path to save the dual-payload sketch file.
/// * `kmer_size` - The length of DNA k-mers to extract (e.g., 21, 31).
///   AA k-mer size is derived automatically as max(3, kmer_size / 3).
/// * `scaled` - The scale factor (e.g., 1000 means keep 1 in 1000 k-mers).
/// * `_threads` - Number of threads (currently unused, reserved for future Rayon integration).
#[pyfunction]
fn sketch(filepath: &str, out_filepath: &str, kmer_size: usize, scaled: u32, nthreads: usize) -> PyResult<()> {
    let aa_kmer_size = std::cmp::max(3, kmer_size / 3);
    let max_hash = (u32::MAX as f64 / scaled as f64) as u32;

    // Collect sequences first: needletail's reader holds a file handle and is not Send,
    // so it cannot cross thread boundaries. We pay one sequential pass to gather owned
    // byte vectors, then fan out the heavy hashing work in parallel.
    let mut sequences: Vec<Vec<u8>> = Vec::new();
    let mut reader = parse_fastx_file(filepath)
        .map_err(|e| PyIOError::new_err(format!("Failed to open {}: {}", filepath, e)))?;
    while let Some(record) = reader.next() {
        let seqrec = record.map_err(|e| PyValueError::new_err(e.to_string()))?;
        sequences.push(seqrec.seq().to_vec());
    }

    // Process sequences in parallel. Each sequence produces its own pair of bitmaps;
    // the reduce step merges them with bitwise OR, which is correct for FracMinHash sets.
    let pool = get_pool(nthreads);
    let (dna_bm, aa_bm) = pool.install(|| {
        sequences.par_iter().map(|seq| {
            let mut dna = RoaringBitmap::new();
            let mut aa  = RoaringBitmap::new();
            let mut upper: Vec<u8> = Vec::new(); // reusable uppercased copy for ntHash

            // DNA k-mer sketching.
            //
            // ntHash only accepts A/C/G/T: any other byte (including lowercase a/c/g/t)
            // makes it panic, and 'N' hashes to a meaningless constant. We therefore
            // uppercase the sequence and roll ntHash over each maximal A/C/G/T stretch.
            // Any k-mer spanning an ambiguous base is skipped, mirroring (and extending)
            // the previous "skip k-mers containing N" behaviour, while soft-masked
            // lowercase bases are folded back in through the uppercasing.
            if seq.len() >= kmer_size {
                upper.clear();
                upper.extend(seq.iter().map(|b| b.to_ascii_uppercase()));
                let len = upper.len();
                let mut start = 0usize;
                for i in 0..=len {
                    let valid = i < len && matches!(upper[i], b'A' | b'C' | b'G' | b'T');
                    if !valid {
                        if i - start >= kmer_size {
                            if let Ok(iter) = NtHashIterator::new(&upper[start..i], kmer_size) {
                                for hash in iter {
                                    let h = hash as u32;
                                    if h <= max_hash {
                                        dna.insert(h);
                                    }
                                }
                            }
                        }
                        start = i + 1;
                    }
                }
            }

            // AA k-mer sketching via 3-frame translation.
            //
            // Translation and hashing are fused: rather than materialising the three
            // translated frames and re-hashing every window, we roll a polynomial hash
            // over the amino-acid stream as codons are decoded — O(1) per residue with
            // no per-sequence allocation. As before, each frame is translated only up to
            // its first stop codon, and any window containing an unknown residue (X) is
            // skipped (an X resets the rolling window).
            if seq.len() >= 3 {
                let high = AA_HASH_BASE.wrapping_pow((aa_kmer_size - 1) as u32);
                let mut ring = vec![0u8; aa_kmer_size]; // last aa_kmer_size residues
                for offset in 0..3 {
                    let mut pos = 0usize;    // ring-buffer write cursor
                    let mut filled = 0usize; // residues accumulated since the last reset
                    let mut h: u64 = 0;      // rolling hash of the current window
                    let mut i = offset;
                    while i + 3 <= seq.len() {
                        let residue = dna_to_aa(&seq[i..i + 3]);
                        i += 3;
                        if residue == b'*' {
                            break; // stop codon ends this frame
                        }
                        if residue == b'X' {
                            // An unknown residue cannot belong to any emitted k-mer:
                            // drop the partial window and start fresh after it.
                            pos = 0;
                            filled = 0;
                            h = 0;
                            continue;
                        }
                        if filled < aa_kmer_size {
                            h = h.wrapping_mul(AA_HASH_BASE).wrapping_add(residue as u64);
                            ring[pos] = residue;
                            filled += 1;
                        } else {
                            let out = ring[pos] as u64;
                            h = h.wrapping_sub(out.wrapping_mul(high));
                            h = h.wrapping_mul(AA_HASH_BASE).wrapping_add(residue as u64);
                            ring[pos] = residue;
                        }
                        pos += 1;
                        if pos == aa_kmer_size {
                            pos = 0;
                        }
                        if filled == aa_kmer_size {
                            let hh = mix64(h) as u32;
                            if hh <= max_hash {
                                aa.insert(hh);
                            }
                        }
                    }
                }
            }

            (dna, aa)
        })
        .reduce(
            || (RoaringBitmap::new(), RoaringBitmap::new()),
            |(mut d1, mut a1), (d2, a2)| { d1 |= d2; a1 |= a2; (d1, a1) },
        )
    });

    write_bitmap_pair(out_filepath, &dna_bm, &aa_bm)
}

/// Compute Containment-based Average Nucleotide/Aminoacid Identity (ANI/AAI)
/// between a focus sketch and a list of target sketches.
///
/// Because we use FracMinHash, we can estimate ANI/AAI purely mathematically 
/// without needing full alignments.
/// The `mode` parameter selects the DNA or AA bitmap from the dual-payload file.
/// When mode is "aa", the effective k-mer size is derived as max(3, kmer_size / 3).
#[pyfunction]
fn containment_ani(focus: &str, targets: Vec<String>, kmer_size: usize, mode: &str) -> PyResult<HashMap<String, f64>> {
    // In AA mode, the effective k-mer size is derived from the DNA k-mer size:
    // a k-mer of N bp translates to N/3 AA residues.
    let eff_kmer = if matches!(mode, "aa" | "AA") {
        std::cmp::max(3, kmer_size / 3)
    } else {
        kmer_size
    };

    let focus_bm = select_bitmap(focus, mode)?;
    let focus_len = focus_bm.len() as f64;
    let mut results = HashMap::new();

    if focus_len == 0.0 {
        for t in &targets { results.insert(t.clone(), 1.0); }
        return Ok(results);
    }

    for target in targets.iter() {
        let target_bm = select_bitmap(target, mode)?;
        let intersection = focus_bm.intersection_len(&target_bm) as f64;
        let containment = intersection / focus_len;

        let ani = if containment > 0.0 {
            1.0 + (1.0 / eff_kmer as f64) * containment.ln()
        } else {
            0.0
        };

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
///
/// The `mode` parameter selects which bitmap (dna/aa) to operate on.
/// In "dna" mode, the parent's AA slot stores the union of children's AA bitmaps
/// (so AA data propagates upward through DNA-based levels for later AA mode sweeps).
/// In "aa" mode, only AA bitmaps are intersected/stripped; DNA passes through unchanged.
#[pyfunction]
fn build_delta_tree(sketches_list: &str, out_tree: &str, is_flat: bool, mode: &str) -> PyResult<()> {
    let list_file = File::open(sketches_list)
        .map_err(|e| PyIOError::new_err(format!("Cannot open sketches list: {}", e)))?;
    let reader = BufReader::new(list_file);

    // Select which bitmap to use for Core/Delta operations: "dna" or "aa"
    let primary = match mode {
        "dna" | "DNA" => true,
        "aa" | "AA" => false,
        _ => return Err(PyValueError::new_err(
            format!("Invalid mode '{}': expected 'dna' or 'aa'", mode)
        )),
    };

    let mut child_paths: Vec<String> = Vec::new();
    let mut child_primary: Vec<RoaringBitmap> = Vec::new(); // DNA or AA (depending on mode)
    let mut child_secondary: Vec<RoaringBitmap> = Vec::new(); // the other bitmap
    let mut core_bitmap = RoaringBitmap::new();
    let mut union_bitmap = RoaringBitmap::new(); // Union of the secondary bitmap (used in DNA mode)
    let mut first = true;

    // PASS 1 (Bottom-Up Step): Load children and compute Core/Union
    for line in reader.lines() {
        let filepath = line.map_err(|e| PyIOError::new_err(e.to_string()))?;
        let filepath = filepath.trim().to_string();
        if filepath.is_empty() { continue; }

        let (dna, aa) = read_bitmap_pair(&filepath)?;
        let (p_bm, s_bm) = if primary { (dna, aa) } else { (aa, dna) };

        if first {
            core_bitmap = p_bm.clone();
            first = false;
        } else {
            core_bitmap &= &p_bm;
        }
        // Union of the secondary bitmap (meaningful in DNA mode)
        union_bitmap |= &s_bm;

        child_paths.push(filepath);
        child_primary.push(p_bm);
        child_secondary.push(s_bm);
    }

    // PASS 2 (Top-Down Step): Strip Core from children's primary bitmap
    if !is_flat {
        for (i, path) in child_paths.iter().enumerate() {
            let mut p_bm = child_primary[i].clone();
            p_bm -= &core_bitmap;

            let s_bm = &child_secondary[i];

            // Write back: primary is stripped, secondary is unchanged
            if primary {
                write_bitmap_pair(path, &p_bm, s_bm)?;
            } else {
                write_bitmap_pair(path, s_bm, &p_bm)?;
            }
        }
    }

    // PASS 3: Save the parent node
    // For DNA mode: Core in DNA slot, Union in AA slot
    // For AA mode: read existing parent (if any) for DNA, use Core for AA
    if primary {
        write_bitmap_pair(out_tree, &core_bitmap, &union_bitmap)?;
    } else {
        // Preserve existing DNA bitmap if the parent file already exists
        let dna_bitmap = if Path::new(out_tree).exists() {
            let (existing_dna, _) = read_bitmap_pair(out_tree)?;
            existing_dna
        } else {
            RoaringBitmap::new()
        };
        write_bitmap_pair(out_tree, &dna_bitmap, &core_bitmap)?;
    }

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
/// ## Pruning
///
/// At each node all of its children are scored upfront.  The children are sorted by
/// distance and only those within `best_distance * (1 + uncertainty/100)` of the closest
/// sibling are enqueued.  This gives logarithmic traversal instead of exhaustive visits.
///
/// `theta` (pruning_threshold) is an additional absolute containment floor: any node whose
/// accumulated containment fraction (accumulated_score / query_len) falls below theta is
/// pruned together with its entire subtree.
///
/// ## Bitmap loading
///
/// Each node's bitmap is loaded exactly once — when its parent evaluates all children.
/// The queue carries the pre-computed accumulated score for the node it points to, so
/// no bitmap is re-read when a node is popped.
///
/// The `mode` parameter selects the DNA or AA bitmap from the dual-payload files.
#[pyfunction]
fn accumulator_search(
    query_sketch: &str,
    tree_root: &str,
    tree_topology: HashMap<String, Vec<(String, String)>>,
    kmer_size: usize,
    theta: f64,
    uncertainty: f64,
    mode: &str,
) -> PyResult<HashMap<String, HashMap<String, f64>>> {

    let eff_kmer = if matches!(mode, "aa" | "AA") {
        std::cmp::max(3, kmer_size / 3)
    } else {
        kmer_size
    };

    let query_bm = select_bitmap(query_sketch, mode)?;
    let query_len = query_bm.len() as f64;

    let mut profiles: HashMap<String, HashMap<String, f64>> = HashMap::new();

    if query_len == 0.0 || !Path::new(tree_root).exists() {
        return Ok(profiles);
    }

    // Inline helper: accumulated k-mer intersection count → ANI distance in [0, 1]
    let score_to_distance = |accumulated: f64| -> f64 {
        let containment = accumulated / query_len;
        let ani = if containment > 0.0 {
            1.0 + (1.0 / eff_kmer as f64) * containment.ln()
        } else {
            0.0
        };
        if ani <= 0.0 { 1.0 } else if ani >= 1.0 { 0.0 } else { 1.0 - ani }
    };

    // Pre-compute the root's accumulated score so every queue entry always carries the
    // final score for the node it represents (parent score + this node's contribution).
    let root_bm = select_bitmap(tree_root, mode)?;
    let root_accumulated = query_bm.intersection_len(&root_bm) as f64;

    // Queue: (node_path, accumulated_score_for_this_node, level_name)
    let mut queue: VecDeque<(String, f64, String)> = VecDeque::new();
    queue.push_back((tree_root.to_string(), root_accumulated, "db".to_string()));

    while let Some((node_path, my_accumulated, level_name)) = queue.pop_front() {
        // Record this node (skip the artificial "db" root above all kingdoms)
        if level_name != "db" {
            if my_accumulated / query_len < theta {
                continue;  // absolute floor — prune this node and don't expand its children
            }
            let distance = score_to_distance(my_accumulated);
            profiles.entry(level_name.clone())
                .or_insert_with(HashMap::new)
                .insert(node_path.clone(), distance);
        }

        let Some(children) = tree_topology.get(&node_path) else { continue; };

        // Score every child to decide which branches are worth pursuing
        let mut candidates: Vec<(String, String, f64, f64)> = Vec::new(); // (path, level, distance, accumulated)

        for (child_path, next_level) in children {
            if !Path::new(child_path).exists() {
                continue;
            }
            let child_bm = select_bitmap(child_path, mode)?;
            let child_accumulated = my_accumulated + query_bm.intersection_len(&child_bm) as f64;

            // Absolute containment floor: prune child and its subtree
            if child_accumulated / query_len < theta {
                continue;
            }

            candidates.push((
                child_path.clone(),
                next_level.clone(),
                score_to_distance(child_accumulated),
                child_accumulated,
            ));
        }

        if candidates.is_empty() {
            continue;
        }

        // Sort by distance ascending so candidates[0] is the closest child
        candidates.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));

        // Uncertainty cutoff: keep every child within best_distance * (1 + uncertainty/100).
        // This is a relative expansion, so a 50% uncertainty keeps all siblings up to 1.5×
        // the closest distance — naturally narrower at lower levels where distances are small.
        let best_distance = candidates[0].2;
        let cutoff = best_distance * (1.0 + uncertainty / 100.0);

        for (child_path, next_level, child_distance, child_accumulated) in candidates {
            if child_distance <= cutoff {
                queue.push_back((child_path, child_accumulated, next_level));
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
/// The `mode` parameter selects the bitmap (dna/aa) for rebalancing.
/// The other bitmap passes through unchanged.
#[pyfunction]
fn update_delta_tree(new_sketch: &str, species_node_path: &str, sibling_sketches: Vec<String>, mode: &str) -> PyResult<()> {
    let primary = matches!(mode, "dna" | "DNA");
    let (new_dna, new_aa) = read_bitmap_pair(new_sketch)?;
    let (species_dna, species_aa) = read_bitmap_pair(species_node_path)?;

    let (mut new_primary, species_primary, species_other, new_other) = if primary {
        (new_dna, species_dna, species_aa, new_aa)
    } else {
        (new_aa, species_aa, species_dna, new_dna)
    };

    let mut new_core = species_primary.clone();
    new_core &= &new_primary;

    let mut lost_core = species_primary.clone();
    lost_core -= &new_core;

    if !lost_core.is_empty() {
        for sibling_path in &sibling_sketches {
            let (sib_dna, sib_aa) = read_bitmap_pair(sibling_path)?;
            let (mut sib_bm, sib_other) = if primary {
                (sib_dna, sib_aa)
            } else {
                (sib_aa, sib_dna)
            };
            sib_bm |= &lost_core;
            if primary {
                write_bitmap_pair(sibling_path, &sib_bm, &sib_other)?;
            } else {
                write_bitmap_pair(sibling_path, &sib_other, &sib_bm)?;
            }
        }
        if primary {
            write_bitmap_pair(species_node_path, &new_core, &species_other)?;
        } else {
            write_bitmap_pair(species_node_path, &species_other, &new_core)?;
        }
    }

    new_primary -= &new_core;

    if primary {
        write_bitmap_pair(new_sketch, &new_primary, &new_other)?;
    } else {
        write_bitmap_pair(new_sketch, &new_other, &new_primary)?;
    }

    Ok(())
}

/// Read a serialized dual-payload sketch file and return the cardinality of the
/// selected bitmap (DNA or AA).
///
/// Cardinality is the exact number of elements (subsampled k-mer hashes) stored
/// in the FracMinHash sketch. In the Delta-SBT architecture this replaces the
/// obsolete "density" metric (ratio of set bits to total bits) that was only
/// meaningful for fixed-size Bloom filters.
#[pyfunction]
fn sketch_cardinality(sketch_path: &str, mode: &str) -> PyResult<u64> {
    let bitmap = select_bitmap(sketch_path, mode)?;
    Ok(bitmap.len())
}

/// The Python Module Definition.
/// 
/// This macro creates the entry points that `pyo3` and `maturin` will compile 
/// into the Python extension module. The name of the function `deltatree` MUST 
/// match the `lib.name` setting in the `Cargo.toml`.
#[pymodule]
fn deltatree(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sketch, m)?)?;
    m.add_function(wrap_pyfunction!(sketch_cardinality, m)?)?;
    m.add_function(wrap_pyfunction!(containment_ani, m)?)?;
    m.add_function(wrap_pyfunction!(accumulator_search, m)?)?;
    m.add_function(wrap_pyfunction!(build_delta_tree, m)?)?;
    m.add_function(wrap_pyfunction!(update_delta_tree, m)?)?;
    Ok(())
}
