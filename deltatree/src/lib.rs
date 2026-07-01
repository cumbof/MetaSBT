//! DeltaTree: High-performance Rust backend for MetaSBT v2.0
//! 
//! This module implements the core algorithms for the MetaSBT Sequence Bloom Tree:
//! every internal node stores the union of its subtree's sub-sampled k-mers. The union bounds
//! which subtrees can contain a query (reachability), while sibling clades are ranked by the query's
//! IDF-weighted (discriminative) union containment: k-mers shared by every sibling carry no signal
//! and are down-weighted to nothing, so a query is not funnelled into the largest clade merely
//! because its union is a denser sample of the universal k-mer space. FracMinHash bounds the hash
//! space, keeping even the root union compact as a Roaring bitmap.
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

/// Target size (in bases) of a DNA sketching chunk. A single-genome `sketch` splits each
/// sequence into overlapping chunks of about this size so that one long contig can be
/// hashed across multiple threads. It is large enough that per-chunk overhead is
/// negligible, yet small enough that a multi-megabase contig still yields many chunks.
const DNA_CHUNK: usize = 262_144;

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
/// genetic code (NCBI translation table 1; identical to the bacterial table 11
/// for every codon, since 11 differs only in alternative start codons, which we
/// do not special-case). Returns b'X' for unknown codons.
///
/// One table is used for every kingdom on purpose. The amino acid sketch is a
/// similarity transform, not a biological annotation: as long as the same code
/// is applied when building the database and when querying, FracMinHash
/// containment between any two genomes stays self-consistent and comparable.
/// The standard code is correct for the vast majority of bacteria, archaea,
/// fungi, other eukaryotes, and host-translated viruses; clades with reassigned
/// codons (e.g. Mycoplasma TGA=W, ciliate TAA/TAG=Q) are merely sketched a
/// little more sparsely, not incorrectly relative to each other.
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
        (b'T', b'G', b'A') => b'*', // stop in the standard genetic code
        (b'T', b'G', b'G') => b'W',
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

/// Sketch the DNA k-mers of `seq` into `dna`, reusing the caller's `upper` buffer.
///
/// ntHash only accepts A/C/G/T: any other byte (including lowercase a/c/g/t) makes it
/// panic, and 'N' hashes to a meaningless constant. We therefore uppercase the sequence
/// and roll ntHash over each maximal A/C/G/T stretch. Any k-mer spanning an ambiguous
/// base is skipped, mirroring (and extending) the previous "skip k-mers containing N"
/// behaviour, while soft-masked lowercase bases are folded back in through uppercasing.
///
/// Because each k-mer's ntHash value depends only on the k-mer itself, this is safe to
/// call on overlapping slices of a longer sequence: as long as consecutive slices overlap
/// by at least `kmer_size - 1`, the union over the slices equals the result for the whole
/// sequence. `sketch` exploits this to chunk a long contig across threads.
fn sketch_dna_into(seq: &[u8], kmer_size: usize, max_hash: u32, dna: &mut RoaringBitmap, upper: &mut Vec<u8>) {
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
}

/// Sketch the amino acid k-mers of `seq` into `aa`, reusing the caller's `ring` buffer
/// (length `aa_kmer_size`) and the precomputed `aa_high` rolling-hash weight.
///
/// Translation and hashing are fused: rather than materialising the three translated
/// frames and re-hashing every window, we roll a polynomial hash over the amino-acid
/// stream as codons are decoded — O(1) per residue with no allocation. Every ORF in each
/// frame is sketched: a stop codon (*) ends the current ORF and resets the rolling window,
/// and any window containing an unknown residue (X) is likewise skipped (an X resets the
/// rolling window). Unlike the DNA path this is position-dependent (frame phase and stop
/// codons), so it must be run over a whole sequence rather than chunked.
fn sketch_aa_into(seq: &[u8], aa_kmer_size: usize, max_hash: u32, aa_high: u64, aa: &mut RoaringBitmap, ring: &mut Vec<u8>) {
    if seq.len() >= 3 {
        for offset in 0..3 {
            let mut pos = 0usize;    // ring-buffer write cursor
            let mut filled = 0usize; // residues accumulated since the last reset
            let mut h: u64 = 0;      // rolling hash of the current window
            let mut i = offset;
            while i + 3 <= seq.len() {
                let residue = dna_to_aa(&seq[i..i + 3]);
                i += 3;
                if residue == b'*' {
                    // Stop codon: end of the current ORF, not the end of the frame.
                    // Reset the rolling window and keep translating the next ORF so every
                    // ORF in the frame is sketched (each emitted k-mer is still
                    // FracMinHash-subsampled by the `max_hash` gate below).
                    pos = 0;
                    filled = 0;
                    h = 0;
                    continue;
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
                    h = h.wrapping_sub(out.wrapping_mul(aa_high));
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
}

/// Sketch one whole sequence into both the DNA and AA bitmaps. Used by `sketch_many`,
/// which walks each genome's sequences serially; `sketch` instead calls the DNA and AA
/// halves separately so it can chunk the DNA pass across threads.
fn sketch_sequence_into(
    seq: &[u8],
    kmer_size: usize,
    aa_kmer_size: usize,
    max_hash: u32,
    aa_high: u64,
    dna: &mut RoaringBitmap,
    aa: &mut RoaringBitmap,
    upper: &mut Vec<u8>,
    ring: &mut Vec<u8>,
) {
    sketch_dna_into(seq, kmer_size, max_hash, dna, upper);
    sketch_aa_into(seq, aa_kmer_size, max_hash, aa_high, aa, ring);
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

    let pool = get_pool(nthreads);
    let high = AA_HASH_BASE.wrapping_pow((aa_kmer_size - 1) as u32);

    // DNA pass. Parallelising over whole sequences leaves a genome that is one huge
    // contig stuck on a single thread, so we instead fan out over overlapping chunks of
    // each sequence. Consecutive chunks overlap by kmer_size-1 bases, which guarantees
    // every k-mer is fully contained in some chunk; since each k-mer's ntHash value
    // depends only on the k-mer, the union over chunks equals the whole-sequence result.
    let overlap = kmer_size.saturating_sub(1);
    let step = DNA_CHUNK.saturating_sub(overlap).max(1);
    let mut chunks: Vec<(usize, usize, usize)> = Vec::new(); // (seq index, start, end)
    for (idx, seq) in sequences.iter().enumerate() {
        if seq.len() < kmer_size {
            continue;
        }
        let mut start = 0usize;
        loop {
            let end = std::cmp::min(start + DNA_CHUNK, seq.len());
            chunks.push((idx, start, end));
            if end == seq.len() {
                break;
            }
            start += step;
        }
    }

    let (dna_bm, aa_bm) = pool.install(|| {
        // DNA: one task per overlapping chunk, merged with bitwise OR.
        let dna_bm = chunks
            .par_iter()
            .map_init(Vec::<u8>::new, |upper, &(idx, s, e)| {
                let mut dna = RoaringBitmap::new();
                sketch_dna_into(&sequences[idx][s..e], kmer_size, max_hash, &mut dna, upper);
                dna
            })
            .reduce(RoaringBitmap::new, |mut a, b| {
                a |= b;
                a
            });

        // AA: one task per sequence (the AA pass is position-dependent and cannot be
        // chunked because frame phase and ORF boundaries depend on absolute position).
        let aa_bm = sequences
            .par_iter()
            .map_init(|| vec![0u8; aa_kmer_size], |ring, seq| {
                let mut aa = RoaringBitmap::new();
                sketch_aa_into(seq, aa_kmer_size, max_hash, high, &mut aa, ring);
                aa
            })
            .reduce(RoaringBitmap::new, |mut a, b| {
                a |= b;
                a
            });

        (dna_bm, aa_bm)
    });

    write_bitmap_pair(out_filepath, &dna_bm, &aa_bm)
}

/// Sketch a batch of genomes, parallelising across genomes rather than within each one.
///
/// `jobs` is a list of `(input_fasta, output_sketch)` pairs. Each genome is read and
/// sketched independently on the Rayon pool: parallelism is across genomes, while each
/// genome's own sequences are processed serially (and its scratch buffers reused). For
/// an `index` run with many genomes this keeps every core busy with far less per-genome
/// overhead than launching a separate Python worker process per genome, and lets Rayon
/// load-balance the whole batch with work stealing. Returns the output paths written.
#[pyfunction]
fn sketch_many(jobs: Vec<(String, String)>, kmer_size: usize, scaled: u32, nthreads: usize) -> PyResult<Vec<String>> {
    let aa_kmer_size = std::cmp::max(3, kmer_size / 3);
    let max_hash = (u32::MAX as f64 / scaled as f64) as u32;
    let high = AA_HASH_BASE.wrapping_pow((aa_kmer_size - 1) as u32);

    let pool = get_pool(nthreads);
    pool.install(|| {
        jobs.par_iter()
            .map(|(in_path, out_path)| {
                // Read this genome's sequences (needletail's reader is not Send, but it
                // lives and dies entirely within this one task, so that is fine here).
                let mut reader = parse_fastx_file(in_path)
                    .map_err(|e| PyIOError::new_err(format!("Failed to open {}: {}", in_path, e)))?;

                let mut dna = RoaringBitmap::new();
                let mut aa = RoaringBitmap::new();
                let mut upper: Vec<u8> = Vec::new();
                let mut ring = vec![0u8; aa_kmer_size];

                while let Some(record) = reader.next() {
                    let seqrec = record.map_err(|e| PyValueError::new_err(e.to_string()))?;
                    sketch_sequence_into(
                        &seqrec.seq(), kmer_size, aa_kmer_size, max_hash, high,
                        &mut dna, &mut aa, &mut upper, &mut ring,
                    );
                }

                write_bitmap_pair(out_path, &dna, &aa)?;
                Ok(out_path.clone())
            })
            .collect()
    })
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

/// Build the Sequence Bloom Tree structure from a list of children sketches.
///
/// Every internal node stores the **union** of its subtree — the set of all sub-sampled
/// k-mer hashes carried by any descendant — exactly as a classic Sequence Bloom Tree does.
/// Because a FracMinHash sketch is a bottom/mod sketch, the union of children sketches is
/// itself a valid FracMinHash sketch of the combined k-mer set, so containment-ANI math
/// stays exact when the node is queried. Unions compose associatively, so a single bottom-up
/// sweep yields the full subtree union at every level: a genus node's union is the union of
/// its species' unions, which are the unions of their genomes, and so on.
///
/// The union (never the intersection) is the correct primitive for membership search. The
/// intersection — the k-mers common to *all* members — shrinks as a clade grows and balloons
/// for singleton clades, which makes `|Q ∩ node| / |Q|` systematically favour singleton/rare
/// lineages and funnels unrelated queries into them. The union has no such cardinality bias:
/// a larger union does not inflate an unrelated query's containment, because FracMinHash only
/// counts shared hashes.
///
/// Storage is bounded and cheap: FracMinHash caps the hash space at `u32::MAX / scaled`, so
/// even the root union (all DB k-mers) fits in a few MB as a Roaring bitmap. There is no need
/// to delta-encode the unions; an "OR down the path" delta only reconstructs a quantity that
/// *grows* toward the leaves (an intersection/Core), and a subtree union shrinks toward the
/// leaves, so that scheme does not apply here.
///
/// The `mode` parameter selects which bitmap (dna/aa) to operate on. In "dna" mode both the
/// DNA and AA unions of the children are written (the AA union propagates upward for the
/// AA-phase search). In "aa" mode only the AA union is (re)written; the DNA slot is preserved.
#[pyfunction]
fn build_delta_tree(sketches_list: &str, out_tree: &str, mode: &str) -> PyResult<()> {
    let list_file = File::open(sketches_list)
        .map_err(|e| PyIOError::new_err(format!("Cannot open sketches list: {}", e)))?;
    let reader = BufReader::new(list_file);

    // Select which bitmap is the primary union target for this sweep: "dna" or "aa"
    let primary = match mode {
        "dna" | "DNA" => true,
        "aa" | "AA" => false,
        _ => return Err(PyValueError::new_err(
            format!("Invalid mode '{}': expected 'dna' or 'aa'", mode)
        )),
    };

    let mut primary_union = RoaringBitmap::new();   // union of the primary bitmap
    let mut secondary_union = RoaringBitmap::new();  // union of the secondary bitmap

    // Bottom-up step: load children and accumulate the union of both payloads
    for line in reader.lines() {
        let filepath = line.map_err(|e| PyIOError::new_err(e.to_string()))?;
        let filepath = filepath.trim().to_string();
        if filepath.is_empty() { continue; }

        let (dna, aa) = read_bitmap_pair(&filepath)?;
        let (p_bm, s_bm) = if primary { (dna, aa) } else { (aa, dna) };

        primary_union |= &p_bm;
        secondary_union |= &s_bm;
    }

    // Save the parent node.
    // DNA mode: DNA union in the DNA slot, AA union in the AA slot.
    // AA mode: rewrite only the AA union; preserve the existing DNA slot.
    if primary {
        write_bitmap_pair(out_tree, &primary_union, &secondary_union)?;
    } else {
        let dna_bitmap = if Path::new(out_tree).exists() {
            let (existing_dna, _) = read_bitmap_pair(out_tree)?;
            existing_dna
        } else {
            RoaringBitmap::new()
        };
        write_bitmap_pair(out_tree, &dna_bitmap, &primary_union)?;
    }

    Ok(())
}

/// The Sequence Bloom Tree search algorithm with IDF-weighted (discriminative) ranking.
///
/// Every node stores the union of its subtree (see `build_delta_tree`): a species node holds
/// the k-mers of all its genomes, a genus node the k-mers of all its species, and so on up to
/// the kingdom. Genome leaves hold their own full sketch.
///
/// ## Two roles for two quantities
///
/// A subtree union is the right primitive for *reachability* but a biased one for *ranking*.
/// Raw union containment `|Q ∩ union| / |Q|` is monotone in clade size: a larger union can only add
/// k-mers, so at the deep, conserved (AA) levels — where every large clade's union saturates the
/// universal k-mer space — the biggest clade tends to win regardless of true membership, and an
/// unrelated query is funnelled into it. This search therefore separates the two roles:
///
/// * **Reachability / pruning** (could the query be in this subtree at all): the raw union is kept
///   as an optimistic upper bound. Because every member's k-mers are a subset of the union, the
///   containment of the query in any single member can never exceed its containment in the union;
///   a union whose containment falls below `theta` therefore cannot hold a good match, so the
///   child and its whole subtree are pruned. This keeps the traversal logarithmic.
///
/// * **Ranking** (which child the query is most *like*): the surviving siblings are ranked by the
///   query's *IDF-weighted* union containment. For the current node's children, each query k-mer `h`
///   is weighted by `idf(h) = ln(N / df(h))`, where `N` is the number of surviving children and
///   `df(h)` is how many of their unions contain `h`. A k-mer shared by every sibling (a universal,
///   conserved k-mer) gets `idf = ln(1) = 0` and contributes nothing; a k-mer specific to one clade
///   gets full weight. A child's score is the sum of `idf` over the query k-mers it contains, divided
///   by the total `idf` mass the query could reach across all siblings. This is size-invariant: a big
///   clade's union wins raw containment only through the universal k-mers, which now weigh nothing, so
///   ranking follows the query's *discriminative* overlap. When there is no signal to separate the
///   siblings (a single survivor, or every shared k-mer is universal so the total mass is zero), the
///   ranking falls back to raw union containment.
///
/// ## Recorded distance
///
/// The distance recorded for each retained node is its *raw* union-containment ANI/AAI (query-to-union),
/// the same quantity as before, so it stays comparable to the cluster boundary for confidence scoring.
/// The IDF weighting only decides *which* siblings are pursued, not the magnitude that is reported.
///
/// ## Pruning / beam
///
/// Surviving children are sorted by their IDF-weighted ranking distance and only those within
/// `best_distance * (1 + uncertainty/100)` of the closest sibling are enqueued. When the closest
/// sibling is an exact match (`best_distance == 0`) the relative window would collapse to zero, so
/// an absolute fallback of `uncertainty/100` keeps near matches in play.
///
/// `theta` (pruning_threshold) is the absolute raw-union-containment floor described above.
///
/// The `mode` parameter selects the DNA or AA bitmap from the dual-payload files (union bound,
/// IDF weighting and recorded distances are all computed in that payload).
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

    // Inline helper: containment in [0, 1] → ANI/AAI distance in [0, 1]
    let dist_from_containment = |containment: f64| -> f64 {
        let ani = if containment > 0.0 {
            1.0 + (1.0 / eff_kmer as f64) * containment.ln()
        } else {
            0.0
        };
        if ani <= 0.0 { 1.0 } else if ani >= 1.0 { 0.0 } else { 1.0 - ani }
    };

    // A surviving child of the node currently being expanded: its query∩union bitmap (reused for
    // the IDF pass) plus the raw union-containment distance that will be reported if it is retained.
    struct Cand {
        path: String,
        level: String,
        inter: RoaringBitmap, // query ∩ child union
        plain_distance: f64,  // raw union-containment ANI/AAI distance (reported)
    }

    // Queue: (node_path, recorded_distance_for_this_node, level_name). The entry node (the
    // search root, or an anchored seed) carries the artificial level "db" and is never recorded;
    // its distance is unused, so a placeholder is pushed for it.
    let mut queue: VecDeque<(String, f64, String)> = VecDeque::new();
    queue.push_back((tree_root.to_string(), 0.0, "db".to_string()));

    while let Some((node_path, node_distance, level_name)) = queue.pop_front() {
        // Record this node (skip the artificial entry root above the first scored level).
        // The distance was already computed (and theta-checked) by the parent that enqueued it.
        if level_name != "db" {
            profiles.entry(level_name.clone())
                .or_insert_with(HashMap::new)
                .insert(node_path.clone(), node_distance);
        }

        let Some(children) = tree_topology.get(&node_path) else { continue; };

        // First pass: load each child's union, keep the query∩union bitmap, and apply the
        // reachability bound. The subtree union is an optimistic upper bound on the containment
        // achievable by any single member, so a union below theta cannot hold a good match — prune
        // the child and its whole subtree.
        let mut cands: Vec<Cand> = Vec::new();
        for (child_path, next_level) in children {
            if !Path::new(child_path).exists() {
                continue;
            }
            let child_union_bm = select_bitmap(child_path, mode)?;
            let inter = &query_bm & &child_union_bm;
            let union_score = inter.len() as f64;
            if union_score / query_len < theta {
                continue;
            }
            cands.push(Cand {
                path: child_path.clone(),
                level: next_level.clone(),
                plain_distance: dist_from_containment(union_score / query_len),
                inter,
            });
        }

        if cands.is_empty() {
            continue;
        }

        // Second pass: inverse-document-frequency weighting over the surviving siblings.
        // df(h) = number of surviving children whose union contains query k-mer h. A k-mer shared
        // by every sibling (universal / conserved) cannot tell them apart, so idf = ln(N/df) sends
        // it to ~0; a k-mer specific to one clade keeps full weight. This strips the clade-size bias
        // out of union containment — a big clade's union wins raw containment only through the
        // universal k-mers, which now weigh nothing — so ranking follows the discriminative overlap.
        let n = cands.len() as f64;
        let mut df: HashMap<u32, u32> = HashMap::new();
        for c in &cands {
            for h in &c.inter {
                *df.entry(h).or_insert(0) += 1;
            }
        }
        let mut idf: HashMap<u32, f64> = HashMap::with_capacity(df.len());
        let mut weighted_total = 0.0f64; // idf mass the query could reach across all siblings
        for (&h, &d) in &df {
            let w = (n / d as f64).ln();
            idf.insert(h, w);
            weighted_total += w;
        }

        // Ranking distance per child. With no discriminative signal (a single survivor, or every
        // shared k-mer is universal so weighted_total is 0) fall back to raw union containment.
        // (path, level, rank_distance, plain_distance)
        let mut ranked: Vec<(String, String, f64, f64)> = Vec::with_capacity(cands.len());
        for c in &cands {
            let rank_distance = if weighted_total > 0.0 {
                let mut w = 0.0f64;
                for h in &c.inter {
                    w += idf[&h];
                }
                dist_from_containment(w / weighted_total)
            } else {
                c.plain_distance
            };
            ranked.push((c.path.clone(), c.level.clone(), rank_distance, c.plain_distance));
        }

        // Sort by ranking distance ascending so ranked[0] is the most discriminative child.
        ranked.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));

        // Uncertainty cutoff: keep every child within best_distance * (1 + uncertainty/100).
        // This is a relative expansion, so a 50% uncertainty keeps all siblings up to 1.5×
        // the closest distance — naturally narrower at lower levels where distances are small.
        // If the best sibling is an exact match the relative window collapses to zero, so fall
        // back to an absolute window of uncertainty/100 to avoid pruning every alternative.
        let best_distance = ranked[0].2;
        let cutoff = if best_distance > 0.0 {
            best_distance * (1.0 + uncertainty / 100.0)
        } else {
            uncertainty / 100.0
        };

        for (child_path, next_level, rank_distance, plain_distance) in ranked {
            if rank_distance <= cutoff {
                // Rank/beam by the discriminative distance, but record the raw union-containment
                // distance so the reported value stays comparable to the level's boundary.
                queue.push_back((child_path, plain_distance, next_level));
            }
        }
    }

    Ok(profiles)
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
    m.add_function(wrap_pyfunction!(sketch_many, m)?)?;
    m.add_function(wrap_pyfunction!(sketch_cardinality, m)?)?;
    m.add_function(wrap_pyfunction!(containment_ani, m)?)?;
    m.add_function(wrap_pyfunction!(accumulator_search, m)?)?;
    m.add_function(wrap_pyfunction!(build_delta_tree, m)?)?;
    Ok(())
}
