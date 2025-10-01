//! Suffix-array helpers using rust-bio (bio 1.6.x).
//!
//! We originally considered using the FM-index (compressed suffix array with
//! Burrows–Wheeler transform) because it provides O(m) substring search time,
//! independent of the reference length n. This is why aligners such as BWA
//! and Bowtie rely on FM-index: they need to search billions of bases quickly.
//!
//! However, in this pipeline we always process **genome bins** of limited size
//! (≈100–200 kb). For these window sizes, suffix array binary search is
//! effectively just as fast: O(m log n), with log n ≈ 18 for n = 200,000.
//! This means only ~18 string comparisons per query, which is negligible
//! compared to the cost of I/O and alignment.
//!
//! ✅ In practice: using suffix arrays directly is simpler, avoids API churn in
//! rust-bio, and does not hurt performance for our bin-based design. FM-index
//! would only offer meaningful gains if we searched across the entire genome
//! at once (~3 Gb). Therefore, suffix array search is the pragmatic choice here.

use bio::data_structures::suffix_array::suffix_array;

/// Workspace holding the reference window and its suffix array.
pub struct SuffixWorkspace {
    pub text: Vec<u8>,   // the reference window sequence
    pub sa: Vec<usize>,  // suffix array over `text`
}

impl SuffixWorkspace {
    /// Build the suffix array over a reference window.
    pub fn build(window: &[u8]) -> Self {
        // Copy the window and add sentinel
        let mut text = window.to_vec();
        text.push(0); // sentinel symbol, smaller than A/C/G/T

        // Build suffix array
        let sa = suffix_array(&text);

        Self { text, sa }
    }


    /// Return all positions where `pat` occurs exactly in `text`.
    /// Positions are 0-based offsets within the window.
    pub fn exact_positions(&self, pat: &[u8]) -> Vec<usize> {
        if pat.is_empty() { return vec![]; }

        let start = self.lower_bound(pat);
        if start >= self.sa.len() { return vec![]; }

        let end = self.upper_bound(pat, start);
        self.sa[start..end]
            .iter()
            .copied()
            .filter(|&pos| pos + pat.len() <= self.text.len() - 1) // exclude sentinel
            .collect()
    }

    /// Binary search for the first suffix >= `pat` (lexicographically by prefix `pat`).
    fn lower_bound(&self, pat: &[u8]) -> usize {
        let (mut lo, mut hi) = (0usize, self.sa.len());
        while lo < hi {
            let mid = (lo + hi) / 2;
            let cmp = prefix_cmp(&self.text[self.sa[mid]..], pat);
            if cmp < 0 {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        lo
    }

    /// Binary search for the first suffix > `pat`, starting from `lo_start`.
    fn upper_bound(&self, pat: &[u8], lo_start: usize) -> usize {
        let (mut lo, mut hi) = (lo_start, self.sa.len());
        while lo < hi {
            let mid = (lo + hi) / 2;
            let cmp = prefix_cmp(&self.text[self.sa[mid]..], pat);
            if cmp <= 0 {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        lo
    }

    /// Check if the suffix starting at `pos` has `pat` as a prefix.
    #[inline]
    fn suffix_has_prefix(&self, pos: usize, pat: &[u8]) -> bool {
        let suf = &self.text[pos..];
        if suf.len() < pat.len() { return false; }
        suf[..pat.len()] == *pat
    }
}

/// Lexicographic comparison of `s` against prefix `pat`.
/// Returns:
///   -1 if `s`'s first `pat.len()` chars are lexicographically smaller than `pat`
///    0 if `s` starts with `pat`
///   +1 if `s`'s first `pat.len()` chars are greater than `pat`
#[inline]
fn prefix_cmp(s: &[u8], pat: &[u8]) -> i32 {
    let n = pat.len().min(s.len());
    for i in 0..n {
        if s[i] < pat[i] { return -1; }
        if s[i] > pat[i] { return 1; }
    }
    if s.len() < pat.len() { -1 } else { 0 }
}
