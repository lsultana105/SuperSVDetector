#[derive(Clone, Debug)]
pub struct Hotspot {
    pub local_pos: usize,
    pub reason: String,
}

pub fn hotspot_near_boundary(
    hotspots: &[Hotspot],
    bin_len: usize,
    margin: usize,
) -> bool {
    hotspots.iter().any(|h| {
        h.local_pos <= margin || h.local_pos >= bin_len.saturating_sub(margin)
    })
}
