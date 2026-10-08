/*! Deterministic random-number mixing shared by search and resampling. */

/// SplitMix64 step used to derive independent deterministic seeds.
pub(crate) fn splitmix64(mut value: u64) -> u64 {
    value = value.wrapping_add(0x9e37_79b9_7f4a_7c15);
    avalanche64(value)
}

/// Finalizes a 64-bit generator state without advancing it.
pub(crate) fn avalanche64(mut value: u64) -> u64 {
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}
