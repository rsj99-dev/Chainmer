#![forbid(unsafe_code)]

//! EvoSeed：面向直系同源候选召回与排序的种子库。
//!
//! Seed Hit 仅产生候选，不构成同源或直系同源判定。

pub mod cli;
pub mod core;
pub mod fasta;
pub mod index;
pub mod scoring;
pub mod search;
pub mod seeds;
pub mod subseqhash2;

pub use core::{anchors, evoseed_keys, evoseed_seeds, EvoParams, PositionedSeed, SeedTrack};
pub use index::SeedIndex;
pub use scoring::{contrastive_weights, specificity_weights, weighted_score};
pub use search::CandidateScore;

/// 包版本；索引元数据与此常量保持一致。
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
/// 兼容历史输出：128-bit 存储宽度，但有效状态最多 64 bit。
pub const HASH_SCHEME: &str = "dual-splitmix64/64-effective-in-u128";
pub const INDEX_FORMAT: &str = "evoseed-index/3";
