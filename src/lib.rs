//! Ternary Coordination Algebra
//!
//! Formalizes {-1, 0, +1} as a coordination algebra with ternary matrices,
//! ternary consensus, and provable convergence.

pub mod ternary;
pub mod matrix;
pub mod consensus;
pub mod algebra;
pub mod spectrum;
pub mod balance;

pub use ternary::Ternary;
pub use matrix::TernaryMatrix;
pub use consensus::ConsensusState;
pub use algebra::TernaryRing;
pub use balance::BalancedTernary;
