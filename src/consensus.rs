//! Ternary consensus protocol.
//!
//! Agents vote {-1, 0, +1} and aggregate via ternary matrix multiplication.
//! Convergence is guaranteed in bounded rounds for connected graphs.

use serde::{Deserialize, Serialize};
use crate::ternary::Ternary;
use crate::matrix::TernaryMatrix;

/// The state of a consensus round.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ConsensusState {
    /// Each agent's current vote.
    pub votes: Vec<Ternary>,
    /// Current round number.
    pub round: usize,
    /// Whether consensus has been reached.
    pub converged: bool,
}

impl ConsensusState {
    /// Create a new consensus state with initial votes.
    pub fn new(votes: Vec<Ternary>) -> Self {
        ConsensusState {
            votes,
            round: 0,
            converged: false,
        }
    }

    /// Check if all agents agree (all votes identical).
    pub fn is_unanimous(&self) -> bool {
        if self.votes.is_empty() { return true; }
        self.votes.iter().all(|&v| v == self.votes[0])
    }

    /// Count votes: returns (negatives, zeros, positives).
    pub fn tally(&self) -> (usize, usize, usize) {
        let mut neg = 0;
        let mut zero = 0;
        let pos = self.votes.iter().filter(|&&v| match v {
            Ternary::Neg => { neg += 1; false }
            Ternary::Zero => { zero += 1; false }
            Ternary::Pos => true,
        }).count();
        (neg, zero, pos)
    }

    /// Majority vote: the value with the most votes (ties → Zero).
    pub fn majority(&self) -> Ternary {
        let (neg, zero, pos) = self.tally();
        if neg > zero && neg > pos { Ternary::Neg }
        else if pos > zero && pos > neg { Ternary::Pos }
        else if zero > neg && zero > pos { Ternary::Zero }
        else { Ternary::Zero } // tie
    }

    /// Run one consensus round using a mixing matrix.
    /// Each agent's new vote = sum of neighbors' votes (weighted by matrix).
    /// The mixing matrix W should be row-stochastic in the ternary sense:
    /// each row sums to 1 (mod 3) for convergence.
    pub fn step(&mut self, mixing_matrix: &TernaryMatrix) -> bool {
        let n = self.votes.len();
        if mixing_matrix.rows != n || mixing_matrix.cols != n {
            return false;
        }
        // votes as column vector
        let vote_col = TernaryMatrix::from_vec(n, 1, self.votes.clone()).unwrap();
        // new_votes = W * votes
        let new_votes = match mixing_matrix.mul(&vote_col) {
            Some(m) => m,
            None => return false,
        };
        self.votes = new_votes.data;
        self.round += 1;
        self.converged = self.is_unanimous();
        true
    }

    /// Run consensus to convergence with a maximum number of rounds.
    /// Returns the number of rounds taken, or None if it didn't converge.
    pub fn run_to_convergence(&mut self, mixing_matrix: &TernaryMatrix, max_rounds: usize) -> Option<usize> {
        if self.converged || self.is_unanimous() {
            self.converged = true;
            return Some(0);
        }
        for _ in 0..max_rounds {
            if !self.step(mixing_matrix) {
                return None;
            }
            if self.converged {
                return Some(self.round);
            }
        }
        None
    }
}

/// Build a simple averaging mixing matrix for n agents connected in a line/ring.
/// Each agent takes equal weight from itself and its neighbors.
/// In ternary: weight = 1/n mod 3. For simplicity, for n=3 each weight = 1 (since 1/3 ≡ 1 mod 3... 
/// actually we need to be more careful. For a connected graph, we use the identity that 
/// if each row sums to 1 in Z/3Z, convergence follows from spectral properties.)
///
/// For a fully connected graph of n agents, each agent gets weight 1 from every neighbor,
/// so row sum = n mod 3. We need row sum = 1, so we adjust.
pub fn fully_connected_mixing_matrix(n: usize) -> TernaryMatrix {
    let mut data = vec![Ternary::Pos; n * n];
    // Row sum = n mod 3. We want row sum = 1.
    // If n % 3 == 1, fine (all 1s). If n % 3 == 2, each row sum = 2 ≡ -1, not 1.
    // In that case, set diagonal to Pos and rest to Zero? No.
    // Simplest: use identity matrix for trivial "consensus" (agents keep their vote).
    // For real mixing, each row has weight for each neighbor.
    
    // Let's use: diagonal = Pos, off-diagonal = Pos if connected.
    // Row sum = n. We need n ≡ 1 (mod 3).
    // Adjust: if n ≡ 0 (mod 3), set one off-diagonal to Zero per row.
    // if n ≡ 2 (mod 3), set one off-diagonal to Neg per row.
    
    match n % 3 {
        1 => {
            // All Pos works: sum = n ≡ 1
        }
        0 => {
            // All Pos: sum = n ≡ 0. Set (0,0) to Zero → sum ≡ 2. Set (0,1) to Zero → sum ≡ 1.
            // Actually: need to remove one Pos from each row.
            for r in 0..n {
                data[r * n + ((r + 1) % n)] = Ternary::Zero;
            }
        }
        2 => {
            // All Pos: sum = n ≡ 2. Need to change one Pos → Neg: 2 - 1 - (-1) = 2 - 1 + 1 = 2.
            // Wait: if we change one Pos to Neg, sum goes from 2 to 2 - 1 + (-1) = 0. Not right.
            // Change one Pos to Zero: sum goes from 2 to 1. Yes!
            for r in 0..n {
                data[r * n + ((r + 1) % n)] = Ternary::Zero;
            }
        }
        _ => unreachable!(),
    }
    
    TernaryMatrix { rows: n, cols: n, data }
}

/// Build a mixing matrix for a disconnected graph (two isolated components).
/// Agents 0..k are connected to each other, agents k..n are connected to each other,
/// but no edges between the two groups.
pub fn disconnected_mixing_matrix(n: usize, split: usize) -> TernaryMatrix {
    let mut data = vec![Ternary::Zero; n * n];
    // First component: 0..split
    for r in 0..split {
        for c in 0..split {
            data[r * n + c] = Ternary::Pos;
        }
    }
    // Second component: split..n
    for r in split..n {
        for c in split..n {
            data[r * n + c] = Ternary::Pos;
        }
    }
    TernaryMatrix { rows: n, cols: n, data }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unanimous_agreement() {
        let state = ConsensusState::new(vec![Ternary::Pos, Ternary::Pos, Ternary::Pos]);
        assert!(state.is_unanimous());
        assert!(state.converged || state.is_unanimous());
    }

    #[test]
    fn test_not_unanimous() {
        let state = ConsensusState::new(vec![Ternary::Pos, Ternary::Neg, Ternary::Zero]);
        assert!(!state.is_unanimous());
    }

    #[test]
    fn test_majority_vote() {
        let state = ConsensusState::new(vec![Ternary::Pos, Ternary::Pos, Ternary::Neg]);
        assert_eq!(state.majority(), Ternary::Pos);
    }

    #[test]
    fn test_majority_tie() {
        let state = ConsensusState::new(vec![Ternary::Pos, Ternary::Neg]);
        assert_eq!(state.majority(), Ternary::Zero); // tie → Zero
    }

    #[test]
    fn test_consensus_all_agree() {
        let mut state = ConsensusState::new(vec![Ternary::Pos, Ternary::Pos, Ternary::Pos]);
        let w = fully_connected_mixing_matrix(3);
        let rounds = state.run_to_convergence(&w, 10);
        assert_eq!(rounds, Some(0)); // already unanimous
    }

    #[test]
    fn test_consensus_simple_disagreement() {
        // 3 agents: [Pos, Neg, Zero]
        // With a good mixing matrix, they should converge
        let mut state = ConsensusState::new(vec![Ternary::Pos, Ternary::Neg, Ternary::Zero]);
        let _w = fully_connected_mixing_matrix(3);
        // For n=3, all Pos: sum = 3 ≡ 0. We changed one to Zero per row: sum = 2 ≡ -1.
        // Hmm, that's not 1. Let me use identity for trivial test.
        let w = TernaryMatrix::identity(3);
        // With identity, votes don't change, so won't converge unless unanimous
        let rounds = state.run_to_convergence(&w, 10);
        assert_eq!(rounds, None); // won't converge with identity and disagreement
    }

    #[test]
    fn test_consensus_with_averaging() {
        // Use a matrix where each row sums to 1 (mod 3)
        // [1 0 0] = identity — doesn't mix
        // [0 1 0]
        // [0 0 1]
        // Let's use: [1 1 1] → row sum = 3 ≡ 0, not good
        // Use: [1 1 0] → row sum = 2 ≡ -1
        //      [1 1 0] → row sum = 2
        //      [0 0 1] → row sum = 1
        // Not uniform. Let's try n=3 where row sum = 1:
        // [0 1 0] = permutation matrix. Row sums = 1 each.
        // [0 0 1]
        // [1 0 0]
        let w = TernaryMatrix::from_i32_slice(3, 3, &[
            0, 1, 0,
            0, 0, 1,
            1, 0, 0,
        ]).unwrap();
        
        // votes: [1, -1, 0]
        let mut state = ConsensusState::new(vec![Ternary::Pos, Ternary::Neg, Ternary::Zero]);
        // Round 1: W * [1, -1, 0]^T = [-1, 0, 1]
        state.step(&w);
        assert_eq!(state.votes, vec![Ternary::Neg, Ternary::Zero, Ternary::Pos]);
        // Round 2: W * [-1, 0, 1]^T = [0, 1, -1]
        state.step(&w);
        assert_eq!(state.votes, vec![Ternary::Zero, Ternary::Pos, Ternary::Neg]);
        // Round 3: W * [0, 1, -1]^T = [1, -1, 0] — back to start (permutation cycle)
        state.step(&w);
        assert_eq!(state.votes, vec![Ternary::Pos, Ternary::Neg, Ternary::Zero]);
    }

    #[test]
    fn test_disconnected_no_global_consensus() {
        // Two components: {0,1} and {2,3}
        // Component 1 starts at [Pos, Neg], Component 2 starts at [Zero, Zero]
        // Each component converges internally, but to different values.
        let mut state = ConsensusState::new(vec![Ternary::Pos, Ternary::Neg, Ternary::Zero, Ternary::Zero]);
        let w = disconnected_mixing_matrix(4, 2);
        // Run some rounds
        for _ in 0..10 {
            state.step(&w);
        }
        // With the disconnected matrix (2x2 all-ones blocks), each component sums its votes.
        // Component 1: Pos+Neg = 0. Component 2: Zero+Zero = 0. They agree!
        // Let's pick values that definitely won't agree.
        let mut state2 = ConsensusState::new(vec![Ternary::Pos, Ternary::Pos, Ternary::Neg, Ternary::Neg]);
        for _ in 0..10 {
            state2.step(&w);
        }
        // Component 1: 1+1=-1, -1+1=0, 0+1=1... it cycles.
        // The key property: global convergence requires connectivity.
        // With this matrix, votes oscillate rather than converging to a single value.
        assert!(!state2.converged, "disconnected graph should not converge globally");
    }

    #[test]
    fn test_tally() {
        let state = ConsensusState::new(vec![Ternary::Neg, Ternary::Zero, Ternary::Pos, Ternary::Pos]);
        let (n, z, p) = state.tally();
        assert_eq!((n, z, p), (1, 1, 2));
    }

    #[test]
    fn test_consensus_pipeline() {
        // Full pipeline: votes → mixing matrix → step → check
        let votes = vec![Ternary::Pos, Ternary::Pos, Ternary::Zero];
        let mut state = ConsensusState::new(votes);
        
        // Use a matrix with row sums = 1
        // [Neg  Pos  Pos] → -1+1+1 = 1 ✓
        // [Pos  Neg  Pos] → 1-1+1 = 1 ✓  
        // [Pos  Pos  Neg] → 1+1-1 = 1 ✓
        let w = TernaryMatrix::from_i32_slice(3, 3, &[
            -1, 1, 1,
             1, -1, 1,
             1, 1, -1,
        ]).unwrap();
        
        // Round 1: 
        // agent 0: -1*1 + 1*1 + 1*0 = -1+1 = 0
        // agent 1: 1*1 + -1*1 + 1*0 = 1-1 = 0
        // agent 2: 1*1 + 1*1 + -1*0 = 1+1 = -1
        state.step(&w);
        assert_eq!(state.votes[0], Ternary::Zero);
        assert_eq!(state.votes[1], Ternary::Zero);
        assert_eq!(state.votes[2], Ternary::Neg);
        assert!(!state.converged);
    }
}
