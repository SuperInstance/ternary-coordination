//! Tutorial: ternary-coordination — Balanced ternary {-1,0,+1} as Z/3Z coordination algebra
//!
//! Ternary matrices, consensus convergence, spectral radius, balanced ternary i32 conversion.
//! Exhaustive axiom verification for ring/field properties.

use ternary_coordination::{
    Ternary, TernaryRing, TernaryMatrix, ConsensusState, BalancedTernary,
    spectrum::spectral_radius,
    consensus::fully_connected_mixing_matrix,
};

fn main() {
    println!("=== Ternary Coordination Tutorial ===\n");

    // Part 1: Ternary values
    println!("Part 1: Balanced ternary values");
    for v in Ternary::all() {
        println!("  {:?} = {}", v, v.to_i32());
    }
    println!();

    // Part 2: Ring arithmetic
    println!("Part 2: Ternary ring arithmetic");
    println!("  (-1) + (+1) = {:?}", TernaryRing::add(Ternary::Neg, Ternary::Pos));
    println!("  (+1) × (+1) = {:?}", TernaryRing::mul(Ternary::Pos, Ternary::Pos));
    println!("  (-1) × (-1) = {:?}", TernaryRing::mul(Ternary::Neg, Ternary::Neg));
    println!("  -(-1) = {:?}", TernaryRing::neg(Ternary::Neg));
    println!();

    // Part 3: Axiom verification
    println!("Part 3: Ring axiom verification");
    println!("  Add associativity: {}", TernaryRing::verify_add_associativity());
    println!("  Mul associativity: {}", TernaryRing::verify_mul_associativity());
    println!("  Add commutativity: {}", TernaryRing::verify_add_commutativity());
    println!("  Mul commutativity: {}", TernaryRing::verify_mul_commutativity());
    println!("  Closure: {}", TernaryRing::verify_closure());
    println!();

    // Part 4: Ternary matrices
    println!("Part 4: Ternary matrices");
    let mat = TernaryMatrix::from_i32_slice(3, 3, &[
         1,  0, -1,
         0,  1,  0,
        -1,  0,  1,
    ]).unwrap();
    println!("  3×3 matrix created");
    let sr = spectral_radius(&mat);
    println!("  Spectral radius: {:.4}", sr);
    println!();

    // Part 5: Consensus
    println!("Part 5: Ternary consensus");
    let votes = vec![Ternary::Pos, Ternary::Pos, Ternary::Zero, Ternary::Neg, Ternary::Pos];
    let mut state = ConsensusState::new(votes);
    let (pos, zero, neg) = state.tally();
    println!("  Votes: {} pos, {} zero, {} neg", pos, zero, neg);
    println!("  Majority: {:?}", state.majority());
    println!("  Is unanimous: {}", state.is_unanimous());
    println!();

    // Part 6: Consensus convergence
    println!("Part 6: Consensus convergence");
    let mixing = fully_connected_mixing_matrix(5);
    let rounds = state.run_to_convergence(&mixing, 100);
    println!("  Converged in: {:?} rounds", rounds);
    println!();

    // Part 7: Balanced ternary ↔ i32
    println!("Part 7: Balanced ternary conversion");
    let bt = BalancedTernary::from_i32(8);
    println!("  8 → balanced ternary → back to i32: {}", bt.to_i32());
    let bt_neg = BalancedTernary::from_i32(-5);
    println!("  -5 → balanced ternary → back to i32: {}", bt_neg.to_i32());
}
