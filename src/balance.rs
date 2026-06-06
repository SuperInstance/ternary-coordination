//! Balanced ternary representation.
//!
//! Convert i32 → balanced ternary digits {-1, 0, +1}.
//! Balanced ternary is a non-standard positional numeral system where digits
//! range from -1 to +1. It was used in the Setun computer (1958, Moscow State University)
//! and studied by Donald Knuth and others.

use serde::{Deserialize, Serialize};
use crate::ternary::Ternary;

/// A balanced ternary number: any i32 represented as Σ digits[i] × 3^i.
/// digits[0] is the least significant trit.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct BalancedTernary {
    pub digits: Vec<Ternary>,
}

impl BalancedTernary {
    /// Create from digits (least significant first).
    pub fn from_digits(digits: Vec<Ternary>) -> Self {
        let mut bt = BalancedTernary { digits };
        bt.trim();
        bt
    }

    /// Create zero.
    pub fn zero() -> Self {
        BalancedTernary { digits: vec![Ternary::Zero] }
    }

    /// Create from an i32 value.
    pub fn from_i32(mut value: i32) -> Self {
        if value == 0 {
            return Self::zero();
        }
        let mut digits = Vec::new();
        while value != 0 {
            let mut r = value % 3;
            value /= 3;
            if r > 1 {
                r -= 3;
                value += 1;
            } else if r < -1 {
                r += 3;
                value -= 1;
            }
            digits.push(Ternary::from_i32(r));
        }
        BalancedTernary { digits }
    }

    /// Convert back to i32.
    pub fn to_i32(&self) -> i32 {
        let mut result: i64 = 0;
        let mut power: i64 = 1;
        for &digit in &self.digits {
            result += digit.to_i32() as i64 * power;
            power *= 3;
        }
        result as i32
    }

    /// Add two balanced ternary numbers using the carry truth table.
    pub fn add(&self, other: &BalancedTernary) -> BalancedTernary {
        let max_len = self.digits.len().max(other.digits.len()) + 2;
        let mut digits = Vec::with_capacity(max_len);
        let mut carry = Ternary::Zero;
        
        for i in 0..max_len {
            let a = self.digits.get(i).copied().unwrap_or(Ternary::Zero);
            let b = other.digits.get(i).copied().unwrap_or(Ternary::Zero);
            let raw = a.to_i32() + b.to_i32() + carry.to_i32();
            
            let (digit, new_carry) = match raw {
                -3 => (Ternary::Zero, Ternary::Neg),
                -2 => (Ternary::Pos, Ternary::Neg),
                -1 => (Ternary::Neg, Ternary::Zero),
                0  => (Ternary::Zero, Ternary::Zero),
                1  => (Ternary::Pos, Ternary::Zero),
                2  => (Ternary::Neg, Ternary::Pos),
                3  => (Ternary::Zero, Ternary::Pos),
                _ => unreachable!("ternary sum out of range: {}", raw),
            };
            digits.push(digit);
            carry = new_carry;
        }
        
        BalancedTernary::from_digits(digits)
    }

    /// Negate a balanced ternary number.
    pub fn neg(&self) -> BalancedTernary {
        BalancedTernary {
            digits: self.digits.iter().map(|&d| d.neg()).collect(),
        }
    }

    /// Subtract two balanced ternary numbers.
    pub fn sub(&self, other: &BalancedTernary) -> BalancedTernary {
        self.add(&other.neg())
    }

    /// Multiply by a single trit.
    pub fn mul_trit(&self, t: Ternary) -> BalancedTernary {
        BalancedTernary {
            digits: self.digits.iter().map(|&d| d.mul(t)).collect(),
        }
    }

    /// Remove trailing zeros (most significant).
    fn trim(&mut self) {
        while self.digits.len() > 1 && self.digits.last() == Some(&Ternary::Zero) {
            self.digits.pop();
        }
    }

    /// Number of trits.
    pub fn trit_count(&self) -> usize {
        self.digits.len()
    }

    /// Format as a string like "+0-" (most significant first).
    pub fn to_string_msb(&self) -> String {
        if self.digits.is_empty() {
            return "0".to_string();
        }
        self.digits.iter().rev().map(|d| match d {
            Ternary::Neg => '-',
            Ternary::Zero => '0',
            Ternary::Pos => '+',
        }).collect()
    }
}

impl std::fmt::Display for BalancedTernary {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_string_msb())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zero() {
        let bt = BalancedTernary::from_i32(0);
        assert_eq!(bt.digits, vec![Ternary::Zero]);
    }

    #[test]
    fn test_one() {
        let bt = BalancedTernary::from_i32(1);
        assert_eq!(bt.digits, vec![Ternary::Pos]);
    }

    #[test]
    fn test_neg_one() {
        let bt = BalancedTernary::from_i32(-1);
        assert_eq!(bt.digits, vec![Ternary::Neg]);
    }

    #[test]
    fn test_two() {
        // 2 = 1*3 + (-1) → digits: [Neg, Pos] (LSB first)
        let bt = BalancedTernary::from_i32(2);
        assert_eq!(bt.digits, vec![Ternary::Neg, Ternary::Pos]);
        assert_eq!(bt.to_i32(), 2);
    }

    #[test]
    fn test_three() {
        // 3 = 0 + 1*3 → digits: [Zero, Pos]
        let bt = BalancedTernary::from_i32(3);
        assert_eq!(bt.digits, vec![Ternary::Zero, Ternary::Pos]);
        assert_eq!(bt.to_i32(), 3);
    }

    #[test]
    fn test_four() {
        // 4 = 1 + 1*3 → digits: [Pos, Pos]
        let bt = BalancedTernary::from_i32(4);
        assert_eq!(bt.digits, vec![Ternary::Pos, Ternary::Pos]);
        assert_eq!(bt.to_i32(), 4);
    }

    #[test]
    fn test_five() {
        // 5: 5/3 = 1 rem 2, 2→-1 carry 1. 1+1=2. 2/3=0 rem 2, 2→-1 carry 1.
        // digits: [-1, -1, 1] = -1 + -3 + 9 = 5 ✓
        let bt = BalancedTernary::from_i32(5);
        assert_eq!(bt.digits, vec![Ternary::Neg, Ternary::Neg, Ternary::Pos]);
        assert_eq!(bt.to_i32(), 5);
    }

    #[test]
    fn test_negative_two() {
        // -2 = 1 + (-1)*3 = 1 - 3 = -2 → digits: [Pos, Neg]
        let bt = BalancedTernary::from_i32(-2);
        assert_eq!(bt.digits, vec![Ternary::Pos, Ternary::Neg]);
        assert_eq!(bt.to_i32(), -2);
    }

    #[test]
    fn test_roundtrip_small() {
        for i in -50i32..=50 {
            let bt = BalancedTernary::from_i32(i);
            assert_eq!(bt.to_i32(), i, "roundtrip failed for {}", i);
        }
    }

    #[test]
    fn test_roundtrip_large() {
        let values = [100, -100, 1000, -1000, 987654, -987654, 715827882, -715827882];
        for &v in &values {
            let bt = BalancedTernary::from_i32(v);
            assert_eq!(bt.to_i32(), v, "roundtrip failed for {}", v);
        }
    }

    #[test]
    fn test_add() {
        let a = BalancedTernary::from_i32(5);
        let b = BalancedTernary::from_i32(-3);
        let sum = a.add(&b);
        assert_eq!(sum.to_i32(), 2);
    }

    #[test]
    fn test_add_various() {
        for a in -20..=20 {
            for b in -20..=20 {
                let ba = BalancedTernary::from_i32(a);
                let bb = BalancedTernary::from_i32(b);
                assert_eq!(ba.add(&bb).to_i32(), a + b, "add failed: {} + {}", a, b);
            }
        }
    }

    #[test]
    fn test_neg() {
        let a = BalancedTernary::from_i32(7);
        let neg_a = a.neg();
        assert_eq!(neg_a.to_i32(), -7);
    }

    #[test]
    fn test_sub() {
        let a = BalancedTernary::from_i32(10);
        let b = BalancedTernary::from_i32(3);
        let diff = a.sub(&b);
        assert_eq!(diff.to_i32(), 7);
    }

    #[test]
    fn test_mul_trit() {
        let a = BalancedTernary::from_i32(5);
        let neg = a.mul_trit(Ternary::Neg);
        assert_eq!(neg.to_i32(), -5);
        
        let zero = a.mul_trit(Ternary::Zero);
        assert_eq!(zero.to_i32(), 0);
    }

    #[test]
    fn test_display() {
        let one = BalancedTernary::from_i32(1);
        assert_eq!(format!("{}", one), "+");
        
        let two = BalancedTernary::from_i32(2);
        assert_eq!(format!("{}", two), "+-");
        
        let zero = BalancedTernary::from_i32(0);
        assert_eq!(format!("{}", zero), "0");
    }

    #[test]
    fn test_add_commutative() {
        let a = BalancedTernary::from_i32(7);
        let b = BalancedTernary::from_i32(-4);
        assert_eq!(a.add(&b).to_i32(), b.add(&a).to_i32());
    }

    #[test]
    fn test_trit_count() {
        assert_eq!(BalancedTernary::from_i32(0).trit_count(), 1);
        assert_eq!(BalancedTernary::from_i32(1).trit_count(), 1);
        assert_eq!(BalancedTernary::from_i32(2).trit_count(), 2);
    }
}
