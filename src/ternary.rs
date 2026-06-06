//! Ternary type {-1, 0, +1} with full arithmetic over Z/3Z.
//!
//! Addition is wrap-around: 1 + 1 = -1 (mod 3), -1 + -1 = 1 (mod 3).
//! Multiplication follows sign rules. This forms a field (Z/3Z, +, ×).

use serde::{Deserialize, Serialize};

/// The ternary digit: negative (-1), zero (0), or positive (+1).
///
/// This type models elements of the field Z/3Z, where:
/// - `Neg` represents -1 ≡ 2 (mod 3)
/// - `Zero` represents 0
/// - `Pos` represents 1
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Ternary {
    Neg = -1,
    Zero = 0,
    Pos = 1,
}

impl Ternary {
    /// Convert from i32 to Ternary (mod 3).
    /// Maps: 0 → Zero, 1 → Pos, 2 → Neg (since 2 ≡ -1 mod 3).
    pub fn from_i32(v: i32) -> Self {
        match ((v % 3) + 3) % 3 {
            0 => Ternary::Zero,
            1 => Ternary::Pos,
            2 => Ternary::Neg,
            _ => unreachable!(),
        }
    }

    /// Convert Ternary to i32.
    pub fn to_i32(self) -> i32 {
        self as i32
    }

    /// Negation: -x. In Z/3Z: -(−1) = 1, −0 = 0, −1 = −1 (i.e., 2 mod 3).
    /// Actually: in Z/3Z, the additive inverse of 1 is 2 ≡ -1, and of 2 is 1.
    pub fn neg(self) -> Self {
        match self {
            Ternary::Neg => Ternary::Pos,
            Ternary::Zero => Ternary::Zero,
            Ternary::Pos => Ternary::Neg,
        }
    }

    /// Addition in Z/3Z (wrap-around).
    pub fn add(self, other: Self) -> Self {
        Self::from_i32(self.to_i32() + other.to_i32())
    }

    /// Subtraction in Z/3Z.
    pub fn sub(self, other: Self) -> Self {
        self.add(other.neg())
    }

    /// Multiplication in Z/3Z.
    pub fn mul(self, other: Self) -> Self {
        Self::from_i32(self.to_i32() * other.to_i32())
    }

    /// Multiplicative inverse. Only Pos and Neg have inverses.
    /// Pos⁻¹ = Pos, Neg⁻¹ = Neg (since (-1)×(-1) = 1).
    pub fn inv(self) -> Option<Self> {
        match self {
            Ternary::Zero => None,
            Ternary::Pos => Some(Ternary::Pos),
            Ternary::Neg => Some(Ternary::Neg),
        }
    }

    /// All variants in order.
    pub fn all() -> [Self; 3] {
        [Ternary::Neg, Ternary::Zero, Ternary::Pos]
    }

    /// Is this zero?
    pub fn is_zero(self) -> bool {
        self == Ternary::Zero
    }
}

impl std::fmt::Display for Ternary {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Ternary::Neg => write!(f, "-1"),
            Ternary::Zero => write!(f, "0"),
            Ternary::Pos => write!(f, "+1"),
        }
    }
}

impl Default for Ternary {
    fn default() -> Self {
        Ternary::Zero
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_all_additions() {
        // All 3×3 = 9 combinations
        // -1 + -1 = 1 (wrap), -1 + 0 = -1, -1 + 1 = 0
        assert_eq!(Ternary::Neg.add(Ternary::Neg), Ternary::Pos); // -2 ≡ 1 mod 3
        assert_eq!(Ternary::Neg.add(Ternary::Zero), Ternary::Neg);
        assert_eq!(Ternary::Neg.add(Ternary::Pos), Ternary::Zero);
        // 0 + anything = anything
        assert_eq!(Ternary::Zero.add(Ternary::Neg), Ternary::Neg);
        assert_eq!(Ternary::Zero.add(Ternary::Zero), Ternary::Zero);
        assert_eq!(Ternary::Zero.add(Ternary::Pos), Ternary::Pos);
        // 1 + -1 = 0, 1 + 0 = 1, 1 + 1 = -1 (wrap)
        assert_eq!(Ternary::Pos.add(Ternary::Neg), Ternary::Zero);
        assert_eq!(Ternary::Pos.add(Ternary::Zero), Ternary::Pos);
        assert_eq!(Ternary::Pos.add(Ternary::Pos), Ternary::Neg); // 2 ≡ -1 mod 3
    }

    #[test]
    fn test_all_multiplications() {
        assert_eq!(Ternary::Neg.mul(Ternary::Neg), Ternary::Pos); // (-1)(-1) = 1
        assert_eq!(Ternary::Neg.mul(Ternary::Zero), Ternary::Zero);
        assert_eq!(Ternary::Neg.mul(Ternary::Pos), Ternary::Neg); // (-1)(1) = -1
        assert_eq!(Ternary::Zero.mul(Ternary::Neg), Ternary::Zero);
        assert_eq!(Ternary::Zero.mul(Ternary::Zero), Ternary::Zero);
        assert_eq!(Ternary::Zero.mul(Ternary::Pos), Ternary::Zero);
        assert_eq!(Ternary::Pos.mul(Ternary::Neg), Ternary::Neg);
        assert_eq!(Ternary::Pos.mul(Ternary::Zero), Ternary::Zero);
        assert_eq!(Ternary::Pos.mul(Ternary::Pos), Ternary::Pos);
    }

    #[test]
    fn test_negation_involution() {
        for t in Ternary::all() {
            assert_eq!(t.neg().neg(), t, "negation of negation should be identity for {:?}", t);
        }
    }

    #[test]
    fn test_zero_additive_identity() {
        for t in Ternary::all() {
            assert_eq!(t.add(Ternary::Zero), t);
            assert_eq!(Ternary::Zero.add(t), t);
        }
    }

    #[test]
    fn test_one_multiplicative_identity() {
        for t in Ternary::all() {
            assert_eq!(t.mul(Ternary::Pos), t);
            assert_eq!(Ternary::Pos.mul(t), t);
        }
    }

    #[test]
    fn test_commutativity() {
        for a in Ternary::all() {
            for b in Ternary::all() {
                assert_eq!(a.add(b), b.add(a));
                assert_eq!(a.mul(b), b.mul(a));
            }
        }
    }

    #[test]
    fn test_associativity() {
        for a in Ternary::all() {
            for b in Ternary::all() {
                for c in Ternary::all() {
                    assert_eq!(a.add(b).add(c), a.add(b.add(c)));
                    assert_eq!(a.mul(b).mul(c), a.mul(b.mul(c)));
                }
            }
        }
    }

    #[test]
    fn test_additive_inverse() {
        for t in Ternary::all() {
            assert_eq!(t.add(t.neg()), Ternary::Zero);
        }
    }

    #[test]
    fn test_multiplicative_inverse() {
        assert_eq!(Ternary::Pos.inv(), Some(Ternary::Pos));
        assert_eq!(Ternary::Neg.inv(), Some(Ternary::Neg));
        assert_eq!(Ternary::Zero.inv(), None);
        // x * x.inv() = 1
        for t in Ternary::all() {
            if let Some(inv) = t.inv() {
                assert_eq!(t.mul(inv), Ternary::Pos);
            }
        }
    }

    #[test]
    fn test_from_i32_roundtrip() {
        for t in Ternary::all() {
            assert_eq!(Ternary::from_i32(t.to_i32()), t);
        }
    }

    #[test]
    fn test_distributivity() {
        for a in Ternary::all() {
            for b in Ternary::all() {
                for c in Ternary::all() {
                    assert_eq!(a.mul(b.add(c)), a.mul(b).add(a.mul(c)));
                }
            }
        }
    }

    #[test]
    fn test_subtraction() {
        for a in Ternary::all() {
            for b in Ternary::all() {
                assert_eq!(a.sub(b), a.add(b.neg()));
            }
        }
    }

    #[test]
    fn test_display() {
        assert_eq!(format!("{}", Ternary::Neg), "-1");
        assert_eq!(format!("{}", Ternary::Zero), "0");
        assert_eq!(format!("{}", Ternary::Pos), "+1");
    }
}
