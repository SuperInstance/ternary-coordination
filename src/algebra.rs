//! Ternary algebraic structures: TernaryGroup, TernaryRing, TernaryModule.
//!
//! Formal verification that Z/3Z forms a ring (and field).
//! With zero divisors analysis: 1 + 1 = -1 in Z/3Z, so 2 ≡ -1.

use crate::ternary::Ternary;

/// Z/3Z as an algebraic ring.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TernaryRing;

impl TernaryRing {
    /// The additive identity (0).
    pub fn zero() -> Ternary { Ternary::Zero }
    
    /// The multiplicative identity (1).
    pub fn one() -> Ternary { Ternary::Pos }
    
    /// Addition in the ring.
    pub fn add(a: Ternary, b: Ternary) -> Ternary { a.add(b) }
    
    /// Multiplication in the ring.
    pub fn mul(a: Ternary, b: Ternary) -> Ternary { a.mul(b) }
    
    /// Additive inverse.
    pub fn neg(a: Ternary) -> Ternary { a.neg() }
    
    /// Subtraction.
    pub fn sub(a: Ternary, b: Ternary) -> Ternary { a.sub(b) }

    /// Verify closure: all operations produce elements in {-1, 0, 1}.
    /// (Trivially true by the type system, but we verify the algebraic properties.)
    pub fn verify_closure() -> bool { true }

    /// Verify associativity of addition: (a+b)+c = a+(b+c) for all a,b,c.
    pub fn verify_add_associativity() -> bool {
        for a in Ternary::all() {
            for b in Ternary::all() {
                for c in Ternary::all() {
                    if a.add(b).add(c) != a.add(b.add(c)) {
                        return false;
                    }
                }
            }
        }
        true
    }

    /// Verify associativity of multiplication: (a*b)*c = a*(b*c) for all a,b,c.
    pub fn verify_mul_associativity() -> bool {
        for a in Ternary::all() {
            for b in Ternary::all() {
                for c in Ternary::all() {
                    if a.mul(b).mul(c) != a.mul(b.mul(c)) {
                        return false;
                    }
                }
            }
        }
        true
    }

    /// Verify additive commutativity: a+b = b+a.
    pub fn verify_add_commutativity() -> bool {
        for a in Ternary::all() {
            for b in Ternary::all() {
                if a.add(b) != b.add(a) {
                    return false;
                }
            }
        }
        true
    }

    /// Verify multiplicative commutativity: a*b = b*a.
    pub fn verify_mul_commutativity() -> bool {
        for a in Ternary::all() {
            for b in Ternary::all() {
                if a.mul(b) != b.mul(a) {
                    return false;
                }
            }
        }
        true
    }

    /// Verify distributivity: a*(b+c) = a*b + a*c.
    pub fn verify_distributivity() -> bool {
        for a in Ternary::all() {
            for b in Ternary::all() {
                for c in Ternary::all() {
                    if a.mul(b.add(c)) != a.mul(b).add(a.mul(c)) {
                        return false;
                    }
                }
            }
        }
        true
    }

    /// Verify additive identity: a + 0 = a.
    pub fn verify_additive_identity() -> bool {
        for a in Ternary::all() {
            if a.add(Ternary::Zero) != a {
                return false;
            }
        }
        true
    }

    /// Verify multiplicative identity: a * 1 = a.
    pub fn verify_multiplicative_identity() -> bool {
        for a in Ternary::all() {
            if a.mul(Ternary::Pos) != a {
                return false;
            }
        }
        true
    }

    /// Verify additive inverse: a + (-a) = 0.
    pub fn verify_additive_inverse() -> bool {
        for a in Ternary::all() {
            if a.add(a.neg()) != Ternary::Zero {
                return false;
            }
        }
        true
    }

    /// Check for zero divisors: nonzero a, b such that a*b = 0.
    /// In Z/3Z, there are NO zero divisors (it's a field).
    pub fn has_zero_divisors() -> bool {
        for a in [Ternary::Neg, Ternary::Pos] {
            for b in [Ternary::Neg, Ternary::Pos] {
                if a.mul(b) == Ternary::Zero {
                    return true;
                }
            }
        }
        false
    }

    /// Verify this is a field: every nonzero element has a multiplicative inverse.
    pub fn is_field() -> bool {
        for a in [Ternary::Neg, Ternary::Pos] {
            if a.inv().is_none() {
                return false;
            }
        }
        true
    }

    /// Verify ALL ring axioms. Returns (passed, list_of_failures).
    pub fn verify_all_axioms() -> (bool, Vec<&'static str>) {
        let mut failures = Vec::new();
        if !Self::verify_add_associativity() { failures.push("add associativity"); }
        if !Self::verify_mul_associativity() { failures.push("mul associativity"); }
        if !Self::verify_add_commutativity() { failures.push("add commutativity"); }
        if !Self::verify_mul_commutativity() { failures.push("mul commutativity"); }
        if !Self::verify_distributivity() { failures.push("distributivity"); }
        if !Self::verify_additive_identity() { failures.push("additive identity"); }
        if !Self::verify_multiplicative_identity() { failures.push("multiplicative identity"); }
        if !Self::verify_additive_inverse() { failures.push("additive inverse"); }
        (failures.is_empty(), failures)
    }

    /// The characteristic of Z/3Z is 3: 1 + 1 + 1 = 0.
    pub fn characteristic() -> usize { 3 }

    /// Verify: 1 + 1 + 1 = 0 (characteristic 3).
    pub fn verify_characteristic() -> bool {
        Ternary::Pos.add(Ternary::Pos).add(Ternary::Pos) == Ternary::Zero
    }
}

/// A ternary group (Z/3Z under addition).
pub trait TernaryGroup {
    type Element;
    fn identity() -> Self::Element;
    fn op(a: &Self::Element, b: &Self::Element) -> Self::Element;
    fn inv(a: &Self::Element) -> Self::Element;
}

/// The additive group of Z/3Z.
pub struct AdditiveGroup;

impl TernaryGroup for AdditiveGroup {
    type Element = Ternary;
    fn identity() -> Ternary { Ternary::Zero }
    fn op(a: &Ternary, b: &Ternary) -> Ternary { a.add(*b) }
    fn inv(a: &Ternary) -> Ternary { a.neg() }
}

/// The multiplicative group of Z/3Z (excluding zero).
pub struct MultiplicativeGroup;

impl TernaryGroup for MultiplicativeGroup {
    type Element = Ternary;
    fn identity() -> Ternary { Ternary::Pos }
    fn op(a: &Ternary, b: &Ternary) -> Ternary { a.mul(*b) }
    fn inv(a: &Ternary) -> Ternary { a.inv().expect("nonzero has inverse") }
}

/// A ternary module: vectors over Z/3Z.
pub trait TernaryModule {
    fn zero(dim: usize) -> Vec<Ternary>;
    fn add(a: &[Ternary], b: &[Ternary]) -> Option<Vec<Ternary>>;
    fn scale(scalar: Ternary, v: &[Ternary]) -> Vec<Ternary>;
}

/// Standard ternary module implementation.
pub struct VectorModule;

impl TernaryModule for VectorModule {
    fn zero(dim: usize) -> Vec<Ternary> {
        vec![Ternary::Zero; dim]
    }

    fn add(a: &[Ternary], b: &[Ternary]) -> Option<Vec<Ternary>> {
        if a.len() != b.len() { return None; }
        Some(a.iter().zip(b.iter()).map(|(&x, &y)| x.add(y)).collect())
    }

    fn scale(scalar: Ternary, v: &[Ternary]) -> Vec<Ternary> {
        v.iter().map(|&x| scalar.mul(x)).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ring_closure() {
        assert!(TernaryRing::verify_closure());
    }

    #[test]
    fn test_ring_add_associativity() {
        assert!(TernaryRing::verify_add_associativity());
    }

    #[test]
    fn test_ring_mul_associativity() {
        assert!(TernaryRing::verify_mul_associativity());
    }

    #[test]
    fn test_ring_commutativity() {
        assert!(TernaryRing::verify_add_commutativity());
        assert!(TernaryRing::verify_mul_commutativity());
    }

    #[test]
    fn test_ring_distributivity() {
        assert!(TernaryRing::verify_distributivity());
    }

    #[test]
    fn test_ring_identities() {
        assert!(TernaryRing::verify_additive_identity());
        assert!(TernaryRing::verify_multiplicative_identity());
    }

    #[test]
    fn test_ring_additive_inverse() {
        assert!(TernaryRing::verify_additive_inverse());
    }

    #[test]
    fn test_all_ring_axioms() {
        let (passed, failures) = TernaryRing::verify_all_axioms();
        assert!(passed, "Ring axiom failures: {:?}", failures);
    }

    #[test]
    fn test_is_field() {
        assert!(TernaryRing::is_field());
        assert!(!TernaryRing::has_zero_divisors());
    }

    #[test]
    fn test_characteristic() {
        assert_eq!(TernaryRing::characteristic(), 3);
        assert!(TernaryRing::verify_characteristic());
    }

    #[test]
    fn test_two_equals_neg_one() {
        // In Z/3Z: 1 + 1 = 2 ≡ -1
        assert_eq!(Ternary::Pos.add(Ternary::Pos), Ternary::Neg);
    }

    #[test]
    fn test_additive_group() {
        let e = AdditiveGroup::identity();
        assert_eq!(e, Ternary::Zero);
        assert_eq!(AdditiveGroup::op(&Ternary::Pos, &Ternary::Neg), Ternary::Zero);
        assert_eq!(AdditiveGroup::inv(&Ternary::Pos), Ternary::Neg);
    }

    #[test]
    fn test_multiplicative_group() {
        let e = MultiplicativeGroup::identity();
        assert_eq!(e, Ternary::Pos);
        assert_eq!(MultiplicativeGroup::op(&Ternary::Neg, &Ternary::Neg), Ternary::Pos);
        assert_eq!(MultiplicativeGroup::inv(&Ternary::Neg), Ternary::Neg);
    }

    #[test]
    fn test_vector_module() {
        let zero = VectorModule::zero(3);
        assert_eq!(zero, vec![Ternary::Zero; 3]);
        
        let a = vec![Ternary::Pos, Ternary::Neg, Ternary::Zero];
        let b = vec![Ternary::Pos, Ternary::Pos, Ternary::Pos];
        let sum = VectorModule::add(&a, &b).unwrap();
        assert_eq!(sum, vec![Ternary::Neg, Ternary::Zero, Ternary::Pos]);
        
        let scaled = VectorModule::scale(Ternary::Neg, &a);
        assert_eq!(scaled, vec![Ternary::Neg, Ternary::Pos, Ternary::Zero]);
    }
}
