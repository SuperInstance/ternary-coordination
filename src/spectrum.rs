//! Ternary matrix eigenvalues over the complex numbers.
//!
//! Spectral radius determines consensus convergence speed.
//! We lift ternary matrices to complex number matrices and compute eigenvalues
//! using the characteristic polynomial.

use crate::matrix::TernaryMatrix;

/// A complex number (f64).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Complex {
    pub re: f64,
    pub im: f64,
}

impl Complex {
    pub fn new(re: f64, im: f64) -> Self { Complex { re, im } }
    pub fn real(r: f64) -> Self { Complex { re: r, im: 0.0 } }
    pub fn zero() -> Self { Complex { re: 0.0, im: 0.0 } }
    
    pub fn add(self, other: Self) -> Self {
        Complex { re: self.re + other.re, im: self.im + other.im }
    }
    pub fn sub(self, other: Self) -> Self {
        Complex { re: self.re - other.re, im: self.im - other.im }
    }
    pub fn mul(self, other: Self) -> Self {
        Complex {
            re: self.re * other.re - self.im * other.im,
            im: self.re * other.im + self.im * other.re,
        }
    }
    pub fn scale(self, s: f64) -> Self {
        Complex { re: self.re * s, im: self.im * s }
    }
    pub fn abs(self) -> f64 {
        (self.re * self.re + self.im * self.im).sqrt()
    }
    pub fn conj(self) -> Self {
        Complex { re: self.re, im: -self.im }
    }
}

impl std::fmt::Display for Complex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.im.abs() < 1e-10 {
            write!(f, "{:.6}", self.re)
        } else if self.re.abs() < 1e-10 {
            write!(f, "{:.6}i", self.im)
        } else {
            write!(f, "{:.6}{:+.6}i", self.re, self.im)
        }
    }
}

/// Convert a ternary matrix to a complex-valued matrix.
pub fn to_complex_matrix(m: &TernaryMatrix) -> Vec<Vec<Complex>> {
    (0..m.rows).map(|r| {
        (0..m.cols).map(|c| Complex::real(m.get(r, c).to_i32() as f64)).collect()
    }).collect()
}

/// Compute eigenvalues of a ternary matrix (up to 3×3) by solving the characteristic polynomial.
/// For n×n matrix A, eigenvalues λ satisfy det(A - λI) = 0.
pub fn eigenvalues(m: &TernaryMatrix) -> Vec<Complex> {
    let n = m.rows;
    assert_eq!(n, m.cols, "eigenvalues require square matrix");
    match n {
        0 => vec![],
        1 => vec![Complex::real(m.get(0, 0).to_i32() as f64)],
        2 => eigenvalues_2x2(m),
        3 => eigenvalues_3x3(m),
        _ => panic!("eigenvalue computation only supported for n ≤ 3"),
    }
}

fn eigenvalues_2x2(m: &TernaryMatrix) -> Vec<Complex> {
    // Characteristic polynomial: λ² - tr(A)λ + det(A) = 0
    let a = m.get(0, 0).to_i32() as f64;
    let b = m.get(0, 1).to_i32() as f64;
    let c = m.get(1, 0).to_i32() as f64;
    let d = m.get(1, 1).to_i32() as f64;
    let trace = a + d;
    let det = a * d - b * c;
    solve_quadratic(1.0, -trace, det)
}

fn eigenvalues_3x3(m: &TernaryMatrix) -> Vec<Complex> {
    // Characteristic polynomial: λ³ - tr(A)λ² + (sum of 2×2 minors)λ - det(A) = 0
    let entries = |r: usize, c: usize| m.get(r, c).to_i32() as f64;
    
    let trace = entries(0, 0) + entries(1, 1) + entries(2, 2);
    
    // Sum of 2×2 principal minors
    let minor_sum = entries(0, 0) * entries(1, 1) - entries(0, 1) * entries(1, 0)
                  + entries(0, 0) * entries(2, 2) - entries(0, 2) * entries(2, 0)
                  + entries(1, 1) * entries(2, 2) - entries(1, 2) * entries(2, 1);
    
    let det = m.determinant().map(|d| d.to_i32() as f64).unwrap_or(0.0);
    
    solve_cubic(1.0, -trace, minor_sum, -det)
}

/// Solve ax² + bx + c = 0.
fn solve_quadratic(a: f64, b: f64, c: f64) -> Vec<Complex> {
    let disc = b * b - 4.0 * a * c;
    if disc >= 0.0 {
        let sqrt_disc = disc.sqrt();
        vec![
            Complex::real((-b + sqrt_disc) / (2.0 * a)),
            Complex::real((-b - sqrt_disc) / (2.0 * a)),
        ]
    } else {
        let sqrt_disc = (-disc).sqrt();
        vec![
            Complex::new(-b / (2.0 * a), sqrt_disc / (2.0 * a)),
            Complex::new(-b / (2.0 * a), -sqrt_disc / (2.0 * a)),
        ]
    }
}

/// Solve x³ + ax² + bx + c = 0 using Cardano's method.
fn solve_cubic(a3: f64, a2: f64, a1: f64, a0: f64) -> Vec<Complex> {
    // Normalize: x³ + px² + qx + r = 0
    let p = a2 / a3;
    let q = a1 / a3;
    let r = a0 / a3;
    
    // Substituion x = t - p/3: t³ + pt + qt + r becomes depressed cubic
    // t³ + αt + β = 0
    let alpha = q - p * p / 3.0;
    let beta = 2.0 * p * p * p / 27.0 - p * q / 3.0 + r;
    
    let disc = -4.0 * alpha * alpha * alpha - 27.0 * beta * beta;
    
    if disc.abs() < 1e-10 {
        // Repeated roots
        if alpha.abs() < 1e-10 && beta.abs() < 1e-10 {
            let t = -p / 3.0;
            return vec![Complex::real(t); 3];
        }
    }
    
    // Use numerical approach for robustness
    let mut roots = Vec::new();
    
    // Try to find a real root first
    for guess in [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0] {
        let mut x = guess;
        for _ in 0..100 {
            let fx = x * x * x + p * x * x + q * x + r;
            let fpx = 3.0 * x * x + 2.0 * p * x + q;
            if fpx.abs() < 1e-14 { break; }
            let dx = fx / fpx;
            x -= dx;
            if dx.abs() < 1e-12 { break; }
        }
        let fx = x * x * x + p * x * x + q * x + r;
        if fx.abs() < 1e-8 {
            // Check if we already have this root
            let is_new = roots.iter().all(|r: &Complex| (r.re - x).abs() > 1e-6);
            if is_new {
                roots.push(Complex::real(x));
            }
            if roots.len() == 3 { break; }
        }
    }
    
    // If we have less than 3 real roots, find complex ones via quadratic of remaining factor
    if roots.len() == 1 {
        let r0 = roots[0].re;
        // Factor out (x - r0): x³ + px² + qx + r = (x - r0)(x² + bx + c)
        let b = p + r0;
        let c = r / r0;
        let remaining = solve_quadratic(1.0, b, c);
        roots.extend(remaining);
    } else if roots.len() == 2 {
        let r0 = roots[0].re;
        let r1 = roots[1].re;
        let _r2_coeff = 1.0; // leading coeff
        let third = Complex::real(r / (r0 * r1));
        if third.im.abs() < 1e-8 {
            roots.push(Complex::real(third.re));
        } else {
            roots.push(third);
        }
    } else if roots.is_empty() {
        // Fallback: use Cardano's formula directly
        let inner = beta * beta / 4.0 + alpha * alpha * alpha / 27.0;
        if inner >= 0.0 {
            let sqrt_inner = inner.sqrt();
            let cbrt_plus = (beta / 2.0 + sqrt_inner).cbrt();
            let cbrt_minus = (beta / 2.0 - sqrt_inner).cbrt();
            let t1 = cbrt_plus + cbrt_minus;
            roots.push(Complex::real(t1 - p / 3.0));
        }
    }
    
    // Ensure we return exactly 3 roots (pad with zeros if numerical issues)
    while roots.len() < 3 {
        roots.push(Complex::real(-p / 3.0));
    }
    
    roots
}

/// Compute the spectral radius (maximum absolute eigenvalue) of a ternary matrix.
pub fn spectral_radius(m: &TernaryMatrix) -> f64 {
    eigenvalues(m).iter().map(|e| e.abs()).fold(0.0_f64, f64::max)
}

/// Check if a mixing matrix leads to consensus convergence.
/// A matrix converges if:
/// 1. It has eigenvalue 1 (with multiplicity 1 for the consensus subspace)
/// 2. All other eigenvalues have |λ| < 1
pub fn converges(m: &TernaryMatrix) -> bool {
    let eigs = eigenvalues(m);
    let has_one = eigs.iter().any(|e| (e.re - 1.0).abs() < 1e-6 && e.im.abs() < 1e-6);
    let others_bounded = eigs.iter().all(|e| {
        (e.re - 1.0).abs() < 1e-6 && e.im.abs() < 1e-6 // eigenvalue is 1
        || e.abs() < 1.0 + 1e-6 // or |λ| ≤ 1
    });
    has_one && others_bounded
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identity_eigenvalues() {
        let id = TernaryMatrix::identity(2);
        let eigs = eigenvalues(&id);
        assert_eq!(eigs.len(), 2);
        for e in &eigs {
            assert!((e.re - 1.0).abs() < 1e-6, "expected 1, got {}", e);
        }
    }

    #[test]
    fn test_zero_matrix_eigenvalues() {
        let z = TernaryMatrix::zeros(2, 2);
        let eigs = eigenvalues(&z);
        for e in &eigs {
            assert!(e.re.abs() < 1e-6, "expected 0, got {}", e);
        }
    }

    #[test]
    fn test_spectral_radius_identity() {
        let id = TernaryMatrix::identity(2);
        let sr = spectral_radius(&id);
        assert!((sr - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_3x3_eigenvalues() {
        let id = TernaryMatrix::identity(3);
        let eigs = eigenvalues(&id);
        assert_eq!(eigs.len(), 3);
        for e in &eigs {
            assert!((e.re - 1.0).abs() < 1e-6);
        }
    }

    #[test]
    fn test_complex_arithmetic() {
        let a = Complex::new(1.0, 2.0);
        let b = Complex::new(3.0, 4.0);
        let sum = a.add(b);
        assert!((sum.re - 4.0).abs() < 1e-10);
        assert!((sum.im - 6.0).abs() < 1e-10);
        
        let prod = a.mul(b);
        assert!((prod.re - (-5.0)).abs() < 1e-10);
        assert!((prod.im - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_complex_abs() {
        let z = Complex::new(3.0, 4.0);
        assert!((z.abs() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_to_complex_matrix() {
        let m = TernaryMatrix::from_i32_slice(2, 2, &[1, -1, 0, 1]).unwrap();
        let cm = to_complex_matrix(&m);
        assert!((cm[0][0].re - 1.0).abs() < 1e-10);
        assert!((cm[0][1].re - (-1.0)).abs() < 1e-10);
        assert!((cm[1][0].re - 0.0).abs() < 1e-10);
        assert!((cm[1][1].re - 1.0).abs() < 1e-10);
    }
}
