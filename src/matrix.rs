//! Ternary matrices over Z/3Z.
//!
//! Supports multiplication, transpose, determinant, inverse (when it exists),
//! and row reduction over the ternary field.

use serde::{Deserialize, Serialize};
use crate::ternary::Ternary;

/// A matrix with ternary entries.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct TernaryMatrix {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<Ternary>,
}

impl TernaryMatrix {
    /// Create a zero matrix of given dimensions.
    pub fn zeros(rows: usize, cols: usize) -> Self {
        TernaryMatrix {
            rows,
            cols,
            data: vec![Ternary::Zero; rows * cols],
        }
    }

    /// Create an identity matrix of given size.
    pub fn identity(n: usize) -> Self {
        let mut m = Self::zeros(n, n);
        for i in 0..n {
            m.data[i * n + i] = Ternary::Pos;
        }
        m
    }

    /// Create from a Vec of Ternary values in row-major order.
    pub fn from_vec(rows: usize, cols: usize, data: Vec<Ternary>) -> Option<Self> {
        if data.len() != rows * cols {
            return None;
        }
        Some(TernaryMatrix { rows, cols, data })
    }

    /// Create from i32 values (converted mod 3).
    pub fn from_i32_slice(rows: usize, cols: usize, vals: &[i32]) -> Option<Self> {
        if vals.len() != rows * cols {
            return None;
        }
        Some(TernaryMatrix {
            rows,
            cols,
            data: vals.iter().map(|&v| Ternary::from_i32(v)).collect(),
        })
    }

    /// Get element at (r, c).
    pub fn get(&self, r: usize, c: usize) -> Ternary {
        self.data[r * self.cols + c]
    }

    /// Set element at (r, c).
    pub fn set(&mut self, r: usize, c: usize, v: Ternary) {
        self.data[r * self.cols + c] = v;
    }

    /// Transpose the matrix.
    pub fn transpose(&self) -> Self {
        let mut result = Self::zeros(self.cols, self.rows);
        for r in 0..self.rows {
            for c in 0..self.cols {
                result.set(c, r, self.get(r, c));
            }
        }
        result
    }

    /// Matrix multiplication over Z/3Z.
    pub fn mul(&self, other: &TernaryMatrix) -> Option<TernaryMatrix> {
        if self.cols != other.rows {
            return None;
        }
        let mut result = Self::zeros(self.rows, other.cols);
        for r in 0..self.rows {
            for c in 0..other.cols {
                let mut sum = Ternary::Zero;
                for k in 0..self.cols {
                    sum = sum.add(self.get(r, k).mul(other.get(k, c)));
                }
                result.set(r, c, sum);
            }
        }
        Some(result)
    }

    /// Scalar multiplication.
    pub fn scale(&self, s: Ternary) -> Self {
        TernaryMatrix {
            rows: self.rows,
            cols: self.cols,
            data: self.data.iter().map(|&v| v.mul(s)).collect(),
        }
    }

    /// Add two matrices.
    pub fn add(&self, other: &TernaryMatrix) -> Option<TernaryMatrix> {
        if self.rows != other.rows || self.cols != other.cols {
            return None;
        }
        let data: Vec<Ternary> = self.data.iter()
            .zip(other.data.iter())
            .map(|(&a, &b)| a.add(b))
            .collect();
        Some(TernaryMatrix { rows: self.rows, cols: self.cols, data })
    }

    /// Determinant for square matrices. Uses cofactor expansion.
    /// Returns None for non-square or matrices larger than we handle.
    pub fn determinant(&self) -> Option<Ternary> {
        if self.rows != self.cols {
            return None;
        }
        Some(self.det_recursive())
    }

    fn det_recursive(&self) -> Ternary {
        let n = self.rows;
        match n {
            0 => Ternary::Pos, // empty determinant = 1 by convention
            1 => self.get(0, 0),
            2 => {
                // ad - bc
                let a = self.get(0, 0);
                let b = self.get(0, 1);
                let c = self.get(1, 0);
                let d = self.get(1, 1);
                a.mul(d).sub(b.mul(c))
            }
            _ => {
                // Cofactor expansion along first row
                let mut det = Ternary::Zero;
                for c in 0..n {
                    let cofactor = self.minor(0, c);
                    let sign = if c % 2 == 0 { Ternary::Pos } else { Ternary::Neg };
                    det = det.add(sign.mul(self.get(0, c)).mul(cofactor.det_recursive()));
                }
                det
            }
        }
    }

    /// Minor matrix: delete row r and column c.
    fn minor(&self, r: usize, c: usize) -> TernaryMatrix {
        let n = self.rows - 1;
        let mut data = Vec::with_capacity(n * n);
        for ri in 0..self.rows {
            if ri == r { continue; }
            for ci in 0..self.cols {
                if ci == c { continue; }
                data.push(self.get(ri, ci));
            }
        }
        TernaryMatrix { rows: n, cols: n, data }
    }

    /// Compute the inverse over Z/3Z, if it exists (determinant must be nonzero).
    pub fn inverse(&self) -> Option<TernaryMatrix> {
        if self.rows != self.cols { return None; }
        let det = self.determinant()?;
        let det_inv = det.inv()?;
        Some(self.adjugate().scale(det_inv))
    }

    /// Adjugate (classical adjoint) matrix.
    fn adjugate(&self) -> TernaryMatrix {
        let n = self.rows;
        let mut data = Vec::with_capacity(n * n);
        for c in 0..n {
            for r in 0..n {
                let sign = if (r + c) % 2 == 0 { Ternary::Pos } else { Ternary::Neg };
                let cofactor = sign.mul(self.minor(r, c).det_recursive());
                data.push(cofactor);
            }
        }
        TernaryMatrix { rows: n, cols: n, data }
    }

    /// Row reduction to row echelon form over Z/3Z.
    /// Returns the reduced matrix and the rank.
    pub fn row_echelon(&self) -> (TernaryMatrix, usize) {
        let mut m = self.clone();
        let mut pivot_row = 0;
        for col in 0..m.cols {
            if pivot_row >= m.rows { break; }
            // Find pivot
            let mut found = None;
            for r in pivot_row..m.rows {
                if !m.get(r, col).is_zero() {
                    found = Some(r);
                    break;
                }
            }
            let Some(pr) = found else { continue; };
            // Swap rows
            if pr != pivot_row {
                for c in 0..m.cols {
                    let tmp = m.get(pivot_row, c);
                    m.set(pivot_row, c, m.get(pr, c));
                    m.set(pr, c, tmp);
                }
            }
            // Scale pivot row so pivot = 1
            let pivot_val = m.get(pivot_row, col);
            let scale = pivot_val.inv().unwrap(); // nonzero, so has inverse
            for c in 0..m.cols {
                m.set(pivot_row, c, m.get(pivot_row, c).mul(scale));
            }
            // Eliminate below
            for r in 0..m.rows {
                if r == pivot_row { continue; }
                let factor = m.get(r, col);
                if factor.is_zero() { continue; }
                for c in 0..m.cols {
                    let val = m.get(r, c).sub(factor.mul(m.get(pivot_row, c)));
                    m.set(r, c, val);
                }
            }
            pivot_row += 1;
        }
        (m, pivot_row)
    }

    /// Check if the matrix is square.
    pub fn is_square(&self) -> bool {
        self.rows == self.cols
    }
}

impl std::fmt::Display for TernaryMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for r in 0..self.rows {
            let row: Vec<String> = (0..self.cols)
                .map(|c| format!("{:>3}", self.get(r, c)))
                .collect();
            writeln!(f, "[{}]", row.join(" "))?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identity_mul() {
        let id = TernaryMatrix::identity(3);
        let a = TernaryMatrix::from_i32_slice(3, 3, &[
            1, -1, 0,
            0, 1, 1,
            -1, 0, 1,
        ]).unwrap();
        assert_eq!(id.mul(&a).unwrap(), a);
        assert_eq!(a.mul(&id).unwrap(), a);
    }

    #[test]
    fn test_matrix_multiply_2x2() {
        let a = TernaryMatrix::from_i32_slice(2, 2, &[1, 1, 0, 1]).unwrap();
        let b = TernaryMatrix::from_i32_slice(2, 2, &[1, 0, -1, 1]).unwrap();
        // [1 1] * [1  0] = [1+(-1)  0+1] = [0 1]
        // [0 1]   [-1 1]   [0+(-1)  0+1]   [-1 1]
        let c = a.mul(&b).unwrap();
        assert_eq!(c.get(0, 0), Ternary::Zero);
        assert_eq!(c.get(0, 1), Ternary::Pos);
        assert_eq!(c.get(1, 0), Ternary::Neg);
        assert_eq!(c.get(1, 1), Ternary::Pos);
    }

    #[test]
    fn test_transpose() {
        let a = TernaryMatrix::from_i32_slice(2, 3, &[
            1, -1, 0,
            0, 1, 1,
        ]).unwrap();
        let t = a.transpose();
        assert_eq!(t.rows, 3);
        assert_eq!(t.cols, 2);
        assert_eq!(t.get(0, 0), Ternary::Pos);
        assert_eq!(t.get(1, 0), Ternary::Neg);
        assert_eq!(t.get(2, 0), Ternary::Zero);
        assert_eq!(t.get(0, 1), Ternary::Zero);
        assert_eq!(t.get(1, 1), Ternary::Pos);
        assert_eq!(t.get(2, 1), Ternary::Pos);
    }

    #[test]
    fn test_double_transpose() {
        let a = TernaryMatrix::from_i32_slice(2, 3, &[
            1, -1, 0,
            0, 1, 1,
        ]).unwrap();
        assert_eq!(a.transpose().transpose(), a);
    }

    #[test]
    fn test_determinant_2x2() {
        // det [1 1] = 1*1 - 1*(-1) = 1+1 = -1 (mod 3)
        //     [-1 1]
        let m = TernaryMatrix::from_i32_slice(2, 2, &[1, 1, -1, 1]).unwrap();
        assert_eq!(m.determinant(), Some(Ternary::Neg)); // 1 - (-1) = 2 ≡ -1

        // det [1 0] = 1
        //     [0 1]
        let id2 = TernaryMatrix::identity(2);
        assert_eq!(id2.determinant(), Some(Ternary::Pos));

        // Singular: det [1 1] = 1*(-1) - 1*1 = -1-1 = 1 (mod 3)... wait
        // [1 1]: det = 1*(-1) - 1*1 = -1-1 = -2 = 1 (mod 3)
        let sing = TernaryMatrix::from_i32_slice(2, 2, &[1, 1, 1, -1]).unwrap();
        assert_eq!(sing.determinant(), Some(Ternary::Pos)); // 1*(-1) - 1*1 = -2 ≡ 1
    }

    #[test]
    fn test_determinant_3x3() {
        let id3 = TernaryMatrix::identity(3);
        assert_eq!(id3.determinant(), Some(Ternary::Pos));

        // [1 1 1]
        // [1 0 1]
        // [0 1 0]
        // det = 1*(0*0 - 1*1) - 1*(1*0 - 1*0) + 1*(1*1 - 0*0)
        //     = 1*(-1) - 1*(0) + 1*(1) = -1 + 0 + 1 = 0
        let m = TernaryMatrix::from_i32_slice(3, 3, &[
            1, 1, 1,
            1, 0, 1,
            0, 1, 0,
        ]).unwrap();
        assert_eq!(m.determinant(), Some(Ternary::Zero));
    }

    #[test]
    fn test_inverse_2x2() {
        // [1 1] has det = -1, inverse of -1 is -1
        // [-1 1]
        let m = TernaryMatrix::from_i32_slice(2, 2, &[1, 1, -1, 1]).unwrap();
        let inv = m.inverse().expect("should be invertible");
        let product = m.mul(&inv).unwrap();
        assert_eq!(product, TernaryMatrix::identity(2));
    }

    #[test]
    fn test_inverse_singular() {
        // All-zeros matrix is singular
        let m = TernaryMatrix::zeros(2, 2);
        assert!(m.inverse().is_none());
    }

    #[test]
    fn test_row_echelon_identity() {
        let id = TernaryMatrix::identity(3);
        let (reduced, rank) = id.row_echelon();
        assert_eq!(rank, 3);
        assert_eq!(reduced, id);
    }

    #[test]
    fn test_row_echelon_reduces() {
        // [1 1] → R2 = R2 - R1 = [0 0]. Rank 1.
        // [1 1]
        // In Z/3Z: R2 - R1 means subtracting, which is same as adding neg.
        let m = TernaryMatrix::from_i32_slice(2, 2, &[1, 1, 1, 1]).unwrap();
        let (reduced, rank) = m.row_echelon();
        assert_eq!(rank, 1);
        // Row 0 should be [1 1] (pivot at col 0, but col 1 = 1, not zero because
        // row echelon doesn't guarantee reduced form for non-trivial row)
        // Actually, our algorithm eliminates BOTH above and below.
        // R1 is [1 1]. Pivot = 1 at (0,0). Scale row by inv(1)=1. 
        // R2 = R2 - 1*R1 = [1-1, 1-1] = [0, 0]. 
        assert_eq!(reduced.get(0, 0), Ternary::Pos);
        // After normalization: R1 has pivot at col 0, no more pivots.
        assert_eq!(reduced.get(1, 0), Ternary::Zero);
        assert_eq!(reduced.get(1, 1), Ternary::Zero);
    }

    #[test]
    fn test_row_echelon_2x3() {
        // [1 0 1] → R2 = R2 - R1
        // [1 0 -1]  → [0 0 1-(-1)] = [0 0 -1] since 1-(-1)=2≡-1
        // But our algorithm also eliminates above: R1 = R1 - 1*R2
        // [1 0 1] → [1 0 1-1*(-1)] = [1 0 1+(-1)] ... wait, the algorithm
        // only eliminates below, not above (row echelon, not RREF).
        // Actually looking at the code, it eliminates both above AND below.
        // So this is actually RREF.
        let m = TernaryMatrix::from_i32_slice(2, 3, &[1, 0, 1, 1, 0, -1]).unwrap();
        let (reduced, rank) = m.row_echelon();
        assert_eq!(rank, 2);
        assert_eq!(reduced.get(0, 0), Ternary::Pos);
        assert_eq!(reduced.get(1, 2), Ternary::Pos);
    }

    #[test]
    fn test_matrix_add() {
        let a = TernaryMatrix::from_i32_slice(2, 2, &[1, 0, -1, 1]).unwrap();
        let b = TernaryMatrix::from_i32_slice(2, 2, &[1, 1, 1, 0]).unwrap();
        let c = a.add(&b).unwrap();
        assert_eq!(c.get(0, 0), Ternary::Neg); // 1+1 = -1
        assert_eq!(c.get(0, 1), Ternary::Pos);
        assert_eq!(c.get(1, 0), Ternary::Zero);
        assert_eq!(c.get(1, 1), Ternary::Pos);
    }

    #[test]
    fn test_scale() {
        let m = TernaryMatrix::from_i32_slice(2, 2, &[1, -1, 0, 1]).unwrap();
        let s = m.scale(Ternary::Neg);
        assert_eq!(s.get(0, 0), Ternary::Neg);
        assert_eq!(s.get(0, 1), Ternary::Pos);
        assert_eq!(s.get(1, 0), Ternary::Zero);
        assert_eq!(s.get(1, 1), Ternary::Neg);
    }
}
