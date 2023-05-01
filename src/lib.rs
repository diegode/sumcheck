use ark_poly::{
    polynomial::multivariate::SparsePolynomial,
    DenseMVPolynomial, Polynomial, multivariate::{SparseTerm, Term},
    univariate::SparsePolynomial as UnivariateSparsePolynomial,
};
use ark_std::{test_rng, UniformRand};
use ark_test_curves::{fp128::Fq, Zero, One};

struct Prover {
    pub g: SparsePolynomial<Fq, SparseTerm>
}

impl Prover {
    fn evaluate_hypercube_sum(&self) -> Fq {
        let mut sum = Fq::zero();
        for i in 0..(1 << self.g.num_vars) {
            let mut point = vec![];
            (0..self.g.num_vars).for_each(|j: usize| {
                point.push(Fq::from((i >> j) & 1));
            });
            sum += self.g.evaluate(&point);
        }
        sum
    }
    
    fn evaluate_hypercube_sum_fixing_point_prefix(&self, point_prefix: &[Fq]) -> UnivariateSparsePolynomial<Fq> {
        let mut sum = UnivariateSparsePolynomial::zero();
        let dimension = self.g.num_vars - point_prefix.len() - 1;
        for i in 0..(1 << dimension) {
            let mut point = point_prefix.to_vec();
            point.push(Fq::one());
            (0..dimension).for_each(|j: usize| {
                point.push(Fq::from((i >> j) & 1));
            });
            sum = sum + self.evaluate_point_excluding_index(&point, point_prefix.len());
        }
        sum
    }
    
    // Assumes that point[k] = 1
    fn evaluate_point_excluding_index(&self, point: &[Fq], k: usize) -> UnivariateSparsePolynomial<Fq> {
        let mut result = UnivariateSparsePolynomial::zero();
        for (coeff, term) in self.g.terms() {
            let result_coeff = coeff * &term.evaluate(point);
            let k_power = degree_in_term(term, k);
            let result_term = UnivariateSparsePolynomial::from_coefficients_vec(vec![(k_power, result_coeff)]);
            result = result + result_term;
        }
        result
    }    
}

pub fn sumcheck_protocol(g: &SparsePolynomial<Fq, SparseTerm>) -> Option<Fq> {
    let prover = Prover { g: g.clone() };
    let hypercube_sum = prover.evaluate_hypercube_sum();
    let rng = &mut test_rng();
    let mut r: Vec<Fq> = vec![];
    let mut previous_gj_eval = Fq::zero();
    for j in 0..g.num_vars {
        let gj = prover.evaluate_hypercube_sum_fixing_point_prefix(&r);
        if gj.degree() != degree(&g, j) {
            return None;
        }
        let expected_gj_sum = if j == 0 { hypercube_sum } else { previous_gj_eval };
        if gj.evaluate(&Fq::zero()) + gj.evaluate(&Fq::one()) != expected_gj_sum {
            return None;
        }
        let rj = Fq::rand(rng);
        r.push(rj);
        previous_gj_eval = gj.evaluate(&rj);
    }
    if previous_gj_eval != g.evaluate(&r) {
        return None;
    }
    Some(hypercube_sum)
}

fn degree_in_term(term: &SparseTerm, k: usize) -> usize {
    term.iter().filter(|(i, _)| *i == k).map(|(_, power)| power).sum()
}

fn degree(poly: &SparsePolynomial<Fq, SparseTerm>, k: usize) -> usize {
    poly.terms().iter().map(|(_, term)| degree_in_term(term, k)).max().unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;
 
    #[test]
    fn thaler_example() {
        // It is 2*x_0^3 + x_0*x_2 + x_1*x_2
        let g = SparsePolynomial::from_coefficients_vec(3,
            vec![
                (Fq::from(2), SparseTerm::new(vec![(0, 3)])),
                (Fq::one(), SparseTerm::new(vec![(0, 1), (2, 1)])),
                (Fq::one(), SparseTerm::new(vec![(1, 1), (2, 1)]))]);    
        assert_eq!(sumcheck_protocol(&g), Some(Fq::from(12)));
    }
}