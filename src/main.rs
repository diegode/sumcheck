use ark_poly::{
    polynomial::multivariate::SparsePolynomial,
    DenseMVPolynomial, Polynomial, multivariate::{SparseTerm, Term},
    univariate::SparsePolynomial as UnivariateSparsePolynomial,
};
use ark_std::{test_rng, UniformRand};
use ark_test_curves::{fp128::Fq, Zero, One};

fn main() {
    let g = build_sample_polynomial();
    let rng = &mut test_rng();
    let mut r: Vec<Fq> = vec![];
    let mut previous_gj = UnivariateSparsePolynomial::zero();
    for j in 0..g.num_vars {
        let gj = evaluate_hypercube_sum_fixing_point_prefix(&g, &r);
        assert_eq!(gj.degree(), degree(&g, j));
        let gj_sum = gj.evaluate(&Fq::zero()) + gj.evaluate(&Fq::one());
        if j == 0 {
            assert_eq!(evaluate_hypercube_sum(&g), gj_sum);
        } else {
            assert_eq!(previous_gj.evaluate(&r[j-1]), gj_sum);
        }
        r.push(Fq::rand(rng));
        previous_gj = gj;
    }
    assert_eq!(previous_gj.evaluate(&r[g.num_vars-1]), g.evaluate(&r)); 
    println!("All good");
}

// Create multivariate polynomial 2*x_0^3 + x_0*x_2 + x_1*x_2
fn build_sample_polynomial() -> SparsePolynomial<Fq, SparseTerm> {
    SparsePolynomial::from_coefficients_vec(
        3,
        vec![
            (Fq::from(2), SparseTerm::new(vec![(0, 3)])),
            (Fq::one(), SparseTerm::new(vec![(0, 1), (2, 1)])),
            (Fq::one(), SparseTerm::new(vec![(1, 1), (2, 1)])),
        ],
    )
}

fn evaluate_hypercube_sum(poly: &SparsePolynomial<Fq, SparseTerm>) -> Fq {
    let mut sum = Fq::zero();
    for i in 0..(1 << poly.num_vars) {
        let mut point = vec![];
        (0..poly.num_vars).for_each(|j: usize| {
            point.push(Fq::from((i >> j) & 1));
        });
        sum += poly.evaluate(&point);
    }
    sum
}

fn evaluate_hypercube_sum_fixing_point_prefix(poly: &SparsePolynomial<Fq, SparseTerm>, point_prefix: &[Fq]) -> UnivariateSparsePolynomial<Fq> {
    let mut sum = UnivariateSparsePolynomial::zero();
    let dimension = poly.num_vars - point_prefix.len() - 1;
    for i in 0..(1 << dimension) {
        let mut point = point_prefix.to_vec();
        point.push(Fq::one());
        (0..dimension).for_each(|j: usize| {
            point.push(Fq::from((i >> j) & 1));
        });
        sum = sum + evaluate_point_excluding_index(poly, &point, point_prefix.len());
    }
    sum
}

// Assumes that point[k] = 1
fn evaluate_point_excluding_index(poly: &SparsePolynomial<Fq, SparseTerm>, point: &[Fq], k: usize) -> UnivariateSparsePolynomial<Fq> {
    let mut result = UnivariateSparsePolynomial::zero();
    for (coeff, term) in poly.terms() {
        let result_coeff = coeff * &term.evaluate(point);
        let k_power = degree_in_term(term, k);
        let result_term = UnivariateSparsePolynomial::from_coefficients_vec(vec![(k_power, result_coeff)]);
        result = result + result_term;
    }
    result
}

fn degree_in_term(term: &SparseTerm, k: usize) -> usize {
    term.iter().filter(|(i, _)| *i == k).map(|(_, power)| power).sum()
}

fn degree(poly: &SparsePolynomial<Fq, SparseTerm>, k: usize) -> usize {
    poly.terms().iter().map(|(_, term)| degree_in_term(term, k)).max().unwrap_or(0)
}