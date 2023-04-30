use ark_poly::{
    polynomial::multivariate::SparsePolynomial,
    DenseMVPolynomial, Polynomial, multivariate::{SparseTerm, Term},
    univariate::SparsePolynomial as UnivariateSparsePolynomial,
};
use ark_test_curves::{fp128::Fq, Zero};

fn main() {
    let g = build_sample_polynomial();
    let c1 = evaluate_hypercube_sum(&g);

    let g1 = evaluate_hypercube_sum_excluding_index(&g, 0);
    if c1 == g1.evaluate(&Fq::zero()) + g1.evaluate(&Fq::from(1)) {
        println!("all good");
    }
}

// Create multivariate polynomial 2*x_0^3 + x_0*x_2 + x_1*x_2
fn build_sample_polynomial() -> SparsePolynomial<Fq, SparseTerm> {
    SparsePolynomial::from_coefficients_vec(
        3,
        vec![
            (Fq::from(2), SparseTerm::new(vec![(0, 3)])),
            (Fq::from(1), SparseTerm::new(vec![(0, 1), (2, 1)])),
            (Fq::from(1), SparseTerm::new(vec![(1, 1), (2, 1)])),
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

fn evaluate_hypercube_sum_excluding_index(poly: &SparsePolynomial<Fq, SparseTerm>, k: usize) -> UnivariateSparsePolynomial<Fq> {
    let mut sum = UnivariateSparsePolynomial::zero();
    for i in 0..(1 << (poly.num_vars - 1)) {
        let mut point = vec![];
        (0..(poly.num_vars - 1)).for_each(|j: usize| {
            point.push(Fq::from((i >> j) & 1));
        });
        point.insert(k, Fq::from(1));
        sum = sum + evaluate_point_excluding_index(poly, &point, k);
    }
    sum
}

// Assumes that point[k] = 1
fn evaluate_point_excluding_index(poly: &SparsePolynomial<Fq, SparseTerm>, point: &[Fq], k: usize) -> UnivariateSparsePolynomial<Fq> {
    let mut result = UnivariateSparsePolynomial::zero();
    for term in poly.terms() {
        let coeff = term.0 * term.1.evaluate(point);
        let k_power = term.1.iter().filter(|(i, _)| *i == k).map(|(_, power)| power).sum();
        let result_term = UnivariateSparsePolynomial::from_coefficients_vec(vec![(k_power, coeff)]);
        result = result + result_term;
    }
    result
}
