use ark_poly::{
    DenseMVPolynomial, Polynomial, multivariate::{SparseTerm, Term},
};
use ark_std::UniformRand;
use ark_test_curves::{fp128::Fq as Field, Zero, One};
use rand::{SeedableRng, rngs::StdRng};

pub type MultivariatePolynomial = ark_poly::polynomial::multivariate::SparsePolynomial<Field, SparseTerm>;
pub type UnivariatePolynomial = ark_poly::univariate::SparsePolynomial<Field>;

struct Prover {
    g: MultivariatePolynomial
}

impl Prover {
    fn evaluate_hypercube_sum(&self) -> Field {
        let mut sum = Field::zero();
        for i in 0..(1 << self.g.num_vars) {
            let mut point = vec![];
            (0..self.g.num_vars).for_each(|j| {
                point.push(Field::from((i >> j) & 1));
            });
            sum += self.g.evaluate(&point);
        }
        sum
    }
    
    fn evaluate_hypercube_sum_fixing_point_prefix(&self, point_prefix: &[Field]) -> UnivariatePolynomial {
        let mut sum = UnivariatePolynomial::zero();
        let dimension = self.g.num_vars - point_prefix.len() - 1;
        for i in 0..(1 << dimension) {
            let mut point = point_prefix.to_vec();
            point.push(Field::one());
            (0..dimension).for_each(|j| {
                point.push(Field::from((i >> j) & 1));
            });
            sum = sum + self.evaluate_point_excluding_index(&point, point_prefix.len());
        }
        sum
    }
    
    // Assumes that point[k] = 1
    fn evaluate_point_excluding_index(&self, point: &[Field], k: usize) -> UnivariatePolynomial {
        let mut result = UnivariatePolynomial::zero();
        for (coeff, term) in self.g.terms() {
            let result_coeff = coeff * &term.evaluate(point);
            let k_power = degree_in_term(term, k);
            let result_term = UnivariatePolynomial::from_coefficients_vec(vec![(k_power, result_coeff)]);
            result = result + result_term;
        }
        result
    }    
}

pub fn sumcheck_protocol(g: &MultivariatePolynomial) -> Option<Field> {
    let prover = Prover { g: g.clone() };
    let hypercube_sum = prover.evaluate_hypercube_sum();
    let mut rng = StdRng::from_entropy();
    let mut r: Vec<Field> = vec![];
    let mut previous_gj_eval = hypercube_sum;
    for j in 0..g.num_vars {
        let gj = prover.evaluate_hypercube_sum_fixing_point_prefix(&r);
        if gj.degree() > degree(g, j) {
            return None;
        }
        if gj.evaluate(&Field::zero()) + gj.evaluate(&Field::one()) != previous_gj_eval {
            return None;
        }
        let rj = Field::rand(&mut rng);
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

fn degree(poly: &MultivariatePolynomial, k: usize) -> usize {
    poly.terms().iter().map(|(_, term)| degree_in_term(term, k)).max().unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;
 
    #[test]
    fn thaler_example() {
        // It is 2*x_0^3 + x_0*x_2 + x_1*x_2
        let g = MultivariatePolynomial::from_coefficients_vec(3,
            vec![
                (Field::from(2), SparseTerm::new(vec![(0, 3)])),
                (Field::one(), SparseTerm::new(vec![(0, 1), (2, 1)])),
                (Field::one(), SparseTerm::new(vec![(1, 1), (2, 1)]))]);
        assert_eq!(sumcheck_protocol(&g), Some(Field::from(12)));
    }
}