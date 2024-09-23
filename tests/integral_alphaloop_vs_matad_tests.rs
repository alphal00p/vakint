mod test_utils;

use symbolica::atom::Atom;
use test_utils::compare_two_evaluations;
use vakint::{EvaluationOrder, VakintSettings};

use crate::test_utils::{convert_test_externals, convert_test_params};

const N_DIGITS_ANLYTICAL_EVALUATION: u32 = 32;
const COMPARISON_REL_THRESHOLD: f64 = 1.0e-25;
// No error on anylic expressions
const MAX_PULL: f64 = 0.0e0;

#[test_log::test]
fn test_integrate_3l_no_numerator() {
    #[rustfmt::skip]
    compare_two_evaluations(
        VakintSettings { n_digits_at_evaluation_time: N_DIGITS_ANLYTICAL_EVALUATION, number_of_terms_in_epsilon_expansion: 3, ..VakintSettings::default()},
        ((&EvaluationOrder::alphaloop_only() ,true),(&EvaluationOrder::matad_only() ,true)),
        Atom::parse(
            "( k(1,11)*k(1,11) )
            *topo(\
                 prop(1,edge(1,2),k(1),muvsq,1)\
                *prop(2,edge(2,3),k(2),muvsq,1)\
                *prop(3,edge(3,1),k(3),muvsq,1)\
                *prop(4,edge(1,4),k(3)-k(1),muvsq,1)\
                *prop(5,edge(2,4),k(1)-k(2),muvsq,1)\
                *prop(6,edge(3,4),k(2)-k(3),muvsq,1)\
            )",
        ).unwrap().as_view(),
        convert_test_params(&[("muvsq".into(), 1.0), ("mursq".into(), 2.0)].iter().cloned().collect(),
        N_DIGITS_ANLYTICAL_EVALUATION),
        convert_test_externals(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            //.map(|i| (i, (17.0*((0) as f64), 4.0*((0) as f64), 3.0*((0) as f64), 12.0*((0) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION),
            COMPARISON_REL_THRESHOLD, MAX_PULL,
        true,
    );
}

#[test_log::test]
fn test_integrate_3l_rank_4() {
    #[rustfmt::skip]
    compare_two_evaluations(
        VakintSettings { n_digits_at_evaluation_time: N_DIGITS_ANLYTICAL_EVALUATION, number_of_terms_in_epsilon_expansion: 3, ..VakintSettings::default()},
        ((&EvaluationOrder::alphaloop_only() ,true),(&EvaluationOrder::matad_only() ,true)),
        Atom::parse(
            "(
                  k(1,11)*k(2,11)*k(1,22)*k(2,22)
                + p(1,11)*k(3,11)*k(3,22)*p(2,22)
                + p(1,11)*p(2,11)*(k(2,22)+k(1,22))*k(2,22) 
             )
            *topo(\
                 prop(1,edge(1,2),k(1),muvsq,1)\
                *prop(2,edge(2,3),k(2),muvsq,1)\
                *prop(3,edge(3,1),k(3),muvsq,1)\
                *prop(4,edge(1,4),k(3)-k(1),muvsq,1)\
                *prop(5,edge(2,4),k(1)-k(2),muvsq,1)\
                *prop(6,edge(3,4),k(2)-k(3),muvsq,1)\
            )",
        ).unwrap().as_view(),
        convert_test_params(&[("muvsq".into(), 1.0), ("mursq".into(), 2.0)].iter().cloned().collect(),
        N_DIGITS_ANLYTICAL_EVALUATION),
        convert_test_externals(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            //.map(|i| (i, (17.0*((0) as f64), 4.0*((0) as f64), 3.0*((0) as f64), 12.0*((0) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION),
            COMPARISON_REL_THRESHOLD, MAX_PULL,
        true,
    );
}
