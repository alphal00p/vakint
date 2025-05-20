mod test_utils;
use vakint::{
    utils::simplify_real, EvaluationMethod, EvaluationOrder, FMFTOptions, LoopNormalizationFactor,
    MATADOptions,
};

use std::vec;

use log::debug;
use std::collections::HashMap;
use symbolica::{
    atom::AtomCore,
    domains::{float::NumericalFloatLike, rational::Fraction},
    parse,
};
use test_utils::{compare_numerical_output, compare_output, get_vakint};
use vakint::{Vakint, VakintSettings};

use crate::test_utils::compare_vakint_evaluation_vs_reference;
use vakint::{externals_from_f64, params_from_f64, vakint_parse};

const N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS: u32 = 32;

#[test_log::test]
fn test_integrate_1l_decorated_indices_alphaloop() {
    #[rustfmt::skip]
    Vakint::initialize_vakint_symbols();
    compare_vakint_evaluation_vs_reference(
        VakintSettings{number_of_terms_in_epsilon_expansion: 3, integral_normalization_factor: LoopNormalizationFactor::MSbar,..VakintSettings::default()},
        EvaluationOrder::alphaloop_only(),
        parse!(
            "( vk::k(3,user_space::mink4(4,11))*vk::k(3,user_space::mink4(4,22)) + vk::k(3,user_space::mink4(4,77))*vk::p(8,user_space::mink4(4,77)) ) \
             * vk::topo( vk::prop(9,vk::edge(66,66),vk::k(3),user_space::MUVsq,1 ))\
            ")
        .unwrap()
        .as_view(),
        params_from_f64(&[
            ("user_space::MUVsq".into(), 1.0), ("vk::mursq".into(), 1.0),
            ("vk::g(user_space::mink4(4,22),user_space::mink4(4,11))".into(), 1.0)
            ].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        HashMap::default(),
        vec![
            (-1, ("0.0".into(), "2.467401100272339654708622749969".into()),),
            (0,  ("0.0".into(), "-5.368458141238928322097907419453".into()),),
            (1,  ("0.0".into(), "9.411704244710969435114178881646".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_1l_decorated_indices_matad() {
    #[rustfmt::skip]
    Vakint::initialize_vakint_symbols();
    compare_vakint_evaluation_vs_reference(
        VakintSettings{number_of_terms_in_epsilon_expansion: 3, integral_normalization_factor: LoopNormalizationFactor::pySecDec,..VakintSettings::default()},
        EvaluationOrder::matad_only(None),
        parse!(
            "( vk::k(3,user_space::mink4(4,11))*vk::k(3,user_space::mink4(4,22)) + vk::k(3,user_space::mink4(4,77))*vk::p(8,user_space::mink4(4,77)) ) \
             * vk::topo( vk::prop(9,vk::edge(66,66),vk::k(3),user_space::MUVsq,1 ))\
            ")
        .unwrap()
        .as_view(),
        params_from_f64(&[
            ("user_space::MUVsq".into(), 1.0), ("vk::mursq".into(), 1.0),
            ("vk::g(user_space::mink4(4,22),user_space::mink4(4,11))".into(), 1.0)
            ].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        HashMap::default(),
        vec![
            (-1, ("0.0".into(), "2.467401100272339654708622749969".into()),),
            (0,  ("0.0".into(), "-5.368458141238928322097907419453".into()),),
            (1,  ("0.0".into(), "9.411704244710969435114178881646".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}
