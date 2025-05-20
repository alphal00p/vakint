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
            "( user_space::sigma(user_space::some_args)*vk::p(1,user_space::mink4(4,33))*vk::p(2,user_space::mink4(4,33))*vk::p(1,user_space::mink4(4,11))*vk::p(2,user_space::mink4(4,22))+vk::k(3,user_space::mink4(4,11))*vk::k(3,user_space::mink4(4,22)) + vk::k(3,user_space::mink4(4,77))*vk::p(1,user_space::mink4(4,77)) ) \
             * vk::topo( vk::prop(9,vk::edge(66,66),vk::k(3),user_space::MUVsq,1 ))\
            ")
        .unwrap()
        .as_view(),
        params_from_f64(&[
            ("user_space::MUVsq".into(), 1.0), ("vk::mursq".into(), 1.0),
            ("vk::g(user_space::mink4(4,22),user_space::mink4(4,11))".into(), 1.0),
            ("vk::p(1,user_space::mink4(4,11))".into(), 1.0),
            ("vk::p(2,user_space::mink4(4,22))".into(), 1.0),
            ("user_space::sigma(user_space::some_args)".into(), 1.0),
            ].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-1, ("0.0".into(), "-36.79980696990178527744327373291".into()),),
            (0,  ("0.0".into(), "99.70093613678094097571666372594".into()),),
            (1,  ("0.0".into(), "-183.0878169087836963770976657747".into()),),
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
            "( user_space::sigma(user_space::some_args)*vk::p(1,user_space::mink4(4,33))*vk::p(2,user_space::mink4(4,33))*vk::p(1,user_space::mink4(4,11))*vk::p(2,user_space::mink4(4,22))+vk::k(3,user_space::mink4(4,11))*vk::k(3,user_space::mink4(4,22)) + vk::k(3,user_space::mink4(4,77))*vk::p(1,user_space::mink4(4,77)) ) \
             * vk::topo( vk::prop(9,vk::edge(66,66),vk::k(3),user_space::MUVsq,1 ))\
            ")
        .unwrap()
        .as_view(),
        params_from_f64(&[
            ("user_space::MUVsq".into(), 1.0), ("vk::mursq".into(), 1.0),
            ("vk::g(user_space::mink4(4,22),user_space::mink4(4,11))".into(), 1.0),
            ("vk::p(1,user_space::mink4(4,11))".into(), 1.0),
            ("vk::p(2,user_space::mink4(4,22))".into(), 1.0),
            ("user_space::sigma(user_space::some_args)".into(), 1.0),
            ].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-1, ("0.0".into(), "-36.79980696990178527744327373291".into()),),
            (0,  ("0.0".into(), "99.70093613678094097571666372594".into()),),
            (1,  ("0.0".into(), "-183.0878169087836963770976657747".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0

    );
}
