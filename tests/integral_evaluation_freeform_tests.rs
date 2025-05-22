mod test_utils;
use vakint::{EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, PySecDecOptions};

use std::vec;

use symbolica::parse;
use vakint::{Vakint, VakintSettings};

use crate::test_utils::compare_vakint_evaluation_vs_reference;
use vakint::{externals_from_f64, params_from_f64};

const N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS: u32 = 32;

const N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS: u32 = 10;
// PySecDec QMC is often very optimistic
const MAX_PULL: f64 = 1.0e99;

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
            (0,  ("0.0".into(), "-35.56610641976561545008896235793".into()),),
            (1,  ("0.0".into(), "-65.21588421381265673942992661242".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_1l_decorated_indices_matad() {
    #[rustfmt::skip]
    Vakint::initialize_vakint_symbols();
    compare_vakint_evaluation_vs_reference(
        VakintSettings{number_of_terms_in_epsilon_expansion: 3, integral_normalization_factor: LoopNormalizationFactor::MSbar,..VakintSettings::default()},
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
            (0,  ("0.0".into(), "-35.56610641976561545008896235793".into()),),
            (1,  ("0.0".into(), "-65.21588421381265673942992661242".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0

    );
}

#[test_log::test]
fn test_integrate_1l_decorated_indices_fmft() {
    #[rustfmt::skip]
    Vakint::initialize_vakint_symbols();
    compare_vakint_evaluation_vs_reference(
        VakintSettings{number_of_terms_in_epsilon_expansion: 5, integral_normalization_factor: LoopNormalizationFactor::MSbar,..VakintSettings::default()},
        EvaluationOrder::fmft_only(None),
        parse!(
            "( user_space::sigma(user_space::some_args)*vk::p(1,user_space::mink4(4,33))*vk::p(2,user_space::mink4(4,33))*vk::p(1,user_space::mink4(4,11))*vk::p(2,user_space::mink4(4,22))+vk::k(3,user_space::mink4(4,11))*vk::k(3,user_space::mink4(4,22)) + vk::k(3,user_space::mink4(4,77))*vk::p(1,user_space::mink4(4,77)) ) \
             * vk::topo( 
                  vk::prop(9,vk::edge(66,66),vk::k(1),user_space::MUVsq,1 )
                * vk::prop(9,vk::edge(66,66),vk::k(2),user_space::MUVsq,1 )
                * vk::prop(9,vk::edge(66,66),vk::k(3),user_space::MUVsq,1 )
                * vk::prop(9,vk::edge(66,66),vk::k(4),user_space::MUVsq,1 )
            )\
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
            (-4, ("-35378.93674652074486878056225484".into(), "0.0".into()),),
            (-3,  ("-140329.6806090741577242311770814".into(), "0.0".into()),),
            (-2,  ("-464844.1053751088101148680385287".into(), "0.0".into()),),
            (-1,  ("-1098012.239402848410636935684308".into(), "0.0".into()),),
            (0,  ("-2358474.482147627644123865090031".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0

    );
}

#[test_log::test]
fn test_integrate_1l_decorated_indices_pysecdec() {
    Vakint::initialize_vakint_symbols();
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{number_of_terms_in_epsilon_expansion: 5, integral_normalization_factor: LoopNormalizationFactor::MSbar, ..VakintSettings::default()},
        EvaluationOrder(vec![EvaluationMethod::PySecDec(PySecDecOptions { reuse_existing_output: Some("./tests_workspace/test_integrate_1l_decorated_indices_pysecdec".into()) ,..PySecDecOptions::default() })]),
        parse!(
            "(user_space::sigma*vk::k(1,user_space::mink4(4,11))*vk::p(1,user_space::mink4(4,11))*vk::k(1,user_space::mink4(4,12))*vk::p(1,user_space::mink4(4,12)))*vk::topo(\
                vk::prop(1,vk::edge(1,1),vk::k(1),user_space::muvsq,2)\
             )"
        )
        .unwrap()
        .as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[
            ("user_space::muvsq".into(), 1.0), 
            ("vk::mursq".into(), 1.0),
            ("user_space::sigma".into(), 1.0)
            ].iter().cloned().collect(),
            N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=1)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS),
        vec![
            (-1, ("0.0".into(),  "-15.4182959953817973987664163098897650903".into()),),
            (0,  ("0.0".into(),  "-15.4182959953817973987664163098897650903".into()),),
            (1,  ("0.0".into(),  "-28.0993361631583491418567542617704801124".into()),),
            (2,  ("0.0".into(),  "-21.9214464510894766034570152742441344730".into()),),
            (3,  ("0.0".into(),  "-31.3082101898603561311809529349519918383".into()),),
        ],
        N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS, MAX_PULL
    );
}
