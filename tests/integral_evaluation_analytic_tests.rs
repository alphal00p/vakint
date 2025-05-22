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
};
use test_utils::{compare_numerical_output, compare_output, get_vakint};
use vakint::{Vakint, VakintSettings};

use crate::test_utils::compare_vakint_evaluation_vs_reference;
use vakint::{externals_from_f64, params_from_f64, vakint_parse};

const N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS: u32 = 32;

#[test_log::test]
fn test_integrate_1l_a() {
    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        use_dot_product_notation: true,
        mu_r_sq_symbol: "mursq".into(),
        integral_normalization_factor: vakint::LoopNormalizationFactor::MSbar,
        run_time_decimal_precision: 16,
        ..VakintSettings::default()
    });

    let mut integral = vakint
        .to_canonical(
            vakint_parse!(
                "(k(1,1)*k(1,2)+k(1,3)*p(1,3))*topo(\
                prop(1,edge(1,1),k(1),muvsq,1)\
            )"
            )
            .unwrap()
            .as_view(),
            true,
        )
        .unwrap();

    integral = vakint.tensor_reduce(integral.as_view()).unwrap();

    let evaluated_integral_res = vakint.evaluate_integral(integral.as_view());
    let evaluated_integral_res_ref = evaluated_integral_res.as_ref();

    let expanded_evaluated_integral_res = evaluated_integral_res_ref.unwrap().expand();
    // let coefs =
    //     expanded_evaluated_integral_res.coefficient_list::<i8>(&[vakint_parse!("vk::ε").unwrap()]);
    // for (v, c) in coefs {
    //     println!("---");
    //     println!("{}: {}", v, c);
    // }
    // panic!();

    let evaluated_integral = compare_output(
        Ok(expanded_evaluated_integral_res.as_view()),
        //evaluated_integral_res_ref.map(|a| a.as_view()),
        simplify_real(
            vakint_parse!(
                "(\
        + ε^-1 * (1/4*𝑖*𝜋^2*muvsq^2*g(1,2))\
        + ε^0  * (3/8*𝑖*𝜋^2*muvsq^2*g(1,2)+1/4*𝑖*𝜋^2*muvsq^2*log(mursq)*g(1,2)-1/4*𝑖*𝜋^2*muvsq^2*log(muvsq)*g(1,2))\
        + ε    * (7/16*𝑖*𝜋^2*muvsq^2*g(1,2)+1/48*𝑖*𝜋^4*muvsq^2*g(1,2)+1/8*𝑖*𝜋^2*muvsq^2*log(mursq)^2*g(1,2)+1/8*𝑖*𝜋^2*muvsq^2*log(muvsq)^2*g(1,2)+3/8*𝑖*𝜋^2*muvsq^2*log(mursq)*g(1,2)-3/8*𝑖*𝜋^2*muvsq^2*log(muvsq)*g(1,2)-1/4*𝑖*𝜋^2*muvsq^2*log(mursq)*log(muvsq)*g(1,2))\
        + ε^2  * (151190863202516241/410199796539607264*𝑖*𝜋^2*muvsq^2*g(1,2)+1/32*𝑖*𝜋^4*muvsq^2*g(1,2)+3/16*𝑖*𝜋^2*muvsq^2*log(mursq)^2*g(1,2)+1/24*𝑖*𝜋^2*muvsq^2*log(mursq)^3*g(1,2)+3/16*𝑖*𝜋^2*muvsq^2*log(muvsq)^2*g(1,2)-1/24*𝑖*𝜋^2*muvsq^2*log(muvsq)^3*g(1,2)+7/16*𝑖*𝜋^2*muvsq^2*log(mursq)*g(1,2)-7/16*𝑖*𝜋^2*muvsq^2*log(muvsq)*g(1,2)+1/48*𝑖*𝜋^4*muvsq^2*log(mursq)*g(1,2)-1/48*𝑖*𝜋^4*muvsq^2*log(muvsq)*g(1,2)-1/8*𝑖*𝜋^2*muvsq^2*log(mursq)^2*log(muvsq)*g(1,2)+1/8*𝑖*𝜋^2*muvsq^2*log(mursq)*log(muvsq)^2*g(1,2)-3/8*𝑖*𝜋^2*muvsq^2*log(mursq)*log(muvsq)*g(1,2))\
            )")
            .unwrap()
            .as_view(),
        ).expand(),
    );

    //debug!("Evaluated integral: {}", evaluated_integral);

    let mut params = HashMap::default();
    params.insert("muvsq".into(), vakint.settings.real_to_prec("1"));
    params.insert("mursq".into(), vakint.settings.real_to_prec("1"));

    let numerical_partial_eval = Vakint::partial_numerical_evaluation(
        &vakint.settings,
        evaluated_integral.as_view(),
        &params,
        None,
    );

    // This test is too unstable as the printout at fixed precision is not accurate enough
    // let numerical_partial_eval_canonical_str = numerical_partial_eval.to_canonical_string();
    // assert_eq!(
    //     numerical_partial_eval_canonical_str,
    //     "-12.0696723514860*g(1,2)*ε^2*𝑖+-5.36845814123893*g(1,2)*𝑖+2.46740110027234*g(1,2)*ε^-1*𝑖+9.41170424471097*g(1,2)*ε*𝑖"
    // );
    assert_eq!(
        numerical_partial_eval.rationalize_coefficients(&Fraction::from(
            0.1_f64.powi((vakint.settings.run_time_decimal_precision - 4) as i32)
        )),
        vakint_parse!("2861033/773022*𝑖*vk::g(1,2)+3640727/573586*𝑖*vk::ε*vk::g(1,2)+1075967/436073*𝑖*vk::ε^-1*vk::g(1,2)+7131281/1067276*𝑖*vk::ε^2*vk::g(1,2)").unwrap()
    );
    let prec = Fraction::from(0.1.pow((vakint.settings.run_time_decimal_precision - 4) as u64));
    _ = compare_output(
        Ok(numerical_partial_eval
            .rationalize_coefficients(&prec)
            .as_view()),
        vakint_parse!(format!(
            "3.701101650408509`{prec}*𝑖*g(1,2)\
                +6.347307988684978`{prec}*𝑖*ε*g(1,2)\
                +2.467401100272340`{prec}*𝑖*ε^-1*g(1,2)\
                +6.681758982674564`{prec}*𝑖*ε^2*g(1,2)",
            prec = vakint.settings.run_time_decimal_precision - 1
        )
        .as_str())
        .unwrap()
        .rationalize_coefficients(&prec),
    );

    params.insert("g(1,2)".into(), vakint.settings.real_to_prec("1"));
    let numerical_full_eval = Vakint::full_numerical_evaluation_without_error(
        &vakint.settings,
        evaluated_integral.as_view(),
        &params,
        None,
    );
    let numerical_full_eval_ref = numerical_full_eval.as_ref();
    debug!(
        "Full eval (metric substituted with 1):\n{}",
        numerical_full_eval_ref.unwrap()
    );
    compare_numerical_output(
        numerical_full_eval_ref,
        vec![
            (-1, ("0.0".into(), "2.467401100272340".into())),
            (0, ("0.0".into(), "3.701101650408509".into())),
            (1, ("0.0".into(), "6.347307988684978".into())),
            (2, ("0.0".into(), "6.681758982674565".into())),
        ],
        vakint.settings.run_time_decimal_precision,
    );
}

#[test_log::test]
fn test_integrate_1l_simple() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{number_of_terms_in_epsilon_expansion: 2, integral_normalization_factor: LoopNormalizationFactor::pySecDec,..VakintSettings::default()},
        EvaluationOrder::analytic_only(),
        vakint_parse!(
            "( 1 )*topo(\
                prop(1,edge(1,1),k(1),muvsq,1)\
            )"
        )
        .unwrap()
        .as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        HashMap::default(),
        vec![
            (-1, ("1.0".into(), "0.0".into()),),
            (0,  ("4.227843350984671393934879099176e-1".into(),  "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_1l_cross_product() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{number_of_terms_in_epsilon_expansion: 5, integral_normalization_factor: LoopNormalizationFactor::MSbar,..VakintSettings::default()},
        EvaluationOrder::analytic_only(),
        vakint_parse!(
            "(k(1,11)*p(1,11)*k(1,12)*p(1,12))*topo(\
                prop(1,edge(1,1),k(1),muvsq,2)\
            )"
        )
        .unwrap()
        .as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=1)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-1, ("0.0".into(), "-15.41829599538179739876641630989".into()),),
            (0,  ("0.0".into(), "-15.41829599538179739876641630989".into()),),
            (1,  ("0.0".into(), "-28.09933616315834914185675426177".into()),),
            (2,  ("0.0".into(), "-21.92144645108947660345701527424".into()),),
            (3,  ("0.0".into(), "-31.30821018986035613118095293495".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_1l_cross_product_with_additional_symbols_numerator() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{number_of_terms_in_epsilon_expansion: 5, integral_normalization_factor: LoopNormalizationFactor::MSbar,..VakintSettings::default()},
        EvaluationOrder::analytic_only(),
        vakint_parse!(
            "(user_space::A*k(1,11)*p(1,11)*k(1,12)*p(1,12)+user_space::B)*topo(\
                prop(1,edge(1,1),k(1),muvsq,2)\
            )"
        )
        .unwrap()
        .as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0), ("user_space::A".into(), 3.0), ("user_space::B".into(), 4.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=1)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-1, ("0.0".into(), "-6.776470381787957720961284930165".into()),),
            (0,  ("0.0".into(), "-46.25488798614539219629924892967".into()),),
            (1,  ("0.0".into(), "-51.82831147814090168009015188908".into()),),
            (2,  ("0.0".into(), "-81.58277415564679875594772935461".into()),),
            (3,  ("0.0".into(), "-69.88990073019845746778737271354".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_1l_dot_product_external() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{integral_normalization_factor: LoopNormalizationFactor::MSbar,..VakintSettings::default()},
        EvaluationOrder::analytic_only(),
        vakint_parse!(
            "(k(1,1)*p(1,1)*k(1,2)*p(2,2))*topo(\
                prop(1,edge(1,1),k(1),muvsq,1)\
            )"
        )
        .unwrap()
        .as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-1, ("0.0".into(), "-7298.572454605580698628106094408".into()),),
            (0,  ("0.0".into(),  "-10947.85868190837104794215914161".into()),),
            (1,  ("0.0".into(), "-18775.33703053016641729482116716".into()),),
            (2,  ("0.0".into(),  "-19764.64307075136314315766281197".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_2l() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{integral_normalization_factor: LoopNormalizationFactor::MSbar,..VakintSettings::default()},
        EvaluationOrder::analytic_only(),
        vakint_parse!(
            "(1)*topo(\
                prop(1,edge(1,2),k(1),muvsq,1)\
                *prop(2,edge(1,2),k(2),muvsq,1)\
                *prop(3,edge(2,1),k(1)+k(2),muvsq,1)\
            )"
        )
        .unwrap()
        .as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=1)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-2, ("-146.1136365510036558546604990331".into(), "0.0".into()),),
            (-1, ("-438.3409096530109675639814970992".into(),  "0.0".into()),),
            (0,  ("-920.6659438677135293597636718376".into(), "0.0".into()),),
            (1,  ("-2358.169673453033244370977417874".into(),  "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_3l() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{integral_normalization_factor: LoopNormalizationFactor::MSbar,..VakintSettings::default()},
        EvaluationOrder::alphaloop_only(),
        vakint_parse!(
            "(1)*topo(\
                 prop(1,edge(1,2),k(1),muvsq,1)\
                *prop(2,edge(2,3),k(2),muvsq,1)\
                *prop(3,edge(3,1),k(3),muvsq,1)\
                *prop(4,edge(1,4),k(3)-k(1),muvsq,1)\
                *prop(5,edge(2,4),k(1)-k(2),muvsq,1)\
                *prop(6,edge(3,4),k(2)-k(3),muvsq,1)\
            )"
        )
        .unwrap()
        .as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=1)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-1, ("0.0".into(), "-2311.289033520460340396770711738".into()),),
            ( 0, ("0.0".into(),  "9647.808284968523284819072033817".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_3l_rank_4() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{integral_normalization_factor: LoopNormalizationFactor::MSbar,..VakintSettings::default()},
        EvaluationOrder::alphaloop_only(),
        vakint_parse!(
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
            )"
        ).unwrap().as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 2.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-3, ("0.0".into(),   "4975068.903103003548576756843470".into()),),
            (-2, ("0.0".into(),  "38893977.28458284777498997027566".into()),),
            (-1, ("0.0".into(),   "172684280.6796404510100040673463".into()),),
            ( 0, ("0.0".into(),  "533601880.7076116546200537059228".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_3l_rank_4_additional_symbols_numerator() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{integral_normalization_factor: LoopNormalizationFactor::MSbar,..VakintSettings::default()},
        EvaluationOrder::alphaloop_only(),
        vakint_parse!(
            "(
                  user_space::A*k(1,11)*k(2,11)*k(1,22)*k(2,22)
                + user_space::B*p(1,11)*k(3,11)*k(3,22)*p(2,22)
                + user_space::C*p(1,11)*p(2,11)*(k(2,22)+k(1,22))*k(2,22) 
             )
            *topo(\
                 prop(1,edge(1,2),k(1),muvsq,1)\
                *prop(2,edge(2,3),k(2),muvsq,1)\
                *prop(3,edge(3,1),k(3),muvsq,1)\
                *prop(4,edge(1,4),k(3)-k(1),muvsq,1)\
                *prop(5,edge(2,4),k(1)-k(2),muvsq,1)\
                *prop(6,edge(3,4),k(2)-k(3),muvsq,1)\
            )"
        ).unwrap().as_view(),
        params_from_f64(&[
            ("muvsq".into(), 1.0), ("mursq".into(), 2.0),
            ("user_space::A".into(), 5.0), ("user_space::B".into(), 7.0), ("user_space::C".into(), 11.0)
            ].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-3, ("0.0".into(),   "51891342.24417464772786998080343".into()),),
            (-2, ("0.0".into(),  "404449675.4942040064883673562925".into()),),
            (-1, ("0.0".into(),   "1791033786.227345613037700936867".into()),),
            ( 0, ("0.0".into(),  "5518874647.019060130141182226994".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_3l_rank_4_matad() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{integral_normalization_factor: LoopNormalizationFactor::MSbar, number_of_terms_in_epsilon_expansion: 5,..VakintSettings::default()},
        EvaluationOrder::matad_only(Some(MATADOptions {direct_numerical_substition: true,..MATADOptions::default()})),
        vakint_parse!(
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
            )"
        ).unwrap().as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 2.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-3, ("0.0".into(),   "4975068.903103003548576756843470".into()),),
            (-2, ("0.0".into(),   "38893977.28458284777498997027566".into()),),
            (-1, ("0.0".into(),   "172684280.6796404510100040673463".into()),),
            ( 0, ("0.0".into(),   "533601880.7076116546200537059229".into()),),
            ( 1, ("0.0".into(),   "1931790768.728074898808769939759".into()),),
            ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_3l_rank_4_matad_additional_symbols_numerator() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{integral_normalization_factor: LoopNormalizationFactor::MSbar, number_of_terms_in_epsilon_expansion: 5,..VakintSettings::default()},
        EvaluationOrder::matad_only(Some(MATADOptions {direct_numerical_substition: true,..MATADOptions::default()})),
        vakint_parse!(
            "(
                  user_space::A*k(1,11)*k(2,11)*k(1,22)*k(2,22)
                + user_space::B*p(1,11)*k(3,11)*k(3,22)*p(2,22)
                + user_space::C*p(1,11)*p(2,11)*(k(2,22)+k(1,22))*k(2,22) 
             )
            *topo(\
                 prop(1,edge(1,2),k(1),muvsq,1)\
                *prop(2,edge(2,3),k(2),muvsq,1)\
                *prop(3,edge(3,1),k(3),muvsq,1)\
                *prop(4,edge(1,4),k(3)-k(1),muvsq,1)\
                *prop(5,edge(2,4),k(1)-k(2),muvsq,1)\
                *prop(6,edge(3,4),k(2)-k(3),muvsq,1)\
            )"
        ).unwrap().as_view(),
        params_from_f64(&[
            ("muvsq".into(), 1.0), ("mursq".into(), 2.0),
            ("user_space::A".into(), 5.0), ("user_space::B".into(), 7.0), ("user_space::C".into(), 11.0)
            ].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-3, ("0.0".into(),   "51891342.24417464772786998080343".into()),),
            (-2, ("0.0".into(),   "404449675.4942040064883673562925".into()),),
            (-1, ("0.0".into(),   "1791033786.227345613037700936867".into()),),
            ( 0, ("0.0".into(),   "5518874647.019060130141182226994".into()),),
            ( 1, ( "0.0".into(),  "19997559685.62917299165116831086".into()),),
            ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_3l_matad() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{integral_normalization_factor: LoopNormalizationFactor::MSbar, number_of_terms_in_epsilon_expansion: 5,..VakintSettings::default()},
        EvaluationOrder::matad_only(None),
        vakint_parse!(
            "(1)*topo(\
                 prop(1,edge(1,2),k(1),muvsq,1)\
                *prop(2,edge(2,3),k(2),muvsq,1)\
                *prop(3,edge(3,1),k(3),muvsq,1)\
                *prop(4,edge(1,4),k(3)-k(1),muvsq,1)\
                *prop(5,edge(2,4),k(1)-k(2),muvsq,1)\
                *prop(6,edge(3,4),k(2)-k(3),muvsq,1)\
            )"
        )
        .unwrap()
        .as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=1)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-1, ("0.0".into(), "-2311.289033520460340396770711738".into()),),
            ( 0, ("0.0".into(),  "9647.808284968523284819072033817".into()),),
            ( 1, ("0.0".into(),  "-40259.80884519899852832255821372".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_4l_h() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { integral_normalization_factor: LoopNormalizationFactor::MSbar, number_of_terms_in_epsilon_expansion: 5, run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "(1)*topo(\
                 prop(1,edge(5,1),k(1),muvsq,1)\
                *prop(2,edge(2,6),k(2),muvsq,1)\
                *prop(3,edge(6,5),k(3),muvsq,1)\
                *prop(4,edge(3,4),k(4),muvsq,1)\
                *prop(5,edge(4,5),k(1)-k(3),muvsq,1)\
                *prop(6,edge(6,3),k(2)-k(3),muvsq,1)\
                *prop(7,edge(4,1),k(3)-k(1)+k(4),muvsq,1)\
                *prop(8,edge(2,3),k(3)-k(2)+k(4),muvsq,1)\
                *prop(9,edge(1,2),k(3)+k(4),muvsq,1)\
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=1)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            // This does not have an analytical expression yet (FMFT not implemented yet)
            (0,  ("-12799.53514305961548130719263292".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_4l_h_rank_4() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { integral_normalization_factor: LoopNormalizationFactor::MSbar, number_of_terms_in_epsilon_expansion: 5, run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "(
                  k(1,11)*k(2,11)*k(1,22)*k(2,22)
                + p(1,11)*k(3,11)*k(3,22)*p(2,22)
                + p(1,11)*p(2,11)*(k(2,22)+k(1,22))*k(2,22) 
             )*topo(\
                 prop(1,edge(5,1),k(1),muvsq,1)\
                *prop(2,edge(2,6),k(2),muvsq,1)\
                *prop(3,edge(6,5),k(3),muvsq,1)\
                *prop(4,edge(3,4),k(4),muvsq,1)\
                *prop(5,edge(4,5),k(1)-k(3),muvsq,1)\
                *prop(6,edge(6,3),k(2)-k(3),muvsq,1)\
                *prop(7,edge(4,1),k(3)-k(1)+k(4),muvsq,1)\
                *prop(8,edge(2,3),k(3)-k(2)+k(4),muvsq,1)\
                *prop(9,edge(1,2),k(3)+k(4),muvsq,1)\
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 3.0), ("mursq".into(), 5.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-4,  ("+10674.59739307939575801964744189".into(), "0.0".into()),),
            (-3,  ("+109876.8599799874165101526102648".into(), "0.0".into()),),
            (-2,  ("+523100.2201190291672501506043102".into(), "0.0".into()),),
            (-1,  ("+2398035.283172097754065659367279".into(), "0.0".into()),),
            ( 0,  ("+45463747.08560375418644783085990".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_4l_h_rank_4_additional_symbols_numerator() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { integral_normalization_factor: LoopNormalizationFactor::MSbar, number_of_terms_in_epsilon_expansion: 5, run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "(
                  user_space::A*k(1,11)*k(2,11)*k(1,22)*k(2,22)
                + user_space::B*p(1,11)*k(3,11)*k(3,22)*p(2,22)
                + user_space::C*p(1,11)*p(2,11)*(k(2,22)+k(1,22))*k(2,22) 
             )*topo(\
                 prop(1,edge(5,1),k(1),muvsq,1)\
                *prop(2,edge(2,6),k(2),muvsq,1)\
                *prop(3,edge(6,5),k(3),muvsq,1)\
                *prop(4,edge(3,4),k(4),muvsq,1)\
                *prop(5,edge(4,5),k(1)-k(3),muvsq,1)\
                *prop(6,edge(6,3),k(2)-k(3),muvsq,1)\
                *prop(7,edge(4,1),k(3)-k(1)+k(4),muvsq,1)\
                *prop(8,edge(2,3),k(3)-k(2)+k(4),muvsq,1)\
                *prop(9,edge(1,2),k(3)+k(4),muvsq,1)\
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 3.0), ("mursq".into(), 5.0), ("user_space::A".into(), 5.0), ("user_space::B".into(), 7.0), ("user_space::C".into(), 11.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-4,  ("+53372.98696539697879009823720947".into(), "0.0".into()),),
            (-3,  ("+549384.2998999370825507630513242".into(), "0.0".into()),),
            (-2,  ("+2207090.697342127753994656146542".into(), "0.0".into()),),
            (-1,  ("+12697906.62436134259734910568075".into(), "0.0".into()),),
            ( 0,  ("+412404431.7658464058084812513050".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
#[allow(non_snake_case)]
fn test_integrate_4l_PR9d_from_H() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD, number_of_terms_in_epsilon_expansion: 5, 
            run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "( 1 )*topo(\
                 prop(1,edge(5,1),k(1),muvsq,1)\
                *prop(2,edge(2,6),k(2),muvsq,1)\
                *prop(3,edge(6,5),k(3),muvsq,0)\
                *prop(4,edge(3,4),k(4),muvsq,2)\
                *prop(5,edge(4,5),k(1)-k(3),muvsq,1)\
                *prop(6,edge(6,3),k(2)-k(3),muvsq,1)\
                *prop(7,edge(4,1),k(3)-k(1)+k(4),muvsq,1)\
                *prop(8,edge(2,3),k(3)-k(2)+k(4),muvsq,1)\
                *prop(9,edge(1,2),k(3)+k(4),muvsq,0)\
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-4,  ("8.333333333333333333333333333333e-2".into(), "0.0".into()),),
            (-3,  ("3.333333333333333333333333333333e-1".into(), "0.0".into()),),
            (-2,  ("-3.144646082033583725553786166618e-1".into(), "0.0".into()),),
            (-1,  ("5.421352941798334340259377610275".into(), "0.0".into()),),
            ( 0,  ("-28.31064373017674207211847384976".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
#[allow(non_snake_case)]
fn test_integrate_4l_PR9d_from_X() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD, number_of_terms_in_epsilon_expansion: 5, 
            run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "( 1 )*topo(\
                prop(1,edge(5,1),k(1),muvsq,1)\
                *prop(2,edge(2,6),k(2),muvsq,1)\
                *prop(3,edge(6,5),k(3),muvsq,0)\
                *prop(4,edge(4,3),k(4),muvsq,2)\
                *prop(5,edge(3,5),k(1)-k(3),muvsq,1)\
                *prop(6,edge(6,4),k(2)-k(3),muvsq,1)\
                *prop(7,edge(3,2),k(3)-k(1)+k(4),muvsq,1)\
                *prop(8,edge(1,4),k(3)-k(2)+k(4),muvsq,1)\
                *prop(9,edge(2,1),k(3)-k(1)-k(2)+k(4),muvsq,0)\
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-4,  ("8.333333333333333333333333333333e-2".into(), "0.0".into()),),
            (-3,  ("3.333333333333333333333333333333e-1".into(), "0.0".into()),),
            (-2,  ("-3.144646082033583725553786166618e-1".into(), "0.0".into()),),
            (-1,  ("5.421352941798334340259377610275".into(), "0.0".into()),),
            ( 0,  ("-28.31064373017674207211847384976".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
#[allow(non_snake_case)]
fn test_integrate_4l_PR9d_from_H_pinch() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD, number_of_terms_in_epsilon_expansion: 5, 
            run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "( 1 )*topo(\
                 prop(1,edge(5,2),k(1),muvsq,1)\
                *prop(2,edge(2,5),k(2),muvsq,1)\
                *prop(4,edge(3,4),k(3),muvsq,2)\
                *prop(5,edge(4,5),k(4),muvsq,1)\
                *prop(6,edge(5,3),k(4)+k(2)-k(1),muvsq,1)\
                *prop(7,edge(4,2),k(3)-k(4),muvsq,1)\
                *prop(8,edge(2,3),k(1)-k(2)+k(3)-k(4),muvsq,1)\
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-4,  ("8.333333333333333333333333333333e-2".into(), "0.0".into()),),
            (-3,  ("3.333333333333333333333333333333e-1".into(), "0.0".into()),),
            (-2,  ("-3.144646082033583725553786166618e-1".into(), "0.0".into()),),
            (-1,  ("5.421352941798334340259377610275".into(), "0.0".into()),),
            ( 0,  ("-28.31064373017674207211847384976".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
#[allow(non_snake_case)]
fn test_integrate_4l_PR9d_from_FG() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD, number_of_terms_in_epsilon_expansion: 5, 
            run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "( 1 )*topo(\
                 prop(1,edge(5,3),k(1),muvsq,1)\
                *prop(2,edge(3,4),k(2),muvsq,2)\
                *prop(3,edge(4,5),k(3),muvsq,1)\
                *prop(4,edge(2,1),k(1)-k(3),muvsq,0)\
                *prop(5,edge(5,1),k(4),muvsq,1)\
                *prop(6,edge(4,2),k(2)-k(3),muvsq,1)\
                *prop(7,edge(1,5),k(1)-k(3)+k(4),muvsq,1)\
                *prop(8,edge(3,2),k(1)-k(2),muvsq,1)\
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-4,  ("8.333333333333333333333333333333e-2".into(), "0.0".into()),),
            (-3,  ("3.333333333333333333333333333333e-1".into(), "0.0".into()),),
            (-2,  ("-3.144646082033583725553786166618e-1".into(), "0.0".into()),),
            (-1,  ("5.421352941798334340259377610275".into(), "0.0".into()),),
            ( 0,  ("-28.31064373017674207211847384976".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
#[allow(non_snake_case)]
fn test_integrate_4l_PR9d_from_FG_pinch() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD, number_of_terms_in_epsilon_expansion: 5, 
                         run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "( 1 )*topo(\
                 prop(1,edge(5,3),k(1),muvsq,1)\
                *prop(2,edge(3,4),k(2),muvsq,2)\
                *prop(3,edge(4,5),k(3),muvsq,1)\
                *prop(5,edge(5,1),k(4),muvsq,1)\
                *prop(6,edge(4,1),k(2)-k(3),muvsq,1)\
                *prop(7,edge(1,5),k(1)-k(3)+k(4),muvsq,1)\
                *prop(8,edge(3,1),k(1)-k(2),muvsq,1)\
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-4,  ("8.333333333333333333333333333333e-2".into(), "0.0".into()),),
            (-3,  ("3.333333333333333333333333333333e-1".into(), "0.0".into()),),
            (-2,  ("-3.144646082033583725553786166618e-1".into(), "0.0".into()),),
            (-1,  ("5.421352941798334340259377610275".into(), "0.0".into()),),
            ( 0,  ("-28.31064373017674207211847384976".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
#[allow(non_snake_case)]
fn test_integrate_4l_PR11d() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { number_of_terms_in_epsilon_expansion: 5, run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "( 1 )*topo(\
                 prop(1,edge(1,2),k(1),muvsq,2)\
                *prop(2,edge(2,5),k(2),muvsq,1)\
                *prop(3,edge(3,4),k(3),muvsq,1)\
                *prop(4,edge(4,5),k(4),muvsq,1)\
                *prop(5,edge(2,3),k(1)-k(2),muvsq,1)\
                *prop(6,edge(4,1),k(3)-k(4),muvsq,1)\
                *prop(7,edge(5,3),k(2)+k(3)-k(1),muvsq,1)\
                *prop(8,edge(1,5),k(3)-k(4)-k(1),muvsq,1)\
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (0,  ("-2.906486288643112641819206002127".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
#[allow(non_snake_case)]
fn test_integrate_4l_clover() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { 
            integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
            number_of_terms_in_epsilon_expansion: 5, run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "( 1 )*topo(
                prop(1, edge(1, 1), k(1), muvsq, 1)*\
                prop(2, edge(1, 1), k(2), muvsq, 1)*\
                prop(3, edge(1, 1), k(3), muvsq, 1)*\
                prop(4, edge(1, 1), k(4), muvsq, 1)
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-4,  ("1.000000000000000000000000000000".into(), "0.0".into()),),
            (-3,  ("4.000000000000000000000000000000".into(), "0.0".into()),),
            (-2,  ("13.28986813369645287294483033329".into(), "0.0".into()),),
            (-1,  ("31.55672999723968577791300378449".into(), "0.0".into()),),
            (0,   ("67.98165058904685502307905531744".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
#[allow(non_snake_case)]
fn test_integrate_4l_clover_with_non_unit_scales() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { number_of_terms_in_epsilon_expansion: 5, run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "( 1 )*topo(
                prop(1, edge(1, 1), k(1), muvsq, 1)*\
                prop(2, edge(1, 1), k(2), muvsq, 1)*\
                prop(3, edge(1, 1), k(3), muvsq, 1)*\
                prop(4, edge(1, 1), k(4), muvsq, 1)
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 3.0), ("mursq".into(), 7.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-4,  ("81.00000000000000000000000000000".into(), "0.0".into()),),
            (-3,  ("-218.9682569565641868485693739496".into(), "0.0".into()),),
            (-2,  ("724.4490568207455247307151297051".into(), "0.0".into()),),
            (-1,  ("-1446.834846767122729283863130264".into(), "0.0".into()),),
            (0,   ("3106.843653546628093141699453745".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
#[allow(non_snake_case)]
fn test_integrate_4l_dotted_clover() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { number_of_terms_in_epsilon_expansion: 5, run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "( 1 )*topo(
                prop(1, edge(1, 1), k(1), muvsq, 2)*\
                prop(2, edge(1, 1), k(2), muvsq, 1)*\
                prop(3, edge(1, 1), k(3), muvsq, 1)*\
                prop(4, edge(1, 1), k(4), muvsq, 1)
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 3.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-4,  ("27.00000000000000000000000000000".into(), "0.0".into()),),
            (-3,  ("-99.98941898552139561618979131653".into(), "0.0".into()),),
            (-2,  ("314.4724379257699038597615012182".into(), "0.0".into()),),
            (-1,  ("-723.7613011959560846715260866565".into(), "0.0".into()),),
            (0,   ("1517.892833437916940808520861336".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
#[allow(non_snake_case)]
fn test_integrate_4l_clover_with_numerator() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { 
            number_of_terms_in_epsilon_expansion: 5, 
            run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS,
            integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD, 
            ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "( 
                  user_space::A * k(1,11)*k(2,11)*k(1,22)*k(2,22)
                + user_space::B * p(1,11)*k(3,11)*k(3,22)*p(2,22)
                + user_space::C * p(1,11)*p(2,11)*(k(2,22)+k(1,22))*k(2,22) 
             )*topo(
                prop(1, edge(1, 1), k(1), muvsq, 2)*\
                prop(2, edge(1, 1), k(2), muvsq, 1)*\
                prop(3, edge(1, 1), k(3), muvsq, 1)*\
                prop(4, edge(1, 1), k(4), muvsq, 1)
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 0.3), ("mursq".into(), 0.7), ("user_space::A".into(), 3.0), ("user_space::B".into(), 4.0), ("user_space::C".into(), 5.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-4,  ("-1.897149599999999855007182247846e-1".into(), "0.0".into()),),
            (-3,  ("-1.495259819655131009380566817668".into(), "0.0".into()),),
            (-2,  ("-6.805240907875078933713325181389".into(), "0.0".into()),),
            (-1,  ("-22.56027900679456203938234477552".into(), "0.0".into()),),
            (0,   ("-60.49337040949871593265194938449".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}
