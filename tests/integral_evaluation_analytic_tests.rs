mod test_utils;
use vakint::{utils::simplify_real, EvaluationOrder};

use std::vec;

use log::debug;
use std::collections::HashMap;
use symbolica::{
    atom::Atom,
    domains::{float::NumericalFloatLike, rational::Fraction},
};
use test_utils::{compare_numerical_output, compare_output, convert_test_params, get_vakint};
use vakint::{Vakint, VakintSettings};

use crate::test_utils::{compare_vakint_evaluation_vs_reference, convert_test_externals};

const N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS: u32 = 32;

#[test_log::test]
fn test_integrate_1l_a() {
    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        use_dot_product_notation: true,
        mu_r_sq_symbol: "mursq".into(),
        integral_normalization_factor: vakint::LoopNormalizationFactor::MSbar,
        n_digits_at_evaluation_time: 16,
        ..VakintSettings::default()
    });

    let mut integral = vakint
        .to_canonical(
            Atom::parse(
                "(k(1,1)*k(1,2)+k(1,3)*p(1,3))*topo(\
                prop(1,edge(1,1),k(1),muvsq,1)\
            )",
            )
            .unwrap()
            .as_view(),
            true,
        )
        .unwrap();

    integral = vakint.tensor_reduce(integral.as_view()).unwrap();

    let evaluated_integral_res = vakint.evaluate_integral(integral.as_view());
    let evaluated_integral_res_ref = evaluated_integral_res.as_ref();

    // let coefs = evaluated_integral_res_ref
    //     .unwrap()
    //     .coefficient_list(State::get_symbol("ε"));
    // for (v, c) in coefs.0 {
    //     println!("{}: {}", v, c);
    // }
    // println!("ε^0: {}", coefs.1);
    let evaluated_integral = compare_output(
        evaluated_integral_res_ref.map(|a| a.as_view()),
        simplify_real(Atom::parse(
            "(\
        + ε^-1 * (1/4*𝑖*𝜋^2*muvsq^2*g(1,2))\
        + ε^-0 * (1/4*muvsq^2*(-𝑖*𝜋^2*log(𝜋)+𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))*g(1,2)+3/8*𝑖*𝜋^2*muvsq^2*g(1,2)+1/4*𝑖*𝜋^2*muvsq^2*log(muvsq^-1)*g(1,2)+1/4*𝑖*𝜋^2*muvsq^2*log(exp(-EulerGamma))*g(1,2))\
        + ε    * ((1/4*muvsq^2*(-𝑖*𝜋^2*log(𝜋)+𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))*g(1,2)+3/8*𝑖*𝜋^2*muvsq^2*g(1,2))*log(exp(-EulerGamma))+(1/4*muvsq^2*(-𝑖*𝜋^2*log(𝜋)+𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))*g(1,2)+3/8*𝑖*𝜋^2*muvsq^2*g(1,2)+1/4*𝑖*𝜋^2*muvsq^2*log(exp(-EulerGamma))*g(1,2))*log(muvsq^-1)+𝑖*𝜋^2*(7/16*muvsq^2*g(1,2)+1/48*𝜋^2*muvsq^2*g(1,2))+3/8*muvsq^2*(-𝑖*𝜋^2*log(𝜋)+𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))*g(1,2)+1/4*muvsq^2*(1/2*𝑖*𝜋^2*log(𝜋)^2+1/2*𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1)^2-𝑖*𝜋^2*log(𝜋)*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))*g(1,2)+1/8*𝑖*𝜋^2*muvsq^2*log(muvsq^-1)^2*g(1,2)+1/8*𝑖*𝜋^2*muvsq^2*log(exp(-EulerGamma))^2*g(1,2))\
        + ε^2  * ((7/16*muvsq^2*g(1,2)+1/48*𝜋^2*muvsq^2*g(1,2))*(-𝑖*𝜋^2*log(𝜋)+𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))+1/2*(1/4*muvsq^2*(-𝑖*𝜋^2*log(𝜋)+𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))*g(1,2)+3/8*𝑖*𝜋^2*muvsq^2*g(1,2))*log(exp(-EulerGamma))^2+(𝑖*𝜋^2*(7/16*muvsq^2*g(1,2)+1/48*𝜋^2*muvsq^2*g(1,2))+3/8*muvsq^2*(-𝑖*𝜋^2*log(𝜋)+𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))*g(1,2)+1/4*muvsq^2*(1/2*𝑖*𝜋^2*log(𝜋)^2+1/2*𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1)^2-𝑖*𝜋^2*log(𝜋)*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))*g(1,2))*log(exp(-EulerGamma))+1/2*(1/4*muvsq^2*(-𝑖*𝜋^2*log(𝜋)+𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))*g(1,2)+3/8*𝑖*𝜋^2*muvsq^2*g(1,2)+1/4*𝑖*𝜋^2*muvsq^2*log(exp(-EulerGamma))*g(1,2))*log(muvsq^-1)^2+((1/4*muvsq^2*(-𝑖*𝜋^2*log(𝜋)+𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))*g(1,2)+3/8*𝑖*𝜋^2*muvsq^2*g(1,2))*log(exp(-EulerGamma))+𝑖*𝜋^2*(7/16*muvsq^2*g(1,2)+1/48*𝜋^2*muvsq^2*g(1,2))+3/8*muvsq^2*(-𝑖*𝜋^2*log(𝜋)+𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))*g(1,2)+1/4*muvsq^2*(1/2*𝑖*𝜋^2*log(𝜋)^2+1/2*𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1)^2-𝑖*𝜋^2*log(𝜋)*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))*g(1,2)+1/8*𝑖*𝜋^2*muvsq^2*log(exp(-EulerGamma))^2*g(1,2))*log(muvsq^-1)+𝑖*𝜋^2*(151190863202516241/410199796539607264*muvsq^2*g(1,2)+1/32*𝜋^2*muvsq^2*g(1,2))+3/8*muvsq^2*(1/2*𝑖*𝜋^2*log(𝜋)^2+1/2*𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1)^2-𝑖*𝜋^2*log(𝜋)*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))*g(1,2)+1/4*muvsq^2*(-1/6*𝑖*𝜋^2*log(𝜋)^3+1/6*𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1)^3+1/2*𝑖*𝜋^2*log(𝜋)^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1)-1/2*𝑖*𝜋^2*log(𝜋)*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1)^2)*g(1,2)+1/24*𝑖*𝜋^2*muvsq^2*log(muvsq^-1)^3*g(1,2)+1/24*𝑖*𝜋^2*muvsq^2*log(exp(-EulerGamma))^3*g(1,2))\
        )",
        )
        .unwrap().as_view()),
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
    debug!(
        "Partial eval: {}",
        numerical_partial_eval.to_canonical_string()
    );

    // This test is too unstable as the printout at fixed precision is not accurate enough
    // let numerical_partial_eval_canonical_str = numerical_partial_eval.to_canonical_string();
    // assert_eq!(
    //     numerical_partial_eval_canonical_str,
    //     "-12.0696723514860*g(1,2)*ε^2*𝑖+-5.36845814123893*g(1,2)*𝑖+2.46740110027234*g(1,2)*ε^-1*𝑖+9.41170424471097*g(1,2)*ε*𝑖"
    // );
    assert_eq!(
        numerical_partial_eval.rationalize_coefficients(&Fraction::from(
            0.1_f64.powi((vakint.settings.n_digits_at_evaluation_time - 4) as i32)
        )),
        Atom::parse("-2879700/536411*𝑖*g(1,2)+3726809/395976*𝑖*ε*g(1,2)+1075967/436073*𝑖*ε^-1*g(1,2)-4041047/334810*𝑖*ε^2*g(1,2)").unwrap()
    );

    let prec = Fraction::from(0.1.pow((vakint.settings.n_digits_at_evaluation_time - 4) as u64));
    compare_output(
        Ok(numerical_partial_eval
            .rationalize_coefficients(&prec)
            .as_view()),
        Atom::parse(
            format!(
                "-5.36845814123893`{prec}*𝑖*g(1,2)\
                +9.4117042447109682`{prec}*𝑖*ε*g(1,2)\
                +2.46740110027234`{prec}*𝑖*ε^-1*g(1,2)\
                -12.0696723514860`{prec}*𝑖*ε^2*g(1,2)",
                prec = vakint.settings.n_digits_at_evaluation_time - 1
            )
            .as_str(),
        )
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
            (-1, ("0.0".into(), "2.46740110027234".into())),
            (0, ("0.0".into(), "-5.36845814123893".into())),
            (1, ("0.0".into(), "9.41170424471097".into())),
            (2, ("0.0".into(), "-12.0696723514860".into())),
        ],
        vakint.settings.n_digits_at_evaluation_time,
    );
}

#[test_log::test]
fn test_integrate_1l_dot_product_external() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings::default(),
        EvaluationOrder::analytic_only(),
        Atom::parse(
            "(k(1,1)*p(1,1)*k(1,2)*p(2,2))*topo(\
                prop(1,edge(1,1),k(1),muvsq,1)\
            )",
        )
        .unwrap()
        .as_view(),
        convert_test_params(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        convert_test_externals(
        &(1..=2)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-1, ("0.0".into(), "-7298.572454605580698628106094408".into()),),
            (0,  ("0.0".into(),  "15879.89918178474997676561014674".into()),),
            (1,  ("0.0".into(), "-27839.82115585504758906774113191".into()),),
            (2,  ("0.0".into(),  "35702.09081569554005452137401203".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_2l() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings::default(),
        EvaluationOrder::analytic_only(),
        Atom::parse(
            "(1)*topo(\
                prop(1,edge(1,2),k(1),muvsq,1)\
                *prop(2,edge(1,2),k(2),muvsq,1)\
                *prop(3,edge(2,1),k(1)+k(2),muvsq,1)\
            )",
        )
        .unwrap()
        .as_view(),
        convert_test_params(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        convert_test_externals(
        &(1..=1)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-2, ("-146.1136365510036558546604990331".into(), "0.0".into()),),
            (-1, ("635.8146971740286947808753047759".into(),  "0.0".into()),),
            (0,  ("-1646.531034471454483109201793220".into(), "0.0".into()),),
            (1,  ("2240.516116133454318298096346441".into(),  "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_3l() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings::default(),
        EvaluationOrder::alphaloop_only(),
        Atom::parse(
            "(1)*topo(\
                 prop(1,edge(1,2),k(1),muvsq,1)\
                *prop(2,edge(2,3),k(2),muvsq,1)\
                *prop(3,edge(3,1),k(3),muvsq,1)\
                *prop(4,edge(1,4),k(3)-k(1),muvsq,1)\
                *prop(5,edge(2,4),k(1)-k(2),muvsq,1)\
                *prop(6,edge(3,4),k(2)-k(3),muvsq,1)\
            )",
        )
        .unwrap()
        .as_view(),
        convert_test_params(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        convert_test_externals(
        &(1..=1)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-1, ("0.0".into(), "-2311.289033520460340396770711738".into()),),
            ( 0, ("0.0".into(),  "35134.99893627257345553503414002".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_3l_rank_4() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings::default(),
        EvaluationOrder::alphaloop_only(),
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
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        convert_test_externals(
        &(1..=2)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-3, ("0.0".into(),   "4975068.903103003548576756843470".into()),),
            (-2, ("0.0".into(),  "-15967412.96033300288621485195252".into()),),
            (-1, ("0.0".into(),   "46275660.33034806160550888444548".into()),),
            ( 0, ("0.0".into(),  "-117731367.8844665539198405383934".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}
