mod test_utils;

use std::vec;

use log::debug;
use std::collections::HashMap;
use symbolica::{
    atom::Atom,
    domains::{float::NumericalFloatLike, rational::Fraction},
};
use test_utils::{compare_numerical_output, compare_output, get_vakint};
use vakint::{Vakint, VakintSettings};

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
        Atom::parse(
            "(\
        + ε^-1 * (1/4*𝑖*𝜋^2*muvsq^2*g(1,2))\
        + ε^-0 * (1/4*muvsq^2*(-𝑖*𝜋^2*log(𝜋)+𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))*g(1,2)+3/8*𝑖*𝜋^2*muvsq^2*g(1,2)+1/4*𝑖*𝜋^2*muvsq^2*log(muvsq^-1)*g(1,2)+1/4*𝑖*𝜋^2*muvsq^2*log(exp(-EulerGamma))*g(1,2))\
        + ε    * ((1/4*muvsq^2*(-𝑖*𝜋^2*log(𝜋)+𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))*g(1,2)+3/8*𝑖*𝜋^2*muvsq^2*g(1,2))*log(exp(-EulerGamma))+(1/4*muvsq^2*(-𝑖*𝜋^2*log(𝜋)+𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))*g(1,2)+3/8*𝑖*𝜋^2*muvsq^2*g(1,2)+1/4*𝑖*𝜋^2*muvsq^2*log(exp(-EulerGamma))*g(1,2))*log(muvsq^-1)+𝑖*𝜋^2*(7/16*muvsq^2*g(1,2)+1/48*𝜋^2*muvsq^2*g(1,2))+3/8*muvsq^2*(-𝑖*𝜋^2*log(𝜋)+𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))*g(1,2)+1/4*muvsq^2*(1/2*𝑖*𝜋^2*log(𝜋)^2+1/2*𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1)^2-𝑖*𝜋^2*log(𝜋)*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))*g(1,2)+1/8*𝑖*𝜋^2*muvsq^2*log(muvsq^-1)^2*g(1,2)+1/8*𝑖*𝜋^2*muvsq^2*log(exp(-EulerGamma))^2*g(1,2))\
        + ε^2  * ((7/16*muvsq^2*g(1,2)+1/48*𝜋^2*muvsq^2*g(1,2))*(-𝑖*𝜋^2*log(𝜋)+𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))+1/2*(1/4*muvsq^2*(-𝑖*𝜋^2*log(𝜋)+𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))*g(1,2)+3/8*𝑖*𝜋^2*muvsq^2*g(1,2))*log(exp(-EulerGamma))^2+(𝑖*𝜋^2*(7/16*muvsq^2*g(1,2)+1/48*𝜋^2*muvsq^2*g(1,2))+3/8*muvsq^2*(-𝑖*𝜋^2*log(𝜋)+𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))*g(1,2)+1/4*muvsq^2*(1/2*𝑖*𝜋^2*log(𝜋)^2+1/2*𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1)^2-𝑖*𝜋^2*log(𝜋)*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))*g(1,2))*log(exp(-EulerGamma))+1/2*(1/4*muvsq^2*(-𝑖*𝜋^2*log(𝜋)+𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))*g(1,2)+3/8*𝑖*𝜋^2*muvsq^2*g(1,2)+1/4*𝑖*𝜋^2*muvsq^2*log(exp(-EulerGamma))*g(1,2))*log(muvsq^-1)^2+((1/4*muvsq^2*(-𝑖*𝜋^2*log(𝜋)+𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))*g(1,2)+3/8*𝑖*𝜋^2*muvsq^2*g(1,2))*log(exp(-EulerGamma))+𝑖*𝜋^2*(7/16*muvsq^2*g(1,2)+1/48*𝜋^2*muvsq^2*g(1,2))+3/8*muvsq^2*(-𝑖*𝜋^2*log(𝜋)+𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))*g(1,2)+1/4*muvsq^2*(1/2*𝑖*𝜋^2*log(𝜋)^2+1/2*𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1)^2-𝑖*𝜋^2*log(𝜋)*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))*g(1,2)+1/8*𝑖*𝜋^2*muvsq^2*log(exp(-EulerGamma))^2*g(1,2))*log(muvsq^-1)+𝑖*𝜋^2*(151190863202516241/410199796539607264*muvsq^2*g(1,2)+1/32*𝜋^2*muvsq^2*g(1,2))+3/8*muvsq^2*(1/2*𝑖*𝜋^2*log(𝜋)^2+1/2*𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1)^2-𝑖*𝜋^2*log(𝜋)*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1))*g(1,2)+1/4*muvsq^2*(-1/6*𝑖*𝜋^2*log(𝜋)^3+1/6*𝑖*𝜋^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1)^3+1/2*𝑖*𝜋^2*log(𝜋)^2*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1)-1/2*𝑖*𝜋^2*log(𝜋)*log(1/4*𝜋^-1*mursq*exp(-EulerGamma)^-1)^2)*g(1,2)+1/24*𝑖*𝜋^2*muvsq^2*log(muvsq^-1)^3*g(1,2)+1/24*𝑖*𝜋^2*muvsq^2*log(exp(-EulerGamma))^3*g(1,2))\
        )",
        )
        .unwrap(),
    );
    debug!("Evaluated integral: {}", evaluated_integral);

    let mut params = HashMap::default();
    params.insert("muvsq".into(), 1.0);
    params.insert("mursq".into(), 1.0);

    let numerical_partial_eval = Vakint::partial_numerical_evaluation(
        &vakint.settings,
        evaluated_integral.as_view(),
        &params,
    );
    debug!(
        "Partial eval: {}",
        numerical_partial_eval.to_canonical_string()
    );
    let numerical_partial_eval_canonical_str =
        numerical_partial_eval.to_canonical_string().replace(
            "(-12.0696723514860*g(1,2)*𝑖+0)*ε^2",
            "-12.0696723514860*g(1,2)*ε^2*𝑖",
        );
    assert_eq!(
        numerical_partial_eval_canonical_str,
        "-12.0696723514860*g(1,2)*ε^2*𝑖+-5.36845814123893*g(1,2)*𝑖+2.46740110027234*g(1,2)*ε^-1*𝑖+9.41170424471097*g(1,2)*ε*𝑖"
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

    params.insert("g(1,2)".into(), 1.0);
    let numerical_full_eval =
        Vakint::full_numerical_evaluation(&vakint.settings, evaluated_integral.as_view(), &params);
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
