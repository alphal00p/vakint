mod test_utils;

use ahash::HashMap;
use log::debug;
use symbolica::atom::Atom;
use test_utils::{compare_numerical_output, compare_output, get_vakint};
use vakint::VakintSettings;

#[test_log::test]
fn test_integrate_1l_a() {
    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        use_dot_product_notation: true,
        integral_normalization: "1".into(),
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

    let evaluated_integral = compare_output(
        evaluated_integral_res_ref
        .map(|a| a.as_view()),
        Atom::parse("\
            + ε^-1 * 1/64*𝑖*𝜋^-2*muvsq^2*g(1,2)\
            + ε^-0 * -1/128*𝑖*𝜋^-2*muvsq^2*log(muvsq*mursq^-1)*g(1,2)\
            + ε    * (1/768*𝑖*muvsq^2*g(1,2)+7/256*𝑖*𝜋^-2*muvsq^2*g(1,2)+1/512*𝑖*𝜋^-2*muvsq^2*log(muvsq*mursq^-1)^2*g(1,2)-3/256*𝑖*𝜋^-2*muvsq^2*log(muvsq*mursq^-1)*g(1,2))\
            + ε^2  * (1/512*𝑖*muvsq^2*g(1,2)+151190863202516241/6563196744633716224*𝑖*𝜋^-2*muvsq^2*g(1,2)-1/1536*𝑖*muvsq^2*log(muvsq*mursq^-1)*g(1,2)+3/1024*𝑖*𝜋^-2*muvsq^2*log(muvsq*mursq^-1)^2*g(1,2)-1/3072*𝑖*𝜋^-2*muvsq^2*log(muvsq*mursq^-1)^3*g(1,2)-7/512*𝑖*𝜋^-2*muvsq^2*log(muvsq*mursq^-1)*g(1,2))+3/128*𝑖*𝜋^-2*muvsq^2*g(1,2)"
        )
        .unwrap(),
    );
    debug!("Evaluated integral: {}", evaluated_integral);

    let mut params = HashMap::default();
    params.insert("muvsq".into(), 1.0);
    params.insert("mursq".into(), 1.0);

    let numerical_partial_eval =
        vakint.partial_numerical_evaluation(evaluated_integral.as_view(), &params);
    debug!(
        "Partial eval: {}",
        numerical_partial_eval.to_canonical_string()
    );
    assert_eq!(
        numerical_partial_eval.to_canonical_string(),
        "1.58314349441153e-3*g(1,2)*ε^-1*𝑖+2.37471524161729e-3*g(1,2)*𝑖+4.07258444855351e-3*g(1,2)*ε*𝑖+4.28717619663842e-3*g(1,2)*ε^2*𝑖"
    );
    /*
    compare_output(
        Ok(numerical_partial_eval.as_view()),
        Atom::parse(
            format!(
                "2.37471524161729e-3`{prec}*𝑖*g(1,2)\
                +4.07258444855351e-3`{prec}*𝑖*ε*g(1,2)\
                +1.58314349441153e-3`{prec}*𝑖*ε^-1*g(1,2)\
                +4.28717619663842e-3`{prec}*𝑖*ε^2*g(1,2)",
                prec = vakint.settings.n_digits_at_evaluation_time - 1
            )
            .as_str(),
        )
        .unwrap(),
    );
    */

    params.insert("g(1,2)".into(), 1.0);
    let numerical_full_eval =
        vakint.full_numerical_evaluation(evaluated_integral.as_view(), &params);
    let numerical_full_eval_ref = numerical_full_eval.as_ref();
    debug!(
        "Full eval (metric substituted with 1):\n{}",
        numerical_full_eval_ref.unwrap()
    );
    compare_numerical_output(
        numerical_full_eval_ref,
        vec![
            (-1, ("0.0".into(), "1.58314349441153e-3".into())),
            (0, ("0.0".into(), "2.37471524161729e-3".into())),
            (1, ("0.0".into(), "4.07258444855351e-3".into())),
            (2, ("0.0".into(), "4.28717619663842e-3".into())),
        ],
        vakint.settings.n_digits_at_evaluation_time,
    );
}
