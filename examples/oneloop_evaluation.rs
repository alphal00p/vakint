use ahash::HashMap;
use vakint::{Vakint, VakintExpression, VakintSettings};

fn main() {
    let vakint = Vakint::new(Some(VakintSettings {
        allow_unknown_integrals: false,
        use_dot_product_notation: true,
        integral_normalization_factor: vakint::LoopNormalizationFactor::MSbar,
        run_time_decimal_precision: 16,
        ..VakintSettings::default()
    }))
    .unwrap();

    // let mut integral = symbolica::try_parse!(
    //     "(
    //           vakint::k(3,mink4(4,11))*vakint::k(3,mink4(4,22))
    //         + vakint::k(3,mink4(4,77))*vakint::p(8,mink4(4,77))
    //         + vakint::p(1,mink4(4,77))*vakint::p(2,mink4(4,77))
    //         + vakint::p(1,mink4(4,99))*vakint::p(2,mink4(4,101))
    //     )*vakint::topo(\
    //     vakint::prop(9,vakint::edge(66,66),vakint::k(3),MUVsq,1)\
    // )"
    // )
    // .unwrap();

    // let mut integral = symbolica::try_parse!(
    //     "(
    //           vakint::k(3,mink4(4,11))*vakint::k(3,mink4(4,22))
    //         + vakint::p(1,mink4(4,11))*vakint::p(2,mink4(4,11))
    //     )*vakint::topo(\
    //     vakint::prop(9,vakint::edge(66,66),vakint::k(3),MUVsq,1)\
    // )"
    // )
    // .unwrap();

    let mut integral = symbolica::try_parse!(
        "(
             ((vakint::k(3,mink4(4,11))+vakint::p(1,mink4(4,11)))*(vakint::k(3,mink4(4,22))+vakint::p(2,mink4(4,22))))^2
        )*vakint::topo(\
        vakint::prop(9,vakint::edge(66,66),vakint::k(3),MUVsq,1)\
    )"
    )
    .unwrap();
    /*
        let mut integral = symbolica::try_parse!(
            "(
                  vakint::k(3,mink4(4,11))*vakint::k(3,mink4(4,11))
                + vakint::p(1,mink4(4,12))*vakint::p(2,mink4(4,12))
            )*vakint::topo(\
            vakint::prop(9,vakint::edge(66,66),vakint::k(3),MUVsq,1)\
        )"
        )
        .unwrap();
    */

    println!(
        "\nInput integral:\n{}\n",
        VakintExpression::try_from(integral.clone()).unwrap()
    );

    integral = vakint.to_canonical(integral.as_view(), true).unwrap();
    println!(
        "Matched integral:\n{}\n",
        VakintExpression::try_from(integral.clone()).unwrap()
    );

    integral = vakint.tensor_reduce(integral.as_view()).unwrap();
    println!(
        "Tensor reduced integral:\n{}\n",
        VakintExpression::try_from(integral.clone()).unwrap()
    );

    integral = vakint.evaluate_integral(integral.as_view()).unwrap();
    println!("Evaluated integral:\n{}\n", integral);

    let mut params = HashMap::default();
    params.insert(
        "oneloop_evaluation::MUVsq".into(),
        vakint.settings.real_to_prec("1.0"),
    );
    params.insert("vakint::mursq".into(), vakint.settings.real_to_prec("1.0"));

    let numerical_partial_eval = Vakint::partial_numerical_evaluation(
        &vakint.settings,
        integral.as_view(),
        &params,
        &HashMap::default(),
        None,
    );
    println!("Partial eval:\n{}\n", numerical_partial_eval);

    let externals = vakint.externals_from_f64(
        &(1..=2)
            .map(|i| {
                (
                    i,
                    (
                        0.17 * ((i + 1) as f64),
                        0.4 * ((i + 2) as f64),
                        0.3 * ((i + 3) as f64),
                        0.12 * ((i + 4) as f64),
                    ),
                )
            })
            .collect(),
    );

    params.insert(
        "vakint::g(oneloop_evaluation::mink4(4,22),oneloop_evaluation::mink4(4,11))".into(),
        vakint.settings.real_to_prec("1.0"),
    );
    params.insert(
        "vakint::p(1,oneloop_evaluation::mink4(4,99))".into(),
        vakint.settings.real_to_prec("1.0"),
    );
    params.insert(
        "vakint::p(2,oneloop_evaluation::mink4(4,101))".into(),
        vakint.settings.real_to_prec("1.0"),
    );
    let numerical_full_eval = Vakint::full_numerical_evaluation_without_error(
        &vakint.settings,
        integral.as_view(),
        &params,
        &HashMap::default(),
        Some(&externals),
    )
    .unwrap();
    println!(
        "Full eval (tensor structures left all substituted with 1):\n{}\n",
        numerical_full_eval
    );
}
