use ahash::HashMap;
use symbolica::atom::Atom;
use vakint::{Vakint, VakintExpression, VakintSettings};

fn main() {
    let vakint = Vakint::new(Some(VakintSettings {
        allow_unknown_integrals: false,
        use_dot_product_notation: true,
        integral_normalization: "1".into(),
        n_digits_at_evaluation_time: 16,
        ..VakintSettings::default()
    }))
    .unwrap();

    let mut integral = Atom::parse(
        "(k(3,11)*k(3,22)+k(3,77)*p(8,77))*topo(\
        prop(9,edge(66,66),k(3),MUVsq,1)\
    )",
    )
    .unwrap();
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
    params.insert("MUVsq".into(), 1.0);
    params.insert("mursq".into(), 1.0);

    let numerical_partial_eval = vakint.partial_numerical_evaluation(integral.as_view(), &params);
    println!("Partial eval:\n{}\n", numerical_partial_eval);

    params.insert("g(11,22)".into(), 1.0);
    let numerical_full_eval = vakint
        .full_numerical_evaluation(integral.as_view(), &params)
        .unwrap();
    println!(
        "Full eval (metric substituted with 1):\n{}\n",
        numerical_full_eval
    );
}
