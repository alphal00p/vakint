mod test_utils;
use symbolica::atom::Atom;
use test_utils::{compare_output, get_vakint};
use vakint::VakintSettings;

#[test_log::test]
fn test_integrate_1l_a() {
    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
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
    integral = vakint.tensor_reduce(integral.as_view(), true).unwrap();

    compare_output(
        vakint.evaluate_integral(integral.as_view()),
        Atom::parse(
            "(\
                -(2*ε-4)^-1*dot(k(1),k(1))*g(1,2)\
            )*topo(I1LA(muvsq,1))",
        )
        .unwrap(),
    );
}
