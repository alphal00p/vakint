mod test_utils;
use log::debug;
use symbolica::atom::Atom;
use test_utils::{compare_output, get_vakint};
use vakint::{VakintError, VakintSettings};

#[test_log::test]
fn test_1l_matching() {
    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        ..VakintSettings::default()
    });

    //println!("Topologies:\n{}", vakint.topologies);

    compare_output(
        vakint
            .to_canonical(
                Atom::parse(
                    "(k(8,2)*k(8,2)+k(8,33)*p(12,33))*topo(\
                    prop(77,edge(42,42),k(8),muvsq,1)\
                )",
                )
                .unwrap()
                .as_view(),
                false,
            )
            .as_ref()
            .map(|a| a.as_view()),
        Atom::parse(
            "(k(1,2)^2+k(1,33)*p(12,33))*topo(\
                prop(1,edge(1,1),k(1),muvsq,1)\
            )",
        )
        .unwrap(),
    );

    compare_output(
        vakint
            .to_canonical(
                Atom::parse(
                    "(k(8,2)*k(8,2)+k(8,33)*p(12,33))*topo(\
                    prop(77,edge(42,42),k(8),muvsq,1)\
                )",
                )
                .unwrap()
                .as_view(),
                true,
            )
            .as_ref()
            .map(|a| a.as_view()),
        Atom::parse("(k(1,2)^2+k(1,33)*p(12,33))*topo(I1LA(muvsq,1))").unwrap(),
    );

    compare_output(
        vakint
            .to_canonical(
                Atom::parse("(k(1,2)^2+k(1,33)*p(12,33))*topo(I1LA(muvsq,-3))")
                    .unwrap()
                    .as_view(),
                false,
            )
            .as_ref()
            .map(|a| a.as_view()),
        Atom::parse("(k(1,2)^2+k(1,33)*p(12,33))*topo(prop(1,edge(1,1),k(1),muvsq,-3))").unwrap(),
    );
}

#[test]
fn test_2l_matching_3prop() {
    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        ..VakintSettings::default()
    });

    //println!("Topologies:\n{}", vakint.topologies);

    compare_output(
        vakint
            .to_canonical(
                Atom::parse(
                    "(k(11,2)*k(11,2)+k(11,77)*k(22,77)+k(22,33)*p(42,33))*topo(\
                            prop(9,edge(7,10),k(11),mUVsq,1)*\
                            prop(33,edge(7,10),k(22),mUVsq,2)*\
                            prop(55,edge(7,10),k(11)+k(22),mUVsq,1)\
                        )",
                )
                .unwrap()
                .as_view(),
                false,
            )
            .as_ref()
            .map(|a| a.as_view()),
        Atom::parse(
            "(k(1,2)*k(1,2)+k(1,77)*k(2,77)+k(2,33)*p(42,33))*topo(\
                        prop(1,edge(1,2),k(1),mUVsq,1)*\
                        prop(2,edge(1,2),k(2),mUVsq,2)*\
                        prop(3,edge(2,1),k(1)+k(2),mUVsq,1)\
                    )",
        )
        .unwrap(),
    );

    compare_output(
        vakint
            .to_canonical(
                Atom::parse(
                    "(k(11,2)*k(11,2)+k(11,77)*k(22,77)+k(22,33)*p(42,33))*topo(\
                            prop(9,edge(7,10),k(11),mUVsq,1)*\
                            prop(33,edge(7,10),k(22),mUVsq,2)*\
                            prop(55,edge(7,10),k(11)+k(22),mUVsq,1)\
                        )",
                )
                .unwrap()
                .as_view(),
                true,
            )
            .as_ref()
            .map(|a| a.as_view()),
        Atom::parse("(k(1,2)*k(1,2)+k(1,77)*k(2,77)+k(2,33)*p(42,33))*topo(I2LA(mUVsq,1,2,1))")
            .unwrap(),
    );

    compare_output(
        vakint.to_canonical(
            Atom::parse("(k(1,2)*k(1,2)+k(1,77)*k(2,77)+k(2,33)*p(42,33))*topo(I2LA(mUVsq,1,2,1))")
                .unwrap()
                .as_view(),
            false,
        ).as_ref().map(|a| a.as_view()),
        Atom::parse("(k(1,2)*k(1,2)+k(1,77)*k(2,77)+k(2,33)*p(42,33))*topo(\
                                                prop(1,edge(1,2),k(1),mUVsq,1)*prop(2,edge(1,2),k(2),mUVsq,2)*\
                                                prop(3,edge(2,1),k(1)+k(2),mUVsq,1)\
                                            )").unwrap(),
    );
}

#[test_log::test]
fn test_2l_matching_pinched() {
    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        ..VakintSettings::default()
    });

    debug!("Topologies:\n{}", vakint.topologies);
    /*
    let mut a = vakint
        .to_canonical(
            Atom::parse(
                "(k(11,2)*k(11,2)+k(11,77)*k(22,77)+k(22,33)*p(42,33))*topo(\
            prop(33,edge(7,10),k(22),mUVsq,2)*\
            prop(55,edge(7,10),k(11),mUVsq,1)\
        )",
            )
            .unwrap()
            .as_view(),
            false,
        )
        .unwrap();
    a = Atom::parse(
        AtomPrinter::new_with_options(a.as_view(), PrintOptions::file())
            .to_string()
            .as_str(),
    )
    .unwrap();
    let b = Atom::parse(
        "(-1*k(1,33)*p(42,33)+k(1,77)*k(2,77)+k(2,2)^2)*topo(prop(1,edge(2,2),k(1),mUVsq,2)*prop(2,edge(2,2),k(2),mUVsq,1))",
    )
    .unwrap();
    println!("BB\n{}", b.to_canonical_string());
    println!(
        "str(a)=str(b) -> {}",
        a.to_canonical_string() == b.to_canonical_string()
    );
    println!("a=b -> {}", a.as_view() == b.as_view());
    */

    compare_output(
        vakint
            .to_canonical(
                Atom::parse(
                    "(k(11,2)*k(11,2)+k(11,77)*k(22,77)+k(22,33)*p(42,33))*topo(\
                            prop(33,edge(7,10),k(22),mUVsq,2)*\
                            prop(55,edge(7,10),k(11),mUVsq,1)\
                        )",
                )
                .unwrap()
                .as_view(),
                false,
            )
            .as_ref()
            .map(|a| a.as_view()),
        Atom::parse(
            "(k(2,2)^2-k(1,33)*p(42,33)+k(1,77)*k(2,77))*topo(\
                        prop(1,edge(2,2),k(1),mUVsq,2)*\
                        prop(2,edge(2,2),k(2),mUVsq,1)\
                    )",
        )
        .unwrap(),
    );

    compare_output(
        vakint
            .to_canonical(
                Atom::parse(
                    "(k(11,2)*k(11,2)+k(11,77)*k(22,77)+k(22,33)*p(42,33))*topo(\
                            prop(33,edge(7,10),k(22),mUVsq,2)*\
                            prop(55,edge(7,10),k(11),mUVsq,1)\
                        )",
                )
                .unwrap()
                .as_view(),
                true,
            )
            .as_ref()
            .map(|a| a.as_view()),
        Atom::parse("(k(2,2)^2-k(1,33)*p(42,33)+k(1,77)*k(2,77))*topo(I2LA(mUVsq,2,1,0))").unwrap(),
    );

    compare_output(
        vakint
            .to_canonical(
                Atom::parse("(k(2,2)^2+k(1,33)*p(42,33)+k(1,77)*k(2,77))*topo(I2LA(mUVsq,2,1,0))")
                    .unwrap()
                    .as_view(),
                false,
            )
            .as_ref()
            .map(|a| a.as_view()),
        Atom::parse(
            "(k(2,2)^2-k(1,33)*p(42,33)+k(1,77)*k(2,77))*topo(\
                        prop(1,edge(2,2),k(1),mUVsq,2)*\
                        prop(2,edge(2,2),k(2),mUVsq,1)\
                )",
        )
        .unwrap(),
    );

    compare_output(
        vakint
            .to_canonical(
                Atom::parse("(k(2,2)^2+k(1,33)*p(42,33)+k(1,77)*k(2,77))*topo(I2LA(mUVsq,2,0,1))")
                    .unwrap()
                    .as_view(),
                false,
            )
            .as_ref()
            .map(|a| a.as_view()),
        Atom::parse(
            "(k(2,2)^2-k(1,33)*p(42,33)+k(1,77)*k(2,77))*topo(\
                        prop(1,edge(2,2),k(1),mUVsq,2)*\
                        prop(3,edge(2,2),k(1)+k(2),mUVsq,1)\
                )",
        )
        .unwrap(),
    );
}

#[test]
fn test_unknown_integrals() {
    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        ..VakintSettings::default()
    });
    let vakint_with_unknown_integrals = get_vakint(VakintSettings {
        allow_unknown_integrals: true,
        ..VakintSettings::default()
    });

    let unknown_integral = Atom::parse(
        "(k(11,2)*k(11,2)+k(11,77)*k(22,77)+k(22,33)*p(42,33))*topo(\
                    prop(33,edge(7,10),k(22),mA,2)*\
                    prop(55,edge(7,10),k(11),mB,1)\
                )",
    )
    .unwrap();

    assert!(matches!(
        vakint.to_canonical(unknown_integral.as_view(), false),
        Err(VakintError::UnreckognizedIntegral(_))
    ));

    // println!("Topologies:\n{}", vakint_with_unknown_integrals.topologies);

    compare_output(
        vakint_with_unknown_integrals
            .to_canonical(unknown_integral.as_view(), false)
            .as_ref()
            .map(|a| a.as_view()),
        Atom::parse(
            "(k(11,2)^2+k(11,77)*k(22,77)+k(22,33)*p(42,33))*topo(UNKNOWN(\
                        prop(33,edge(7,10),k(22),mA,2)*\
                        prop(55,edge(7,10),k(11),mB,1))\
                    )",
        )
        .unwrap(),
    );
}
