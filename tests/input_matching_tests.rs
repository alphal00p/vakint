mod test_utils;
use log::debug;
use symbolica::try_parse;
use test_utils::{compare_output, get_vakint};
use vakint::{VakintError, VakintSettings};

#[test_log::test]
fn test_1l_matching() {
    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        ..VakintSettings::default()
    });

    debug!("Topologies:\n{}", vakint.topologies);
    _ = compare_output(
        vakint
            .to_canonical(
                try_parse!(
                    "(lm(8,2)*lm(8,2)+lm(8,33)*pext(12,33))*vk::topo(\
                    vk::prop(77,vk::edge(42,42),lm(8),muvsq,1)\
                )"
                )
                .unwrap()
                .as_view(),
                false,
            )
            .as_ref()
            .map(|a| a.as_view()),
        try_parse!(
            "(vk::k(1,2)^2+vk::k(1,33)*pext(12,33))*vk::topo(\
                vk::prop(1,vk::edge(1,1),vk::k(1),muvsq,1)\
            )"
        )
        .unwrap(),
    );

    _ = compare_output(
        vakint
            .to_canonical(
                try_parse!(
                    "(lm(8,2)*lm(8,2)+lm(8,33)*pext(12,33))*vk::topo(\
                    vk::prop(77,vk::edge(42,42),lm(8),muvsq,1)\
                )"
                )
                .unwrap()
                .as_view(),
                true,
            )
            .as_ref()
            .map(|a| a.as_view()),
        try_parse!("(vk::k(1,2)^2+vk::k(1,33)*pext(12,33))*vk::topo(vk::I1L(muvsq,1))").unwrap(),
    );

    _ = compare_output(
        vakint
            .to_canonical(
                try_parse!("(lm(1,2)^2+lm(1,33)*pext(12,33))*vk::topo(vk::I1L(muvsq,-3))")
                    .unwrap()
                    .as_view(),
                false,
            )
            .as_ref()
            .map(|a| a.as_view()),
        try_parse!("(lm(1,2)^2+lm(1,33)*pext(12,33))*vk::topo(vk::prop(1,vk::edge(1,1),vk::k(1),muvsq,-3))")
            .unwrap(),
    );
}

#[test]
fn test_2l_matching_3prop() {
    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        ..VakintSettings::default()
    });

    //println!("Topologies:\n{}", vakint.topologies);

    _ = compare_output(
        vakint
            .to_canonical(
                try_parse!(
                    "(k(11,2)*k(11,2)+k(11,77)*k(22,77)+k(22,33)*p(42,33))*vk::topo(\
                            vk::prop(9,vk::edge(7,10),k(11),mUVsq,1)*\
                            vk::prop(33,vk::edge(7,10),k(22),mUVsq,2)*\
                            vk::prop(55,vk::edge(7,10),k(11)+k(22),mUVsq,1)\
                        )"
                )
                .unwrap()
                .as_view(),
                false,
            )
            .as_ref()
            .map(|a| a.as_view()),
        try_parse!(
            "(vk::k(1,2)*vk::k(1,2)+vk::k(1,77)*vk::k(2,77)+vk::k(2,33)*p(42,33))*vk::topo(\
                        vk::prop(1,vk::edge(1,2),vk::k(1),mUVsq,1)*\
                        vk::prop(2,vk::edge(1,2),vk::k(2),mUVsq,2)*\
                        vk::prop(3,vk::edge(2,1),vk::k(1)+vk::k(2),mUVsq,1)\
                    )"
        )
        .unwrap(),
    );

    _ = compare_output(
        vakint
            .to_canonical(
                try_parse!(
                    "(vk::k(11,2)*vk::k(11,2)+vk::k(11,77)*vk::k(22,77)+vk::k(22,33)*p(42,33))*vk::topo(\
                            vk::prop(9,vk::edge(7,10),vk::k(11),mUVsq,1)*\
                            vk::prop(33,vk::edge(7,10),vk::k(22),mUVsq,2)*\
                            vk::prop(55,vk::edge(7,10),vk::k(11)+vk::k(22),mUVsq,1)\
                        )"
                )
                .unwrap()
                .as_view(),
                true,
            )
            .as_ref()
            .map(|a| a.as_view()),
        try_parse!("(vk::k(1,2)*vk::k(1,2)+vk::k(1,77)*vk::k(2,77)+vk::k(2,33)*p(42,33))*vk::topo(vk::I2L(mUVsq,1,2,1))")
            .unwrap(),
    );

    _ = compare_output(
        vakint.to_canonical(
            try_parse!("(k(1,2)*k(1,2)+k(1,77)*k(2,77)+k(2,33)*p(42,33))*vk::topo(vk::I2L(mUVsq,1,2,1))")
                .unwrap()
                .as_view(),
            false,
        ).as_ref().map(|a| a.as_view()),
        try_parse!("(k(1,2)*k(1,2)+k(1,77)*k(2,77)+k(2,33)*p(42,33))*vk::topo(\
                                                vk::prop(1,vk::edge(1,2),vk::k(1),mUVsq,1)*vk::prop(2,vk::edge(1,2),vk::k(2),mUVsq,2)*\
                                                vk::prop(3,vk::edge(2,1),vk::k(1)+vk::k(2),mUVsq,1)\
                                            )").unwrap(),
    );
}

#[test_log::test]
fn test_2l_matching_pinched() {
    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        ..VakintSettings::default()
    });

    _ = compare_output(
        vakint
            .to_canonical(
                try_parse!(
                    "(k(11,2)*k(11,2)+k(11,77)*k(22,77)+k(22,33)*p(42,33))*vk::topo(\
                            vk::prop(33,vk::edge(10,10),k(22),mUVsq,2)*\
                            vk::prop(55,vk::edge(10,10),k(11),mUVsq,1)\
                        )"
                )
                .unwrap()
                .as_view(),
                false,
            )
            .as_ref()
            .map(|a| a.as_view()),
        try_parse!(
            "(vk::k(2,2)^2+vk::k(1,33)*p(42,33)+vk::k(1,77)*vk::k(2,77))*vk::topo(\
                        vk::prop(1,vk::edge(1,1),vk::k(1),mUVsq,2)*\
                        vk::prop(2,vk::edge(1,1),vk::k(2),mUVsq,1)
                    )"
        )
        .unwrap(),
    );

    _ = compare_output(
        vakint
            .to_canonical(
                try_parse!(
                    "(k(11,2)*k(11,2)+k(11,77)*k(22,77)+k(22,33)*p(42,33))*vk::topo(\
                            vk::prop(33,vk::edge(10,10),k(22),mUVsq,2)*\
                            vk::prop(55,vk::edge(10,10),k(11),mUVsq,1)\
                        )"
                )
                .unwrap()
                .as_view(),
                true,
            )
            .as_ref()
            .map(|a| a.as_view()),
        try_parse!("(vk::k(2,2)^2+vk::k(1,33)*p(42,33)+vk::k(1,77)*vk::k(2,77))*vk::topo(vk::I2L_pinch_3(mUVsq,2,1,0))")
            .unwrap(),
    );

    _ = compare_output(
        vakint
            .to_canonical(
                try_parse!(
                    "(vk::k(2,2)^2+vk::k(1,33)*p(42,33)+vk::k(1,77)*vk::k(2,77))*vk::topo(vk::I2L(mUVsq,2,1,0))"
                )
                .unwrap()
                .as_view(),
                false,
            )
            .as_ref()
            .map(|a| a.as_view()),
        try_parse!(
            "(vk::k(2,2)^2+vk::k(1,33)*p(42,33)+vk::k(1,77)*vk::k(2,77))*vk::topo(\
                        vk::prop(1,vk::edge(1,2),vk::k(1),mUVsq,2)*\
                        vk::prop(2,vk::edge(1,2),vk::k(2),mUVsq,1)*\
                        vk::prop(3,vk::edge(2,1),vk::k(1)+vk::k(2),mUVsq,0)
                )"
        )
        .unwrap(),
    );

    _ = compare_output(
        vakint
            .to_canonical(
                try_parse!(
                    "(vk::k(2,2)^2+vk::k(1,33)*p(42,33)+vk::k(1,77)*vk::k(2,77))*vk::topo(vk::I2L_pinch_3(mUVsq,2,1,0))"
                )
                .unwrap()
                .as_view(),
                false,
            )
            .as_ref()
            .map(|a| a.as_view()),
        try_parse!(
            "(vk::k(2,2)^2+vk::k(1,33)*p(42,33)+vk::k(1,77)*vk::k(2,77))*vk::topo(\
                        vk::prop(1,vk::edge(1,1),vk::k(1),mUVsq,2)*\
                        vk::prop(2,vk::edge(1,1),vk::k(2),mUVsq,1)
                )"
        )
        .unwrap(),
    );
}

#[test_log::test]
fn test_3l_matching_with_zero_powers_in_short_form() {
    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        ..VakintSettings::default()
    });

    debug!("Topologies:\n{}", vakint.topologies);

    _ = compare_output(
        vakint
            .to_canonical(
                try_parse!(
                    "( 1 )*vk::topo(\
                vk::prop(1,vk::edge(1,2),k(1),muvsq,1)\
              * vk::prop(2,vk::edge(1,2),k(2),muvsq,1)\
              * vk::prop(3,vk::edge(1,2),k(3),muvsq,1)\
              * vk::prop(4,vk::edge(2,1),k(1)+k(2)+k(3),muvsq,2)\
            )"
                )
                .unwrap()
                .as_view(),
                true,
            )
            .as_ref()
            .map(|a| a.as_view()),
        try_parse!("vk::topo(vk::I3L_pinch_1_6(muvsq,0,1,1,1,2,0))").unwrap(),
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

    let unknown_integral = try_parse!(
        "(vk::k(1,2)*vk::k(1,2)+vk::k(1,77)*vk::k(2,77)+vk::k(2,33)*p(42,33))*vk::topo(\
                    vk::prop(1,vk::edge(7,10),vk::k(2),mA,2)*\
                    vk::prop(2,vk::edge(7,10),vk::k(1),mB,1)\
                )"
    )
    .unwrap();

    assert!(matches!(
        vakint.to_canonical(unknown_integral.as_view(), false),
        Err(VakintError::UnreckognizedIntegral(_))
    ));

    // println!("Topologies:\n{}", vakint_with_unknown_integrals.topologies);

    _ = compare_output(
        vakint_with_unknown_integrals
            .to_canonical(unknown_integral.as_view(), false)
            .as_ref()
            .map(|a| a.as_view()),
        try_parse!(
            "(vk::k(1,2)^2+vk::k(1,77)*vk::k(2,77)+vk::k(2,33)*p(42,33))*vk::topo(vk::UNKNOWN(\
                        vk::prop(1,vk::edge(7,10),vk::k(2),mA,2)*\
                        vk::prop(2,vk::edge(7,10),vk::k(1),mB,1))\
                    )"
        )
        .unwrap(),
    );
}

#[test]
fn test_2l_pinched_matching() {
    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        ..VakintSettings::default()
    });

    // println!("Topologies:\n{}", vakint_with_unknown_integrals.topologies);

    _ = compare_output(
        vakint
            .to_canonical(
                try_parse!("vk::topo(vk::prop(1,vk::edge(1,1),k(1),muvsq,1)*vk::prop(2,vk::edge(1,1),k(2),muvsq,1))")
                    .unwrap()
                    .as_view(),
                true,
            )
            .as_ref()
            .map(|a| a.as_view()),
        try_parse!("vk::topo(vk::I2L_pinch_3(muvsq,1,1,0))").unwrap(),
    );
}
