use symbolica::try_parse;
use vakint::{Vakint, VakintExpression, VakintSettings};

fn main() {
    let vakint = Vakint::new(Some(VakintSettings {
        allow_unknown_integrals: false,
        ..VakintSettings::default()
    }))
    .unwrap();

    //println!("Supported topologies:\n{}", vakint.topologies);

    let input = try_parse!(
        "(vk::k(11,2)*vk::k(11,2)+vk::k(11,77)*vk::k(22,77)+vk::k(22,33)*p(42,33))*vk::topo(\
                    vk::prop(9,vk::edge(7,10),vk::k(11),mUVsqA,1)*\
                    vk::prop(33,vk::edge(7,10),vk::k(22),mUVsqA,2)*\
                    vk::prop(55,vk::edge(7,10),vk::k(11)+vk::k(22),mUVsqA,1)\
                )+\
                (2*vk::k(33,2)*vk::k(33,2)+17*vk::k(44,33)*p(17,33))*vk::topo(\
                    vk::prop(7,vk::edge(9,21),vk::k(33),mUVsqB,1)*\
                    vk::prop(13,vk::edge(9,21),vk::k(44),mUVsqB,2)*\
                    vk::prop(17,vk::edge(9,21),vk::k(33)+vk::k(44),mUVsqB,1)\
                )"
    )
    .unwrap();

    let output = vakint.to_canonical(input.as_view(), true).unwrap();

    println!(
        "\nInput:\n\n{}\n\nhas been matched to\n\n{}\n",
        VakintExpression::try_from(input).unwrap(),
        VakintExpression::try_from(output).unwrap()
    );
}
