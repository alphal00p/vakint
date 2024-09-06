mod test_utils;

use symbolica::atom::Atom;
use test_utils::compare_analytical_vs_pysecdec;

const COMPARISON_WITH_PYSECDEC_REL_THRESHOLD: f64 = 1.0e-7;
// QMC pySecDec error is very optimistic
const MAX_PULL: f64 = 100.0;

#[test_log::test]
fn test_integrate_1l_pysecdec() {
    compare_analytical_vs_pysecdec(
        Atom::parse(
            "(1)*topo(\
            prop(1,edge(1,1),k(1),muvsq,1)\
        )",
        )
        .unwrap()
        .as_view(),
        [("muvsq".into(), 1.0), ("mursq".into(), 1.0)]
            .iter()
            .cloned()
            .collect(),
        (1..=1)
            .map(|i| (format!("p{}", i), (13.0, 4.0, 3.0, 12.0)))
            .collect(),
        COMPARISON_WITH_PYSECDEC_REL_THRESHOLD,
        MAX_PULL,
        true,
    );
}

#[test_log::test]
fn test_integrate_1l_pysecdec_non_unit_mass() {
    compare_analytical_vs_pysecdec(
        Atom::parse(
            "(1)*topo(\
            prop(1,edge(1,1),k(1),muvsq,1)\
        )",
        )
        .unwrap()
        .as_view(),
        [("muvsq".into(), 2.0), ("mursq".into(), 1.0)]
            .iter()
            .cloned()
            .collect(),
        (1..=1)
            .map(|i| (format!("p{}", i), (13.0, 4.0, 3.0, 12.0)))
            .collect(),
        COMPARISON_WITH_PYSECDEC_REL_THRESHOLD,
        MAX_PULL,
        true,
    );
}

#[test_log::test]
fn test_integrate_1l_pysecdec_non_unit_scale() {
    compare_analytical_vs_pysecdec(
        Atom::parse(
            "(1)*topo(\
            prop(1,edge(1,1),k(1),muvsq,1)\
        )",
        )
        .unwrap()
        .as_view(),
        [("muvsq".into(), 1.0), ("mursq".into(), 2.0)]
            .iter()
            .cloned()
            .collect(),
        (1..=1)
            .map(|i| (format!("p{}", i), (13.0, 4.0, 3.0, 12.0)))
            .collect(),
        COMPARISON_WITH_PYSECDEC_REL_THRESHOLD,
        MAX_PULL,
        true,
    );
}

#[test_log::test]
fn test_integrate_1l_pysecdec_num_rank_two() {
    compare_analytical_vs_pysecdec(
        Atom::parse(
            "((k(1,33)*k(1,33))^2+k(1,55)*p(1,55))*topo(\
            prop(1,edge(1,1),k(1),muvsq,1)\
        )",
        )
        .unwrap()
        .as_view(),
        [("muvsq".into(), 1.0), ("mursq".into(), 1.0)]
            .iter()
            .cloned()
            .collect(),
        (1..=1)
            .map(|i| (format!("p{}", i), (13.0, 4.0, 3.0, 12.0)))
            .collect(),
        COMPARISON_WITH_PYSECDEC_REL_THRESHOLD,
        MAX_PULL,
        true,
    );
}
