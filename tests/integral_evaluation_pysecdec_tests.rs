mod test_utils;
use vakint::{EvaluationMethod, EvaluationOrder, PySecDecOptions};

use std::vec;

use symbolica::atom::Atom;
use test_utils::convert_test_params;

use crate::test_utils::{compare_vakint_evaluation_vs_reference, convert_test_externals};

const N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS: u32 = 10;
// PySecDec QMC is often very optimistic
const MAX_PULL: f64 = 1.0e5;

#[test_log::test]
fn test_integrate_2l_different_masses() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        EvaluationOrder(vec![EvaluationMethod::PySecDec(PySecDecOptions::default())]),
        Atom::parse(
            "(1)*topo(\
                prop(1,edge(1,2),k(1),muvsqA,1)\
                *prop(2,edge(1,2),k(2),muvsqB,1)\
                *prop(3,edge(2,1),k(1)+k(2),muvsqC,1)\
            )",
        )
        .unwrap()
        .as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        convert_test_params(&[("muvsqA".into(), 1.0), ("muvsqB".into(), 1.0), ("muvsqC".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS),
        convert_test_externals(
        &(1..=1)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS),
        vec![
            (-2, ("-146.1136365510036558546604990331".into(), "0.0".into()),),
            (-1, ("635.8146971740286947808753047759".into(),  "0.0".into()),),
            (0,  ("-1646.531034471454483109201793220".into(), "0.0".into()),),
            (1,  ("2240.516116133454318298096346441".into(),  "0.0".into()),),
        ],
        N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS, MAX_PULL
    );
}
