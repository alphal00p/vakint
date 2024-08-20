use std::f64::consts::LOG2_10;

use colored::Colorize;
use symbolica::{
    atom::{Atom, AtomView},
    domains::float::{Complex, Float},
    printer::{AtomPrinter, PrintOptions},
};
use vakint::{NumericalEvaluationResult, Vakint, VakintError, VakintSettings};

pub fn get_vakint(vakint_settings: VakintSettings) -> Vakint {
    match Vakint::new(Some(vakint_settings)) {
        Ok(r) => r,
        Err(err) => panic!("Failed to initialize vakint: {}", err),
    }
}

pub fn compare_output(output: Result<AtomView, &VakintError>, expected_output: Atom) -> Atom {
    match output {
        Ok(r) => {
            let r_processed = Atom::parse(
                AtomPrinter::new_with_options(r, PrintOptions::file())
                    .to_string()
                    .as_str(),
            )
            .unwrap();
            let r_processed_view = r_processed.as_view();
            if r_processed_view != expected_output.as_view() {
                println!(
                    "Output does not match expected output:\n{}\n!=\n{}",
                    format!("{}", r_processed_view).red(),
                    format!("{}", expected_output).green()
                );
            }
            assert_eq!(r_processed_view, expected_output.as_view());
            r_processed_view.to_owned()
        }
        Err(err) => panic!("Error: {}", err),
    }
}

#[allow(unused)]
pub fn compare_numerical_output(
    output: Result<&NumericalEvaluationResult, &VakintError>,
    expected_output: Vec<(i64, (String, String))>,
    prec: u32,
) {
    let binary_prec: u32 = ((prec as f64) * LOG2_10).floor() as u32;

    match output {
        Ok(numerical_result) => {
            let r = numerical_result.get_epsilon_coefficients();
            if r.len() != expected_output.len() {
                println!(
                    "Output does not match expected output: length mismatch: {} != {}",
                    r.len(),
                    expected_output.len()
                );
            }
            assert!(r.len() == expected_output.len());
            for ((o_pwr, o), (e_pwr, (e_real, e_cmplx))) in r.iter().zip(expected_output) {
                if *o_pwr != e_pwr {
                    println!("Power mismatch: {} != {}", o_pwr, e_pwr);
                    assert_eq!(*o_pwr, e_pwr);
                }
                let trgt = Complex::new(
                    Float::parse(e_real.as_str(), Some(binary_prec)).unwrap(),
                    Float::parse(e_cmplx.as_str(), Some(binary_prec)).unwrap(),
                );
                let scale = if trgt.norm_squared() > Float::with_val(binary_prec, 0.0) {
                    trgt.norm_squared()
                } else {
                    Float::with_val(binary_prec, 0.0)
                };
                let o_prec = (o.clone() - trgt.clone()).norm_squared() / scale;
                let trgt_prec = Float::with_val(binary_prec, (10.0_f64).powi(-(prec as i32)));
                if o_prec > trgt_prec {
                    println!(
                        "Output does not match expected output:\n{}\n!=\n{} (error: {} > target precision {})",
                        format!("{}", o).red(),
                        format!("{}", trgt).green(),
                        o_prec, trgt_prec
                    );
                }
                assert!(o_prec < trgt_prec);
            }
        }
        Err(err) => panic!("Error: {}", err),
    }
}
