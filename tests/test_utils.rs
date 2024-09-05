use std::f64::consts::LOG2_10;

use colored::Colorize;
use log::debug;
use std::vec;
use symbolica::{
    atom::{Atom, AtomView},
    domains::float::{Complex, Float},
    printer::{AtomPrinter, PrintOptions},
};

use std::collections::HashMap;
use vakint::{EvaluationMethod, LoopNormalizationFactor, PySecDecOptions, VakintSettings};
use vakint::{NumericalEvaluationResult, Vakint, VakintError};

pub fn get_vakint(vakint_settings: VakintSettings) -> Vakint {
    match Vakint::new(Some(vakint_settings)) {
        Ok(r) => r,
        Err(err) => panic!("Failed to initialize vakint: {}", err),
    }
}

#[allow(unused)]
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
            if r != expected_output.as_view() {
                println!(
                    "Output does not match expected output:\n{}\n!=\n{}",
                    format!("{}", r_processed_view).red(),
                    format!("{}", expected_output).green()
                );
            }
            assert_eq!(r, expected_output.as_view());
            r.to_owned()
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

#[allow(unused)]
pub fn compare_analytical_vs_pysecdec(
    integra_view: AtomView,
    numerical_masses: HashMap<String, f64, ahash::RandomState>,
    numerical_external_momenta: HashMap<String, (f64, f64, f64, f64), ahash::RandomState>,
    rel_threshold: f64,
    max_pull: f64,
    quiet: bool,
) {
    // First perform the analytic evaluation

    let mut vakint_analytic_settings = VakintSettings {
        allow_unknown_integrals: false,
        use_dot_product_notation: true,
        mu_r_sq_symbol: "mursq".into(),
        integral_normalization_factor: LoopNormalizationFactor::MSbar,
        n_digits_at_evaluation_time: 16,
        ..VakintSettings::default()
    };
    vakint_analytic_settings.evaluation_order = vakint_analytic_settings
        .evaluation_order
        .iter()
        .filter(|&method| !matches!(method, EvaluationMethod::PySecDec(_)))
        .cloned()
        .collect::<Vec<_>>();
    let mut vakint = get_vakint(vakint_analytic_settings);

    let mut eval_params = HashMap::default();
    eval_params.insert(
        "muvsq".into(),
        *numerical_masses
            .get("muvsq")
            .unwrap_or_else(|| panic!("muvsq not found in numerical_masses")),
    );
    eval_params.insert(
        "mursq".into(),
        *numerical_masses
            .get("mursq")
            .unwrap_or_else(|| panic!("mursq not found in numerical_masses")),
    );

    let integral = vakint.to_canonical(integra_view, true).unwrap();

    let integral_reduced = vakint.tensor_reduce(integral.as_view()).unwrap();
    let benchmark_evaluated_integral = vakint
        .evaluate_integral(integral_reduced.as_view())
        .unwrap();

    let benchmark = Vakint::full_numerical_evaluation(
        &vakint.settings,
        benchmark_evaluated_integral.as_view(),
        &eval_params,
    )
    .unwrap();

    // Now perform the pySecDec evaluation
    vakint.settings.evaluation_order = vec![EvaluationMethod::PySecDec(PySecDecOptions {
        quiet,
        relative_precision: rel_threshold * 1.0e-2,
        numerical_masses,
        numerical_external_momenta,
    })];

    let pysec_dec_eval = match vakint.evaluate_integral(integral.as_view()) {
        Ok(eval) => eval,
        Err(e) => {
            panic!("Error during pySecDec evaluation: {}", e);
        }
    };
    let (pysecdec_central, pysecdec_error) = match Vakint::full_numerical_evaluation_with_error(
        &vakint.settings,
        pysec_dec_eval.as_view(),
        &eval_params,
    ) {
        Ok(eval) => eval,
        Err(e) => {
            panic!("Error during parsing of numerical pySecDec result: {}", e);
        }
    };

    let (matches, msg) = benchmark.does_approx_match(
        &pysecdec_central,
        Some(&pysecdec_error),
        rel_threshold,
        max_pull,
    );
    if !matches || !quiet {
        println!("Benchmark        :\n{}", benchmark);
        println!("PySecDec central :\n{}", pysecdec_central);
        println!("PySecDec error   :\n{}", pysecdec_error);
        println!("{}", msg)
    } else {
        debug!("Benchmark        :\n{}", benchmark);
        debug!("PySecDec central :\n{}", pysecdec_central);
        debug!("PySecDec error   :\n{}", pysecdec_error);
        debug!("{}", msg)
    }
    assert!(matches, "Benchmark and numerical result do not match");
}
