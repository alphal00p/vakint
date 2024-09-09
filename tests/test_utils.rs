use std::f64::consts::LOG2_10;

use colored::Colorize;
use log::{debug, info};
use std::vec;
use symbolica::{
    atom::{Atom, AtomView},
    domains::float::{Complex, Float, Real, RealNumberLike},
    printer::{AtomPrinter, PrintOptions},
};

use std::collections::HashMap;
use vakint::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, Momentum, PySecDecOptions,
    VakintSettings,
};
use vakint::{NumericalEvaluationResult, Vakint, VakintError};

pub fn get_vakint(vakint_settings: VakintSettings) -> Vakint {
    match Vakint::new(Some(vakint_settings)) {
        Ok(r) => r,
        Err(err) => panic!("Failed to initialize vakint: {}", err),
    }
}

#[allow(unused)]
pub fn convert_test_params(
    params: &HashMap<String, f64, ahash::RandomState>,
    decimal_prec: u32,
) -> HashMap<String, Complex<Float>, ahash::RandomState> {
    let binary_prec: u32 = ((decimal_prec as f64) * LOG2_10).floor() as u32;
    params
        .iter()
        .map(|(k, v)| {
            (
                k.clone(),
                Complex::new(
                    Float::with_val(binary_prec, v),
                    Float::with_val(binary_prec, 0),
                ),
            )
        })
        .collect()
}

#[allow(unused)]
pub fn convert_test_externals(
    externals: &HashMap<usize, (f64, f64, f64, f64), ahash::RandomState>,
    decimal_prec: u32,
) -> HashMap<usize, Momentum, ahash::RandomState> {
    let binary_prec: u32 = ((decimal_prec as f64) * LOG2_10).floor() as u32;
    externals
        .iter()
        .map(|(&k, v)| {
            (
                k,
                (
                    Complex::new(
                        Float::with_val(binary_prec, v.0),
                        Float::with_val(binary_prec, 0.0),
                    ),
                    Complex::new(
                        Float::with_val(binary_prec, v.1),
                        Float::with_val(binary_prec, 0.0),
                    ),
                    Complex::new(
                        Float::with_val(binary_prec, v.2),
                        Float::with_val(binary_prec, 0.0),
                    ),
                    Complex::new(
                        Float::with_val(binary_prec, v.3),
                        Float::with_val(binary_prec, 0.0),
                    ),
                ),
            )
        })
        .collect()
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
    numerical_masses: HashMap<String, Complex<Float>, ahash::RandomState>,
    numerical_external_momenta: HashMap<usize, Momentum, ahash::RandomState>,
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
        evaluation_order: EvaluationOrder::analytic_only(),
        ..VakintSettings::default()
    };
    let mut vakint = get_vakint(vakint_analytic_settings);

    let mut eval_params = HashMap::default();
    eval_params.insert(
        "muvsq".into(),
        numerical_masses
            .get("muvsq")
            .unwrap_or_else(|| panic!("muvsq not found in numerical_masses"))
            .to_owned(),
    );
    eval_params.insert(
        "mursq".into(),
        numerical_masses
            .get("mursq")
            .unwrap_or_else(|| panic!("mursq not found in numerical_masses"))
            .to_owned(),
    );

    let integral = vakint.to_canonical(integra_view, true).unwrap();

    let integral_reduced = vakint.tensor_reduce(integral.as_view()).unwrap();
    let benchmark_evaluated_integral = vakint
        .evaluate_integral(integral_reduced.as_view())
        .unwrap();

    let benchmark = Vakint::full_numerical_evaluation_without_error(
        &vakint.settings,
        benchmark_evaluated_integral.as_view(),
        &eval_params,
        Some(&numerical_external_momenta),
    )
    .unwrap();

    // Now perform the pySecDec evaluation
    let f64_numerical_masses = numerical_masses
        .iter()
        .map(|(k, v)| (k.clone(), v.norm().re.to_f64()))
        .collect();
    let f64_numerical_external_momenta = numerical_external_momenta
        .iter()
        .map(|(k, v)| {
            (
                format!("p{}", k),
                (
                    v.0.norm().re.to_f64(),
                    v.1.norm().re.to_f64(),
                    v.2.norm().re.to_f64(),
                    v.3.norm().re.to_f64(),
                ),
            )
        })
        .collect();
    vakint.settings.evaluation_order =
        EvaluationOrder(vec![EvaluationMethod::PySecDec(PySecDecOptions {
            quiet,
            relative_precision: rel_threshold * 1.0e-2,
            numerical_masses: f64_numerical_masses,
            numerical_external_momenta: f64_numerical_external_momenta,
        })]);

    let pysec_dec_eval = match vakint.evaluate_integral(integral.as_view()) {
        Ok(eval) => eval,
        Err(e) => {
            panic!("Error during pySecDec evaluation: {}", e);
        }
    };
    let (pysecdec_central, pysecdec_error) = match Vakint::full_numerical_evaluation(
        &vakint.settings,
        pysec_dec_eval.as_view(),
        &eval_params,
        Some(&numerical_external_momenta),
    ) {
        Ok(eval) => eval,
        Err(e) => {
            panic!("Error during parsing of numerical pySecDec result: {}", e);
        }
    };

    let (matches, msg) = benchmark.does_approx_match(
        &pysecdec_central,
        pysecdec_error.as_ref(),
        rel_threshold,
        max_pull,
    );
    if !matches || !quiet {
        println!("Analytic         :\n{}", benchmark);
        println!("PySecDec central :\n{}", pysecdec_central);
        println!("PySecDec error   :\n{}", pysecdec_error.unwrap());
        println!("{}", msg)
    } else {
        debug!("Analytic         :\n{}", benchmark);
        debug!("PySecDec central :\n{}", pysecdec_central);
        debug!("PySecDec error   :\n{}", pysecdec_error.unwrap());
        debug!("{}", msg)
    }
    assert!(matches, "Benchmark and numerical result do not match.");
}

#[allow(unused)]
pub fn compare_vakint_evaluation_vs_reference(
    evaluation_order: EvaluationOrder,
    integra_view: AtomView,
    numerical_masses: HashMap<String, Complex<Float>, ahash::RandomState>,
    numerical_external_momenta: HashMap<usize, Momentum, ahash::RandomState>,
    expected_output: Vec<(i64, (String, String))>,
    prec: u32,
    max_pull: f64,
) {
    // First perform Vakint evaluation

    let mut vakint_settings = VakintSettings {
        allow_unknown_integrals: true,
        use_dot_product_notation: true,
        mu_r_sq_symbol: "mursq".into(),
        integral_normalization_factor: LoopNormalizationFactor::MSbar,
        n_digits_at_evaluation_time: prec,
        evaluation_order,
        ..VakintSettings::default()
    };

    // Insert parameters in pySecDec evaluation too if present
    for eo in vakint_settings.evaluation_order.0.iter_mut() {
        if let EvaluationMethod::PySecDec(pysecdec_options) = eo {
            let f64_numerical_masses = numerical_masses
                .iter()
                .map(|(k, v)| (k.clone(), v.norm().re.to_f64()))
                .collect();
            let f64_numerical_external_momenta = numerical_external_momenta
                .iter()
                .map(|(k, v)| {
                    (
                        format!("p{}", k),
                        (
                            v.0.norm().re.to_f64(),
                            v.1.norm().re.to_f64(),
                            v.2.norm().re.to_f64(),
                            v.3.norm().re.to_f64(),
                        ),
                    )
                })
                .collect();
            pysecdec_options.relative_precision = 10.0_f64.powi(-(prec as i32));
            pysecdec_options.numerical_masses = f64_numerical_masses;
            pysecdec_options.numerical_external_momenta = f64_numerical_external_momenta;
        }
    }

    let mut vakint = get_vakint(vakint_settings);

    let integral = vakint.to_canonical(integra_view, true).unwrap();

    let integral = vakint.tensor_reduce(integral.as_view()).unwrap();

    let integral = vakint
        .evaluate_integral(integral.as_view())
        .unwrap_or_else(|op| panic!("Failed to evaluate integral: {}", op));

    let (result, error) = Vakint::full_numerical_evaluation(
        &vakint.settings,
        integral.as_view(),
        &numerical_masses,
        Some(&numerical_external_momenta),
    )
    .unwrap();

    let binary_prec: u32 = ((prec.max(16) as f64) * LOG2_10).floor() as u32;

    let reference = NumericalEvaluationResult(
        expected_output
            .iter()
            .map(|(eps_pwr, (re, im))| {
                (
                    *eps_pwr,
                    Complex::new(
                        Float::parse(re.as_str(), Some(binary_prec)).unwrap(),
                        Float::parse(im.as_str(), Some(binary_prec)).unwrap(),
                    ),
                )
            })
            .collect::<Vec<_>>(),
    );
    let (matches, msg) = result.does_approx_match(
        &reference,
        error.as_ref(),
        0.1_f64.powi((prec as i32) - 2),
        max_pull,
    );
    if !matches {
        info!("Result     :\n{}", result);
        if error.is_some() {
            println!("Error      :\n{}", error.unwrap());
        }
        info!("Reference  :\n{}", reference);
        info!("{}", msg)
    } else {
        debug!("Result        :\n{}", result);
        if error.is_some() {
            println!("Error      :\n{}", error.unwrap());
        }
        debug!("Reference     :\n{}", reference);
        debug!("{}", msg)
    }
    assert!(matches, "Vakint result and reference result do not match");
}
