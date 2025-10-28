use std::f64::consts::LOG2_10;

use colored::Colorize;
use log::{debug, info};
use symbolica::{
    atom::{Atom, AtomView},
    domains::float::{Complex, Float},
    printer::{AtomPrinter, PrintOptions},
    try_parse,
};

use std::collections::HashMap;
use vakint::{EvaluationOrder, LoopNormalizationFactor, Momentum, VakintSettings};
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
            let r_processed = try_parse!(
                AtomPrinter::new_with_options(
                    r,
                    PrintOptions {
                        hide_namespace: Some("tests"),
                        ..PrintOptions::file()
                    }
                )
                .to_string()
                .as_str()
            )
            .unwrap();
            let r_processed_view = r_processed.as_view();
            if r != expected_output.as_view() {
                println!(
                    "Output does not match expected output:\n{}\n!=\n{}",
                    format!(
                        "{}",
                        AtomPrinter::new_with_options(r_processed_view, PrintOptions::file())
                    )
                    .red(),
                    format!(
                        "{}",
                        AtomPrinter::new_with_options(
                            expected_output.as_view(),
                            PrintOptions::file()
                        )
                    )
                    .green()
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
                        o_prec,
                        trgt_prec
                    );
                }
                assert!(o_prec < trgt_prec);
            }
        }
        Err(err) => panic!("Error: {}", err),
    }
}

// EvaluationOrder::analytic_only()

#[allow(unused, clippy::too_many_arguments)]
pub fn compare_two_evaluations(
    vakint_default_settings: VakintSettings,
    evaluation_orders: ((&EvaluationOrder, bool), (&EvaluationOrder, bool)),
    integra_view: AtomView,
    numerical_masses: HashMap<String, Float, ahash::RandomState>,
    numerical_external_momenta: HashMap<usize, Momentum, ahash::RandomState>,
    rel_threshold: f64,
    max_pull: f64,
    quiet: bool,
) {
    let mut mod_evaluation_order_a = evaluation_orders.0.0.clone();
    let mut mod_evaluation_order_b = evaluation_orders.1.0.clone();
    for eval_order in [&mut mod_evaluation_order_a, &mut mod_evaluation_order_b] {
        eval_order.adjust(
            Some(quiet),
            rel_threshold * 1.0e-2,
            &numerical_masses,
            &HashMap::default(),
            &numerical_external_momenta,
        );
    }

    // First perform the first evaluation type

    let mut vakint_analytic_settings = VakintSettings {
        allow_unknown_integrals: false,
        use_dot_product_notation: true,
        mu_r_sq_symbol: "mursq".into(),
        integral_normalization_factor: LoopNormalizationFactor::pySecDec,
        evaluation_order: mod_evaluation_order_a.clone(),
        ..vakint_default_settings
    };
    let mut vakint = get_vakint(vakint_analytic_settings);

    let mut eval_params = HashMap::default();
    if numerical_masses.contains_key("user_space::muv") {
        eval_params.insert(
            "user_space::muv".into(),
            numerical_masses
                .get("user_space::muv")
                .unwrap_or_else(|| panic!("user_space::muv not found in numerical_masses"))
                .to_owned(),
        );
    } else {
        eval_params.insert(
            "muvsq".into(),
            numerical_masses
                .get("muvsq")
                .unwrap_or_else(|| panic!("muvsq not found in numerical_masses"))
                .to_owned(),
        );
    }
    eval_params.insert(
        "mursq".into(),
        numerical_masses
            .get("mursq")
            .unwrap_or_else(|| panic!("mursq not found in numerical_masses"))
            .to_owned(),
    );

    let integral = vakint.to_canonical(integra_view, true).unwrap();

    let integral_reduced = if [evaluation_orders.0, evaluation_orders.1]
        .iter()
        .any(|(_a, do_reduction)| *do_reduction)
    {
        vakint.tensor_reduce(integral.as_view()).unwrap()
    } else {
        integral.clone()
    };

    debug!(
        "Evaluating integral with evaluation order: {}",
        format!("{}", mod_evaluation_order_a).green(),
    );
    let benchmark_evaluated_integral = match vakint.evaluate_integral(if evaluation_orders.0.1 {
        integral_reduced.as_view()
    } else {
        integral.as_view()
    }) {
        Ok(eval) => eval,
        Err(e) => {
            panic!(
                "Error during benchmark integral evaluation with {} :: error :\n{}",
                format!("{}", mod_evaluation_order_a).red(),
                e
            );
        }
    };

    let (benchmark_central, benchmark_error) = match Vakint::full_numerical_evaluation(
        &vakint.settings,
        benchmark_evaluated_integral.as_view(),
        &eval_params,
        &HashMap::default(),
        Some(&numerical_external_momenta),
    ) {
        Ok(eval) => eval,
        Err(e) => {
            panic!(
                "Error during numerical evaluation of benchmark result with {} :: error:\n{}",
                format!("{}", mod_evaluation_order_a).red(),
                e
            );
        }
    };
    debug!(
        "Benchmark {} :: central :\n{}",
        format!("{}", mod_evaluation_order_a).green(),
        benchmark_central.clone()
    );
    if benchmark_error.is_some() {
        debug!(
            "Benchmark {} :: error   :\n{}",
            format!("{}", mod_evaluation_order_a).green(),
            benchmark_error.as_ref().unwrap().clone()
        );
    }

    // Now perform the second evaluation
    vakint.settings.evaluation_order = mod_evaluation_order_b.clone();
    debug!(
        "Evaluating integral with evaluation order: {}",
        format!("{}", mod_evaluation_order_b).green(),
    );
    let comparison_eval = match vakint.evaluate_integral(if evaluation_orders.1.1 {
        integral_reduced.as_view()
    } else {
        integral.as_view()
    }) {
        Ok(eval) => eval,
        Err(e) => {
            panic!(
                "Error during comparison valuation with: {} :: error :\n{}",
                format!("{}", mod_evaluation_order_b).red(),
                e
            );
        }
    };
    let (tested_central, tested_error) = match Vakint::full_numerical_evaluation(
        &vakint.settings,
        comparison_eval.as_view(),
        &eval_params,
        &HashMap::default(),
        Some(&numerical_external_momenta),
    ) {
        Ok(eval) => eval,
        Err(e) => {
            panic!(
                "Error during numerical evaluation of benchmark result with {} :: error:\n{}",
                format!("{}", mod_evaluation_order_b).red(),
                e
            );
        }
    };
    debug!(
        "Tested {} :: central :\n{}",
        format!("{}", mod_evaluation_order_b).green(),
        tested_central.clone()
    );
    if tested_error.is_some() {
        debug!(
            "Tested {} :: error   :\n{}",
            format!("{}", mod_evaluation_order_b).green(),
            tested_error.as_ref().unwrap().clone()
        );
    }

    let mut combined_error = match (&benchmark_error, &tested_error) {
        (Some(b), Some(t)) => Some(b.aggregate_errors(t)),
        (Some(b), None) => Some(b.clone()),
        (None, Some(t)) => Some(t.clone()),
        (None, None) => None,
    };
    let (matches, msg) = benchmark_central.does_approx_match(
        &tested_central,
        combined_error.as_ref(),
        rel_threshold,
        max_pull,
    );
    if !matches || !quiet {
        println!("\n{}\n", "<><><><><>".green());
        println!(
            "Benchmark {} :: central :\n{}",
            format!("{}", mod_evaluation_order_a).green(),
            benchmark_central
        );
        if benchmark_error.is_some() {
            println!(
                "Benchmark {} :: error   :\n{}",
                format!("{}", mod_evaluation_order_a).green(),
                benchmark_error.unwrap()
            );
        }
        println!(
            "Tested {} :: central :\n{}",
            format!("{}", mod_evaluation_order_b).green(),
            tested_central
        );
        if tested_error.is_some() {
            println!(
                "Tested {} :: error   :\n{}",
                format!("{}", mod_evaluation_order_b).green(),
                tested_error.unwrap()
            );
        }
        println!("{}", msg);
        println!("\n{}\n", "<><><><><>".green());
    } else {
        println!("\n{}\n", "<><><><><>".green());
        info!(
            "Benchmark {} :: central :\n{}",
            format!("{}", mod_evaluation_order_a).green(),
            benchmark_central
        );
        if benchmark_error.is_some() {
            info!(
                "Benchmark {} :: error   :\n{}",
                format!("{}", mod_evaluation_order_a).green(),
                benchmark_error.unwrap()
            );
        }
        info!(
            "Tested {} :: central :\n{}",
            format!("{}", mod_evaluation_order_b).green(),
            tested_central
        );
        if tested_error.is_some() {
            info!(
                "Tested {} :: error   :\n{}",
                format!("{}", mod_evaluation_order_b).green(),
                tested_error.unwrap()
            );
        }
        info!("{}", msg);
        println!("\n{}\n", "<><><><><>".green());
    }
    assert!(matches, "Benchmark and numerical result do not match.");
}

#[allow(unused, clippy::too_many_arguments)]
pub fn compare_vakint_evaluation_vs_reference(
    vakint_default_settings: VakintSettings,
    evaluation_order: EvaluationOrder,
    integra_view: AtomView,
    numerical_masses: HashMap<String, Float, ahash::RandomState>,
    numerical_external_momenta: HashMap<usize, Momentum, ahash::RandomState>,
    expected_output: Vec<(i64, (String, String))>,
    prec: u32,
    max_pull: f64,
) {
    // Adjust evaluation method options
    let mut mod_evaluation_order = evaluation_order.clone();
    mod_evaluation_order.adjust(
        None,
        10.0_f64.powi(-(prec as i32)),
        &numerical_masses,
        &HashMap::default(),
        &numerical_external_momenta,
    );

    // First perform Vakint evaluation
    let mut vakint_settings = VakintSettings {
        allow_unknown_integrals: true,
        use_dot_product_notation: true,
        mu_r_sq_symbol: "vakint::mursq".into(),
        run_time_decimal_precision: prec,
        evaluation_order: mod_evaluation_order.clone(),
        ..vakint_default_settings
    };

    let mut vakint = get_vakint(vakint_settings);

    let integral = vakint.to_canonical(integra_view, true).unwrap();

    let integral = vakint.tensor_reduce(integral.as_view()).unwrap();

    let integral = vakint
        .evaluate_integral(integral.as_view())
        .unwrap_or_else(|op| panic!("Failed to evaluate integral: {}", op));

    debug!("Result before numerical evaluation:\n{}", integral);
    let (result, error) = Vakint::full_numerical_evaluation(
        &vakint.settings,
        integral.as_view(),
        &numerical_masses,
        &HashMap::default(),
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
        println!(
            "Result from {}:\n{}",
            format!("{}", mod_evaluation_order).green(),
            result
        );
        if error.is_some() {
            println!("Error:\n{}", error.unwrap());
        }
        println!("Reference:\n{}", reference);
        println!("{}", msg)
    } else {
        debug!(
            "Result from {}:\n{}",
            format!("{}", mod_evaluation_order).green(),
            result
        );
        if error.is_some() {
            debug!("Error:\n{}", error.unwrap());
        }
        debug!("Reference:\n{}", reference);
        debug!("{}", msg)
    }
    assert!(matches, "Vakint result and reference result do not match");
}
