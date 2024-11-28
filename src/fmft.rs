use std::collections::HashMap;

use crate::fmft_numerics::{MASTERS_EXPANSION, MASTERS_NUMERIC_SUBSTITUTIONS};
use crate::utils::set_precision_in_float_atom;
use crate::utils::{self, set_precision_in_polynomial_atom};
use crate::{
    fmft_numerics::{ADDITIONAL_CONSTANTS, POLY_GAMMA_SUBSTITUTIONS},
    gt_condition,
};
use colored::Colorize;
use log::debug;
use regex::Regex;
use string_template_plus::{Render, RenderOptions, Template};
use symbolica::{
    atom::{Atom, AtomView},
    domains::{integer::Integer, rational::Rational},
    fun,
    id::{Condition, MatchSettings, Pattern, PatternOrMap},
    state::State,
};

use crate::{
    get_integer_from_match, number_condition, symbol_condition, symbols::S, FMFTOptions, TEMPLATES,
};

use crate::{ReplacementRules, Vakint, VakintError, VakintSettings};

pub struct FMFT {
    settings: VakintSettings,
}

impl FMFT {
    pub fn with_settings(settings: VakintSettings) -> Self {
        FMFT { settings }
    }

    pub fn process_fmft_form_output(
        &self,
        processed_form_output: Atom,
    ) -> Result<Atom, VakintError> {
        let mut res = processed_form_output.to_owned();
        res = Pattern::parse("d").unwrap().replace_all(
            res.as_view(),
            &Pattern::parse("4-2*ep").unwrap().into(),
            None,
            None,
        );
        // Temporarily work with the variable "ep" instead of the UTF-8 symbol for epsilon
        res = Atom::new_var(State::get_symbol(self.settings.epsilon_symbol.as_str()))
            .into_pattern()
            .replace_all(
                res.as_view(),
                &Pattern::parse("ep").unwrap().into(),
                None,
                None,
            );
        res = fun!(S.vkdot, fun!(S.p, S.id1_a), fun!(S.p, S.id2_a))
            .into_pattern()
            .replace_all(
                res.as_view(),
                &fun!(S.dot, fun!(S.p, S.id1_a), fun!(S.p, S.id2_a))
                    .into_pattern()
                    .into(),
                Some(
                    &(Condition::from((S.id1_, number_condition()))
                        & Condition::from((S.id2_, number_condition()))),
                ),
                None,
            );
        Ok(res)
    }

    pub fn substitute_gam_functions(&self, result: AtomView) -> Atom {
        let mut res = result.to_owned();
        res = Pattern::parse("Gam(x_,y_)").unwrap().replace_all(
            res.as_view(),
            &Pattern::parse("exp(ep*y_*EulerGamma)*Gamma(x_+ep*y_)")
                .unwrap()
                .into(),
            None,
            None,
        );
        res = Pattern::parse("iGam(x_,y_)").unwrap().replace_all(
            res.as_view(),
            &Pattern::parse("exp(-ep*y_*EulerGamma)/Gamma(x_+ep*y_)")
                .unwrap()
                .into(),
            None,
            None,
        );
        res
    }

    pub fn expand_masters(&self, result: AtomView) -> Result<Atom, VakintError> {
        let mut r = result.to_owned();
        r.repeat_map(Box::new(move |av: AtomView| {
            let mut res = av.to_owned();
            res = self.substitute_gam_functions(res.as_view());
            for (src, (trgt, restriction)) in MASTERS_EXPANSION.iter() {
                res = src.into_pattern().replace_all(
                    res.as_view(),
                    &trgt.into_pattern().into(),
                    Some(restriction),
                    None,
                );
            }
            res
        }));
        if let Some(m) = Pattern::parse("GammaArgs(x_,y_)")
            .unwrap()
            .pattern_match(
                r.as_view(),
                &Condition::default(),
                &MatchSettings::default(),
            )
            .next()
        {
            return Err(VakintError::FMFTError(
            format!("FMFT result contains a Gamma function whose numerical evaluation is not implemented in vakint: Gamma({}+{}*{})",
                m.match_stack.get(S.x_).unwrap().to_atom(),
                m.match_stack.get(S.y_).unwrap().to_atom(),
                self.settings.epsilon_symbol
            ),
        ));
        }
        Ok(r)
    }

    pub fn substitute_masters(&self, result: AtomView) -> Result<Atom, VakintError> {
        let processed_constants = MASTERS_NUMERIC_SUBSTITUTIONS
            .iter()
            .map(|(src, (trgt, condition))| {
                (
                    src,
                    (
                        set_precision_in_polynomial_atom(
                            trgt.as_view(),
                            State::get_symbol("ep"),
                            &self.settings,
                        ),
                        condition.clone(),
                    ),
                )
            })
            .collect::<Vec<_>>();
        // for (a, (b, _c)) in processed_constants.clone() {
        //     println!("{} -> {}", a, b);
        // }
        let mut r = result.to_owned();
        r.repeat_map(Box::new(move |av: AtomView| {
            let mut res = av.to_owned();
            for (src, (trgt, matching_condition)) in processed_constants.iter() {
                res = src.into_pattern().replace_all(
                    res.as_view(),
                    &trgt.into_pattern().into(),
                    Some(matching_condition),
                    None,
                );
            }
            res
        }));
        // println!("DONE! {}", r);
        Ok(r)
    }

    pub fn substitute_poly_gamma(&self, result: AtomView) -> Result<Atom, VakintError> {
        let processed_constants = POLY_GAMMA_SUBSTITUTIONS
            .iter()
            .map(|(src, trgt)| {
                (
                    src,
                    set_precision_in_float_atom(trgt.as_view(), &self.settings),
                )
            })
            .collect::<Vec<_>>();
        let mut r = result.to_owned();

        r = Pattern::parse("Gamma(n_)").unwrap().replace_all(
            r.as_view(),
            &PatternOrMap::Map(Box::new(move |match_in| {
                let n = get_integer_from_match(match_in.get(S.n_).unwrap()).unwrap() as u32;
                Atom::new_num(Integer::factorial(n - 1))
            })),
            Some(&Condition::from((S.n_, gt_condition(0)))),
            None,
        );
        r.repeat_map(Box::new(move |av: AtomView| {
            let mut res = av.to_owned();
            for (src, trgt) in processed_constants.iter() {
                res = src.into_pattern().replace_all(
                    res.as_view(),
                    &trgt.into_pattern().into(),
                    None,
                    None,
                );
            }
            res
        }));

        if let Some(m) = Pattern::parse("PolyGamma(x_,y_)")
            .unwrap()
            .pattern_match(
                r.as_view(),
                &Condition::default(),
                &MatchSettings::default(),
            )
            .next()
        {
            return Err(VakintError::FMFTError(
            format!("FMFT result contains a PolyGamma function whose numerical evaluation is not implemented in vakint: PolyGamma({},{})",
                m.match_stack.get(S.x_).unwrap().to_atom(),
                m.match_stack.get(S.y_).unwrap().to_atom()
            ),
        ));
        }
        Ok(r)
    }

    pub fn substitute_additional_constants(&self, result: AtomView) -> Result<Atom, VakintError> {
        let processed_constants = ADDITIONAL_CONSTANTS
            .iter()
            .map(|(src, trgt)| {
                (
                    src,
                    set_precision_in_float_atom(trgt.as_view(), &self.settings),
                )
            })
            .collect::<Vec<_>>();
        let mut r = result.to_owned();
        r.repeat_map(Box::new(move |av: AtomView| {
            let mut res = av.to_owned();
            for (src, trgt) in processed_constants.iter() {
                res = src.into_pattern().replace_all(
                    res.as_view(),
                    &trgt.into_pattern().into(),
                    None,
                    None,
                );
            }
            res
        }));

        Ok(r)
    }
}

impl Vakint {
    pub fn fmft_evaluate(
        &self,
        vakint: &Vakint,
        input_numerator: AtomView,
        integral_specs: &ReplacementRules,
        options: &FMFTOptions,
    ) -> Result<Atom, VakintError> {
        let integral = integral_specs.canonical_topology.get_integral();

        debug!(
            "Processing the following integral with {}:\n{}",
            "FMFT".green(),
            integral
        );

        let err = Err(VakintError::InvalidGenericExpression(format!(
            "Could not find the shorthand integral name in the FMFT expression: {}",
            integral
                .short_expression
                .as_ref()
                .map(|a| a.to_canonical_string())
                .unwrap_or("None".to_string())
        )));

        let fmft = FMFT::with_settings(vakint.settings.clone());

        let integral_name = if let Some(short_expression) = integral.short_expression.as_ref() {
            if let Some(m) = fun!(S.fun_, S.any_a___)
                .into_pattern()
                .pattern_match(
                    short_expression.as_view(),
                    &Condition::from((S.fun_, symbol_condition())),
                    &MatchSettings::default(),
                )
                .next()
            {
                m.match_stack
                    .get(S.fun_)
                    .unwrap()
                    .to_atom()
                    .to_canonical_string()
            } else {
                return err;
            }
        } else {
            return err;
        };

        let muv_sq_symbol = if let Some(m) = fun!(S.fun_, S.x_a, S.any_a___)
            .into_pattern()
            .pattern_match(
                integral.short_expression.as_ref().unwrap().as_view(),
                &(Condition::from((S.x_, symbol_condition()))
                    & Condition::from((S.fun_, symbol_condition()))),
                &MatchSettings::default(),
            )
            .next()
        {
            match m.match_stack.get(S.x_).unwrap().to_atom() {
                Atom::Var(s) => s.get_symbol(),
                _ => {
                    return Err(VakintError::MalformedGraph(format!(
                        "Could not find muV in graph:\n{}",
                        integral.short_expression.as_ref().unwrap()
                    )));
                }
            }
        } else {
            return Err(VakintError::MalformedGraph(format!(
                "Could not find muV in graph:\n{}",
                integral.short_expression.as_ref().unwrap()
            )));
        };

        // Here we map the propagators in the correct order for the definition of the topology in FMFT
        let vakint_to_fmft_edge_map = match integral_name.as_str().split("_pinch_").next().unwrap()
        {
            "I4L_H" => vec![2, 3, 4, 5, 6, 7, 8, 9, 1],
            "I4L_X" => vec![2, 3, 4, 5, 6, 7, 8, 9, 10],
            "I4L_BMW" => vec![3, 4, 5, 6, 7, 8, 9, 10],
            "I4L_FG" => vec![1, 3, 4, 5, 6, 7, 8, 9],
            _ => {
                return Err(VakintError::InvalidGenericExpression(format!(
                    "Integral {} is not supported by FMFT.",
                    integral_name
                )))
            }
        };

        let mut numerator = Vakint::convert_to_dot_notation(input_numerator);

        let mut numerator_for_scaling_analysis = numerator.clone();
        // Rescale dot products
        numerator_for_scaling_analysis = fun!(S.dot, fun!(S.x_, S.id1_a), fun!(S.y_, S.id2_a))
            .into_pattern()
            .replace_all(
                numerator_for_scaling_analysis.as_view(),
                &(&S.lambda_a * fun!(S.vkdot, fun!(S.x_, S.id1_a), fun!(S.y_, S.id2_a)))
                    .into_pattern()
                    .into(),
                Some(
                    &(Condition::from((S.id1_, number_condition()))
                        & Condition::from((S.id2_, number_condition()))),
                ),
                None,
            );
        // Rescale the UV mass
        numerator_for_scaling_analysis = Atom::new_var(muv_sq_symbol).into_pattern().replace_all(
            numerator_for_scaling_analysis.as_view(),
            &(&S.lambda_a * Atom::new_var(muv_sq_symbol))
                .into_pattern()
                .into(),
            Some(
                &(Condition::from((S.id1_, number_condition()))
                    & Condition::from((S.id2_, number_condition()))),
            ),
            None,
        );
        let (scaling_coefficients, remainder) =
            numerator_for_scaling_analysis.coefficient_list(S.lambda);

        let mut muv_sq_dimension = if remainder.is_zero() && scaling_coefficients.len() == 1 {
            let (coef, _) = scaling_coefficients[0].clone();
            if let Some(m) = Pattern::parse(format!("{}^{}", S.lambda, S.pow_).as_str())
                .unwrap()
                .pattern_match(
                    coef.as_view(),
                    &Condition::default(),
                    &MatchSettings::default(),
                )
                .next()
            {
                get_integer_from_match(m.match_stack.get(S.pow_).unwrap()).unwrap()
            } else {
                return Err(VakintError::FMFTError(format!(
                    "The numerator does not appear to have uniform dimensionality: {}",
                    numerator_for_scaling_analysis
                )));
            }
        } else if scaling_coefficients.is_empty() {
            0
        } else {
            return Err(VakintError::FMFTError(format!(
                "The numerator does not appear to have uniform dimensionality: {}",
                numerator_for_scaling_analysis
            )));
        };

        // Set the numerator to one and restore it using uniform dimensionality later
        numerator = Atom::new_var(muv_sq_symbol).into_pattern().replace_all(
            numerator.as_view(),
            &Pattern::parse("1").unwrap().into(),
            None,
            None,
        );

        if utils::could_match(
            &fun!(S.dot, fun!(S.p, S.id1_a), fun!(S.k, S.id2_a)).into_pattern(),
            numerator.as_view(),
        ) {
            return Err(VakintError::InvalidNumerator(
                format!("Make sure the numerator has been tensor-reduced before being processed by FMFT : {}", numerator)
            ));
        }

        // Now map all exterior dot products into a special function `vkdot` so that it does not interfere with FMFT
        // Do not forget to normalize by the dimensionality, which is muv^2 in this case.
        numerator = fun!(S.dot, fun!(S.p, S.id1_a), fun!(S.p, S.id2_a))
            .into_pattern()
            .replace_all(
                numerator.as_view(),
                &(fun!(S.vkdot, fun!(S.p, S.id1_a), fun!(S.p, S.id2_a))
                    / (Atom::new_var(muv_sq_symbol)))
                .into_pattern()
                .into(),
                Some(
                    &(Condition::from((S.id1_, number_condition()))
                        & Condition::from((S.id2_, number_condition()))),
                ),
                None,
            );

        // And finally map all interior products into the form p<i>.p<i> expected by fmft, where <i> is the edge id carrying momentum k<i>.
        let momenta = integral_specs.get_propagator_property_list("q_");
        let mut lmb_prop_indices = vec![];
        for i_loop in 1..=integral.n_loops {
            if let Some((i_edge, _)) = momenta.iter().find(|(_i_prop, k)| {
                k.as_view() == fun!(S.k, Atom::new_num(i_loop as i64)).as_view()
            }) {
                lmb_prop_indices.push(*i_edge as isize);
            } else {
                lmb_prop_indices.push(-1);
            }
        }

        let vakint_to_fmft_edge_map_copy = vakint_to_fmft_edge_map.clone();
        numerator = fun!(S.dot, fun!(S.k, S.id1_a), fun!(S.k, S.id2_a))
            .into_pattern()
            .replace_all(
                numerator.as_view(),
                &PatternOrMap::Map(Box::new(move |match_in| {
                    let id1 = lmb_prop_indices[(get_integer_from_match(match_in.get(S.id1_).unwrap()).unwrap()-1) as usize];
                    if id1 < 0 {
                        panic!(
                            "Could not find LMB edge for momentum k({}) in a topology supported by FMFT and used in numerator.",
                            get_integer_from_match(match_in.get(S.id1_).unwrap()).unwrap()
                        );
                    }
                    let id2 = lmb_prop_indices[(get_integer_from_match(match_in.get(S.id2_).unwrap()).unwrap()-1) as usize];
                    if id2 < 0 {
                        panic!(
                            "Could not find LMB edge for momentum k({}) in a topology supported by FMFT and used in numerator.",
                            get_integer_from_match(match_in.get(S.id2_).unwrap()).unwrap()
                        );
                    }
                    let i_edge1 =
                        vakint_to_fmft_edge_map_copy[(id1 as usize) - 1];
                    let i_edge2 =
                        vakint_to_fmft_edge_map_copy[(id2 as usize) - 1];
                    // in FMFT, the loop momenta dot products need to be written p<i>.p<i>
                    // The outter dot will be converted to an inner dot in the next step
                    // Again we must normalize by the dimensionality muv^2
                    Atom::parse(format!("dot(p{},p{})", i_edge1, i_edge2).as_str()).unwrap()
                        / (Atom::new_var(muv_sq_symbol))
                })),
                Some(
                    &(Condition::from((S.id1_, number_condition()))
                        & Condition::from((S.id2_, number_condition()))),
                ),
                None,
            );

        // Substitute eps by (4-d)/2
        numerator = Pattern::parse(&vakint.settings.epsilon_symbol)
            .unwrap()
            .replace_all(
                numerator.as_view(),
                &Pattern::parse("(4-d)/2").unwrap().into(),
                None,
                None,
            );

        let dot_produce_replacer =
            Regex::new(r"dot\((?<vecA>[\w|\d]+),(?<vecB>[\w|\d]+)\)").unwrap();
        let numerator_string = dot_produce_replacer
            .replace_all(&numerator.to_canonical_string(), "($vecA.$vecB)")
            .to_string();

        let powers = (1..=integral.n_props)
            .map(|i_prop| {
                integral_specs
                    .canonical_expression_substitutions
                    .get(&fun!(S.pow, Atom::new_num(i_prop as i64)))
                    .map(|a| a.try_into().unwrap())
                    .unwrap_or(0)
            })
            .collect::<Vec<_>>();

        let integral_string = powers
            .iter()
            .zip(vakint_to_fmft_edge_map)
            .map(|(pwr, fmft_edge_index)| format!("d{}^{}", fmft_edge_index, -pwr))
            .collect::<Vec<_>>()
            .join("*");

        // Offset dimensionality by the denominators
        muv_sq_dimension -= powers.iter().sum::<i64>();
        // Include the scaling of the measure
        muv_sq_dimension += 2 * integral.n_loops as i64;
        debug!(
            "{}: The integral was detected to be of overall dimension GeV^{} ...",
            "FMFT".green(),
            2 * muv_sq_dimension
        );

        //println!("FMFT input string: {}", format!("({})*({})", numerator_string, integral_string));

        // Replace functions with 1 and get all remaining symbols
        let mut numerator_additional_symbols = Pattern::parse("f_(args__)")
            .unwrap()
            .replace_all(
                input_numerator,
                &Atom::parse("1").unwrap().into_pattern().into(),
                None,
                None,
            )
            .get_all_symbols(false);
        let eps_symbol = State::get_symbol(vakint.settings.epsilon_symbol.clone());
        numerator_additional_symbols.retain(|&s| s != eps_symbol);

        let template = Template::parse_template(TEMPLATES.get("run_fmft.txt").unwrap()).unwrap();
        let mut vars: HashMap<String, String> = HashMap::new();
        vars.insert("numerator".into(), numerator_string);
        vars.insert("integral".into(), integral_string);
        vars.insert("symbols".into(), muv_sq_symbol.to_string());
        if !numerator_additional_symbols.is_empty() {
            vars.insert(
                "additional_symbols".into(),
                format!(
                    "Auto S {};",
                    numerator_additional_symbols
                        .iter()
                        .map(|item| item.to_string())
                        .collect::<Vec<_>>()
                        .join(", "),
                ),
            );
        } else {
            vars.insert("additional_symbols".into(), "".into());
        }
        let rendered = template
            .render(&RenderOptions {
                variables: vars,
                ..Default::default()
            })
            .unwrap();

        let form_result = vakint.run_form(
            &["fmft.frm".into()],
            ("run_fmft.frm".into(), rendered),
            vec![],
            vakint.settings.clean_tmp_dir,
            vakint.settings.temporary_directory.clone(),
        )?;

        let processed_form_result = self.process_form_output(form_result)?;
        let mut evaluated_integral = fmft.process_fmft_form_output(processed_form_result)?;
        debug!(
            "{}: raw result from FORM:\n{}",
            "FMFT".green(),
            evaluated_integral
        );
        // Restore dimensionality now
        evaluated_integral = evaluated_integral
            * Atom::parse(format!("{}^{}", muv_sq_symbol, muv_sq_dimension).as_str()).unwrap();

        let mut fmft_normalization_correction = Atom::parse(
            format!(
                "( 
                    (𝑖*(𝜋^((4-2*{eps})/2)))\
                  * (exp(-EulerGamma))^({eps})\
                  * (exp(-logmUVmu-log_mu_sq))^({eps})\
                 )^{n_loops}",
                eps = self.settings.epsilon_symbol,
                n_loops = integral.n_loops
            )
            .as_str(),
        )
        .unwrap();

        // Since FMFT uses euclidean denominator, we must adjust the overall sign by (-1) per quadratic denominator with power one.
        fmft_normalization_correction = fmft_normalization_correction
            * Atom::parse(format!("((-1)^{})", powers.iter().sum::<i64>()).as_str()).unwrap();

        // Adjust normalization factor
        let mut complete_normalization = fmft_normalization_correction
            * S.n_loops.into_pattern().replace_all(
                vakint
                    .settings
                    .get_integral_normalization_factor_atom()?
                    .as_view(),
                &Atom::new_num(integral.n_loops as i64).into_pattern().into(),
                None,
                None,
            );
        complete_normalization =
            Atom::new_var(State::get_symbol(self.settings.epsilon_symbol.as_str()))
                .into_pattern()
                .replace_all(
                    complete_normalization.as_view(),
                    &Pattern::parse("ep").unwrap().into(),
                    None,
                    None,
                );

        evaluated_integral = evaluated_integral * complete_normalization;

        if options.expand_masters {
            let expansion_depth = vakint.settings.number_of_terms_in_epsilon_expansion
                - (integral.n_loops as i64)
                - 1;
            debug!(
                "{}: Expanding master integrals with terms up to and including {}^{} ...",
                "FMFT".green(),
                self.settings.epsilon_symbol,
                expansion_depth
            );
            evaluated_integral = fmft.expand_masters(evaluated_integral.as_view())?;

            debug!(
                "{}: Series expansion of the result up to and including terms of order {}^{} ...",
                "FMFT".green(),
                self.settings.epsilon_symbol,
                expansion_depth
            );
            evaluated_integral = match evaluated_integral.series(
                State::get_symbol("ep"),
                Atom::Zero.as_view(),
                Rational::from(expansion_depth),
                true,
            ) {
                Ok(a) => a,
                Err(e) => return Err(VakintError::SymbolicaError(e.to_string())),
            }
            .to_atom();

            // Sanity check
            if let Some(m) = Pattern::parse("Oep(x_,y_)")
                .unwrap()
                .pattern_match(
                    evaluated_integral.as_view(),
                    &Condition::default(),
                    &MatchSettings::default(),
                )
                .next()
            {
                return Err(VakintError::FMFTError(format!(
                    "FMFT expansion yielded terms beyond expansion depth supported: Oep({},{})",
                    m.match_stack.get(S.x_).unwrap().to_atom(),
                    m.match_stack.get(S.y_).unwrap().to_atom(),
                )));
            }

            if options.susbstitute_masters {
                debug!("{}: Substituting master integrals coefficient with their numerical evaluations...", "FMFT".green());
                evaluated_integral = fmft.substitute_masters(evaluated_integral.as_view())?;
                debug!(
                    "{}: Substituting PolyGamma and period constants...",
                    "FMFT".green()
                );
                evaluated_integral = evaluated_integral.expand();
                evaluated_integral = fmft.substitute_poly_gamma(evaluated_integral.as_view())?;
                evaluated_integral =
                    fmft.substitute_additional_constants(evaluated_integral.as_view())?;
                // Sanity check
                if let Some(m) = Pattern::parse("Oep(x_,y_)")
                    .unwrap()
                    .pattern_match(
                        evaluated_integral.as_view(),
                        &Condition::default(),
                        &MatchSettings::default(),
                    )
                    .next()
                {
                    return Err(VakintError::FMFTError(format!(
                        "FMFT expansion yielded terms beyond expansion depth supported: Oep({},{})",
                        m.match_stack.get(S.x_).unwrap().to_atom(),
                        m.match_stack.get(S.y_).unwrap().to_atom(),
                    )));
                }
            }
        }

        evaluated_integral = Pattern::parse("ep").unwrap().replace_all(
            evaluated_integral.as_view(),
            &Atom::new_var(State::get_symbol(self.settings.epsilon_symbol.as_str()))
                .into_pattern()
                .into(),
            None,
            None,
        );

        if !vakint.settings.use_dot_product_notation {
            evaluated_integral =
                Vakint::convert_from_dot_notation(evaluated_integral.as_view(), false);
        }

        let log_muv_mu_sq = fun!(
            State::LOG,
            Atom::new_var(muv_sq_symbol)
                / Atom::new_var(State::get_symbol(vakint.settings.mu_r_sq_symbol.as_str()))
        );

        let log_mu_sq = fun!(
            State::LOG,
            Atom::new_var(State::get_symbol(vakint.settings.mu_r_sq_symbol.as_str()))
        );

        evaluated_integral = Pattern::parse("logmUVmu").unwrap().replace_all(
            evaluated_integral.as_view(),
            &(log_muv_mu_sq).into_pattern().into(),
            None,
            None,
        );
        evaluated_integral = Pattern::parse("log_mu_sq").unwrap().replace_all(
            evaluated_integral.as_view(),
            &(log_mu_sq).into_pattern().into(),
            None,
            None,
        );

        // println!(
        //     "evaluated_integral: {}",
        //     evaluated_integral.to_canonical_string()
        // );

        Ok(evaluated_integral)
    }
}
