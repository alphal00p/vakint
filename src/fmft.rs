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
use symbolica::printer::{AtomPrinter, PrintOptions};
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    domains::{integer::Integer, rational::Rational},
    function,
    id::{Condition, PatternOrMap},
};

use crate::{
    get_integer_from_atom, number_condition, symbol_condition, symbols::S, FMFTOptions, TEMPLATES,
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
        res = res.replace_all(
            &vk_parse!("d").unwrap().to_pattern(),
            vk_parse!("4-2*ep").unwrap().to_pattern(),
            None,
            None,
        );
        // Temporarily work with the variable "ep" instead of the UTF-8 symbol for epsilon
        res = res.replace_all(
            &Atom::new_var(vk_symbol!(self.settings.epsilon_symbol.as_str())).to_pattern(),
            vk_parse!("ep").unwrap().to_pattern(),
            None,
            None,
        );
        res = res.replace_all(
            &function!(S.vkdot, function!(S.p, S.id1_a), function!(S.p, S.id2_a)).to_pattern(),
            function!(S.dot, function!(S.p, S.id1_a), function!(S.p, S.id2_a)).to_pattern(),
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
        res = res.replace_all(
            &vk_parse!("Gam(x_,y_)").unwrap().to_pattern(),
            vk_parse!("exp(ep*y_*EulerGamma)*Gamma(x_+ep*y_)").unwrap().to_pattern(),
            None,
            None,
        );
        res = res.replace_all(
            &vk_parse!("iGam(x_,y_)").unwrap().to_pattern(),
            vk_parse!("exp(-ep*y_*EulerGamma)/Gamma(x_+ep*y_)").unwrap().to_pattern(),
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
                res = res.replace_all(
                    &src.to_pattern(),
                    trgt.to_pattern(),
                    Some(restriction),
                    None,
                );
            }
            res
        }));
        if let Some(m) = r
            .pattern_match(&vk_parse!("GammaArgs(x_,y_)").unwrap().into(), None, None)
            .next()
        {
            return Err(VakintError::FMFTError(
            format!("FMFT result contains a Gamma function whose numerical evaluation is not implemented in vakint: Gamma({}+{}*{})",
                m.get(&S.x_).unwrap(),
                m.get(&S.y_).unwrap(),
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
                            vk_symbol!("ep"),
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
                res = res.replace_all(
                    &src.to_pattern(),
                    trgt.to_pattern(),
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

        r = r.replace_all(
            &vk_parse!("Gamma(n_)").unwrap().into(),
            PatternOrMap::Map(Box::new(move |match_in| {
                let n = get_integer_from_atom(match_in.get(S.n_).unwrap().to_atom().as_view())
                    .unwrap() as u32;
                Atom::new_num(Integer::factorial(n - 1))
            })),
            Some(&Condition::from((S.n_, gt_condition(0)))),
            None,
        );
        r.repeat_map(Box::new(move |av: AtomView| {
            let mut res = av.to_owned();
            for (src, trgt) in processed_constants.iter() {
                res = res.replace_all(&src.to_pattern(), trgt.to_pattern(), None, None);
            }
            res
        }));

        if let Some(m) = r
            .pattern_match(&vk_parse!("PolyGamma(x_,y_)").unwrap().into(), None, None)
            .next()
        {
            return Err(VakintError::FMFTError(
            format!("FMFT result contains a PolyGamma function whose numerical evaluation is not implemented in vakint: PolyGamma({},{})",
                m.get(&S.x_).unwrap(),
                m.get(&S.y_).unwrap()
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
                res = res.replace_all(&src.to_pattern(), trgt.to_pattern(), None, None);
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
            if let Some(m) = short_expression
                .pattern_match(
                    &function!(S.fun_, S.any_a___).to_pattern(),
                    Some(&Condition::from((S.fun_, symbol_condition()))),
                    None,
                )
                .next()
            {
                AtomPrinter::new_with_options(m.get(&S.fun_).unwrap().as_atom_view(), PrintOptions::file_no_namespace()).to_string()
            } else {
                return err;
            }
        } else {
            return err;
        };

        let muv_sq_symbol = if let Some(m) = integral
            .short_expression
            .as_ref()
            .unwrap()
            .pattern_match(
                &function!(S.fun_, S.x_a, S.any_a___).to_pattern(),
                Some(
                    &(Condition::from((S.x_, symbol_condition()))
                        & Condition::from((S.fun_, symbol_condition()))),
                ),
                None,
            )
            .next()
        {
            match m.get(&S.x_).unwrap() {
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

        if utils::could_match(
            &function!(S.dot, function!(S.p, S.id1_a), function!(S.k, S.id2_a)).to_pattern(),
            numerator.as_view(),
        ) {
            return Err(VakintError::InvalidNumerator(
                format!("Make sure the numerator has been tensor-reduced before being processed by FMFT : {}", numerator)
            ));
        }

        // Now map all exterior dot products into a special function `vkdot` so that it does not interfere with FMFT
        // Do not forget to normalize by the dimensionality, which is muv^2 in this case.
        numerator = numerator.replace_all(
            &function!(S.dot, function!(S.p, S.id1_a), function!(S.p, S.id2_a)).to_pattern(),
            function!(S.vkdot, function!(S.p, S.id1_a), function!(S.p, S.id2_a)).to_pattern(),
            Some(
                &(Condition::from((S.id1_, number_condition()))
                    & Condition::from((S.id2_, number_condition()))),
            ),
            None,
        );

        // And finally map all interior products into the form p<i>.p<i> expected by fmft, where <i> is the edge id carrying momentum k<i>.
        let momenta: HashMap<usize, Atom> = integral_specs.get_propagator_property_list("q_");
        let mut lmb_prop_indices = vec![];
        for i_loop in 1..=integral.n_loops {
            if let Some((i_edge, _)) = momenta.iter().find(|(_i_prop, k)| {
                k.as_view() == function!(S.k, Atom::new_num(i_loop as i64)).as_view()
            }) {
                lmb_prop_indices.push(*i_edge as isize);
            } else {
                lmb_prop_indices.push(-1);
            }
        }

        let vakint_to_fmft_edge_map_copy = vakint_to_fmft_edge_map.clone();
        numerator = 
        numerator.replace_all(
                &function!(S.dot, function!(S.k, S.id1_a), function!(S.k, S.id2_a))
            .to_pattern(),
                PatternOrMap::Map(Box::new(move |match_in| {
                    let id1 = lmb_prop_indices[(get_integer_from_atom(match_in.get(S.id1_).unwrap().to_atom().as_view()).unwrap()-1) as usize];
                    if id1 < 0 {
                        panic!(
                            "Could not find LMB edge for momentum k({}) in a topology supported by FMFT and used in numerator.",
                            get_integer_from_atom(match_in.get(S.id1_).unwrap().to_atom().as_view()).unwrap()
                        );
                    }
                    let id2 = lmb_prop_indices[(get_integer_from_atom(match_in.get(S.id2_).unwrap().to_atom().as_view()).unwrap()-1) as usize];
                    if id2 < 0 {
                        panic!(
                            "Could not find LMB edge for momentum k({}) in a topology supported by FMFT and used in numerator.",
                            get_integer_from_atom(match_in.get(S.id2_).unwrap().to_atom().as_view()).unwrap()
                        );
                    }
                    let i_edge1 =
                        vakint_to_fmft_edge_map_copy[(id1 as usize) - 1];
                    let i_edge2 =
                        vakint_to_fmft_edge_map_copy[(id2 as usize) - 1];
                    // in FMFT, the loop momenta dot products need to be written p<i>.p<i>
                    // The outter dot will be converted to an inner dot in the next step
                    // Again we must normalize by the dimensionality muv^2
                    vk_parse!(format!("dot(p{},p{})", i_edge1, i_edge2).as_str()).unwrap()
                        * (Atom::new_var(muv_sq_symbol))
                })),
                Some(
                    &(Condition::from((S.id1_, number_condition()))
                        & Condition::from((S.id2_, number_condition()))),
                ),
                None,
            );

        // Substitute eps by (4-d)/2
        numerator = 
        numerator.replace_all(
                &vk_parse!(&vakint.settings.epsilon_symbol).unwrap().into(),
                vk_parse!("(4-d)/2").unwrap().to_pattern(),
                None,
                None,
            );

        let dot_produce_replacer =
            Regex::new(r"dot\((?<vecA>[\w|\d]+),(?<vecB>[\w|\d]+)\)").unwrap();
        let numerator_string = dot_produce_replacer
            .replace_all(
                &AtomPrinter::new_with_options(
                    numerator.as_view(),
                    PrintOptions::file_no_namespace(),
                )
                .to_string(), "($vecA.$vecB)")
            .to_string();

        let powers = (1..=integral.n_props)
            .map(|i_prop| {
                integral_specs
                    .canonical_expression_substitutions
                    .get(&function!(S.pow, Atom::new_num(i_prop as i64)))
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

        // Replace functions with 1 and get all remaining symbols
        let mut numerator_additional_symbols = 
            input_numerator.replace_all(
                &vk_parse!("f_(args__)").unwrap().into(),
                vk_parse!("1").unwrap().to_pattern(),
                None,
                None,
            )
            .get_all_symbols(false);
        let eps_symbol = vk_symbol!(vakint.settings.epsilon_symbol.clone());
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
                        .map(|item| item.get_stripped_name())
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

        // Offset dimensionality by the denominators
        let muv_sq_dimension = 2 * (integral.n_loops as i64) - powers.iter().sum::<i64>();

        evaluated_integral = evaluated_integral
            * vk_parse!(format!("{}^{}", muv_sq_symbol, muv_sq_dimension).as_str()).unwrap();

        let fmft_normalization_correction = vk_parse!(format!(
            "( 
                (𝑖*(𝜋^((4-2*{eps})/2)))\
              * (exp(-EulerGamma))^({eps})\
              * (exp(-logmUVmu-log_mu_sq))^({eps})\
             )^{n_loops}",
            eps = self.settings.epsilon_symbol,
            n_loops = integral.n_loops
        )
        .as_str())
        .unwrap();

        // Adjust normalization factor
        let mut complete_normalization = fmft_normalization_correction
            * vakint
            .settings
            .get_integral_normalization_factor_atom()?.replace_all(
                &S.n_loops.to_pattern(),
                Atom::new_num(integral.n_loops as i64).to_pattern(),
                None,
                None,
            );
        complete_normalization = 
        complete_normalization.replace_all(
                &Atom::new_var(vk_symbol!(self.settings.epsilon_symbol.as_str()))
            .to_pattern(),
                vk_parse!("ep").unwrap().to_pattern(),
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
                vk_symbol!("ep"),
                Atom::Zero.as_view(),
                Rational::from(expansion_depth),
                true,
            ) {
                Ok(a) => a,
                Err(e) => return Err(VakintError::SymbolicaError(e.to_string())),
            }
            .to_atom();

            // Sanity check
            if let Some(m) = 
                evaluated_integral.pattern_match(
                    &vk_parse!("Oep(x_,y_)")
                    .unwrap().to_pattern(),
                    None,
                    None,
                )
                .next()
            {
                return Err(VakintError::FMFTError(format!(
                    "FMFT expansion yielded terms beyond expansion depth supported: Oep({},{})",
                    m.get(&S.x_).unwrap(),
                    m.get(&S.y_).unwrap(),
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
                if let Some(m) = 
                    evaluated_integral.pattern_match(
                        &vk_parse!("Oep(x_,y_)")
                        .unwrap().to_pattern(),
                        None,
                        None,
                    )
                    .next()
                {
                    return Err(VakintError::FMFTError(format!(
                        "FMFT expansion yielded terms beyond expansion depth supported: Oep({},{})",
                        m.get(&S.x_).unwrap(),
                        m.get(&S.y_).unwrap(),
                    )));
                }
            }
        }

        evaluated_integral = evaluated_integral.replace_all(
            &vk_parse!("ep").unwrap().into(),
            Atom::new_var(vk_symbol!(self.settings.epsilon_symbol.as_str()))
                .to_pattern(),
            None,
            None,
        );

        if !vakint.settings.use_dot_product_notation {
            evaluated_integral =
                Vakint::convert_from_dot_notation(evaluated_integral.as_view(), false);
        }

        let log_muv_mu_sq = function!(
            Atom::LOG,
            Atom::new_var(muv_sq_symbol)
                / Atom::new_var(vk_symbol!(vakint.settings.mu_r_sq_symbol.as_str()))
        );

        let log_mu_sq = function!(
            Atom::LOG,
            Atom::new_var(vk_symbol!(vakint.settings.mu_r_sq_symbol.as_str()))
        );

        evaluated_integral = evaluated_integral.replace_all(
            &vk_parse!("logmUVmu").unwrap().into(),
            (log_muv_mu_sq).to_pattern(),
            None,
            None,
        );
        evaluated_integral = evaluated_integral.replace_all(
            &vk_parse!("log_mu_sq").unwrap().into(),
            (log_mu_sq).to_pattern(),
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
