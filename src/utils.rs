use ahash::HashMap;
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    coefficient::CoefficientView,
    id::{Condition, MatchSettings, Pattern, PatternOrMap, PatternRestriction},
    symb,
};

use crate::{eq_condition, VakintSettings};

pub fn replace_until_stable(
    target: AtomView<'_>,
    lhs: &Pattern,
    rhs: &PatternOrMap,
    conditions: Option<&Condition<PatternRestriction>>,
    settings: Option<&MatchSettings>,
) -> Atom {
    let mut old_res = target.replace_all(lhs, rhs, conditions, settings);
    loop {
        let res = old_res.replace_all(lhs, rhs, conditions, settings);
        if res == old_res {
            break;
        }
        old_res = res;
    }
    old_res
}

pub fn simplify_real(input: AtomView) -> Atom {
    let mut res = replace_until_stable(
        input,
        &Pattern::parse("x_^y_").unwrap(),
        &Pattern::parse("1").unwrap().into(),
        Some(&Condition::from((symb!("y_"), eq_condition(0)))),
        None,
    );
    res = replace_until_stable(
        res.as_view(),
        &Pattern::parse("log(exp(x_))").unwrap(),
        &Pattern::parse("x_").unwrap().into(),
        None,
        None,
    );
    res = replace_until_stable(
        res.as_view(),
        &Pattern::parse("exp(log(x_))").unwrap(),
        &Pattern::parse("x_").unwrap().into(),
        None,
        None,
    );
    res = replace_until_stable(
        res.as_view(),
        &Pattern::parse("log(x_^p_)").unwrap(),
        &Pattern::parse("p_*log(x_)").unwrap().into(),
        None,
        None,
    );
    res = replace_until_stable(
        res.as_view(),
        &Pattern::parse("exp(x_)*exp(y_)").unwrap(),
        &Pattern::parse("exp(x_*y_)").unwrap().into(),
        None,
        None,
    );
    res
}

pub fn could_match(pattern: &Pattern, target: AtomView) -> bool {
    return target
        .pattern_match(
            pattern,
            Some(&Condition::default()),
            Some(&MatchSettings::default()),
        )
        .next()
        .is_some();
}

pub fn set_precision_in_float_atom(input: AtomView, settings: &VakintSettings) -> Atom {
    let binary_prec = settings.get_binary_precision();
    if let AtomView::Num(fl_view) = input {
        match fl_view.get_coeff_view() {
            CoefficientView::Float(fl) => {
                let mut truncated_float = fl.to_float();
                truncated_float.set_prec(binary_prec);
                Atom::new_num(truncated_float)
            }
            _ => input.to_owned(),
        }
    } else {
        input.to_owned()
    }
}

pub fn set_precision_in_polynomial_atom(
    input: AtomView,
    variable: Symbol,
    settings: &VakintSettings,
) -> Atom {
    let mut res = Atom::Zero;
    let coeffs = input.coefficient_list::<i8, _>(&[Atom::new_var(variable)]);
    for (var_with_power, coeff) in coeffs {
        let new_coeff = set_precision_in_float_atom(coeff.as_view(), settings);
        res = res + new_coeff * var_with_power;
    }

    res
}

pub fn split_linear_atom(a: AtomView, variable: AtomView) -> (Atom, Atom) {
    let split_atom: HashMap<_, _> = a
        .coefficient_list::<i8, _>(&[variable])
        .into_iter()
        .collect();
    if split_atom.len() > 2 {
        panic!(
            "Could not split linear atom '{}' into at most two parts using variable '{}'. Split atom:\n{}",
            a, variable, split_atom.iter().map(|(k, v)| format!("{}: {}", k, v)).collect::<Vec<_>>().join("\n")
        );
    } else {
        let res = (
            split_atom
                .get(&variable.to_owned())
                .unwrap_or(&Atom::Zero)
                .to_owned(),
            split_atom
                .get(&Atom::new_num(1))
                .unwrap_or(&Atom::Zero)
                .to_owned(),
        );
        res
    }
}

#[test]
fn test_set_precision_in_float_atom() {
    println!(
        "Res={}",
        set_precision_in_float_atom(
            Atom::parse("1.23456789113245423`100").unwrap().as_view(),
            &VakintSettings {
                run_time_decimal_precision: 30,
                ..VakintSettings::default()
            }
        )
    );
}
