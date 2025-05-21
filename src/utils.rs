use ahash::HashMap;
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    coefficient::CoefficientView,
    id::{Condition, MatchSettings, Pattern, PatternRestriction},
};
use vakint_macros::vk_parse;

use crate::{eq_condition, VakintSettings};

pub(crate) mod vakint_macros {
    macro_rules! vk_parse {
        ($s: expr) => {{
            symbolica::parse!($s, crate::NAMESPACE)
        }};
        ($s: expr, $ns: expr) => {{
            symbolica::parse!($sm, $ns)
        }};
    }
    macro_rules! vk_symbol {
        ($s: expr) => {{
            if format!("{}", $s).starts_with(format!("{}::", crate::NAMESPACE).as_str()) {
                symbolica::symbol!($s)
            } else {
                symbolica::symbol!(format!("{}::{}", crate::NAMESPACE, $s))
            }
        }};
    }

    pub(crate) use {vk_parse, vk_symbol};
}

#[macro_export]
macro_rules! vakint_parse {
    ($s: expr) => {{
        symbolica::parse!($s, $crate::NAMESPACE)
    }};
    ($s: expr, $ns: expr) => {{
        symbolica::parse!($sm, $ns)
    }};
}

#[macro_export]
macro_rules! vakint_symbol {
    ($s: expr) => {{
        symbolica::symbol!(format!("{}::{}", $crate::NAMESPACE, $s))
    }};
}

pub fn replace_until_stable(
    target: AtomView<'_>,
    lhs: &Pattern,
    rhs: &Pattern,
    conditions: Option<&Condition<PatternRestriction>>,
    settings: Option<&MatchSettings>,
) -> Atom {
    // Specifiying settings is not possible anymore given new Symbolica v0.16 API
    assert!(settings.is_none());

    // We do not want to use `repeat()` below because it does not use the same termination condition
    // which we want in our case to be that the atom is stable.

    // let mut res = target.replace(lhs);
    // if let Some(c) = conditions {
    //     res = res.when(c);
    // }
    // let a = res.repeat().with(rhs);
    // a
    let mut replacer = target.replace(lhs);
    if let Some(c) = conditions {
        replacer = replacer.when(c);
    }
    let mut old_res = replacer.with(rhs);
    loop {
        let mut replacer = old_res.replace(lhs);
        if let Some(c) = conditions {
            replacer = replacer.when(c);
        }
        let res = replacer.with(rhs);
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
        &vk_parse!("x_^y_").unwrap().to_pattern(),
        &vk_parse!("1").unwrap().to_pattern(),
        Some(&Condition::from((crate::vk_symbol!("y_"), eq_condition(0)))),
        None,
    );
    res = replace_until_stable(
        res.as_view(),
        &vk_parse!("log(exp(x_))").unwrap().to_pattern(),
        &vk_parse!("x_").unwrap().to_pattern(),
        None,
        None,
    );
    res = replace_until_stable(
        res.as_view(),
        &vk_parse!("exp(log(x_))").unwrap().to_pattern(),
        &vk_parse!("x_").unwrap().to_pattern(),
        None,
        None,
    );
    res = replace_until_stable(
        res.as_view(),
        &vk_parse!("log(x_^p_)").unwrap().to_pattern(),
        &vk_parse!("p_*log(x_)").unwrap().to_pattern(),
        None,
        None,
    );
    res = replace_until_stable(
        res.as_view(),
        &vk_parse!("exp(x_)*exp(y_)").unwrap().to_pattern(),
        &vk_parse!("exp(x_*y_)").unwrap().to_pattern(),
        None,
        None,
    );
    res
}

pub fn could_match(pattern: &Pattern, target: AtomView) -> bool {
    target
        .pattern_match(
            pattern,
            Some(&Condition::default()),
            Some(&MatchSettings::default()),
        )
        .next()
        .is_some()
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
    let coeffs = input.coefficient_list::<i8>(&[Atom::new_var(variable)]);
    for (var_with_power, coeff) in coeffs {
        let new_coeff = set_precision_in_float_atom(coeff.as_view(), settings);
        res = res + new_coeff * var_with_power;
    }

    res
}

pub fn split_linear_atom(a: AtomView, variable: AtomView) -> (Atom, Atom) {
    let split_atom: HashMap<_, _> = a.coefficient_list::<i8>(&[variable]).into_iter().collect();
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
            vk_parse!("1.23456789113245423`100").unwrap().as_view(),
            &VakintSettings {
                run_time_decimal_precision: 30,
                ..VakintSettings::default()
            }
        )
    );
}
