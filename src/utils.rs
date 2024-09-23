use symbolica::{
    atom::{Atom, AtomView},
    coefficient::CoefficientView,
    id::{Condition, MatchSettings, Pattern, PatternOrMap, PatternRestriction},
};

use crate::VakintSettings;

pub fn replace_until_stable(
    target: AtomView<'_>,
    lhs: &Pattern,
    rhs: &PatternOrMap,
    conditions: Option<&Condition<PatternRestriction>>,
    settings: Option<&MatchSettings>,
) -> Atom {
    let mut old_res = lhs.replace_all(target, rhs, conditions, settings);
    loop {
        let res = lhs.replace_all(old_res.as_view(), rhs, conditions, settings);
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
    return pattern
        .pattern_match(target, &Condition::default(), &MatchSettings::default())
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

#[test]
fn test_set_precision_in_float_atom() {
    println!(
        "Res={}",
        set_precision_in_float_atom(
            Atom::parse("1.23456789113245423`100").unwrap().as_view(),
            &VakintSettings {
                n_digits_at_evaluation_time: 30,
                ..VakintSettings::default()
            }
        )
    );
}
