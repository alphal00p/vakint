use symbolica::{
    atom::{Atom, AtomView},
    id::{Condition, MatchSettings, Pattern, PatternOrMap, WildcardAndRestriction},
};

pub fn replace_until_stable(
    target: AtomView<'_>,
    lhs: &Pattern,
    rhs: &PatternOrMap,
    conditions: Option<&Condition<WildcardAndRestriction>>,
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
