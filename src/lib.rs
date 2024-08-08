use anyhow::Result;
use colored::Colorize;
use std::{
    collections::{HashMap, HashSet},
    fmt,
};
use symbolica::{
    atom::{representation::InlineNum, AsAtomView, Atom, AtomView, SliceType, Symbol},
    coefficient::CoefficientView,
    fun,
    id::{Condition, Match, MatchSettings, Pattern, PatternRestriction, WildcardAndRestriction},
    state::{FunctionAttribute, State},
};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum VakintError {
    #[error("invalid input integral found: '{0}'")]
    InvalidIntegralFormat(String),
    #[error("invalid generic expression in supported integral: {0}")]
    InvalidGenericExpression(String),
    #[error("invalid expression for the momentum of an edge: {0}")]
    InvalidMomentumExpression(String),
    #[error("invalid short expression supplied for integral: {0}")]
    InvalidShortExpression(String),
    #[error(
        "the following integral could not be identified using any of the supported topologies: {0}"
    )]
    UnreckognizedIntegral(String),
    #[error(
        "not all parts of the numerator have been identified w.r.t the canonical \
        expression (make sure to use the name '{0}' for momenta external to the topology).\n\
        This check can be disabled by setting `verify_numerator_identification` to false.\
        \nLeft-over: {1}"
    )]
    NumeratorNotReplaced(String, String),
    #[error("unknown vakint error")]
    Unknown,
}

fn propagators_condition() -> PatternRestriction {
    PatternRestriction::Filter(Box::new(move |m| {
        let props = match m {
            Match::Single(a) => vec![*a],
            Match::Multiple(SliceType::Mul, atoms) => atoms.clone(),
            _ => return false,
        };

        let pattern = Pattern::parse("prop(propID_,uedge(nl_,nr_),q_,mUVsq_,pow_)").unwrap();
        let number_node_condition = Condition::from((State::get_symbol("nl_"), gt_condition(0)))
            & Condition::from((State::get_symbol("nr_"), gt_condition(0)))
            & Condition::from((State::get_symbol("propID_"), gt_condition(0)))
            & Condition::from((State::get_symbol("mUVsq_"), symbol_or_number()))
            & Condition::from((State::get_symbol("pow_"), symbol_or_number()));
        for p in props {
            if pattern
                .pattern_match(p, &number_node_condition, &MatchSettings::default())
                .next()
                .is_none()
            {
                return false;
            }
        }
        true
    }))
}

fn gt_condition(value: i64) -> PatternRestriction {
    PatternRestriction::Filter(Box::new(move |m| {
        if let Match::Single(AtomView::Num(a)) = m {
            a.get_coeff_view() > InlineNum::new(value, 1).as_num_view().get_coeff_view()
        } else {
            false
        }
    }))
}

#[allow(unused)]
fn lt_condition(value: i64) -> PatternRestriction {
    PatternRestriction::Filter(Box::new(move |m| {
        if let Match::Single(AtomView::Num(a)) = m {
            a.get_coeff_view() < InlineNum::new(value, 1).as_num_view().get_coeff_view()
        } else {
            false
        }
    }))
}

fn range_condition(min: i64, max: i64) -> PatternRestriction {
    PatternRestriction::Filter(Box::new(move |m| {
        if let Match::Single(AtomView::Num(a)) = m {
            a.get_coeff_view() <= InlineNum::new(max, 1).as_num_view().get_coeff_view()
                && a.get_coeff_view() >= InlineNum::new(min, 1).as_num_view().get_coeff_view()
        } else {
            false
        }
    }))
}

fn number_condition() -> PatternRestriction {
    PatternRestriction::Filter(Box::new(move |m| {
        matches!(m, Match::Single(AtomView::Num(_)))
    }))
}

fn symbol_or_number() -> PatternRestriction {
    PatternRestriction::Filter(Box::new(move |m| {
        matches!(
            m,
            Match::Single(AtomView::Num(_))
                | Match::Single(AtomView::Var(_))
                | Match::FunctionName(_)
        )
    }))
}

fn symbol_condition() -> PatternRestriction {
    PatternRestriction::Filter(Box::new(move |m| {
        matches!(m, Match::Single(AtomView::Var(_)) | Match::FunctionName(_))
    }))
}

fn apply_restriction_to_symbols(
    symbols: Vec<Symbol>,
    restriction: &PatternRestriction,
) -> Condition<(Symbol, PatternRestriction)> {
    symbols[1..].iter().fold(
        Condition::from((symbols[0], restriction.clone())),
        |acc, &s| acc & Condition::from((s, restriction.clone())),
    )
}

#[derive(Debug, Clone)]
pub struct ReplacementRules {
    canonical_topology: Topology,
    edge_ids_canonical_to_input_map: HashMap<usize, usize>,
    canonical_expression_substitutions: HashMap<Atom, Atom>,
    numerator_substitutions: HashMap<Atom, Atom>,
}

impl fmt::Display for ReplacementRules {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "ReplacementRules {{")?;
        writeln!(f, "  canonical_topology: {}", self.canonical_topology)?;
        writeln!(
            f,
            "  edge_ids_canonical_to_input_map: {{ {} }}",
            self.edge_ids_canonical_to_input_map
                .iter()
                .map(|(k, v)| format!("{} -> {}", k, v))
                .collect::<Vec<String>>()
                .join(", ")
        )?;
        writeln!(
            f,
            "  canonical_expression_substitutions: {{ {} }}",
            self.canonical_expression_substitutions
                .iter()
                .map(|(k, v)| format!("{} -> {}", k, v))
                .collect::<Vec<String>>()
                .join(", ")
        )?;
        writeln!(
            f,
            "  numerator_substitutions: {{ {} }}",
            self.numerator_substitutions
                .iter()
                .map(|(k, v)| format!("{} -> {}", k, v))
                .collect::<Vec<String>>()
                .join(", ")
        )?;
        write!(f, "}}")
    }
}
impl Default for ReplacementRules {
    fn default() -> Self {
        ReplacementRules {
            canonical_topology: Topology::Unknown(Integral::new(0, None, None).unwrap()),
            edge_ids_canonical_to_input_map: HashMap::new(),
            canonical_expression_substitutions: HashMap::new(),
            numerator_substitutions: HashMap::new(),
        }
    }
}

pub struct Topologies(Vec<Topology>);

impl fmt::Display for Topologies {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            self.0
                .iter()
                .enumerate()
                .map(|(i, t)| format!("{}) {}", i + 1, t))
                .collect::<Vec<_>>()
                .join("\n")
        )
    }
}

impl Topologies {
    pub fn generate_topologies_with_contractions(
        n_tot_props: usize,
        canonical_expression: AtomView,
        short_expression: AtomView,
        contractions: Vec<Vec<usize>>,
    ) -> Result<Self, VakintError> {
        let mut topologies = vec![];
        for contracted_prop_indices in contractions {
            let mut contracted_canonical_expression = canonical_expression.to_owned().clone();
            let mut contracted_short_expression = short_expression.to_owned().clone();
            for prop_id in contracted_prop_indices.iter() {
                contracted_canonical_expression =
                    Pattern::parse(format!("prop({},args__)", prop_id).as_str())
                        .unwrap()
                        .replace_all(
                            contracted_canonical_expression.as_view(),
                            &Pattern::parse("1").unwrap(),
                            None,
                            None,
                        );
                contracted_short_expression = Pattern::parse(format!("pow({})", prop_id).as_str())
                    .unwrap()
                    .replace_all(
                        contracted_short_expression.as_view(),
                        &Pattern::parse("0").unwrap(),
                        None,
                        None,
                    );
            }
            /*
            if !contracted_prop_indices.is_empty() {
                let short_integral_symbol = if let Some(m) = Pattern::parse("fn_(args__)")
                    .unwrap()
                    .pattern_match(
                        contracted_short_expression.as_view(),
                        &Condition::default(),
                        &MatchSettings::default(),
                    )
                    .next()
                {
                    if let Match::FunctionName(s) =
                        m.match_stack.get(State::get_symbol("fn_")).unwrap()
                    {
                        s.to_string()
                    } else {
                        return Err(VakintError::InvalidShortExpression(format!(
                            "{}",
                            short_expression
                        )));
                    }
                } else {
                    return Err(VakintError::InvalidShortExpression(format!(
                        "{}",
                        short_expression
                    )));
                };
                contracted_short_expression =
                    Pattern::parse(format!("{}(args__)", short_integral_symbol).as_str())
                        .unwrap()
                        .replace_all(
                            contracted_short_expression.as_view(),
                            &Pattern::parse(
                                format!(
                                    "{}_{}(args__)",
                                    short_integral_symbol,
                                    contracted_prop_indices
                                        .iter()
                                        .map(|c| format!("{}", c))
                                        .collect::<Vec<_>>()
                                        .join("_")
                                )
                                .as_str(),
                            )
                            .unwrap(),
                            None,
                            None,
                        );
            }
            */
            topologies.push(
                Integral::new(
                    n_tot_props,
                    Some(contracted_canonical_expression),
                    Some(contracted_short_expression),
                )?
                .into(),
            );
        }
        Ok(Topologies(topologies))
    }

    pub fn generate_topologies(settings: &VakintSettings) -> Result<Self, VakintError> {
        // One-loop topology
        let mut topologies = Topologies(vec![Integral::new(
            1,
            Some(Atom::parse("topo(prop(1,edge(1,1),k(1),msq(1),pow(1)))").unwrap()),
            Some(Atom::parse("I1LA(msq(1),pow(1))").unwrap()),
        )?
        .into()]);
        // Two-loop topologies
        topologies.0.extend(
            Topologies::generate_topologies_with_contractions(
                3,
                Atom::parse(
                    "topo(\
                        prop(1,edge(1,2),k(1),msq(1),pow(1))\
                        *prop(2,edge(1,2),k(2),msq(1),pow(2))\
                        *prop(3,edge(2,1),k(1)+k(2),msq(1),pow(3))\
                    )",
                )
                .unwrap()
                .as_view(),
                Atom::parse("I2LA(msq(1),pow(1),pow(2),pow(3))")
                    .unwrap()
                    .as_view(),
                vec![vec![], vec![3]],
            )?
            .0,
        );

        if settings.allow_unknown_integrals {
            topologies
                .0
                .push(Topology::Unknown(Integral::new(0, None, None)?));
        }
        Ok(topologies)
    }

    fn match_topologies_to_user_input(
        &self,
        input: AtomView,
    ) -> Result<Option<ReplacementRules>, VakintError> {
        for topology in self.0.iter() {
            if let Some(replacement_rules) = topology.match_topology_to_user_input(input)? {
                return Ok(Some(replacement_rules));
            }
        }
        Ok(None)
    }
}

#[allow(unused)]
#[derive(Debug, Clone)]
enum Topology {
    OneLoop(Integral),
    TwoLoop(Integral),
    ThreeLoop(Integral),
    FourLoop(Integral),
    Unknown(Integral),
}

impl fmt::Display for Topology {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Topology::OneLoop(i) => write!(f, "{}{}{}", "OneLoop(".magenta(), i, ")".magenta()),
            Topology::TwoLoop(i) => write!(f, "{}{}{}", "TwoLoop(".magenta(), i, ")".magenta()),
            Topology::ThreeLoop(i) => write!(f, "{}{}{}", "ThreeLoop(".magenta(), i, ")".magenta()),
            Topology::FourLoop(i) => write!(f, "{}{}{}", "FourLoop(".magenta(), i, ")".magenta()),
            Topology::Unknown(i) => write!(f, "{}{}{}", "Unknown(".red(), i, ")".red()),
        }
    }
}

impl From<Integral> for Topology {
    fn from(i: Integral) -> Self {
        match i.n_loops {
            1 => Topology::OneLoop(i),
            2 => Topology::TwoLoop(i),
            3 => Topology::ThreeLoop(i),
            4 => Topology::FourLoop(i),
            _ => Topology::Unknown(i),
        }
    }
}

impl Topology {
    fn get_integral(&self) -> &Integral {
        match self {
            Topology::OneLoop(i)
            | Topology::TwoLoop(i)
            | Topology::ThreeLoop(i)
            | Topology::FourLoop(i)
            | Topology::Unknown(i) => i,
        }
    }

    fn to_canonical(
        &self,
        integral: AtomView,
        replacement_rules: &ReplacementRules,
        short_form: bool,
    ) -> Atom {
        match self {
            Topology::Unknown(_) => Pattern::parse("topo(props_)").unwrap().replace_all(
                integral,
                &Pattern::parse("topo(UNKNOWN(props_))").unwrap(),
                None,
                None,
            ),
            t => t.get_integral().to_canonical(replacement_rules, short_form),
        }
    }

    fn match_topology_to_user_input(
        &self,
        input: AtomView,
    ) -> Result<Option<ReplacementRules>, VakintError> {
        match self {
            Topology::Unknown(_) => {
                // println!(
                //     "\n>>>Trying to match with the unknown integral pattern {}\n<<<",
                //     self
                // );
                let undirected_input = Pattern::parse("edge(x_,y_)").unwrap().replace_all(
                    input,
                    &Pattern::parse("uedge(x_,y_)").unwrap(),
                    None,
                    None,
                );
                let unknown_integral = self.get_integral();
                if unknown_integral
                    .generic_pattern
                    .pattern
                    .pattern_match(
                        undirected_input.as_view(),
                        &unknown_integral.generic_pattern.conditions,
                        &unknown_integral.generic_pattern.match_settings,
                    )
                    .next()
                    .is_some()
                {
                    //println!("Found a match!");
                    Ok(Some(ReplacementRules::default()))
                } else {
                    //println!("Does not match!");
                    Ok(None)
                }
            }
            t => {
                // println!("\n>>>Trying to match: {}\n<<<", self);
                if let Some(mut replacement_rule) =
                    t.get_integral().match_integral_to_user_input(input)?
                {
                    replacement_rule.canonical_topology = t.clone();
                    // println!("Found a match! ->\n{}", replacement_rule);
                    Ok(Some(replacement_rule))
                } else {
                    Ok(None)
                }
            }
        }
    }
}

#[derive(Debug, Clone)]
struct Integral {
    n_loops: usize,
    n_props: usize,
    name: String,
    generic_pattern: FullPattern,
    canonical_expression: Option<Atom>,
    short_expression: Option<Atom>,
    short_expression_pattern: Option<Pattern>,
}

impl fmt::Display for Integral {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "name='{}', n_loops={}, n_props_top_topo={}\n   | canonical_expression='{}',\n   | canonical_pattern='{}',\n   | short_expression='{}'",
            self.name,
            self.n_loops,
            self.n_props,
            self.canonical_expression
                .as_ref()
                .map(|e| format!("{}", e))
                .unwrap_or("N/A".into()).blue(),
            format!("{}",self.generic_pattern.pattern.to_atom().unwrap()).blue(),
            self.short_expression
                .as_ref()
                .map(|e| format!("{}", e))
                .unwrap_or("N/A".into()).green()
        )
    }
}

fn get_integer_from_match(m: &Match<'_>) -> Option<i64> {
    if let Match::Single(AtomView::Num(a)) = m {
        match a.get_coeff_view() {
            CoefficientView::Natural(n, 1) => Some(n),
            _ => None,
        }
    } else {
        None
    }
}

fn get_node_ids(
    match_stack: &HashMap<symbolica::atom::Symbol, symbolica::id::Match>,
) -> Result<(usize, usize), VakintError> {
    let id_node_left = if let Some(id) =
        get_integer_from_match(match_stack.get(&State::get_symbol("nl_")).unwrap())
    {
        id
    } else {
        return Err(VakintError::InvalidGenericExpression(format!(
            "Left node must be an integer: {}",
            match_stack.get(&State::get_symbol("nl_")).unwrap()
        )));
    };
    let id_node_right = if let Some(id) =
        get_integer_from_match(match_stack.get(&State::get_symbol("nr_")).unwrap())
    {
        id
    } else {
        return Err(VakintError::InvalidGenericExpression(format!(
            "Right node must be an integer: {}",
            match_stack.get(&State::get_symbol("nl_")).unwrap()
        )));
    };

    Ok((id_node_left as usize, id_node_right as usize))
}

#[allow(clippy::type_complexity)]
fn get_individual_momenta(momentum: &Match) -> Result<Vec<(Symbol, (Atom, i64))>, VakintError> {
    let momentum_atom =
        match momentum {
            Match::Single(a) => a.to_owned(),
            Match::Multiple(SliceType::Add, atoms) => Atom::parse(
                atoms
                    .iter()
                    .map(|a| a.to_string())
                    .collect::<Vec<_>>()
                    .join("+")
                    .as_str(),
            )
            .unwrap(),
            _ => return Err(VakintError::InvalidMomentumExpression(
                "Invalid expression for edge momentum. It is not in the form '<symbol>_(<int>)'."
                    .into(),
            )),
        };

    let err = Err(VakintError::InvalidMomentumExpression(format!(
        "Edge momentum {} does not contain only <symbol>_(<int>).",
        momentum_atom
    )));
    let mut res = vec![];
    let mut test = momentum_atom.to_owned();

    while let Some(m) = Pattern::parse("q_(lmbID_)")
        .unwrap()
        .pattern_match(
            test.as_view(),
            &(Condition::from((State::get_symbol("q_"), symbol_condition()))
                & Condition::from((State::get_symbol("lmbID_"), number_condition()))),
            &MatchSettings::default(),
        )
        .next()
    {
        let a_match = m.match_stack.get(State::get_symbol("lmbID_")).unwrap();
        let (mom_symbol, (atom_id, mom_id)) = (
            if let Match::FunctionName(s) = m.match_stack.get(State::get_symbol("q_")).unwrap() {
                *s
            } else {
                return err;
            },
            if let Match::Single(a) = a_match {
                if let Some(id) = get_integer_from_match(a_match) {
                    (a.to_owned(), id)
                } else {
                    return err;
                }
            } else {
                return err;
            },
        );
        test = fun!(mom_symbol, &atom_id).into_pattern().replace_all(
            test.as_view(),
            &Atom::Zero.into_pattern(),
            None,
            None,
        );

        res.push((mom_symbol, (atom_id, mom_id)));
    }

    if !test.is_zero() {
        return err;
    }

    Ok(res)
}

fn get_prop_with_id(
    expression: AtomView,
    prop_id: usize,
) -> Option<HashMap<symbolica::atom::Symbol, symbolica::id::Match>> {
    let no_settings = MatchSettings::default();
    let pattern =
        Pattern::parse(format!("prop({},edge(nl_,nr_),q_,mUVsq_,pow_)", prop_id).as_str()).unwrap();
    let number_node_condition = Condition::from((State::get_symbol("nl_"), gt_condition(0)))
        & Condition::from((State::get_symbol("nr_"), gt_condition(0)));
    if let Some(m) = pattern
        .pattern_match(expression, &number_node_condition, &no_settings)
        .next()
    {
        Some(
            m.match_stack
                .get_matches()
                .iter()
                .cloned()
                .collect::<HashMap<_, _>>(),
        )
    } else {
        None
    }
}

impl Integral {
    fn new(
        // This refers to the *total* number of propagators in the top-level topology,
        // i.e. the number of entries in the corresponding short_expression
        tot_n_props: usize,
        canonical_expression: Option<Atom>,
        short_expression: Option<Atom>,
    ) -> Result<Integral, VakintError> {
        if canonical_expression.is_none() {
            // This is an unknown topology
            let all_accepting_pattern = FullPattern {
                pattern: Pattern::parse("topo(props_)").unwrap(),
                conditions: Condition::from((State::get_symbol("props_"), propagators_condition())),
                match_settings: MatchSettings::default(),
            };
            let short_expression_pattern = short_expression.as_ref().map(|a| a.into_pattern());
            return Ok(Integral {
                name: "Unknown".into(),
                n_loops: 0,
                n_props: 0,
                generic_pattern: all_accepting_pattern,
                canonical_expression,
                short_expression,
                short_expression_pattern,
            });
        }
        let e = canonical_expression.clone().unwrap();

        let mut loop_mom_indices = HashSet::<i64>::default();
        let mut next_atom = e.clone();
        let mut old_atom = e.clone();

        let mut generic_expression = Atom::parse("1").unwrap();
        let mut generic_condition = Condition::default();

        for i_prop in 1..=tot_n_props {
            if let Some(m) = get_prop_with_id(e.as_view(), i_prop) {
                // Format check for the momenta
                let momenta = get_individual_momenta(m.get(&State::get_symbol("q_")).unwrap())?;

                let mass_symbol_string = if let Match::Single(a) =
                    m.get(&State::get_symbol("mUVsq_")).unwrap()
                {
                    if let Some(m3) = Pattern::parse("msq(mid_)")
                        .unwrap()
                        .pattern_match(
                            a.as_atom_view(),
                            &Condition::from((
                                State::get_symbol("mid_"),
                                range_condition(1, tot_n_props as i64),
                            )),
                            &MatchSettings::default(),
                        )
                        .next()
                    {
                        format!(
                            "msq{}_",
                            get_integer_from_match(
                                m3.match_stack.get(State::get_symbol("mid_")).unwrap(),
                            )
                            .unwrap()
                        )
                    } else {
                        return Err(VakintError::InvalidGenericExpression(
                            format!("Generic expression does not have masses formatted as msq(integer in [1,n_props]): {}",a),
                        ));
                    }
                } else {
                    return Err(VakintError::InvalidGenericExpression(
                        "Generic expression does not have masses formatted as msq(<m_integer_id>)."
                            .into(),
                    ));
                };

                // If the power is a pattern, then match it
                let power_match = m.get(&State::get_symbol("pow_")).unwrap();
                let pow_symbol_string = if let Some(pow) = get_integer_from_match(power_match) {
                    format!("{}", pow)
                } else if let Match::Single(a) = power_match {
                    if let Some(m3) = Pattern::parse("pow(pid_)")
                        .unwrap()
                        .pattern_match(
                            a.as_atom_view(),
                            &Condition::from((
                                State::get_symbol("pid_"),
                                range_condition(1, tot_n_props as i64),
                            )),
                            &MatchSettings::default(),
                        )
                        .next()
                    {
                        format!(
                            "pow{}_",
                            get_integer_from_match(
                                m3.match_stack.get(State::get_symbol("pid_")).unwrap(),
                            )
                            .unwrap()
                        )
                    } else {
                        return Err(VakintError::InvalidGenericExpression(
                            format!("Generic expression does not have powers formatted as pow(integer in [1,n_props]): {}",a),
                        ));
                    }
                } else {
                    return Err(VakintError::InvalidGenericExpression(
                        "Generic expression does not have powers formatted as pow(<m_integer_id>)."
                            .into(),
                    ));
                };

                for (mom_symbol, (_mom_atom_id, _mom_id)) in momenta {
                    if mom_symbol != State::get_symbol("k") {
                        return Err(VakintError::InvalidGenericExpression(
                            "Generic expression does not have momenta involving only expressions of the type k(<integer>)."
                                .to_string(),
                        ));
                    }
                    loop_mom_indices.insert(_mom_id);
                }

                let (id_node_left, id_node_right) = get_node_ids(&m)?;

                generic_expression = generic_expression * Atom::parse(
                    format!("prop(id{id}_,uedge(n{id_l_node}l_,n{id_r_node}r_),q{id}__,{mass},{power})",
                    id=i_prop,id_l_node=id_node_left, id_r_node=id_node_right,
                    mass=mass_symbol_string, power=pow_symbol_string).as_str()).unwrap();
                generic_condition = generic_condition
                    & Condition::from((State::get_symbol(mass_symbol_string), symbol_or_number()))
                    & Condition::from((
                        State::get_symbol(format!("pow{}_", i_prop).as_str()),
                        number_condition(),
                    ))
                    & Condition::from((
                        State::get_symbol(format!("id{}_", i_prop).as_str()),
                        number_condition(),
                    ))
                    & apply_restriction_to_symbols(
                        vec![
                            State::get_symbol(format!("n{}l_", id_node_left).as_str()),
                            State::get_symbol(format!("n{}r_", id_node_right).as_str()),
                        ],
                        &symbol_or_number(),
                    );
                next_atom = Pattern::parse(format!("prop({},args__)", i_prop).as_str())
                    .unwrap()
                    .replace_all(
                        old_atom.as_view(),
                        &Pattern::parse("1").unwrap(),
                        None,
                        None,
                    );
                old_atom = next_atom.clone();
            }
        }

        if next_atom != Atom::parse("topo(1)").unwrap() {
            return Err(VakintError::InvalidGenericExpression(format!(
                "Not all propagators of the generic expression supplied have been successfully identified. Left-over: {}",next_atom)
            ));
        }
        generic_expression = fun!(State::get_symbol("topo"), &generic_expression);
        let generic_pattern = FullPattern {
            pattern: generic_expression.into_pattern(),
            conditions: generic_condition,
            match_settings: MatchSettings::default(),
        };

        let name = if let Some(m) = Pattern::parse("f_(args__)")
            .unwrap()
            .pattern_match(
                short_expression.as_ref().unwrap().as_view(),
                &Condition::default(),
                &MatchSettings::default(),
            )
            .next()
        {
            if let Some(Match::FunctionName(s)) = m.match_stack.get(State::get_symbol("f_")) {
                s.to_string()
            } else {
                return Err(VakintError::InvalidShortExpression(
                    short_expression.unwrap().to_string(),
                ));
            }
        } else {
            return Err(VakintError::InvalidShortExpression(
                short_expression.unwrap().to_string(),
            ));
        };

        let mut short_expression_pattern = short_expression.clone().unwrap();
        for i_prop in 1..=tot_n_props {
            short_expression_pattern = Pattern::parse(format!("pow({})", i_prop).as_str())
                .unwrap()
                .replace_all(
                    short_expression_pattern.as_view(),
                    &Pattern::parse(format!("pow{}_", i_prop).as_str()).unwrap(),
                    Some(&Condition::default()),
                    Some(&MatchSettings {
                        allow_new_wildcards_on_rhs: true,
                        ..MatchSettings::default()
                    }),
                );
            short_expression_pattern = Pattern::parse(format!("msq({})", i_prop).as_str())
                .unwrap()
                .replace_all(
                    short_expression_pattern.as_view(),
                    &Pattern::parse(format!("msq{}_", i_prop).as_str()).unwrap(),
                    Some(&Condition::default()),
                    Some(&MatchSettings {
                        allow_new_wildcards_on_rhs: true,
                        ..MatchSettings::default()
                    }),
                );
        }

        Ok(Integral {
            name,
            n_loops: loop_mom_indices.len(),
            n_props: tot_n_props,
            generic_pattern,
            canonical_expression,
            short_expression,
            short_expression_pattern: Some(short_expression_pattern.into_pattern()),
        })
    }

    fn match_integral_to_short_user_input(
        &self,
        input: AtomView,
    ) -> Result<Option<ReplacementRules>, VakintError> {
        let unwrapped_input = Pattern::parse("topo(integral_)").unwrap().replace_all(
            input,
            &Pattern::parse("integral_").unwrap(),
            None,
            None,
        );

        if let Some(short_expression_pattern) = self.short_expression_pattern.as_ref() {
            if let Some(m1) = short_expression_pattern
                .pattern_match(
                    unwrapped_input.as_view(),
                    &apply_restriction_to_symbols(
                        (1..=self.n_props)
                            .map(|i_prop| State::get_symbol(format!("pow{}_", i_prop)))
                            .collect(),
                        &symbol_or_number(),
                    ),
                    &MatchSettings::default(),
                )
                .next()
            {
                let mut replacement_rules = ReplacementRules::default();
                for i_prop in 1..=self.n_props {
                    if let Some(Match::Single(a)) = m1
                        .match_stack
                        .get(State::get_symbol(format!("pow{}_", i_prop).as_str()))
                    {
                        replacement_rules.canonical_expression_substitutions.insert(
                            Atom::parse(format!("pow({})", i_prop).as_str()).unwrap(),
                            a.to_owned(),
                        );
                    }
                }
                // Dummy substitutions for the numerator in this case
                for i_loop in 1..=self.n_loops {
                    replacement_rules.numerator_substitutions.insert(
                        Atom::parse(format!("k({},idx_)", i_loop).as_str()).unwrap(),
                        Atom::parse(format!("k({},idx_)", i_loop).as_str()).unwrap(),
                    );
                }
                Ok(Some(replacement_rules))
            } else {
                // Make sure the user did not intend to pass a short form expression
                if let Some(m) = Pattern::parse("fn_(args__)")
                    .unwrap()
                    .pattern_match(input, &Condition::default(), &MatchSettings::default())
                    .next()
                {
                    if let Match::FunctionName(s) =
                        m.match_stack.get(State::get_symbol("fn_")).unwrap()
                    {
                        if *s == State::get_symbol(self.name.clone()) {
                            Err(VakintError::InvalidShortExpression(format!("{}", input)))
                        } else {
                            Ok(None)
                        }
                    } else {
                        Ok(None)
                    }
                } else {
                    Ok(None)
                }
            }
        } else {
            Ok(None)
        }
    }

    fn match_integral_to_user_input(
        &self,
        input: AtomView,
    ) -> Result<Option<ReplacementRules>, VakintError> {
        // Check if the input is a short expression
        if let Ok(Some(a)) = self.match_integral_to_short_user_input(input) {
            return Ok(Some(a));
        }

        let undirected_input = Pattern::parse("edge(x_,y_)").unwrap().replace_all(
            input,
            &Pattern::parse("uedge(x_,y_)").unwrap(),
            None,
            None,
        );

        // println!("input: {}", undirected_input);
        // println!(
        //     "pattern: {}",
        //     self.generic_pattern.pattern.to_atom().unwrap()
        // );
        if let Some(m1) = self
            .generic_pattern
            .pattern
            .pattern_match(
                undirected_input.as_view(),
                &self.generic_pattern.conditions,
                &self.generic_pattern.match_settings,
            )
            .next()
        {
            let mut replacement_rules = ReplacementRules::default();

            for prop_id in 1..=self.n_props {
                if let Some(canonical_prop_match) = get_prop_with_id(
                    self.canonical_expression.as_ref().unwrap().as_view(),
                    prop_id,
                ) {
                    if let Some(Match::Single(a)) = m1
                        .match_stack
                        .get(State::get_symbol(format!("pow{}_", prop_id).as_str()))
                    {
                        replacement_rules.canonical_expression_substitutions.insert(
                            Atom::parse(format!("pow({})", prop_id).as_str()).unwrap(),
                            a.to_owned(),
                        );
                    }

                    if let Some(Match::Single(a)) = m1
                        .match_stack
                        .get(State::get_symbol(format!("msq{}_", prop_id)))
                    {
                        replacement_rules.canonical_expression_substitutions.insert(
                            Atom::parse(format!("msq({})", prop_id).as_str()).unwrap(),
                            a.to_owned(),
                        );
                    }

                    let input_prop_id = if let Some(i) = get_integer_from_match(
                        m1.match_stack
                            .get(State::get_symbol(format!("id{}_", prop_id).as_str()))
                            .unwrap(),
                    ) {
                        i as usize
                    } else {
                        panic!("Match id from input expression must be an integer.")
                    };

                    replacement_rules
                        .edge_ids_canonical_to_input_map
                        .insert(prop_id, input_prop_id);

                    // get the node ids in user's input for this prop as well as the one in the canonical expression
                    let input_prop_match = get_prop_with_id(input, input_prop_id).unwrap();
                    let (input_id_node_left, input_id_node_right) =
                        get_node_ids(&input_prop_match).unwrap();

                    let canonical_ids = get_node_ids(&canonical_prop_match).unwrap();

                    let (canonical_id_node_left, canonical_id_node_right) = (
                        get_integer_from_match(
                            m1.match_stack
                                .get(State::get_symbol(
                                    format!("n{}l_", canonical_ids.0).as_str(),
                                ))
                                .unwrap(),
                        )
                        .unwrap() as usize,
                        get_integer_from_match(
                            m1.match_stack
                                .get(State::get_symbol(
                                    format!("n{}r_", canonical_ids.1).as_str(),
                                ))
                                .unwrap(),
                        )
                        .unwrap() as usize,
                    );

                    let is_edge_flipped = if (input_id_node_left, input_id_node_right)
                        == (canonical_id_node_left, canonical_id_node_right)
                    {
                        false
                    } else if (input_id_node_right, input_id_node_left)
                        == (canonical_id_node_left, canonical_id_node_right)
                    {
                        true
                    } else {
                        unreachable!(
                            "Nodes IDs should have been matches: ({},{})!=({},{})=(n{}l_,n{}r_)",
                            input_id_node_left,
                            input_id_node_right,
                            canonical_id_node_left,
                            canonical_id_node_right,
                            canonical_ids.0,
                            canonical_ids.1
                        )
                    };

                    let potential_lmb_match =
                        input_prop_match.get(&State::get_symbol("q_")).unwrap();
                    // This is an LMB that'll need replacement in the numerator
                    if let Match::Single(AtomView::Fun(_)) = potential_lmb_match {
                        let (input_momentum_symbol, (input_momentum_atom_id, _input_momentum_id)) =
                            get_individual_momenta(potential_lmb_match)?.pop().unwrap();
                        let input_lmb_pattern = fun!(
                            input_momentum_symbol,
                            &input_momentum_atom_id,
                            &Atom::parse("idx_").unwrap()
                        );

                        let mut canonical_momenta_atom_for_pattern =
                            if let Match::Single(canonical_momenta) =
                                canonical_prop_match.get(&State::get_symbol("q_")).unwrap()
                            {
                                Pattern::parse("k(ilmb_)").unwrap().replace_all(
                                    *canonical_momenta,
                                    &Pattern::parse("k(ilmb_,idx_)").unwrap(),
                                    Some(&Condition::from((
                                        State::get_symbol("ilmb_"),
                                        number_condition(),
                                    ))),
                                    Some(&MatchSettings {
                                        allow_new_wildcards_on_rhs: true,
                                        ..MatchSettings::default()
                                    }),
                                )
                            } else {
                                unreachable!()
                            };

                        if is_edge_flipped {
                            canonical_momenta_atom_for_pattern =
                                canonical_momenta_atom_for_pattern * -1
                        }

                        replacement_rules
                            .numerator_substitutions
                            .insert(input_lmb_pattern, canonical_momenta_atom_for_pattern);
                    }
                } else {
                    replacement_rules.canonical_expression_substitutions.insert(
                        Atom::parse(format!("pow({})", prop_id).as_str()).unwrap(),
                        Atom::Zero,
                    );
                }
            }
            Ok(Some(replacement_rules))
        } else {
            Ok(None)
        }
    }

    fn to_canonical(&self, replacement_rules: &ReplacementRules, short_form: bool) -> Atom {
        let mut new_expression = if short_form {
            fun!(
                State::get_symbol("topo"),
                self.short_expression.as_ref().unwrap().clone()
            )
        } else {
            self.canonical_expression.as_ref().unwrap().clone()
        };
        for (source, target) in replacement_rules.canonical_expression_substitutions.iter() {
            new_expression = source.into_pattern().replace_all(
                new_expression.as_view(),
                &target.into_pattern(),
                None,
                None,
            );
        }

        // Remove propagators with zero powers
        new_expression = Pattern::parse("prop(propID_,edge(nl_,nr_),q_,mUVsq_,0)")
            .unwrap()
            .replace_all(
                new_expression.as_view(),
                &Pattern::parse("1").unwrap(),
                None,
                None,
            );

        new_expression
    }
}

pub struct VakintSettings {
    #[allow(unused)]
    pub use_pysecdec: bool,
    pub external_momenta_symbol: String,
    pub verify_numerator_identification: bool,
    pub allow_unknown_integrals: bool,
}

#[allow(clippy::derivable_impls)]
impl Default for VakintSettings {
    fn default() -> Self {
        VakintSettings {
            use_pysecdec: false,
            external_momenta_symbol: "p".into(),
            verify_numerator_identification: true,
            allow_unknown_integrals: true,
        }
    }
}

#[derive(Debug, Clone)]
struct FullPattern {
    pattern: Pattern,
    conditions: Condition<WildcardAndRestriction>,
    match_settings: MatchSettings,
}

impl From<Pattern> for FullPattern {
    fn from(pattern: Pattern) -> FullPattern {
        FullPattern {
            pattern,
            conditions: Condition::default(),
            match_settings: MatchSettings::default(),
        }
    }
}

pub struct Vakint {
    #[allow(unused)]
    pub settings: VakintSettings,
    pub topologies: Topologies,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VakintTerm {
    pub integral: Atom,
    pub numerator: Atom,
}

impl From<VakintTerm> for Atom {
    fn from(vakint_term: VakintTerm) -> Atom {
        vakint_term.integral * vakint_term.numerator
    }
}

impl fmt::Display for VakintTerm {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "({}) x {}",
            format!("{}", self.numerator).cyan(),
            format!("{}", self.integral).green(),
        )
    }
}

impl VakintTerm {
    pub fn canonicalize(
        &mut self,
        vakint: &Vakint,
        replacement_rules: &ReplacementRules,
        short_form: bool,
    ) -> Result<(), VakintError> {
        self.integral = replacement_rules.canonical_topology.to_canonical(
            self.integral.as_view(),
            replacement_rules,
            short_form,
        );
        if matches!(replacement_rules.canonical_topology, Topology::Unknown(_),) {
            return Ok(()); // No further canonicalization possible
        }

        let mut new_numerator = self.numerator.clone();
        let mut test = self.numerator.clone();
        for (source, target) in replacement_rules.numerator_substitutions.iter() {
            new_numerator = source.into_pattern().replace_all(
                new_numerator.as_view(),
                &target.into_pattern(),
                None,
                None,
            );
            test = source.into_pattern().replace_all(
                test.as_view(),
                &Atom::parse("1").unwrap().into_pattern(),
                None,
                None,
            );
        }

        // Make sure to also set all externals to zero for the test
        test = Pattern::parse(
            format!("{}(momID_,idx_)", vakint.settings.external_momenta_symbol).as_str(),
        )
        .unwrap()
        .replace_all(
            test.as_view(),
            &Atom::parse("1").unwrap().into_pattern(),
            None,
            None,
        );

        if vakint.settings.verify_numerator_identification && !matches!(test, Atom::Num(_)) {
            return Err(VakintError::NumeratorNotReplaced(
                vakint.settings.external_momenta_symbol.clone(),
                test.to_string(),
            ));
        }
        self.numerator = new_numerator;

        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VakintExpression(Vec<VakintTerm>);

impl fmt::Display for VakintExpression {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}{}",
            if self.0.len() > 1 { " " } else { "" },
            self.0
                .iter()
                .map(|term| format!("{}", term))
                .collect::<Vec<_>>()
                .join("\n+")
        )
    }
}

impl VakintExpression {
    pub fn split_integrals(input: AtomView) -> Result<Vec<VakintTerm>, VakintError> {
        let (integrals, remainder) = input.to_owned().coefficient_list(State::get_symbol("topo"));

        if remainder != Atom::parse("0").unwrap() {
            return Err(VakintError::InvalidIntegralFormat(format!("{}", remainder)));
        }

        let mut res = vec![];
        for (integral, numerator) in integrals.iter() {
            res.push(VakintTerm {
                integral: integral.clone(),
                numerator: numerator.clone(),
            });
        }
        Ok(res)
    }

    pub fn canonicalize(
        &mut self,
        vakint: &Vakint,
        topologies: &Topologies,
        short_form: bool,
    ) -> Result<(), VakintError> {
        for term in self.0.iter_mut() {
            if let Some(replacement_rules) =
                topologies.match_topologies_to_user_input(term.integral.as_view())?
            {
                //println!("replacement_rules = {}", replacement_rules,);
                term.canonicalize(vakint, &replacement_rules, short_form)?;
            } else {
                return Err(VakintError::UnreckognizedIntegral(
                    term.integral.to_string(),
                ));
            }
        }
        Ok(())
    }

    #[allow(dead_code)]
    pub fn map<F>(&mut self, f: F) -> VakintExpression
    where
        F: Fn(&VakintTerm) -> VakintTerm,
    {
        VakintExpression(self.0.iter().map(f).collect())
    }

    #[allow(dead_code)]
    pub fn map_numerator<F>(&self, f: F) -> VakintExpression
    where
        F: Fn(AtomView) -> Atom,
    {
        VakintExpression(
            self.0
                .iter()
                .map(|term| VakintTerm {
                    integral: term.integral.clone(),
                    numerator: f(term.numerator.as_view()),
                })
                .collect(),
        )
    }

    #[allow(dead_code)]
    pub fn map_integrals<F>(&self, f: F) -> VakintExpression
    where
        F: Fn(AtomView) -> Atom,
    {
        VakintExpression(
            self.0
                .iter()
                .map(|term| VakintTerm {
                    integral: f(term.integral.as_view()),
                    numerator: term.numerator.clone(),
                })
                .collect(),
        )
    }
}

impl TryFrom<Atom> for VakintExpression {
    type Error = VakintError;
    fn try_from(atom: Atom) -> Result<Self, VakintError> {
        Ok(VakintExpression(VakintExpression::split_integrals(
            atom.as_view(),
        )?))
    }
}

impl TryFrom<AtomView<'_>> for VakintExpression {
    type Error = VakintError;
    fn try_from(atom_view: AtomView) -> Result<Self, VakintError> {
        Ok(VakintExpression(VakintExpression::split_integrals(
            atom_view,
        )?))
    }
}

impl From<VakintExpression> for Atom {
    fn from(vakint_expr: VakintExpression) -> Atom {
        let mut res = Atom::Zero;
        for term in vakint_expr.0.iter() {
            let t: Atom = VakintTerm::into(term.clone());
            res = res + t;
        }
        res
    }
}

impl Vakint {
    fn initialize_symbolica_symbols() {
        _ = State::get_symbol_with_attributes("uedge", &[FunctionAttribute::Symmetric]);
    }
    pub fn new(settings: Option<VakintSettings>) -> Result<Self, VakintError> {
        Vakint::initialize_symbolica_symbols();
        let vakint_settings = settings.unwrap_or_default();
        let topologies = Topologies::generate_topologies(&vakint_settings)?;
        let vakint = Vakint {
            settings: vakint_settings,
            topologies,
        };
        Ok(vakint)
    }

    pub fn evaluate(&self, input: AtomView) -> Result<Atom, VakintError> {
        let mut vakint_expr = VakintExpression::try_from(input)?;
        vakint_expr.canonicalize(self, &self.topologies, false)?;
        Ok(vakint_expr.into())
    }

    pub fn to_canonical(&self, input: AtomView, short_form: bool) -> Result<Atom, VakintError> {
        let mut vakint_expr = VakintExpression::try_from(input)?;
        vakint_expr.canonicalize(self, &self.topologies, short_form)?;
        Ok(vakint_expr.into())
    }
}
