use std::fmt;

use colored::Colorize;
use symbolica::{
    atom::{Atom, AtomView},
    id::Pattern,
};

use crate::{
    get_node_ids, get_prop_with_id, symbols::S, Integral, ReplacementRules, VakintError,
    VakintSettings,
};

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
    pub fn generate_topologies(settings: &VakintSettings) -> Result<Self, VakintError> {
        // One-loop topology
        let mut topologies = Topologies(vec![Integral::new(
            1,
            Some(Atom::parse("topo(prop(1,edge(1,1),k(1),msq(1),pow(1)))").unwrap()),
            Some(Atom::parse("I1LA(msq(1),pow(1))").unwrap()),
            Some(Atom::parse("1").unwrap()),
            None,
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
                Some(Atom::parse("1").unwrap().as_view()),
                None,
            )?
            .0,
        );

        if settings.allow_unknown_integrals {
            topologies
                .0
                .push(Topology::Unknown(Integral::new(0, None, None, None, None)?));
        }
        Ok(topologies)
    }

    pub fn generate_topologies_with_contractions(
        n_tot_props: usize,
        canonical_expression: AtomView,
        short_expression: AtomView,
        contractions: Vec<Vec<usize>>,
        matad_expression: Option<AtomView>,
        fmft_expression: Option<AtomView>,
    ) -> Result<Self, VakintError> {
        let mut topologies = vec![];
        for contracted_prop_indices in contractions {
            let mut contracted_canonical_expression = canonical_expression.to_owned();
            let mut contracted_short_expression = short_expression.to_owned();
            let mut contracted_matad_expression = if let Some(av) = matad_expression {
                av.to_owned()
            } else {
                Atom::Zero
            };
            let mut contracted_fmft_expression = if let Some(av) = fmft_expression {
                av.to_owned()
            } else {
                Atom::Zero
            };
            let mut nodes_to_merge = vec![];
            for prop_id in contracted_prop_indices.iter() {
                if let Some(m) =
                    get_prop_with_id(contracted_canonical_expression.as_view(), *prop_id)
                {
                    let (left_node_id, right_node_id) = get_node_ids(&m).unwrap();
                    nodes_to_merge.push((right_node_id, left_node_id));
                } else {
                    return Err(VakintError::InvalidIntegralFormat(format!(
                        "Cannot contract propagatore id {} as it is not found in integral: {}",
                        prop_id, contracted_canonical_expression
                    )));
                }
                contracted_canonical_expression =
                    Pattern::parse(format!("prop({},args__)", prop_id).as_str())
                        .unwrap()
                        .replace_all(
                            contracted_canonical_expression.as_view(),
                            &S.one.into_pattern().into(),
                            None,
                            None,
                        );

                for a in [
                    &mut contracted_short_expression,
                    &mut contracted_matad_expression,
                    &mut contracted_fmft_expression,
                ] {
                    if !a.is_zero() {
                        *a = Pattern::parse(format!("pow({})", prop_id).as_str())
                            .unwrap()
                            .replace_all(a.as_view(), &S.zero.into_pattern().into(), None, None);
                    }
                }
            }
            let mut old_contracted_canonical_expression = contracted_canonical_expression.clone();
            'replace_contracted_nodes: loop {
                for (old_node_id, new_node_id) in nodes_to_merge.iter() {
                    for (lhs, rhs) in [
                        (
                            format!("edge({},nr_)", old_node_id),
                            format!("edge({},nr_)", new_node_id),
                        ),
                        (
                            format!("edge(nl_,{})", old_node_id),
                            format!("edge(nl_,{})", new_node_id),
                        ),
                    ] {
                        contracted_canonical_expression =
                            Pattern::parse(lhs.as_str()).unwrap().replace_all(
                                contracted_canonical_expression.as_view(),
                                &Pattern::parse(rhs.as_str()).unwrap().into(),
                                None,
                                None,
                            );
                    }
                }
                if old_contracted_canonical_expression == contracted_canonical_expression {
                    break 'replace_contracted_nodes;
                }
                old_contracted_canonical_expression = contracted_canonical_expression.clone();
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
                    if contracted_matad_expression.is_zero() {
                        None
                    } else {
                        Some(contracted_matad_expression)
                    },
                    if contracted_fmft_expression.is_zero() {
                        None
                    } else {
                        Some(contracted_fmft_expression)
                    },
                )?
                .into(),
            );
        }
        Ok(Topologies(topologies))
    }

    pub fn match_topologies_to_user_input(
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
pub enum Topology {
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
    pub fn get_integral(&self) -> &Integral {
        match self {
            Topology::OneLoop(i)
            | Topology::TwoLoop(i)
            | Topology::ThreeLoop(i)
            | Topology::FourLoop(i)
            | Topology::Unknown(i) => i,
        }
    }

    pub fn get_integral_mut(&mut self) -> &mut Integral {
        match self {
            Topology::OneLoop(i)
            | Topology::TwoLoop(i)
            | Topology::ThreeLoop(i)
            | Topology::FourLoop(i)
            | Topology::Unknown(i) => i,
        }
    }

    pub fn to_canonical(
        &self,
        integral: AtomView,
        replacement_rules: &ReplacementRules,
        short_form: bool,
    ) -> Atom {
        match self {
            Topology::Unknown(_) => Pattern::parse("topo(props_)").unwrap().replace_all(
                integral,
                &Pattern::parse("topo(UNKNOWN(props_))").unwrap().into(),
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
                    &Pattern::parse("uedge(x_,y_)").unwrap().into(),
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
