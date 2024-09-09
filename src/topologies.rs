use std::{collections::HashSet, fmt};

use colored::Colorize;
use symbolica::{
    atom::{Atom, AtomView},
    id::{Condition, MatchSettings, Pattern},
    state::State,
};

use crate::{
    get_individual_momenta, get_node_ids, get_prop_with_id, symbols::S, Integral, ReplacementRules,
    VakintError, VakintSettings,
};
pub enum TopologyContractions {
    Custom(Vec<Vec<usize>>),
    Automatic,
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
    pub fn generate_topologies(settings: &VakintSettings) -> Result<Self, VakintError> {
        // ==
        // One-loop topology
        // ==
        let mut topologies = Topologies(vec![Integral::new(
            1,
            Some(Atom::parse("topo(prop(1,edge(1,1),k(1),msq(1),pow(1)))").unwrap()),
            Some(Atom::parse("I1L(msq(1),pow(1))").unwrap()),
            Some(Atom::parse("1").unwrap()),
            None,
        )?
        .into()]);
        // ==
        // Two-loop topologies
        // ==
        topologies.0.extend(
            Topologies::generate_topologies_with_contractions(
                Atom::parse(
                    "topo(\
                         prop(1,edge(1,2),k(1),msq(1),pow(1))\
                        *prop(2,edge(1,2),k(2),msq(1),pow(2))\
                        *prop(3,edge(2,1),k(1)+k(2),msq(1),pow(3))\
                    )",
                )
                .unwrap()
                .as_view(),
                Atom::parse("I2L(msq(1),pow(1),pow(2),pow(3))")
                    .unwrap()
                    .as_view(),
                TopologyContractions::Custom(vec![vec![], vec![3]]),
                Some(Atom::parse("1").unwrap().as_view()),
                None,
            )?
            .0,
        );
        // ==
        // Three-loop topologies
        // ==
        topologies.0.extend(
            Topologies::generate_topologies_with_contractions(
                Atom::parse(
                    "topo(\
                         prop(1,edge(1,2),k(1),msq(1),pow(1))\
                        *prop(2,edge(2,3),k(2),msq(1),pow(2))\
                        *prop(3,edge(3,1),k(3),msq(1),pow(3))\
                        *prop(4,edge(1,4),k(3)-k(1),msq(1),pow(4))\
                        *prop(5,edge(2,4),k(1)-k(2),msq(1),pow(5))\
                        *prop(6,edge(3,4),k(2)-k(3),msq(1),pow(6))\
                    )",
                )
                .unwrap()
                .as_view(),
                Atom::parse("I3L(msq(1),pow(1),pow(2),pow(3),pow(4),pow(5),pow(6))")
                    .unwrap()
                    .as_view(),
                TopologyContractions::Automatic,
                Some(Atom::parse("1").unwrap().as_view()),
                None,
            )?
            .0,
        );
        // ==
        // Four-loop topologies
        // ==
        // H topology in FMFT
        // in FMFT paper, FMFT edge ids are these edge ids + 1, except for p_1 which is edge #9 here.
        topologies.0.extend(
            Topologies::generate_topologies_with_contractions(
                Atom::parse(
                    "topo(\
                         prop(1,edge(5,1),k(1),msq(1),pow(1))\
                        *prop(2,edge(2,6),k(2),msq(1),pow(2))\
                        *prop(3,edge(6,5),k(3),msq(1),pow(3))\
                        *prop(4,edge(4,3),k(4),msq(1),pow(4))\
                        *prop(5,edge(3,6),k(1)-k(3),msq(1),pow(5))\
                        *prop(6,edge(6,4),k(2)-k(3),msq(1),pow(6))\
                        *prop(7,edge(3,1),k(3)-k(1)+k(4),msq(1),pow(7))\
                        *prop(8,edge(2,4),k(3)-k(2)+k(4),msq(1),pow(8))\
                        *prop(9,edge(1,2),k(3)+k(4),msq(1),pow(9))\
                    )",
                )
                .unwrap()
                .as_view(),
                Atom::parse(
                    "I4L_H(msq(1),pow(1),pow(2),pow(3),pow(4),pow(5),pow(6),pow(7),pow(8),pow(9))",
                )
                .unwrap()
                .as_view(),
                TopologyContractions::Custom(vec![vec![]]),
                None,
                Some(Atom::parse("FMFT_H").unwrap().as_view()),
            )?
            .0,
        );
        // X topology in FMFT
        // in FMFT paper, FMFT edge ids are these edge ids + 1.
        topologies.0.extend(
            Topologies::generate_topologies_with_contractions(
                Atom::parse(
                    "topo(\
                         prop(1,edge(5,1),k(1),msq(1),pow(1))\
                        *prop(2,edge(2,6),k(2),msq(1),pow(2))\
                        *prop(3,edge(6,5),k(3),msq(1),pow(3))\
                        *prop(4,edge(4,3),k(4),msq(1),pow(4))\
                        *prop(5,edge(3,5),k(1)-k(3),msq(1),pow(5))\
                        *prop(6,edge(6,4),k(2)-k(3),msq(1),pow(6))\
                        *prop(7,edge(3,2),k(3)-k(1)+k(4),msq(1),pow(7))\
                        *prop(8,edge(1,4),k(3)-k(2)+k(4),msq(1),pow(8))\
                        *prop(9,edge(2,1),k(3)-k(1)-k(2)+k(4),msq(1),pow(9))\
                    )",
                )
                .unwrap()
                .as_view(),
                Atom::parse(
                    "I4L_X(msq(1),pow(1),pow(2),pow(3),pow(4),pow(5),pow(6),pow(7),pow(8),pow(9))",
                )
                .unwrap()
                .as_view(),
                TopologyContractions::Custom(vec![vec![]]),
                None,
                Some(Atom::parse("FMFT_X").unwrap().as_view()),
            )?
            .0,
        );
        // BMW topology in FMFT
        // When sorting both these edges and the FMT ones according to their ID, both lists match.
        topologies.0.extend(
            Topologies::generate_topologies_with_contractions(
                Atom::parse(
                    "topo(\
                         prop(1,edge(1,2),k(1),msq(1),pow(1))\
                        *prop(2,edge(2,5),k(2),msq(1),pow(2))\
                        *prop(3,edge(3,4),k(3),msq(1),pow(3))\
                        *prop(4,edge(4,5),k(4),msq(1),pow(4))\
                        *prop(5,edge(2,3),k(1)-k(2),msq(1),pow(5))\
                        *prop(6,edge(4,1),k(3)-k(4),msq(1),pow(6))\
                        *prop(7,edge(5,3),k(2)+k(3)-k(1),msq(1),pow(7))\
                        *prop(8,edge(1,5),k(3)-k(4)-k(1),msq(1),pow(8))\
                    )",
                )
                .unwrap()
                .as_view(),
                Atom::parse(
                    "I4L_BMW(msq(1),pow(1),pow(2),pow(3),pow(4),pow(5),pow(6),pow(7),pow(8))",
                )
                .unwrap()
                .as_view(),
                TopologyContractions::Custom(vec![vec![]]),
                None,
                Some(Atom::parse("FMFT_BMW").unwrap().as_view()),
            )?
            .0,
        );
        // FG topology in FMFT
        // When sorting both these edges and the FMT ones according to their ID, both lists match.
        topologies.0.extend(
            Topologies::generate_topologies_with_contractions(
                Atom::parse(
                    "topo(\
                         prop(1,edge(5,3),k(1),msq(1),pow(1))\
                        *prop(2,edge(3,4),k(2),msq(1),pow(2))\
                        *prop(3,edge(4,5),k(3),msq(1),pow(3))\
                        *prop(4,edge(2,1),k(1)-k(3),msq(1),pow(4))\
                        *prop(5,edge(5,1),k(4),msq(1),pow(5))\
                        *prop(6,edge(4,2),k(2)-k(3),msq(1),pow(6))\
                        *prop(7,edge(1,5),k(1)-k(3)+k(4),msq(1),pow(7))\
                        *prop(8,edge(3,2),k(1)-k(2),msq(1),pow(8))\
                    )",
                )
                .unwrap()
                .as_view(),
                Atom::parse(
                    "I4L_FG(msq(1),pow(1),pow(2),pow(3),pow(4),pow(5),pow(6),pow(7),pow(8))",
                )
                .unwrap()
                .as_view(),
                TopologyContractions::Automatic,
                None,
                Some(Atom::parse("FMFT_FG").unwrap().as_view()),
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

    pub fn count_propagators_in_integral(integral: AtomView) -> usize {
        let mut n_props: usize = 0;
        let mut prop_checker = integral.to_owned();
        while Pattern::parse("prop(args__)")
            .unwrap()
            .replace_iter(
                prop_checker.to_owned().as_view(),
                &S.one.into_pattern().into(),
                &Condition::default(),
                &MatchSettings::default(),
            )
            .next(&mut prop_checker)
            .is_some()
        {
            n_props += 1;
        }
        n_props
    }

    pub fn generate_topologies_with_contractions(
        canonical_expression: AtomView,
        short_expression: AtomView,
        contractions: TopologyContractions,
        matad_expression: Option<AtomView>,
        fmft_expression: Option<AtomView>,
    ) -> Result<Self, VakintError> {
        let n_tot_props = Topologies::count_propagators_in_integral(canonical_expression);
        if n_tot_props == 0 {
            return Err(VakintError::InvalidIntegralFormat(format!(
                "No propagators found in integral '{}'",
                canonical_expression
            )));
        }
        let mut topologies = vec![];

        let unique_contractions = match contractions {
            TopologyContractions::Custom(c) => c,
            TopologyContractions::Automatic => {
                let master_topology = Topology::generate_topology_with_contraction(
                    n_tot_props,
                    canonical_expression,
                    short_expression,
                    vec![],
                    matad_expression,
                    fmft_expression,
                )?;
                let cs = master_topology
                    .get_integral()
                    .graph
                    .find_unique_contractions();
                topologies.push(master_topology);
                cs
            }
        };

        for contracted_prop_indices in unique_contractions {
            topologies.push(Topology::generate_topology_with_contraction(
                n_tot_props,
                canonical_expression,
                short_expression,
                contracted_prop_indices,
                matad_expression,
                fmft_expression,
            )?);
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
    pub fn generate_topology_with_contraction(
        n_tot_props: usize,
        canonical_expression: AtomView,
        short_expression: AtomView,
        contraction: Vec<usize>,
        matad_expression: Option<AtomView>,
        fmft_expression: Option<AtomView>,
    ) -> Result<Topology, VakintError> {
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
        for prop_id in contraction.iter() {
            if let Some(m) = get_prop_with_id(contracted_canonical_expression.as_view(), *prop_id) {
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
        Ok(Integral::new(
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
        .into())
    }

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
                    // Make sure all propagators have IDs starting at 1 until n_props
                    let mut test_input = input.to_owned();
                    let mut prop_id = 1;
                    let mut prop_pattern =
                        Pattern::parse(format!("prop({},args__)", prop_id).as_str()).unwrap();
                    while prop_pattern
                        .pattern_match(
                            test_input.as_view(),
                            &Condition::default(),
                            &MatchSettings::default(),
                        )
                        .next()
                        .is_some()
                    {
                        test_input = prop_pattern.replace_all(
                            test_input.as_view(),
                            &Atom::new_num(1).into_pattern().into(),
                            None,
                            None,
                        );
                        prop_id += 1;
                        prop_pattern =
                            Pattern::parse(format!("prop({},args__)", prop_id).as_str()).unwrap();
                    }
                    test_input = Pattern::parse("UNKNOWN(args__)").unwrap().replace_all(
                        test_input.as_view(),
                        &Pattern::parse("args__").unwrap().into(),
                        None,
                        None,
                    );
                    if test_input != Atom::parse("topo(1)").unwrap() {
                        return Err(VakintError::InvalidIntegralFormat(format!(
                            "UNKNOWN integrals must have propagator ids ranging from 1 to their maximal number of propagators: {}",
                            test_input
                        )));
                    }

                    let mut loop_momenta_ids = HashSet::new();
                    for i_prop in 1..=prop_id - 1 {
                        for (mom_symbol, (_, mom_id)) in get_individual_momenta(
                            get_prop_with_id(input, i_prop)
                                .unwrap()
                                .get(&State::get_symbol("q_"))
                                .unwrap(),
                        )?
                        .iter()
                        {
                            if *mom_symbol != S.k {
                                return Err(VakintError::InvalidIntegralFormat(format!(
                                    "Unknown integrals must have only loop momenta with symbol {}: {}",
                                    S.k,mom_symbol
                                )));
                            }
                            if *mom_id < 1 {
                                return Err(VakintError::InvalidIntegralFormat(format!(
                                    "Unknown integrals must have loop momenta ids starting from 1, not {}",
                                    mom_id
                                )));
                            }
                            loop_momenta_ids.insert(*mom_id);
                        }
                    }
                    let mut sorted_loop_momenta_ids = loop_momenta_ids
                        .iter()
                        .map(|mom_id| *mom_id as usize)
                        .collect::<Vec<_>>();
                    sorted_loop_momenta_ids.sort();
                    if sorted_loop_momenta_ids
                        != (1..=sorted_loop_momenta_ids.len()).collect::<Vec<_>>()
                    {
                        return Err(VakintError::InvalidIntegralFormat(format!(
                            "Unknown integrals must have loop momenta ids ranging from 1 to their maximal number of loop momenta, not {:?}",
                            sorted_loop_momenta_ids
                        )));
                    }
                    let n_loops = *sorted_loop_momenta_ids.last().unwrap();
                    //println!("Found a match!");
                    let mut replacement_rules = ReplacementRules::default();
                    replacement_rules.canonical_expression_substitutions.insert(
                        Atom::parse("integral").unwrap(),
                        Pattern::parse("UNKNOWN(int_)").unwrap().replace_all(
                            input,
                            &Pattern::parse("int_").unwrap().into(),
                            None,
                            None,
                        ),
                    );
                    replacement_rules.canonical_expression_substitutions.insert(
                        Atom::parse("n_props").unwrap(),
                        Atom::new_num((prop_id - 1) as i64),
                    );
                    replacement_rules.canonical_expression_substitutions.insert(
                        Atom::parse("n_loops").unwrap(),
                        Atom::new_num(n_loops as i64),
                    );
                    Ok(Some(replacement_rules))
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
