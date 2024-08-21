use anyhow::Result;
use colored::Colorize;
#[allow(unused)]
use log::{debug, info, warn};

use regex::Regex;
use rug::float::Constant;
use std::{
    collections::{HashMap, HashSet},
    env,
    f64::consts::LOG2_10,
    fmt, fs,
    process::{Command, ExitStatus, Stdio},
};
use string_template_plus::{Render, RenderOptions, Template};
use symbolica::{
    atom::{
        representation::InlineNum, AsAtomView, Atom, AtomView, FunctionBuilder, SliceType, Symbol,
        Var,
    },
    domains::{
        float::{Complex, Float},
        integer::IntegerRing,
        rational::{Fraction, Rational},
    },
    fun,
    id::{Condition, Match, MatchSettings, Pattern, PatternRestriction, WildcardAndRestriction},
    printer::{AtomPrinter, PrintOptions},
    state::{FunctionAttribute, State},
};
use thiserror::Error;
use version_compare::{compare_to, Cmp};

use phf::phf_map;

static MINIMAL_FORM_VERSION: &str = "4.2.1";

#[allow(unused)]
static METRIC_SYMBOL: &str = "g";
static LOOP_MOMENTUM_SYMBOL: &str = "k";
static EXTERNAL_MOMENTUM_SYMBOL: &str = "p";

static FORM_SRC: phf::Map<&'static str, &'static str> = phf_map! {
    "integrateduv.frm" =>  include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/form_src/integrateduv.frm"
    )),
    "tensorreduce.frm" =>  include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/form_src/tensorreduce.frm"
    )),
    "pvtab10.h" =>  include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/form_src/pvtab10.h"
    )),
    "fmft.frm" =>  include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/form_src/fmft.frm"
    )),
};

static TEMPLATES: phf::Map<&'static str, &'static str> = phf_map! {
    "run_tensor_reduction.txt" =>  include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/templates/run_tensor_reduction.txt"
    )),
    "run_alphaloop_integral_evaluation.txt" => include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/templates/run_alphaloop_integral_evaluation.txt"
    ))
};

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
    #[error("invalid numerator expression: {0}")]
    InvalidNumerator(String),
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
    #[error("FORM run crashed with the following error:\nstderr: {0}\nYou can rerun the script using:\n{1}")]
    FormError(String, String, String),
    #[error("{0}")]
    FormVersion(String),
    #[error("FORM is not installed in your system and required for vakint to work.")]
    FormUnavailable,
    #[error("Could not find FORM output file 'out.txt':\nstderr: {0}\nYou can rerun the script using:\n{1}")]
    MissingFormOutput(String, String, String),
    #[error("Symbolica could not parse FORM output:\n{0}\nError:{1}")]
    FormOutputParsingError(String, String),
    #[error(transparent)]
    IoError(#[from] std::io::Error), // Add this variant to convert std::io::Error to VakintError
    #[error("{0}")]
    MalformedGraph(String),
    #[error("Symbolica error: {0}")]
    SymbolicaError(String),
    #[error("Numerical evaluation error: {0}")]
    EvaluationError(String),
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

fn even_condition() -> PatternRestriction {
    PatternRestriction::Filter(Box::new(move |m| {
        #[allow(warnings)]
        if let Match::Single(AtomView::Num(_)) = m {
            get_integer_from_match(&m).unwrap() % 2 == 0
        } else {
            false
        }
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
            canonical_topology: Topology::Unknown(Integral::default()),
            edge_ids_canonical_to_input_map: HashMap::new(),
            canonical_expression_substitutions: HashMap::new(),
            numerator_substitutions: HashMap::new(),
        }
    }
}

impl ReplacementRules {
    fn apply_replacement_rules(&mut self) {
        let integral = self.canonical_topology.get_integral_mut();
        for (source, target) in self.canonical_expression_substitutions.iter() {
            for expr in [
                integral.canonical_expression.as_mut(),
                integral.short_expression.as_mut(),
                integral.alphaloop_expression.as_mut(),
                integral.matad_expression.as_mut(),
                integral.fmft_expression.as_mut(),
            ]
            .into_iter()
            .flatten()
            {
                *expr = source.into_pattern().replace_all(
                    expr.as_view(),
                    &target.into_pattern(),
                    None,
                    None,
                );
            }
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
                            &Pattern::parse("1").unwrap(),
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
                            .replace_all(a.as_view(), &Pattern::parse("0").unwrap(), None, None);
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
                                &Pattern::parse(rhs.as_str()).unwrap(),
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

#[derive(Debug, Clone)]
struct Edge {
    id: usize,
    left_node_id: usize,
    right_node_id: usize,
    momentum: Atom,
}

impl fmt::Display for Edge {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "(#{}|{}->{}|{})",
            self.id, self.left_node_id, self.right_node_id, self.momentum
        )
    }
}

#[derive(Debug, Clone)]
struct Node {
    id: usize,
    edges: Vec<(usize, EdgeDirection)>,
}

impl fmt::Display for Node {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "(#{}|[{}])",
            self.id,
            self.edges
                .iter()
                .map(|(e_id, dir)| format!("{}@{}", dir, e_id))
                .collect::<Vec<_>>()
                .join(",")
        )
    }
}

#[derive(Debug, Clone)]
enum EdgeDirection {
    Incoming,
    Outgoing,
}

impl EdgeDirection {
    fn is_incoming(&self) -> bool {
        match self {
            EdgeDirection::Incoming => true,
            EdgeDirection::Outgoing => false,
        }
    }
    #[allow(unused)]
    fn is_outgoing(&self) -> bool {
        match self {
            EdgeDirection::Incoming => false,
            EdgeDirection::Outgoing => true,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct Graph {
    edges: HashMap<usize, Edge>,
    nodes: HashMap<usize, Node>,
}

impl fmt::Display for Graph {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut sorted_edges = self.edges.iter().collect::<Vec<_>>();
        sorted_edges.sort_by(|(e1_id, _), (e2_id, _)| e1_id.partial_cmp(e2_id).unwrap());
        let mut sorted_nodes = self.nodes.iter().collect::<Vec<_>>();
        sorted_nodes.sort_by(|(n1_id, _), (n2_id, _)| n1_id.partial_cmp(n2_id).unwrap());
        write!(
            f,
            "Edges: {}\nNodes: {}",
            sorted_edges
                .iter()
                .map(|(_e_id, e)| format!("{}", e))
                .collect::<Vec<_>>()
                .join(" "),
            sorted_nodes
                .iter()
                .map(|(_n_id, n)| format!("{}", n))
                .collect::<Vec<_>>()
                .join(" ")
        )
    }
}

impl fmt::Display for EdgeDirection {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            EdgeDirection::Incoming => write!(f, "IN"),
            EdgeDirection::Outgoing => write!(f, "OUT"),
        }
    }
}

impl Graph {
    pub fn to_graphviz(&self) -> String {
        format!(
            "digraph G {{\n{}\n}}",
            self.edges
                .values()
                .map(|e| format!(
                    "  {} -> {} [label=\"{}|{}\"]",
                    e.left_node_id, e.right_node_id, e.id, e.momentum
                ))
                .collect::<Vec<_>>()
                .join("\n")
        )
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
    fn get_integral(&self) -> &Integral {
        match self {
            Topology::OneLoop(i)
            | Topology::TwoLoop(i)
            | Topology::ThreeLoop(i)
            | Topology::FourLoop(i)
            | Topology::Unknown(i) => i,
        }
    }

    fn get_integral_mut(&mut self) -> &mut Integral {
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
pub struct Integral {
    n_loops: usize,
    n_props: usize,
    name: String,
    generic_pattern: FullPattern,
    canonical_expression: Option<Atom>,
    short_expression: Option<Atom>,
    short_expression_pattern: Option<Pattern>,
    alphaloop_expression: Option<Atom>,
    fmft_expression: Option<Atom>,
    matad_expression: Option<Atom>,
    #[allow(unused)]
    graph: Graph,
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
    let n = m.to_atom();
    match n.try_into() {
        Ok(res) => Some(res),
        Err(_) => None,
    }
    /*
    if let Match::Single(AtomView::Num(a)) = m {
        match a.get_coeff_view() {
            CoefficientView::Natural(n, 1) => Some(n),
            _ => None,
        }
    } else {
        None
    }
    */
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

impl Default for Integral {
    fn default() -> Self {
        Integral::new(0, None, None, None, None).unwrap()
    }
}

impl Integral {
    pub fn get_graph_from_expression(
        integral_expression: AtomView,
        tot_n_props: usize,
    ) -> Result<Graph, VakintError> {
        let mut graph = Graph::default();
        for i_prop in 1..=tot_n_props {
            if let Some(m) = get_prop_with_id(integral_expression, i_prop) {
                // Format check for the momenta
                let momentum = m.get(&State::get_symbol("q_")).unwrap().to_atom();

                let (left_node_id, right_node_id) = get_node_ids(&m)?;

                graph.edges.insert(
                    i_prop,
                    Edge {
                        id: i_prop,
                        momentum,
                        left_node_id,
                        right_node_id,
                    },
                );
            }
        }

        for (&e_id, edge) in graph.edges.iter() {
            if let Some(node) = graph.nodes.get_mut(&edge.left_node_id) {
                node.edges.push((e_id, EdgeDirection::Outgoing));
            } else {
                graph.nodes.insert(
                    edge.left_node_id,
                    Node {
                        id: edge.left_node_id,
                        edges: vec![(e_id, EdgeDirection::Outgoing)],
                    },
                );
            }
            if let Some(node) = graph.nodes.get_mut(&edge.right_node_id) {
                node.edges.push((e_id, EdgeDirection::Incoming));
            } else {
                graph.nodes.insert(
                    edge.right_node_id,
                    Node {
                        id: edge.right_node_id,
                        edges: vec![(e_id, EdgeDirection::Incoming)],
                    },
                );
            }
        }

        for (n_id, nodes) in graph.nodes.iter() {
            if nodes.edges.len() <= 1 {
                return Err(VakintError::MalformedGraph(format!("Node {} is connected to only {} edges, this cannot be for a vaccuum graph. Graph:\n{}",
                    n_id, nodes.edges.len(), integral_expression
                )));
            }
        }

        Ok(graph)
    }

    pub fn new(
        // This refers to the *total* number of propagators in the top-level topology,
        // i.e. the number of entries in the corresponding short_expression
        tot_n_props: usize,
        canonical_expression: Option<Atom>,
        short_expression: Option<Atom>,
        matad_expression: Option<Atom>,
        fmft_expression: Option<Atom>,
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
                canonical_expression: None,
                short_expression,
                short_expression_pattern,
                alphaloop_expression: None,
                fmft_expression: None,
                matad_expression: None,
                graph: Graph::default(),
            });
        }
        let e = canonical_expression.clone().unwrap();

        let graph = Integral::get_graph_from_expression(e.as_view(), tot_n_props)?;

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

        let mut alphaloop_expression = Atom::new_num(1);

        for (_n_id, node) in graph.nodes.iter() {
            let mut fb = FunctionBuilder::new(State::get_symbol("vxs"));
            for (e_id, dir) in &node.edges {
                fb = fb.add_arg(
                    &(graph.edges.get(e_id).unwrap().momentum.clone()
                        * if dir.is_incoming() { 1 } else { -1 }),
                );
            }
            alphaloop_expression = alphaloop_expression * fb.finish();
        }
        for (&e_id, edge) in graph.edges.iter() {
            alphaloop_expression = alphaloop_expression
                * fun!(
                    State::get_symbol("uvprop"),
                    &edge.momentum,
                    fun!(State::get_symbol("pow"), Atom::new_num(e_id as i64))
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
            alphaloop_expression: Some(alphaloop_expression),
            fmft_expression,
            matad_expression,
            graph,
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
                        // We do not want to match propagators with zero powers in the short form,
                        // as these should be matched to the pinched version with a hardcoded zero power
                        if a.is_zero() {
                            return Ok(None);
                        }
                        replacement_rules.canonical_expression_substitutions.insert(
                            Atom::parse(format!("pow({})", i_prop).as_str()).unwrap(),
                            a.to_owned(),
                        );
                    }
                    if let Some(Match::Single(a)) = m1
                        .match_stack
                        .get(State::get_symbol(format!("msq{}_", i_prop)))
                    {
                        replacement_rules.canonical_expression_substitutions.insert(
                            Atom::parse(format!("msq({})", i_prop).as_str()).unwrap(),
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

pub enum EvaluationApproach {
    AlphaLoop,
    MATAD,
    FMFT,
    PySecDec,
}

impl fmt::Display for EvaluationApproach {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            EvaluationApproach::AlphaLoop => write!(f, "AlphaLoop"),
            EvaluationApproach::MATAD => write!(f, "MATTAD"),
            EvaluationApproach::FMFT => write!(f, "FMFT"),
            EvaluationApproach::PySecDec => write!(f, "PySecDec"),
        }
    }
}

impl EvaluationApproach {
    pub fn supports(&self, vakint: &Vakint, topology: &Topology) -> bool {
        match self {
            EvaluationApproach::AlphaLoop => {
                topology.get_integral().alphaloop_expression.is_some()
                    && vakint.settings.number_of_terms_in_epsilon_expansion <= 3
                    && topology.get_integral().n_loops <= 3
            }
            EvaluationApproach::MATAD => {
                topology.get_integral().matad_expression.is_some()
                    && topology.get_integral().n_loops <= 3
            }
            EvaluationApproach::FMFT => {
                topology.get_integral().fmft_expression.is_some()
                    && topology.get_integral().n_loops == 4
            }
            EvaluationApproach::PySecDec => true,
        }
    }

    pub fn evaluate_integral(
        &self,
        vakint: &Vakint,
        numerator: AtomView,
        integral_specs: &ReplacementRules,
    ) -> Result<Atom, VakintError> {
        match self {
            EvaluationApproach::AlphaLoop => {
                vakint.alpha_loop_evaluate(vakint, numerator, integral_specs)
            }
            EvaluationApproach::MATAD => unimplemented!("MATAD evaluations not yet implemented"),
            EvaluationApproach::FMFT => unimplemented!("FMFT evaluations not yet implemented"),
            EvaluationApproach::PySecDec => {
                unimplemented!("PySecDec evaluations not yet implemented")
            }
        }
    }
}

pub struct VakintSettings {
    #[allow(unused)]
    pub use_pysecdec: bool,
    pub epsilon_symbol: String,
    pub mu_r_sq_symbol: String,
    pub form_exe_path: String,
    pub verify_numerator_identification: bool,
    pub integral_normalization: String,
    // TODO optimize for 16 digits -> f64
    pub n_digits_at_evaluation_time: u32,
    pub allow_unknown_integrals: bool,
    pub evaluation_order: Vec<EvaluationApproach>,
    // This quantity is typically set equal to *one plus the maximum loop count* of the UV regularisation problem considered.
    // For example when considering a 2-loop problem, then:
    //   a) for the nested one-loop integrals appearing, the single pole, finite term *and* order-epsilon term will need to be considered.
    //   b) for the two-loop integrals, the double pole, single pole and finite terms will be needed, so again three terms
    pub number_of_terms_in_epsilon_expansion: usize,
    pub use_dot_product_notation: bool,
}

#[allow(clippy::derivable_impls)]
impl Default for VakintSettings {
    fn default() -> Self {
        VakintSettings {
            use_pysecdec: false,
            epsilon_symbol: "ε".into(),
            mu_r_sq_symbol: "mursq".into(),
            form_exe_path: env::var("FORM_PATH").unwrap_or("form".into()),
            verify_numerator_identification: true,
            n_digits_at_evaluation_time: 16,
            integral_normalization: "1".into(),
            allow_unknown_integrals: true,
            evaluation_order: vec![
                EvaluationApproach::AlphaLoop,
                EvaluationApproach::MATAD,
                EvaluationApproach::FMFT,
                EvaluationApproach::PySecDec,
            ],
            // Default to a three-loop UV subtraction problem, for which alphaLoop implementation can be used.
            number_of_terms_in_epsilon_expansion: 3,
            use_dot_product_notation: false,
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
    pub vectors: Vec<(String, i64)>,
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
    pub fn evaluate_integral(&mut self, vakint: &Vakint) -> Result<(), VakintError> {
        let mut integral_specs = if let Some(replacement_rules) = vakint
            .topologies
            .match_topologies_to_user_input(self.integral.as_view())?
        {
            replacement_rules
        } else {
            return Err(VakintError::UnreckognizedIntegral(
                self.integral.to_string(),
            ));
        };

        integral_specs.apply_replacement_rules();

        'eval: for evaluation_approach in vakint.settings.evaluation_order.iter() {
            if evaluation_approach.supports(vakint, &integral_specs.canonical_topology) {
                self.numerator = evaluation_approach.evaluate_integral(
                    vakint,
                    self.numerator.as_atom_view(),
                    &integral_specs,
                )?;
                self.integral = Atom::new_num(1);
                break 'eval;
            }
        }
        Ok(())
    }

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
        test = Pattern::parse(format!("{}(momID_,idx_)", EXTERNAL_MOMENTUM_SYMBOL).as_str())
            .unwrap()
            .replace_all(
                test.as_view(),
                &Atom::parse("1").unwrap().into_pattern(),
                None,
                None,
            );

        // Substitute metric in as well
        test = Pattern::parse("g(idx1_,idx2_)").unwrap().replace_all(
            test.as_view(),
            &Atom::parse("1").unwrap().into_pattern(),
            None,
            None,
        );

        if vakint.settings.verify_numerator_identification && !matches!(test, Atom::Num(_)) {
            return Err(VakintError::NumeratorNotReplaced(
                EXTERNAL_MOMENTUM_SYMBOL.into(),
                test.to_string(),
            ));
        }
        self.numerator = new_numerator;

        self.vectors = VakintTerm::identify_vectors_in_numerator(self.numerator.as_view())?;

        Ok(())
    }

    pub fn identify_vectors_in_numerator(
        numerator: AtomView,
    ) -> Result<Vec<(String, i64)>, VakintError> {
        let mut vectors = HashSet::new();

        let vector_matcher_pattern = Pattern::parse("vec_(id_,idx_)").unwrap();
        let vector_conditions = Condition::from((State::get_symbol("vec_"), symbol_condition()))
            & Condition::from((State::get_symbol("id_"), number_condition()))
            & Condition::from((State::get_symbol("idx_"), number_condition()));
        let vector_match_settings = MatchSettings::default();
        let mut vector_matcher = vector_matcher_pattern.pattern_match(
            numerator,
            &vector_conditions,
            &vector_match_settings,
        );

        while let Some(m) = vector_matcher.next() {
            if let Match::FunctionName(vec_symbol) =
                m.match_stack.get(State::get_symbol("vec_")).unwrap()
            {
                if *vec_symbol == State::get_symbol(LOOP_MOMENTUM_SYMBOL)
                    || *vec_symbol == State::get_symbol(EXTERNAL_MOMENTUM_SYMBOL)
                {
                    vectors.insert((
                        vec_symbol.to_string(),
                        get_integer_from_match(
                            m.match_stack.get(State::get_symbol("id_")).unwrap(),
                        )
                        .unwrap(),
                    ));
                }
            } else {
                unreachable!("Vector name should be a symbol.")
            }
        }
        Ok(vectors.iter().cloned().collect::<Vec<_>>())
    }

    pub fn tensor_reduce(&mut self, vakint: &Vakint) -> Result<(), VakintError> {
        let mut form_numerator = self.numerator.clone();

        // Make sure to undo the dot product notation.
        // If it was not used, the command below will do nothing.
        form_numerator = Vakint::convert_from_dot_notation(form_numerator.as_view());
        for (vec, id) in self.vectors.iter() {
            form_numerator = Pattern::parse(format!("{}({},idx_)", vec, id).as_str())
                .unwrap()
                .replace_all(
                    form_numerator.as_view(),
                    &Pattern::parse(
                        format!(
                            "vec{}({}{},idx_)",
                            if *vec == EXTERNAL_MOMENTUM_SYMBOL {
                                "1"
                            } else {
                                ""
                            },
                            vec,
                            id
                        )
                        .as_str(),
                    )
                    .unwrap(),
                    Some(&Condition::from((
                        State::get_symbol("idx_"),
                        number_condition(),
                    ))),
                    None,
                )
        }

        let template =
            Template::parse_template(TEMPLATES.get("run_tensor_reduction.txt").unwrap()).unwrap();
        let mut vars: HashMap<String, String> = HashMap::new();
        vars.insert(
            "numerator".into(),
            vakint.prepare_expression_for_form(form_numerator)?,
        );
        let rendered = template
            .render(&RenderOptions {
                variables: vars,
                ..Default::default()
            })
            .unwrap();
        let form_result = vakint.run_form(
            &["tensorreduce.frm".into(), "pvtab10.h".into()],
            ("run_tensor_reduction.frm".into(), rendered),
            vec![],
        )?;
        let mut reduced_numerator = vakint.process_form_output(form_result)?;
        for (vec, id) in self.vectors.iter() {
            reduced_numerator = Pattern::parse(format!("{}{}", vec, id).as_str())
                .unwrap()
                .replace_all(
                    reduced_numerator.as_view(),
                    &Pattern::parse(format!("{}({})", vec, id).as_str()).unwrap(),
                    None,
                    None,
                );
        }

        if !vakint.settings.use_dot_product_notation {
            reduced_numerator = Vakint::convert_from_dot_notation(reduced_numerator.as_view());
        }

        self.numerator = reduced_numerator;
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
                vectors: VakintTerm::identify_vectors_in_numerator(numerator.as_view())?,
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

    pub fn tensor_reduce(&mut self, vakint: &Vakint) -> Result<(), VakintError> {
        for term in self.0.iter_mut() {
            term.tensor_reduce(vakint)?;
        }
        Ok(())
    }

    pub fn evaluate_integral(&mut self, vakint: &Vakint) -> Result<(), VakintError> {
        for term in self.0.iter_mut() {
            term.evaluate_integral(vakint)?;
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
                    vectors: term.vectors.clone(),
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
                    vectors: term.vectors.clone(),
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

pub struct NumericalEvaluationResult(Vec<(i64, Complex<Float>)>);

impl fmt::Display for NumericalEvaluationResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            self.0
                .iter()
                .map(|(power, float_eval)| format!(
                    "{} : {}",
                    format!(
                        "ε^{}{}",
                        match *power {
                            power if power > 0 => "+",
                            0 => " ",
                            power if power < 0 => "",
                            _ => unreachable!(),
                        },
                        power
                    )
                    .green(),
                    format!("{}", float_eval).blue()
                ))
                .collect::<Vec<_>>()
                .join("\n")
        )
    }
}

impl NumericalEvaluationResult {
    pub fn get_epsilon_coefficients(&self) -> Vec<(i64, Complex<Float>)> {
        self.0.clone()
    }
}

impl Vakint {
    fn initialize_symbolica_symbols() {
        // unoriented edge
        _ = State::get_symbol_with_attributes("uedge", &[FunctionAttribute::Symmetric]);
        // dot product
        _ = State::get_symbol_with_attributes("dot", &[FunctionAttribute::Symmetric]);
        // metric tensor
        _ = State::get_symbol_with_attributes("g", &[FunctionAttribute::Symmetric]);
    }

    pub fn get_coeff_map(prec: u32) -> impl Fn(&Fraction<IntegerRing>) -> Complex<Float> {
        move |x| Complex::new(x.to_multi_prec_float(prec), Float::with_val(prec, 0.))
    }

    pub fn get_constants_map(
        &self,
        params: &HashMap<String, f64, ahash::RandomState>,
    ) -> Result<HashMap<Atom, Complex<Float>, ahash::random_state::RandomState>, VakintError> {
        let mut const_map: HashMap<Atom, Complex<Float>, ahash::random_state::RandomState> =
            HashMap::default();

        let binary_prec: u32 =
            ((self.settings.n_digits_at_evaluation_time as f64) * LOG2_10).floor() as u32;
        const_map.insert(
            Atom::from(Var::new(State::get_symbol("𝜋"))),
            Complex::new(
                Float::with_val(binary_prec, Constant::Pi),
                Float::with_val(binary_prec, 0),
            ),
        );

        const_map.insert(
            Atom::from(Var::new(State::get_symbol("𝑖"))),
            Complex::new(
                Float::with_val(binary_prec, 0),
                Float::with_val(binary_prec, 1),
            ),
        );

        const_map.insert(
            fun!(State::LOG, Atom::new_num(2)),
            Complex::new(
                Float::with_val(binary_prec, Constant::Log2),
                Float::new(binary_prec),
            ),
        );

        for (symb, value) in params.iter() {
            const_map.insert(
                Atom::parse(symb.as_str()).unwrap(),
                Complex::new(
                    Float::with_val(binary_prec, value),
                    Float::with_val(binary_prec, 0),
                ),
            );
        }

        Ok(const_map)
    }

    pub fn partial_numerical_evaluation(
        &self,
        integral: AtomView,
        params: &HashMap<String, f64, ahash::RandomState>,
    ) -> Atom {
        let const_map = self.get_constants_map(params).unwrap();

        let mut res = integral.to_owned();
        for (src, trgt) in const_map.iter() {
            res = src.into_pattern().replace_all(
                res.as_view(),
                &((Atom::new_num(trgt.re.clone())
                    + Atom::new_var(State::get_symbol("𝑖")) * Atom::new_num(trgt.im.clone()))
                .into_pattern()),
                None,
                None,
            );
        }

        res
    }

    pub fn full_numerical_evaluation(
        &self,
        integral: AtomView,
        params: &HashMap<String, f64, ahash::RandomState>,
    ) -> Result<NumericalEvaluationResult, VakintError> {
        let epsilon_coeffs = integral.coefficient_list(State::get_symbol("ε"));

        let mut epsilon_coeffs_vec = epsilon_coeffs
            .0
            .iter()
            .map(|(eps_atom, coeff)| {
                if let Some(m) = Pattern::parse("ε^n_")
                    .unwrap()
                    .pattern_match(
                        eps_atom.as_view(),
                        &Condition::from((State::get_symbol("n_"), number_condition())),
                        &MatchSettings::default(),
                    )
                    .next()
                {
                    (
                        get_integer_from_match(m.match_stack.get(State::get_symbol("n_")).unwrap())
                            .unwrap(),
                        coeff,
                    )
                } else if *eps_atom == Atom::parse("ε").unwrap() {
                    (1, coeff)
                } else {
                    panic!("Epsilon atom should be of the form ε^n_")
                }
            })
            .collect::<Vec<_>>();

        epsilon_coeffs_vec.push((0, &epsilon_coeffs.1));

        let map = self.get_constants_map(params).unwrap();
        let map_view = map.iter().map(|(k, v)| (k.as_view(), v.clone())).collect();
        let binary_prec: u32 =
            ((self.settings.n_digits_at_evaluation_time as f64) * LOG2_10).floor() as u32;
        let mut epsilon_coeffs_vec_floats = vec![];
        for (i64, coeff) in epsilon_coeffs_vec.iter() {
            epsilon_coeffs_vec_floats.push((
                *i64,
                match coeff.evaluate(
                    &Vakint::get_coeff_map(binary_prec),
                    &map_view,
                    &HashMap::default(),
                    &mut HashMap::default(),
                ) {
                    Ok(x) => x,
                    Err(e) => {
                        return Err(VakintError::EvaluationError(format!(
                            "Is some tensor structure left?: {}",
                            e
                        )));
                    }
                },
            ));
        }

        epsilon_coeffs_vec_floats.sort_by(|(i1, _), (i2, _)| i1.cmp(i2));
        Ok(NumericalEvaluationResult(epsilon_coeffs_vec_floats))
    }

    fn alpha_loop_evaluate(
        &self,
        vakint: &Vakint,
        numerator: AtomView,
        integral_specs: &ReplacementRules,
    ) -> Result<Atom, VakintError> {
        let integral = integral_specs.canonical_topology.get_integral();
        let alphaloop_expression = integral.alphaloop_expression.as_ref().unwrap().as_view();

        /*
        println!("Numerator : {}", numerator);
        println!("Evaluating AlphaLoop : {}", alphaloop_expression);
        println!("Graph:\n{}", integral.graph.to_graphviz());
        */

        // Make sure to undo the dot product notation.
        // If it was not used, the command below will do nothing.
        let mut form_expression = numerator.to_owned() * alphaloop_expression.to_owned();
        form_expression = Vakint::convert_to_dot_notation(form_expression.as_view());
        while let Some(m) = Pattern::parse("k(id_)")
            .unwrap()
            .pattern_match(
                form_expression.as_view(),
                &Condition::from((State::get_symbol("id_"), number_condition())),
                &MatchSettings::default(),
            )
            .next()
        {
            let k_id = get_integer_from_match(m.match_stack.get(State::get_symbol("id_")).unwrap())
                .unwrap();
            form_expression = Pattern::parse(format!("k({})", k_id).as_str())
                .unwrap()
                .replace_all(
                    form_expression.as_view(),
                    &Pattern::parse(format!("k{}", k_id).as_str()).unwrap(),
                    None,
                    None,
                );
        }

        let template = Template::parse_template(
            TEMPLATES
                .get("run_alphaloop_integral_evaluation.txt")
                .unwrap(),
        )
        .unwrap();
        let mut vars: HashMap<String, String> = HashMap::new();
        vars.insert(
            "numerator".into(),
            vakint.prepare_expression_for_form(form_expression)?,
        );
        let rendered = template
            .render(&RenderOptions {
                variables: vars,
                ..Default::default()
            })
            .unwrap();
        let form_result = vakint.run_form(
            &["integrateduv.frm".into()],
            ("run_alphaloop_integral_evaluation.frm".into(), rendered),
            vec![
                "-D".into(),
                format!("MAXPOLE={}", integral.n_loops),
                "-D".into(),
                format!(
                    "SELECTEDEPSILONORDER={}",
                    vakint.settings.number_of_terms_in_epsilon_expansion - integral.n_loops
                ),
            ],
        )?;

        let mut evaluated_integral = vakint.process_form_output(form_result)?;

        if vakint.settings.use_dot_product_notation {
            evaluated_integral = Vakint::convert_from_dot_notation(evaluated_integral.as_view());
        }

        let muv_sq_symbol = if let Some(m) = Pattern::parse("prop(args__,m_,pow_)")
            .unwrap()
            .pattern_match(
                integral.canonical_expression.as_ref().unwrap().as_view(),
                &Condition::default(),
                &MatchSettings::default(),
            )
            .next()
        {
            m.match_stack
                .get(State::get_symbol("m_"))
                .unwrap()
                .to_atom()
        } else {
            return Err(VakintError::MalformedGraph(format!(
                "Could not find muV in graph:\n{}",
                integral.canonical_expression.as_ref().unwrap()
            )));
        };

        evaluated_integral = Pattern::parse("mUV").unwrap().replace_all(
            evaluated_integral.as_view(),
            &muv_sq_symbol
                .pow((Atom::new_num(1) / Atom::new_num(2)).as_atom_view())
                .into_pattern(),
            None,
            None,
        );

        evaluated_integral = Pattern::parse("logmUVmu").unwrap().replace_all(
            evaluated_integral.as_view(),
            &((Atom::new_num(1) / Atom::new_num(2))
                * fun!(
                    State::LOG,
                    muv_sq_symbol
                        / Atom::new_var(State::get_symbol(vakint.settings.mu_r_sq_symbol.as_str()))
                ))
            .into_pattern(),
            None,
            None,
        );

        evaluated_integral = evaluated_integral
            * Atom::parse(vakint.settings.integral_normalization.as_str()).unwrap();

        let expanded_evaluation = match evaluated_integral.series(
            State::get_symbol(vakint.settings.epsilon_symbol.as_str()),
            Atom::Zero.as_atom_view(),
            Rational::from(vakint.settings.number_of_terms_in_epsilon_expansion - integral.n_loops),
            true,
        ) {
            Ok(a) => a,
            Err(e) => return Err(VakintError::SymbolicaError(e.to_string())),
        };

        /*
        println!("exp res: {}", expanded_evaluation);
        println!("exp res a: {}", expanded_evaluation.to_atom());
        let tt = expanded_evaluation
            .to_atom()
            .coefficient_list(State::get_symbol(vakint.settings.epsilon_symbol.as_str()));
        println!("miaou {}", tt.1);
        for (els, l) in tt.0 {
            println!("exp res coefs a: {} {}", els, l);
        }
        */

        Ok(expanded_evaluation.to_atom())
    }

    fn get_form_version(&self) -> Result<String, VakintError> {
        let mut cmd = Command::new(self.settings.form_exe_path.as_str());
        cmd.arg("-version");
        let output = if let Ok(o) = cmd.output() {
            o
        } else {
            return Err(VakintError::FormUnavailable);
        };

        if !ExitStatus::success(&output.status) {
            return Err(VakintError::FormUnavailable);
        }
        let output_str = String::from_utf8_lossy(&output.stdout).into_owned();
        let re = Regex::new(r"FORM\s([\.|\d]+)").unwrap();
        let mut versions = vec![];
        for (_, [version]) in re.captures_iter(output_str.as_str()).map(|c| c.extract()) {
            versions.push(version);
        }
        if versions.is_empty() {
            return Err(VakintError::FormVersion(format!(
                "Could not obtain form version from command:\n{:?}\nOutput was:\n{}",
                cmd, output_str
            )));
        }
        Ok(versions[0].into())
    }

    pub fn convert_from_dot_notation(atom: AtomView) -> Atom {
        let mut expr = atom.to_owned();
        let mut running_dummy_index = 1;

        while let Some(m) = Pattern::parse("dot(v1_(id1_),v2_(id2_))")
            .unwrap()
            .pattern_match(
                expr.as_view(),
                &(Condition::from((State::get_symbol("v1_"), symbol_condition()))
                    & Condition::from((State::get_symbol("v2_"), symbol_condition()))
                    & Condition::from((State::get_symbol("id1_"), number_condition()))
                    & Condition::from((State::get_symbol("id2_"), number_condition()))),
                &MatchSettings::default(),
            )
            .next()
        {
            let (id1, id2) = (
                get_integer_from_match(m.match_stack.get(State::get_symbol("id1_")).unwrap())
                    .unwrap(),
                get_integer_from_match(m.match_stack.get(State::get_symbol("id2_")).unwrap())
                    .unwrap(),
            );
            let (v1, v2) = match (
                m.match_stack.get(State::get_symbol("v1_")),
                m.match_stack.get(State::get_symbol("v2_")),
            ) {
                (Some(Match::FunctionName(s1)), Some(Match::FunctionName(s2))) => (s1, s2),
                _ => unreachable!("Vectors have to be symbols"),
            };

            expr = m.target.into_pattern().replace_all(
                expr.as_view(),
                &Pattern::parse(
                    format!(
                        "{}({},{})*{}({},{})",
                        v1, id1, running_dummy_index, v2, id2, running_dummy_index
                    )
                    .as_str(),
                )
                .unwrap(),
                None,
                None,
            );
            running_dummy_index += 1;
        }
        expr
    }

    pub fn convert_to_dot_notation(atom: AtomView) -> Atom {
        let mut expr = atom.to_owned();
        expr = Pattern::parse("v_(id_,idx_)^n_").unwrap().replace_all(
            expr.as_view(),
            &Pattern::parse("dot(id_,id_)^(n_/2)").unwrap(),
            Some(
                &(Condition::from((State::get_symbol("v_"), symbol_condition()))
                    & Condition::from((State::get_symbol("id_"), number_condition()))
                    & Condition::from((State::get_symbol("idx_"), number_condition()))
                    & Condition::from((State::get_symbol("n_"), even_condition()))),
            ),
            None,
        );

        // dot products
        expr = Pattern::parse("v1_(id1_,idx_)*v2_(id2_,idx_)")
            .unwrap()
            .replace_all(
                expr.as_view(),
                &Pattern::parse("dot(v1_(id1_),v2_(id2_))").unwrap(),
                Some(
                    &(Condition::from((State::get_symbol("v1_"), symbol_condition()))
                        & Condition::from((State::get_symbol("v2_"), symbol_condition()))
                        & Condition::from((State::get_symbol("id1_"), number_condition()))
                        & Condition::from((State::get_symbol("id2_"), number_condition()))
                        & Condition::from((State::get_symbol("idx_"), number_condition()))),
                ),
                None,
            );

        // metric contraction
        expr = Pattern::parse("g(idx1_,idx2_)*v1_(id1_,idx1_)*v2_(id2_,idx2_)")
            .unwrap()
            .replace_all(
                expr.as_view(),
                &Pattern::parse("dot(v1_(id1_),v2_(id2_))").unwrap(),
                Some(
                    &(Condition::from((State::get_symbol("v1_"), symbol_condition()))
                        & Condition::from((State::get_symbol("v2_"), symbol_condition()))
                        & Condition::from((State::get_symbol("id1_"), number_condition()))
                        & Condition::from((State::get_symbol("id2_"), number_condition()))
                        & Condition::from((State::get_symbol("idx1_"), number_condition()))
                        & Condition::from((State::get_symbol("idx2_"), number_condition()))),
                ),
                None,
            );

        expr
    }

    pub fn new(settings: Option<VakintSettings>) -> Result<Self, VakintError> {
        Vakint::initialize_symbolica_symbols();
        let vakint_settings = settings.unwrap_or_default();
        let topologies = Topologies::generate_topologies(&vakint_settings)?;
        let vakint = Vakint {
            settings: vakint_settings,
            topologies,
        };
        let form_version = vakint.get_form_version()?;
        match compare_to(form_version.clone(), MINIMAL_FORM_VERSION, Cmp::Ge) {
            Ok(valid) => {
                if valid {
                    debug!(
                        "FORM successfully detected with version '{}'.",
                        form_version
                    );
                } else {
                    return Err(VakintError::FormVersion(format!(
                        "FORM version installed on your system does not meet minimal requirements: {}<{}.",
                        form_version, MINIMAL_FORM_VERSION
                    )));
                }
            }
            Err(_) => {
                return Err(VakintError::FormVersion(format!(
                    "Could not parse FORM version '{}'.",
                    form_version
                )))
            }
        };
        Ok(vakint)
    }

    pub fn prepare_expression_for_form(&self, expression: Atom) -> Result<String, VakintError> {
        let processed = expression.clone();
        let processed = Pattern::parse(&self.settings.epsilon_symbol)
            .unwrap()
            .replace_all(
                processed.as_view(),
                &Pattern::parse("ep").unwrap(),
                None,
                None,
            );
        Ok(AtomPrinter::new_with_options(processed.as_view(), PrintOptions::file()).to_string())
    }

    pub fn process_form_output(&self, form_output: String) -> Result<Atom, VakintError> {
        match Atom::parse(form_output.replace("i_", "𝑖").as_str()) {
            Ok(mut processed) => {
                processed = Pattern::parse("rat(x_,y_)").unwrap().replace_all(
                    processed.as_view(),
                    &Pattern::parse("(x_/y_)").unwrap(),
                    None,
                    None,
                );
                processed = Pattern::parse("rat(x_)").unwrap().replace_all(
                    processed.as_view(),
                    &Pattern::parse("x_").unwrap(),
                    None,
                    None,
                );
                processed = Pattern::parse("ep").unwrap().replace_all(
                    processed.as_view(),
                    &Pattern::parse(&self.settings.epsilon_symbol).unwrap(),
                    None,
                    None,
                );
                processed = Pattern::parse("pi").unwrap().replace_all(
                    processed.as_view(),
                    &Pattern::parse("𝜋").unwrap(),
                    None,
                    None,
                );
                processed = Pattern::parse("g(idx1_,idx2_)").unwrap().replace_all(
                    processed.as_view(),
                    &Pattern::parse(format!("{}(idx1_,idx2_)", METRIC_SYMBOL).as_str()).unwrap(),
                    Some(
                        &(Condition::from((State::get_symbol("idx1_"), number_condition()))
                            & Condition::from((State::get_symbol("idx2_"), number_condition()))),
                    ),
                    None,
                );
                processed = Pattern::parse("g(v1_,v2_)").unwrap().replace_all(
                    processed.as_view(),
                    &Pattern::parse("dot(v1_,v2_)").unwrap(),
                    Some(
                        &(Condition::from((State::get_symbol("v1_"), symbol_condition()))
                            & Condition::from((State::get_symbol("v2_"), symbol_condition()))),
                    ),
                    None,
                );
                processed = Pattern::parse("g(v1_(args1_),v2_(args2_))")
                    .unwrap()
                    .replace_all(
                        processed.as_view(),
                        &Pattern::parse("dot(v1_(args1_),v2_(args2_))").unwrap(),
                        Some(
                            &(Condition::from((State::get_symbol("v1_"), symbol_condition()))
                                & Condition::from((State::get_symbol("v2_"), symbol_condition()))),
                        ),
                        None,
                    );
                Ok(processed)
            }
            Err(err) => Err(VakintError::FormOutputParsingError(form_output, err)),
        }
    }

    pub fn run_form(
        &self,
        resources: &[String],
        input: (String, String),
        options: Vec<String>,
    ) -> Result<String, VakintError> {
        let tmp_dir = &env::temp_dir().join("vakint_temp");
        if tmp_dir.exists() {
            fs::remove_dir_all(tmp_dir)?;
        }
        fs::create_dir(tmp_dir)?;

        for resource in resources.iter() {
            fs::write(tmp_dir.join(resource), FORM_SRC.get(resource).unwrap())?;
        }
        fs::write(tmp_dir.join(&input.0).to_str().unwrap(), &input.1)?;
        let mut cmd = Command::new(self.settings.form_exe_path.as_str());
        for opt in options {
            cmd.arg(opt);
        }
        cmd.arg(input.0);
        let output = cmd
            .current_dir(tmp_dir)
            .stderr(Stdio::inherit())
            .stdout(Stdio::null())
            .output()?;
        if !ExitStatus::success(&output.status) {
            return Err(VakintError::FormError(
                String::from_utf8_lossy(&output.stderr).into(),
                format!("{:?}", cmd),
                tmp_dir.to_str().unwrap().into(),
            ));
        }
        if !tmp_dir.join("out.txt").exists() {
            return Err(VakintError::MissingFormOutput(
                String::from_utf8_lossy(&output.stderr).into(),
                format!("{:?}", cmd),
                tmp_dir.to_str().unwrap().into(),
            ));
        }
        let result = fs::read_to_string(tmp_dir.join("out.txt"))?;
        fs::remove_dir_all(tmp_dir)?;
        Ok(result)
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

    pub fn tensor_reduce(&self, input: AtomView) -> Result<Atom, VakintError> {
        let mut vakint_expr = VakintExpression::try_from(input)?;
        vakint_expr.tensor_reduce(self)?;
        Ok(vakint_expr.into())
    }

    pub fn evaluate_integral(&self, input: AtomView) -> Result<Atom, VakintError> {
        let mut vakint_expr = VakintExpression::try_from(input)?;
        vakint_expr.evaluate_integral(self)?;
        Ok(vakint_expr.into())
    }
}
