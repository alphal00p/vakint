pub mod fmft;
pub mod fmft_numerics;
pub mod graph;
pub mod matad;
pub mod matad_numerics;
pub mod symbols;
pub mod topologies;
pub mod utils;

use ahash::RandomState;
use anyhow::Result;
use colored::Colorize;
use graph::Graph;
#[allow(unused)]
use log::{debug, info, warn};

use regex::Regex;
use rug::float::Constant;
use std::{
    collections::{hash_map::Entry, HashMap, HashSet},
    env,
    f64::consts::LOG2_10,
    fmt,
    fs::{self, File},
    io::{BufRead, BufReader, Write},
    ops::Div,
    path::PathBuf,
    process::{Command, ExitStatus, Stdio},
    sync::{Arc, LazyLock, Mutex},
};
use string_template_plus::{Render, RenderOptions, Template};
#[allow(unused)]
use symbolica::id::PatternOrMap;
use symbolica::{
    atom::{
        representation::InlineNum, AsAtomView, Atom, AtomView, FunctionBuilder, SliceType, Symbol,
        Var,
    },
    domains::{
        atom::AtomField,
        float::{Complex, Float, NumericalFloatLike, Real, RealNumberLike, SingleFloat},
        integer::IntegerRing,
        rational::{Fraction, Rational},
    },
    fun,
    id::{
        Condition, Match, MatchSettings, Pattern, PatternRestriction, Replacement,
        WildcardRestriction,
    },
    poly::series::Series,
    printer::{AtomPrinter, PrintOptions},
    state::State,
    transformer::Transformer,
};
use utils::simplify_real;

use thiserror::Error;
use topologies::{Topologies, Topology};
use version_compare::{compare_to, Cmp};

use symbols::{EXTERNAL_MOMENTUM_SYMBOL, LOOP_MOMENTUM_SYMBOL, METRIC_SYMBOL, S};

use phf::phf_map;

pub type Momentum = (
    Complex<Float>,
    Complex<Float>,
    Complex<Float>,
    Complex<Float>,
);

static MINIMAL_FORM_VERSION: &str = "4.2.1";
static MINIMAL_PYSECDEC_VERSION: &str = "1.6.4";

static FORM_SRC: phf::Map<&'static str, &'static str> = phf_map! {
    "integrateduv.frm" =>  include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/form_src/alphaloop/integrateduv.frm"
    )),
    "tensorreduce.frm" =>  include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/form_src/alphaloop/tensorreduce.frm"
    )),
    "pvtab10.h" =>  include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/form_src/alphaloop/pvtab10.h"
    )),
    "fmft.frm" =>  include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/form_src/fmft/fmft.frm"
    )),
    "matad-ng.hh" =>  include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/form_src/matad/matad-ng.hh"
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
    )),
    "run_pySecDec_template.txt" => include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/templates/run_pySecDec_template.txt"
    )),
    "run_matad.txt" => include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/templates/run_matad.txt"
    )),
    "run_fmft.txt" => include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/templates/run_fmft.txt"
    )),
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
    #[error("Could not find a method suitable for evaluating this integral up to 𝒪(ε^{1}): {0}")]
    NoEvaluationMethodFound(String, i64),
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
    #[error("FORM run crashed with the following error:\nstderr: {0}\nstdout: {1}\nYou can rerun the script using:\n{2}")]
    FormError(String, String, String, String),
    #[error("{0}")]
    FormVersion(String),
    #[error("FORM is not installed in your system and required for vakint to work.")]
    FormUnavailable,
    #[error("PySecDec error: {0}")]
    PySecDecError(String),
    #[error("{0}")]
    PySecDecVersion(String),
    #[error("PySecDec is not installed in your system and required for vakint to evaluate with PySecDec. Install it with 'pip install pysecdec'")]
    PySecDecUnavailable,
    #[error("Could not find FORM output file 'out.txt':\nstderr: {0}\nYou can rerun the script using:\n{1}")]
    MissingFormOutput(String, String, String),
    #[error("Symbolica could not parse PySecDec output:\n{0}\nError:{1}")]
    PySecDecOutputParsingError(String, String),
    #[error("Symbolica could not parse FORM output:\n{0}\nError:{1}")]
    FormOutputParsingError(String, String),
    #[error(transparent)]
    IoError(#[from] std::io::Error), // Add this variant to convert std::io::Error to VakintError
    #[error("{0}")]
    MalformedGraph(String),
    #[error("Invalid loop normalization factor specified: {0}.\nError: {1}\nNote that the following symbols are allowed in the expression: {2}")]
    InvalidLoopNormalization(String, String, String),
    #[error("Symbolica error: {0}")]
    SymbolicaError(String),
    #[error("MATAD error: {0}")]
    MATADError(String),
    #[error("FMFT error: {0}")]
    FMFTError(String),
    #[error("Numerical evaluation error: {0}")]
    EvaluationError(String),
    #[error("unknown vakint error")]
    Unknown,
}

pub fn params_from_f64(
    params: &HashMap<String, f64, ahash::RandomState>,
    decimal_prec: u32,
) -> HashMap<String, Complex<Float>, ahash::RandomState> {
    let binary_prec: u32 = ((decimal_prec as f64) * LOG2_10).floor() as u32;
    params
        .iter()
        .map(|(k, v)| {
            (
                k.clone(),
                Complex::new(
                    Float::with_val(binary_prec, v),
                    Float::with_val(binary_prec, 0),
                ),
            )
        })
        .collect()
}

pub fn externals_from_f64(
    externals: &HashMap<usize, (f64, f64, f64, f64), ahash::RandomState>,
    decimal_prec: u32,
) -> HashMap<usize, Momentum, ahash::RandomState> {
    let binary_prec: u32 = ((decimal_prec as f64) * LOG2_10).floor() as u32;
    externals
        .iter()
        .map(|(&k, v)| {
            (
                k,
                (
                    Complex::new(
                        Float::with_val(binary_prec, v.0),
                        Float::with_val(binary_prec, 0.0),
                    ),
                    Complex::new(
                        Float::with_val(binary_prec, v.1),
                        Float::with_val(binary_prec, 0.0),
                    ),
                    Complex::new(
                        Float::with_val(binary_prec, v.2),
                        Float::with_val(binary_prec, 0.0),
                    ),
                    Complex::new(
                        Float::with_val(binary_prec, v.3),
                        Float::with_val(binary_prec, 0.0),
                    ),
                ),
            )
        })
        .collect()
}

fn propagators_condition() -> WildcardRestriction {
    WildcardRestriction::Filter(Box::new(move |m| {
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

fn gt_condition(value: i64) -> WildcardRestriction {
    WildcardRestriction::Filter(Box::new(move |m| {
        if let Match::Single(AtomView::Num(a)) = m {
            a.get_coeff_view() > InlineNum::new(value, 1).as_num_view().get_coeff_view()
        } else {
            false
        }
    }))
}

fn eq_condition(value: i64) -> WildcardRestriction {
    WildcardRestriction::Filter(Box::new(move |m| {
        if let Match::Single(AtomView::Num(a)) = m {
            a.get_coeff_view() == InlineNum::new(value, 1).as_num_view().get_coeff_view()
        } else {
            false
        }
    }))
}

#[allow(unused)]
fn lt_condition(value: i64) -> WildcardRestriction {
    WildcardRestriction::Filter(Box::new(move |m| {
        if let Match::Single(AtomView::Num(a)) = m {
            a.get_coeff_view() < InlineNum::new(value, 1).as_num_view().get_coeff_view()
        } else {
            false
        }
    }))
}

fn range_condition(min: i64, max: i64) -> WildcardRestriction {
    WildcardRestriction::Filter(Box::new(move |m| {
        if let Match::Single(AtomView::Num(a)) = m {
            a.get_coeff_view() <= InlineNum::new(max, 1).as_num_view().get_coeff_view()
                && a.get_coeff_view() >= InlineNum::new(min, 1).as_num_view().get_coeff_view()
        } else {
            false
        }
    }))
}

fn number_condition() -> WildcardRestriction {
    WildcardRestriction::Filter(Box::new(move |m| {
        matches!(m, Match::Single(AtomView::Num(_)))
    }))
}

fn even_condition() -> WildcardRestriction {
    WildcardRestriction::Filter(Box::new(move |m| {
        #[allow(warnings)]
        if let Match::Single(AtomView::Num(_)) = m {
            get_integer_from_match(&m).unwrap() % 2 == 0
        } else {
            false
        }
    }))
}

fn symbol_or_number() -> WildcardRestriction {
    WildcardRestriction::Filter(Box::new(move |m| {
        matches!(
            m,
            Match::Single(AtomView::Num(_))
                | Match::Single(AtomView::Var(_))
                | Match::FunctionName(_)
        )
    }))
}

fn symbol_condition() -> WildcardRestriction {
    WildcardRestriction::Filter(Box::new(move |m| {
        matches!(m, Match::Single(AtomView::Var(_)) | Match::FunctionName(_))
    }))
}

fn apply_restriction_to_symbols(
    symbols: Vec<Symbol>,
    restriction: &WildcardRestriction,
) -> Condition<PatternRestriction> {
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
    fn apply_replacement_rules(&mut self) -> Result<(), VakintError> {
        if matches!(self.canonical_topology, Topology::Unknown(_)) {
            let integral = self.canonical_topology.get_integral_mut();
            integral.canonical_expression = Some(
                self.canonical_expression_substitutions
                    .get(&Atom::parse("integral").unwrap())
                    .unwrap()
                    .to_owned(),
            );
            let n_props: i64 = self
                .canonical_expression_substitutions
                .get(&Atom::parse("n_props").unwrap())
                .unwrap()
                .try_into()
                .unwrap();
            integral.n_props = n_props as usize;
            let n_loops: i64 = self
                .canonical_expression_substitutions
                .get(&Atom::parse("n_loops").unwrap())
                .unwrap()
                .try_into()
                .unwrap();
            integral.n_loops = n_loops as usize;
            integral.graph = Graph::new_from_atom(
                integral.canonical_expression.as_ref().unwrap().as_view(),
                integral.n_props,
            )?;
            return Ok(());
        }

        let integral = self.canonical_topology.get_integral_mut();

        for (source, target) in self.canonical_expression_substitutions.iter() {
            for expr in [
                integral.canonical_expression.as_mut(),
                integral.short_expression.as_mut(),
                integral.alphaloop_expression.as_mut(),
            ]
            .into_iter()
            .flatten()
            {
                *expr = source.into_pattern().replace_all(
                    expr.as_view(),
                    &target.into_pattern().into(),
                    None,
                    None,
                );
            }
        }
        Ok(())
    }

    fn get_propagator_property_list(&self, property: &str) -> HashMap<usize, Atom> {
        let mut property_list = HashMap::new();
        let integral = self.canonical_topology.get_integral();
        let canonical_expression_view = integral.canonical_expression.as_ref().unwrap().as_view();
        for prop_id in 1..=integral.n_props {
            if let Some(m) = get_prop_with_id(canonical_expression_view, prop_id) {
                // println!("property {}", property);
                // println!(
                //     "key {}",
                //     &m.get(&State::get_symbol(property)).unwrap().to_atom()
                // );
                // println!("self.canonical_expression_substitution {}", self);
                let v_a = m.get(&State::get_symbol(property)).unwrap().to_atom();
                property_list.insert(prop_id, v_a);
            }
        }
        property_list
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
    applicable_evaluation_methods: EvaluationOrder,
    graph: Graph,
    node_pairs: HashMap<Symbol, HashSet<Symbol>>,
    unoriented_generic_pattern: FullPattern,
}

impl fmt::Display for Integral {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "name='{}', n_loops={}, n_props_top_topo={}\n   | canonical_expression='{}',\n   | canonical_pattern='{}',\n   | short_expression='{}'",
            if self.name=="UNKNOWN" {self.name.red()} else {self.name.green()},
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
fn get_individual_momenta_from_atom(
    momentum: AtomView,
) -> Result<Vec<(Symbol, (Atom, i64))>, VakintError> {
    let wrapped_momentum = fun!(State::get_symbol("mom"), momentum.to_owned());
    if let Some(m) = Pattern::parse("mom(q_)")
        .unwrap()
        .pattern_match(
            wrapped_momentum.as_view(),
            &Condition::default(),
            &MatchSettings::default(),
        )
        .next()
    {
        get_individual_momenta(m.match_stack.get(State::get_symbol("q_")).unwrap())
    } else {
        Err(VakintError::InvalidMomentumExpression(format!(
            "Edge momentum {} does not contain only 'q_'.",
            momentum
        )))
    }
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
            &Atom::Zero.into_pattern().into(),
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
        Integral::new(0, None, None, EvaluationOrder::empty()).unwrap()
    }
}

impl Integral {
    pub fn new(
        // This refers to the *total* number of propagators in the top-level topology,
        // i.e. the number of entries in the corresponding short_expression
        tot_n_props: usize,
        canonical_expression: Option<Atom>,
        short_expression: Option<Atom>,
        applicable_evaluation_methods: EvaluationOrder,
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
                generic_pattern: all_accepting_pattern.clone(),
                canonical_expression: None,
                short_expression,
                short_expression_pattern,
                alphaloop_expression: None,
                applicable_evaluation_methods,
                graph: Graph::default(),
                unoriented_generic_pattern: all_accepting_pattern,
                node_pairs: HashMap::new(),
            });
        }
        let e = canonical_expression.clone().unwrap();

        let graph = Graph::new_from_atom(e.as_view(), tot_n_props)?;

        let mut loop_mom_indices = HashSet::<i64>::default();
        let mut next_atom = e.clone();
        let mut old_atom = e.clone();

        let mut generic_expression = Atom::parse("1").unwrap();
        let mut generic_condition: Condition<PatternRestriction> = Condition::default();
        let mut unoriented_generic_expression = Atom::parse("1").unwrap();
        let mut unoriented_generic_condition: Condition<PatternRestriction> = Condition::default();
        let mut node_pairs = HashMap::<Symbol, HashSet<Symbol>>::default();

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
                unoriented_generic_expression = unoriented_generic_expression * Atom::parse(
                        format!("prop(id{id}_,uedge(n{id_l_node}_,n{id_r_node}_),q{id}__,{mass},{power})",
                        id=i_prop,id_l_node=id_node_left, id_r_node=id_node_right,
                        mass=mass_symbol_string, power=pow_symbol_string).as_str()).unwrap();

                let new_conditions =
                    Condition::from((State::get_symbol(mass_symbol_string), symbol_or_number()))
                        & Condition::from((
                            State::get_symbol(format!("pow{}_", i_prop).as_str()),
                            number_condition(),
                        ))
                        & Condition::from((
                            State::get_symbol(format!("id{}_", i_prop).as_str()),
                            number_condition(),
                        ));
                generic_condition = generic_condition
                    & new_conditions.clone()
                    & apply_restriction_to_symbols(
                        vec![
                            State::get_symbol(format!("n{}l_", id_node_left).as_str()),
                            State::get_symbol(format!("n{}r_", id_node_right).as_str()),
                        ],
                        &symbol_or_number(),
                    );
                unoriented_generic_condition = unoriented_generic_condition
                    & new_conditions
                    & apply_restriction_to_symbols(
                        vec![
                            State::get_symbol(format!("n{}_", id_node_left).as_str()),
                            State::get_symbol(format!("n{}_", id_node_right).as_str()),
                        ],
                        &symbol_or_number(),
                    );
                for (id, side) in [(id_node_left, "l"), (id_node_right, "r")] {
                    let new_entry = State::get_symbol(format!("n{}{}_", id, side));
                    node_pairs
                        .entry(State::get_symbol(format!("n{}_", id)))
                        .and_modify(|v| {
                            v.insert(new_entry);
                        })
                        .or_insert(HashSet::from([new_entry]));
                }
                next_atom = Pattern::parse(format!("prop({},args__)", i_prop).as_str())
                    .unwrap()
                    .replace_all(
                        old_atom.as_view(),
                        &Pattern::parse("1").unwrap().into(),
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
        let unoriented_generic_pattern = FullPattern {
            pattern: unoriented_generic_expression.into_pattern(),
            conditions: unoriented_generic_condition,
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
                    &Pattern::parse(format!("pow{}_", i_prop).as_str())
                        .unwrap()
                        .into(),
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
                    &Pattern::parse(format!("msq{}_", i_prop).as_str())
                        .unwrap()
                        .into(),
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
            applicable_evaluation_methods,
            graph,
            unoriented_generic_pattern,
            node_pairs,
        })
    }

    fn match_integral_to_short_user_input(
        &self,
        input: AtomView,
    ) -> Result<Option<ReplacementRules>, VakintError> {
        let unwrapped_input = Pattern::parse("topo(integral_)").unwrap().replace_all(
            input,
            &Pattern::parse("integral_").unwrap().into(),
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
                        // NO: We do not want to match propagators with zero powers in the short form,
                        // as these should be matched to the pinched version with a hardcoded zero power
                        // ZERO_POWERS: Better to keep them to allow use to match to higher-level inputs
                        // if a.is_zero() {
                        //     return Ok(None);
                        // }
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
            &Pattern::parse("uedge(x_,y_)").unwrap().into(),
            None,
            None,
        );

        // Make sure that the unoriented topology could at least match, otherwise abort immediately
        if self
            .unoriented_generic_pattern
            .pattern
            .pattern_match(
                undirected_input.as_view(),
                &self.unoriented_generic_pattern.conditions,
                &MatchSettings::default(),
            )
            .next()
            .is_none()
        {
            return Ok(None);
        }

        let mut topology_matcher = self.generic_pattern.pattern.pattern_match(
            undirected_input.as_view(),
            &self.generic_pattern.conditions,
            &self.generic_pattern.match_settings,
        );
        'outer: while let Some(m1) = topology_matcher.next() {
            // println!("VHDEBUG match result:");
            // for (k, v) in m1.match_stack.get_matches() {
            //     println!("{} -> {}", k, v.to_atom());
            // }
            // Make sure that all matched nodes are distinct
            let mut node_matches: HashMap<Symbol, HashSet<Atom>> = HashMap::default();
            for (unoriented_node, oriented_nodes) in self.node_pairs.iter() {
                node_matches
                    .entry(*unoriented_node)
                    .or_insert(HashSet::from_iter(
                        oriented_nodes
                            .iter()
                            .map(|&on| m1.match_stack.get(on).unwrap().to_atom()),
                    ));
            }

            // This check that all oriented node point to the same value
            if node_matches.values().any(|nodes| nodes.len() != 1) {
                continue 'outer;
            }
            // This check that all unoriented nodes point to different values
            if HashSet::<Atom>::from_iter(
                node_matches
                    .values()
                    .map(|nodes| nodes.iter().next().unwrap())
                    .cloned(),
            )
            .len()
                != node_matches.len()
            {
                // We do not need a continue here because there is not enough distinct nodes w.r.t canonical expression I believe
                break 'outer;
            }

            let mut replacement_rules = ReplacementRules::default();

            for prop_id in 1..=self.n_props {
                if let Some(canonical_prop_match) = get_prop_with_id(
                    self.canonical_expression.as_ref().unwrap().as_view(),
                    prop_id,
                ) {
                    for var_prop_id in 1..=self.n_props {
                        if let Some(mtmp) = m1
                            .match_stack
                            .get(State::get_symbol(format!("pow{}_", var_prop_id).as_str()))
                        {
                            if let Entry::Vacant(e) =
                                replacement_rules.canonical_expression_substitutions.entry(
                                    Atom::parse(format!("pow({})", var_prop_id).as_str()).unwrap(),
                                )
                            {
                                e.insert(mtmp.to_atom());
                            }
                        }
                        if let Some(mtmp) = m1
                            .match_stack
                            .get(State::get_symbol(format!("msq{}_", var_prop_id)))
                        {
                            if let Entry::Vacant(e) =
                                replacement_rules.canonical_expression_substitutions.entry(
                                    Atom::parse(format!("msq({})", var_prop_id).as_str()).unwrap(),
                                )
                            {
                                e.insert(mtmp.to_atom());
                            }
                        }
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
                            Pattern::parse("k(ilmb_)").unwrap().replace_all(
                                canonical_prop_match
                                    .get(&State::get_symbol("q_"))
                                    .unwrap()
                                    .to_atom()
                                    .as_view(),
                                &Pattern::parse("k(ilmb_,idx_)").unwrap().into(),
                                Some(&Condition::from((
                                    State::get_symbol("ilmb_"),
                                    number_condition(),
                                ))),
                                Some(&MatchSettings {
                                    allow_new_wildcards_on_rhs: true,
                                    ..MatchSettings::default()
                                }),
                            );
                        if is_edge_flipped {
                            canonical_momenta_atom_for_pattern =
                                (canonical_momenta_atom_for_pattern * -1).expand()
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
            // println!("VHDEBUG FROM integral: {}", self);
            // println!("VHDEBUG Found match: {}", replacement_rules);
            // panic!("VHDEBUG STOP");
            return Ok(Some(replacement_rules));
        }

        // let msg = [
        //     format!("User topology: {}", undirected_input),
        //     format!(
        //         "Unoriented pattern: {}",
        //         self.unoriented_generic_pattern.pattern.to_atom().unwrap()
        //     ),
        //     format!(
        //         "Oriented pattern: {}",
        //         self.generic_pattern.pattern.to_atom().unwrap()
        //     ),
        // ];
        // unreachable!("There should have been a match in Vakint at this stage. This is a logic error in vakint.\n{}",msg.join("\n"));

        Ok(None)
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
                &target.into_pattern().into(),
                None,
                None,
            );
        }

        // NO: Remove propagators with zero powers
        // ZERO_POWERS: Better to keep them to allow use to match to higher-level inputs
        // new_expression = Pattern::parse("prop(propID_,edge(nl_,nr_),q_,mUVsq_,0)")
        //     .unwrap()
        //     .replace_all(
        //         new_expression.as_view(),
        //         &Pattern::parse("1").unwrap().into(),
        //         None,
        //         None,
        //     );

        new_expression
    }
}

#[derive(Debug, Clone)]
pub struct PySecDecOptions {
    pub quiet: bool,
    pub relative_precision: f64,
    pub numerical_masses: HashMap<String, f64>,
    pub numerical_external_momenta: HashMap<String, (f64, f64, f64, f64)>,
    pub min_n_evals: u64,
    pub max_n_evals: u64,
    pub reuse_existing_output: Option<String>,
}

impl Default for PySecDecOptions {
    fn default() -> Self {
        // Give some random values to the external and masses. The user is expected to change these.
        let mut numerical_masses = HashMap::default();
        numerical_masses.insert("muvsq".into(), 1.0);
        let mut numerical_external_momenta = HashMap::default();
        for i in 1..=10 {
            numerical_external_momenta.insert(format!("p{}", i), (13.0, 4.0, 3.0, 12.0));
        }
        PySecDecOptions {
            quiet: true,
            relative_precision: 1.0e-7,
            numerical_masses,
            numerical_external_momenta,
            min_n_evals: 10_000,
            max_n_evals: 1_000_000_000_000,
            reuse_existing_output: None,
        }
    }
}

impl fmt::Display for PySecDecOptions {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "quiet={}, relative_precision={:.e}, min_n_evals={}, max_n_evals={}",
            self.quiet, self.relative_precision, self.min_n_evals, self.max_n_evals
        )
    }
}

#[derive(Debug, Clone)]
pub struct FMFTOptions {
    pub expand_masters: bool,
    pub susbstitute_masters: bool,
}

impl Default for FMFTOptions {
    fn default() -> Self {
        FMFTOptions {
            expand_masters: true,
            susbstitute_masters: true,
        }
    }
}

impl fmt::Display for FMFTOptions {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "expand_masters={}, susbstitute_masters={}",
            self.expand_masters, self.susbstitute_masters
        )
    }
}

#[derive(Debug, Clone)]
pub struct MATADOptions {
    pub expand_masters: bool,
    pub susbstitute_masters: bool,
    pub substitute_hpls: bool,
    pub direct_numerical_substition: bool,
}

impl Default for MATADOptions {
    fn default() -> Self {
        MATADOptions {
            expand_masters: true,
            susbstitute_masters: true,
            substitute_hpls: true,
            direct_numerical_substition: true,
        }
    }
}

impl fmt::Display for MATADOptions {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "expand_masters={}, susbstitute_masters={}, substitute_hpls={}, direct_numerical_substition={}",
            self.expand_masters, self.susbstitute_masters, self.substitute_hpls, self.direct_numerical_substition
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VakintDependency {
    FORM,
    PySecDec,
}

#[derive(Debug, Clone)]
pub enum EvaluationMethod {
    AlphaLoop,
    MATAD(MATADOptions),
    FMFT(FMFTOptions),
    PySecDec(PySecDecOptions),
}

impl fmt::Display for EvaluationMethod {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            EvaluationMethod::AlphaLoop => write!(f, "AlphaLoop"),
            EvaluationMethod::MATAD(opts) => write!(f, "MATAD ({})", opts),
            EvaluationMethod::FMFT(opts) => write!(f, "FMFT ({})", opts),
            EvaluationMethod::PySecDec(opts) => {
                write!(f, "PySecDec ({})", opts)
                //write!(f, "PySecDec")
            }
        }
    }
}

impl EvaluationMethod {
    pub fn supports(&self, vakint: &Vakint, topology: &Topology) -> bool {
        match self {
            EvaluationMethod::AlphaLoop => {
                topology
                    .get_integral()
                    .applicable_evaluation_methods
                    .0
                    .iter()
                    .any(|m| matches!(m, EvaluationMethod::AlphaLoop))
                    && topology.get_integral().alphaloop_expression.is_some()
                    && vakint.settings.number_of_terms_in_epsilon_expansion <= 4
                    && topology.get_integral().n_loops <= 3
            }
            EvaluationMethod::MATAD(_) => {
                topology
                    .get_integral()
                    .applicable_evaluation_methods
                    .0
                    .iter()
                    .any(|m| matches!(m, EvaluationMethod::MATAD(_)))
                    && topology.get_integral().n_loops <= 3
            }
            EvaluationMethod::FMFT(_) => {
                topology
                    .get_integral()
                    .applicable_evaluation_methods
                    .0
                    .iter()
                    .any(|m| matches!(m, EvaluationMethod::FMFT(_)))
                    && topology.get_integral().n_loops == 4
            }
            EvaluationMethod::PySecDec(_) => topology
                .get_integral()
                .applicable_evaluation_methods
                .0
                .iter()
                .any(|m| matches!(m, EvaluationMethod::PySecDec(_))),
        }
    }

    pub fn dependencies(&self) -> Vec<VakintDependency> {
        match self {
            EvaluationMethod::AlphaLoop => {
                vec![VakintDependency::FORM]
            }
            EvaluationMethod::MATAD(_) => {
                vec![VakintDependency::FORM]
            }
            EvaluationMethod::FMFT(_) => {
                vec![VakintDependency::FORM]
            }
            EvaluationMethod::PySecDec(_) => {
                vec![VakintDependency::FORM, VakintDependency::PySecDec]
            }
        }
    }

    pub fn adjust(
        &mut self,
        quiet: Option<bool>,
        relative_precision: f64,
        numerical_masses: &HashMap<String, Complex<Float>, RandomState>,
        numerical_external_momenta: &HashMap<usize, Momentum, RandomState>,
    ) {
        match self {
            EvaluationMethod::AlphaLoop => {}
            EvaluationMethod::MATAD(_) => {}
            EvaluationMethod::FMFT(_) => {}
            EvaluationMethod::PySecDec(opts) => {
                let f64_numerical_masses = numerical_masses
                    .iter()
                    .map(|(k, v)| (k.clone(), v.norm().re.to_f64()))
                    .collect();
                let f64_numerical_external_momenta = numerical_external_momenta
                    .iter()
                    .map(|(k, v)| {
                        (
                            format!("p{}", k),
                            (
                                v.0.norm().re.to_f64(),
                                v.1.norm().re.to_f64(),
                                v.2.norm().re.to_f64(),
                                v.3.norm().re.to_f64(),
                            ),
                        )
                    })
                    .collect();
                if let Some(is_quiet) = quiet {
                    opts.quiet = is_quiet;
                }
                opts.relative_precision = relative_precision;
                opts.numerical_masses = f64_numerical_masses;
                opts.numerical_external_momenta = f64_numerical_external_momenta;
            }
        }
    }

    pub fn evaluate_integral(
        &self,
        vakint: &Vakint,
        numerator: AtomView,
        integral_specs: &ReplacementRules,
    ) -> Result<Atom, VakintError> {
        let result = match self {
            EvaluationMethod::AlphaLoop => {
                vakint.alphaloop_evaluate(vakint, numerator, integral_specs)
            }
            EvaluationMethod::MATAD(opts) => {
                vakint.matad_evaluate(vakint, numerator, integral_specs, opts)
            }
            EvaluationMethod::FMFT(opts) => {
                vakint.fmft_evaluate(vakint, numerator, integral_specs, opts)
            }
            EvaluationMethod::PySecDec(opts) => {
                vakint.pysecdec_evaluate(vakint, numerator, integral_specs, opts)
            }
        }?;
        // Simplify logarithms and zero powers knowing that all arguments are real
        Ok(simplify_real(result.as_view()))
    }
}

#[derive(Debug, Clone)]
pub struct EvaluationOrder(pub Vec<EvaluationMethod>);

impl Default for EvaluationOrder {
    fn default() -> Self {
        EvaluationOrder(vec![
            EvaluationMethod::AlphaLoop,
            EvaluationMethod::MATAD(MATADOptions::default()),
            EvaluationMethod::FMFT(FMFTOptions::default()),
            EvaluationMethod::PySecDec(PySecDecOptions::default()),
        ])
    }
}

impl fmt::Display for EvaluationOrder {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "[{}]",
            self.0
                .iter()
                .map(|m| format!("{}", m))
                .collect::<Vec<_>>()
                .join("|")
        )
    }
}

impl EvaluationOrder {
    pub fn numerical_only() -> Self {
        EvaluationOrder(vec![EvaluationMethod::PySecDec(PySecDecOptions::default())])
    }
    pub fn analytic_only() -> Self {
        EvaluationOrder(vec![
            EvaluationMethod::AlphaLoop,
            EvaluationMethod::MATAD(MATADOptions::default()),
            EvaluationMethod::FMFT(FMFTOptions::default()),
        ])
    }
    pub fn alphaloop_only() -> Self {
        EvaluationOrder(vec![EvaluationMethod::AlphaLoop])
    }
    pub fn matad_only(matad_options: Option<MATADOptions>) -> Self {
        EvaluationOrder(vec![EvaluationMethod::MATAD(
            matad_options.unwrap_or_default(),
        )])
    }
    pub fn fmft_only(fmft_options: Option<FMFTOptions>) -> Self {
        EvaluationOrder(vec![EvaluationMethod::FMFT(
            fmft_options.unwrap_or_default(),
        )])
    }
    pub fn pysecdec_only(pysecdec_options: Option<PySecDecOptions>) -> Self {
        EvaluationOrder(vec![EvaluationMethod::PySecDec(
            pysecdec_options.unwrap_or_default(),
        )])
    }
    pub fn empty() -> Self {
        EvaluationOrder(vec![])
    }
    pub fn all() -> Self {
        EvaluationOrder(vec![
            EvaluationMethod::AlphaLoop,
            EvaluationMethod::MATAD(MATADOptions::default()),
            EvaluationMethod::FMFT(FMFTOptions::default()),
            EvaluationMethod::PySecDec(PySecDecOptions::default()),
        ])
    }
    pub fn all_but_fmft() -> Self {
        EvaluationOrder(vec![
            EvaluationMethod::AlphaLoop,
            EvaluationMethod::MATAD(MATADOptions::default()),
            EvaluationMethod::PySecDec(PySecDecOptions::default()),
        ])
    }

    pub fn adjust(
        &mut self,
        quiet: Option<bool>,
        relative_precision: f64,
        numerical_masses: &HashMap<String, Complex<Float>, RandomState>,
        numerical_external_momenta: &HashMap<usize, Momentum, RandomState>,
    ) {
        for method in self.0.iter_mut() {
            method.adjust(
                quiet,
                relative_precision,
                numerical_masses,
                numerical_external_momenta,
            );
        }
    }
}

#[derive(Debug, Clone)]
pub struct VakintSettings {
    #[allow(unused)]
    pub epsilon_symbol: String,
    pub mu_r_sq_symbol: String,
    pub form_exe_path: String,
    pub python_exe_path: String,
    pub verify_numerator_identification: bool,
    pub integral_normalization_factor: LoopNormalizationFactor,
    pub run_time_decimal_precision: u32,
    pub allow_unknown_integrals: bool,
    pub clean_tmp_dir: bool,
    pub evaluation_order: EvaluationOrder,
    // This quantity is typically set equal to *one plus the maximum loop count* of the UV regularisation problem considered.
    // For example when considering a 2-loop problem, then:
    //   a) for the nested one-loop integrals appearing, the single pole, finite term *and* order-epsilon term will need to be considered.
    //   b) for the two-loop integrals, the double pole, single pole and finite terms will be needed, so again three terms
    pub number_of_terms_in_epsilon_expansion: i64,
    pub use_dot_product_notation: bool,
    pub temporary_directory: Option<String>,
}

impl VakintSettings {
    pub fn get_integral_normalization_factor_atom(&self) -> Result<Atom, VakintError> {
        self.integral_normalization_factor.to_atom(self)
    }

    pub fn get_binary_precision(&self) -> u32 {
        ((self.run_time_decimal_precision.max(17) as f64) * LOG2_10).floor() as u32
    }

    pub fn real_to_prec(&self, re: &str) -> Complex<Float> {
        let prec = self.get_binary_precision();
        Complex::new(
            Float::parse(re, Some(prec)).unwrap(),
            Float::parse("0", Some(prec)).unwrap(),
        )
    }

    pub fn complex_to_prec(&self, re: &str, im: &str) -> Complex<Float> {
        let prec = self.get_binary_precision();
        Complex::new(
            Float::parse(re, Some(prec)).unwrap(),
            Float::parse(im, Some(prec)).unwrap(),
        )
    }
}

#[derive(Debug, Clone)]
pub enum LoopNormalizationFactor {
    #[allow(non_camel_case_types)]
    pySecDec,
    MSbar,
    FMFTandMATAD,
    Custom(String),
}

impl LoopNormalizationFactor {
    pub fn to_expression(&self) -> String {
        match self {
            LoopNormalizationFactor::pySecDec => "(𝑖*(𝜋^((4-2*eps)/2)))^(-n_loops)".into(),
            LoopNormalizationFactor::FMFTandMATAD => {
                "( 𝑖*(𝜋^((4-2*eps)/2)) * (exp(-EulerGamma))^(eps) )^(-n_loops)".into()
            }
            LoopNormalizationFactor::MSbar => {
                "(exp(log_mu_sq)/(4*𝜋*exp(-EulerGamma)))^(eps*n_loops)".into()
            }
            LoopNormalizationFactor::Custom(s) => s.clone(),
        }
    }

    pub fn static_allowed_symbols() -> Vec<String> {
        vec![
            "eps".into(),
            "log_mu_sq".into(),
            "EulerGamma".into(),
            "𝑖".into(),
            "I".into(),
            "𝜋".into(),
            "pi".into(),
        ]
    }

    pub fn allowed_symbols(settings: &VakintSettings) -> Vec<String> {
        let mut allowed_symbols = LoopNormalizationFactor::static_allowed_symbols();
        allowed_symbols.push(settings.epsilon_symbol.clone());
        allowed_symbols
    }

    pub fn to_atom(&self, settings: &VakintSettings) -> Result<Atom, VakintError> {
        let mut a = Atom::try_from(self)?;
        a = Pattern::parse("eps").unwrap().replace_all(
            a.as_view(),
            &Pattern::parse(&settings.epsilon_symbol).unwrap().into(),
            None,
            None,
        );
        Ok(a)
    }

    pub fn validate(
        &self,
        settings: &VakintSettings,
    ) -> Result<(Atom, Series<AtomField>, NumericalEvaluationResult), VakintError> {
        let expr: Atom = self.to_atom(settings)?;
        let expanded_expr = Pattern::parse("n_loops").unwrap().replace_all(
            expr.as_view(),
            &Atom::new_num(1).into_pattern().into(),
            None,
            None,
        );

        let expanded_expr = match expanded_expr.series(
            State::get_symbol(settings.epsilon_symbol.as_str()),
            Atom::Zero.as_atom_view(),
            Rational::from(settings.number_of_terms_in_epsilon_expansion - 1),
            true,
        ) {
            Ok(a) => a,
            Err(e) => return Err(VakintError::SymbolicaError(e.to_string())),
        };

        let mut expanded_expr_atom = expanded_expr.to_atom();
        let log_mu_sq = fun!(
            State::LOG,
            Atom::new_var(State::get_symbol(settings.mu_r_sq_symbol.as_str()))
        );
        expanded_expr_atom = Pattern::parse("log_mu_sq").unwrap().replace_all(
            expanded_expr_atom.as_view(),
            &(log_mu_sq).into_pattern().into(),
            None,
            None,
        );

        let mut params = HashMap::default();
        params.insert("mursq".into(), settings.real_to_prec("1"));
        let num_res = match Vakint::full_numerical_evaluation_without_error(
            settings,
            expanded_expr_atom.as_view(),
            &params,
            None,
        ) {
            Ok(r) => r,
            Err(e) => {
                return Err(VakintError::InvalidLoopNormalization(
                    format!("{}", expr),
                    e.to_string(),
                    LoopNormalizationFactor::allowed_symbols(settings).join(","),
                ))
            }
        };
        Ok((expr, expanded_expr, num_res))
    }
}

impl fmt::Display for LoopNormalizationFactor {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            LoopNormalizationFactor::pySecDec => {
                write!(
                    f,
                    "{} convention [ {} ]",
                    String::from("PySecDec").green(),
                    self.to_expression()
                )
            }
            LoopNormalizationFactor::MSbar => {
                write!(
                    f,
                    "{} convention [ {} ]",
                    String::from("MSbar").green(),
                    self.to_expression()
                )
            }
            LoopNormalizationFactor::FMFTandMATAD => {
                write!(
                    f,
                    "{} convention [ {} ]",
                    String::from("FMFTandMATAD").green(),
                    self.to_expression()
                )
            }
            LoopNormalizationFactor::Custom(_) => {
                write!(
                    f,
                    "{} convention [ {} ]",
                    String::from("Custom").green(),
                    self.to_expression()
                )
            }
        }
    }
}

impl TryFrom<&LoopNormalizationFactor> for Atom {
    type Error = VakintError;
    fn try_from(expr: &LoopNormalizationFactor) -> Result<Self, VakintError> {
        match Atom::parse(&expr.to_expression()) {
            Ok(a) => {
                let mut processed_a = Pattern::parse("I").unwrap().replace_all(
                    a.as_view(),
                    &Pattern::parse("𝑖").unwrap().into(),
                    None,
                    None,
                );
                processed_a = Pattern::parse("pi").unwrap().replace_all(
                    processed_a.as_view(),
                    &Pattern::parse("𝜋").unwrap().into(),
                    None,
                    None,
                );
                Ok(processed_a)
            }
            Err(e) => Err(VakintError::InvalidLoopNormalization(
                format!("{}", expr),
                e.to_string(),
                LoopNormalizationFactor::static_allowed_symbols().join(","),
            )),
        }
    }
}

#[allow(clippy::derivable_impls)]
impl Default for VakintSettings {
    fn default() -> Self {
        VakintSettings {
            epsilon_symbol: "ε".into(),
            mu_r_sq_symbol: "mursq".into(),
            form_exe_path: env::var("FORM_PATH").unwrap_or("form".into()),
            python_exe_path: env::var("PYTHON_BIN_PATH").unwrap_or("python3".into()),
            verify_numerator_identification: false,
            run_time_decimal_precision: 32,
            integral_normalization_factor: LoopNormalizationFactor::pySecDec,
            allow_unknown_integrals: true,
            clean_tmp_dir: env::var("VAKINT_NO_CLEAN_TMP_DIR").is_err(),
            evaluation_order: EvaluationOrder::default(),
            // Default to a three-loop UV subtraction problem, for which alphaLoop implementation can be used.
            number_of_terms_in_epsilon_expansion: 4,
            use_dot_product_notation: false,
            temporary_directory: None,
        }
    }
}

#[derive(Debug, Clone)]
struct FullPattern {
    pattern: Pattern,
    conditions: Condition<PatternRestriction>,
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

#[derive(Debug, Clone)]
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

        integral_specs.apply_replacement_rules()?;
        self.apply_numerator_replacement_rules(&integral_specs, &vakint.settings)?;

        let mut could_evaluate_integral = false;
        'eval: for evaluation_approach in vakint.settings.evaluation_order.0.iter() {
            if evaluation_approach.supports(vakint, &integral_specs.canonical_topology) {
                let evaluated_integral = evaluation_approach.evaluate_integral(
                    vakint,
                    self.numerator.as_atom_view(),
                    &integral_specs,
                )?;
                self.numerator = simplify_real(evaluated_integral.as_view());
                self.integral = Atom::new_num(1);
                could_evaluate_integral = true;
                break 'eval;
            }
        }
        if !could_evaluate_integral {
            return Err(VakintError::NoEvaluationMethodFound(
                self.integral.to_string(),
                vakint.settings.number_of_terms_in_epsilon_expansion
                    - (integral_specs.canonical_topology.get_integral().n_loops as i64)
                    - 1,
            ));
        }
        Ok(())
    }

    pub fn apply_numerator_replacement_rules(
        &mut self,
        replacement_rules: &ReplacementRules,
        vakint_settings: &VakintSettings,
    ) -> Result<(), VakintError> {
        let mut new_numerator = self.numerator.clone();
        let mut test = self.numerator.clone();

        // Here we must be careful that substituations can be of the form:
        // {a-> b, b->c}
        // so that we must not apply the substitutions to the outcome of the previous substitution.
        // We will use replace_all_multiple for this
        let casted_replacement_rules = replacement_rules
            .numerator_substitutions
            .iter()
            .map(|(source, target)| {
                (
                    source.clone().into_pattern(),
                    target.clone().as_view().into_pattern().into(),
                )
            })
            .collect::<Vec<_>>();
        let one_substitution_pattern = Atom::parse("1").unwrap().into_pattern().into();
        new_numerator = new_numerator.replace_all_multiple(
            casted_replacement_rules
                .iter()
                .map(|(source, target)| Replacement::new(source, target))
                .collect::<Vec<_>>()
                .as_slice(),
        );
        test = test.replace_all_multiple(
            casted_replacement_rules
                .iter()
                .map(|(source, _target)| Replacement::new(source, &one_substitution_pattern))
                .collect::<Vec<_>>()
                .as_slice(),
        );

        // Make sure to also set all externals to zero for the test
        test = Pattern::parse(format!("{}(momID_,idx_)", EXTERNAL_MOMENTUM_SYMBOL).as_str())
            .unwrap()
            .replace_all(
                test.as_view(),
                &Atom::parse("1").unwrap().into_pattern().into(),
                None,
                None,
            );

        // Substitute metric in as well
        test = Pattern::parse("g(idx1_,idx2_)").unwrap().replace_all(
            test.as_view(),
            &Atom::parse("1").unwrap().into_pattern().into(),
            None,
            None,
        );

        // Substitute epsilon regulator
        test = Pattern::parse(&vakint_settings.epsilon_symbol)
            .unwrap()
            .replace_all(
                test.as_view(),
                &Atom::parse("1").unwrap().into_pattern().into(),
                None,
                None,
            );

        if vakint_settings.verify_numerator_identification && !matches!(test, Atom::Num(_)) {
            return Err(VakintError::NumeratorNotReplaced(
                EXTERNAL_MOMENTUM_SYMBOL.into(),
                test.to_string(),
            ));
        }
        self.numerator = new_numerator;

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

        self.apply_numerator_replacement_rules(replacement_rules, &vakint.settings)?;

        Ok(())
    }

    pub fn identify_vectors_in_numerator(
        numerator: AtomView,
    ) -> Result<Vec<(String, i64)>, VakintError> {
        let mut vectors = HashSet::new();
        // make sure the numerator is in the form of vec_(id_,idx_)
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
        form_numerator = Vakint::convert_from_dot_notation(form_numerator.as_view(), false);
        let vectors = Self::identify_vectors_in_numerator(form_numerator.as_view())?;
        for (vec, id) in vectors.iter() {
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
                    .unwrap()
                    .into(),
                    Some(&Condition::from((
                        State::get_symbol("idx_"),
                        number_condition(),
                    ))),
                    None,
                )
        }

        let template =
            Template::parse_template(TEMPLATES.get("run_tensor_reduction.txt").unwrap()).unwrap();

        // Replace functions with 1 and get all remaining symbols
        let mut numerator_additional_symbols = Pattern::parse("f_(args__)")
            .unwrap()
            .replace_all(
                form_numerator.as_view(),
                &Atom::parse("1").unwrap().into_pattern().into(),
                None,
                None,
            )
            .get_all_symbols(false);
        let eps_symbol = State::get_symbol(vakint.settings.epsilon_symbol.clone());
        numerator_additional_symbols.retain(|&s| s != eps_symbol);

        let mut vars: HashMap<String, String> = HashMap::new();

        vars.insert(
            "numerator".into(),
            vakint.prepare_expression_for_form(form_numerator)?,
        );
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
            &["tensorreduce.frm".into(), "pvtab10.h".into()],
            ("run_tensor_reduction.frm".into(), rendered),
            vec![],
            vakint.settings.clean_tmp_dir,
            vakint.settings.temporary_directory.clone(),
        )?;

        let mut reduced_numerator = vakint.process_form_output(form_result)?;

        for (vec, id) in vectors.iter() {
            reduced_numerator = Pattern::parse(format!("{}{}", vec, id).as_str())
                .unwrap()
                .replace_all(
                    reduced_numerator.as_view(),
                    &Pattern::parse(format!("{}({})", vec, id).as_str())
                        .unwrap()
                        .into(),
                    None,
                    None,
                );
        }

        if !vakint.settings.use_dot_product_notation {
            reduced_numerator =
                Vakint::convert_from_dot_notation(reduced_numerator.as_view(), false);
        }

        self.numerator = reduced_numerator;
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VakintExpression(pub Vec<VakintTerm>);

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

#[derive(Debug, Clone, Default)]
pub struct NumericalEvaluationResult(pub Vec<(i64, Complex<Float>)>);

impl fmt::Display for NumericalEvaluationResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.0.is_empty() {
            write!(f, "{}", "Empty result".green())
        } else {
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
}

impl NumericalEvaluationResult {
    pub fn to_atom(&self, epsilon_symbol: Symbol) -> Atom {
        let mut res = Atom::Zero;
        for (exp, coeff) in self.get_epsilon_coefficients() {
            res = res
                + (Atom::new_num(coeff.re) + Atom::new_var(State::I) * Atom::new_num(coeff.im))
                    * Atom::new_var(epsilon_symbol).pow(&Atom::new_num(exp))
        }
        res
    }

    pub fn from_atom(
        input: AtomView,
        epsilon_symbol: Symbol,
        settings: &VakintSettings,
    ) -> Result<Self, VakintError> {
        let epsilon_coeffs = input.coefficient_list(epsilon_symbol);

        let mut epsilon_coeffs_vec = epsilon_coeffs
            .0
            .iter()
            .map(|(eps_atom, coeff)| {
                if let Some(m) = Pattern::parse(format!("{}^n_", epsilon_symbol).as_str())
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
                } else if *eps_atom == Atom::new_var(epsilon_symbol) {
                    (1, coeff)
                } else {
                    panic!("Epsilon atom should be of the form {}^n_", epsilon_symbol)
                }
            })
            .collect::<Vec<_>>();

        if !epsilon_coeffs.1.is_zero() {
            epsilon_coeffs_vec.push((0, &epsilon_coeffs.1));
        }

        let binary_prec = settings.get_binary_precision();
        let mut epsilon_coeffs_vec_floats = vec![];
        for (i64, coeff) in epsilon_coeffs_vec.iter() {
            epsilon_coeffs_vec_floats.push((
                *i64,
                match coeff.evaluate(
                    &Vakint::get_coeff_map(binary_prec),
                    &HashMap::default(),
                    &HashMap::default(),
                    &mut HashMap::default(),
                ) {
                    Ok(x) => x,
                    Err(e) => {
                        return Err(VakintError::EvaluationError(format!(
                            "Coefficients of the Laurent series are not numbers. Expression: '{}' | Error: '{}'",
                            coeff.to_canonical_string(),
                            e
                        )));
                    }
                },
            ));
        }

        epsilon_coeffs_vec_floats.sort_by(|(i1, _), (i2, _)| i1.cmp(i2));
        Ok(NumericalEvaluationResult(epsilon_coeffs_vec_floats))
    }

    pub fn from_vec(input: Vec<(i64, (String, String))>, settings: &VakintSettings) -> Self {
        let binary_prec = settings.get_binary_precision();
        NumericalEvaluationResult(
            input
                .iter()
                .map(|(eps_pwr, (re, im))| {
                    (
                        *eps_pwr,
                        Complex::new(
                            Float::parse(re.as_str(), Some(binary_prec)).unwrap(),
                            Float::parse(im.as_str(), Some(binary_prec)).unwrap(),
                        ),
                    )
                })
                .collect::<Vec<_>>(),
        )
    }

    pub fn get_epsilon_coefficients(&self) -> Vec<(i64, Complex<Float>)> {
        self.0.clone()
    }
    pub fn get_epsilon_coefficient(&self, power: i64) -> Complex<Float> {
        match self.0.iter().find(|(i, _f)| *i == power) {
            Some((_, f)) => f.clone(),
            None => {
                if self.0.is_empty() {
                    Complex::new(Float::with_val(53, 0.), Float::with_val(53, 0.))
                } else {
                    Complex::new(self.0[0].1.re.zero(), self.0[0].1.re.zero())
                }
            }
        }
    }

    pub fn is_zero(&self) -> bool {
        self.0.iter().all(|(_, f)| f.is_zero())
    }

    pub fn aggregate_errors(&self, other: &NumericalEvaluationResult) -> NumericalEvaluationResult {
        let mut res = NumericalEvaluationResult::default();
        let mut orders = vec![];
        for (o, _val) in self.0.iter() {
            orders.push(*o);
        }
        for (o, _val) in other.0.iter() {
            orders.push(*o);
        }
        orders.sort();
        for o in orders {
            let (a, b) = (
                &self.get_epsilon_coefficient(o),
                &other.get_epsilon_coefficient(o),
            );

            res.0.push((0, (a * a + b * b).sqrt()));
        }
        res
    }

    pub fn does_approx_match(
        &self,
        other: &NumericalEvaluationResult,
        error: Option<&NumericalEvaluationResult>,
        threshold: f64,
        max_pull: f64,
    ) -> (bool, String) {
        let mut powers = HashSet::new();
        for (power, _) in self.0.iter() {
            powers.insert(*power);
        }
        for (power, _) in other.0.iter() {
            powers.insert(*power);
        }
        let mut power_vec = powers.iter().cloned().collect::<Vec<_>>();
        power_vec.sort();

        for power in power_vec {
            let (self_val, other_val) = (
                self.get_epsilon_coefficient(power),
                other.get_epsilon_coefficient(power),
            );
            let comparisons = vec![
                (
                    "real",
                    self_val.re.clone(),
                    other_val.re.clone(),
                    error.map(|e| e.get_epsilon_coefficient(power).re.clone()),
                ),
                (
                    "imaginary",
                    self_val.im.clone(),
                    other_val.im.clone(),
                    error.map(|e| e.get_epsilon_coefficient(power).im.clone()),
                ),
            ];
            for (part, r, o, e) in comparisons.iter() {
                let f_max_pull = Float::with_val(r.prec(), max_pull);
                let f_threshold = Float::with_val(r.prec(), threshold);
                let delta = (r.clone() - o).norm();
                if let Some(err) = e {
                    if delta > f_max_pull * err.norm() {
                        return (
                            false,
                            format!(
                                "{} part of ε^{} coefficient does not match within max pull: {} != {} (pull = {})",
                                part, power, r, o, delta/err.norm()
                            ),
                        );
                    }
                }
                let scale = (r.norm() + o.norm()) / Float::with_val(r.prec(), 2.0);
                if scale.is_zero() {
                    if delta > f_threshold {
                        return (
                            false,
                            format!(
                                "{} part of ε^{} coefficient does not match within abs. error required: {} != {} (abs. error = {})",
                                part, power, r, o, delta
                            ),
                        );
                    }
                } else {
                    let rel_error = delta.div(scale);
                    if rel_error > f_threshold {
                        return (
                            false,
                            format!(
                                "{} part of ε^{} coefficient does not match within rel. error required: {} != {} (rel. error = {})",
                                part, power, r, o, rel_error
                            ),
                        );
                    }
                }
            }
        }
        (true, "matches!".into())
    }
}

impl Vakint {
    pub fn new(settings: Option<VakintSettings>) -> Result<Self, VakintError> {
        // Force initialization of symbolica symbols with proper attributes
        LazyLock::force(&S);
        let vakint_settings = settings.unwrap_or_default();

        // Verify that the chosen normalisation only contains the expected symbols
        let (_full_atom, expanded, evaluated) = vakint_settings
            .integral_normalization_factor
            .validate(&vakint_settings)?;

        debug!(
            "Loop normalisation factor considered:\nFull                          : {}\nExpanded (n_loops=1)          : {}\nEvaluated (n_loops=1, mu_r=1) :\n{}",
            vakint_settings.integral_normalization_factor, expanded, evaluated
        );

        let topologies = Topologies::generate_topologies(&vakint_settings)?;
        let vakint = Vakint {
            settings: vakint_settings,
            topologies,
        };

        if vakint
            .settings
            .evaluation_order
            .0
            .iter()
            .any(|em| em.dependencies().contains(&VakintDependency::FORM))
        {
            let form_version = vakint.get_form_version()?;
            match compare_to(form_version.clone(), MINIMAL_FORM_VERSION, Cmp::Ge) {
                Ok(valid) => {
                    if valid {
                        debug!(
                            "{} successfully detected with version '{}'",
                            "FORM".green(),
                            form_version.green()
                        );
                    } else {
                        return Err(VakintError::FormVersion(format!(
                        "{} version installed on your system does not meet minimal requirements: {}<{}",
                        "FORM".red(),
                        form_version.red(), MINIMAL_FORM_VERSION
                    )));
                    }
                }
                Err(_) => {
                    return Err(VakintError::FormVersion(format!(
                        "Could not parse {} version '{}'.",
                        "FORM".red(),
                        form_version.red()
                    )))
                }
            };
        }

        if vakint
            .settings
            .evaluation_order
            .0
            .iter()
            .any(|em| em.dependencies().contains(&VakintDependency::PySecDec))
        {
            let pysecdec_version = vakint.get_pysecdec_version()?;
            match compare_to(pysecdec_version.clone(), MINIMAL_PYSECDEC_VERSION, Cmp::Ge) {
                Ok(valid) => {
                    if valid {
                        debug!(
                            "{} successfully detected with version '{}'.",
                            "PySecDec".green(),
                            pysecdec_version.green()
                        );
                    } else {
                        return Err(VakintError::FormVersion(format!(
                            "{} version installed on your system does not meet minimal requirements: {}<{}.",
                            "PySecDec".red(),
                            pysecdec_version.red(), MINIMAL_PYSECDEC_VERSION
                        )));
                    }
                }
                Err(_) => {
                    return Err(VakintError::FormVersion(format!(
                        "Could not parse {} version '{}'.",
                        "PySecDec".red(),
                        pysecdec_version.red()
                    )))
                }
            };
        }
        //println!("Topologies generated:\n{}", vakint.topologies);

        Ok(vakint)
    }

    pub fn params_from_f64(
        &self,
        params: &HashMap<String, f64>,
    ) -> HashMap<String, Complex<Float>, ahash::RandomState> {
        let hm = HashMap::<String, f64, ahash::RandomState>::from_iter(
            params.iter().map(|(k, v)| (k.clone(), *v)),
        );
        params_from_f64(&hm, self.settings.run_time_decimal_precision)
    }

    pub fn externals_from_f64(
        &self,
        externals: &HashMap<usize, (f64, f64, f64, f64)>,
    ) -> HashMap<usize, Momentum, ahash::RandomState> {
        let em = HashMap::<usize, (f64, f64, f64, f64), ahash::RandomState>::from_iter(
            externals.iter().map(|(k, v)| (*k, *v)),
        );
        externals_from_f64(&em, self.settings.run_time_decimal_precision)
    }

    pub fn numerical_evaluation(
        &self,
        expression: AtomView,
        params: &HashMap<String, Complex<Float>, ahash::RandomState>,
        externals: Option<&HashMap<usize, Momentum, ahash::RandomState>>,
    ) -> Result<(NumericalEvaluationResult, Option<NumericalEvaluationResult>), VakintError> {
        Vakint::full_numerical_evaluation(&self.settings, expression, params, externals)
    }

    pub fn get_coeff_map(prec: u32) -> impl Fn(&Fraction<IntegerRing>) -> Complex<Float> {
        move |x| Complex::new(x.to_multi_prec_float(prec), Float::with_val(prec, 0.))
    }

    pub fn get_constants_map(
        settings: &VakintSettings,
        params: &HashMap<String, Complex<Float>, ahash::RandomState>,
        externals: Option<&HashMap<usize, Momentum, ahash::RandomState>>,
    ) -> Result<HashMap<Atom, Complex<Float>, ahash::random_state::RandomState>, VakintError> {
        let mut const_map: HashMap<Atom, Complex<Float>, ahash::random_state::RandomState> =
            HashMap::default();

        let binary_prec = settings.get_binary_precision();
        const_map.insert(
            Atom::from(Var::new(State::get_symbol("𝜋"))),
            Complex::new(
                Float::with_val(binary_prec, Constant::Pi),
                Float::with_val(binary_prec, 0),
            ),
        );

        if let Some(ext_p) = externals.as_ref() {
            for (&i, val1) in ext_p.iter() {
                for (&j, val2) in ext_p.iter() {
                    if i > j {
                        continue;
                    }
                    let dot_product = val1.0.to_owned() * val2.0.to_owned()
                        - val1.1.to_owned() * val2.1.to_owned()
                        - val1.2.to_owned() * val2.2.to_owned()
                        - val1.3.to_owned() * val2.3.to_owned();
                    const_map.insert(
                        fun!(
                            S.dot,
                            fun!(S.p, Atom::new_num(i as i64)),
                            fun!(S.p, Atom::new_num(j as i64))
                        ),
                        dot_product.to_owned(),
                    );
                }
            }
        }

        const_map.insert(
            Atom::from(Var::new(State::get_symbol("EulerGamma"))),
            Complex::new(
                Float::with_val(binary_prec, Constant::Euler),
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
            const_map.insert(Atom::parse(symb.as_str()).unwrap(), value.clone());
        }

        Ok(const_map)
    }

    pub fn partial_numerical_evaluation(
        settings: &VakintSettings,
        integral: AtomView,
        params: &HashMap<String, Complex<Float>, ahash::RandomState>,
        externals: Option<&HashMap<usize, Momentum, ahash::RandomState>>,
    ) -> Atom {
        let const_map = Vakint::get_constants_map(settings, params, externals).unwrap();

        let mut res = Vakint::convert_to_dot_notation(integral);
        for (src, trgt) in const_map.iter() {
            res = src.into_pattern().replace_all(
                res.as_view(),
                &((Atom::new_num(trgt.re.clone())
                    + Atom::new_var(State::get_symbol("𝑖")) * Atom::new_num(trgt.im.clone()))
                .into_pattern())
                .into(),
                None,
                None,
            );
        }

        res
    }

    pub fn full_numerical_evaluation(
        settings: &VakintSettings,
        integral: AtomView,
        params: &HashMap<String, Complex<Float>, ahash::RandomState>,
        externals: Option<&HashMap<usize, Momentum, ahash::RandomState>>,
    ) -> Result<(NumericalEvaluationResult, Option<NumericalEvaluationResult>), VakintError> {
        if S.error_flag
            .into_pattern()
            .pattern_match(integral, &Condition::default(), &MatchSettings::default())
            .next()
            .is_none()
        {
            return Ok((
                Vakint::full_numerical_evaluation_without_error(
                    settings, integral, params, externals,
                )?,
                None,
            ));
        }
        let (error_atom, integral_atom) = integral.coefficient_list(S.error_flag_symbol);
        assert!(error_atom.len() == 1);
        let mut central = Vakint::full_numerical_evaluation_without_error(
            settings,
            integral_atom.as_view(),
            params,
            externals,
        )?;
        let mut error = Vakint::full_numerical_evaluation_without_error(
            settings,
            error_atom[0].1.as_view(),
            params,
            externals,
        )?;

        for (i, eval) in central.0.iter() {
            if !error.0.iter().any(|(j, _)| i == j) {
                error
                    .0
                    .push((*i, Complex::new(eval.re.zero(), eval.re.zero())));
            }
        }
        for (i, eval) in error.0.iter() {
            if !central.0.iter().any(|(j, _)| i == j) {
                central
                    .0
                    .push((*i, Complex::new(eval.re.zero(), eval.re.zero())));
            }
        }
        Ok((central, Some(error)))
    }

    pub fn full_numerical_evaluation_without_error(
        settings: &VakintSettings,
        integral: AtomView,
        params: &HashMap<String, Complex<Float>, ahash::RandomState>,
        externals: Option<&HashMap<usize, Momentum, ahash::RandomState>>,
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

        if !epsilon_coeffs.1.is_zero() {
            epsilon_coeffs_vec.push((0, &epsilon_coeffs.1));
        }

        let map = Vakint::get_constants_map(settings, params, externals).unwrap();
        let map_view: HashMap<AtomView<'_>, Complex<Float>, RandomState> =
            map.iter().map(|(k, v)| (k.as_view(), v.clone())).collect();

        let binary_prec = settings.get_binary_precision();
        let mut epsilon_coeffs_vec_floats = vec![];
        for (i64, coeff) in epsilon_coeffs_vec.iter() {
            let coeff_processed = Vakint::convert_to_dot_notation(coeff.as_view());
            // for (k, v) in map_view.iter() {
            //     println!("{} -> {}", k, v);
            // }
            // println!("coeff_processed={}", coeff_processed);

            epsilon_coeffs_vec_floats.push((
                *i64,
                match coeff_processed.evaluate(
                    &Vakint::get_coeff_map(binary_prec),
                    &map_view,
                    &HashMap::default(),
                    &mut HashMap::default(),
                ) {
                    Ok(x) => x,
                    Err(e) => {
                        return Err(VakintError::EvaluationError(format!(
                            "Is some tensor structure left? | Expression: '{}' | Error: '{}'",
                            coeff.to_canonical_string(),
                            e
                        )));
                    }
                },
            ));
        }

        epsilon_coeffs_vec_floats.sort_by(|(i1, _), (i2, _)| i1.cmp(i2));
        Ok(NumericalEvaluationResult(epsilon_coeffs_vec_floats))
    }

    fn pysecdec_evaluate(
        &self,
        vakint: &Vakint,
        input_numerator: AtomView,
        integral_specs: &ReplacementRules,
        options: &PySecDecOptions,
    ) -> Result<Atom, VakintError> {
        let integral = integral_specs.canonical_topology.get_integral();

        let pysecdec_inputs = if options.reuse_existing_output.is_some()
            && PathBuf::from(options.reuse_existing_output.as_ref().unwrap()).exists()
        {
            vec![("run_pySecDec.py".into(), "".into())]
        } else {
            debug!(
                "Processing the following integral with {}:\n{}",
                "PySecDec".green(),
                integral
            );
            let dot_product_numerator = Vakint::convert_from_dot_notation(input_numerator, false);
            let vectors =
                VakintTerm::identify_vectors_in_numerator(dot_product_numerator.as_view())?;
            let numerator_atom = Vakint::convert_to_dot_notation(input_numerator);
            let numerator = numerator_atom.as_view();
            let mut processed_numerator = Vakint::convert_to_dot_notation(numerator);

            // Make sure there is no open index left
            if Pattern::parse("s_(id_,idx_)")
                .unwrap()
                .pattern_match(
                    processed_numerator.as_view(),
                    &(Condition::from((State::get_symbol("id_"), number_condition()))
                        & Condition::from((State::get_symbol("idx_"), number_condition()))
                        & Condition::from((State::get_symbol("s_"), symbol_condition()))),
                    &MatchSettings::default(),
                )
                .next()
                .is_some()
            {
                return Err(VakintError::InvalidNumerator(
                format!("PySecDec can only handle scalar numerator. If you have open indices, make sure they are contracted with external momenta: {}",processed_numerator),
            ));
            }

            // Check if numerator contains additional symbols
            // First, replace functions with 1 and get all remaining symbols
            let mut numerator_additional_symbols = Pattern::parse("f_(args__)")
                .unwrap()
                .replace_all(
                    input_numerator,
                    &Atom::parse("1").unwrap().into_pattern().into(),
                    None,
                    None,
                )
                .get_all_symbols(false);
            let eps_symbol: Symbol = State::get_symbol(vakint.settings.epsilon_symbol.clone());
            numerator_additional_symbols.retain(|&s| s != eps_symbol);

            // Convert back from dot notation
            processed_numerator =
                Vakint::convert_from_dot_notation(processed_numerator.as_view(), true);

            let mut lorentz_indices = HashSet::<String>::new();
            let arc_mutex_lorentz_indices = Arc::new(Mutex::new(HashSet::new()));

            let dot_product_matcher = Pattern::parse("v1_(id1_,idx_)*v2_(id2_,idx_)").unwrap();
            // Powers higher than two cannot occur as different dummy indices would have been used in
            // the call 'processed_numerator = Vakint::convert_from_dot_notation(processed_numerator.as_view(), true)'
            let square_matcher = Pattern::parse("v1_(id1_,idx_)^2").unwrap();
            let mut old_processed_numerator = processed_numerator.clone();
            loop {
                let arc_mutex_lorentz_indices_sent = arc_mutex_lorentz_indices.clone();
                let dot_product_transformer =
                    Transformer::Map(Box::new(move |a_in: AtomView, a_out: &mut Atom| {
                        if let AtomView::Fun(s) = a_in {
                            let a_in = s.to_slice();
                            arc_mutex_lorentz_indices_sent
                                .lock()
                                .unwrap()
                                .insert(format!("mu{}", a_in.get(4).to_canonical_string()));
                            *a_out = Atom::parse(
                                format!(
                                    "{}{}(mu{})*{}{}(mu{})",
                                    a_in.get(0).to_canonical_string(),
                                    a_in.get(1).to_canonical_string(),
                                    a_in.get(4).to_canonical_string(),
                                    a_in.get(2).to_canonical_string(),
                                    a_in.get(3).to_canonical_string(),
                                    a_in.get(4).to_canonical_string(),
                                )
                                .as_str(),
                            )
                            .unwrap();
                        };

                        Ok(())
                    }));

                processed_numerator = dot_product_matcher.replace_all(
                    old_processed_numerator.as_view(),
                    &Pattern::Transformer(Box::new((
                        Some(Pattern::parse("arg(v1_,id1_,v2_,id2_,idx_)").unwrap()),
                        vec![dot_product_transformer.clone()],
                    )))
                    .into(),
                    None,
                    Some(&MatchSettings {
                        rhs_cache_size: 0,
                        ..Default::default()
                    }),
                );

                processed_numerator = square_matcher.replace_all(
                    processed_numerator.as_view(),
                    &Pattern::Transformer(Box::new((
                        Some(Pattern::parse("arg(v1_,id1_,v1_,id1_,idx_)").unwrap()),
                        vec![dot_product_transformer],
                    )))
                    .into(),
                    None,
                    Some(&MatchSettings {
                        rhs_cache_size: 0,
                        ..Default::default()
                    }),
                );
                lorentz_indices.extend(arc_mutex_lorentz_indices.lock().unwrap().clone());
                if old_processed_numerator == processed_numerator {
                    break;
                } else {
                    old_processed_numerator = processed_numerator.clone();
                }
            }

            let mut m = Atom::new();
            let power_list_map = integral_specs.get_propagator_property_list("pow_");

            let mass_list_map = integral_specs.get_propagator_property_list("mUVsq_");
            let mut masses = HashSet::new();
            let mut power_list: Vec<i64> = vec![];

            let mut loop_momenta_ids_found = HashSet::new();
            let pysecdec_propagators = format!(
                "[{}]",
                integral
                    .graph
                    .edges
                    .iter()
                    .map(|(e_id, e)| format!(
                        "'({})**2{}'",
                        {
                            m = e.momentum.clone();
                            power_list.push(
                                power_list_map
                                    .get(e_id)
                                    .unwrap()
                                    .to_owned()
                                    .try_into()
                                    .unwrap(),
                            );
                            for i_loop in 1..=integral.n_loops {
                                let p = Pattern::parse(format!("k({})", i_loop).as_str()).unwrap();
                                if utils::could_match(&p, m.as_view()) {
                                    loop_momenta_ids_found.insert(i_loop);
                                }
                                m = p.replace_all(
                                    m.as_view(),
                                    &Pattern::parse(format!("k{}", i_loop).as_str())
                                        .unwrap()
                                        .into(),
                                    None,
                                    None,
                                );
                            }
                            AtomPrinter::new_with_options(m.as_view(), PrintOptions::file())
                        },
                        {
                            let m = mass_list_map.get(e_id).unwrap().to_owned();
                            if m.is_zero() {
                                "".into()
                            } else {
                                masses.insert(m.get_symbol().unwrap().to_string());
                                format!("-{}", m.get_symbol().unwrap())
                            }
                        }
                    ))
                    .collect::<Vec<_>>()
                    .join(", "),
            );
            let n_loops_in_topology = loop_momenta_ids_found.len();

            // This should always hold, because pinched higher loops that reduce the loop counts must be captured by lower loop count topologies.
            // It is important that this hold otherwise MSbar normalization factor would be off.
            assert!(n_loops_in_topology == integral.n_loops);

            let mut vars: HashMap<String, String> = HashMap::new();

            vars.insert("graph_name".into(), "pySecDecRun".into());

            vars.insert("propagators".into(), pysecdec_propagators);
            let mut sorted_lorentz_indices = lorentz_indices.iter().cloned().collect::<Vec<_>>();
            sorted_lorentz_indices.sort();
            vars.insert(
                "lorentz_indices".into(),
                format!(
                    "[{}]",
                    sorted_lorentz_indices
                        .iter()
                        .map(|li| format!("'{}'", li))
                        .collect::<Vec<_>>()
                        .join(",")
                ),
            );
            vars.insert(
                "loop_momenta".into(),
                format!(
                    "[{}]",
                    (1..=integral.n_loops)
                        .map(|i| format!("'k{}'", i))
                        .collect::<Vec<_>>()
                        .join(",")
                ),
            );
            let mut external_momenta = vectors
                .iter()
                .filter_map(|(v, id)| {
                    if v == EXTERNAL_MOMENTUM_SYMBOL {
                        Some(format!("{}{}", v, id))
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>();
            external_momenta.sort();
            vars.insert(
                "external_momenta".into(),
                format!(
                    "[{}]",
                    external_momenta
                        .iter()
                        .map(|em| format!("'{}'", em))
                        .collect::<Vec<_>>()
                        .join(",")
                ),
            );
            vars.insert(
                "power_list".into(),
                format!(
                    "tuple([{}])",
                    power_list
                        .iter()
                        .map(|p| format!("{}", p))
                        .collect::<Vec<_>>()
                        .join(",")
                ),
            );
            let numerator_path = String::from("numerator.txt");
            vars.insert("numerator_path".into(), numerator_path.clone());
            // Expand the numerator around epsilon=0 to make sure it is polynomial
            processed_numerator = processed_numerator
                .series(
                    State::get_symbol(vakint.settings.epsilon_symbol.as_str()),
                    Atom::Zero.as_atom_view(),
                    Rational::from(
                        vakint.settings.number_of_terms_in_epsilon_expansion
                            - (integral.n_loops as i64),
                    ),
                    true,
                )
                .unwrap()
                .to_atom();
            let mut numerator_string = processed_numerator
                .to_canonical_string()
                .replace(vakint.settings.epsilon_symbol.as_str(), "eps");
            // Powers higher than two cannot occur as different dummy indices would have been used in
            // the call 'processed_numerator = Vakint::convert_from_dot_notation(processed_numerator.as_view(), true)'
            let remove_squares_re =
                Regex::new(r"(?<vec>[\w|-|\d]+)\((?<idx>mu-?\d+)\)\^2").unwrap();
            numerator_string = remove_squares_re
                .replace_all(numerator_string.as_str(), "$vec$id($idx)*$vec$id($idx)")
                .to_string();

            let description_path = "description.txt".into();
            let description_string = format!(
                "Integral:\n{}\nNumerator:\n{}",
                integral_specs,
                numerator.to_owned()
            );

            vars.insert("n_loops".into(), format!("{}", n_loops_in_topology));
            vars.insert("loop_additional_prefactor".into(), "1".into());
            vars.insert("additional_overall_factor".into(), "1.0".into());

            let mut masses_vec = masses.iter().cloned().collect::<Vec<_>>();
            masses_vec.sort();
            let mut real_parameters = vec![];
            let mut replacement_rules = vec![];

            // Add masses
            for m in masses_vec.iter() {
                real_parameters.push(format!("'{}'", m.clone()));
            }

            //Potential extra parameters
            let mut sorted_additional_numerator_symbols =
                numerator_additional_symbols.iter().collect::<Vec<_>>();
            sorted_additional_numerator_symbols.sort();
            for additional_param in sorted_additional_numerator_symbols.iter() {
                real_parameters.push(format!("'{}'", additional_param));
            }

            // And the external momenta
            for iv1 in 0..external_momenta.len() {
                for iv2 in iv1..external_momenta.len() {
                    let v1 = external_momenta[iv1].clone();
                    let v2 = external_momenta[iv2].clone();
                    replacement_rules.push((format!("'{}*{}'", v1, v2), format!("'{}{}'", v1, v2)));
                    real_parameters.push(format!("'{}{}'", v1, v2));
                }
            }

            vars.insert(
                "replacement_rules".into(),
                format!(
                    "[{}]",
                    replacement_rules
                        .iter()
                        .map(|(k, v)| format!("({},{})", k, v))
                        .collect::<Vec<_>>()
                        .join(",\n")
                ),
            );

            vars.insert(
                "real_parameters".into(),
                format!("[{}]", real_parameters.join(",")),
            );

            vars.insert("complex_parameters_input".into(), "[]".into());
            vars.insert("real_parameters_input".into(), "[]".into());
            let mut default_external_momenta = vec![];
            for p in external_momenta {
                if let Some(f) = options.numerical_external_momenta.get(&p) {
                    default_external_momenta.push((p, f));
                } else {
                    return Err(VakintError::EvaluationError(format!("Missing specification of numerical value for external momentum '{}'. Specify it in the PySecDecOptions of Vakint.", p)));
                }
            }
            vars.insert(
                "default_externals".into(),
                format!(
                    r"[{}\n]",
                    default_external_momenta
                        .iter()
                        .map(|(pi, f)| format!(
                            "({:.e},{:.e},{:.e},{:.e}) #{}",
                            f.0, f.1, f.2, f.3, pi
                        ))
                        .collect::<Vec<_>>()
                        .join(r"\n,")
                ),
            );
            let mut default_masses = vec![];
            for m in masses_vec {
                if let Some(num_m) = options.numerical_masses.get(&m) {
                    default_masses.push((m, num_m));
                } else {
                    return Err(VakintError::EvaluationError(format!("Missing specification of numerical value for mass '{}'. Specify it in the PySecDecOptions of Vakint.", m)));
                }
            }
            for additional_param in sorted_additional_numerator_symbols.iter() {
                if let Some(num_additional_param) =
                    options.numerical_masses.get(&additional_param.to_string())
                {
                    default_masses.push((additional_param.to_string(), num_additional_param));
                } else {
                    return Err(VakintError::EvaluationError(format!("Missing specification of numerical value for additional numerator symbol '{}'. Specify it in the PySecDecOptions of Vakint.", additional_param)));
                }
            }
            vars.insert(
                "default_masses".into(),
                format!(
                    r"[{}\n]",
                    default_masses
                        .iter()
                        .map(|(m, num_m)| format!("{:.e} #{}", num_m, m))
                        .collect::<Vec<_>>()
                        .join(r"\n,")
                ),
            );

            vars.insert(
                "max_epsilon_order".into(),
                format!(
                    "{}",
                    vakint.settings.number_of_terms_in_epsilon_expansion
                        - (integral.n_loops as i64)
                        - 1
                ),
            );
            vars.insert("contour_deformation".into(), "False".into());

            vars.insert(
                "drawing_input_power_list".into(),
                vars.get("power_list").unwrap().clone(),
            );
            vars.insert("drawing_input_external_lines".into(), "[]".into());
            let internal_lines_str = format!(
                "[{}]",
                integral
                    .graph
                    .edges
                    .iter()
                    .map(|(e_id, e)| {
                        format!("['e{}',[{},{}]]", e_id, e.left_node_id, e.right_node_id)
                    })
                    .collect::<Vec<_>>()
                    .join(",")
            );
            vars.insert("drawing_input_internal_lines".into(), internal_lines_str);

            vars.insert("couplings_prefactor".into(), "'1'".into());
            vars.insert("couplings_values".into(), "{}".into());

            let template =
                Template::parse_template(TEMPLATES.get("run_pySecDec_template.txt").unwrap())
                    .unwrap();

            let rendered = template
                .render(&RenderOptions {
                    variables: vars,
                    ..Default::default()
                })
                .unwrap_or_else(|e| {
                    panic!("Error while rendering template: {}", e);
                });

            vec![
                ("run_pySecDec.py".into(), rendered),
                (numerator_path, numerator_string),
                (description_path, description_string),
            ]
        };

        let mut pysecdec_options = vec!["-i".into(), "qmc".into()];
        pysecdec_options.push("--eps_rel".into());
        pysecdec_options.push(format!("{:.e}", options.relative_precision));
        pysecdec_options.push("--min_n".into());
        pysecdec_options.push(format!("{}", options.min_n_evals));
        pysecdec_options.push("--max_eval".into());
        pysecdec_options.push(format!("{}", options.max_n_evals));
        if options.quiet {
            pysecdec_options.push("-q".into());
        }
        let pysecdec_output = vakint.run_pysecdec(
            &pysecdec_inputs,
            pysecdec_options,
            vakint.settings.clean_tmp_dir && options.reuse_existing_output.is_none(),
            options.reuse_existing_output.clone(),
            vakint.settings.temporary_directory.clone(),
        )?;

        if !options.quiet {
            info!(
                "PySecDec raw output:\ncentral: {}\nerror : {}",
                pysecdec_output[0], pysecdec_output[1]
            );
        }

        /*
        let pysecdec_output = [
            String::from("ep^(-1)*(2.0000000000000000e+00+(0.0000000000000000e+00)*I)+ep^(0)*(-5.4072569092295630e-01+(0.0000000000000000e+00)*I)+ep^(1)*(2.7180301350542537e+00+(0.0000000000000000e+00)*I)+ep^(2)*(-8.5638398947489147e-01+(0.0000000000000000e+00)*I)"),
            String::from("ep^(-1)*(7.9539672483176226e-17+(0.0000000000000000e+00)*I)+ep^(0)*(5.9466474629934134e-17+(0.0000000000000000e+00)*I)+ep^(1)*(1.1644579020645423e-16+(0.0000000000000000e+00)*I)+ep^(2)*(8.0813279157421749e-17+(0.0000000000000000e+00)*I)")
        ];
        */

        let log_mu_sq = fun!(
            State::LOG,
            Atom::new_var(State::get_symbol(vakint.settings.mu_r_sq_symbol.as_str()))
        );

        let pysecdec_normalization_correction = Atom::parse(
            format!(
                "(  𝑖*(𝜋^((4-2*{eps})/2))\
                )^{n_loops}",
                eps = self.settings.epsilon_symbol,
                n_loops = integral.n_loops
            )
            .as_str(),
        )
        .unwrap();

        let evaluated_integral = pysecdec_output
            .iter()
            .map(|out| match Atom::parse(out.replace("I", "𝑖").as_str()) {
                Ok(mut processed) => {
                    processed = Pattern::parse("ep").unwrap().replace_all(
                        processed.as_view(),
                        &Pattern::parse(&self.settings.epsilon_symbol)
                            .unwrap()
                            .into(),
                        None,
                        None,
                    );
                    processed = processed
                        * pysecdec_normalization_correction.to_owned()
                        * S.n_loops.into_pattern().replace_all(
                            vakint
                                .settings
                                .get_integral_normalization_factor_atom()?
                                .as_view(),
                            &Atom::new_num(integral.n_loops as i64).into_pattern().into(),
                            None,
                            None,
                        );

                    let expanded_evaluation = match processed.series(
                        State::get_symbol(vakint.settings.epsilon_symbol.as_str()),
                        Atom::Zero.as_atom_view(),
                        Rational::from(
                            vakint.settings.number_of_terms_in_epsilon_expansion
                                - (integral.n_loops as i64)
                                - 1,
                        ),
                        true,
                    ) {
                        Ok(a) => a,
                        Err(e) => return Err(VakintError::SymbolicaError(e.to_string())),
                    };

                    let mut evaluated_i = expanded_evaluation.to_atom();
                    evaluated_i = Pattern::parse("log_mu_sq").unwrap().replace_all(
                        evaluated_i.as_view(),
                        &(log_mu_sq).into_pattern().into(),
                        None,
                        None,
                    );

                    Ok(evaluated_i)
                }
                Err(err) => Err(VakintError::PySecDecOutputParsingError(out.clone(), err)),
            })
            .collect::<Result<Vec<_>, VakintError>>()?;

        Ok(evaluated_integral[0].to_owned()
            + evaluated_integral[1].to_owned() * S.error_flag.to_owned())
    }

    fn alphaloop_evaluate(
        &self,
        vakint: &Vakint,
        numerator: AtomView,
        integral_specs: &ReplacementRules,
    ) -> Result<Atom, VakintError> {
        let integral = integral_specs.canonical_topology.get_integral();

        debug!(
            "Processing the following integral with {}:\n{}",
            "AlphaLoop".green(),
            integral
        );

        let muv_sq_symbol = if let Some(m) = Pattern::parse("prop(args__,m_,pow_)")
            .unwrap()
            .pattern_match(
                integral.canonical_expression.as_ref().unwrap().as_view(),
                &Condition::default(),
                &MatchSettings::default(),
            )
            .next()
        {
            match m
                .match_stack
                .get(State::get_symbol("m_"))
                .unwrap()
                .to_atom()
            {
                Atom::Var(s) => s.get_symbol(),
                _ => {
                    return Err(VakintError::MalformedGraph(format!(
                        "Could not find muV in graph:\n{}",
                        integral.canonical_expression.as_ref().unwrap()
                    )));
                }
            }
        } else {
            return Err(VakintError::MalformedGraph(format!(
                "Could not find muV in graph:\n{}",
                integral.canonical_expression.as_ref().unwrap()
            )));
        };

        let alphaloop_expression = integral.alphaloop_expression.as_ref().unwrap().as_view();

        let mut numerator = numerator.to_owned();
        numerator = Atom::new_var(muv_sq_symbol).into_pattern().replace_all(
            numerator.as_view(),
            &Pattern::parse("mUV^2").unwrap().into(),
            None,
            None,
        );
        // println!("Numerator : {}", numerator);
        // println!("Evaluating AlphaLoop : {}", alphaloop_expression);
        // println!("Graph:\n{}", integral.graph.to_graphviz());

        // Make sure to undo the dot product notation.
        // If it was not used, the command below will do nothing.
        let mut form_expression = numerator * alphaloop_expression.to_owned();
        form_expression = Vakint::convert_to_dot_notation(form_expression.as_view());
        // println!("Input expression with dot products : {}", form_expression);

        let mut vector_mapping: HashMap<Atom, Atom> = HashMap::new();
        while let Some(m) = Pattern::parse("v_(id_)")
            .unwrap()
            .pattern_match(
                form_expression.as_view(),
                &(Condition::from((State::get_symbol("id_"), number_condition()))
                    & Condition::from((State::get_symbol("v_"), symbol_condition()))),
                &MatchSettings::default(),
            )
            .next()
        {
            let v = match m
                .match_stack
                .get(State::get_symbol("v_"))
                .unwrap()
                .to_atom()
            {
                Atom::Var(s) => s.get_symbol(),
                _ => {
                    return Err(VakintError::MalformedGraph(format!(
                        "Could not find v in graph:\n{}",
                        integral.canonical_expression.as_ref().unwrap()
                    )))
                }
            };
            let v_id = get_integer_from_match(m.match_stack.get(State::get_symbol("id_")).unwrap())
                .unwrap();

            form_expression = Pattern::parse(format!("{}({})", v, v_id).as_str())
                .unwrap()
                .replace_all(
                    form_expression.as_view(),
                    &Pattern::parse(format!("{}{}", v, v_id).as_str())
                        .unwrap()
                        .into(),
                    None,
                    None,
                );
            vector_mapping.insert(
                Atom::parse(format!("{}{}", v, v_id).as_str()).unwrap(),
                Atom::parse(format!("{}({})", v, v_id).as_str()).unwrap(),
            );
        }
        // println!("Input expression for FORM : {}", form_expression);

        let template = Template::parse_template(
            TEMPLATES
                .get("run_alphaloop_integral_evaluation.txt")
                .unwrap(),
        )
        .unwrap();

        // Replace functions with 1 and get all remaining symbols
        let mut numerator_additional_symbols = Pattern::parse("f_(args__)")
            .unwrap()
            .replace_all(
                form_expression.as_view(),
                &Atom::parse("1").unwrap().into_pattern().into(),
                None,
                None,
            )
            .get_all_symbols(false);
        let eps_symbol = State::get_symbol(vakint.settings.epsilon_symbol.clone());
        numerator_additional_symbols.retain(|&s| s != eps_symbol);

        let mut vars: HashMap<String, String> = HashMap::new();
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
                    vakint.settings.number_of_terms_in_epsilon_expansion
                        - (integral.n_loops as i64)
                        - 1
                ),
            ],
            vakint.settings.clean_tmp_dir,
            vakint.settings.temporary_directory.clone(),
        )?;

        let mut evaluated_integral = vakint.process_form_output(form_result)?;
        debug!(
            "{}: raw result from FORM:\n{}",
            "AlphaLoop".green(),
            evaluated_integral
        );

        // Convert vectors back from pi(j) notation to p(i,j) notation
        for (s, t) in vector_mapping.iter() {
            evaluated_integral = s.into_pattern().replace_all(
                evaluated_integral.as_view(),
                &t.into_pattern().into(),
                None,
                None,
            );
        }

        if !vakint.settings.use_dot_product_notation {
            evaluated_integral =
                Vakint::convert_from_dot_notation(evaluated_integral.as_view(), false);
        }

        evaluated_integral = Pattern::parse("mUV").unwrap().replace_all(
            evaluated_integral.as_view(),
            &Atom::new_var(muv_sq_symbol)
                .pow((Atom::new_num(1) / Atom::new_num(2)).as_atom_view())
                .into_pattern()
                .into(),
            None,
            None,
        );

        let log_muv_mu_sq = fun!(
            State::LOG,
            Atom::new_var(muv_sq_symbol)
                / Atom::new_var(State::get_symbol(vakint.settings.mu_r_sq_symbol.as_str()))
        );

        let log_mu_sq = fun!(
            State::LOG,
            Atom::new_var(State::get_symbol(vakint.settings.mu_r_sq_symbol.as_str()))
        );

        //*((2*𝜋)^(2*{eps}))\
        //*((4*𝜋*exp(EulerGamma)))^{eps}
        //*((2*𝜋)^(4-2*{eps}))
        //*𝑖*(𝜋^((4-2*{eps})/2))
        // Make sure to split off the logarithmic terms with one term showing explicitely a ratio of scales and the other
        // having just a logarithm of the renormalization scale so that cancellations are symbolic when using `log_mu_sq`
        // in the normalization choice.
        // We must keep the name logmUVmu as it is reserved in the alphaloop implementation and corresponds to log(mUV^2/mu^2)
        // This is also the reason we do not simplify the expression exp(-logmUVmu+log_mu_sq)
        let alphaloop_normalization_correction = Atom::parse(
            format!(
                "(\
                    𝑖*(𝜋^((4-2*{eps})/2))\
                 * (exp(-EulerGamma))^({eps})\
                 * (exp(-logmUVmu-log_mu_sq))^({eps})\
                 )^{n_loops}",
                eps = self.settings.epsilon_symbol,
                n_loops = integral.n_loops
            )
            .as_str(),
        )
        .unwrap();

        evaluated_integral = evaluated_integral
            * alphaloop_normalization_correction
            * S.n_loops.into_pattern().replace_all(
                vakint
                    .settings
                    .get_integral_normalization_factor_atom()?
                    .as_view(),
                &Atom::new_num(integral.n_loops as i64).into_pattern().into(),
                None,
                None,
            );

        let expanded_evaluation = match evaluated_integral.series(
            State::get_symbol(vakint.settings.epsilon_symbol.as_str()),
            Atom::Zero.as_atom_view(),
            Rational::from(
                vakint.settings.number_of_terms_in_epsilon_expansion
                    - (integral.n_loops as i64)
                    - 1,
            ),
            true,
        ) {
            Ok(a) => a,
            Err(e) => return Err(VakintError::SymbolicaError(e.to_string())),
        };
        evaluated_integral = expanded_evaluation.to_atom();

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

        Ok(evaluated_integral)
    }

    fn get_pysecdec_version(&self) -> Result<String, VakintError> {
        let mut cmd = Command::new(self.settings.python_exe_path.as_str());
        cmd.arg("-c");
        cmd.arg("import pySecDec; print(pySecDec.__version__)");
        let output = if let Ok(o) = cmd.output() {
            o
        } else {
            return Err(VakintError::PySecDecUnavailable);
        };
        if !ExitStatus::success(&output.status) {
            return Err(VakintError::PySecDecUnavailable);
        }
        let output_str = String::from_utf8_lossy(&output.stdout).into_owned();
        let re = Regex::new(r"([\.|\d]+)").unwrap();
        let mut versions = vec![];
        for (_, [version]) in re.captures_iter(output_str.as_str()).map(|ci| ci.extract()) {
            versions.push(version);
        }
        if versions.is_empty() {
            return Err(VakintError::PySecDecVersion(format!(
                "Could not obtain PySecDec version from command:\n{:?}\nOutput was:\n{}",
                cmd, output_str
            )));
        }
        Ok(versions[0].into())
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
        for (_, [version]) in re.captures_iter(output_str.as_str()).map(|ci| ci.extract()) {
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

    pub fn convert_from_dot_notation(atom: AtomView, no_power: bool) -> Atom {
        let mut expr = atom.to_owned().expand();
        let mut running_dummy_index = 1;

        if no_power {
            while let Some(m) = Pattern::parse("dot(v1_(id1_),v2_(id2_))^n_")
                .unwrap()
                .pattern_match(
                    expr.as_view(),
                    &(Condition::from((State::get_symbol("v1_"), symbol_condition()))
                        & Condition::from((State::get_symbol("v2_"), symbol_condition()))
                        & Condition::from((State::get_symbol("id1_"), number_condition()))
                        & Condition::from((State::get_symbol("id2_"), number_condition()))
                        & Condition::from((State::get_symbol("n_"), number_condition()))),
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
                let pow =
                    get_integer_from_match(m.match_stack.get(State::get_symbol("n_")).unwrap())
                        .unwrap() as usize;
                let new_expression = (0..pow)
                    .map(|_i| {
                        let r = format!(
                            "{}({},{})*{}({},{})",
                            v1, id1, running_dummy_index, v2, id2, running_dummy_index
                        );
                        running_dummy_index += 1;
                        r
                    })
                    .collect::<Vec<_>>()
                    .join("*");

                expr = m.target.into_pattern().replace_all(
                    expr.as_view(),
                    &Pattern::parse(new_expression.as_str()).unwrap().into(),
                    None,
                    None,
                );
            }
        }

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
                .unwrap()
                .into(),
                None,
                None,
            );
            running_dummy_index += 1;
        }

        expr
    }

    pub fn convert_to_dot_notation(atom: AtomView) -> Atom {
        let mut old_expr = atom.to_owned().expand();

        loop {
            let mut expr = Pattern::parse("v_(id_,idx_)^n_").unwrap().replace_all(
                old_expr.as_view(),
                &Pattern::parse("dot(v_(id_),v_(id_))^(n_/2)")
                    .unwrap()
                    .into(),
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
                    &Pattern::parse("dot(v1_(id1_),v2_(id2_))").unwrap().into(),
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
                    &Pattern::parse("dot(v1_(id1_),v2_(id2_))").unwrap().into(),
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

            if expr == old_expr {
                break;
            } else {
                old_expr = expr.to_owned();
            }
        }
        old_expr
    }

    pub fn prepare_expression_for_form(&self, expression: Atom) -> Result<String, VakintError> {
        let processed = expression.clone();
        let processed = Pattern::parse(&self.settings.epsilon_symbol)
            .unwrap()
            .replace_all(
                processed.as_view(),
                &Pattern::parse("ep").unwrap().into(),
                None,
                None,
            );
        Ok(AtomPrinter::new_with_options(processed.as_view(), PrintOptions::file()).to_string())
    }

    pub fn process_form_output(&self, form_output: String) -> Result<Atom, VakintError> {
        let processed_form_str = form_output
            .replace("i_", "𝑖")
            .replace("\\\n", "\n")
            .split("\n")
            .map(|s| s.trim())
            .collect::<Vec<_>>()
            .join("");
        match Atom::parse(processed_form_str.as_str()) {
            Ok(mut processed) => {
                processed = Pattern::parse("rat(x_,y_)").unwrap().replace_all(
                    processed.as_view(),
                    &Pattern::parse("(x_/y_)").unwrap().into(),
                    None,
                    None,
                );
                processed = Pattern::parse("rat(x_)").unwrap().replace_all(
                    processed.as_view(),
                    &Pattern::parse("x_").unwrap().into(),
                    None,
                    None,
                );
                processed = Pattern::parse("ep").unwrap().replace_all(
                    processed.as_view(),
                    &Pattern::parse(&self.settings.epsilon_symbol)
                        .unwrap()
                        .into(),
                    None,
                    None,
                );
                processed = Pattern::parse("pi").unwrap().replace_all(
                    processed.as_view(),
                    &Pattern::parse("𝜋").unwrap().into(),
                    None,
                    None,
                );

                processed = Pattern::parse("g(idx1_,idx2_)").unwrap().replace_all(
                    processed.as_view(),
                    &Pattern::parse(format!("{}(idx1_,idx2_)", METRIC_SYMBOL).as_str())
                        .unwrap()
                        .into(),
                    Some(
                        &(Condition::from((State::get_symbol("idx1_"), number_condition()))
                            & Condition::from((State::get_symbol("idx2_"), number_condition()))),
                    ),
                    None,
                );
                processed = Pattern::parse("g(v1_,v2_)").unwrap().replace_all(
                    processed.as_view(),
                    &Pattern::parse("dot(v1_,v2_)").unwrap().into(),
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
                        &Pattern::parse("dot(v1_(args1_),v2_(args2_))")
                            .unwrap()
                            .into(),
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

    pub fn run_pysecdec(
        &self,
        input: &[(String, String)],
        options: Vec<String>,
        clean: bool,
        reused_path: Option<String>,
        temporary_directory: Option<String>,
    ) -> Result<Vec<String>, VakintError> {
        let mut generate_pysecdec_sources = true;
        let tmp_dir = &(if let Some(reused_path_specified) = reused_path.as_ref() {
            let specified_dir = PathBuf::from(reused_path_specified);
            if !specified_dir.exists() {
                if fs::create_dir(specified_dir.clone()).is_err() {
                    return Err(VakintError::PySecDecError(format!(
                        "Could not create user-specified directory to be reused for pySecDec run '{}'.",
                        reused_path_specified.green()
                    )));
                } else {
                    warn!("User-specified directory '{}' not found, and was instead created and pySecDec sources will be regenerated.", reused_path_specified);
                }
            } else {
                generate_pysecdec_sources = false;
                warn!("{}",format!("User requested to re-use existing directory '{}' for the pysecdec run.\nThis is of course potentially unsafe and should be used for debugging only.\nRemove that directory to start clean.", reused_path_specified).red());
            }
            specified_dir
        } else {
            let tmp_directory = if let Some(temp_path) = temporary_directory {
                &PathBuf::from(temp_path).join("vakint_temp_pysecdec")
            } else {
                &env::temp_dir().join("vakint_temp_pysecdec")
            };
            tmp_directory.clone()
        });

        if generate_pysecdec_sources {
            if tmp_dir.exists() {
                fs::remove_dir_all(tmp_dir)?;
            }
            fs::create_dir(tmp_dir)?;
            for input in input.iter() {
                fs::write(tmp_dir.join(&input.0).to_str().unwrap(), &input.1)?;
            }
        };

        let mut cmd = Command::new(self.settings.python_exe_path.as_str());
        cmd.arg(input[0].clone().0);
        for opt in options {
            cmd.arg(opt);
        }
        cmd.current_dir(tmp_dir);

        if !clean {
            info!("Running {} with command: {:?}", "PySecDec".green(), cmd);
            info!("You can follow the run live with 'tail -f follow_run.txt' in that temporary directory");
        } else {
            debug!("Running {} with command: {:?}", "PySecDec".green(), cmd);
            debug!("You can follow the run live with 'tail -f follow_run.txt' in that temporary directory");
        }

        let mut child = cmd.stderr(Stdio::piped()).stdout(Stdio::piped()).spawn()?;

        let stdout = child.stdout.take().unwrap();

        let reader = BufReader::new(stdout);
        let mut follow_file = File::create(tmp_dir.join("follow_run.txt"))?;

        for line in reader.lines() {
            let line_with_new_line = format!("{}\n", line?);
            follow_file.write_all(line_with_new_line.as_bytes())?;
            follow_file.flush()?;
        }

        let status = child.wait()?;

        if !ExitStatus::success(&status) {
            return Err(VakintError::FormError(
                "N/A".into(),
                "N/A".into(),
                format!("{:?}", cmd),
                tmp_dir.to_str().unwrap().into(),
            ));
        }
        if !tmp_dir.join("out.txt").exists() {
            return Err(VakintError::MissingFormOutput(
                "N/A".into(),
                format!("{:?}", cmd),
                tmp_dir.to_str().unwrap().into(),
            ));
        }
        let result = fs::read_to_string(tmp_dir.join("out.txt"))?;
        if clean {
            fs::remove_dir_all(tmp_dir)?;
        }
        Ok(result.split('\n').map(|s| s.into()).collect::<Vec<_>>())
    }

    pub fn run_form(
        &self,
        resources: &[String],
        input: (String, String),
        options: Vec<String>,
        clean: bool,
        temporary_directory: Option<String>,
    ) -> Result<String, VakintError> {
        let tmp_dir = if let Some(temp_path) = temporary_directory {
            &PathBuf::from(temp_path).join("vakint_temp")
        } else {
            &env::temp_dir().join("vakint_temp")
        };

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
        cmd.current_dir(tmp_dir);
        if !clean {
            info!("Running {} with command: {:?}", "FORM".green(), cmd);
        } else {
            debug!("Running {} with command: {:?}", "FORM".green(), cmd);
        }
        let output = cmd.stderr(Stdio::piped()).stdout(Stdio::piped()).output()?;
        if !ExitStatus::success(&output.status) {
            return Err(VakintError::FormError(
                String::from_utf8_lossy(&output.stderr).into(),
                String::from_utf8_lossy(&output.stdout).into(),
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
        if clean {
            fs::remove_dir_all(tmp_dir)?;
        }
        Ok(result)
    }

    pub fn evaluate(&self, input: AtomView) -> Result<Atom, VakintError> {
        let mut vakint_expr = VakintExpression::try_from(input)?;
        vakint_expr.canonicalize(self, &self.topologies, false)?;
        vakint_expr.tensor_reduce(self)?;
        vakint_expr.evaluate_integral(self)?;
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
