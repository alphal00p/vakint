use std::sync::LazyLock;

use symbolica::{
    atom::{Atom, FunctionAttribute, Symbol},
    symb,
};

pub static METRIC_SYMBOL: &str = "g";
pub static LOOP_MOMENTUM_SYMBOL: &str = "k";
pub static EXTERNAL_MOMENTUM_SYMBOL: &str = "p";

#[allow(dead_code)]
pub struct VakintSymbols {
    pub uedge: Symbol,
    pub dot: Symbol,
    pub vkdot: Symbol,
    pub g: Symbol,
    pub x: Symbol,
    pub y: Symbol,
    pub xa: Atom,
    pub ya: Atom,
    pub x_: Symbol,
    pub y_: Symbol,
    pub x_a: Atom,
    pub y_a: Atom,
    pub one: Atom,
    pub zero: Atom,
    pub error_flag_symbol: Symbol,
    pub error_flag: Atom,
    pub n_loops: Atom,
    pub p: Symbol,
    pub k: Symbol,
    pub id1_: Symbol,
    pub id2_: Symbol,
    pub id1_a: Atom,
    pub id2_a: Atom,
    pub pow: Symbol,
    pub pow_: Symbol,
    pub fun_: Symbol,
    pub fun_a: Atom,
    pub any___: Symbol,
    pub any_a___: Atom,
    pub n_: Symbol,
    pub a_: Symbol,
    pub b_: Symbol,
    pub lambda: Symbol,
    pub lambda_a: Atom,
}

pub static S: LazyLock<VakintSymbols> = LazyLock::new(|| VakintSymbols {
    uedge: Symbol::new_with_attributes("uedge", &[FunctionAttribute::Symmetric]).unwrap(),
    dot: Symbol::new_with_attributes(
        "dot",
        &[FunctionAttribute::Symmetric, FunctionAttribute::Linear],
    )
    .unwrap(),
    vkdot: Symbol::new_with_attributes(
        "vkdot",
        &[FunctionAttribute::Symmetric, FunctionAttribute::Linear],
    )
    .unwrap(),
    g: Symbol::new_with_attributes("METRIC_SYMBOL", &[FunctionAttribute::Symmetric]).unwrap(),
    x: symb!("x"),
    y: symb!("y"),
    xa: Atom::new_var(symb!("xa")),
    ya: Atom::new_var(symb!("ya")),
    x_: symb!("x_"),
    y_: symb!("y_"),
    x_a: Atom::new_var(symb!("x_")),
    y_a: Atom::new_var(symb!("y_")),
    one: Atom::new_num(1),
    zero: Atom::Zero,
    error_flag_symbol: symb!("ERROR"),
    error_flag: Atom::new_var(symb!("ERROR")),
    n_loops: Atom::new_var(symb!("n_loops")),
    p: symb!(EXTERNAL_MOMENTUM_SYMBOL),
    k: symb!(LOOP_MOMENTUM_SYMBOL),
    id1_: symb!("id1_"),
    id2_: symb!("id2_"),
    id1_a: Atom::new_var(symb!("id1_")),
    id2_a: Atom::new_var(symb!("id2_")),
    pow: symb!("pow"),
    pow_: symb!("pow_"),
    fun_: symb!("fun_"),
    fun_a: Atom::new_var(symb!("fun_")),
    any___: symb!("any___"),
    any_a___: Atom::new_var(symb!("any___")),
    n_: symb!("n_"),
    a_: symb!("a_"),
    b_: symb!("b_"),
    lambda: symb!("VakintLambdaScalingAnalysis"),
    lambda_a: Atom::new_var(symb!("VakintLambdaScalingAnalysis")),
});
