use std::sync::LazyLock;

use symbolica::{
    atom::{Atom, Symbol},
    state::{FunctionAttribute, State},
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
}

pub static S: LazyLock<VakintSymbols> = LazyLock::new(|| VakintSymbols {
    uedge: State::get_symbol_with_attributes("uedge", &[FunctionAttribute::Symmetric]).unwrap(),
    dot: State::get_symbol_with_attributes(
        "dot",
        &[FunctionAttribute::Symmetric, FunctionAttribute::Linear],
    )
    .unwrap(),
    vkdot: State::get_symbol_with_attributes(
        "vkdot",
        &[FunctionAttribute::Symmetric, FunctionAttribute::Linear],
    )
    .unwrap(),
    g: State::get_symbol_with_attributes("METRIC_SYMBOL", &[FunctionAttribute::Symmetric]).unwrap(),
    x: State::get_symbol("x"),
    y: State::get_symbol("y"),
    xa: Atom::new_var(State::get_symbol("xa")),
    ya: Atom::new_var(State::get_symbol("ya")),
    x_: State::get_symbol("x_"),
    y_: State::get_symbol("y_"),
    x_a: Atom::new_var(State::get_symbol("x_")),
    y_a: Atom::new_var(State::get_symbol("y_")),
    one: Atom::new_num(1),
    zero: Atom::Zero,
    error_flag_symbol: State::get_symbol("ERROR"),
    error_flag: Atom::new_var(State::get_symbol("ERROR")),
    n_loops: Atom::new_var(State::get_symbol("n_loops")),
    p: State::get_symbol(EXTERNAL_MOMENTUM_SYMBOL),
    k: State::get_symbol(LOOP_MOMENTUM_SYMBOL),
    id1_: State::get_symbol("id1_"),
    id2_: State::get_symbol("id2_"),
    id1_a: Atom::new_var(State::get_symbol("id1_")),
    id2_a: Atom::new_var(State::get_symbol("id2_")),
    pow: State::get_symbol("pow"),
    pow_: State::get_symbol("pow_"),
    fun_: State::get_symbol("fun_"),
    fun_a: Atom::new_var(State::get_symbol("fun_")),
    any___: State::get_symbol("any___"),
    any_a___: Atom::new_var(State::get_symbol("any___")),
    n_: State::get_symbol("n_"),
});
