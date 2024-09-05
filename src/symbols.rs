use std::sync::LazyLock;

use symbolica::{
    atom::{Atom, Symbol},
    state::{FunctionAttribute, State},
};

#[allow(dead_code)]
pub struct VakintSymbols {
    pub uedge: Symbol,
    pub dot: Symbol,
    pub g: Symbol,
    pub x: Symbol,
    pub y: Symbol,
    pub xa: Atom,
    pub ya: Atom,
    pub one: Atom,
    pub zero: Atom,
    pub error_flag_symbol: Symbol,
    pub error_flag: Atom,
    pub n_loops: Atom,
}

#[allow(dead_code)]
pub static S: LazyLock<VakintSymbols> = LazyLock::new(|| VakintSymbols {
    uedge: State::get_symbol_with_attributes("uedge", &[FunctionAttribute::Symmetric]).unwrap(),
    dot: State::get_symbol_with_attributes(
        "dot",
        &[FunctionAttribute::Symmetric, FunctionAttribute::Linear],
    )
    .unwrap(),
    g: State::get_symbol_with_attributes("g", &[FunctionAttribute::Symmetric]).unwrap(),
    x: State::get_symbol("x"),
    y: State::get_symbol("y"),
    xa: Atom::new_var(State::get_symbol("xa")),
    ya: Atom::new_var(State::get_symbol("ya")),
    one: Atom::new_num(1),
    zero: Atom::Zero,
    error_flag_symbol: State::get_symbol("ERROR"),
    error_flag: Atom::new_var(State::get_symbol("ERROR")),
    n_loops: Atom::new_var(State::get_symbol("n_loops")),
});
