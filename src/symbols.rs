use std::sync::LazyLock;

use crate::utils::vakint_macros::vk_symbol;
use symbolica::{
    atom::{Atom, Symbol},
    symbol,
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
    pub g_form: Symbol,
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
    uedge: symbol!(format!("{}::uedge",crate::NAMESPACE); Symmetric).unwrap(),
    dot: symbol!(
        format!("{}::dot",crate::NAMESPACE); Symmetric,Linear
    )
    .unwrap(),
    vkdot: symbol!(
        format!("{}::vkdot",crate::NAMESPACE);  Symmetric, Linear
    )
    .unwrap(),
    g: symbol!(format!("{}::{}",crate::NAMESPACE,METRIC_SYMBOL); Symmetric).unwrap(),
    g_form: symbol!(format!("{}::g",crate::NAMESPACE); Symmetric).unwrap(),
    x: vk_symbol!("x"),
    y: vk_symbol!("y"),
    xa: Atom::new_var(vk_symbol!("xa")),
    ya: Atom::new_var(vk_symbol!("ya")),
    x_: vk_symbol!("x_"),
    y_: vk_symbol!("y_"),
    x_a: Atom::new_var(vk_symbol!("x_")),
    y_a: Atom::new_var(vk_symbol!("y_")),
    one: Atom::new_num(1),
    zero: Atom::Zero,
    error_flag_symbol: vk_symbol!("ERROR"),
    error_flag: Atom::new_var(vk_symbol!("ERROR")),
    n_loops: Atom::new_var(vk_symbol!("n_loops")),
    p: vk_symbol!(EXTERNAL_MOMENTUM_SYMBOL),
    k: vk_symbol!(LOOP_MOMENTUM_SYMBOL),
    id1_: vk_symbol!("id1_"),
    id2_: vk_symbol!("id2_"),
    id1_a: Atom::new_var(vk_symbol!("id1_")),
    id2_a: Atom::new_var(vk_symbol!("id2_")),
    pow: vk_symbol!("pow"),
    pow_: vk_symbol!("pow_"),
    fun_: vk_symbol!("fun_"),
    fun_a: Atom::new_var(vk_symbol!("fun_")),
    any___: vk_symbol!("any___"),
    any_a___: Atom::new_var(vk_symbol!("any___")),
    n_: vk_symbol!("n_"),
    a_: vk_symbol!("a_"),
    b_: vk_symbol!("b_"),
    lambda: vk_symbol!("VakintLambdaScalingAnalysis"),
    lambda_a: Atom::new_var(vk_symbol!("VakintLambdaScalingAnalysis")),
});
