use colored::Colorize;
use log::debug;
use symbolica::atom::{Atom, AtomView};

use crate::{ReplacementRules, Vakint, VakintError};

#[allow(unused)]
impl Vakint {
    pub fn matad_evaluate(
        &self,
        vakint: &Vakint,
        input_numerator: AtomView,
        integral_specs: &ReplacementRules,
    ) -> Result<Atom, VakintError> {
        let integral = integral_specs.canonical_topology.get_integral();

        debug!(
            "Processing the following integral with {}:\n{}",
            "MATAD".green(),
            integral
        );

        println!("integral_specs=\n{}", integral_specs);
        println!("Numerator=\n{}", input_numerator);

        let mut numerator = Vakint::convert_to_dot_notation(input_numerator);

        Ok(Atom::Zero)
    }
}
