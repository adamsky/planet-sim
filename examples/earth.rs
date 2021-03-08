use plasim_rs::{Planet, Sim};

fn main() {
    let mut planet = Sim::new_planet(Planet::Earth);
    planet.step();
}
