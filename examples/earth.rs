use planet_sim::{Planet, Sim};

fn main() {
    let mut planet = Sim::new_planet(Planet::Earth);
    std::thread::sleep(std::time::Duration::from_secs(2));
    // planet.step();
}
