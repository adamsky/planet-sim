use planet_sim::{Planet, Sim};

fn main() {
    let mut planet = Sim::new_planet(Planet::Earth);
    planet.step();
    // std::thread::sleep(std::time::Duration::from_secs(2));
}
