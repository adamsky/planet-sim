mod config;
mod constants;
mod sim;

pub type Float = f32;
pub type Int = i32;
pub type Vec2d<T> = Vec<Vec<T>>;
pub type Vec3d<T> = Vec<Vec<Vec<T>>>;
pub type Vec4d<T> = Vec<Vec<Vec<Vec<T>>>>;

pub use sim::Planet;
pub use sim::Sim;
