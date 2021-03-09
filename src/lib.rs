mod config;
mod constants;
mod sim;

pub type FloatNum = f32;
pub type Int = i32;
pub type Vec2d<T> = Vec<Vec<T>>;
pub type Vec3d<T> = Vec<Vec<Vec<T>>>;
pub type Vec4d<T> = Vec<Vec<Vec<Vec<T>>>>;

pub use sim::Planet;
pub use sim::Sim;

pub fn min(first: FloatNum, second: FloatNum) -> FloatNum {
    first.min(second)
}

pub fn max(first: FloatNum, second: FloatNum) -> FloatNum {
    first.max(second)
}
