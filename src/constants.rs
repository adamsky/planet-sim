//! Global constants.

use crate::{Float, Int};

/// Number of tracers, 1st. reserved for q
pub const NTRACE: usize = 1;
/// Number of processes
pub const NPRO: usize = 4;

/// Number of latitudes
pub const NLAT: usize = 64;
/// Number of levels
pub const NLEV: usize = 10;

/// Number of longitudes
pub const NLON: usize = NLAT + NLAT;
/// Triangular truncation
pub const NTRU: usize = (NLON - 1) / 3;
/// Latitudes per process
pub const NLPP: usize = NLAT / NPRO;
/// Horizontal part
pub const NHOR: usize = NLON * NLPP;
/// Number of gridpoints
pub const NUGP: usize = NLON * NLAT;
/// Dimension of packed fields
pub const NPGP: usize = NLON * NLAT / 2;
/// Levels - 1
pub const NLEM: usize = NLEV - 1;
/// Levels + 1
pub const NLEP: usize = NLEV + 1;
/// Levels squared
pub const NLSQ: usize = NLEV * NLEV;
/// Truncation + 1
pub const NTP1: usize = NTRU + 1;
/// No of real global modes
pub const NRSP: usize = (NTRU + 1) * (NTRU + 2);
/// No of complex global modes
pub const NCSP: usize = NRSP / 2;
/// Modes per process
pub const NSPP: usize = (NRSP + NPRO - 1) / NPRO;
/// Dim of spectral fields
pub const NESP: usize = NSPP * NPRO;
/// Dim of Vert. Coord. Tab
pub const NVCT: usize = 2 * (NLEV + 1);
/// Dim for zonal mean diagnostics
pub const NZOM: usize = 2 * NTP1;
/// Master node
pub const NROOT: usize = 0;

/// ez = 1 / sqrt(3/8)
pub const EZ: Float = 1.63299310207;
/// Pi
pub const PI: Float = 3.14159265359;
/// 2 Pi
pub const TWOPI: Float = PI + PI;
/// Gas constant for water vapour
pub const RV: Float = 461.51;
/// Specific heat for water vapour
pub const ACPV: Float = 1870.;
/// Melting point (CO2) - for Mars
pub const TMELT_CO2: Float = 148.0;

/// Stefan-Bolzman constant
pub const SBK: Float = 5.67E-8;
