use serde::{Deserialize, Serialize};

use crate::{FloatNum, Int};

#[derive(Deserialize, Serialize)]
pub struct Config {
    planet: PlanetConfig,
}

#[derive(Deserialize, Serialize)]
pub struct PlanetConfig {
    // astronomy
    /// 1 = fix the planetary orbit
    pub nfixorb: Int,
    /// Eccentricity for fixed orbits (AMIP-II value)
    pub eccen: FloatNum,
    /// Longitude of perihelion
    pub mvelp: FloatNum,
    /// Obliquity (deg) (AMIP-II)
    pub obliq: FloatNum,
    /// Rotation speed (factor)
    pub rotspd: FloatNum,
    /// Sidereal day length
    pub sidereal_day: FloatNum,
    /// Solar day length
    pub solar_day: FloatNum,
    /// Sidereal year length
    pub sidereal_year: FloatNum,
    /// Tropical year length
    pub tropical_year: FloatNum,

    // atmosphere
    /// Kappa (Poisson constant R/Cp)
    pub akap: FloatNum,
    /// Lapse rate
    pub alr: FloatNum,
    /// Gas constant
    pub gascon: FloatNum,
    // Parameters for Magnus-Teten-Formula
    // for saturation vapor pressure over liquid water
    /// Parameter for Magnus-Teten-Formula
    pub ra1: FloatNum,
    /// Parameter for Magnus-Teten-Formula
    pub ra2: FloatNum,
    /// Parameter for Magnus-Teten-Formula
    pub ra4: FloatNum,

    // numerics
    /// Time filter constant
    pub pnu: FloatNum,

    // physics
    /// Gravity (mean on NN)
    pub ga: FloatNum,
    /// Radius
    pub plarad: FloatNum,

    // radiation
    /// Solar constant
    pub gsol0: FloatNum,

    // radmod
    /// Mass (10^24 kg)
    pub p_mass: FloatNum,
    /// Volume (10^10 km3)
    pub p_volume: FloatNum,
    /// Equatorial radius
    pub p_radius_eq: FloatNum,
    /// Polar radius
    pub p_radius_po: FloatNum,
    /// Mean radius
    pub p_radius_me: FloatNum,
    /// Ellipticity
    pub p_ellipticity: FloatNum,
    /// Density (kg/m3)
    pub p_density: FloatNum,
    /// Bond albedo
    pub p_albedo: FloatNum,
    /// Black body temperature
    pub p_blackt: FloatNum,
    /// Sidereal rotation period
    pub p_sidrot: FloatNum,
    /// Perihelion (10^6 km)
    pub p_perihelion: FloatNum,
    /// Aphelion (10^6 km)
    pub p_aphelion: FloatNum,
    // Sidereal orbit period
    // p_sidorbit: Float, // = sidereal_year / sidereal_day
}

impl Default for PlanetConfig {
    fn default() -> Self {
        Self::earth()
    }
}

impl PlanetConfig {
    pub fn earth() -> Self {
        Self {
            nfixorb: 0,
            eccen: 0.016715,
            mvelp: 102.7,
            obliq: 23.441,
            rotspd: 1.0,
            // 23h 56m 04s
            sidereal_day: 86164.0916,
            // 24h 00m 00s
            solar_day: 86400.0,
            // 365d 06h 09m 09s
            sidereal_year: 31558149.0,
            // 365d 05h 49m 16s
            tropical_year: 31556956.0,
            akap: 0.286,
            alr: 0.0065,
            gascon: 287.0,
            ra1: 610.78,
            ra2: 17.2693882,
            ra4: 35.86,
            pnu: 0.1,
            ga: 9.80665,
            plarad: 6371220.0,
            gsol0: 1365.0,

            p_mass: 5.9736,
            p_volume: 108.321,
            p_radius_eq: 6378.0,
            p_radius_po: 6356.0,
            p_radius_me: 6371220.0 / 1000.0, // plarad / 1000.0,
            p_ellipticity: 0.0034,
            p_density: 5520.0,
            p_albedo: 0.385,
            p_blackt: 247.3,
            p_sidrot: 86164.0916 / 3600.0, // sidereal_day / 3600.0,
            p_perihelion: 147.1,
            p_aphelion: 152.1,
        }
    }

    pub fn mars() -> Self {
        Self {
            nfixorb: 0,
            eccen: 0.016715,
            mvelp: 102.7,
            obliq: 23.441,
            rotspd: 1.0,
            // 23h 56m 04s
            sidereal_day: 86164.0916,
            // 24h 00m 00s
            solar_day: 86400.0,
            // 365d 06h 09m 09s
            sidereal_year: 31558149.0,
            // 365d 05h 49m 16s
            tropical_year: 31556956.0,
            akap: 0.286,
            alr: 0.0065,
            gascon: 287.0,
            ra1: 610.78,
            ra2: 17.2693882,
            ra4: 35.86,
            pnu: 0.1,
            ga: 9.80665,
            plarad: 6371220.0,
            gsol0: 1365.0,

            p_mass: 0.6419,
            p_volume: 16.318,
            p_radius_eq: 3393.0,
            p_radius_po: 3373.0,
            p_radius_me: 3390.0,
            p_ellipticity: 0.0065,
            p_density: 3933.0,
            p_albedo: 0.16,
            p_blackt: 216.6,
            p_sidrot: 24.6229,
            p_perihelion: 206.6,
            p_aphelion: 249.2,
        }
    }
}
