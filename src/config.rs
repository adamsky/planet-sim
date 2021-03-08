use serde::{Deserialize, Serialize};

use crate::{Float, Int};

#[derive(Deserialize, Serialize)]
pub struct Config {
    planet: PlanetConfig,
}

#[derive(Deserialize, Serialize)]
pub struct PlanetConfig {
    // astronomy
    /// 1 = fix the planetary orbit
    nfixorb: Int,
    /// Eccentricity for fixed orbits (AMIP-II value)
    eccen: Float,
    /// Longitude of perihelion
    mvelp: Float,
    /// Obliquity (deg) (AMIP-II)
    obliq: Float,
    /// Rotation speed (factor)
    rotspd: Float,
    /// Sidereal day length
    sidereal_day: Float,
    /// Solar day length
    solar_day: Float,
    /// Sidereal year length
    sidereal_year: Float,
    /// Tropical year length
    tropical_year: Float,

    // atmosphere
    /// Kappa (Poisson constant R/Cp)
    akap: Float,
    /// Lapse rate
    alr: Float,
    /// Gas constant
    gascon: Float,
    // Parameters for Magnus-Teten-Formula
    // for saturation vapor pressure over liquid water
    /// Parameter for Magnus-Teten-Formula
    ra1: Float,
    /// Parameter for Magnus-Teten-Formula
    ra2: Float,
    /// Parameter for Magnus-Teten-Formula
    ra4: Float,

    // numerics
    /// Time filter constant
    pnu: Float,

    // physics
    /// Gravity (mean on NN)
    ga: Float,
    /// Radius
    plarad: Float,

    // radiation
    /// Solar constant
    gsol0: Float,

    // radmod
    /// Mass (10^24 kg)
    p_mass: Float,
    /// Volume (10^10 km3)
    p_volume: Float,
    /// Equatorial radius
    p_radius_eq: Float,
    /// Polar radius
    p_radius_po: Float,
    /// Mean radius
    p_radius_me: Float,
    /// Ellipticity
    p_ellipticity: Float,
    /// Density (kg/m3)
    p_density: Float,
    /// Bond albedo
    p_albedo: Float,
    /// Black body temperature
    p_blackt: Float,
    /// Sidereal rotation period
    p_sidrot: Float,
    /// Perihelion (10^6 km)
    p_perihelion: Float,
    /// Aphelion (10^6 km)
    p_aphelion: Float,
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
