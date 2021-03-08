pub mod radiation;

use crate::config::PlanetConfig;
use crate::constants::{NESP, NHOR, NLAT, NLEP, NLEV, NLON, NLPP, NPRO, NSPP, NTRACE, NTRU};
use crate::{Float, Int};

#[derive(Default)]
pub struct Sim {
    // resolution: Resolution,
    datetime: DateTime,

    int_scalars: IntScalars,

    spec_arrays: SpectralArrays,
    grid_arrays_undim: GridArraysUndim,
    grid_arrays_dim: GridArraysDim,
    basic_radiation: BasicRadiation,
    rad: radiation::Radiation,
    surface: Surface,
    water: Water,
    biome: Biome,
    soil: Soil,

    diagno_arrays: DiagnosticArrays,
    lat_arrays: LatitudeArrays,
    lev_arrays: LevelArrays,

    seed: RandomSeed,
    planet_vars: PlanetVars,
}

impl Sim {
    pub fn new(planet_cfg: PlanetConfig) -> Self {
        let sim: Sim = Default::default();
        sim
    }

    pub fn new_planet(planet: Planet) -> Self {
        let cfg = match planet {
            Planet::Earth => PlanetConfig::earth(),
            Planet::Mars => PlanetConfig::mars(),
            // _ => unimplemented!(),
        };
        Sim::new(cfg)
    }

    pub fn step(&mut self) {
        unimplemented!()
    }
}

pub enum Planet {
    Earth,
    Mars,
}

// struct Resolution {
//     /// Number of latitudes
//     n_latitudes: u8,
//     /// Number of levels
//     n_levels: u8,
//     // /// Number of processes
//     // NPRO_ATM: u8
// }
//
// impl Default for Resolution {
//     fn default() -> Self {
//         Self {
//             n_latitudes: 64,
//             n_levels: 10,
//         }
//     }
// }

struct IntScalars {
    /// Add noise for kick > 0
    pub kick: Int,
    /// Global switch for planet mars
    pub mars: Int,
    /// master switch for output: 0=no output
    pub noutput: Int,
    /// write data interval: 0 = once per day
    pub nafter: Int,
    /// 1: switch to aqua planet mode
    pub naqua: Int,
    /// 1: run vegetation module
    pub nveg: Int,
    /// number of modes to print
    pub ncoeff: Int,
    /// write diagnostics interval 0 = every 10th. day
    pub ndiag: Int,
    /// 1: run with GUI
    pub ngui: Int,

    /// index of longitude for column mode
    pub sellon: Int,
    /// number of initial timesteps
    pub nkits: Int,
    /// 1 for true, 0 for false
    pub nrestart: Int,
    /// switches radiation off/on  1/0
    pub nrad: Int,
    /// vertical diffusion 1/0
    pub nflux: Int,
    /// advection 1/0=(y/n)
    pub nadv: Int,
    /// horizontal diffusion 1/0=(y/n)
    pub nhordif: Int,
    /// equidistant sigma levels (1/0)=(y/n)
    pub neqsig: Int,
    /// comprehensive print out (only for checks!)
    pub nprint: Int,
    /// grid point for print out (only for checks!)
    pub nprhor: Int,
    /// pack spectral fields on output
    pub npacksp: Int,
    /// pack gridpoint fields on output
    pub npackgp: Int,
    /// accumulation counter for diagnistics
    pub naccuout: Int,
    /// switch for frank's gp-diagnostic arrays
    pub ndiaggp: Int,
    /// switch for frank's sp-diagnostic arrays
    pub ndiagsp: Int,
    /// switch for cloud forcing diagnostic
    pub ndiagcf: Int,

    /// number of additional 2-d gp-diagnostic arrays
    pub ndiaggp2d: Int,
    /// number of additional 3-d gp-diagnostic arrays
    pub ndiaggp3d: Int,
    /// number of additional 2-d sp-diagnostic arrays
    pub ndiagsp2d: Int,
    /// number of additional 3-d sp-diagnostic arrays
    pub ndiagsp3d: Int,
    /// divergence damping countdown
    pub ndivdamp: Int,
    /// critical wavenumber for horizontal diffusion
    pub nhdiff: Int,
    /// switch for time use diagnostics
    pub ntime: Int,
    /// radiation day for perpetual integration
    pub nperpetual: Int,
    /// number of sea points on grid
    pub n_sea_points: Int,
    /// switch for entropy diagnostics
    pub nentropy: Int,
    /// switch for 3d entropy diagnostics
    pub nentro3d: Int,
    /// switch for energy diagnostics
    pub nenergy: Int,
    /// switch for 3d energy diagnostics
    pub nener3d: Int,
    /// switch for heating due to momentum dissipation
    pub ndheat: Int,
    /// length of random seed (set by lib call)
    pub nseedlen: Int,
    /// enable (1) or disable (0) Semi Lagrangian Advection
    pub nsela: Int,
    /// switch for LnPs initialization
    pub nspinit: Int,
    /// switch for top sponge layer
    pub nsponge: Int,
    /// 1: spectral q   0: gridpoint q (semi-Langrangian)
    pub nqspec: Int,
}

impl Default for IntScalars {
    fn default() -> Self {
        Self {
            kick: 1,
            mars: 0,
            noutput: 1,
            nafter: 0,
            naqua: 0,
            nveg: 1,
            ncoeff: 0,
            ndiag: 0,
            ngui: 0,

            sellon: 1,
            nkits: 3,
            nrestart: 0,
            nrad: 1,
            nflux: 1,
            nadv: 1,
            nhordif: 1,
            neqsig: 0,
            nprint: 0,
            nprhor: 0,
            npacksp: 0,
            npackgp: 0,
            naccuout: 0,
            ndiaggp: 0,
            ndiagsp: 0,
            ndiagcf: 0,

            ndiaggp2d: 0,
            ndiaggp3d: 0,
            ndiagsp2d: 0,
            ndiagsp3d: 0,
            ndivdamp: 0,
            nhdiff: 15,
            ntime: 0,
            nperpetual: 0,
            n_sea_points: 0,
            nentropy: 0,
            nentro3d: 0,
            nenergy: 0,
            nener3d: 0,
            ndheat: 1,
            nseedlen: 0,
            nsela: 0,
            nspinit: 0,
            nsponge: 0,
            nqspec: 1,
            // ..Default::default()
        }
    }
}

struct SpectralArrays {
    /// Spectral Divergence
    pub sd: [[Float; NESP]; NLEV],
    /// Spectral Temperature
    pub st: [[Float; NESP]; NLEV],
    /// Spectral Vorticity
    pub sz: [[Float; NESP]; NLEV],
    /// Spectral Specific Humidity
    pub sq: [[Float; NESP]; NLEV],
    /// Spectral Pressure (ln Ps)
    pub sp: [Float; NESP],
    /// Spectral Orography
    pub so: [Float; NESP],
    /// Spectral Restoration Temperature
    pub sr: [[Float; NESP]; NLEV],

    /// Spectral Divergence  Partial
    pub sdp: [[Float; NSPP]; NLEV],
    /// Spectral Temperature Partial
    pub stp: [[Float; NSPP]; NLEV],
    /// Spectral Vorticity   Partial
    pub szp: [[Float; NSPP]; NLEV],
    /// Spectral S.Humidity  Partial
    pub sqp: [[Float; NSPP]; NLEV],
    /// Spectral Pressure    Partial
    pub spp: [Float; NSPP],
    /// Spectral Orography   Partial
    pub sop: [Float; NSPP],
    /// Spectral Restoration Partial
    pub srp: [[Float; NSPP]; NLEV],

    /// Spectral Divergence  Tendency
    pub sdt: [[Float; NSPP]; NLEV],
    /// Spectral Temperature Tendency
    pub stt: [[Float; NSPP]; NLEV],
    /// Spectral Vorticity   Tendency
    pub szt: [[Float; NSPP]; NLEV],
    /// Spectral S.Humidity  Tendency
    pub sqt: [[Float; NSPP]; NLEV],
    /// Spectral Pressure    Tendency
    pub spt: [Float; NSPP],

    /// Spectral Divergence  Minus
    pub sdm: [[Float; NSPP]; NLEV],
    /// Spectral Temperature Minus
    pub stm: [[Float; NSPP]; NLEV],
    /// Spectral Vorticity   Minus
    pub szm: [[Float; NSPP]; NLEV],
    /// Spectral S.Humidity  Minus
    pub sqm: [[Float; NSPP]; NLEV],
    /// Spectral Pressure    Minus
    pub spm: [Float; NSPP],

    /// horizontal diffusion
    pub sak: [[Float; NESP]; NLEV],
    /// horizontal diffusion partial
    pub sakpp: [[Float; NSPP]; NLEV],
    /// specific humidity for output
    pub sqout: [[Float; NESP]; NLEV],
    /// Factors for output normalization
    pub spnorm: [Float; NESP],

    /// Holds wavenumber
    pub nindex: [Int; NESP],
    /// Holds wavenumber
    pub nscatsp: [Int; NPRO],
    /// Holds wavenumber
    pub ndel: [Int; NLEV],

    /// Difference between instances
    pub sdd: Vec<Vec<Float>>,
    /// Difference between instances
    pub std: Vec<Vec<Float>>,
    /// Difference between instances
    pub szd: Vec<Vec<Float>>,
    /// Difference between instances
    pub spd: Vec<Float>,
}

impl Default for SpectralArrays {
    fn default() -> Self {
        Self {
            nindex: [NTRU as Int; NESP],
            nscatsp: [NSPP as Int; NPRO],
            ndel: [2; NLEV],
            ..Default::default()
        }
    }
}

/// Global Gridpoint Arrays (un-dimensionalized)
struct GridArraysUndim {
    /// Divergence
    pub gd: [[Float; NHOR]; NLEV],
    /// Temperature (-t0)
    pub gt: [[Float; NHOR]; NLEV],
    /// Absolute vorticity
    pub gz: [[Float; NHOR]; NLEV],
    /// Spec. humidity
    pub gq: [[Float; NHOR]; NLEV],
    /// Zonal wind (*cos(phi))
    pub gu: [[Float; NHOR]; NLEV],
    /// Meridional wind (*cos(phi))
    pub gv: [[Float; NHOR]; NLEV],
    /// t-tendency
    pub gtdt: [[Float; NHOR]; NLEV],
    /// q-tendency
    pub gqdt: [[Float; NHOR]; NLEV],
    /// u-tendency
    pub gudt: [[Float; NHOR]; NLEV],
    /// v-tendency
    pub gvdt: [[Float; NHOR]; NLEV],
    /// Surface pressure or ln(ps)
    pub gp: [Float; NHOR],
    /// dln(ps)/dphi
    pub gpj: [Float; NHOR],

    /// 1/cos(phi)**2
    pub rcsq: [Float; NHOR],
}

impl Default for GridArraysUndim {
    fn default() -> Self {
        Self {
            gd: [[0.; NHOR]; NLEV],
            gt: [[0.; NHOR]; NLEV],
            gz: [[0.; NHOR]; NLEV],
            gq: [[0.; NHOR]; NLEV],
            gu: [[0.; NHOR]; NLEV],
            gv: [[0.; NHOR]; NLEV],
            gtdt: [[0.; NHOR]; NLEV],
            gqdt: [[0.; NHOR]; NLEV],
            gudt: [[0.; NHOR]; NLEV],
            gvdt: [[0.; NHOR]; NLEV],
            gp: [0.; NHOR],
            gpj: [0.; NHOR],
            rcsq: [0.; NHOR],
        }
    }
}

/// Global Gridpoint Arrays (dimensionalized)
struct GridArraysDim {
    /// Temperature
    pub dt: [[Float; NHOR]; NLEP],
    /// Spec. humidity
    pub dq: [[Float; NHOR]; NLEP],
    /// Zonal wind [m/s]
    pub du: [[Float; NHOR]; NLEP],
    /// Meridional wind [m/s]
    pub dv: [[Float; NHOR]; NLEP],
    /// Surface pressure
    pub dp: [Float; NHOR],

    /// Saturation humidity
    pub dqsat: [[Float; NHOR]; NLEP],
    /// Adiabatic q-tendencies (for eg kuo)
    pub dqt: [[Float; NHOR]; NLEP],
    /// Cloud cover
    pub dcc: [[Float; NHOR]; NLEP],
    /// Liquid water content
    pub dql: [[Float; NHOR]; NLEP],
    /// Vertical velocity (dp/dt)
    pub dw: [[Float; NHOR]; NLEV],
    /// t-tendency
    pub dtdt: [[Float; NHOR]; NLEP],
    /// q-tendency
    pub dqdt: [[Float; NHOR]; NLEP],
    /// u-tendency
    pub dudt: [[Float; NHOR]; NLEP],
    /// v-tendency
    pub dvdt: [[Float; NHOR]; NLEP],

    /// Surface pressure at time t
    pub dp0: [Float; NHOR],
    /// Zonal wind at time t
    pub du0: [[Float; NHOR]; NLEP],
    /// Meridional wind at time t
    pub dv0: [[Float; NHOR]; NLEP],
    /// Trace array
    pub dtrace: [[[[Float; NLON]; NLAT]; NLEV]; NTRACE],
}

impl Default for GridArraysDim {
    fn default() -> Self {
        Self {
            dtrace: [[[[1.; NLON]; NLAT]; NLEV]; NTRACE],
            ..Default::default()
        }
    }
}

struct BasicRadiation {
    /// Albedo
    pub dalb: [Float; NHOR],
    /// Net solar radiation
    pub dswfl: [[Float; NHOR]; NLEP],
    /// Net thermal radiation
    pub dlwfl: [[Float; NHOR]; NLEP],
    /// Net radiation (SW + LW)
    pub dflux: [[Float; NHOR]; NLEP],
    /// Solar radiation upward
    pub dfu: [[Float; NHOR]; NLEP],
    /// Solar radiation downward
    pub dfd: [[Float; NHOR]; NLEP],
    /// Thermal radiation upward
    pub dftu: [[Float; NHOR]; NLEP],
    /// Thermal radiation downward
    pub dftd: [[Float; NHOR]; NLEP],
}

impl Default for BasicRadiation {
    fn default() -> Self {
        Self {
            dalb: [0.; NHOR],
            dswfl: [[0.; NHOR]; NLEP],
            dlwfl: [[0.; NHOR]; NLEP],
            dflux: [[0.; NHOR]; NLEP],
            dfu: [[0.; NHOR]; NLEP],
            dfd: [[0.; NHOR]; NLEP],
            dftu: [[0.; NHOR]; NLEP],
            dftd: [[0.; NHOR]; NLEP],
        }
    }
}

struct Surface {
    /// Surface wetness
    pub drhs: [Float; NHOR],
    /// Land(1)/sea(0) mask
    pub dls: [Float; NHOR],
    /// Rougthness length
    pub dz0: [Float; NHOR],
    /// Ice thickness
    pub diced: [Float; NHOR],
    /// Ice cover
    pub dicec: [Float; NHOR],
    /// x-surface wind stress
    pub dtaux: [Float; NHOR],
    /// y-surface wind stress
    pub dtauy: [Float; NHOR],
    /// u-star**3 (needed eg. for coupling)
    pub dust3: [Float; NHOR],
    /// Surface sensible heat flx
    pub dshfl: [Float; NHOR],
    /// Surface latent heat flx
    pub dlhfl: [Float; NHOR],
    /// Surface evaporation
    pub devap: [Float; NHOR],
    /// Surface air temperature
    pub dtsa: [Float; NHOR],
    /// Mixed-layer depth (output from ocean)
    pub dmld: [Float; NHOR],

    // TODO what is this?
    pub dshdt: [Float; NHOR],
    // TODO what is this?
    pub dlhdt: [Float; NHOR],
}

impl Default for Surface {
    fn default() -> Self {
        Self {
            dls: [1.; NHOR],
            ..Default::default()
        }
    }
}

struct Water {
    /// Convective Precip (m/s)
    pub dprc: [Float; NHOR],
    /// Large Scale Precip (m/s)
    pub dprl: [Float; NHOR],
    /// Snow Fall (m/s)
    pub dprs: [Float; NHOR],
    /// Vertical integrated specific humidity (kg/m**2)
    pub dqvi: [Float; NHOR],
}

impl Default for Water {
    fn default() -> Self {
        Self {
            dprc: [0.; NHOR],
            dprl: [0.; NHOR],
            dprs: [0.; NHOR],
            dqvi: [0.; NHOR],
        }
    }
}

struct Biome {
    /// Forest cover (fract.)
    pub dforest: [Float; NHOR],
    /// Field capacity (m)
    pub dwmax: [Float; NHOR],
}

impl Default for Biome {
    fn default() -> Self {
        Self {
            dforest: [0.5; NHOR],
            ..Default::default()
        }
    }
}

struct Soil {
    /// Soil wetness (m)
    pub dwatc: [Float; NHOR],
    /// Surface runoff (m/s)
    pub drunoff: [Float; NHOR],
    /// Snow depth (m)
    pub dsnow: [Float; NHOR],
    /// Snow melt (m/s water eq.)
    pub dsmelt: [Float; NHOR],
    /// Snow depth change (m/s water eq.)
    pub dsndch: [Float; NHOR],
    /// Soil temperature uppermost level (K)
    pub dtsoil: [Float; NHOR],
    /// Soil temperature level 2 (K)
    pub dtd2: [Float; NHOR],
    /// Soil temperature level 3 (K)
    pub dtd3: [Float; NHOR],
    /// Soil temperature level 4 (K)
    pub dtd4: [Float; NHOR],
    /// Soil temperature lowermost level (K)
    pub dtd5: [Float; NHOR],
    /// Glacier mask (0.,1.)
    pub dglac: [Float; NHOR],
}

impl Default for Soil {
    fn default() -> Self {
        Self {
            dwatc: [0.; NHOR],
            drunoff: [0.; NHOR],
            dsnow: [0.; NHOR],
            dsmelt: [0.; NHOR],
            dsndch: [0.; NHOR],
            dtsoil: [0.; NHOR],
            dtd2: [0.; NHOR],
            dtd3: [0.; NHOR],
            dtd4: [0.; NHOR],
            dtd5: [0.; NHOR],
            dglac: [0.; NHOR],
        }
    }
}

struct DiagnosticArrays {
    // TODO what is this?
    pub ndl: [Int; NLEV],

    pub csu: [[Float; NLAT]; NLEV],
    pub csv: [[Float; NLAT]; NLEV],
    pub cst: [[Float; NLAT]; NLEV],
    pub csm: [[Float; NLAT]; NLEV],
    pub ccc: [[Float; NLAT]; NLEV],
    pub span: [Float; NESP],

    /// 2-d diagnostics
    pub dgp2d: Vec<Vec<Float>>,
    pub dsp2d: Vec<Vec<Float>>,
    /// 3-d diagnostics
    pub dgp3d: Vec<Vec<Vec<Float>>>,
    pub dsp3d: Vec<Vec<Vec<Float>>>,

    /// Cloud forcing diagnostics
    pub dclforc: Vec<Vec<Float>>,
    /// Entropy diagnostics
    pub dentropy: Vec<Vec<Float>>,
    /// Entropy diagnostics 3d
    pub dentro3d: Vec<Vec<Vec<Float>>>,
    /// Energy diagnostics
    pub denergy: Vec<Vec<Float>>,
    /// Energy diagnostics 3d
    pub dener3d: Vec<Vec<Vec<Float>>>,

    /// ps for entropy diagnostics
    pub dentrop: Vec<Float>,
    /// t for entropy diagnostics
    pub dentrot: Vec<Vec<Float>>,
    /// q for entropy diagnostics
    pub dentroq: Vec<Vec<Float>>,
    /// 2d entropy for diagnostics
    pub dentro: Vec<Float>,

    // accumulated output
    /// Acculumated evaporation
    pub aevap: [Float; NHOR],
    /// Acculumated lage scale precip
    pub aprl: [Float; NHOR],
    /// Acculumated convective precip
    pub aprc: [Float; NHOR],
    /// Acculumated snow fall
    pub aprs: [Float; NHOR],
    /// Acculumated sensible heat flux
    pub ashfl: [Float; NHOR],
    /// Acculumated latent heat flux
    pub alhfl: [Float; NHOR],
    /// Acculumated surface runoff
    pub aroff: [Float; NHOR],
    /// Acculumated snow melt
    pub asmelt: [Float; NHOR],
    /// Acculumated snow depth change
    pub asndch: [Float; NHOR],
    /// Acculumated total cloud cover
    pub acc: [Float; NHOR],
    /// Acculumated surface solar radiation
    pub assol: [Float; NHOR],
    /// Acculumated surface thermal radiation
    pub asthr: [Float; NHOR],
    /// Acculumated top solar radiation
    pub atsol: [Float; NHOR],
    /// Acculumated top thermal radiation
    pub atthr: [Float; NHOR],
    /// Acculumated surface solar radiation upward
    pub assolu: [Float; NHOR],
    /// Acculumated surface thermal radiation upward
    pub asthru: [Float; NHOR],
    /// Acculumated top solar radiation upward
    pub atsolu: [Float; NHOR],
    /// Acculumated zonal wind stress
    pub ataux: [Float; NHOR],
    /// Acculumated meridional wind stress
    pub atauy: [Float; NHOR],
    /// Acculumated vertical integrated q
    pub aqvi: [Float; NHOR],
    /// Acculumated surface air temperature
    pub atsa: [Float; NHOR],
    /// Maximum surface air temperature
    pub atsama: [Float; NHOR],
    /// Minimum surface air temperature
    pub atsami: [Float; NHOR],
    /// Acculumated surface temperature
    pub ats0: [Float; NHOR],
}

impl Default for DiagnosticArrays {
    fn default() -> Self {
        Self {
            ndl: [0; NLEV],
            csu: [[0.; NLAT]; NLEV],
            csv: [[0.; NLAT]; NLEV],
            cst: [[0.; NLAT]; NLEV],
            csm: [[0.; NLAT]; NLEV],
            ccc: [[0.; NLAT]; NLEV],
            span: [0.; NESP],
            dgp2d: vec![],
            dsp2d: vec![],
            dgp3d: vec![],
            dsp3d: vec![],
            dclforc: vec![],
            dentropy: vec![],
            dentro3d: vec![],
            denergy: vec![],
            dener3d: vec![],
            dentrop: vec![],
            dentrot: vec![],
            dentroq: vec![],
            dentro: vec![],
            aevap: [0.; NHOR],
            aprl: [0.; NHOR],
            aprc: [0.; NHOR],
            aprs: [0.; NHOR],
            ashfl: [0.; NHOR],
            alhfl: [0.; NHOR],
            aroff: [0.; NHOR],
            asmelt: [0.; NHOR],
            asndch: [0.; NHOR],
            acc: [0.; NHOR],
            assol: [0.; NHOR],
            asthr: [0.; NHOR],
            atsol: [0.; NHOR],
            atthr: [0.; NHOR],
            assolu: [0.; NHOR],
            asthru: [0.; NHOR],
            atsolu: [0.; NHOR],
            ataux: [0.; NHOR],
            atauy: [0.; NHOR],
            aqvi: [0.; NHOR],
            atsa: [0.; NHOR],
            atsama: [0.; NHOR],
            atsami: [0.; NHOR],
            ats0: [0.; NHOR],
        }
    }
}

struct LatitudeArrays {
    pub chlat: [[char; 3]; NLAT],
    //TODO "real kind=8"
    /// sin(phi)
    pub sid: [Float; NLAT],
    //TODO "real kind=8"
    /// Gaussian weights
    pub gwd: [Float; NLAT],
    /// cos(phi)**2
    pub csq: [Float; NLAT],
    /// cos(phi)
    pub cola: [Float; NLAT],
    /// 1 / cos(phi)
    pub rcs: [Float; NLAT],
    /// latitude in degrees
    pub deglat: [Float; NLPP],
}

impl Default for LatitudeArrays {
    fn default() -> Self {
        Self {
            chlat: [['a'; 3]; NLAT],
            sid: [0.; NLAT],
            gwd: [0.; NLAT],
            csq: [0.; NLAT],
            cola: [0.; NLAT],
            rcs: [0.; NLAT],
            deglat: [0.; NLPP],
        }
    }
}

struct LevelArrays {
    /// Diffusion time scale for divergence (days)
    pub tdissd: [Float; NLEV],
    /// Diffusion time scale for vorticity (days)
    pub tdissz: [Float; NLEV],
    /// Diffusion time scale for temperature (days)
    pub tdisst: [Float; NLEV],
    /// Diffusion time scale for sp. humidity (days)
    pub tdissq: [Float; NLEV],

    pub restim: [Float; NLEV],
    pub t0: [Float; NLEV],
    pub tfrc: [Float; NLEV],
    pub sigh: [Float; NLEV],
    pub damp: [Float; NLEV],
    pub dsigma: [Float; NLEV],
    pub sigma: [Float; NLEV],
    pub sigmah: [Float; NLEV],
    pub t01s2: [Float; NLEV],
    pub tkp: [Float; NLEV],

    pub c: [[Float; NLEV]; NLEV],
    pub g: [[Float; NLEV]; NLEV],
    pub tau: [[Float; NLEV]; NLEV],
    pub bm1: [[[Float; NLEV]; NLEV]; NTRU],
}

impl Default for LevelArrays {
    fn default() -> Self {
        Self {
            tdissd: [0.20; NLEV],
            tdissz: [1.10; NLEV],
            tdisst: [5.60; NLEV],
            tdissq: [0.1; NLEV],
            t0: [250.0; NLEV],
            ..Default::default()
        }
    }
}

#[derive(Default)]
struct RandomSeed {
    // TODO replacing namelists with monolithic config?
    /// Settable in namelist
    pub seed: [Int; 8],
    /// Machine dependent seed
    pub meed: Vec<Int>,
}

struct PlanetVars {
    /// Planet name
    pub yplanet: String,
    /// Global switch to fix orbit
    nfixorb: Int,

    /// Kappa
    akap: Float,
    /// Lapse rate
    alr: Float,
    /// Gravity
    ga: Float,
    /// Gas constant for dry air
    gascon: Float,
    /// Planet radius
    plarad: Float,
    /// Planet radius
    pnu: Float,

    /// Length of sidereal day (sec)
    sidereal_day: Float,
    /// Length of solar day (sec)
    solar_day: Float,
    /// Length of sidereal year (sec)
    sidereal_year: Float,
    /// Length of tropical year (sec)
    tropical_year: Float,

    /// Omega used for scaling
    ww: Float,
    /// Orography scaling
    oroscale: Float,

    ra1: Float,
    ra2: Float,
    ra4: Float,

    /// Specific heat for dry air
    /// `acpd = gascon / akap`
    acpd: Float,
    /// `acpv / acpd - 1.0`
    adv: Float,
    /// `cv = plarad * ww`
    cv: Float,
    /// `ct = CV * CV / gascon`
    ct: Float,
    /// Time filter 2
    /// `pnu21 = 1.0 - 2.0 * pnu `
    pnu21: Float,
    /// Rd/Rv
    /// `rdbrv = gascon / RV`
    rdbrv: Float,
    /// Rotation speed (factor)
    rotspd: Float,
    /// Eccentricity of Orbit
    eccen: Float,
    /// Obliquity of Orbit
    obliq: Float,
    /// Longitude of moving vernal equinox
    mvelp: Float,
}

impl Default for PlanetVars {
    fn default() -> Self {
        Self {
            oroscale: 1.0,
            rotspd: 1.0,
            ..Default::default()
        }
    }
}

struct DateTime {
    /// Current timestep
    pub nstep: Int,
    /// Start timestep for this run
    pub nstep1: Int,
    /// Start timestep for this run
    pub mstep: Int,
    /// Start timestep for this run
    pub mocd: Int,

    /// Start year
    pub n_start_year: Int,
    /// Start month
    pub n_start_month: Int,
    /// Start step since 1-Jan-0000
    pub n_start_step: Int,
    /// Needed for time interpolation
    pub n_days_per_month: Int,
    /// Set to 365 for real calendar
    pub n_days_per_year: Int,

    /// Years to run
    pub n_run_years: Int,
    /// Months to run
    pub n_run_months: Int,
    /// Days  to run (debugging)
    pub n_run_days: Int,
    /// Steps to run (debugging)
    pub n_run_steps: Int,

    /// minutes/timestep = 1day/ntspd
    pub mpstep: Int,
    /// Number of timesteps per day
    pub ntspd: Int,
    /// Number of writes per day
    pub nwpd: Int,
    /// Date & time array
    pub ndatim: [Int; 7],
    /// Start of run
    pub tmstart: Float,
}

impl Default for DateTime {
    fn default() -> Self {
        Self {
            n_start_year: 1,
            n_start_month: 1,
            n_days_per_month: 30,
            n_days_per_year: 360,
            ndatim: [-1; 7],
            ..Default::default()
        }
    }
}
