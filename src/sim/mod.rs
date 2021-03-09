pub mod cal;
pub mod radiation;

use crate::config::PlanetConfig;
use crate::constants::{EZ, NESP, NHOR, NLAT, NLEP, NLEV, NLON, NLPP, NPRO, NSPP, NTRACE, NTRU};
use crate::{FloatNum, Int, Vec2d, Vec3d, Vec4d};

#[derive(Default, Clone)]
pub struct Sim {
    // resolution: Resolution,
    datetime: DateTime,

    int_scalars: IntScalars,
    real_scalars: RealScalars,

    spec_arrays: SpectralArrays,
    grid_arrays_undim: GridArraysUndim,
    grid_arrays_dim: GridArraysDim,
    rad_basic: BasicRadiation,
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
        let mut sim: Sim = Default::default();
        sim.planet_vars.eccen = planet_cfg.eccen;
        sim.planet_vars.obliq = planet_cfg.obliq;
        sim.planet_vars.mvelp = planet_cfg.mvelp;
        sim.planet_vars.nfixorb = planet_cfg.nfixorb;
        println!("init start");
        sim.rad.initialize(
            sim.int_scalars.ndheat,
            sim.int_scalars.neqsig,
            sim.planet_vars.nfixorb,
            sim.datetime.n_start_year,
            sim.real_scalars.co2,
            sim.planet_vars.eccen,
            sim.planet_vars.obliq,
            sim.planet_vars.mvelp,
        );
        println!("init done");
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
        radiation::step(self);
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

#[derive(Clone)]
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
        println!("IntScalars::default()");
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

#[derive(Clone)]
struct RealScalars {
    /// Latent heat of sublimation
    pub als: FloatNum,
    /// Latent heat of vaporization
    pub alv: FloatNum,
    /// planetary vorticity
    pub plavor: FloatNum,
    /// angle threshold for solar radiation
    pub dawn: FloatNum,
    /// timestep (sec)
    pub deltsec: FloatNum,
    /// timestep (sec) * 2
    pub deltsec2: FloatNum,
    /// deltsec * Omega (ww)
    pub delt: FloatNum,
    /// 2 * delt
    pub delt2: FloatNum,

    pub dtep: FloatNum,
    pub dtns: FloatNum,
    pub dtrop: FloatNum,
    pub dttrp: FloatNum,

    /// Temperature ground in mean profile
    pub tgr: FloatNum,
    /// global mean surface pressure
    pub psurf: FloatNum,
    /// start time (for performance estimates)
    pub time0: FloatNum,
    /// atm. co2 concentration (ppmv)
    pub co2: FloatNum,
    /// diagnostic U max
    pub umax: FloatNum,
    /// diagnostic T2m mean
    pub t2mean: FloatNum,
    /// Melting point (H2O)
    pub tmelt: FloatNum,
    /// diagnostic precipitation mean
    pub precip: FloatNum,
    /// diagnostic evaporation
    pub evap: FloatNum,
    /// Outgoing longwave radiation
    pub olr: FloatNum,
    /// damping time (days) for sponge layer
    pub dampsp: FloatNum,
}

impl Default for RealScalars {
    fn default() -> Self {
        Self {
            als: 2.8345E6,
            alv: 2.5008E6,
            plavor: EZ,
            dawn: 0.0,
            deltsec: 0.0,
            deltsec2: 0.0,
            delt: 0.0,
            delt2: 0.0,
            dtep: 0.0,
            dtns: 0.0,
            dtrop: 12000.0,
            dttrp: 2.0,
            tgr: 288.0,
            psurf: 101100.0,
            time0: 0.0,
            co2: 360.0,
            umax: 0.0,
            t2mean: 0.0,
            tmelt: 273.16,
            precip: 0.0,
            evap: 0.0,
            olr: 0.0,
            dampsp: 0.0,
        }
    }
}

#[derive(Clone)]
struct SpectralArrays {
    /// Spectral Divergence
    pub sd: Vec2d<FloatNum>, // (NESP, NLEV)
    /// Spectral Temperature
    pub st: Vec2d<FloatNum>, // (NESP, NLEV)
    /// Spectral Vorticity
    pub sz: Vec2d<FloatNum>, // (NESP, NLEV)
    /// Spectral Specific Humidity
    pub sq: Vec2d<FloatNum>, // (NESP, NLEV)
    /// Spectral Pressure (ln Ps)
    pub sp: Vec<FloatNum>, // (NESP)
    /// Spectral Orography
    pub so: Vec<FloatNum>, // (NESP)
    /// Spectral Restoration Temperature
    pub sr: Vec2d<FloatNum>, // (NESP, NLEV)

    /// Spectral Divergence  Partial
    pub sdp: Vec2d<FloatNum>, // (NSPP, NLEV)
    /// Spectral Temperature Partial
    pub stp: Vec2d<FloatNum>, // (NSPP, NLEV)
    /// Spectral Vorticity   Partial
    pub szp: Vec2d<FloatNum>, // (NSPP, NLEV)
    /// Spectral S.Humidity  Partial
    pub sqp: Vec2d<FloatNum>, // (NSPP, NLEV)
    /// Spectral Pressure    Partial
    pub spp: Vec<FloatNum>, // (NSPP)
    /// Spectral Orography   Partial
    pub sop: Vec<FloatNum>, // (NSPP)
    /// Spectral Restoration Partial
    pub srp: Vec2d<FloatNum>, // (NSPP, NLEV)

    /// Spectral Divergence  Tendency
    pub sdt: Vec2d<FloatNum>, // (NSPP, NLEV)
    /// Spectral Temperature Tendency
    pub stt: Vec2d<FloatNum>, // (NSPP, NLEV)
    /// Spectral Vorticity   Tendency
    pub szt: Vec2d<FloatNum>, // (NSPP, NLEV)
    /// Spectral S.Humidity  Tendency
    pub sqt: Vec2d<FloatNum>, // (NSPP, NLEV)
    /// Spectral Pressure    Tendency
    pub spt: Vec<FloatNum>, // (NSPP)

    /// Spectral Divergence  Minus
    pub sdm: Vec2d<FloatNum>, // (NSPP, NLEV)
    /// Spectral Temperature Minus
    pub stm: Vec2d<FloatNum>, // (NSPP, NLEV)
    /// Spectral Vorticity   Minus
    pub szm: Vec2d<FloatNum>, // (NSPP, NLEV)
    /// Spectral S.Humidity  Minus
    pub sqm: Vec2d<FloatNum>, // (NSPP, NLEV)
    /// Spectral Pressure    Minus
    pub spm: Vec<FloatNum>, // (NSPP)

    /// horizontal diffusion
    pub sak: Vec2d<FloatNum>, // (NESP, NLEV)
    /// horizontal diffusion partial
    pub sakpp: Vec2d<FloatNum>, // (NSPP, NLEV)
    /// specific humidity for output
    pub sqout: Vec2d<FloatNum>, // (NESP, NLEV)
    /// Factors for output normalization
    pub spnorm: Vec<FloatNum>, // (NESP)

    /// Holds wavenumber
    pub nindex: Vec<Int>, // (NESP)
    /// Holds wavenumber
    pub nscatsp: Vec<Int>, // (NPRO)
    /// Holds wavenumber
    pub ndel: Vec<Int>, // (NLEV)

    /// Difference between instances
    pub sdd: Vec<Vec<FloatNum>>,
    /// Difference between instances
    pub std: Vec<Vec<FloatNum>>,
    /// Difference between instances
    pub szd: Vec<Vec<FloatNum>>,
    /// Difference between instances
    pub spd: Vec<FloatNum>,
}

impl Default for SpectralArrays {
    fn default() -> Self {
        Self {
            sd: vec![],
            st: vec![],
            sz: vec![],
            sq: vec![],
            sp: vec![],
            so: vec![],
            sr: vec![],
            sdp: vec![],
            stp: vec![],
            szp: vec![],
            sqp: vec![],
            spp: vec![],
            sop: vec![],
            srp: vec![],
            sdt: vec![],
            stt: vec![],
            szt: vec![],
            sqt: vec![],
            spt: vec![],
            sdm: vec![],
            stm: vec![],
            szm: vec![],
            sqm: vec![],
            spm: vec![],
            sak: vec![],
            sakpp: vec![],
            sqout: vec![],
            spnorm: vec![],
            nindex: vec![NTRU as Int; NESP],
            nscatsp: vec![NSPP as Int; NPRO],
            ndel: vec![2; NLEV],
            sdd: vec![],
            std: vec![],
            szd: vec![],
            spd: vec![],
        }
    }
}

/// Global Gridpoint Arrays (un-dimensionalized)
#[derive(Default, Clone)]
struct GridArraysUndim {
    /// Divergence
    pub gd: Vec2d<FloatNum>, // (NHOR, NLEV)
    /// Temperature (-t0)
    pub gt: Vec2d<FloatNum>, // (NHOR, NLEV)
    /// Absolute vorticity
    pub gz: Vec2d<FloatNum>, // (NHOR, NLEV)
    /// Spec. humidity
    pub gq: Vec2d<FloatNum>, // (NHOR, NLEV)
    /// Zonal wind (*cos(phi))
    pub gu: Vec2d<FloatNum>, // (NHOR, NLEV)
    /// Meridional wind (*cos(phi))
    pub gv: Vec2d<FloatNum>, // (NHOR, NLEV)
    /// t-tendency
    pub gtdt: Vec2d<FloatNum>, // (NHOR, NLEV)
    /// q-tendency
    pub gqdt: Vec2d<FloatNum>, // (NHOR, NLEV)
    /// u-tendency
    pub gudt: Vec2d<FloatNum>, // (NHOR, NLEV)
    /// v-tendency
    pub gvdt: Vec2d<FloatNum>, // (NHOR, NLEV)
    /// Surface pressure or ln(ps)
    pub gp: Vec<FloatNum>, // (NHOR)
    /// dln(ps)/dphi
    pub gpj: Vec<FloatNum>, // (NHOR)

    /// 1/cos(phi)**2
    pub rcsq: Vec<FloatNum>, // (NHOR)
}

// impl Default for GridArraysUndim {
//     fn default() -> Self {
//         Self {
//             gd: vec![vec![0.; NLEV]; NHOR],
//             gt: vec![vec![0.; NLEV]; NHOR],
//             gz: vec![vec![0.; NLEV]; NHOR],
//             gq: vec![vec![0.; NLEV]; NHOR],
//             gu: vec![vec![0.; NLEV]; NHOR],
//             gv: vec![vec![0.; NLEV]; NHOR],
//             gtdt: vec![vec![0.; NLEV]; NHOR],
//             gqdt: vec![vec![0.; NLEV]; NHOR],
//             gudt: vec![vec![0.; NLEV]; NHOR],
//             gvdt: vec![vec![0.; NLEV]; NHOR],
//             gp: vec![0.; NHOR],
//             gpj: vec![0.; NHOR],
//             rcsq: vec![0.; NHOR],
//         }
//     }
// }

/// Global Gridpoint Arrays (dimensionalized)
#[derive(Clone)]
struct GridArraysDim {
    /// Temperature
    pub dt: Vec2d<FloatNum>, // (NHOR, NLEP)
    /// Spec. humidity
    pub dq: Vec2d<FloatNum>, // (NHOR, NLEP)
    /// Zonal wind [m/s]
    pub du: Vec2d<FloatNum>, // (NHOR, NLEP)
    /// Meridional wind [m/s]
    pub dv: Vec2d<FloatNum>, // (NHOR, NLEP)
    /// Surface pressure
    pub dp: Vec<FloatNum>, // (NHOR)

    /// Saturation humidity
    pub dqsat: Vec2d<FloatNum>, // (NHOR, NLEP)
    /// Adiabatic q-tendencies (for eg kuo)
    pub dqt: Vec2d<FloatNum>, // (NHOR, NLEP)
    /// Cloud cover
    pub dcc: Vec2d<FloatNum>, // (NHOR, NLEP)
    /// Liquid water content
    pub dql: Vec2d<FloatNum>, // (NHOR, NLEP)
    /// Vertical velocity (dp/dt)
    pub dw: Vec2d<FloatNum>, // (NHOR, NLEV)
    /// t-tendency
    pub dtdt: Vec2d<FloatNum>, // (NHOR, NLEP)
    /// q-tendency
    pub dqdt: Vec2d<FloatNum>, // (NHOR, NLEP)
    /// u-tendency
    pub dudt: Vec2d<FloatNum>, // (NHOR, NLEP)
    /// v-tendency
    pub dvdt: Vec2d<FloatNum>, // (NHOR, NLEP)

    /// Surface pressure at time t
    pub dp0: Vec<FloatNum>, // (NHOR)
    /// Zonal wind at time t
    pub du0: Vec2d<FloatNum>, // (NHOR, NLEP)
    /// Meridional wind at time t
    pub dv0: Vec2d<FloatNum>, // (NHOR, NLEP)
    /// Trace array
    pub dtrace: Vec4d<FloatNum>, // (NLON, NLAT, NLEV, NTRACE)
}

impl Default for GridArraysDim {
    fn default() -> Self {
        Self {
            dt: vec![vec![0.0; NLEP]; NHOR],
            dq: vec![vec![0.0; NLEP]; NHOR],
            du: vec![vec![0.0; NLEP]; NHOR],
            dv: vec![vec![0.0; NLEP]; NHOR],
            dp: vec![0.0; NHOR],
            dqsat: vec![],
            dqt: vec![],
            dcc: vec![],
            dql: vec![],
            dw: vec![],
            dtdt: vec![],
            dqdt: vec![],
            dudt: vec![],
            dvdt: vec![],
            dp0: vec![],
            du0: vec![],
            dv0: vec![],
            dtrace: vec![vec![vec![vec![1.; NTRACE]; NLEV]; NLAT]; NLON],
        }
    }
}

#[derive(Clone)]
struct BasicRadiation {
    /// Albedo
    pub dalb: Vec<FloatNum>, // (NHOR)
    /// Net solar radiation
    pub dswfl: Vec2d<FloatNum>, // (NHOR, NLEP)
    /// Net thermal radiation
    pub dlwfl: Vec2d<FloatNum>, // (NHOR, NLEP)
    /// Net radiation (SW + LW)
    pub dflux: Vec2d<FloatNum>, // (NHOR, NLEP)
    /// Solar radiation upward
    pub dfu: Vec2d<FloatNum>, // (NHOR, NLEP)
    /// Solar radiation downward
    pub dfd: Vec2d<FloatNum>, // (NHOR, NLEP)
    /// Thermal radiation upward
    pub dftu: Vec2d<FloatNum>, // (NHOR, NLEP)
    /// Thermal radiation downward
    pub dftd: Vec2d<FloatNum>, // (NHOR, NLEP)
}

impl Default for BasicRadiation {
    fn default() -> Self {
        Self {
            dalb: vec![0.; NHOR],
            dswfl: vec![vec![0.; NLEP]; NHOR],
            dlwfl: vec![vec![0.; NLEP]; NHOR],
            dflux: vec![vec![0.; NLEP]; NHOR],
            dfu: vec![vec![0.; NLEP]; NHOR],
            dfd: vec![vec![0.; NLEP]; NHOR],
            dftu: vec![vec![0.; NLEP]; NHOR],
            dftd: vec![vec![0.; NLEP]; NHOR],
        }
    }
}

#[derive(Clone)]
struct Surface {
    /// Surface wetness
    pub drhs: Vec<FloatNum>, // (NHOR)
    /// Land(1)/sea(0) mask
    pub dls: Vec<FloatNum>, // (NHOR)
    /// Rougthness length
    pub dz0: Vec<FloatNum>, // (NHOR)
    /// Ice thickness
    pub diced: Vec<FloatNum>, // (NHOR)
    /// Ice cover
    pub dicec: Vec<FloatNum>, // (NHOR)
    /// x-surface wind stress
    pub dtaux: Vec<FloatNum>, // (NHOR)
    /// y-surface wind stress
    pub dtauy: Vec<FloatNum>, // (NHOR)
    /// u-star**3 (needed eg. for coupling)
    pub dust3: Vec<FloatNum>, // (NHOR)
    /// Surface sensible heat flx
    pub dshfl: Vec<FloatNum>, // (NHOR)
    /// Surface latent heat flx
    pub dlhfl: Vec<FloatNum>, // (NHOR)
    /// Surface evaporation
    pub devap: Vec<FloatNum>, // (NHOR)
    /// Surface air temperature
    pub dtsa: Vec<FloatNum>, // (NHOR)
    /// Mixed-layer depth (output from ocean)
    pub dmld: Vec<FloatNum>, // (NHOR)

    // TODO what is this?
    pub dshdt: Vec<FloatNum>, // (NHOR)
    // TODO what is this?
    pub dlhdt: Vec<FloatNum>, // (NHOR)
}

impl Default for Surface {
    fn default() -> Self {
        Self {
            drhs: vec![],
            dls: vec![1.; NHOR],
            dz0: vec![],
            diced: vec![],
            dicec: vec![],
            dtaux: vec![],
            dtauy: vec![],
            dust3: vec![],
            dshfl: vec![],
            dlhfl: vec![],
            devap: vec![],
            dtsa: vec![],
            dmld: vec![],
            dshdt: vec![],
            dlhdt: vec![],
        }
    }
}

#[derive(Clone)]
struct Water {
    /// Convective Precip (m/s)
    pub dprc: Vec<FloatNum>, // (NHOR)
    /// Large Scale Precip (m/s)
    pub dprl: Vec<FloatNum>, // (NHOR)
    /// Snow Fall (m/s)
    pub dprs: Vec<FloatNum>, // (NHOR)
    /// Vertical integrated specific humidity (kg/m**2)
    pub dqvi: Vec<FloatNum>, // (NHOR)
}

impl Default for Water {
    fn default() -> Self {
        Self {
            dprc: vec![0.; NHOR],
            dprl: vec![0.; NHOR],
            dprs: vec![0.; NHOR],
            dqvi: vec![0.; NHOR],
        }
    }
}

#[derive(Clone)]
struct Biome {
    /// Forest cover (fract.)
    pub dforest: Vec<FloatNum>, // (NHOR)
    /// Field capacity (m)
    pub dwmax: Vec<FloatNum>, // (NHOR)
}

impl Default for Biome {
    fn default() -> Self {
        Self {
            dforest: vec![0.5; NHOR],
            dwmax: vec![],
        }
    }
}

#[derive(Clone)]
struct Soil {
    /// Soil wetness (m)
    pub dwatc: Vec<FloatNum>, // (NHOR)
    /// Surface runoff (m/s)
    pub drunoff: Vec<FloatNum>, // (NHOR)
    /// Snow depth (m)
    pub dsnow: Vec<FloatNum>, // (NHOR)
    /// Snow melt (m/s water eq.)
    pub dsmelt: Vec<FloatNum>, // (NHOR)
    /// Snow depth change (m/s water eq.)
    pub dsndch: Vec<FloatNum>, // (NHOR)
    /// Soil temperature uppermost level (K)
    pub dtsoil: Vec<FloatNum>, // (NHOR)
    /// Soil temperature level 2 (K)
    pub dtd2: Vec<FloatNum>, // (NHOR)
    /// Soil temperature level 3 (K)
    pub dtd3: Vec<FloatNum>, // (NHOR)
    /// Soil temperature level 4 (K)
    pub dtd4: Vec<FloatNum>, // (NHOR)
    /// Soil temperature lowermost level (K)
    pub dtd5: Vec<FloatNum>, // (NHOR)
    /// Glacier mask (0.,1.)
    pub dglac: Vec<FloatNum>, // (NHOR)
}

impl Default for Soil {
    fn default() -> Self {
        Self {
            dwatc: vec![0.; NHOR],
            drunoff: vec![0.; NHOR],
            dsnow: vec![0.; NHOR],
            dsmelt: vec![0.; NHOR],
            dsndch: vec![0.; NHOR],
            dtsoil: vec![0.; NHOR],
            dtd2: vec![0.; NHOR],
            dtd3: vec![0.; NHOR],
            dtd4: vec![0.; NHOR],
            dtd5: vec![0.; NHOR],
            dglac: vec![0.; NHOR],
        }
    }
}

#[derive(Clone)]
struct DiagnosticArrays {
    // TODO what is this?
    pub ndl: Vec<Int>, // (NLEV)

    pub csu: Vec2d<FloatNum>, // (NLAT, NLEV)
    pub csv: Vec2d<FloatNum>, // (NLAT, NLEV)
    pub cst: Vec2d<FloatNum>, // (NLAT, NLEV)
    pub csm: Vec2d<FloatNum>, // (NLAT, NLEV)
    pub ccc: Vec2d<FloatNum>, // (NLAT, NLEV)
    pub span: Vec<FloatNum>,  // (NESP)

    /// 2-d diagnostics
    pub dgp2d: Vec<Vec<FloatNum>>,
    pub dsp2d: Vec<Vec<FloatNum>>,
    /// 3-d diagnostics
    pub dgp3d: Vec<Vec<Vec<FloatNum>>>,
    pub dsp3d: Vec<Vec<Vec<FloatNum>>>,

    /// Cloud forcing diagnostics
    pub dclforc: Vec<Vec<FloatNum>>,
    /// Entropy diagnostics
    pub dentropy: Vec<Vec<FloatNum>>,
    /// Entropy diagnostics 3d
    pub dentro3d: Vec<Vec<Vec<FloatNum>>>,
    /// Energy diagnostics
    pub denergy: Vec<Vec<FloatNum>>,
    /// Energy diagnostics 3d
    pub dener3d: Vec<Vec<Vec<FloatNum>>>,

    /// ps for entropy diagnostics
    pub dentrop: Vec<FloatNum>,
    /// t for entropy diagnostics
    pub dentrot: Vec<Vec<FloatNum>>,
    /// q for entropy diagnostics
    pub dentroq: Vec<Vec<FloatNum>>,
    /// 2d entropy for diagnostics
    pub dentro: Vec<FloatNum>,

    // accumulated output
    /// Acculumated evaporation
    pub aevap: Vec<FloatNum>, // (NHOR)
    /// Acculumated lage scale precip
    pub aprl: Vec<FloatNum>, // (NHOR)
    /// Acculumated convective precip
    pub aprc: Vec<FloatNum>, // (NHOR)
    /// Acculumated snow fall
    pub aprs: Vec<FloatNum>, // (NHOR)
    /// Acculumated sensible heat flux
    pub ashfl: Vec<FloatNum>, // (NHOR)
    /// Acculumated latent heat flux
    pub alhfl: Vec<FloatNum>, // (NHOR)
    /// Acculumated surface runoff
    pub aroff: Vec<FloatNum>, // (NHOR)
    /// Acculumated snow melt
    pub asmelt: Vec<FloatNum>, // (NHOR)
    /// Acculumated snow depth change
    pub asndch: Vec<FloatNum>, // (NHOR)
    /// Acculumated total cloud cover
    pub acc: Vec<FloatNum>, // (NHOR)
    /// Acculumated surface solar radiation
    pub assol: Vec<FloatNum>, // (NHOR)
    /// Acculumated surface thermal radiation
    pub asthr: Vec<FloatNum>, // (NHOR)
    /// Acculumated top solar radiation
    pub atsol: Vec<FloatNum>, // (NHOR)
    /// Acculumated top thermal radiation
    pub atthr: Vec<FloatNum>, // (NHOR)
    /// Acculumated surface solar radiation upward
    pub assolu: Vec<FloatNum>, // (NHOR)
    /// Acculumated surface thermal radiation upward
    pub asthru: Vec<FloatNum>, // (NHOR)
    /// Acculumated top solar radiation upward
    pub atsolu: Vec<FloatNum>, // (NHOR)
    /// Acculumated zonal wind stress
    pub ataux: Vec<FloatNum>, // (NHOR)
    /// Acculumated meridional wind stress
    pub atauy: Vec<FloatNum>, // (NHOR)
    /// Acculumated vertical integrated q
    pub aqvi: Vec<FloatNum>, // (NHOR)
    /// Acculumated surface air temperature
    pub atsa: Vec<FloatNum>, // (NHOR)
    /// Maximum surface air temperature
    pub atsama: Vec<FloatNum>, // (NHOR)
    /// Minimum surface air temperature
    pub atsami: Vec<FloatNum>, // (NHOR)
    /// Acculumated surface temperature
    pub ats0: Vec<FloatNum>, // (NHOR)
}

impl Default for DiagnosticArrays {
    fn default() -> Self {
        Self {
            ndl: vec![0; NLEV],
            csu: vec![vec![0.; NLEV]; NLAT],
            csv: vec![vec![0.; NLEV]; NLAT],
            cst: vec![vec![0.; NLEV]; NLAT],
            csm: vec![vec![0.; NLEV]; NLAT],
            ccc: vec![vec![0.; NLEV]; NLAT],
            span: vec![0.; NESP],
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
            aevap: vec![0.; NHOR],
            aprl: vec![0.; NHOR],
            aprc: vec![0.; NHOR],
            aprs: vec![0.; NHOR],
            ashfl: vec![0.; NHOR],
            alhfl: vec![0.; NHOR],
            aroff: vec![0.; NHOR],
            asmelt: vec![0.; NHOR],
            asndch: vec![0.; NHOR],
            acc: vec![0.; NHOR],
            assol: vec![0.; NHOR],
            asthr: vec![0.; NHOR],
            atsol: vec![0.; NHOR],
            atthr: vec![0.; NHOR],
            assolu: vec![0.; NHOR],
            asthru: vec![0.; NHOR],
            atsolu: vec![0.; NHOR],
            ataux: vec![0.; NHOR],
            atauy: vec![0.; NHOR],
            aqvi: vec![0.; NHOR],
            atsa: vec![0.; NHOR],
            atsama: vec![0.; NHOR],
            atsami: vec![0.; NHOR],
            ats0: vec![0.; NHOR],
        }
    }
}

#[derive(Clone)]
struct LatitudeArrays {
    pub chlat: Vec<String>, // [[char; 3]; NLAT],
    //TODO "real kind=8"
    /// sin(phi)
    pub sid: Vec<FloatNum>, // (NLAT)
    //TODO "real kind=8"
    /// Gaussian weights
    pub gwd: Vec<FloatNum>, // (NLAT)
    /// cos(phi)**2
    pub csq: Vec<FloatNum>, // (NLAT)
    /// cos(phi)
    pub cola: Vec<FloatNum>, // (NLAT)
    /// 1 / cos(phi)
    pub rcs: Vec<FloatNum>, // (NLAT)
    /// latitude in degrees
    pub deglat: Vec<FloatNum>, // (NLPP)
}

impl Default for LatitudeArrays {
    fn default() -> Self {
        Self {
            chlat: vec!["abc".to_string(); NLAT],
            sid: vec![0.; NLAT],
            gwd: vec![0.; NLAT],
            csq: vec![0.; NLAT],
            cola: vec![0.; NLAT],
            rcs: vec![0.; NLAT],
            deglat: vec![0.; NLPP],
        }
    }
}

#[derive(Clone)]
struct LevelArrays {
    /// Diffusion time scale for divergence (days)
    pub tdissd: Vec<FloatNum>, // (NLEV)
    /// Diffusion time scale for vorticity (days)
    pub tdissz: Vec<FloatNum>, // (NLEV)
    /// Diffusion time scale for temperature (days)
    pub tdisst: Vec<FloatNum>, // (NLEV)
    /// Diffusion time scale for sp. humidity (days)
    pub tdissq: Vec<FloatNum>, // (NLEV)

    pub restim: Vec<FloatNum>, // (NLEV)
    pub t0: Vec<FloatNum>,     // (NLEV)
    pub tfrc: Vec<FloatNum>,   // (NLEV)
    pub sigh: Vec<FloatNum>,   // (NLEV)
    pub damp: Vec<FloatNum>,   // (NLEV)
    pub dsigma: Vec<FloatNum>, // (NLEV)
    pub sigma: Vec<FloatNum>,  // (NLEV)
    pub sigmah: Vec<FloatNum>, // (NLEV)
    pub t01s2: Vec<FloatNum>,  // (NLEV)
    pub tkp: Vec<FloatNum>,    // (NLEV)

    pub c: Vec2d<FloatNum>,   // (NLEV, NLEV)
    pub g: Vec2d<FloatNum>,   // (NLEV, NLEV)
    pub tau: Vec2d<FloatNum>, // (NLEV, NLEV)
    pub bm1: Vec3d<FloatNum>, // (NLEV, NLEV, NTRU)
}

impl Default for LevelArrays {
    fn default() -> Self {
        Self {
            tdissd: vec![0.20; NLEV],
            tdissz: vec![1.10; NLEV],
            tdisst: vec![5.60; NLEV],
            tdissq: vec![0.1; NLEV],
            restim: vec![0.0; NLEV],
            t0: vec![250.0; NLEV],
            tfrc: vec![0.0; NLEV],
            sigh: vec![0.0; NLEV],
            damp: vec![0.0; NLEV],
            dsigma: vec![0.0; NLEV],
            sigma: vec![0.0; NLEV],
            sigmah: vec![0.0; NLEV],
            t01s2: vec![0.0; NLEV],
            tkp: vec![0.0; NLEV],
            c: vec![],
            g: vec![],
            tau: vec![],
            bm1: vec![],
        }
    }
}

#[derive(Default, Clone)]
struct RandomSeed {
    // TODO replacing namelists with monolithic config?
    /// Settable in namelist
    pub seed: [Int; 8],
    /// Machine dependent seed
    pub meed: Vec<Int>,
}

#[derive(Clone)]
struct PlanetVars {
    /// Planet name
    pub yplanet: String,
    /// Global switch to fix orbit
    nfixorb: Int,

    /// Kappa
    akap: FloatNum,
    /// Lapse rate
    alr: FloatNum,
    /// Gravity
    ga: FloatNum,
    /// Gas constant for dry air
    gascon: FloatNum,
    /// Planet radius
    plarad: FloatNum,
    /// Planet radius
    pnu: FloatNum,

    /// Length of sidereal day (sec)
    sidereal_day: FloatNum,
    /// Length of solar day (sec)
    solar_day: FloatNum,
    /// Length of sidereal year (sec)
    sidereal_year: FloatNum,
    /// Length of tropical year (sec)
    tropical_year: FloatNum,

    /// Omega used for scaling
    ww: FloatNum,
    /// Orography scaling
    oroscale: FloatNum,

    ra1: FloatNum,
    ra2: FloatNum,
    ra4: FloatNum,

    /// Specific heat for dry air
    /// `acpd = gascon / akap`
    acpd: FloatNum,
    /// `acpv / acpd - 1.0`
    adv: FloatNum,
    /// `cv = plarad * ww`
    cv: FloatNum,
    /// `ct = CV * CV / gascon`
    ct: FloatNum,
    /// Time filter 2
    /// `pnu21 = 1.0 - 2.0 * pnu `
    pnu21: FloatNum,
    /// Rd/Rv
    /// `rdbrv = gascon / RV`
    rdbrv: FloatNum,
    /// Rotation speed (factor)
    rotspd: FloatNum,
    /// Eccentricity of Orbit
    eccen: FloatNum,
    /// Obliquity of Orbit
    obliq: FloatNum,
    /// Longitude of moving vernal equinox
    mvelp: FloatNum,
}

impl Default for PlanetVars {
    fn default() -> Self {
        Self {
            yplanet: "".to_string(),
            nfixorb: 0,
            akap: 0.0,
            alr: 0.0,
            ga: 0.0,
            gascon: 0.0,
            plarad: 0.0,
            pnu: 0.0,
            sidereal_day: 0.0,
            solar_day: 0.0,
            sidereal_year: 0.0,
            tropical_year: 0.0,
            ww: 0.0,
            oroscale: 1.0,
            ra1: 0.0,
            ra2: 0.0,
            ra4: 0.0,
            acpd: 0.0,
            adv: 0.0,
            cv: 0.0,
            ct: 0.0,
            pnu21: 0.0,
            rdbrv: 0.0,
            rotspd: 1.0,
            eccen: 0.0,
            obliq: 0.0,
            mvelp: 0.0,
        }
    }
}

#[derive(Copy, Clone)]
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
    pub tmstart: FloatNum,
}

impl Default for DateTime {
    fn default() -> Self {
        println!("DateTime::default()");
        Self {
            nstep: 0,
            nstep1: 0,
            mstep: 0,
            mocd: 0,
            n_start_year: 1,
            n_start_month: 1,
            n_start_step: 0,
            n_days_per_month: 30,
            n_days_per_year: 360,
            n_run_years: 0,
            n_run_months: 0,
            n_run_days: 0,
            n_run_steps: 0,
            mpstep: 0,
            ntspd: 0,
            nwpd: 0,
            ndatim: [-1; 7],
            tmstart: 0.0,
        }
    }
}
