pub mod radiation;

use crate::config::PlanetConfig;
use crate::constants::{NESP, NHOR, NLAT, NLEP, NLEV, NLON, NLPP, NPRO, NSPP, NTRACE, NTRU};
use crate::{Float, Int, Vec2d, Vec3d, Vec4d};

#[derive(Default, Clone)]
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

#[derive(Copy, Clone)]
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
struct SpectralArrays {
    /// Spectral Divergence
    pub sd: Vec2d<Float>, // (NESP, NLEV)
    /// Spectral Temperature
    pub st: Vec2d<Float>, // (NESP, NLEV)
    /// Spectral Vorticity
    pub sz: Vec2d<Float>, // (NESP, NLEV)
    /// Spectral Specific Humidity
    pub sq: Vec2d<Float>, // (NESP, NLEV)
    /// Spectral Pressure (ln Ps)
    pub sp: Vec<Float>, // (NESP)
    /// Spectral Orography
    pub so: Vec<Float>, // (NESP)
    /// Spectral Restoration Temperature
    pub sr: Vec2d<Float>, // (NESP, NLEV)

    /// Spectral Divergence  Partial
    pub sdp: Vec2d<Float>, // (NSPP, NLEV)
    /// Spectral Temperature Partial
    pub stp: Vec2d<Float>, // (NSPP, NLEV)
    /// Spectral Vorticity   Partial
    pub szp: Vec2d<Float>, // (NSPP, NLEV)
    /// Spectral S.Humidity  Partial
    pub sqp: Vec2d<Float>, // (NSPP, NLEV)
    /// Spectral Pressure    Partial
    pub spp: Vec<Float>, // (NSPP)
    /// Spectral Orography   Partial
    pub sop: Vec<Float>, // (NSPP)
    /// Spectral Restoration Partial
    pub srp: Vec2d<Float>, // (NSPP, NLEV)

    /// Spectral Divergence  Tendency
    pub sdt: Vec2d<Float>, // (NSPP, NLEV)
    /// Spectral Temperature Tendency
    pub stt: Vec2d<Float>, // (NSPP, NLEV)
    /// Spectral Vorticity   Tendency
    pub szt: Vec2d<Float>, // (NSPP, NLEV)
    /// Spectral S.Humidity  Tendency
    pub sqt: Vec2d<Float>, // (NSPP, NLEV)
    /// Spectral Pressure    Tendency
    pub spt: Vec<Float>, // (NSPP)

    /// Spectral Divergence  Minus
    pub sdm: Vec2d<Float>, // (NSPP, NLEV)
    /// Spectral Temperature Minus
    pub stm: Vec2d<Float>, // (NSPP, NLEV)
    /// Spectral Vorticity   Minus
    pub szm: Vec2d<Float>, // (NSPP, NLEV)
    /// Spectral S.Humidity  Minus
    pub sqm: Vec2d<Float>, // (NSPP, NLEV)
    /// Spectral Pressure    Minus
    pub spm: Vec<Float>, // (NSPP)

    /// horizontal diffusion
    pub sak: Vec2d<Float>, // (NESP, NLEV)
    /// horizontal diffusion partial
    pub sakpp: Vec2d<Float>, // (NSPP, NLEV)
    /// specific humidity for output
    pub sqout: Vec2d<Float>, // (NESP, NLEV)
    /// Factors for output normalization
    pub spnorm: Vec<Float>, // (NESP)

    /// Holds wavenumber
    pub nindex: Vec<Int>, // (NESP)
    /// Holds wavenumber
    pub nscatsp: Vec<Int>, // (NPRO)
    /// Holds wavenumber
    pub ndel: Vec<Int>, // (NLEV)

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
    pub gd: Vec2d<Float>, // (NHOR, NLEV)
    /// Temperature (-t0)
    pub gt: Vec2d<Float>, // (NHOR, NLEV)
    /// Absolute vorticity
    pub gz: Vec2d<Float>, // (NHOR, NLEV)
    /// Spec. humidity
    pub gq: Vec2d<Float>, // (NHOR, NLEV)
    /// Zonal wind (*cos(phi))
    pub gu: Vec2d<Float>, // (NHOR, NLEV)
    /// Meridional wind (*cos(phi))
    pub gv: Vec2d<Float>, // (NHOR, NLEV)
    /// t-tendency
    pub gtdt: Vec2d<Float>, // (NHOR, NLEV)
    /// q-tendency
    pub gqdt: Vec2d<Float>, // (NHOR, NLEV)
    /// u-tendency
    pub gudt: Vec2d<Float>, // (NHOR, NLEV)
    /// v-tendency
    pub gvdt: Vec2d<Float>, // (NHOR, NLEV)
    /// Surface pressure or ln(ps)
    pub gp: Vec<Float>, // (NHOR)
    /// dln(ps)/dphi
    pub gpj: Vec<Float>, // (NHOR)

    /// 1/cos(phi)**2
    pub rcsq: Vec<Float>, // (NHOR)
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
    pub dt: Vec2d<Float>, // (NHOR, NLEP)
    /// Spec. humidity
    pub dq: Vec2d<Float>, // (NHOR, NLEP)
    /// Zonal wind [m/s]
    pub du: Vec2d<Float>, // (NHOR, NLEP)
    /// Meridional wind [m/s]
    pub dv: Vec2d<Float>, // (NHOR, NLEP)
    /// Surface pressure
    pub dp: Vec<Float>, // (NHOR)

    /// Saturation humidity
    pub dqsat: Vec2d<Float>, // (NHOR, NLEP)
    /// Adiabatic q-tendencies (for eg kuo)
    pub dqt: Vec2d<Float>, // (NHOR, NLEP)
    /// Cloud cover
    pub dcc: Vec2d<Float>, // (NHOR, NLEP)
    /// Liquid water content
    pub dql: Vec2d<Float>, // (NHOR, NLEP)
    /// Vertical velocity (dp/dt)
    pub dw: Vec2d<Float>, // (NHOR, NLEV)
    /// t-tendency
    pub dtdt: Vec2d<Float>, // (NHOR, NLEP)
    /// q-tendency
    pub dqdt: Vec2d<Float>, // (NHOR, NLEP)
    /// u-tendency
    pub dudt: Vec2d<Float>, // (NHOR, NLEP)
    /// v-tendency
    pub dvdt: Vec2d<Float>, // (NHOR, NLEP)

    /// Surface pressure at time t
    pub dp0: Vec<Float>, // (NHOR)
    /// Zonal wind at time t
    pub du0: Vec2d<Float>, // (NHOR, NLEP)
    /// Meridional wind at time t
    pub dv0: Vec2d<Float>, // (NHOR, NLEP)
    /// Trace array
    pub dtrace: Vec4d<Float>, // (NLON, NLAT, NLEV, NTRACE)
}

impl Default for GridArraysDim {
    fn default() -> Self {
        Self {
            dt: vec![],
            dq: vec![],
            du: vec![],
            dv: vec![],
            dp: vec![],
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
    pub dalb: Vec<Float>, // (NHOR)
    /// Net solar radiation
    pub dswfl: Vec2d<Float>, // (NHOR, NLEP)
    /// Net thermal radiation
    pub dlwfl: Vec2d<Float>, // (NHOR, NLEP)
    /// Net radiation (SW + LW)
    pub dflux: Vec2d<Float>, // (NHOR, NLEP)
    /// Solar radiation upward
    pub dfu: Vec2d<Float>, // (NHOR, NLEP)
    /// Solar radiation downward
    pub dfd: Vec2d<Float>, // (NHOR, NLEP)
    /// Thermal radiation upward
    pub dftu: Vec2d<Float>, // (NHOR, NLEP)
    /// Thermal radiation downward
    pub dftd: Vec2d<Float>, // (NHOR, NLEP)
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
    pub drhs: Vec<Float>, // (NHOR)
    /// Land(1)/sea(0) mask
    pub dls: Vec<Float>, // (NHOR)
    /// Rougthness length
    pub dz0: Vec<Float>, // (NHOR)
    /// Ice thickness
    pub diced: Vec<Float>, // (NHOR)
    /// Ice cover
    pub dicec: Vec<Float>, // (NHOR)
    /// x-surface wind stress
    pub dtaux: Vec<Float>, // (NHOR)
    /// y-surface wind stress
    pub dtauy: Vec<Float>, // (NHOR)
    /// u-star**3 (needed eg. for coupling)
    pub dust3: Vec<Float>, // (NHOR)
    /// Surface sensible heat flx
    pub dshfl: Vec<Float>, // (NHOR)
    /// Surface latent heat flx
    pub dlhfl: Vec<Float>, // (NHOR)
    /// Surface evaporation
    pub devap: Vec<Float>, // (NHOR)
    /// Surface air temperature
    pub dtsa: Vec<Float>, // (NHOR)
    /// Mixed-layer depth (output from ocean)
    pub dmld: Vec<Float>, // (NHOR)

    // TODO what is this?
    pub dshdt: Vec<Float>, // (NHOR)
    // TODO what is this?
    pub dlhdt: Vec<Float>, // (NHOR)
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
    pub dprc: Vec<Float>, // (NHOR)
    /// Large Scale Precip (m/s)
    pub dprl: Vec<Float>, // (NHOR)
    /// Snow Fall (m/s)
    pub dprs: Vec<Float>, // (NHOR)
    /// Vertical integrated specific humidity (kg/m**2)
    pub dqvi: Vec<Float>, // (NHOR)
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
    pub dforest: Vec<Float>, // (NHOR)
    /// Field capacity (m)
    pub dwmax: Vec<Float>, // (NHOR)
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
    pub dwatc: Vec<Float>, // (NHOR)
    /// Surface runoff (m/s)
    pub drunoff: Vec<Float>, // (NHOR)
    /// Snow depth (m)
    pub dsnow: Vec<Float>, // (NHOR)
    /// Snow melt (m/s water eq.)
    pub dsmelt: Vec<Float>, // (NHOR)
    /// Snow depth change (m/s water eq.)
    pub dsndch: Vec<Float>, // (NHOR)
    /// Soil temperature uppermost level (K)
    pub dtsoil: Vec<Float>, // (NHOR)
    /// Soil temperature level 2 (K)
    pub dtd2: Vec<Float>, // (NHOR)
    /// Soil temperature level 3 (K)
    pub dtd3: Vec<Float>, // (NHOR)
    /// Soil temperature level 4 (K)
    pub dtd4: Vec<Float>, // (NHOR)
    /// Soil temperature lowermost level (K)
    pub dtd5: Vec<Float>, // (NHOR)
    /// Glacier mask (0.,1.)
    pub dglac: Vec<Float>, // (NHOR)
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

    pub csu: Vec2d<Float>, // (NLAT, NLEV)
    pub csv: Vec2d<Float>, // (NLAT, NLEV)
    pub cst: Vec2d<Float>, // (NLAT, NLEV)
    pub csm: Vec2d<Float>, // (NLAT, NLEV)
    pub ccc: Vec2d<Float>, // (NLAT, NLEV)
    pub span: Vec<Float>,  // (NESP)

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
    pub aevap: Vec<Float>, // (NHOR)
    /// Acculumated lage scale precip
    pub aprl: Vec<Float>, // (NHOR)
    /// Acculumated convective precip
    pub aprc: Vec<Float>, // (NHOR)
    /// Acculumated snow fall
    pub aprs: Vec<Float>, // (NHOR)
    /// Acculumated sensible heat flux
    pub ashfl: Vec<Float>, // (NHOR)
    /// Acculumated latent heat flux
    pub alhfl: Vec<Float>, // (NHOR)
    /// Acculumated surface runoff
    pub aroff: Vec<Float>, // (NHOR)
    /// Acculumated snow melt
    pub asmelt: Vec<Float>, // (NHOR)
    /// Acculumated snow depth change
    pub asndch: Vec<Float>, // (NHOR)
    /// Acculumated total cloud cover
    pub acc: Vec<Float>, // (NHOR)
    /// Acculumated surface solar radiation
    pub assol: Vec<Float>, // (NHOR)
    /// Acculumated surface thermal radiation
    pub asthr: Vec<Float>, // (NHOR)
    /// Acculumated top solar radiation
    pub atsol: Vec<Float>, // (NHOR)
    /// Acculumated top thermal radiation
    pub atthr: Vec<Float>, // (NHOR)
    /// Acculumated surface solar radiation upward
    pub assolu: Vec<Float>, // (NHOR)
    /// Acculumated surface thermal radiation upward
    pub asthru: Vec<Float>, // (NHOR)
    /// Acculumated top solar radiation upward
    pub atsolu: Vec<Float>, // (NHOR)
    /// Acculumated zonal wind stress
    pub ataux: Vec<Float>, // (NHOR)
    /// Acculumated meridional wind stress
    pub atauy: Vec<Float>, // (NHOR)
    /// Acculumated vertical integrated q
    pub aqvi: Vec<Float>, // (NHOR)
    /// Acculumated surface air temperature
    pub atsa: Vec<Float>, // (NHOR)
    /// Maximum surface air temperature
    pub atsama: Vec<Float>, // (NHOR)
    /// Minimum surface air temperature
    pub atsami: Vec<Float>, // (NHOR)
    /// Acculumated surface temperature
    pub ats0: Vec<Float>, // (NHOR)
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
    pub sid: Vec<Float>, // (NLAT)
    //TODO "real kind=8"
    /// Gaussian weights
    pub gwd: Vec<Float>, // (NLAT)
    /// cos(phi)**2
    pub csq: Vec<Float>, // (NLAT)
    /// cos(phi)
    pub cola: Vec<Float>, // (NLAT)
    /// 1 / cos(phi)
    pub rcs: Vec<Float>, // (NLAT)
    /// latitude in degrees
    pub deglat: Vec<Float>, // (NLPP)
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
    pub tdissd: Vec<Float>, // (NLEV)
    /// Diffusion time scale for vorticity (days)
    pub tdissz: Vec<Float>, // (NLEV)
    /// Diffusion time scale for temperature (days)
    pub tdisst: Vec<Float>, // (NLEV)
    /// Diffusion time scale for sp. humidity (days)
    pub tdissq: Vec<Float>, // (NLEV)

    pub restim: Vec<Float>, // (NLEV)
    pub t0: Vec<Float>,     // (NLEV)
    pub tfrc: Vec<Float>,   // (NLEV)
    pub sigh: Vec<Float>,   // (NLEV)
    pub damp: Vec<Float>,   // (NLEV)
    pub dsigma: Vec<Float>, // (NLEV)
    pub sigma: Vec<Float>,  // (NLEV)
    pub sigmah: Vec<Float>, // (NLEV)
    pub t01s2: Vec<Float>,  // (NLEV)
    pub tkp: Vec<Float>,    // (NLEV)

    pub c: Vec2d<Float>,   // (NLEV, NLEV)
    pub g: Vec2d<Float>,   // (NLEV, NLEV)
    pub tau: Vec2d<Float>, // (NLEV, NLEV)
    pub bm1: Vec3d<Float>, // (NLEV, NLEV, NTRU)
}

impl Default for LevelArrays {
    fn default() -> Self {
        Self {
            tdissd: vec![0.20; NLEV],
            tdissz: vec![1.10; NLEV],
            tdisst: vec![5.60; NLEV],
            tdissq: vec![0.1; NLEV],
            restim: vec![],
            t0: vec![250.0; NLEV],
            tfrc: vec![],
            sigh: vec![],
            damp: vec![],
            dsigma: vec![],
            sigma: vec![],
            sigmah: vec![],
            t01s2: vec![],
            tkp: vec![],
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
    pub tmstart: Float,
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
