//! Radiation module, equivalent of radmod.

use crate::constants::{NHOR, NLEP, NLEV, NLON, NLPP, NTRU};
use crate::{Float, Int, Sim, Vec2d};

/// Minimum value for eccen
pub const ORB_ECCEN_MIN: Float = 0.0;
/// Maximum value for eccen
pub const ORB_ECCEN_MAX: Float = 0.1;
/// Minimum value for obliq
pub const ORB_OBLIQ_MIN: Float = -90.0;
/// Maximum value for obliq
pub const ORB_OBLIQ_MAX: Float = 90.0;
/// Minimum value for mvelp
pub const ORB_MVELP_MIN: Float = 0.0;
/// Maximum value for mvelp
pub const ORB_MVELP_MAX: Float = 360.0;
/// undefined/unset/invalid value
pub const ORB_UNDEF_REAL: Float = 1.; // 1.e36;
/// flag to use default orbit
pub const ORB_DEFAULT: Float = ORB_UNDEF_REAL;
/// undefined/unset/invalid value
pub const ORB_UNDEF_INT: Int = 2000000000;
/// flag to not use input year
pub const ORB_NOT_YEAR_BASED: Int = ORB_UNDEF_INT;

#[derive(Clone)]
pub struct Radiation {
    // configurables
    /// Solar constant (set in planet module)
    pub gsol0: Float,
    /// Cos of lat of insolation if ncstsol=1
    pub solclat: Float,
    /// Cos of dec of insolation if ncstsol=1
    pub solcdec: Float,
    /// Cloud grayness (-1 = computed)
    pub clgray: Float,
    /// Absorption coefficient h2o continuum (lwr)
    pub th2oc: Float,
    /// Tuning of cloud albedo range1
    pub tswr1: Float,
    /// Tuning of cloud back scattering c. range2
    pub tswr2: Float,
    /// Tuning of cloud s. scattering alb. range2
    pub tswr3: Float,

    /// Solar constant (set in planet module)
    pub tpofmt: Float,
    /// Solar constant (set in planet module)
    pub acllwr: Float,

    /// Parameter to define o3 profile
    pub a0o3: Float,
    /// Parameter to define o3 profile
    pub a1o3: Float,
    /// Parameter to define o3 profile
    pub aco3: Float,
    /// Parameter to define o3 profile
    pub bo3: Float,
    /// Parameter to define o3 profile
    pub co3: Float,
    /// Parameter to define o3 profile
    pub toffo3: Float,
    /// Scale o3 concentration
    pub o3scale: Float,
    /// Switch for ozon (0=no,1=yes,2=datafile)
    pub no3: Int,
    /// Switch for solang (1/0=yes/no)
    pub nsol: Int,

    /// Switch for swr (1/0=yes/no)
    pub nswr: Int,
    /// Switch for lwr (1/0=yes/no)
    pub nlwr: Int,
    /// Switch for computed cloud props.(1/0=y/n)
    pub nswrcl: Int,
    /// Switch for rayleigh scat. (1/0=yes/no)
    pub nrscat: Int,
    /// Switch for daily cycle of insolation (0 = daily mean insolation)
    pub ndcycle: Int,
    /// Switch to set constant insolation on the whole planet (0/1)=(off/on)
    pub ncstsol: Int,
    /// Year before present (1950 AD), default = 2000 AD
    pub iyrbp: Int,

    /// Cloud albedos spectral range 1
    pub rcl1: [Float; 3],
    /// Cloud albedos spectral range 2
    pub rcl2: [Float; 3],
    /// Cloud absorptivities spectral range 2
    pub acl2: [Float; 3],

    // arrays
    /// Cosine of solar zenit angle
    pub gmu0: Vec<Float>, // NHOR
    /// Cosine of solar zenit angle
    pub gmu1: Vec<Float>, // NHOR
    /// Ozon concentration (kg/kg)
    pub dqo3: Vec2d<Float>, // (NHOR, NLEV)
    /// co2 concentration (ppmv)
    pub dqco2: Vec2d<Float>, // (NHOR, NLEV)
    /// lwr temperature tendencies
    pub dtdtlwr: Vec2d<Float>, // (NHOR, NLEV)
    /// swr temperature tendencies
    pub dtdtswr: Vec2d<Float>, // (NHOR, NLEV)
    /// Climatological O3 (used if NO3=2)
    pub dqo3cl: Vec<Vec<Vec<Float>>>,

    // scalars
    /// Earth-sun distance factor ( i.e. (1/r)**2 )
    pub gdist2: Float,
    /// CPU time for radiation
    pub time4rad: Float,
    /// CPU time for short wave radiation
    pub time4swr: Float,
    /// CPU time for long wave radiation
    pub time4lwr: Float,

    // orbital parameters
    pub ORB_UNDEF_INT: Int,
    /// Earth's obliquity in radians
    pub obliqr: Float,
    /// Mean longitude of perihelion at the vernal equinox (radians)
    pub lambm0: Float,
    /// Earth's moving vernal equinox longitude of perihelion plus pi (radians)
    pub mvelpp: Float,
    /// Earth-sun distance factor ( i.e. (1/r)**2 )
    pub eccf: Float,
    /// Year AD to calculate orbit for
    pub iyrad: Int,
    /// Flag to print-out status information or not.
    /// (This turns off ALL status printing including error messages)
    pub log_print: bool,

    // extended entropy/energy diagnostics
    pub dftde1: Vec2d<Float>, // (NHOR, NLEP)
    pub dftde2: Vec2d<Float>, // (NHOR, NLEP)
    pub dftue1: Vec2d<Float>, // (NHOR, NLEP)
    pub dftue2: Vec2d<Float>, // (NHOR, NLEP)
    pub dftu0: Vec2d<Float>,  // (NHOR, NLEP)
    pub dft: Vec2d<Float>,    // (NHOR, NLEP)

    // auxiliary variables for solar zenit angle calculations
    /// cos(lat)*cos(decl)
    pub solclatcdec: Float,
    /// sin(lat)
    pub solslat: Float,
    /// sin(decl)
    pub solsdec: Float,
    /// sin(lat)*sin(decl)
    pub solslatsdec: Float,
    /// temporary zenit angle
    pub zmuz: Float,
}

impl Default for Radiation {
    // values for planet earth
    fn default() -> Self {
        Self {
            gsol0: 1367.0,
            solclat: 1.0,
            solcdec: 1.0,
            clgray: -1.0,
            th2oc: 0.024,
            tswr1: 0.077,
            tswr2: 0.065,
            tswr3: 0.0055,
            tpofmt: 1.00,
            acllwr: 0.100,
            a0o3: 0.25,
            a1o3: 0.11,
            aco3: 0.08,
            bo3: 20000.,
            co3: 5000.,
            toffo3: 90.,
            o3scale: 1.0,
            no3: 1,
            nsol: 1,
            nswr: 1,
            nlwr: 1,
            nswrcl: 1,
            nrscat: 1,
            ndcycle: 0,
            ncstsol: 0,
            iyrbp: -50,
            rcl1: [0.15, 0.30, 0.60],
            rcl2: [0.15, 0.30, 0.60],
            acl2: [0.05, 0.10, 0.20],

            gmu0: vec![],
            gmu1: vec![],
            dqo3: vec![],
            dqco2: vec![],
            dtdtlwr: vec![],
            dtdtswr: vec![],
            dqo3cl: vec![],
            gdist2: 1.0,

            time4rad: 0.0,
            time4swr: 0.0,
            time4lwr: 0.0,
            ORB_UNDEF_INT: 2000000000,
            obliqr: 0.0,
            lambm0: 0.0,
            mvelpp: 0.0,
            eccf: 0.0,
            iyrad: 0,
            log_print: true,

            dftde1: vec![],
            dftde2: vec![],
            dftue1: vec![],
            dftue2: vec![],
            dftu0: vec![],
            dft: vec![],
            solclatcdec: 0.0,
            solslat: 0.0,
            solsdec: 0.0,
            solslatsdec: 0.0,
            zmuz: 0.0,
        }
    }
}

impl Radiation {
    pub fn initialize(
        &mut self,
        ndheat: Int,
        neqsig: Int,
        nfixorb: Int,
        n_start_year: Int,
        // planet vars
        eccen: Float,
        obliq: Float,
        mvelp: Float,
    ) {
        // select radiation setup based on resolution
        let mut jtune = 0;
        if ndheat > 0 {
            if NTRU == 21 || NTRU == 1 {
                if NLEV == 5 {
                    if self.ndcycle == 1 {
                        jtune = 0;
                    } else if neqsig == 1 {
                        jtune = 0;
                    } else {
                        self.tswr1 = 0.02;
                        self.tswr2 = 0.065;
                        self.tswr3 = 0.004;
                        self.th2oc = 0.024;
                        jtune = 1;
                    }
                }
                //
                else if NLEV == 10 {
                    if self.ndcycle == 1 {
                        jtune = 0;
                    } else if neqsig == 1 {
                        jtune = 0;
                    } else {
                        self.th2oc = 0.024;
                        self.tswr1 = 0.077;
                        self.tswr2 = 0.065;
                        self.tswr3 = 0.0055;
                        jtune = 1;
                    }
                }
            }
            //
            else if NTRU == 31 {
                if NLEV == 10 {
                    if self.ndcycle == 1 {
                        jtune = 0;
                    } else if neqsig == 1 {
                        jtune = 0;
                    } else {
                        self.tswr1 = 0.077;
                        self.tswr2 = 0.067;
                        self.tswr3 = 0.0055;
                        self.th2oc = 0.024;
                        jtune = 1;
                    }
                }
            }
            //
            else if NTRU == 42 {
                if NLEV == 1 {
                    jtune = 0;
                } else if neqsig == 1 {
                    jtune = 0;
                } else {
                    self.tswr1 = 0.089;
                    self.tswr2 = 0.06;
                    self.tswr3 = 0.0048;
                    self.th2oc = 0.0285;
                    jtune = 1;
                }
            }
        }

        if jtune == 0 {
            println!(
                "No radiation setup for this resolution: NTRU: {}, NLEV: {}",
                NTRU, NLEV
            );
            println!("Using default setup. You may need to tune the radiation")
        }

        // TODO read config

        //
        // determine orbital parameters
        //

        self.iyrbp = 1950 - n_start_year;
        self.iyrad = 1950 - self.iyrbp;

        if nfixorb == 1 {
            // fixed orbital params (default AMIP II)
            self.iyrad = self.ORB_UNDEF_INT;
        }

        orb_params(
            self.iyrad,
            eccen,
            obliq,
            mvelp,
            self.log_print,
            &mut self.obliqr,
            &mut self.lambm0,
            &mut self.mvelpp,
        )

        // if (nfixorb == 1) then ! fixed orbital params (default AMIP II)
        // iyrad = ORB_UNDEF_INT
        // endif
        // call orb_params(iyrad, eccen, obliq, mvelp                          &
        //     &               ,obliqr, lambm0, mvelpp, log_print, mypid, nroot,nud)
    }
}

/// Calculate short wave radiation fluxes
///
/// this parameterization of sw-radiation bases on transmissivities
/// from Lacis & Hansen (1974) for clear sky (H2O,O3,Rayleigh)
/// and Stephens (1978) + Stephens et al. (1984) for clouds.
/// for the verical integration, the adding method is used.
/// some aspects of the realisation are taken from
/// 'a simple radiation parameterization for use in mesoscale models'
/// by S. Bakan (Max-Planck Institut fuer Meteorologie)
/// (unfortunately neither published nor finished)
///
/// no PUMA *subs* are used
///
/// the following PUMA variables are used/modified:
///
/// ga               : gravity acceleration (m/s2) (used)
/// sigma(NLEV)      : full level sigma  (used)
/// sigmah(NLEV)     : half level sigma  (used)
/// dsigma(NLEV)     : delta sigma (half level)  (used)
/// dp(NHOR)         : surface pressure (Pa) (used)
/// dalb(NHOR)       : surface albedo (used)
/// dq(NHOR,NLEP)    : specific humidity (kg/kg) (used)
/// dql(NHOR,NLEP)   : cloud liquid water content (kg/kg) (used)
/// dcc(NHOR,NLEP)   : cloud cover (frac.) (used)
/// dswfl(NHOR,NLEP) : short wave radiation (W/m2)  (modified)
/// dfu(NHOR,NLEP)   : short wave radiation upward (W/m2) (modified)
/// dfd(NHOR,NLEP)   : short wave radiation downward (W/m2) (modified)
fn swr(mut sim: &mut Sim) {
    // 0) define local parameters and arrays

    const zero: Float = 0.000001; // if insolation < zero : fluxes=0.
    const zsolar1: Float = 0.517; // spectral partitioning 1 (wl < 0.75mue)
    const zsolar2: Float = 0.483; // spectral partitioning 2 (wl > 0.75mue)
    const zbetta: Float = 1.66; // magnification factor water vapour
    const zmbar: Float = 1.9; // magnification factor ozon
    const zro3: Float = 2.14; // ozon density (kg/m**3 STP)
    const zfo3: Float = 100. / zro3; // transfere o3 to cm STP

    // transmissivities 1-l
    let zt1: Vec<Vec<Float>>; // (NHOR,NLEP);
    let zt2: Vec<Vec<Float>>; // (NHOR,NLEP);
                              // reflexivities l-1 (scattered)
    let zr1s: Vec<Vec<Float>>; // (NHOR,NLEP);
    let zr2s: Vec<Vec<Float>>; // (NHOR,NLEP);
                               // reflexivities l-NL (direct)
    let zrl1: Vec<Vec<Float>>; // (NHOR,NLEP);
    let zrl2: Vec<Vec<Float>>; // (NHOR,NLEP);
                               // reflexivities l-NL (scattered)
    let zrl1s: Vec<Vec<Float>>; // (NHOR,NLEP);
    let zrl2s: Vec<Vec<Float>>; // (NHOR,NLEP);

    // layer transmissivity (down)
    let ztb1: Vec<Vec<Float>>; // (NHOR,NLEV);
    let ztb2: Vec<Vec<Float>>; // (NHOR,NLEV);
                               // layer transmissivity (up)
    let ztb1u: Vec<Vec<Float>>; // (NHOR,NLEV);
    let ztb2u: Vec<Vec<Float>>; // (NHOR,NLEV);
                                //  layer reflexivity (direct)
    let zrb1: Vec<Vec<Float>>; // (NHOR,NLEV);
    let zrb2: Vec<Vec<Float>>; // (NHOR,NLEV);
                               //  layer reflexibity (scattered)
    let zrb1s: Vec<Vec<Float>>; // (NHOR,NLEV);
    let zrb2s: Vec<Vec<Float>>; // (NHOR,NLEV);

    let zo3l: Vec<Vec<Float>>; // (NHOR,NLEV); // ozon amount (top-l)
    let zxo3l: Vec<Vec<Float>>; // (NHOR,NLEV); // effective ozon amount (top-l)
    let zwvl: Vec<Vec<Float>>; // (NHOR,NLEV); // water vapor amount (top-l)
    let zywvl: Vec<Vec<Float>>; // (NHOR,NLEV); // effective water vapor amount (top-l)
    let zrcs: Vec<Vec<Float>>; // (NHOR,NLEV); // clear sky reflexivity (downward beam)
    let zrcsu: Vec<Vec<Float>>; // (NHOR,NLEV); // clear sky reflexivity (upward beam)

    // top solar radiation
    let zftop1: Vec<Float>; // (NHOR);
    let zftop2: Vec<Float>; // (NHOR);
                            // upward fluxes
    let zfu1: Vec<Float>; // (NHOR);
    let zfu2: Vec<Float>; // (NHOR);
                          // downward fluxes
    let zfd1: Vec<Float>; // (NHOR);
    let zfd2: Vec<Float>; // (NHOR);

    let mut zmu0: Vec<Float> = vec![0.; NHOR]; // (NHOR); zenit angle
    let mut zmu1: Vec<Float> = vec![0.; NHOR]; // (NHOR); zenit angle
    let zcs: Vec<Float>; // (NHOR); clear sky part
    let zm: Vec<Float>; // (NHOR); magnification factor
    let zo3: Vec<Float>; // (NHOR); ozon amount
    let zo3t: Vec<Float>; // (NHOR); total ozon amount (top-sfc)
    let zxo3t: Vec<Float>; // (NHOR); effective total ozon amount (top-sfc)

    // ozon transmissivity (downward/upward beam)
    let zto3: Vec<Float>; // (NHOR); zenit angle
    let zto3u: Vec<Float>; // (NHOR); zenit angle
                           // total ozon transmissivities (d/u)
    let zto3t: Vec<Float>; // (NHOR);
    let zto3tu: Vec<Float>; // (NHOR);

    let zwv: Vec<Float>; // (NHOR); water vapor amount
    let zwvt: Vec<Float>; // (NHOR); total water vapor amount (top-sfc)
    let zywvt: Vec<Float>; // (NHOR); total effective water vapor amount (top-sfc)

    // water vapor trasmissivity (d/u)
    let ztwv: Vec<Float>; // (NHOR);
    let ztwvu: Vec<Float>; // (NHOR);
                           // total water vapor transmissivities (d/u)
    let ztwvt: Vec<Float>; // (NHOR);
    let ztwvtu: Vec<Float>; //

    // reflexivities combined layer (direct)
    let zra1: Vec<Float>; // (NHOR);
    let zra2: Vec<Float>; // (NHOR);
                          // reflexivities combined layer (scatterd)
    let zra1s: Vec<Float>; // (NHOR);
    let zra2s: Vec<Float>; // (NHOR);
                           // transmissivities combined layer (di)
    let zta1: Vec<Float>; // (NHOR);
    let zta2: Vec<Float>; // (NHOR);
                          // transmissivities combined layer (sc)
    let zta1s: Vec<Float>; // (NHOR);
    let zta2s: Vec<Float>; // (NHOR);

    let z1mrabr: Vec<Float>; // (NHOR); 1/(1.-rb*ra(*))

    // cloud reflexivities (direct)
    let mut zrcl1: Vec<Vec<Float>> = vec![vec![0.; NLEV]; NHOR]; // (NHOR,NLEV);
    let zrcl2: Vec<Vec<Float>>; // (NHOR,NLEV);
                                // cloud reflexivities (scattered)
    let mut zrcl1s: Vec<Vec<Float>> = vec![vec![0.; NLEV]; NHOR]; // (NHOR,NLEV);
    let zrcl2s: Vec<Vec<Float>>; // (NHOR,NLEV);
                                 // cloud transmissivities
    let ztcl2: Vec<Vec<Float>>; // (NHOR,NLEV);
    let ztcl2s: Vec<Vec<Float>>; // (NHOR,NLEV);

    // arrays for diagnostic cloud properties
    let zlwp: Vec<Float>; // (NHOR);
    let ztau: Vec<Float>; // (NHOR);
    let zlog: Vec<Float>; // (NHOR);
    let zb2: Vec<Float>; // (NHOR);
    let zom0: Vec<Float>; // (NHOR);
    let zuz: Vec<Float>; // (NHOR);
    let zun: Vec<Float>; // (NHOR);
    let zr: Vec<Float>; // (NHOR);
    let zexp: Vec<Float>; // (NHOR);
    let zu: Vec<Float>; // (NHOR);
    let zb1: Vec<Float>; // (NHOR);

    let losun: Vec<bool>; // (NHOR) flag for gridpoints with insolation

    // cosine of zenith angle

    if sim.rad.ndcycle == 0 {
        let mut js = 1;
        let mut je = NLON;
        for jlat in 1..NLPP {
            // TODO is this correctly translated
            // icnt = count(gmu0(js:je) > 0.0)
            let icnt = sim.rad.gmu0[js..je].iter().filter(|n| **n > 0.).count();
            if icnt > 0 {
                let zsum: Float = sim.rad.gmu0[js..je].iter().sum();
                zmu0[js..je].iter_mut().map(|n| *n = zsum / icnt as Float); // used for clouds
                zmu1[js..je].iter_mut().map(|n| *n = zsum / NLON as Float); // used for insolation
            } else {
                zmu0[js..je].iter_mut().map(|n| *n = 0.0);
                zmu1[js..je].iter_mut().map(|n| *n = 0.0);
            }
            js = js + NLON;
            je = je + NLON;
        } // jlat
    } else {
        zmu0 = sim.rad.gmu0.clone();
        zmu1 = sim.rad.gmu0.clone();
    } // (ndcycle == 0)

    // top solar radiation downward
    // zftop1(:) = zsolar1 * gsol0 * gdist2 * zmu1(:)
    // zftop2(:) = zsolar2 * gsol0 * gdist2 * zmu1(:)
    zftop1 = zmu1
        .iter_mut()
        .map(|mut n| *n * (zsolar1 * sim.rad.gsol0 * sim.rad.gdist2))
        .collect();
    zftop2 = zmu1
        .iter_mut()
        .map(|mut n| *n * (zsolar2 * sim.rad.gsol0 * sim.rad.gdist2))
        .collect();

    // from this point on, all computations are made only for
    // points with solar insolation > zero

    // losun(:) = (zftop1(:) + zftop2(:) > zero)
    losun = vec![false; zftop1.len()]
        .iter_mut()
        .enumerate()
        .map(|(i, n)| (zftop1[i] + zftop2[i]) > 0.)
        .collect();

    // cloud properites

    zcs = vec![1.0; NHOR]; // Clear sky fraction (1.0 = clear sky)
    let zmu00: Float = 0.5;
    let zb3 = sim.rad.tswr1 * zmu00.sqrt() / zmu00;
    let zb4 = sim.rad.tswr2 * zmu00.sqrt();
    let zb5 = sim.rad.tswr3 * zmu00 * zmu00;

    // prescribed

    if sim.rad.nswrcl == 0 {
        for jlev in 1..NLEV {
            if sim.lev_arrays.sigma[jlev] <= 1. / 3. {
                zrcl1s
                    .iter_mut()
                    .map(|n| n[jlev] = sim.rad.rcl1[1] / (sim.rad.rcl1[1] + zmu00));
                zrcl1.iter_mut().enumerate().map(|(i, n)| {
                    n[jlev] = zcs[i] * sim.rad.rcl1[1] / (sim.rad.rcl1[1] + zmu0[i])
                        + (1. - zcs[i]) * zrcl1s[i][jlev]
                });
            }
        }
    }

    //       if (nswrcl == 0) then
    //        do jlev=1,NLEV
    //         if(sigma(jlev) <= 1./3.) then
    //          zrcl1s(:,jlev)=rcl1(1)/(rcl1(1)+zmu00)
    //          zrcl1(:,jlev)=zcs(:)*rcl1(1)/(rcl1(1)+zmu0(:))                 &
    //      &                +(1.-zcs(:))*zrcl1s(:,jlev)
    //          zrcl2s(:,jlev)=AMIN1(1.-acl2(1),rcl2(1)/(rcl2(1)+zmu00))
    //          zrcl2(:,jlev)=AMIN1(1.-acl2(1),zcs(:)*rcl2(1)/(rcl2(1)+zmu0(:))&
    //      &                                +(1.-zcs(:))*zrcl2s(:,jlev))
    //          ztcl2s(:,jlev)=1.-zrcl2s(:,jlev)-acl2(1)
    //          ztcl2(:,jlev)=1.-zrcl2(:,jlev)-acl2(1)
    //         elseif(sigma(jlev) > 1./3. .and. sigma(jlev) <= 2./3.) then
    //          zrcl1s(:,jlev)=rcl1(2)/(rcl1(2)+zmu00)
    //          zrcl1(:,jlev)=zcs(:)*rcl1(2)/(rcl1(2)+zmu0(:))                 &
    //      &                +(1.-zcs(:))*zrcl1s(:,jlev)
    //          zrcl2s(:,jlev)=AMIN1(1.-acl2(2),rcl2(2)/(rcl2(2)+zmu00))
    //          zrcl2(:,jlev)=AMIN1(1.-acl2(2),zcs(:)*rcl2(2)/(rcl2(2)+zmu0(:))&
    //      &                                 +(1.-zcs(:))*zrcl2s(:,jlev))
    //          ztcl2s(:,jlev)=1.-zrcl2s(:,jlev)-acl2(2)
    //          ztcl2(:,jlev)=1.-zrcl2(:,jlev)-acl2(2)
    //         else
    //          zrcl1s(:,jlev)=rcl1(3)/(rcl1(3)+zmu00)
    //          zrcl1(:,jlev)=zcs(:)*rcl1(3)/(rcl1(3)+zmu0(:))                 &
    //      &                +(1.-zcs(:))*zrcl1s(:,jlev)
    //          zrcl2s(:,jlev)=AMIN1(1.-acl2(3),rcl2(3)/(rcl2(3)+zmu00))
    //          zrcl2(:,jlev)=AMIN1(1.-acl2(3),zcs(:)*rcl2(3)/(rcl2(3)+zmu0(:))&
    //      &                                 +(1.-zcs(:))*zrcl2s(:,jlev))
    //          ztcl2s(:,jlev)=1.-zrcl2s(:,jlev)-acl2(3)
    //          ztcl2(:,jlev)=1.-zrcl2(:,jlev)-acl2(3)
    //         endif
    //         zcs(:)=zcs(:)*(1.-dcc(:,jlev))
    //        enddo
    //       else
    //        zrcl1(:,:)=0.0
    //        zrcl2(:,:)=0.0
    //        ztcl2(:,:)=1.0
    //        zrcl1s(:,:)=0.0
    //        zrcl2s(:,:)=0.0
    //        ztcl2s(:,:)=1.0
    //        do jlev=1,NLEV
    //         where(losun(:) .and. (dcc(:,jlev) > 0.))
    //          zlwp(:) = min(1000.0,1000.*dql(:,jlev)*dp(:)/ga*dsigma(jlev))
    //          ztau(:) = 2.0 * ALOG10(zlwp(:)+1.5)**3.9
    //          zlog(:) = log(1000.0 / ztau(:))
    //          zb2(:)  = zb4 / ALOG(3.+0.1*ztau(:))
    //          zom0(:) = min(0.9999,1.0 - zb5 * zlog(:))
    //          zun(:)  = 1.0 - zom0(:)
    //          zuz(:)  = zun(:) + 2.0 * zb2(:) * zom0(:)
    //          zu(:)   = SQRT(zuz(:)/zun(:))
    //          zexp(:) = exp(min(25.0,ztau(:)*SQRT(zuz(:)*zun(:))/zmu00))
    //          zr(:)   = (zu(:)+1.)*(zu(:)+1.)*zexp(:)                      &
    //      &           - (zu(:)-1.)*(zu(:)-1.)/zexp(:)
    //          zrcl1s(:,jlev)=1.-1./(1.+zb3*ztau(:))
    //          ztcl2s(:,jlev)=4.*zu(:)/zr(:)
    //          zrcl2s(:,jlev)=(zu(:)*zu(:)-1.)/zr(:)*(zexp(:)-1./zexp(:))
    //
    //          zb1(:)  = tswr1*SQRT(zmu0(:))
    //          zb2(:)  = tswr2*SQRT(zmu0(:))/ALOG(3.+0.1*ztau(:))
    //          zom0(:) = min(0.9999,1.-tswr3*zmu0(:)*zmu0(:)*zlog(:))
    //          zun(:)  = 1.0 - zom0(:)
    //          zuz(:)  = zun(:) + 2.0 * zb2(:) * zom0(:)
    //          zu(:)   = SQRT(zuz(:)/zun(:))
    //          zexp(:) = exp(min(25.0,ztau(:)*SQRT(zuz(:)*zun(:))/zmu0(:)))
    //          zr(:)   = (zu(:)+1.)*(zu(:)+1.)*zexp(:)                      &
    //      &           - (zu(:)-1.)*(zu(:)-1.)/zexp(:)
    //          zrcl1(:,jlev)=1.-1./(1.+zb1(:)*ztau(:)/zmu0(:))
    //          ztcl2(:,jlev)=4.*zu(:)/zr(:)
    //          zrcl2(:,jlev)=(zu(:)*zu(:)-1.)/zr(:)*(zexp(:)-1./zexp(:))
    //          zrcl1(:,jlev)=zcs(:)*zrcl1(:,jlev)+(1.-zcs(:))*zrcl1s(:,jlev)
    //          ztcl2(:,jlev)=zcs(:)*ztcl2(:,jlev)+(1.-zcs(:))*ztcl2s(:,jlev)
    //          zrcl2(:,jlev)=zcs(:)*zrcl2(:,jlev)+(1.-zcs(:))*zrcl2s(:,jlev)
    //         endwhere
    //         zcs(:)=zcs(:)*(1.-dcc(:,jlev))
    //        enddo ! jlev
    //       endif ! (nswrcl == 0)
    //
    //     // magnification factor
    //
    //       where(losun(:))
    //        zm(:)=35./SQRT(1.+1224.*zmu0(:)*zmu0(:))
    //
    //     // absorber amount and clear sky fraction
    //
    //        zcs(:)=1.
    //        zo3t(:)=0.
    //        zxo3t(:)=0.
    //        zwvt(:)=0.
    //        zywvt(:)=0.
    //       endwhere
    //       do jlev=1,NLEV
    //        where(losun(:))
    //         zo3(:)=zfo3*dsigma(jlev)*dp(:)*dqo3(:,jlev)/ga
    //         zo3t(:)=zo3t(:)+zo3(:)
    //         zxo3t(:)=zcs(:)*(zxo3t(:)+zm(:)*zo3(:))                         &
    //      &          +(1.-zcs(:))*(zxo3t(:)+zmbar*zo3(:))
    //         zo3l(:,jlev)=zo3t(:)
    //         zxo3l(:,jlev)=zxo3t(:)
    //         zwv(:)=0.1*dsigma(jlev)*dq(:,jlev)*dp(:)/ga                     &
    //      &        *SQRT(273./dt(:,jlev))*sigma(jlev)*dp(:)/100000.
    //         zwvt(:)=zwvt(:)+zwv(:)
    //         zywvt(:)=zcs(:)*(zywvt(:)+zm(:)*zwv(:))                         &
    //      &          +(1.-zcs(:))*(zywvt(:)+zbetta*zwv(:))
    //         zwvl(:,jlev)=zwvt(:)
    //         zywvl(:,jlev)=zywvt(:)
    //         zcs(:)=zcs(:)*(1.-dcc(:,jlev))
    //         zrcs(:,jlev)=0.
    //         zrcsu(:,jlev)=0.
    //        endwhere
    //       end do
    //
    //     // compute optical properties
    //     //
    //     // downward loop
    //     //
    //     // preset
    //
    //       where(losun(:))
    //        zta1(:)=1.
    //        zta1s(:)=1.
    //        zra1(:)=0.
    //        zra1s(:)=0.
    //        zta2(:)=1.
    //        zta2s(:)=1.
    //        zra2(:)=0.
    //        zra2s(:)=0.
    // !
    //        zto3t(:)=1.
    //        zo3(:)=zxo3t(:)+zmbar*zo3t(:)
    //        zto3tu(:)=1.                                                     &
    //      &          -(0.02118*zo3(:)/(1.+0.042*zo3(:)+0.000323*zo3(:)**2)   &
    //      &           +1.082*zo3(:)/((1.+138.6*zo3(:))**0.805)               &
    //      &           +0.0658*zo3(:)/(1.+(103.6*zo3(:))**3))/zsolar1
    //        ztwvt(:)=1.
    //        zwv(:)=zywvt(:)+zbetta*zwvt(:)
    //        ztwvtu(:)=1.-2.9*zwv(:)/((1.+141.5*zwv(:))**0.635+5.925*zwv(:))  &
    //      &            /zsolar2
    // !
    // !     clear sky scattering (Rayleigh scatterin lower most level only)
    // !
    //        zrcs(:,NLEV)=(0.219/(1.+0.816*zmu0(:))*zcs(:)                    &
    //      &              +0.144*(1.-zcs(:)-dcc(:,NLEV)))*nrscat
    //        zrcsu(:,NLEV)=0.144*(1.-dcc(:,NLEV))*nrscat
    //       endwhere
    // !
    //       do jlev=1,NLEV
    //        where(losun(:))
    //         zt1(:,jlev)=zta1(:)
    //         zt2(:,jlev)=zta2(:)
    //         zr1s(:,jlev)=zra1s(:)
    //         zr2s(:,jlev)=zra2s(:)
    // !
    // !     set single layer R and T:
    // !
    // !     1. spectral range 1:
    // !
    // !     a) R
    // !     clear part: rayleigh scattering (only lowermost level)
    // !     cloudy part: cloud albedo
    // !
    //         zrb1(:,jlev)=zrcs(:,jlev)+zrcl1(:,jlev)*dcc(:,jlev)
    //         zrb1s(:,jlev)=zrcsu(:,jlev)+zrcl1s(:,jlev)*dcc(:,jlev)
    // !
    // !     b) T
    // !
    // !     ozon absorption
    // !
    // !     downward beam
    // !
    //         zo3(:)=zxo3l(:,jlev)
    //         zto3(:)=(1.                                                     &
    //      &          -(0.02118*zo3(:)/(1.+0.042*zo3(:)+0.000323*zo3(:)**2)   &
    //      &           +1.082*zo3(:)/((1.+138.6*zo3(:))**0.805)               &
    //      &           +0.0658*zo3(:)/(1.+(103.6*zo3(:))**3))/zsolar1)        &
    //      &         /zto3t(:)
    //         zto3t(:)=zto3t(:)*zto3(:)
    // !
    // !     upward scattered beam
    // !
    //         zo3(:)=zxo3t(:)+zmbar*(zo3t(:)-zo3l(:,jlev))
    //         zto3u(:)=zto3tu(:)                                              &
    //      &         /(1.-(0.02118*zo3(:)/(1.+0.042*zo3(:)+0.000323*zo3(:)**2)&
    //      &              +1.082*zo3(:)/((1.+138.6*zo3(:))**0.805)            &
    //      &              +0.0658*zo3(:)/(1.+(103.6*zo3(:))**3))/zsolar1)
    //         zto3tu(:)=zto3tu(:)/zto3u(:)
    // !
    // !     total T = 1-(A(ozon)+R(rayl.))*(1-dcc)-R(cloud)*dcc
    // !
    //         ztb1(:,jlev)=1.-(1.-zto3(:))*(1.-dcc(:,jlev))-zrb1(:,jlev)
    //         ztb1u(:,jlev)=1.-(1.-zto3u(:))*(1.-dcc(:,jlev))-zrb1s(:,jlev)
    // !
    // !     make combined layer R_ab, R_abs, T_ab and T_abs
    // !
    //         z1mrabr(:)=1./(1.-zra1s(:)*zrb1s(:,jlev))
    //         zra1(:)=zra1(:)+zta1(:)*zrb1(:,jlev)*zta1s(:)*z1mrabr(:)
    //         zta1(:)=zta1(:)*ztb1(:,jlev)*z1mrabr(:)
    //         zra1s(:)=zrb1s(:,jlev)+ztb1u(:,jlev)*zra1s(:)*ztb1(:,jlev)      &
    //      &                        *z1mrabr(:)
    //         zta1s(:)=ztb1u(:,jlev)*zta1s(:)*z1mrabr(:)
    // !
    // !     2. spectral range 2:
    // !
    // !     a) R
    // !
    // !     cloud albedo
    // !
    //         zrb2(:,jlev)=zrcl2(:,jlev)*dcc(:,jlev)
    //         zrb2s(:,jlev)=zrcl2s(:,jlev)*dcc(:,jlev)
    // !
    // !     b) T
    // !
    // !     water vapor absorption
    // !
    // !     downward beam
    // !
    //        zwv(:)=zywvl(:,jlev)
    //        ztwv(:)=(1.-2.9*zwv(:)/((1.+141.5*zwv(:))**0.635+5.925*zwv(:))   &
    //      &            /zsolar2)                                             &
    //      &        /ztwvt(:)
    //        ztwvt(:)=ztwvt(:)*ztwv(:)
    // !
    // !     upward scattered beam
    // !
    //        zwv(:)=zywvt(:)+zbetta*(zwvt(:)-zwvl(:,jlev))
    //        ztwvu(:)=ztwvtu(:)                                               &
    //      &         /(1.-2.9*zwv(:)/((1.+141.5*zwv(:))**0.635+5.925*zwv(:))  &
    //      &            /zsolar2)
    //        ztwvtu(:)=ztwvtu(:)/ztwvu(:)
    // !
    // !     total T = 1-A(water vapor)*(1.-dcc)-(A(cloud)+R(cloud))*dcc
    // !
    //         ztb2(:,jlev)=1.-(1.-ztwv(:))*(1.-dcc(:,jlev))                   &
    //      &              -(1.-ztcl2(:,jlev))*dcc(:,jlev)
    //         ztb2u(:,jlev)=1.-(1.-ztwvu(:))*(1.-dcc(:,jlev))                 &
    //      &               -(1.-ztcl2s(:,jlev))*dcc(:,jlev)
    // !
    // !     make combined layer R_ab, R_abs, T_ab and T_abs
    // !
    //         z1mrabr(:)=1./(1.-zra2s(:)*zrb2s(:,jlev))
    //         zra2(:)=zra2(:)+zta2(:)*zrb2(:,jlev)*zta2s(:)*z1mrabr(:)
    //         zta2(:)=zta2(:)*ztb2(:,jlev)*z1mrabr(:)
    //         zra2s(:)=zrb2s(:,jlev)+ztb2u(:,jlev)*zra2s(:)*ztb2(:,jlev)      &
    //      &                       *z1mrabr(:)
    //         zta2s(:)=ztb2u(:,jlev)*zta2s(:)*z1mrabr(:)
    //        endwhere
    //       enddo
    //       where(losun(:))
    //        zt1(:,NLEP)=zta1(:)
    //        zt2(:,NLEP)=zta2(:)
    //        zr1s(:,NLEP)=zra1s(:)
    //        zr2s(:,NLEP)=zra2s(:)
    // !
    // !     upward loop
    // !
    // !     make upward R
    // !
    //        zra1s(:)=dalb(:)
    //        zra2s(:)=dalb(:)
    // !
    // !      set albedo for the direct beam (for ocean use ECHAM3 param)
    // !
    //        dalb(:)=dls(:)*dalb(:)+(1.-dls(:))*dicec(:)*dalb(:)              &
    //      &        +(1.-dls(:))*(1.-dicec(:))*AMIN1(0.05/(zmu0(:)+0.15),0.15)
    //        zra1(:)=dalb(:)
    //        zra2(:)=dalb(:)
    //       endwhere
    //       do jlev=NLEV,1,-1
    //        where(losun(:))
    //         zrl1(:,jlev+1)=zra1(:)
    //         zrl2(:,jlev+1)=zra2(:)
    //         zrl1s(:,jlev+1)=zra1s(:)
    //         zrl2s(:,jlev+1)=zra2s(:)
    //         zra1(:)=zrb1(:,jlev)+ztb1(:,jlev)*zra1(:)*ztb1u(:,jlev)         &
    //      &                      /(1.-zra1s(:)*zrb1s(:,jlev))
    //         zra1s(:)=zrb1s(:,jlev)+ztb1u(:,jlev)*zra1s(:)*ztb1u(:,jlev)     &
    //      &                        /(1.-zra1s(:)*zrb1s(:,jlev))
    //         zra2(:)=zrb2(:,jlev)+ztb2(:,jlev)*zra2(:)*ztb2u(:,jlev)         &
    //      &                      /(1.-zra2s(:)*zrb2s(:,jlev))
    //         zra2s(:)=zrb2s(:,jlev)+ztb2u(:,jlev)*zra2s(:)*ztb2u(:,jlev)     &
    //      &                        /(1.-zra2s(:)*zrb2s(:,jlev))
    //        endwhere
    //       enddo
    //       where(losun(:))
    //        zrl1(:,1)=zra1(:)
    //        zrl2(:,1)=zra2(:)
    //        zrl1s(:,1)=zra1s(:)
    //        zrl2s(:,1)=zra2s(:)
    //       endwhere
    // !
    // !     fluxes at layer interfaces
    // !
    //       do jlev=1,NLEP
    //        where(losun(:))
    //         z1mrabr(:)=1./(1.-zr1s(:,jlev)*zrl1s(:,jlev))
    //         zfd1(:)=zt1(:,jlev)*z1mrabr(:)
    //         zfu1(:)=-zt1(:,jlev)*zrl1(:,jlev)*z1mrabr(:)
    //         z1mrabr(:)=1./(1.-zr2s(:,jlev)*zrl2s(:,jlev))
    //         zfd2(:)=zt2(:,jlev)*z1mrabr(:)
    //         zfu2(:)=-zt2(:,jlev)*zrl2(:,jlev)*z1mrabr(:)
    //         dfu(:,jlev)=zfu1(:)*zftop1(:)+zfu2(:)*zftop2(:)
    //         dfd(:,jlev)=zfd1(:)*zftop1(:)+zfd2(:)*zftop2(:)
    //         dswfl(:,jlev)=dfu(:,jlev)+dfd(:,jlev)
    //        endwhere
    //       enddo
}

/// Compute long wave radiation
///
/// o3-, co2-, h2o- and cloud-absorption is considered
///
/// clear sky absorptivities from Sasamori 1968
/// (J. Applied Meteorology, 7, 721-729)
///
/// no PUMA *subs* are used
///
/// the following PUMA variables are used/modified:
///
/// ga               : gravity accelleration (m/s2) (used)
/// sigma(NLEV)      : full level sigma  (used)
/// sigmah(NLEV)     : half level sigma  (used)
/// dsigma(NLEV)     : delta sigma (half level)  (used)
/// dp(NHOR)         : surface pressure (Pa) (used)
/// dq(NHOR,NLEP)    : specific humidity (kg/kg) (used)
/// dt(NHOR,NLEP)    : temperature (K) (used)
/// dcc(NHOR,NLEP)   : cloud cover (frac.) (used)
/// dlwfl(NHOR,NLEP) : long wave radiation (W/m2)  (modified)
/// dftu(NHOR,NLEP)  : long wave radiation upward (W/m2) (modified)
/// dftd(NHOR,NLEP)  : long wave radiation downward (W/m2) (modified)
fn lwr() {
    //
    // 0) define local parameters
    //

    const zmmair: Float = 0.0289644; // molecular weight air (kg/mol)
    const zmmco2: Float = 0.0440098; // molecular weight co2 (kg/mol)
    const zpv2pm: Float = zmmco2 / zmmair; // transfere co2 ppvol to ppmass
    const zrco2: Float = 1.9635; // co2 density (kg/m3 stp)
    const zro3: Float = 2.14; // o3 density (kg/m3 stp)
    const zttop: Float = 0.; // t at top of atmosphere

    // scaling factors for transmissivities
    // uh2o in g/cm**2 (zfh2o=0.1=1000/(100*100) *kg/m**2)
    // uco2 in cm-STP  (zfco2=100./rco2 * kg/m**2)
    // uo3  in cm-STP  (zfo3=100./ro3 * kg/m**2)
    const zfh2o: Float = 0.1;
    const zfco2: Float = 100. / zrco2;
    const zfo3: Float = 100. / zro3;
    const zt0: Float = 295.;

    // local arrays
    // TODO make sure how this syntax works: "array(0:VAR)"
    // real zbu(NHOR,0:NLEP)
    let zbu = vec![vec![0.; NHOR]; NLEP + 1]; // effective SBK*T**4 for upward radiation
    let zbd = vec![vec![0.; NHOR]; NLEP + 1]; // effective SBK*T**4 for downward radiation
    let zst4h = vec![vec![0.; NHOR]; NLEP]; // SBK*T**4  on half levels
    let zst4 = vec![vec![0.; NHOR]; NLEP]; // SBK*T**4  on ull levels
    let ztau = vec![vec![0.; NHOR]; NLEV]; // total transmissivity
    let zq = vec![vec![0.; NHOR]; NLEV]; // modified water vapour
    let zqo3 = vec![vec![0.; NHOR]; NLEV]; // modified ozon
    let zqco2 = vec![vec![0.; NHOR]; NLEV]; // modified co2
    let ztausf = vec![vec![0.; NHOR]; NLEV]; // total transmissivity to surface
    let ztaucs = vec![vec![0.; NHOR]; NLEV]; // clear sky transmissivity
    let ztaucc0 = vec![vec![0.; NHOR]; NLEV]; // layer transmissivity cloud
    let ztaucc = vec![0.; NHOR]; // cloud transmissivity
    let ztau0 = vec![0.; NHOR]; // approx. layer transmissivity
    let zsumwv = vec![0.; NHOR]; // effective water vapor amount
    let zsumo3 = vec![0.; NHOR]; // effective o3 amount
    let zsumco2 = vec![0.; NHOR]; // effective co2 amount
    let zsfac = vec![0; NHOR]; // scaling factor
    let zah2o = vec![0; NHOR]; // water vapor absorptivity
    let zaco2 = vec![0.; NHOR]; // co2 absorptivity
    let zao3 = vec![0.; NHOR]; // o3 absorptivity
    let zth2o = vec![0.; NHOR]; // water vapor - co2 overlap transmissivity
    let zbdl = vec![0.; NHOR]; // layer evective downward rad.
    let zeps = vec![0.; NHOR]; // surface emissivity
    let zps2 = vec![0.; NHOR]; // ps**2
    let mut zsigh2 = vec![0.; NLEP]; // sigmah**2

    // entropy
    let zbue1 = vec![vec![0; NHOR]; NLEP + 1];
    let zbue2 = vec![vec![0; NHOR]; NLEP + 1];

    //
    // 1) set some necessary (helpful) bits
    //

    let zero = 0.000001; // small number

    let zao30 = 0.209 * (0.00007 as Float).powf(0.436); // to get a(o3)=0 for o3=0
    let zco20 = 0.0676 * (0.01022 as Float).powf(0.421); // to get a(co2)=0 for co2=0
    let zh2o0a = 0.846 * (0.0000359 as Float).powf(0.243); // to get a(h2o)=0 for h2o=0
    let zh2o0 = 0.832 * ((0.0286 as Float).powf(0.26)); // to get t(h2o)=1 for h2o=0

    // to make a(o3) continues at 0.01cm:
    let zao3c =
        0.209 * ((0.01 + 0.00007) as Float).powf(0.436) - zao30 - 0.0212 * (0.01 as Float).log10();

    // to make a(co2) continues at 1cm:
    let zaco2c = 0.0676 * (1.01022 as Float).powf(0.421) - zco20;

    // to make a(h2o) continues at 0.01gm:
    let zah2oc =
        0.846 * (((0.01 + 3.59E-5) as Float).powf(0.243)) - zh2o0a - 0.24 * (0.02 as Float).log10();

    // to make t(h2o) continues at 2gm :
    let zth2oc = 1. - (0.832 * (((2. + 0.0286) as Float).powf(0.26)) - zh2o0)
        + 0.1196 * ((2. - 0.6931) as Float).ln();

    zsigh2[1] = 0.;

    // TODO translate this
    //zsigh2(2: NLEP) = sigmah(1: NLEV) * *2;
    // zps2(:)=dp(:)*dp(:)

    // 2) calc. stb*t**4 + preset fluxes

    // TODO
    // dftu(:,:)=0.;
    // dftd(:,:)=0.;

    // stb*t**4 on full and half levels

    // full levels and surface

    // TODO
    // zst4(:,1:NLEP)=SBK*dt(:,1:NLEP)**4;

    // half level (incl. toa and near surface)

    // TODO
    // zst4h(:,1)=zst4(:,1)-(zst4(:,1)-zst4(:,2))*sigma(1)/(sigma(1)-sigma(2));

    for jlev in 2..NLEV {
        let jlem = jlev - 1;
        // TODO
        // zst4h(:,jlev)=(zst4(:,jlev)*(sigma(jlem)-sigmah(jlem))
        //                    +zst4(:,jlem)*(sigmah(jlem)-sigma(jlev)))
        //                   /(sigma(jlem)-sigma(jlev));
    }

    // TODO
    // where((zst4(:,NLEV)-zst4h(:,NLEV))                                &
    //     &     *(zst4(:,NLEV)-zst4(:,NLEP)) > 0.)
    // zst4h(:,NLEP)=zst4(:,NLEP)
    // elsewhere
    // zst4h(:,NLEP)=zst4(:,NLEV)                                       &
    //     &              +(zst4(:,NLEV)-zst4(:,NLEM))*(1.-sigma(NLEV))       &
    //     &              /(sigma(NLEV)-sigma(NLEM))
    // endwhere

    // TODO
    // // top downward flux, surface grayness and surface upward flux
    // zbd(:,0)=SBK*zttop**4
    // zeps(:)=dls(:)+0.98*(1.-dls(:))
    // zbu(:,NLEP)=zeps(:)*zst4(:,NLEP)
    // zbue1(:,NLEP)=0.
    // zbue2(:,NLEP)=zbu(:,NLEP)

    // TODO more
}

/// Calculate earth's orbital parameters using Dave Threshers
/// formula which came from Berger, Andre.  1978
/// "A Simple Algorithm to Compute Long-Term Variations
/// of Daily Insolation".  Contribution 18, Institute of Astronomy and
/// Geophysics, Universite Catholique de Louvain, Louvain-la-Neuve,
/// Belgium.
///
/// Original Author: Erik Kluzek
/// Date:            Oct/97
fn orb_params(
    iyear_AD: Int,
    eccen: Float,
    obliq: Float,
    mvelp: Float,
    log_print: bool,
    mut obliqr: &mut Float,
    mut lambm0: &mut Float,
    mut mvelpp: &mut Float,
) {
    // Parameters for calculating earth's orbital characteristics

    const POBLEN: usize = 47; // number of elements in the series to calc obliquity
    const PECCLEN: usize = 19; // number of elements in the series to calc eccentricity
    const PMVELEN: usize = 78; // number of elements in the series to calc vernal equinox

    let degrad: Float; // degrees to radians conversion factor
    let obamp: [Float; POBLEN]; // amplitudes for obliquity cosine series
    let obrate: [Float; POBLEN]; // rates for obliquity cosine series
    let obphas: [Float; POBLEN]; // phases for obliquity cosine series
    let ecamp: [Float; PECCLEN]; // amplitudes for eccentricity/fvelp cosine/sine series
    let ecrate: [Float; PECCLEN]; // rates for eccentricity/fvelp cosine/sine series
    let ecphas: [Float; PECCLEN]; // phases for eccentricity/fvelp cosine/sine series
    let mvamp: [Float; PMVELEN]; // amplitudes for mvelp sine series
    let mvrate: [Float; PMVELEN]; // rates for mvelp sine series
    let mvphas: [Float; PMVELEN]; // phases for mvelp sine series
    let mut yb4_1950AD: Float = 0.; // number of years before 1950 AD

    let psecdeg: Float = 1. / 3600.; // arc seconds to degrees conversion

    // Cosine series data for computation of obliquity:
    // amplitude (arc seconds), rate (arc seconds/year), phase (degrees).

    obamp = [
        -2462.2214466,
        -857.3232075,
        -629.3231835,
        -414.2804924,
        -311.7632587,
        308.9408604,
        -162.5533601,
        -116.1077911,
        101.1189923,
        -67.6856209,
        24.9079067,
        22.5811241,
        -21.1648355,
        -15.6549876,
        15.3936813,
        14.6660938,
        -11.7273029,
        10.2742696,
        6.4914588,
        5.8539148,
        -5.4872205,
        -5.4290191,
        5.1609570,
        5.0786314,
        -4.0735782,
        3.7227167,
        3.3971932,
        -2.8347004,
        -2.6550721,
        -2.5717867,
        -2.4712188,
        2.4625410,
        2.2464112,
        -2.0755511,
        -1.9713669,
        -1.8813061,
        -1.8468785,
        1.8186742,
        1.7601888,
        -1.5428851,
        1.4738838,
        -1.4593669,
        1.4192259,
        -1.1818980,
        1.1756474,
        -1.1316126,
        1.0896928,
    ];

    obrate = [
        31.609974, 32.620504, 24.172203, 31.983787, 44.828336, 30.973257, 43.668246, 32.246691,
        30.599444, 42.681324, 43.836462, 47.439436, 63.219948, 64.230478, 1.010530, 7.437771,
        55.782177, 0.373813, 13.218362, 62.583231, 63.593761, 76.438310, 45.815258, 8.448301,
        56.792707, 49.747842, 12.058272, 75.278220, 65.241008, 64.604291, 1.647247, 7.811584,
        12.207832, 63.856665, 56.155990, 77.448840, 6.801054, 62.209418, 20.656133, 48.344406,
        55.145460, 69.000539, 11.071350, 74.291298, 11.047742, 0.636717, 12.844549,
    ];

    obphas = [
        251.9025, 280.8325, 128.3057, 292.7252, 15.3747, 263.7951, 308.4258, 240.0099, 222.9725,
        268.7809, 316.7998, 319.6024, 143.8050, 172.7351, 28.9300, 123.5968, 20.2082, 40.8226,
        123.4722, 155.6977, 184.6277, 267.2772, 55.0196, 152.5268, 49.1382, 204.6609, 56.5233,
        200.3284, 201.6651, 213.5577, 17.0374, 164.4194, 94.5422, 131.9124, 61.0309, 296.2073,
        135.4894, 114.8750, 247.0691, 256.6114, 32.1008, 143.6804, 16.8784, 160.6835, 27.5932,
        348.1074, 82.6496,
    ];

    // Cosine/sine series data for computation of eccentricity and
    // fixed vernal equinox longitude of perihelion (fvelp):
    // amplitude, rate (arc seconds/year), phase (degrees).

    ecamp = [
        0.01860798,
        0.01627522,
        -0.01300660,
        0.00988829,
        -0.00336700,
        0.00333077,
        -0.00235400,
        0.00140015,
        0.00100700,
        0.00085700,
        0.00064990,
        0.00059900,
        0.00037800,
        -0.00033700,
        0.00027600,
        0.00018200,
        -0.00017400,
        -0.00012400,
        0.00001250,
    ];

    ecrate = [
        4.2072050, 7.3460910, 17.8572630, 17.2205460, 16.8467330, 5.1990790, 18.2310760,
        26.2167580, 6.3591690, 16.2100160, 3.0651810, 16.5838290, 18.4939800, 6.1909530,
        18.8677930, 17.4255670, 6.1860010, 18.4174410, 0.6678630,
    ];

    ecphas = [
        28.620089, 193.788772, 308.307024, 320.199637, 279.376984, 87.195000, 349.129677,
        128.443387, 154.143880, 291.269597, 114.860583, 332.092251, 296.414411, 145.769910,
        337.237063, 152.092288, 126.839891, 210.667199, 72.108838,
    ];

    // Sine series data for computation of moving vernal equinox longitude of perihelion:
    //  amplitude (arc seconds), rate (arc seconds/year), phase (degrees).

    mvamp = [
        7391.0225890,
        2555.1526947,
        2022.7629188,
        -1973.6517951,
        1240.2321818,
        953.8679112,
        -931.7537108,
        872.3795383,
        606.3544732,
        -496.0274038,
        456.9608039,
        346.9462320,
        -305.8412902,
        249.6173246,
        -199.1027200,
        191.0560889,
        -175.2936572,
        165.9068833,
        161.1285917,
        139.7878093,
        -133.5228399,
        117.0673811,
        104.6907281,
        95.3227476,
        86.7824524,
        86.0857729,
        70.5893698,
        -69.9719343,
        -62.5817473,
        61.5450059,
        -57.9364011,
        57.1899832,
        -57.0236109,
        -54.2119253,
        53.2834147,
        52.1223575,
        -49.0059908,
        -48.3118757,
        -45.4191685,
        -42.2357920,
        -34.7971099,
        34.4623613,
        -33.8356643,
        33.6689362,
        -31.2521586,
        -30.8798701,
        28.4640769,
        -27.1960802,
        27.0860736,
        -26.3437456,
        24.7253740,
        24.6732126,
        24.4272733,
        24.0127327,
        21.7150294,
        -21.5375347,
        18.1148363,
        -16.9603104,
        -16.1765215,
        15.5567653,
        15.4846529,
        15.2150632,
        14.5047426,
        -14.3873316,
        13.1351419,
        12.8776311,
        11.9867234,
        11.9385578,
        11.7030822,
        11.6018181,
        -11.2617293,
        -10.4664199,
        10.4333970,
        -10.2377466,
        10.1934446,
        -10.1280191,
        10.0289441,
        -10.0034259,
    ];

    mvrate = [
        31.609974, 32.620504, 24.172203, 0.636717, 31.983787, 3.138886, 30.973257, 44.828336,
        0.991874, 0.373813, 43.668246, 32.246691, 30.599444, 2.147012, 10.511172, 42.681324,
        13.650058, 0.986922, 9.874455, 13.013341, 0.262904, 0.004952, 1.142024, 63.219948,
        0.205021, 2.151964, 64.230478, 43.836462, 47.439436, 1.384343, 7.437771, 18.829299,
        9.500642, 0.431696, 1.160090, 55.782177, 12.639528, 1.155138, 0.168216, 1.647247,
        10.884985, 5.610937, 12.658184, 1.010530, 1.983748, 14.023871, 0.560178, 1.273434,
        12.021467, 62.583231, 63.593761, 76.438310, 4.280910, 13.218362, 17.818769, 8.359495,
        56.792707, 8.448301, 1.978796, 8.863925, 0.186365, 8.996212, 6.771027, 45.815258,
        12.002811, 75.278220, 65.241008, 18.870667, 22.009553, 64.604291, 11.498094, 0.578834,
        9.237738, 49.747842, 2.147012, 1.196895, 2.133898, 0.173168,
    ];

    mvphas = [
        251.9025, 280.8325, 128.3057, 348.1074, 292.7252, 165.1686, 263.7951, 15.3747, 58.5749,
        40.8226, 308.4258, 240.0099, 222.9725, 106.5937, 114.5182, 268.7809, 279.6869, 39.6448,
        126.4108, 291.5795, 307.2848, 18.9300, 273.7596, 143.8050, 191.8927, 125.5237, 172.7351,
        316.7998, 319.6024, 69.7526, 123.5968, 217.6432, 85.5882, 156.2147, 66.9489, 20.2082,
        250.7568, 48.0188, 8.3739, 17.0374, 155.3409, 94.1709, 221.1120, 28.9300, 117.1498,
        320.5095, 262.3602, 336.2148, 233.0046, 155.6977, 184.6277, 267.2772, 78.9281, 123.4722,
        188.7132, 180.1364, 49.1382, 152.5268, 98.2198, 97.4808, 221.5376, 168.2438, 161.1199,
        55.0196, 262.6495, 200.3284, 201.6651, 294.6547, 99.8233, 213.5577, 154.1631, 232.7153,
        138.3034, 204.6609, 106.5938, 250.4676, 332.3345, 27.3039,
    ];

    // Local variables

    let i: Int; // Index for series summations
    let mut obsum: Float; // Obliquity series summation
    let mut cossum: Float; // Cosine series summation for eccentricity/fvelp
    let mut sinsum: Float; // Sine series summation for eccentricity/fvelp
    let mut fvelp: Float = 0.0; // Fixed vernal equinox longitude of perihelion
    let mut mvsum: Float; // mvelp series summation
    let beta: Float; // Intermediate argument for lambm0
    let years: Float; // Years to time of interest (negative = past; positive = future)

    let mut eccen2: Float; // eccentricity squared
    let mut eccen3: Float; // eccentricity cubed
    let pi: Float = 4. * (1. as Float).atan() as Float; // pi

    // radinp and algorithms below will need a degrees to radians conversion factor.
    degrad = pi / 180.;

    // Check for flag to use input orbit parameters
    if iyear_AD == ORB_NOT_YEAR_BASED {
        // Check input obliq, eccen, and mvelp to ensure reasonable
        if obliq == ORB_UNDEF_REAL {
            if log_print {
                // TODO print
                // write(nud,*)'(orb_params) Have to specify orbital parameters:'
                // write(nud,*) 'Either set: '                                   &
                //     &                ,'iyear_AD, OR [obliq, eccen, and mvelp]:'
                // write(nud,*)'iyear_AD is the year to simulate the orbit for ' &
                //     &                ,'(ie. 1950): '
                // write(nud,*)'obliq, eccen, mvelp specify the orbit directly:'
                // write(nud,*)'The AMIP II settings (for a 1995 orbit) are: '
                // write(nud,*)' obliq = 23.4441'
                // write(nud,*)' eccen = 0.016715'
                // write(nud,*)' mvelp = 102.7'
            }
            panic!();
        } else if log_print {
            // TODO print
            // write(nud,*)'(orb_params) Use input orbital parameters: '
        }
        if obliq.lt(&ORB_OBLIQ_MIN) || obliq.gt(&ORB_OBLIQ_MAX) {
            if log_print {
                // TODO print
                // write(nud,*) '(orb_params): Input obliquity unreasonable: '  &
                //     &                  ,obliq
            }
            panic!();
        }
        if eccen.lt(&ORB_ECCEN_MIN) || eccen.gt(&ORB_ECCEN_MAX) {
            if log_print {
                // TODO print
                // write(nud,*) '(orb_params): Input eccentricity unreasonable: '&
                //     &                 ,eccen
            }
            panic!();
        }

        if mvelp.lt(&ORB_MVELP_MIN) || mvelp.gt(&ORB_MVELP_MAX) {
            if log_print {
                // TODO print
                // write(nud,*)'(orb_params): Input mvelp unreasonable: ', mvelp
            }
            panic!();
        }
        eccen2 = eccen * eccen;
        eccen3 = eccen2 * eccen;
    }
    // Otherwise calculate based on years before present
    else {
        yb4_1950AD = 1950.0 - iyear_AD as Float;
        if yb4_1950AD.abs().gt(&1000000.0) {
            if log_print {
                // TODO print
                // write(nud,*)'(orb_params) orbit only valid for years+-1000000'
                // write(nud,*)'(orb_params) Relative to 1950 AD'
                // write(nud,*)'(orb_params) # of years before 1950: ',yb4_1950AD
                // write(nud,*)'(orb_params) Year to simulate was  : ',iyear_AD
            }
            panic!();
        }

        // The following calculates the earth's obliquity, orbital eccentricity
        // (and various powers of it) and vernal equinox mean longitude of
        // perihelion for years in the past (future = negative of years past),
        // using constants (see parameter section) given in the program of:
        //
        // Berger, Andre.  1978  A Simple Algorithm to Compute Long-Term Variations
        // of Daily Insolation.  Contribution 18, Institute of Astronomy and
        // Geophysics, Universite Catholique de Louvain, Louvain-la-Neuve, Belgium.
        //
        // and formulas given in the paper (where less precise constants are also
        // given):
        //
        // Berger, Andre.  1978.  Long-Term Variations of Daily Insolation and
        // Quaternary Climatic Changes.  J. of the Atmo. Sci. 35:2362-2367
        //
        // The algorithm is valid only to 1,000,000 years past or hence.
        // For a solution valid to 5-10 million years past see the above author.
        // Algorithm below is better for years closer to present than is the
        // 5-10 million year solution.
        //
        // Years to time of interest must be negative of years before present
        // (1950) in formulas that follow.
        //
        years = -yb4_1950AD;
        //
        // In the summations below, cosine or sine arguments, which end up in
        // degrees, must be converted to radians via multiplication by degrad.
        //
        // Summation of cosine series for obliquity (epsilon in Berger 1978) in
        // degrees. Convert the amplitudes and rates, which are in arc seconds, into
        // degrees via multiplication by psecdeg (arc seconds to degrees conversion
        // factor).  For obliq, first term is Berger 1978's epsilon star; second
        // term is series summation in degrees.

        obsum = 0.0;

        for i in 0..POBLEN {
            obsum = obsum
                + obamp[i] * psecdeg * ((obrate[i] * psecdeg * years + obphas[i]) * degrad).cos();
        }
        let obliq = 23.320556 + obsum;

        // Summation of cosine and sine series for computation of eccentricity
        // (eccen; e in Berger 1978) and fixed vernal equinox longitude of perihelion
        // (fvelp; pi in Berger 1978), which is used for computation of moving vernal
        // equinox longitude of perihelion.  Convert the rates, which are in arc
        // seconds, into degrees via multiplication by psecdeg.

        cossum = 0.0;

        for i in 0..PECCLEN {
            cossum = cossum + ecamp[i] * ((ecrate[i] * psecdeg * years + ecphas[i]) * degrad).cos();
        }

        sinsum = 0.0;

        for i in 0..PECCLEN {
            sinsum = sinsum + ecamp[i] * ((ecrate[i] * psecdeg * years + ecphas[i]) * degrad).sin();
        }

        // Use summations to calculate eccentricity

        eccen2 = cossum * cossum + sinsum * sinsum;
        let eccen = eccen2.sqrt();
        eccen3 = eccen2 * eccen;

        // A series of cases for fvelp, which is in radians.

        if cossum.abs() <= 1.0E-8 {
            if sinsum == 0.0 {
                fvelp = 0.0;
            } else if sinsum < 0.0 {
                fvelp = 1.5 * pi;
            } else if sinsum > 0.0 {
                fvelp = 0.5 * pi;
            }
        } else if cossum < 0.0 {
            fvelp = (sinsum / cossum).atan() + pi;
        } else if cossum > 0.0 {
            if sinsum < 0.0 {
                fvelp = (sinsum / cossum).atan() + 2.0 * pi;
            } else {
                fvelp = (sinsum / cossum).atan();
            }
        }

        // Summation of sine series for computation of moving vernal equinox longitude
        // of perihelion (mvelp; omega bar in Berger 1978) in degrees.  For mvelp,
        // first term is fvelp in degrees; second term is Berger 1978's psi bar times
        // years and in degrees; third term is Berger 1978's zeta; fourth term is
        // series summation in degrees.  Convert the amplitudes and rates, which are
        // in arc seconds, into degrees via multiplication by psecdeg.  Series summation
        // plus second and third terms constitute Berger 1978's psi, which is the
        // general precession.

        mvsum = 0.0;

        for i in 0..PMVELEN {
            mvsum = mvsum
                + mvamp[i] * psecdeg * ((mvrate[i] * psecdeg * years + mvphas[i]) * degrad).sin();
        }
        let mut mvelp = fvelp / degrad + 50.439273 * psecdeg * years + 3.392506 + mvsum;

        // Cases to make sure mvelp is between 0 and 360.

        while mvelp < 0.0 {
            mvelp = mvelp + 360.0;
        }
        while mvelp >= 360.0 {
            mvelp = mvelp - 360.0;
        }
    } // end of test on whether to calculate or use input orbital params

    // Orbit needs the obliquity in radians

    *obliqr = obliq * degrad;

    // 180 degrees must be added to mvelp since observations are made from the
    // earth and the sun is considered (wrongly for the algorithm) to go around
    // the earth. For a more graphic explanation see Appendix B in:
    //
    // A. Berger, M. Loutre and C. Tricot. 1993.  Insolation and Earth's Orbital
    // Periods.  J. of Geophysical Research 98:10,341-10,362.
    //
    // Additionally, orbit will need this value in radians. So mvelp becomes
    // mvelpp (mvelp plus pi)
    //
    *mvelpp = (mvelp + 180.) * degrad;

    // Set up an argument used several times in lambm0 calculation ahead.

    beta = (1. - eccen2).sqrt();

    // The mean longitude at the vernal equinox (lambda m nought in Berger
    // 1978; in radians) is calculated from the following formula given in
    // Berger 1978.  At the vernal equinox the true longitude (lambda in Berger
    // 1978) is 0.

    *lambm0 = 2.
        * ((0.5 * eccen + 0.125 * eccen3) * (1. + beta) * mvelpp.sin()
            - 0.25 * eccen2 * (0.5 + beta) * (2. * *mvelpp).sin()
            + 0.125 * eccen3 * (1. / 3. + beta) * (3. * *mvelpp).sin());

    // if ( log_print ) then
    //  if(mypid==nroot) then
    //   write(nud,'(/," *****************************************")')
    //   write(nud,'(" *     Computed Orbital Parameters       *")')
    //   write(nud,'(" *****************************************")')
    //   write(nud,'(" * Year AD           =  ",i16  ," *")') iyear_AD
    //   write(nud,'(" * Eccentricity      =  ",f16.6," *")') eccen
    //   write(nud,'(" * Obliquity (deg)   =  ",f16.6," *")') obliq
    //   write(nud,'(" * Obliquity (rad)   =  ",f16.6," *")') obliqr
    //   write(nud,'(" * Long of perh(deg) =  ",f16.6," *")') mvelp
    //   write(nud,'(" * Long of perh(rad) =  ",f16.6," *")') mvelpp
    //   write(nud,'(" * Long at v.e.(rad) =  ",f16.6," *")') lambm0
    //   write(nud,'(" *****************************************")')
    //  end if
    // end if
}

/// Compute earth/orbit parameters using formula suggested by
/// Duane Thresher.
///
/// Original version:  Erik Kluzek
/// Date:              Oct/1997
///
/// Modification: 22-Feb-2006 (ek) - get days/yr from pumamod
fn orb_decl(
    calday: Float,
    eccen: Float,
    mvelpp: Float,
    lambm0: Float,
    obliqr: Float,
    n_days_per_year: Int,
    ndatim: Vec<Int>,
    mut delta: &mut Float,
    mut eccf: &mut Float,
) {
    // Calday of vernal equinox, correct for Jan 1 = calday 1
    let ve: Float = 80.5;
    let pie = std::f64::consts::PI as Float;

    // Lambda m, earth's mean longitude of perihelion (radians)
    let mut lambm: Float;
    // Intermediate argument involving lambm
    let mut lmm: Float;
    // Lambda, the earth's longitude of perihelion
    let lamb: Float;
    // Inverse normalized sun/earth distance
    let mut invrho: Float;
    // Sine of lmm
    let sinl: Float;

    // Compute eccentricity factor and solar declination using
    // day value where a round day (such as 213.0) refers to 0z at
    // Greenwich longitude.
    //
    // Use formulas from Berger, Andre 1978: Long-Term Variations of Daily
    // Insolation and Quaternary Climatic Changes. J. of the Atmo. Sci.
    // 35:2362-2367.
    //
    // To get the earth's true longitude (position in orbit; lambda in Berger 1978),
    // which is necessary to find the eccentricity factor and declination,
    // must first calculate the mean longitude (lambda m in Berger 1978) at
    // the present day.  This is done by adding to lambm0 (the mean longitude
    // at the vernal equinox, set as March 21 at noon, when lambda = 0; in radians)
    // an increment (delta lambda m in Berger 1978) that is the number of
    // days past or before (a negative increment) the vernal equinox divided by
    // the days in a model year times the 2*pi radians in a complete orbit.

    lambm = lambm0 + (calday - ve) * 2. * pie / (n_days_per_year + ndatim[6]) as Float; // ndatim(7) = leap year
    lmm = lambm - mvelpp;

    // The earth's true longitude, in radians, is then found from
    // the formula in Berger 1978:

    sinl = lmm.sin();
    lamb = lambm
        + eccen
            * (2. * sinl
                + eccen
                    * (1.25 * (2. * lmm).sin()
                        + eccen * ((13.0 / 12.0) * (3. * lmm).sin() - 0.25 * sinl)));

    // Using the obliquity, eccentricity, moving vernal equinox longitude of
    // perihelion (plus), and earth's true longitude, the declination (delta)
    // and the normalized earth/sun distance (rho in Berger 1978; actually inverse
    // rho will be used), and thus the eccentricity factor (eccf), can be calculated
    // from formulas given in Berger 1978.

    invrho = (1. + eccen * (lamb - mvelpp).cos()) / (1. - eccen * eccen);

    // Set solar declination and eccentricity factor

    *delta = (obliqr.sin() * lamb.sin()).asin();
    *eccf = invrho * invrho;
}

/// Print out the information on the input orbital characteristics
///
/// Original version:  Erik Kluzek
/// Date:              Oct/1997
fn orb_print(iyear_AD: Int, eccen: Float, obliq: Float, mvelp: Float) {
    if iyear_AD == ORB_NOT_YEAR_BASED {
        if obliq == ORB_UNDEF_REAL {
            println!("Orbit parameters not set!");
        } else {
            println!("Orbital parameters:");
            println!("Obliquity (degree): {}", obliq);
            println!("Eccentricity (unitless): {}", eccen);
            println!("Long. of moving Perhelion (deg): {}", mvelp);
        }
    } else if iyear_AD > 0 {
        println!(
            "Orbital parameters calculated for given year: {} AD",
            iyear_AD
        );
    } else {
        println!(
            "Orbital parameters calculated for given year: {} BC",
            iyear_AD
        );
    }
}
