//! Radiation module, equivalent of radmod.

use crate::constants::{NHOR, NLEM, NLEP, NLEV, NLON, NLPP, NTRU, SBK, TWOPI};
use crate::{max, min, FloatNum, Int, Sim, Vec2d};

use num_traits::Float;

/// Minimum value for eccen
pub const ORB_ECCEN_MIN: FloatNum = 0.0;
/// Maximum value for eccen
pub const ORB_ECCEN_MAX: FloatNum = 0.1;
/// Minimum value for obliq
pub const ORB_OBLIQ_MIN: FloatNum = -90.0;
/// Maximum value for obliq
pub const ORB_OBLIQ_MAX: FloatNum = 90.0;
/// Minimum value for mvelp
pub const ORB_MVELP_MIN: FloatNum = 0.0;
/// Maximum value for mvelp
pub const ORB_MVELP_MAX: FloatNum = 360.0;
/// undefined/unset/invalid value
pub const ORB_UNDEF_REAL: FloatNum = 1.; // 1.e36;
/// flag to use default orbit
pub const ORB_DEFAULT: FloatNum = ORB_UNDEF_REAL;
/// undefined/unset/invalid value
pub const ORB_UNDEF_INT: Int = 2000000000;
/// flag to not use input year
pub const ORB_NOT_YEAR_BASED: Int = ORB_UNDEF_INT;

#[derive(Clone)]
pub struct Radiation {
    // configurables
    /// Solar constant (set in planet module)
    pub gsol0: FloatNum,
    /// Cos of lat of insolation if ncstsol=1
    pub solclat: FloatNum,
    /// Cos of dec of insolation if ncstsol=1
    pub solcdec: FloatNum,
    /// Cloud grayness (-1 = computed)
    pub clgray: FloatNum,
    /// Absorption coefficient h2o continuum (lwr)
    pub th2oc: FloatNum,
    /// Tuning of cloud albedo range1
    pub tswr1: FloatNum,
    /// Tuning of cloud back scattering c. range2
    pub tswr2: FloatNum,
    /// Tuning of cloud s. scattering alb. range2
    pub tswr3: FloatNum,

    /// Solar constant (set in planet module)
    pub tpofmt: FloatNum,
    /// Solar constant (set in planet module)
    pub acllwr: FloatNum,

    /// Parameter to define o3 profile
    pub a0o3: FloatNum,
    /// Parameter to define o3 profile
    pub a1o3: FloatNum,
    /// Parameter to define o3 profile
    pub aco3: FloatNum,
    /// Parameter to define o3 profile
    pub bo3: FloatNum,
    /// Parameter to define o3 profile
    pub co3: FloatNum,
    /// Parameter to define o3 profile
    pub toffo3: FloatNum,
    /// Scale o3 concentration
    pub o3scale: FloatNum,
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
    pub rcl1: [FloatNum; 3],
    /// Cloud albedos spectral range 2
    pub rcl2: [FloatNum; 3],
    /// Cloud absorptivities spectral range 2
    pub acl2: [FloatNum; 3],

    // arrays
    /// Cosine of solar zenit angle
    pub gmu0: Vec<FloatNum>, // NHOR
    /// Cosine of solar zenit angle
    pub gmu1: Vec<FloatNum>, // NHOR
    /// Ozon concentration (kg/kg)
    pub dqo3: Vec2d<FloatNum>, // (NHOR, NLEV)
    /// co2 concentration (ppmv)
    pub dqco2: Vec2d<FloatNum>, // (NHOR, NLEV)
    /// lwr temperature tendencies
    pub dtdtlwr: Vec2d<FloatNum>, // (NHOR, NLEV)
    /// swr temperature tendencies
    pub dtdtswr: Vec2d<FloatNum>, // (NHOR, NLEV)
    /// Climatological O3 (used if NO3=2)
    pub dqo3cl: Vec<Vec<Vec<FloatNum>>>,

    // scalars
    /// Earth-sun distance factor ( i.e. (1/r)**2 )
    pub gdist2: FloatNum,
    /// CPU time for radiation
    pub time4rad: FloatNum,
    /// CPU time for short wave radiation
    pub time4swr: FloatNum,
    /// CPU time for long wave radiation
    pub time4lwr: FloatNum,

    // orbital parameters
    pub ORB_UNDEF_INT: Int,
    /// Earth's obliquity in radians
    pub obliqr: FloatNum,
    /// Mean longitude of perihelion at the vernal equinox (radians)
    pub lambm0: FloatNum,
    /// Earth's moving vernal equinox longitude of perihelion plus pi (radians)
    pub mvelpp: FloatNum,
    /// Earth-sun distance factor ( i.e. (1/r)**2 )
    pub eccf: FloatNum,
    /// Year AD to calculate orbit for
    pub iyrad: Int,
    /// Flag to print-out status information or not.
    /// (This turns off ALL status printing including error messages)
    pub log_print: bool,

    // extended entropy/energy diagnostics
    pub dftde1: Vec2d<FloatNum>, // (NHOR, NLEP)
    pub dftde2: Vec2d<FloatNum>, // (NHOR, NLEP)
    pub dftue1: Vec2d<FloatNum>, // (NHOR, NLEP)
    pub dftue2: Vec2d<FloatNum>, // (NHOR, NLEP)
    pub dftu0: Vec2d<FloatNum>,  // (NHOR, NLEP)
    pub dftd0: Vec2d<FloatNum>,  // (NHOR, NLEP)

    // auxiliary variables for solar zenit angle calculations
    /// cos(lat)*cos(decl)
    pub solclatcdec: FloatNum,
    /// sin(lat)
    pub solslat: FloatNum,
    /// sin(decl)
    pub solsdec: FloatNum,
    /// sin(lat)*sin(decl)
    pub solslatsdec: FloatNum,
    /// temporary zenit angle
    pub zmuz: FloatNum,
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

            gmu0: vec![0.0; NHOR],
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
            dftd0: vec![],
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
        co2: FloatNum,
        // planet vars
        eccen: FloatNum,
        obliq: FloatNum,
        mvelp: FloatNum,
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
        );

        // read climatological ozone
        if self.no3 == 2 {
            self.dqo3cl = vec![vec![vec![0.0; 14]; NLEV]; NHOR]; // allocate(dqo3cl(NHOR,NLEV,0:13))

            // TODO probably not needed
            // call mpsurfgp('dqo3cl',dqo3cl,NHOR,NLEV*14)
        }

        // set co2 3d-field (enable external co2 by if statement)
        if co2 > 0. {
            self.dqco2 = vec![vec![co2; NLEV]; NHOR];
        }
    }
}

/// do the radiation calculations
/// this *sub* is called by PUMA (PUMA-interface)
///
/// no PUMA *subs* are used
///
/// the following PUMA variables are used/modified:
///
/// ga               : gravity accelleration (m/s2) (used)
/// acpd             : specific heat of dry air (J/kgK) (used)
/// ADV              : ACPV/acpd - 1  (used)
/// sigma(NLEV)      : sigma of T-levels (used)
/// dp(NHOR)         : surface pressure (Pa) (used)
/// dq(NHOR,NLEP)    : specific humidity (kg/kg) (used)
/// dtdt(NHOR,NLEP)  : temperature tendencies (K/s) (modified)
/// dswfl(NHOR,NLEP) : short wave radiation (W/m2)  (modified)
/// dlwfl(NHOR,NLEP) : long wave radiation (W/m2)   (modified)
/// dfu(NHOR,NLEP)   : short wave radiation upward (W/m2) (modified)
/// dfd(NHOR,NLEP)   : short wave radiation downward (W/m2) (modified)
/// dftu(NHOR,NLEP)  : long wave radiation upward (W/m2) (modified)
/// dftd(NHOR,NLEP)  : long wave radiation downward (W/m2) (modified)
/// dflux(NHOR,NLEP) : total radiation (W/m2) (modified)
///
/// the following radiation *subs* are called:
///
/// solang           : calc. cosine of solar zenit angle
/// mko3             : calc. ozon distribution
/// swr              : calc. short wave radiation fluxes
/// wr              : calc. long wave radiation fluxes
pub fn step(mut sim: &mut Sim) {
    //  0) define local arrays

    // temperature tendency due to rad (K/s)
    let mut zdtdt: Vec2d<FloatNum> = vec![vec![0.0; NLEV]; NHOR]; // (NHOR,NLEV)

    // allocatable arrays for diagnostic
    let zprf1: Vec2d<FloatNum>;
    let zprf2: Vec2d<FloatNum>;
    let zprf3: Vec2d<FloatNum>;
    let zprf4: Vec2d<FloatNum>;
    let zprf5: Vec2d<FloatNum>;
    let zprf6: Vec2d<FloatNum>;
    let zprf7: Vec2d<FloatNum>;
    let zprf8: Vec2d<FloatNum>;
    let zprf9: Vec2d<FloatNum>;
    let zprf10: Vec2d<FloatNum>;
    let zprf11: Vec2d<FloatNum>;
    let zprf12: Vec2d<FloatNum>;
    let zcc: Vec2d<FloatNum>;
    let zalb: Vec2d<FloatNum>;
    let zdtdte: Vec2d<FloatNum>;

    // cpu time estimates
    if sim.int_scalars.ntime == 1 {
        // TODO
        // mksecond(zsec, 0.);
    }

    // 1) set all fluxes to zero
    sim.rad_basic.dfu = vec![vec![0.0; NLEP]; NHOR]; // short wave radiation upward
    sim.rad_basic.dfd = vec![vec![0.0; NLEP]; NHOR]; // short wave radiation downward
    sim.rad_basic.dftu = vec![vec![0.0; NLEP]; NHOR]; // long wave radiation upward
    sim.rad_basic.dftd = vec![vec![0.0; NLEP]; NHOR]; // long wave radiation downward
    sim.rad_basic.dswfl = vec![vec![0.0; NLEP]; NHOR]; // total short wave radiation
    sim.rad_basic.dlwfl = vec![vec![0.0; NLEP]; NHOR]; // total long wave radiation
    sim.rad.dftue1 = vec![vec![0.0; NLEP]; NHOR]; // entropy
    sim.rad.dftue2 = vec![vec![0.0; NLEP]; NHOR]; // entropy

    // 2) compute cosine of solar zenit angle for each gridpoint
    if sim.rad.nsol == 1 {
        // TODO
        // solang();
    }

    // 3) compute ozon distribution
    if sim.rad.no3 < 3 {
        mko3(sim);
    }

    // 4) short wave radiation
    //
    // a) if clear sky diagnostic is switched on:
    if sim.int_scalars.ndiagcf > 0 {
        // TODO
        // allocate(zcc(NHOR,NLEP))
        // allocate(zalb(NHOR))
        // zcc(:,:)=dcc(:,:)
        // zalb(:)=dalb(:)
        // dcc(:,:)=0.
        // if(nswr==1) call swr
        // dclforc(:,1)=dswfl(:,NLEP)
        // dclforc(:,3)=dswfl(:,1)
        // dclforc(:,5)=dfu(:,1)
        // dclforc(:,6)=dfu(:,NLEP)
        // dcc(:,:)=zcc(:,:)
        // dalb(:)=zalb(:)
        // deallocate(zalb)
    }
    // b) normal computation
    if sim.int_scalars.ntime == 1 {
        // TODO
        // mksecond(zsec1, 0.);
    }
    if sim.rad.nswr == 1 {
        swr(sim);
    }
    if sim.int_scalars.ntime == 1 {
        // TODO
        // mksecond(zsec1, zsec1);
        // time4swr = time4swr + zsec1;
    }

    // 5) long wave radiation
    //
    // a) if clear sky diagnostic is switched on:
    if sim.int_scalars.ndiagcf > 0 {
        // TODO
        // zcc(:,:)=dcc(:,:)
        // dcc(:,:)=0.
        // if(nlwr==1) call lwr
        // dclforc(:,2)=dlwfl(:,NLEP)
        // dclforc(:,4)=dlwfl(:,1)
        // dclforc(:,7)=dftu(:,NLEP)
        // dcc(:,:)=zcc(:,:)
        // deallocate(zcc)
    }
    // b) normal computation
    if sim.int_scalars.ntime == 1 {
        // TODO
        // mksecond(zsec1, 0.);
    }
    if sim.rad.nlwr == 1 {
        lwr(sim);
    }
    if sim.int_scalars.ntime == 1 {
        // TODO
        // mksecond(zsec1, zsec1);
        // time4lwr = time4lwr + zsec1;
    }

    // 6) Total flux
    for (i, n) in sim.rad_basic.dflux.iter_mut().enumerate() {
        for (j, m) in n.iter_mut().enumerate() {
            *m = sim.rad_basic.dlwfl[i][j] + sim.rad_basic.dswfl[i][j];
        }
    }

    // 7) compute tendencies and add them to PUMA dtdt
    for jlev in 0..NLEV {
        let jlep = jlev + 1;
        for (i, n) in zdtdt.iter_mut().enumerate() {
            n[jlev] = -sim.planet_vars.ga
                * (sim.rad_basic.dflux[i][jlep] - sim.rad_basic.dflux[i][jlev])
                / (sim.lev_arrays.dsigma[jlev]
                    * sim.grid_arrays_dim.dp[i]
                    * sim.planet_vars.acpd
                    * (1. + sim.planet_vars.adv * sim.grid_arrays_dim.dq[i][jlev]));
        }
        for (i, n) in sim.grid_arrays_dim.dtdt.iter_mut().enumerate() {
            n[jlev] = n[jlev] + zdtdt[i][jlev];
        }
        for (i, n) in sim.rad.dtdtswr.iter_mut().enumerate() {
            n[jlev] = sim.planet_vars.ga
                * (sim.rad_basic.dswfl[i][jlep] - sim.rad_basic.dswfl[i][jlev])
                / (sim.lev_arrays.dsigma[jlev]
                    * sim.grid_arrays_dim.dp[i]
                    * sim.planet_vars.acpd
                    * (1. + sim.planet_vars.adv * sim.grid_arrays_dim.dq[i][jlev]));
        }
        for (i, n) in sim.rad.dtdtlwr.iter_mut().enumerate() {
            n[jlev] = -sim.planet_vars.ga
                * (sim.rad_basic.dlwfl[i][jlep] - sim.rad_basic.dlwfl[i][jlev])
                / (sim.lev_arrays.dsigma[jlev]
                    * sim.grid_arrays_dim.dp[i]
                    * sim.planet_vars.acpd
                    * (1. + sim.planet_vars.adv * sim.grid_arrays_dim.dq[i][jlev]));
        }
    }

    // !
    // !**   8) dbug printout if nprint=2 (see pumamod)
    // !
    //
    //       if (nprint==2) then
    //        allocate(zprf1(NLON*NLAT,NLEP))
    //        allocate(zprf2(NLON*NLAT,NLEP))
    //        allocate(zprf3(NLON*NLAT,NLEP))
    //        allocate(zprf4(NLON*NLAT,NLEP))
    //        allocate(zprf5(NLON*NLAT,NLEP))
    //        allocate(zprf6(NLON*NLAT,NLEP))
    //        allocate(zprf7(NLON*NLAT,NLEP))
    //        allocate(zprf8(NLON*NLAT,NLEV))
    //        allocate(zprf9(NHOR,NLEV))
    //        allocate(zprf10(NHOR,NLEV))
    //        allocate(zprf11(NLON*NLAT,NLEV))
    //        allocate(zprf12(NLON*NLAT,NLEV))
    //        call mpgagp(zprf1,dfd,NLEP)
    //        call mpgagp(zprf2,dfu,NLEP)
    //        call mpgagp(zprf3,dswfl,NLEP)
    //        call mpgagp(zprf4,dftd,NLEP)
    //        call mpgagp(zprf5,dftu,NLEP)
    //        call mpgagp(zprf6,dlwfl,NLEP)
    //        call mpgagp(zprf7,dflux,NLEP)
    //        call mpgagp(zprf8,zdtdt,NLEV)
    //        do jlev = 1 , NLEV
    //         jlep=jlev+1
    //         zprf9(:,jlev)=-ga*(dswfl(:,jlep)-dswfl(:,jlev))                 &
    //      &               /(dsigma(jlev)*dp(:)*acpd*(1.+ADV*dq(:,jlev)))
    //         zprf10(:,jlev)=-ga*(dlwfl(:,jlep)-dlwfl(:,jlev))                &
    //      &                /(dsigma(jlev)*dp(:)*acpd*(1.+ADV*dq(:,jlev)))
    //        enddo
    //        call mpgagp(zprf11,zprf9,NLEV)
    //        call mpgagp(zprf12,zprf10,NLEV)
    //        if(mypid==NROOT) then
    //         do jlev=1,NLEP
    //          write(nud,*)'L= ',jlev,' swd= ',zprf1(nprhor,jlev)                  &
    //      &                    ,' swu= ',zprf2(nprhor,jlev)                  &
    //      &                    ,' swt= ',zprf3(nprhor,jlev)
    //          write(nud,*)'L= ',jlev,' lwd= ',zprf4(nprhor,jlev)                  &
    //      &                    ,' lwu= ',zprf5(nprhor,jlev)                  &
    //      &                    ,' lwt= ',zprf6(nprhor,jlev)
    //          write(nud,*)'L= ',jlev,' totalflux= ',zprf7(nprhor,jlev)
    //         enddo
    //         do jlev=1,NLEV
    //          write(nud,*)'L= ',jlev,' dtdt= ',zprf8(nprhor,jlev)                 &
    //      &                    ,' dtsw= ',zprf11(nprhor,jlev)                &
    //      &                    ,' dtlw= ',zprf12(nprhor,jlev)
    //         enddo
    //        endif
    //        deallocate(zprf1)
    //        deallocate(zprf2)
    //        deallocate(zprf3)
    //        deallocate(zprf4)
    //        deallocate(zprf5)
    //        deallocate(zprf6)
    //        deallocate(zprf7)
    //        deallocate(zprf8)
    //        deallocate(zprf9)
    //        deallocate(zprf10)
    //        deallocate(zprf11)
    //        deallocate(zprf12)
    //       endif
    //
    // !
    // !     franks dbug
    // !
    //
    //       if(ndiaggp==1) then
    //        do jlev = 1 , NLEV
    //         jlep=jlev+1
    //         dgp3d(:,jlev,5)=-ga*(dlwfl(:,jlep)-dlwfl(:,jlev))               &
    //      &               /(dsigma(jlev)*dp(:)*acpd*(1.+ADV*dq(:,jlev)))
    //         dgp3d(:,jlev,6)=-ga*(dswfl(:,jlep)-dswfl(:,jlev))               &
    //      &               /(dsigma(jlev)*dp(:)*acpd*(1.+ADV*dq(:,jlev)))
    //        enddo
    //       end if
    // !
    // !     entropy/energy diagnostics
    // !
    //       if(nentropy > 0) then
    //        dentropy(:,16)=dlwfl(:,NLEP)/dt(:,NLEP)
    //        dentropy(:,17)=dswfl(:,NLEP)/dt(:,NLEP)
    //        dentropy(:,27)=dftd(:,NLEP)/dt(:,NLEP)
    //        dentropy(:,28)=dftu(:,NLEP)/dt(:,NLEP)
    //        dentropy(:,9)=0.
    //        dentropy(:,10)=0.
    //        dentropy(:,21)=0.
    //        dentropy(:,22)=0.
    //        dentropy(:,23)=0.
    //        dentropy(:,24)=0.
    //        dentropy(:,26)=0.
    //        dentropy(:,29)=0.
    //        dentropy(:,30)=0.
    //        do jlev=1,NLEV
    //         jlep=jlev+1
    //         dentro(:)=dftu0(:,jlev)/dentrot(:,jlev)
    //         dentropy(:,29)=dentropy(:,29)+dentro(:)
    //         if(nentro3d > 0) dentro3d(:,jlev,19)=dentro(:)
    //         dentro(:)=dftd0(:,jlev)/dentrot(:,jlev)
    //         dentropy(:,30)=dentropy(:,30)+dentro(:)
    //         if(nentro3d > 0) dentro3d(:,jlev,20)=dentro(:)
    //         dentro(:)=-ga*(dlwfl(:,jlep)-dlwfl(:,jlev))                     &
    //      &           /(dsigma(jlev)*dp(:)*acpd*(1.+ADV*dq(:,jlev)))         &
    //      &         *acpd*(1.+adv*dentroq(:,jlev))*dentrop(:)/ga*dsigma(jlev)&
    //      &         /dentrot(:,jlev)
    //         dentropy(:,9)=dentropy(:,9)+dentro(:)
    //         if(nentro3d > 0) dentro3d(:,jlev,9)=dentro(:)
    //         dentro(:)=-ga*(dswfl(:,jlep)-dswfl(:,jlev))                     &
    //      &           /(dsigma(jlev)*dp(:)*acpd*(1.+ADV*dq(:,jlev)))         &
    //      &         *acpd*(1.+adv*dentroq(:,jlev))*dentrop(:)/ga*dsigma(jlev)&
    //      &         /dentrot(:,jlev)
    //         dentropy(:,10)=dentropy(:,10)+dentro(:)
    //         if(nentro3d > 0) dentro3d(:,jlev,10)=dentro(:)
    //         dentro(:)=-ga*(dftd(:,jlep)-dftd(:,jlev))                       &
    //      &           /(dsigma(jlev)*dp(:)*acpd*(1.+ADV*dq(:,jlev)))         &
    //      &         *acpd*(1.+adv*dentroq(:,jlev))*dentrop(:)/ga*dsigma(jlev)&
    //      &         /dentrot(:,jlev)
    //         dentropy(:,21)=dentropy(:,21)+dentro(:)
    //         if(nentro3d > 0) dentro3d(:,jlev,15)=dentro(:)
    //         dentro(:)=-ga*(dftu(:,jlep)-dftu(:,jlev))                       &
    //      &           /(dsigma(jlev)*dp(:)*acpd*(1.+ADV*dq(:,jlev)))         &
    //      &         *acpd*(1.+adv*dentroq(:,jlev))*dentrop(:)/ga*dsigma(jlev)&
    //      &         /dentrot(:,jlev)
    //         dentropy(:,22)=dentropy(:,22)+dentro(:)
    //         if(nentro3d > 0) dentro3d(:,jlev,16)=dentro(:)
    //         dentro(:)=-ga*(dftue1(:,jlep)-dftue1(:,jlev))                   &
    //      &           /(dsigma(jlev)*dp(:)*acpd*(1.+ADV*dq(:,jlev)))         &
    //      &         *acpd*(1.+adv*dentroq(:,jlev))*dentrop(:)/ga*dsigma(jlev)&
    //      &         /dentrot(:,jlev)
    //         dentropy(:,23)=dentropy(:,23)+dentro(:)
    //         if(nentro3d > 0) dentro3d(:,jlev,17)=dentro(:)
    //         dentro(:)=-ga*(dftue2(:,jlep)-dftue2(:,jlev))                   &
    //      &           /(dsigma(jlev)*dp(:)*acpd*(1.+ADV*dq(:,jlev)))         &
    //      &         *acpd*(1.+adv*dentroq(:,jlev))*dentrop(:)/ga*dsigma(jlev)&
    //      &         /dentrot(:,jlev)
    //         dentropy(:,24)=dentropy(:,24)+dentro(:)
    //         if(nentro3d > 0) dentro3d(:,jlev,18)=dentro(:)
    //         dentropy(:,26)=dentropy(:,26)+dentro(:)*dentrot(:,jlev)
    //        enddo
    //        dentropy(:,25)=(dftu(:,NLEP)+dentropy(:,26))/dt(:,NLEP)
    //        dentropy(:,26)=-dentropy(:,26)/dt(:,NLEP)
    //       endif
    //       if(nenergy > 0) then
    //        allocate(zdtdte(NHOR,NLEV))
    //        denergy(:,9)=0.
    //        denergy(:,10)=0.
    //        denergy(:,17)=0.
    //        denergy(:,18)=0.
    //        denergy(:,19)=0.
    //        denergy(:,20)=0.
    //        denergy(:,28)=0.
    //        do jlev=1,NLEV
    //         jlep=jlev+1
    //         denergy(:,9)=denergy(:,9)-(dlwfl(:,jlep)-dlwfl(:,jlev))
    //         denergy(:,10)=denergy(:,10)-(dswfl(:,jlep)-dswfl(:,jlev))
    //         denergy(:,17)=denergy(:,17)-(dftd(:,jlep)-dftd(:,jlev))
    //         denergy(:,18)=denergy(:,18)-(dftu(:,jlep)-dftu(:,jlev))
    //         denergy(:,19)=denergy(:,19)-(dftue1(:,jlep)-dftue1(:,jlev))
    //         denergy(:,20)=denergy(:,20)-(dftue2(:,jlep)-dftue2(:,jlev))
    //         denergy(:,28)=denergy(:,28)+dt(:,jlev)*dsigma(jlev)
    //         if(nener3d > 0) then
    //          dener3d(:,jlev,9)=-(dlwfl(:,jlep)-dlwfl(:,jlev))
    //          dener3d(:,jlev,10)=-(dswfl(:,jlep)-dswfl(:,jlev))
    //          dener3d(:,jlev,17)=-(dftd(:,jlep)-dftd(:,jlev))
    //          dener3d(:,jlev,18)=-(dftu(:,jlep)-dftu(:,jlev))
    //          dener3d(:,jlev,19)=-(dftue1(:,jlep)-dftue1(:,jlev))
    //          dener3d(:,jlev,20)=-(dftue2(:,jlep)-dftue2(:,jlev))
    //          dener3d(:,jlev,28)=dt(:,jlev)*dsigma(jlev)
    //         endif
    //        enddo
    //        deallocate(zdtdte)
    //       endif
    //
    //       if(ntime == 1) then
    //        call mksecond(zsec,zsec)
    //        time4rad=time4rad+zsec
    //       endif
}

/// compute ozon distribution
///
/// used subroutines:
///
/// ndayofyear : module calmod - compute day of year
///
/// the following variables from pumamod are used:
///
/// TWOPI    : 2*PI
/// nstep    : PLASIM time step
/// sid(NLAT): double precision sines of gaussian latitudes
fn mko3(mut sim: &mut Sim) {
    // TODO what is this exactly?
    // interface
    //      integer function ndayofyear(kstep)
    //         integer, intent(in) :: kstep
    //      end function ndayofyear
    //   end interface

    let mut jlat: Int; // latitude index
    let mut jlev: Int; // level index
    let mut jh1: Int;
    let mut jh2: Int;
    let mut imonth: Int = 0; // current month
    let mut jm: Int = 0; // index to next or previous month

    const zt0: FloatNum = 255.;
    const zro3: FloatNum = 2.14;
    const zfo3: FloatNum = 100.0 / zro3;

    let mut zcday: FloatNum = 0.; // current day
    let zconst: FloatNum;
    let zw: FloatNum = 0.; // interpolation weight

    let mut za: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR)
    let mut zh: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR)
    let mut zo3t: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR)
    let mut zo3: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR)

    // compute synthetic ozone distribution
    if sim.rad.no3 == 1 {
        // TODO create and use calendar interface here to get ndayofyear
        // zcday = ndayofyear(nstep);

        for jlat in 0..NLPP {
            // horizontal index for end of latitude
            jh2 = (jlat * NLON) as Int;
            // horizontal index for start of latitude
            jh1 = jh2 - NLON as Int + 1;

            // ndatim(6) = leap year
            for i in jh1 as usize..jh2 as usize + 1 {
                za[i] = sim.rad.a0o3
                    + sim.rad.a1o3 * (sim.lat_arrays.sid[jlat]).abs()
                    + sim.rad.aco3
                        * sim.lat_arrays.sid[jlat]
                        * (TWOPI * (zcday - sim.rad.toffo3)
                            / (sim.datetime.n_days_per_year + sim.datetime.ndatim[6]) as FloatNum)
                            .cos();
            }
        } // jlat

        zconst = (-sim.rad.bo3 / sim.rad.co3).exp();
        for (i, n) in zo3t.iter_mut().enumerate() {
            *n = za[i];
        }
        for n in zh.iter_mut() {
            *n = 0.0;
        }

        for jlev in (1..NLEV).rev() {
            for (i, n) in zh.iter_mut().enumerate() {
                *n = *n
                    - sim.grid_arrays_dim.dt[i][jlev] * sim.planet_vars.gascon / sim.planet_vars.ga
                        * (sim.lev_arrays.sigmah[jlev - 1] / sim.lev_arrays.sigmah[jlev]).ln();
            }
            for (i, n) in zo3.iter_mut().enumerate() {
                *n = -(za[i] + za[i] * zconst) / (1. + ((zh[i] - sim.rad.bo3) / sim.rad.co3)).exp()
                    + zo3t[i];
            }
            for (i, n) in sim.rad.dqo3.iter_mut().enumerate() {
                n[jlev] = zo3[i] * sim.planet_vars.ga
                    / (zfo3 * sim.lev_arrays.dsigma[jlev] * sim.grid_arrays_dim.dp[i]);
            }
            for (i, n) in zo3t.iter_mut().enumerate() {
                *n = *n - zo3[i];
            }
        }

        for (i, n) in sim.rad.dqo3.iter_mut().enumerate() {
            n[0] = zo3t[i] * sim.planet_vars.ga
                / (zfo3 * sim.lev_arrays.dsigma[0] * sim.grid_arrays_dim.dp[i]);
        }
    }
    // else interpolate from climatological ozone
    else if sim.rad.no3 == 2 {
        for (i, n) in sim.rad.dqo3.iter_mut().enumerate() {
            for (j, m) in n.iter_mut().enumerate() {
                *m = (1.0 - zw) * sim.rad.dqo3cl[i][j][imonth as usize]
                    + zw * sim.rad.dqo3cl[i][j][jm as usize];
            }
        }
    }

    if sim.rad.o3scale != 1.0 {
        for (i, n) in sim.rad.dqo3.iter_mut().enumerate() {
            for (j, m) in n.iter_mut().enumerate() {
                *m = sim.rad.o3scale * *m;
            }
        }
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

    const zero: FloatNum = 0.000001; // if insolation < zero : fluxes=0.
    const zsolar1: FloatNum = 0.517; // spectral partitioning 1 (wl < 0.75mue)
    const zsolar2: FloatNum = 0.483; // spectral partitioning 2 (wl > 0.75mue)
    const zbetta: FloatNum = 1.66; // magnification factor water vapour
    const zmbar: FloatNum = 1.9; // magnification factor ozon
    const zro3: FloatNum = 2.14; // ozon density (kg/m**3 STP)
    const zfo3: FloatNum = 100. / zro3; // transfere o3 to cm STP

    // transmissivities 1-l
    let mut zt1: Vec<Vec<FloatNum>> = vec![vec![0.; NLEP]; NHOR]; // (NHOR,NLEP);
    let mut zt2: Vec<Vec<FloatNum>> = vec![vec![0.; NLEP]; NHOR]; // (NHOR,NLEP);

    // reflexivities l-1 (scattered)
    let mut zr1s: Vec<Vec<FloatNum>> = vec![vec![0.; NLEP]; NHOR]; // (NHOR,NLEP);
    let mut zr2s: Vec<Vec<FloatNum>> = vec![vec![0.; NLEP]; NHOR]; // (NHOR,NLEP);

    // reflexivities l-NL (direct)
    let mut zrl1: Vec<Vec<FloatNum>> = vec![vec![0.; NLEP]; NHOR]; // (NHOR,NLEP);
    let mut zrl2: Vec<Vec<FloatNum>> = vec![vec![0.; NLEP]; NHOR]; // (NHOR,NLEP);

    // reflexivities l-NL (scattered)
    let mut zrl1s: Vec<Vec<FloatNum>> = vec![vec![0.; NLEP]; NHOR]; // (NHOR,NLEP);
    let mut zrl2s: Vec<Vec<FloatNum>> = vec![vec![0.; NLEP]; NHOR]; // (NHOR,NLEP);

    // layer transmissivity (down)
    let mut ztb1: Vec<Vec<FloatNum>> = vec![vec![0.; NLEV]; NHOR]; // (NHOR,NLEV);
    let mut ztb2: Vec<Vec<FloatNum>> = vec![vec![0.; NLEV]; NHOR]; // (NHOR,NLEV);

    // layer transmissivity (up)
    let mut ztb1u: Vec<Vec<FloatNum>> = vec![vec![0.; NLEV]; NHOR]; // (NHOR,NLEV);
    let mut ztb2u: Vec<Vec<FloatNum>> = vec![vec![0.; NLEV]; NHOR]; // (NHOR,NLEV);

    //  layer reflexivity (direct)
    let mut zrb1: Vec<Vec<FloatNum>> = vec![vec![0.; NLEV]; NHOR]; // (NHOR,NLEV);
    let mut zrb2: Vec<Vec<FloatNum>> = vec![vec![0.; NLEV]; NHOR]; // (NHOR,NLEV);

    //  layer reflexibity (scattered)
    let mut zrb1s: Vec<Vec<FloatNum>> = vec![vec![0.; NLEV]; NHOR]; // (NHOR,NLEV);
    let mut zrb2s: Vec<Vec<FloatNum>> = vec![vec![0.; NLEV]; NHOR]; // (NHOR,NLEV);

    let mut zo3l: Vec<Vec<FloatNum>> = vec![vec![0.; NLEV]; NHOR]; // (NHOR,NLEV); // ozon amount (top-l)
    let mut zxo3l: Vec<Vec<FloatNum>> = vec![vec![0.; NLEV]; NHOR]; // (NHOR,NLEV); // effective ozon amount (top-l)
    let mut zwvl: Vec<Vec<FloatNum>> = vec![vec![0.; NLEV]; NHOR]; // (NHOR,NLEV); // water vapor amount (top-l)
    let mut zywvl: Vec<Vec<FloatNum>> = vec![vec![0.; NLEV]; NHOR]; // (NHOR,NLEV); // effective water vapor amount (top-l)
    let mut zrcs: Vec<Vec<FloatNum>> = vec![vec![0.; NLEV]; NHOR]; // (NHOR,NLEV); // clear sky reflexivity (downward beam)
    let mut zrcsu: Vec<Vec<FloatNum>> = vec![vec![0.; NLEV]; NHOR]; // (NHOR,NLEV); // clear sky reflexivity (upward beam)

    // top solar radiation
    let mut zftop1: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);
    let mut zftop2: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);

    // upward fluxes
    let mut zfu1: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);
    let mut zfu2: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);

    // downward fluxes
    let mut zfd1: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);
    let mut zfd2: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);

    let mut zmu0: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR); zenit angle
    let mut zmu1: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR); zenit angle
    let mut zcs: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR); clear sky part
    let mut zm: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR); magnification factor
    let mut zo3: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR); ozon amount
    let mut zo3t: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR); total ozon amount (top-sfc)
    let mut zxo3t: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR); effective total ozon amount (top-sfc)

    // ozon transmissivity (downward/upward beam)
    let mut zto3: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR); zenit angle
    let mut zto3u: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR); zenit angle

    // total ozon transmissivities (d/u)
    let mut zto3t: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);
    let mut zto3tu: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);

    let mut zwv: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR); water vapor amount
    let mut zwvt: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR); total water vapor amount (top-sfc)
    let mut zywvt: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR); total effective water vapor amount (top-sfc)

    // water vapor trasmissivity (d/u)
    let mut ztwv: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);
    let mut ztwvu: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);

    // total water vapor transmissivities (d/u)
    let mut ztwvt: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);
    let mut ztwvtu: Vec<FloatNum> = vec![0.; NHOR]; //

    // reflexivities combined layer (direct)
    let mut zra1: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);
    let mut zra2: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);

    // reflexivities combined layer (scatterd)
    let mut zra1s: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);
    let mut zra2s: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);

    // transmissivities combined layer (di)
    let mut zta1: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);
    let mut zta2: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);

    // transmissivities combined layer (sc)
    let mut zta1s: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);
    let mut zta2s: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);

    let mut z1mrabr: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR); 1/(1.-rb*ra(*))

    // cloud reflexivities (direct)
    let mut zrcl1: Vec<Vec<FloatNum>> = vec![vec![0.; NLEV]; NHOR]; // (NHOR,NLEV);
    let mut zrcl2: Vec<Vec<FloatNum>> = vec![vec![0.; NLEV]; NHOR]; // (NHOR,NLEV);
                                                                    // cloud reflexivities (scattered)
    let mut zrcl1s: Vec<Vec<FloatNum>> = vec![vec![0.; NLEV]; NHOR]; // (NHOR,NLEV);
    let mut zrcl2s: Vec<Vec<FloatNum>> = vec![vec![0.; NLEV]; NHOR]; // (NHOR,NLEV);

    // cloud transmissivities
    let mut ztcl2: Vec<Vec<FloatNum>> = vec![vec![0.; NLEV]; NHOR]; // (NHOR,NLEV);
    let mut ztcl2s: Vec<Vec<FloatNum>> = vec![vec![0.; NLEV]; NHOR]; // (NHOR,NLEV);

    // arrays for diagnostic cloud properties
    let mut zlwp: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);
    let mut ztau: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);
    let mut zlog: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);
    let mut zb2: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);
    let mut zom0: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);
    let mut zuz: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);
    let mut zun: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);
    let mut zr: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);
    let mut zexp: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);
    let mut zu: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);
    let mut zb1: Vec<FloatNum> = vec![0.; NHOR]; // (NHOR);

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
                let zsum: FloatNum = sim.rad.gmu0[js..je].iter().sum();
                zmu0[js..je]
                    .iter_mut()
                    .map(|n| *n = zsum / icnt as FloatNum); // used for clouds
                zmu1[js..je]
                    .iter_mut()
                    .map(|n| *n = zsum / NLON as FloatNum); // used for insolation
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
    let zmu00: FloatNum = 0.5;
    let zb3 = sim.rad.tswr1 * zmu00.sqrt() / zmu00;
    let zb4 = sim.rad.tswr2 * zmu00.sqrt();
    let zb5 = sim.rad.tswr3 * zmu00 * zmu00;

    // prescribed

    if sim.rad.nswrcl == 0 {
        for jlev in 0..NLEV {
            if sim.lev_arrays.sigma[jlev] <= 1. / 3. {
                zrcl1s
                    .iter_mut()
                    .map(|n| n[jlev] = sim.rad.rcl1[0] / (sim.rad.rcl1[0] + zmu00));
                zrcl1.iter_mut().enumerate().map(|(i, n)| {
                    n[jlev] = zcs[i] * sim.rad.rcl1[0] / (sim.rad.rcl1[0] + zmu0[i])
                        + (1. - zcs[i]) * zrcl1s[i][jlev]
                });
                zrcl2s.iter_mut().map(|n| {
                    n[jlev] =
                        (1. - sim.rad.acl2[0]).min(sim.rad.rcl1[0] / (sim.rad.rcl2[0] + zmu00))
                });
                zrcl2.iter_mut().enumerate().map(|(i, n)| {
                    n[jlev] = (1. - sim.rad.acl2[0]).min(
                        zcs[i] * sim.rad.rcl2[0] / (sim.rad.rcl2[0] + zmu0[i])
                            + (1. - zcs[i]) * zrcl2s[i][jlev],
                    )
                });
                ztcl2s
                    .iter_mut()
                    .enumerate()
                    .map(|(i, n)| n[jlev] = 1. - zrcl2s[i][jlev] - sim.rad.acl2[0]);
                ztcl2
                    .iter_mut()
                    .enumerate()
                    .map(|(i, n)| n[jlev] = 1. - zrcl2[i][jlev] - sim.rad.acl2[0]);
            } else if sim.lev_arrays.sigma[jlev] > 1. / 3. && sim.lev_arrays.sigma[jlev] <= 2. / 3.
            {
                zrcl1s
                    .iter_mut()
                    .enumerate()
                    .map(|(i, n)| n[jlev] = sim.rad.rcl1[1] / (sim.rad.rcl1[1] + zmu00));
                zrcl1.iter_mut().enumerate().map(|(i, n)| {
                    n[jlev] = zcs[i] * sim.rad.rcl1[1] / (sim.rad.rcl1[1] + zmu0[i])
                        + (1. - zcs[i]) * zrcl1s[i][jlev]
                });
                zrcl2s
                    .iter_mut()
                    .enumerate()
                    .map(|(i, n)| n[jlev] = (1. - sim.rad.acl2[1]).min(sim.rad.rcl2[1] + zmu00));
                zrcl2.iter_mut().enumerate().map(|(i, n)| {
                    n[jlev] = min(
                        1. - sim.rad.acl2[1],
                        zcs[i] * sim.rad.rcl2[1] / (sim.rad.rcl2[1] + zmu0[i])
                            + (1. - zcs[i]) * zrcl2s[i][jlev],
                    )
                });
                ztcl2s
                    .iter_mut()
                    .enumerate()
                    .map(|(i, n)| n[jlev] = 1. - zrcl2s[i][jlev] - sim.rad.acl2[1]);
                ztcl2
                    .iter_mut()
                    .enumerate()
                    .map(|(i, n)| n[jlev] = 1. - zrcl2[i][jlev] - sim.rad.acl2[1]);
            } else {
                zrcl1s
                    .iter_mut()
                    .enumerate()
                    .map(|(i, n)| n[jlev] = sim.rad.rcl1[2] / (sim.rad.rcl1[2] + zmu00));
                zrcl1.iter_mut().enumerate().map(|(i, n)| {
                    n[jlev] = zcs[i] * sim.rad.rcl1[2] / (sim.rad.rcl1[2] + zmu0[i])
                        + (1. - zcs[i]) * zrcl1s[i][jlev]
                });
                zrcl2s.iter_mut().enumerate().map(|(i, n)| {
                    n[jlev] = min(
                        1. - sim.rad.acl2[2],
                        sim.rad.rcl2[2] / (sim.rad.rcl2[2] + zmu00),
                    )
                });
                zrcl2.iter_mut().enumerate().map(|(i, n)| {
                    n[jlev] = min(
                        1. - sim.rad.acl2[2],
                        zcs[i] * sim.rad.rcl2[2] / (sim.rad.rcl2[2] + zmu0[i])
                            + (1. - zcs[i]) * zrcl2s[i][jlev],
                    )
                });
                ztcl2s
                    .iter_mut()
                    .enumerate()
                    .map(|(i, n)| n[jlev] = 1. - zrcl2s[i][jlev] - sim.rad.acl2[2]);
                ztcl2
                    .iter_mut()
                    .enumerate()
                    .map(|(i, n)| n[jlev] = 1. - zrcl2[i][jlev] - sim.rad.acl2[2]);
            }
            zcs.iter_mut()
                .enumerate()
                .map(|(i, n)| *n = *n * (1. - sim.grid_arrays_dim.dcc[i][jlev]));
        }
    } else {
        zrcl1 = vec![vec![0.; NLEV]; NHOR];
        zrcl2 = vec![vec![0.; NLEV]; NHOR];
        ztcl2 = vec![vec![1.; NLEV]; NHOR];
        zrcl1s = vec![vec![0.; NLEV]; NHOR];
        zrcl2s = vec![vec![0.; NLEV]; NHOR];
        ztcl2s = vec![vec![1.; NLEV]; NHOR];

        for jlev in 0..NLEV {
            // TODO this sort of iteration is probably not the most efficient
            for nh in 0..NHOR {
                if losun[nh] && sim.grid_arrays_dim.dcc[nh][jlev] > 0. {
                    zlwp[nh] = min(
                        1000.0,
                        1000. * sim.grid_arrays_dim.dql[nh][jlev] * sim.grid_arrays_dim.dp[nh]
                            / sim.planet_vars.ga
                            * sim.lev_arrays.dsigma[jlev],
                    );
                    ztau[nh] = 2.0 * ((zlwp[nh] + 1.5).log10()).powf(3.9);
                    zlog[nh] = (1000.0 / ztau[nh]).ln();
                    zb2[nh] = zb4 / (3. + 0.1 * ztau[nh]).ln();
                    zom0[nh] = (0.9999 as FloatNum).min(1.0 - zb5 * zlog[nh]);
                    zun[nh] = 1.0 - zom0[nh];
                    zuz[nh] = zun[nh] + 2.0 * zb2[nh] * zom0[nh];
                    zu[nh] = (zuz[nh] / zun[nh]).sqrt();
                    zexp[nh] = (min(25.0, ztau[nh] * (zuz[nh] * zun[nh]).sqrt()) / zmu00).exp();
                    zr[nh] = (zu[nh] + 1.) * (zu[nh] + 1.) * zexp[nh]
                        - (zu[nh] - 1.) * (zu[nh] - 1.) / zexp[nh];
                    zrcl1s[nh][jlev] = 1. - 1. / (1. + zb3 * ztau[nh]);
                    ztcl2s[nh][jlev] = 4. * zu[nh] / zr[nh];
                    zrcl2s[nh][jlev] = (zu[nh] * zu[nh] - 1.) / zr[nh] * (zexp[nh] - 1. / zexp[nh]);

                    zb1[nh] = sim.rad.tswr1 * zmu0[nh].sqrt();
                    zb2[nh] = sim.rad.tswr2 * zmu0[nh].sqrt() / (3. + 0.1 * ztau[nh]).ln();
                    zom0[nh] = min(0.9999, 1. - sim.rad.tswr3 * zmu0[nh] * zmu0[nh] * zlog[nh]);
                    zun[nh] = 1.0 - zom0[nh];
                    zuz[nh] = zun[nh] + 2.0 * zb2[nh] * zom0[nh];
                    zu[nh] = (zuz[nh] / zun[nh]).sqrt();
                    zexp[nh] = (min(25.0, ztau[nh] * (zuz[nh] * zun[nh]).sqrt() / zmu0[nh])).exp();
                    zr[nh] = (zu[nh] + 1.) * (zu[nh] + 1.) * zexp[nh]
                        - (zu[nh] - 1.) * (zu[nh] - 1.) / zexp[nh];
                    zrcl1[nh][jlev] = 1. - 1. / (1. + zb1[nh] * ztau[nh] / zmu0[nh]);
                    ztcl2[nh][jlev] = 4. * zu[nh] / zr[nh];
                    zrcl2[nh][jlev] = (zu[nh] * zu[nh] - 1.) / zr[nh] * (zexp[nh] - 1. / zexp[nh]);
                    zrcl1[nh][jlev] = zcs[nh] * zrcl1[nh][jlev] + (1. - zcs[nh]) * zrcl1s[nh][jlev];
                    ztcl2[nh][jlev] = zcs[nh] * ztcl2[nh][jlev] + (1. - zcs[nh]) * ztcl2s[nh][jlev];
                    zrcl2[nh][jlev] = zcs[nh] * zrcl2[nh][jlev] + (1. - zcs[nh]) * zrcl2s[nh][jlev];
                }
            }
            zcs.iter_mut()
                .enumerate()
                .map(|(i, n)| *n = *n * (1. - sim.grid_arrays_dim.dcc[i][jlev]));
        } // jlev
    } // (nswrcl == 0)

    // magnification factor

    for nh in 0..NHOR {
        if losun[nh] {
            zm[nh] = 35. / (1. + 1224. * zmu0[nh] * zmu0[nh]).sqrt();

            // absorber amount and clear sky fraction
            zcs[nh] = 1.;
            zo3t[nh] = 0.;
            zxo3t[nh] = 0.;
            zwvt[nh] = 0.;
            zywvt[nh] = 0.;
        }
    }

    for jlev in 0..NLEV {
        for nh in 0..NHOR {
            if losun[nh] {
                zo3[nh] = zfo3
                    * sim.lev_arrays.dsigma[jlev]
                    * sim.grid_arrays_dim.dp[nh]
                    * sim.rad.dqo3[nh][jlev]
                    / sim.planet_vars.ga;
                zo3t[nh] = zo3t[nh] + zo3[nh];
                zxo3t[nh] = zcs[nh] * (zxo3t[nh] + zm[nh] * zo3[nh])
                    + (1. - zcs[nh]) * (zxo3t[nh] + zmbar * zo3[nh]);
                zo3l[nh][jlev] = zo3t[nh];
                zxo3l[nh][jlev] = zxo3t[nh];
                zwv[nh] = 0.1
                    * sim.lev_arrays.dsigma[jlev]
                    * sim.grid_arrays_dim.dq[nh][jlev]
                    * sim.grid_arrays_dim.dp[nh]
                    / sim.planet_vars.ga
                    * (273. / sim.grid_arrays_dim.dt[nh][jlev]).sqrt()
                    * sim.lev_arrays.sigma[jlev]
                    * sim.grid_arrays_dim.dp[nh]
                    / 100000.;
                zwvt[nh] = zwvt[nh] + zwv[nh];
                zywvt[nh] = zcs[nh] * (zywvt[nh] + zm[nh] * zwv[nh])
                    + (1. - zcs[nh]) * (zywvt[nh] + zbetta * zwv[nh]);
                zwvl[nh][jlev] = zwvt[nh];
                zywvl[nh][jlev] = zywvt[nh];
                zcs[nh] = zcs[nh] * (1. - sim.grid_arrays_dim.dcc[nh][jlev]);
                zrcs[nh][jlev] = 0.;
                zrcsu[nh][jlev] = 0.;
            }
        }
    }

    // compute optical properties
    //
    // downward loop
    //
    // preset

    for nh in 0..NHOR {
        if losun[nh] {
            zta1[nh] = 1.;
            zta1s[nh] = 1.;
            zra1[nh] = 0.;
            zra1s[nh] = 0.;
            zta2[nh] = 1.;
            zta2s[nh] = 1.;
            zra2[nh] = 0.;
            zra2s[nh] = 0.;

            zto3t[nh] = 1.;
            zo3[nh] = zxo3t[nh] + zmbar * zo3t[nh];
            zto3tu[nh] = 1.
                - (0.02118 * zo3[nh] / (1. + 0.042 * zo3[nh] + 0.000323 * (zo3[nh]).powf(2.))
                    + 1.082 * zo3[nh] / ((1. + 138.6 * zo3[nh]).powf(0.805))
                    + 0.0658 * zo3[nh] / (1. + (103.6 * zo3[nh]).powf(3.)))
                    / zsolar1;
            ztwvt[nh] = 1.;
            zwv[nh] = zywvt[nh] + zbetta * zwvt[nh];
            ztwvtu[nh] = 1.
                - 2.9 * zwv[nh] / ((1. + 141.5 * zwv[nh]).powf(0.635) + 5.925 * zwv[nh]) / zsolar2;

            //clear sky scattering (Rayleigh scatterin lower most level only)
            zrcs[nh][NLEV - 1] = (0.219 / (1. + 0.816 * zmu0[nh]) * zcs[nh]
                + 0.144 * (1. - zcs[nh] - sim.grid_arrays_dim.dcc[nh][NLEV - 1]))
                * sim.rad.nrscat as FloatNum;
            zrcsu[nh][NLEV - 1] =
                0.144 * (1. - sim.grid_arrays_dim.dcc[nh][NLEV - 1]) * sim.rad.nrscat as FloatNum;
        }
    }

    for jlev in 0..NLEV {
        for nh in 0..NHOR {
            if losun[nh] {
                zt1[nh][jlev] = zta1[nh];
                zt2[nh][jlev] = zta2[nh];
                zr1s[nh][jlev] = zra1s[nh];
                zr2s[nh][jlev] = zra2s[nh];

                // set single layer R and T:
                //
                // 1. spectral range 1:
                //
                // a) R
                //
                // clear part: rayleigh scattering (only lowermost level)
                // cloudy part: cloud albedo
                zrb1[nh][jlev] =
                    zrcs[nh][jlev] + zrcl1[nh][jlev] * sim.grid_arrays_dim.dcc[nh][jlev];
                zrb1s[nh][jlev] =
                    zrcsu[nh][jlev] + zrcl1s[nh][jlev] * sim.grid_arrays_dim.dcc[nh][jlev];

                // b) T
                //
                // ozon absorption
                //
                // downward beam
                zo3[nh] = zxo3l[nh][jlev];
                zto3[nh] = (1.
                    - (0.02118 * zo3[nh] / (1. + 0.042 * zo3[nh] + 0.000323 * zo3[nh].powf(2.))
                        + 1.082 * zo3[nh] / ((1. + 138.6 * zo3[nh]).powf(0.805))
                        + 0.0658 * zo3[nh] / (1. + (103.6 * zo3[nh]).powf(3.)))
                        / zsolar1)
                    / zto3t[nh];
                zto3t[nh] = zto3t[nh] * zto3[nh];

                // upward scattered beam
                zo3[nh] = zxo3t[nh] + zmbar * (zo3t[nh] - zo3l[nh][jlev]);
                zto3u[nh] = zto3tu[nh]
                    / (1.
                        - (0.02118 * zo3[nh]
                            / (1. + 0.042 * zo3[nh] + 0.000323 * zo3[nh].powf(2.))
                            + 1.082 * zo3[nh] / ((1. + 138.6 * zo3[nh]).powf(0.805))
                            + 0.0658 * zo3[nh] / 1.
                            + (103.6 * zo3[nh]).powf(3.))
                            / zsolar1);
                zto3tu[nh] = zto3tu[nh] / zto3u[nh];

                // total T = 1-(A(ozon)+R(rayl.))*(1-dcc)-R(cloud)*dcc
                ztb1[nh][jlev] = 1.
                    - (1. - zto3[nh]) * (1. - sim.grid_arrays_dim.dcc[nh][jlev])
                    - zrb1[nh][jlev];
                ztb1u[nh][jlev] = 1.
                    - (1. - zto3u[nh]) * (1. - sim.grid_arrays_dim.dcc[nh][jlev])
                    - zrb1s[nh][jlev];

                // make combined layer R_ab, R_abs, T_ab and T_abs
                z1mrabr[nh] = 1. / (1. - zra1s[nh] * zrb1s[nh][jlev]);
                zra1[nh] = zra1[nh] + zta1[nh] * zrb1[nh][jlev] * zta1s[nh] * z1mrabr[nh];
                zta1[nh] = zta1[nh] * ztb1[nh][jlev] * z1mrabr[nh];
                zra1s[nh] =
                    zrb1s[nh][jlev] + ztb1u[nh][jlev] * zra1s[nh] * ztb1[nh][jlev] * z1mrabr[nh];
                zta1s[nh] = ztb1u[nh][jlev] * zta1s[nh] * z1mrabr[nh];

                // 2. spectral range 2:
                //
                // a) R
                //
                //  cloud albedo
                zrb2[nh][jlev] = zrcl2[nh][jlev] * sim.grid_arrays_dim.dcc[nh][jlev];
                zrb2s[nh][jlev] = zrcl2s[nh][jlev] * sim.grid_arrays_dim.dcc[nh][jlev];

                // b) T
                //
                // water vapor absorption
                //
                // downward beam
                zwv[nh] = zywvl[nh][jlev];
                ztwv[nh] = (1.
                    - 2.9 * zwv[nh]
                        / ((1. + 141.5 * zwv[nh]).powf(0.635) + 5.925 * zwv[nh])
                        / zsolar2)
                    / ztwvt[nh];
                ztwvt[nh] = ztwvt[nh] * ztwv[nh];

                // upward scattered beam
                zwv[nh] = zywvt[nh] + zbetta * (zwvt[nh] - zwvl[nh][jlev]);
                ztwvu[nh] = ztwvtu[nh]
                    / (1.
                        - 2.9 * zwv[nh]
                            / ((1. + 141.5 * zwv[nh]).powf(0.635) + 5.925 * zwv[nh])
                            / zsolar2);
                ztwvtu[nh] = ztwvtu[nh] / ztwvu[nh];

                // total T = 1-A(water vapor)*(1.-dcc)-(A(cloud)+R(cloud))*dcc
                ztb2[nh][jlev] = 1.
                    - (1. - ztwv[nh]) * (1. - sim.grid_arrays_dim.dcc[nh][jlev])
                    - (1. - ztcl2[nh][jlev]) * sim.grid_arrays_dim.dcc[nh][jlev];
                ztb2u[nh][jlev] = 1.
                    - (1. - ztwvu[nh]) * (1. - sim.grid_arrays_dim.dcc[nh][jlev])
                    - (1. - ztcl2s[nh][jlev]) * sim.grid_arrays_dim.dcc[nh][jlev];

                // make combined layer R_ab, R_abs, T_ab and T_abs
                z1mrabr[nh] = 1. / (1. - zra2s[nh] * zrb2s[nh][jlev]);
                zra2[nh] = zra2[nh] + zta2[nh] * zrb2[nh][jlev] * zta2s[nh] * z1mrabr[nh];
                zta2[nh] = zta2[nh] * ztb2[nh][jlev] * z1mrabr[nh];
                zra2s[nh] =
                    zrb2s[nh][jlev] + ztb2u[nh][jlev] * zra2s[nh] * ztb2[nh][jlev] * z1mrabr[nh];
                zta2s[nh] = ztb2u[nh][jlev] * zta2s[nh] * z1mrabr[nh];
            }
        }
    }

    for nh in 0..NHOR {
        if losun[nh] {
            zt1[nh][NLEP - 1] = zta1[nh];
            zt2[nh][NLEP - 1] = zta2[nh];
            zr1s[nh][NLEP - 1] = zra1s[nh];
            zr2s[nh][NLEP - 1] = zra2s[nh];

            // upward loop
            //
            // make upward R
            zra1s[nh] = sim.rad_basic.dalb[nh];
            zra2s[nh] = sim.rad_basic.dalb[nh];

            // set albedo for the direct beam (for ocean use ECHAM3 param)
            sim.rad_basic.dalb[nh] = sim.surface.dls[nh] * sim.rad_basic.dalb[nh]
                + (1. - sim.surface.dls[nh]) * sim.surface.dicec[nh] * sim.rad_basic.dalb[nh]
                + (1. - sim.surface.dls[nh])
                    * (1. - sim.surface.dicec[nh])
                    * min(0.05 / (zmu0[nh] + 0.15), 0.15);
            zra1[nh] = sim.rad_basic.dalb[nh];
            zra2[nh] = sim.rad_basic.dalb[nh];
        }
    }

    for jlev in (0..NLEV - 1).rev() {
        for nh in 0..NHOR {
            if losun[nh] {
                zrl1[nh][jlev + 1] = zra1[nh];
                zrl2[nh][jlev + 1] = zra2[nh];
                zrl1s[nh][jlev + 1] = zra1s[nh];
                zrl2s[nh][jlev + 1] = zra2s[nh];
                zra1[nh] = zrb1[nh][jlev]
                    + ztb1[nh][jlev] * zra1[nh] * ztb1u[nh][jlev]
                        / (1. - zra1s[nh] * zrb1s[nh][jlev]);
                zra1s[nh] = zrb1s[nh][jlev]
                    + ztb1u[nh][jlev] * zra1s[nh] * ztb1u[nh][jlev]
                        / (1. - zra1s[nh] * zrb1s[nh][jlev]);
                zra2[nh] = zrb2[nh][jlev]
                    + ztb2[nh][jlev] * zra2[nh] * ztb2u[nh][jlev]
                        / (1. - zra2s[nh] * zrb2s[nh][jlev]);
                zra2s[nh] = zrb2s[nh][jlev]
                    + ztb2u[nh][jlev] * zra2s[nh] * ztb2u[nh][jlev]
                        / (1. - zra2s[nh] * zrb2s[nh][jlev]);
            }
        }
    }
    for nh in 0..NHOR {
        if losun[nh] {
            zrl1[nh][0] = zra1[nh];
            zrl2[nh][0] = zra2[nh];
            zrl1s[nh][0] = zra1s[nh];
            zrl2s[nh][0] = zra2s[nh];
        }
    }

    // fluxes at layer interfaces

    for jlev in 0..NLEP {
        for nh in 0..NHOR {
            if losun[nh] {
                z1mrabr[nh] = 1. / (1. - zr1s[nh][jlev] * zrl1s[nh][jlev]);
                zfd1[nh] = zt1[nh][jlev] * z1mrabr[nh];
                zfu1[nh] = -zt1[nh][jlev] * zrl1[nh][jlev] * z1mrabr[nh];
                z1mrabr[nh] = 1. / (1. - zr2s[nh][jlev] * zrl2s[nh][jlev]);
                zfd2[nh] = zt2[nh][jlev] * z1mrabr[nh];
                zfu2[nh] = -zt2[nh][jlev] * zrl2[nh][jlev] * z1mrabr[nh];
                sim.rad_basic.dfu[nh][jlev] = zfu1[nh] * zftop1[nh] + zfu2[nh] * zftop2[nh];
                sim.rad_basic.dfd[nh][jlev] = zfd1[nh] * zftop1[nh] + zfd2[nh] * zftop2[nh];
                sim.rad_basic.dswfl[nh][jlev] =
                    sim.rad_basic.dfu[nh][jlev] + sim.rad_basic.dfd[nh][jlev];
            }
        }
    }
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
fn lwr(mut sim: &mut Sim) {
    //
    // 0) define local parameters
    //

    const zmmair: FloatNum = 0.0289644; // molecular weight air (kg/mol)
    const zmmco2: FloatNum = 0.0440098; // molecular weight co2 (kg/mol)
    const zpv2pm: FloatNum = zmmco2 / zmmair; // transfere co2 ppvol to ppmass
    const zrco2: FloatNum = 1.9635; // co2 density (kg/m3 stp)
    const zro3: FloatNum = 2.14; // o3 density (kg/m3 stp)
    const zttop: FloatNum = 0.; // t at top of atmosphere

    // scaling factors for transmissivities
    // uh2o in g/cm**2 (zfh2o=0.1=1000/(100*100) *kg/m**2)
    // uco2 in cm-STP  (zfco2=100./rco2 * kg/m**2)
    // uo3  in cm-STP  (zfo3=100./ro3 * kg/m**2)
    const zfh2o: FloatNum = 0.1;
    const zfco2: FloatNum = 100. / zrco2;
    const zfo3: FloatNum = 100. / zro3;
    const zt0: FloatNum = 295.;

    // local arrays
    // TODO make sure how this syntax works: "array(0:VAR)"
    let mut zbu = vec![vec![0.; NLEP + 1]; NHOR]; // effective SBK*T**4 for upward radiation (NHOR,0:NLEP)
    let mut zbd = vec![vec![0.; NLEP + 1]; NHOR]; // effective SBK*T**4 for downward radiation (NHOR,0:NLEP)
    let mut zst4h = vec![vec![0.; NLEP]; NHOR]; // SBK*T**4  on half levels (NHOR,NLEP)
    let mut zst4 = vec![vec![0.; NLEP]; NHOR]; // SBK*T**4  on ull levels (NHOR,NLEP)
    let mut ztau = vec![vec![0.; NLEV]; NHOR]; // total transmissivity (NHOR,NLEV)
    let mut zq = vec![vec![0.; NLEV]; NHOR]; // modified water vapour (NHOR,NLEV)
    let mut zqo3 = vec![vec![0.; NLEV]; NHOR]; // modified ozon (NHOR,NLEV)
    let mut zqco2 = vec![vec![0.; NLEV]; NHOR]; // modified co2 (NHOR,NLEV)
    let mut ztausf = vec![vec![0.; NLEV]; NHOR]; // total transmissivity to surface (NHOR,NLEV)
    let mut ztaucs = vec![vec![0.; NLEV]; NHOR]; // clear sky transmissivity (NHOR,NLEV)
    let mut ztaucc0 = vec![vec![0.; NLEV]; NHOR]; // layer transmissivity cloud (NHOR,NLEV)

    let mut ztaucc = vec![0.; NHOR]; // cloud transmissivity (NHOR)
    let mut ztau0 = vec![0.; NHOR]; // approx. layer transmissivity (NHOR)
    let mut zsumwv = vec![0.; NHOR]; // effective water vapor amount (NHOR)
    let mut zsumo3 = vec![0.; NHOR]; // effective o3 amount (NHOR)
    let mut zsumco2 = vec![0.; NHOR]; // effective co2 amount (NHOR)
    let mut zsfac = vec![0.; NHOR]; // scaling factor (NHOR)
    let mut zah2o = vec![0.; NHOR]; // water vapor absorptivity (NHOR)
    let mut zaco2 = vec![0.; NHOR]; // co2 absorptivity (NHOR)
    let mut zao3 = vec![0.; NHOR]; // o3 absorptivity (NHOR)
    let mut zth2o = vec![0.; NHOR]; // water vapor - co2 overlap transmissivity (NHOR)
    let mut zbdl = vec![0.; NHOR]; // layer evective downward rad. (NHOR)
    let mut zeps = vec![0.; NHOR]; // surface emissivity (NHOR)
    let mut zps2 = vec![0.; NHOR]; // ps**2 (NHOR)
    let mut zsigh2 = vec![0.; NLEP]; // sigmah**2 (NLEP)

    // entropy
    let mut zbue1 = vec![vec![0.; NLEP + 1]; NHOR]; // (NHOR,0:NLEP)
    let mut zbue2 = vec![vec![0.; NLEP + 1]; NHOR]; // (NHOR,0:NLEP)

    //
    // 1) set some necessary (helpful) bits
    //

    let zero = 0.000001; // small number

    let zao30 = 0.209 * (0.00007 as FloatNum).powf(0.436); // to get a(o3)=0 for o3=0
    let zco20 = 0.0676 * (0.01022 as FloatNum).powf(0.421); // to get a(co2)=0 for co2=0
    let zh2o0a = 0.846 * (0.0000359 as FloatNum).powf(0.243); // to get a(h2o)=0 for h2o=0
    let zh2o0 = 0.832 * ((0.0286 as FloatNum).powf(0.26)); // to get t(h2o)=1 for h2o=0

    // to make a(o3) continues at 0.01cm:
    let zao3c = 0.209 * ((0.01 + 0.00007) as FloatNum).powf(0.436)
        - zao30
        - 0.0212 * (0.01 as FloatNum).log10();

    // to make a(co2) continues at 1cm:
    let zaco2c = 0.0676 * (1.01022 as FloatNum).powf(0.421) - zco20;

    // to make a(h2o) continues at 0.01gm:
    let zah2oc = 0.846 * (((0.01 + 3.59E-5) as FloatNum).powf(0.243))
        - zh2o0a
        - 0.24 * (0.02 as FloatNum).log10();

    // to make t(h2o) continues at 2gm :
    let zth2oc = 1. - (0.832 * (((2. + 0.0286) as FloatNum).powf(0.26)) - zh2o0)
        + 0.1196 * ((2. - 0.6931) as FloatNum).ln();

    zsigh2[0] = 0.;

    //zsigh2(2:NLEP) = sigmah(1: NLEV)**2;
    zsigh2
        .iter_mut()
        .skip(1)
        .enumerate()
        .map(|(i, n)| *n = sim.lev_arrays.sigmah[i - 1].powf(2.));

    // zps2(:)=dp(:)*dp(:)
    zps2.iter_mut()
        .enumerate()
        .map(|(i, n)| *n = sim.grid_arrays_dim.dp[i] * sim.grid_arrays_dim.dp[i]);

    // 2) calc. stb*t**4 + preset fluxes

    // TODO this has already been initialized elsewhere
    sim.rad_basic.dftu = vec![vec![0.; NLEP]; NHOR];
    sim.rad_basic.dftd = vec![vec![0.; NLEP]; NHOR];

    // stb*t**4 on full and half levels

    // full levels and surface

    // zst4(:,1:NLEP)=SBK*dt(:,1:NLEP)**4;
    for (i, n) in zst4.iter_mut().enumerate() {
        for (j, m) in n.iter_mut().enumerate() {
            *m = SBK * sim.grid_arrays_dim.dt[i][j].powf(4.);
        }
    }

    // half level (incl. toa and near surface)

    // zst4h(:,1)=zst4(:,1)-(zst4(:,1)-zst4(:,2))*sigma(1)/(sigma(1)-sigma(2));
    zst4h.iter_mut().enumerate().map(|(i, n)| {
        n[0] = zst4[i][0]
            - (zst4[i][0] - zst4[i][1]) * sim.lev_arrays.sigma[0]
            - sim.lev_arrays.sigma[1]
    });

    for jlev in 1..NLEV {
        let jlem = jlev - 1;
        zst4h.iter_mut().enumerate().map(|(i, n)| {
            n[jlev] = (zst4[i][jlev] * (sim.lev_arrays.sigma[jlem] - sim.lev_arrays.sigmah[jlem])
                + zst4[i][jlem] * (sim.lev_arrays.sigmah[jlem])
                - (sim.lev_arrays.sigma[jlev]))
                / (sim.lev_arrays.sigma[jlem] - sim.lev_arrays.sigma[jlev])
        });
    }

    for nh in 0..NHOR {
        if (zst4[nh][NLEV - 1] - zst4h[nh][NLEV - 1]) * (zst4[nh][NLEV - 1] - zst4[nh][NLEP - 1])
            > 0.
        {
            zst4h[nh][NLEP - 1] = zst4[nh][NLEP - 1];
        } else {
            zst4h[nh][NLEP - 1] = zst4[nh][NLEV - 1]
                + (zst4[nh][NLEV - 1] - zst4[nh][NLEM - 1]) * (1. - sim.lev_arrays.sigma[NLEV - 1])
                    / (sim.lev_arrays.sigma[NLEV - 1])
                - sim.lev_arrays.sigma[NLEM - 1];
        }
    }

    // top downward flux, surface grayness and surface upward flux

    zbd.iter_mut().map(|n| n[0] = SBK * zttop.powf(4.));
    zeps.iter_mut()
        .enumerate()
        .map(|(i, n)| *n = sim.surface.dls[i] + 0.98 * (1. - sim.surface.dls[i]));
    zbu.iter_mut()
        .enumerate()
        .map(|(i, n)| n[NLEP - 1] = zeps[i] * zst4[i][NLEP - 1]);
    zbue1.iter_mut().map(|n| n[NLEP - 1] = 0.);
    zbue2
        .iter_mut()
        .enumerate()
        .map(|(i, n)| n[NLEP - 1] = zbu[i][NLEP - 1]);

    // 3) vertical loops
    //
    // a) modified absorber amounts and cloud transmissivities
    //

    for jlev in 0..NLEV {
        let jlep = jlev + 1;
        let zzf1 =
            sim.lev_arrays.sigma[jlev] * sim.lev_arrays.dsigma[jlev] / sim.planet_vars.ga / 100000.;
        let zzf2 = zpv2pm * 1E-6;
        let zzf3 =
            -1.66 * sim.rad.acllwr * 1000. * sim.lev_arrays.dsigma[jlev] / sim.planet_vars.ga;
        zsfac
            .iter_mut()
            .enumerate()
            .map(|(i, n)| *n = zzf1 * zps2[i]);
        zq.iter_mut()
            .enumerate()
            .map(|(i, n)| n[jlev] = zfh2o * zsfac[i] * sim.grid_arrays_dim.dq[i][jlev]);
        zqo3.iter_mut()
            .enumerate()
            .map(|(i, n)| n[jlev] = zfo3 * zsfac[i] * sim.rad.dqo3[i][jlev]);
        zqco2
            .iter_mut()
            .enumerate()
            .map(|(i, n)| n[jlev] = zfco2 * zzf2 * zsfac[i] * sim.rad.dqco2[i][jlev]);

        if sim.rad.clgray > 0. {
            ztaucc0
                .iter_mut()
                .enumerate()
                .map(|(i, n)| n[jlev] = 1. - sim.grid_arrays_dim.dcc[i][jlev] * sim.rad.clgray);
        } else {
            ztaucc0.iter_mut().enumerate().map(|(i, n)| {
                n[jlev] = 1.
                    - sim.grid_arrays_dim.dcc[i][jlev]
                        * (1.
                            - (zzf3 * sim.grid_arrays_dim.dql[i][jlev] * sim.grid_arrays_dim.dp[i])
                                .exp())
            });
        }
    }

    // b) transmissivities, effective radiations and fluxes

    for jlev in 0..NLEV {
        // TODO how to handle these weirdly -1 indexed arrays
        // (there are two as far as I can tell: zbu, zbd)
        // let jlem = jlev as Int - 1;
        let jlem = jlev;
        ztaucc = vec![1.; NHOR];
        zsumwv = vec![0.; NHOR];
        zsumo3 = vec![0.; NHOR];
        zsumco2 = vec![0.; NHOR];

        // transmissivities
        for jlev2 in jlev..NLEV {
            let jlep2 = jlev2 + 1;
            for (i, n) in zsumwv.iter_mut().enumerate() {
                *n = *n + zq[i][jlev2];
            }
            for (i, n) in zsumo3.iter_mut().enumerate() {
                *n = *n + zqo3[i][jlev2];
            }
            for (i, n) in zsumco2.iter_mut().enumerate() {
                *n = *n + zqco2[i][jlev2];
            }

            // clear sky transmisivity
            //
            // h2o absorption:
            //
            // a) 6.3mu
            //
            for nh in 0..NHOR {
                if zsumwv[nh] <= 0.01 {
                    zah2o[nh] = 0.846 * (zsumwv[nh] + 3.59E-5).powf(0.243) - zh2o0a;
                } else {
                    zah2o[nh] = 0.24 * (zsumwv[nh] + 0.01).log10() + zah2oc;
                }
            }

            // b) continuum
            //
            if sim.rad.th2oc > 0. {
                for (i, n) in zah2o.iter_mut().enumerate() {
                    *n = min(*n + (1. - (-sim.rad.th2oc * zsumwv[i]).exp()), 1.)
                }
            }

            // co2 absorption:
            //
            for nh in 0..NHOR {
                if zsumco2[nh] <= 1.0 {
                    zaco2[nh] = 0.0676 * (zsumco2[nh] + 0.01022).powf(0.421) - zco20;
                } else {
                    zaco2[nh] = 0.0546 * (zsumco2[nh]).log10() + zaco2c;
                }
            }

            //  Boer et al. (1984) scheme for t(h2o) at co2 overlapp
            //
            for nh in 0..NHOR {
                if zsumwv[nh] <= 2. {
                    zth2o[nh] = 1. - (0.832 * (zsumwv[nh] + 0.0286).powf(0.26) - zh2o0);
                } else {
                    zth2o[nh] = max(0., zth2oc - 0.1196 * (zsumwv[nh] - 0.6931).ln());
                }
            }

            // o3 absorption:
            //
            for nh in 0..NHOR {
                if zsumo3[nh] <= 0.01 {
                    zao3[nh] = 0.209 * (zsumo3[nh] + 7E-5).powf(0.436) - zao30;
                } else {
                    zao3[nh] = 0.0212 * (zsumo3[nh]).log10() + zao3c;
                }
            }

            // total clear sky transmissivity
            for (i, n) in ztaucs.iter_mut().enumerate() {
                n[jlev2] = 1. - zah2o[i] - zao3[i] - zaco2[i] * zth2o[i];
            }

            // bound transmissivity:
            for (i, n) in ztaucs.iter_mut().enumerate() {
                n[jlev2] = min(1. - zero, max(zero, n[jlev2]));
            }

            // cloud transmisivity assuming random overlap
            for (i, n) in ztaucc.iter_mut().enumerate() {
                *n = *n * ztaucc0[i][jlev2];
            }

            // total transmissivity
            for (i, n) in ztau.iter_mut().enumerate() {
                n[jlev2] = ztaucs[i][jlev2] * ztaucc[i];
            }
        }

        // upward and downward effective SBK*T**4
        //

        if jlev == 1 {
            for n in ztau0.iter_mut() {
                *n = 1. - zero;
            }

            for jlev2 in 0..NLEV {
                let jlep2 = jlev2 + 1;
                for (i, n) in ztau0.iter_mut().enumerate() {
                    *n = ztaucs[i][jlev2] / *n * ztaucc0[i][jlev2] * sim.rad.tpofmt;
                }
                for (i, n) in ztau0.iter_mut().enumerate() {
                    *n = min(1. - zero, max(zero, *n));
                }

                for nh in 0..NHOR {
                    if (zst4[nh][jlev2] - zst4h[nh][jlev2]) * (zst4[nh][jlev2] - zst4h[nh][jlep2])
                        > 0.
                    {
                        zbd[nh][jlev2] =
                            0.5 * zst4[nh][jlev2] + 0.25 * (zst4h[nh][jlev2] + zst4h[nh][jlep2]);
                        zbu[nh][jlev2] = zbd[nh][jlev2];
                        zbue1[nh][jlev2] = zbu[nh][jlev2];
                        zbue2[nh][jlev2] = 0.;
                    } else {
                        zbd[nh][jlev2] = (zst4h[nh][jlep2] - ztau0[nh] * zst4h[nh][jlev2])
                            / (1. - ztau0[nh])
                            - (zst4h[nh][jlev2] - zst4h[nh][jlep2]) / ztau0[nh].ln();
                        zbu[nh][jlev2] = zst4h[nh][jlev2] + zst4h[nh][jlep2] - zbd[nh][jlev2];
                        zbue1[nh][jlev2] = zbu[nh][jlev2];
                        zbue2[nh][jlev2] = 0.;
                    }
                }

                for (i, n) in ztau0.iter_mut().enumerate() {
                    *n = ztaucs[i][jlev2];
                }
            }
        }

        // fluxes
        //

        for (i, n) in sim.rad_basic.dftu.iter_mut().enumerate() {
            n[jlev] = n[jlev] - zbu[i][jlev];
        }
        for (i, n) in sim.rad_basic.dftd.iter_mut().enumerate() {
            n[jlev] = n[jlev] + zbd[i][jlem];
        }
        for (i, n) in sim.rad.dftue1.iter_mut().enumerate() {
            n[jlev] = n[jlev] - zbue1[i][jlev];
        }
        for (i, n) in sim.rad.dftue2.iter_mut().enumerate() {
            n[jlev] = n[jlev] - zbue2[i][jlev];
        }
        for (i, n) in zbdl.iter_mut().enumerate() {
            *n = zbd[i][jlem] - zbd[i][jlev];
        }

        for jlev2 in jlev..NLEV {
            let jlep2 = jlev2 + 1;

            for (i, n) in sim.rad_basic.dftu.iter_mut().enumerate() {
                n[jlev] = n[jlev] - (zbu[i][jlep2] - zbu[i][jlev2]) * ztau[i][jlev2];
            }
            for (i, n) in sim.rad.dftue1.iter_mut().enumerate() {
                n[jlev] = n[jlev] - (zbue1[i][jlep2] - zbue1[i][jlev2]) * ztau[i][jlev2];
            }
            for (i, n) in sim.rad.dftue2.iter_mut().enumerate() {
                n[jlev] = n[jlev] - (zbue2[i][jlep2] - zbue2[i][jlev2]) * ztau[i][jlev2];
            }
            for (i, n) in sim.rad_basic.dftd.iter_mut().enumerate() {
                n[jlep2] = n[jlep2] + zbdl[i] * ztau[i][jlev2];
            }
        }

        if jlev == 1 {
            for (i, n) in sim.rad.dftu0.iter_mut().enumerate() {
                n[0] = -zbu[i][0] * (1. - ztau[i][0]);
            }
            for jlev2 in 1..NLEV {
                let jlem = jlev2 - 1;
                for (i, n) in sim.rad.dftu0.iter_mut().enumerate() {
                    n[jlev2] = -zbu[i][jlev2] * (ztau[i][jlem] - ztau[i][jlev2]);
                }
            }
        }

        // collect transmissivity to surface
        for (i, n) in ztausf.iter_mut().enumerate() {
            n[jlev] = ztau[i][NLEV - 1];
        }
    }

    for jlev in 0..NLEV - 1 {
        let jlep = jlev + 1;
        for (i, n) in sim.rad.dftd0.iter_mut().enumerate() {
            n[jlev] = -zbd[i][jlev] * (ztausf[i][jlep] - ztausf[i][jlev]);
        }
    }

    // TODO jlev is not available here, source code faulty?
    // for (i, n) in sim.rad.dftd0.iter_mut().enumerate() {
    //     n[NLEV - 1] = -zbd[i][jlev] * (1. - ztausf[i][NLEV - 1]);
    // }

    // complite surface lwr
    //
    for (i, n) in sim.rad_basic.dftu.iter_mut().enumerate() {
        n[NLEP - 1] = n[NLEP - 1] - zbu[i][NLEP - 1];
    }
    for (i, n) in sim.rad_basic.dftd.iter_mut().enumerate() {
        n[NLEP - 1] = n[NLEP - 1] + zbd[i][NLEV - 1];
    }
    for (i, n) in sim.rad.dftue1.iter_mut().enumerate() {
        n[NLEP - 1] = n[NLEP - 1] - zbue1[i][NLEP - 1];
    }
    for (i, n) in sim.rad.dftue2.iter_mut().enumerate() {
        n[NLEP - 1] = n[NLEP - 1] - zbue2[i][NLEP - 1];
    }

    // correct for non black suface
    //
    for (i, n) in zeps.iter_mut().enumerate() {
        *n = (1. - *n) * sim.rad_basic.dftd[i][NLEP - 1];
    }
    for jlev in 0..NLEV {
        for (i, n) in sim.rad_basic.dftu.iter_mut().enumerate() {
            n[jlev] = n[jlev] - ztausf[i][jlev] * zeps[i];
        }
        for (i, n) in sim.rad.dftue2.iter_mut().enumerate() {
            n[jlev] = n[jlev] - ztausf[i][jlev] * zeps[i];
        }
    }

    // total longwave radiation
    for (i, n) in sim.rad_basic.dlwfl.iter_mut().enumerate() {
        for (j, m) in n.iter_mut().enumerate() {
            *m = sim.rad_basic.dftu[i][j] + sim.rad_basic.dftd[i][j];
        }
    }
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
    eccen: FloatNum,
    obliq: FloatNum,
    mvelp: FloatNum,
    log_print: bool,
    mut obliqr: &mut FloatNum,
    mut lambm0: &mut FloatNum,
    mut mvelpp: &mut FloatNum,
) {
    // Parameters for calculating earth's orbital characteristics

    const POBLEN: usize = 47; // number of elements in the series to calc obliquity
    const PECCLEN: usize = 19; // number of elements in the series to calc eccentricity
    const PMVELEN: usize = 78; // number of elements in the series to calc vernal equinox

    let degrad: FloatNum; // degrees to radians conversion factor
    let obamp: [FloatNum; POBLEN]; // amplitudes for obliquity cosine series
    let obrate: [FloatNum; POBLEN]; // rates for obliquity cosine series
    let obphas: [FloatNum; POBLEN]; // phases for obliquity cosine series
    let ecamp: [FloatNum; PECCLEN]; // amplitudes for eccentricity/fvelp cosine/sine series
    let ecrate: [FloatNum; PECCLEN]; // rates for eccentricity/fvelp cosine/sine series
    let ecphas: [FloatNum; PECCLEN]; // phases for eccentricity/fvelp cosine/sine series
    let mvamp: [FloatNum; PMVELEN]; // amplitudes for mvelp sine series
    let mvrate: [FloatNum; PMVELEN]; // rates for mvelp sine series
    let mvphas: [FloatNum; PMVELEN]; // phases for mvelp sine series
    let mut yb4_1950AD: FloatNum = 0.; // number of years before 1950 AD

    let psecdeg: FloatNum = 1. / 3600.; // arc seconds to degrees conversion

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
    let mut obsum: FloatNum; // Obliquity series summation
    let mut cossum: FloatNum; // Cosine series summation for eccentricity/fvelp
    let mut sinsum: FloatNum; // Sine series summation for eccentricity/fvelp
    let mut fvelp: FloatNum = 0.0; // Fixed vernal equinox longitude of perihelion
    let mut mvsum: FloatNum; // mvelp series summation
    let beta: FloatNum; // Intermediate argument for lambm0
    let years: FloatNum; // Years to time of interest (negative = past; positive = future)

    let mut eccen2: FloatNum; // eccentricity squared
    let mut eccen3: FloatNum; // eccentricity cubed
    let pi: FloatNum = 4. * (1. as FloatNum).atan() as FloatNum; // pi

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
        yb4_1950AD = 1950.0 - iyear_AD as FloatNum;
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
    calday: FloatNum,
    eccen: FloatNum,
    mvelpp: FloatNum,
    lambm0: FloatNum,
    obliqr: FloatNum,
    n_days_per_year: Int,
    ndatim: Vec<Int>,
    mut delta: &mut FloatNum,
    mut eccf: &mut FloatNum,
) {
    // Calday of vernal equinox, correct for Jan 1 = calday 1
    let ve: FloatNum = 80.5;
    let pie = std::f64::consts::PI as FloatNum;

    // Lambda m, earth's mean longitude of perihelion (radians)
    let mut lambm: FloatNum;
    // Intermediate argument involving lambm
    let mut lmm: FloatNum;
    // Lambda, the earth's longitude of perihelion
    let lamb: FloatNum;
    // Inverse normalized sun/earth distance
    let mut invrho: FloatNum;
    // Sine of lmm
    let sinl: FloatNum;

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

    lambm = lambm0 + (calday - ve) * 2. * pie / (n_days_per_year + ndatim[6]) as FloatNum; // ndatim(7) = leap year
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
fn orb_print(iyear_AD: Int, eccen: FloatNum, obliq: FloatNum, mvelp: FloatNum) {
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
