use crate::{FloatNum, Int};
use std::ops::Rem;

#[derive(Clone)]
pub struct Calendar {
    pub mondays: [Int; 13],
    pub mona365: [Int; 13],
    pub monaccu: [Int; 13],

    pub ny400d: Int,
    pub ny100d: Int,
    pub ny004d: Int,
    pub ny001d: Int,
    pub nud: Int,

    // TODO ?
    // These values are copied from pumamod in subroutine calini
    pub n_days_per_month: Int,
    pub n_days_per_year: Int,
    pub n_start_step: Int,
    pub ntspd: Int,
    pub solar_day: FloatNum,
}

impl Default for Calendar {
    fn default() -> Self {
        Self {
            mondays: [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
            mona365: [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365],
            monaccu: [0; 13],
            ny400d: 400 * 365 + 97,
            ny100d: 100 * 365 + 24,
            ny004d: 4 * 365 + 1,
            ny001d: 365,
            nud: 6,
            n_days_per_month: 30,
            n_days_per_year: 360,
            n_start_step: 0,
            ntspd: 1,
            solar_day: 86400.0, // sec
        }
    }
}

fn yday2mmdd(cal: &Calendar, mut kyday: &mut Int, mut kmon: &mut Int, mut kday: &mut Int) {
    if cal.n_days_per_year == 365 {
        *kmon = 1;
        while *kyday > cal.mona365[*kmon as usize] {
            *kmon = *kmon + 1;
        }
        *kday = *kyday - cal.monaccu[*kmon as usize - 1];
    } else {
        *kmon = (*kyday - 1) / cal.n_days_per_month;
        *kday = *kyday - cal.n_days_per_month * *kmon;
    }
}

fn nweekday(kday: Int) -> Int {
    (kday + 5).rem(7) as Int
}

pub fn ndayofyear(kstep: Int, cal: &mut Calendar) -> Int {
    let mut idatim = [0; 7];

    if cal.n_days_per_year == 365 {
        step2cal(kstep, cal.ntspd, &mut idatim, cal);
        return idatim[2] + cal.monaccu[idatim[1] as usize - 1];
    } else {
        step2cal30(kstep, cal.ntspd, &mut idatim, cal);
        return idatim[2] + cal.n_days_per_month * (idatim[1] - 1);
    }

    0
}

//  kstep     ! time step since simulation start
//  ktspd     ! time steps per day
//  kdatim(7) ! year,month,day,hour,min,weekday,leapyear
fn step2cal(kstep: Int, ktspd: Int, kdatim: &mut [Int; 7], mut cal: &mut Calendar) {
    let mut iyea: Int; // current year of simulation
    let mut imon: Int; // current month of simulation
    let mut iday: Int; // current day of simulation
    let mut ihou: Int; // current hour of simulation
    let mut imin: Int; // current minute of simulation

    let mut idall: Int;
    let istp: Int;

    let iy400: Int;
    let id400: Int;
    let iy100: Int;
    let id100: Int;
    let iy004: Int;
    let id004: Int;
    let iy001: Int;
    let mut id001: Int;
    let jmon: Int;

    let leap: bool;

    idall = kstep / ktspd;

    iy400 = idall / cal.ny400d; // segment (of 400 years)
    id400 = idall.rem(cal.ny400d);

    if id400 <= cal.ny100d {
        // century year is leap year
        iy100 = 0; // century in segment [0]
        id100 = id400;
        iy004 = id100 / cal.ny004d; // tetrade in century [0..24]
        id004 = id100.rem(cal.ny004d);
        leap = id004 <= cal.ny001d;
        if leap {
            iy001 = 0; // year in tetrade [0]
            id001 = id004;
        } else {
            iy001 = (id004 - 1) / cal.ny001d; // year in tetrade [1,2,3]
            id001 = (id004 - 1).rem(cal.ny001d);
        }
    } else {
        // century year is not leap year
        iy100 = (id400 - 1) / cal.ny100d; // century in segment [1,2,3]
        id100 = (id400 - 1).rem(cal.ny100d);

        if id100 < cal.ny004d - 1 {
            iy004 = 0; // tetrade in century [0]
            id004 = id100;
            leap = false;
            iy001 = id004 / cal.ny001d; // year in tetrade [1,2,3]
            id001 = id004.rem(cal.ny001d);
        } else {
            iy004 = (id100 + 1) / cal.ny004d; // tetrade in century [1..24]
            id004 = (id100 + 1).rem(cal.ny004d);
            leap = id004 <= cal.ny001d;
            if leap {
                iy001 = 0; // year in tetrade [0]
                id001 = id004;
            } else {
                iy001 = (id004 - 1) / cal.ny001d;
                id001 = (id004 - 1).rem(cal.ny001d);
            }
        }
    }

    iyea = iy400 * 400 + iy100 * 100 + iy004 + 4 + iy001;

    cal.monaccu[0] = cal.mondays[0];
    cal.monaccu[1] = cal.mondays[1];
    cal.monaccu[2] = cal.mondays[1] + cal.mondays[2];
    if leap {
        cal.monaccu[2] = cal.monaccu[2] + 1;
        for jmon in 3..13 {
            cal.monaccu[jmon] = cal.monaccu[jmon - 1] + cal.mondays[jmon];
        }
        imon = 1;
        id001 = id001 + 1;
        while id001 > cal.monaccu[imon as usize] {
            imon = imon + 1;
        }
        iday = id001 - cal.monaccu[imon as usize - 1];

        istp = kstep.rem(ktspd);
        imin = (istp * 1440) / ktspd;
        ihou = imin / 60;
        imin = imin.rem(60);

        kdatim[0] = iyea;
        kdatim[1] = imon;
        kdatim[2] = iday;
        kdatim[3] = ihou;
        kdatim[4] = imin;
        kdatim[5] = (kstep / ktspd + 5).rem(7); // day of week
        if leap {
            kdatim[6] = 1;
        } else {
            kdatim[6] = 0;
        }
    }
}

fn cal2step(
    ktspd: Int,          // time steps per day
    kyea: Int,           // current year of simulation
    kmon: Int,           // current month of simulation
    kday: Int,           // current day of simulation
    khou: Int,           // current hour of simulation
    kmin: Int,           // current minute of simulation
    mut kstep: &mut Int, //time step since simulation start
    mut cal: &mut Calendar,
) {
    let mut idall: Int;
    let mut ilp: Int;

    let mut iy400: Int;
    let mut id400: Int;
    let mut iy100: Int;
    let mut id100: Int;
    let mut iy004: Int;
    let mut id004: Int;
    let mut jmon: Int;

    let leap: bool;

    // simplified calendar
    if cal.n_days_per_year != 365 {
        *kstep =
            ktspd * (kyea * cal.n_days_per_year + (kmon - 1) * cal.n_days_per_month + kday - 1);
        return;
    }

    iy400 = kyea / 400; // segment [400]
    id400 = kyea.rem(400); // year in segment [0..399]
    iy100 = id400 / 100; // century [0,1,2,3]
    id100 = id400.rem(100); // year in century [0..99]
    iy004 = id100 / 4; // tetrade [0..24]
    id004 = id100.rem(4); // year in tetrade [0,1,2,3]

    leap = (id004 == 0 && (id100 != 0 || id400 == 0));

    ilp = -1;
    if id004 > 0 {
        ilp = ilp + 1;
    }
    if iy100 > 0 && id100 == 0 {
        ilp = ilp + 1;
    }

    cal.monaccu[0] = cal.mondays[0];
    cal.monaccu[1] = cal.mondays[1];
    cal.monaccu[2] = cal.mondays[1];
    if leap {
        cal.monaccu[2] = cal.monaccu[2] + 1;
    }
    for jmon in 3..13 {
        cal.monaccu[jmon] = cal.monaccu[jmon - 1] + cal.mondays[jmon];
    }

    idall = iy400 * cal.ny400d
        + iy100 * cal.ny100d
        + iy004 * cal.ny004d
        + id004 * cal.ny001d
        + cal.monaccu[kmon as usize - 1]
        + kday
        + ilp;
    *kstep = ktspd * idall + (ktspd * (khou * 60 + kmin)) / 1440;
}

// kstep     ! time step since simulation start
// ktspd     ! time steps per day
// kdatim(7) ! year,month,day,hour,min,weekday,leapyear
fn step2cal30(kstep: Int, ktspd: Int, mut kdatim: &mut [Int; 7], cal: &mut Calendar) {
    let mut iyea: Int; // current year of simulation
    let mut imon: Int; // current month of simulation
    let mut iday: Int; // current day of simulation
    let mut ihou: Int; // current hour of simulation
    let mut imin: Int; // current minute of simulation
    let mut idall: Int;
    let mut istp: Int;

    idall = kstep / ktspd;
    iyea = idall / cal.n_days_per_year;
    idall = idall.rem(cal.n_days_per_year);
    imon = idall / cal.n_days_per_month + 1;
    iday = idall.rem(cal.n_days_per_month) + 1;
    istp = kstep.rem(ktspd);
    imin = ((istp as FloatNum * cal.solar_day) / (ktspd * 60) as FloatNum) as Int;
    ihou = imin / 60;
    imin = imin.rem(60);

    kdatim[0] = iyea;
    kdatim[1] = imon;
    kdatim[2] = iday;
    kdatim[3] = ihou;
    kdatim[4] = imin;
    kdatim[5] = 0; // day of week
    kdatim[6] = 0; // leap year
}

pub fn ntomin(
    kstep: Int,
    mut kmin: &mut Int,
    mut khou: &mut Int,
    mut kday: &mut Int,
    mut kmon: &mut Int,
    mut kyea: &mut Int,
    mut cal: &mut Calendar,
) {
    let mut idatim = [0; 7];

    if cal.n_days_per_year == 365 {
        step2cal(kstep, cal.ntspd, &mut idatim, cal);
    } else {
        step2cal30(kstep, cal.ntspd, &mut idatim, cal);
    }

    *kyea = idatim[0];
    *kmon = idatim[1];
    *kday = idatim[2];
    *khou = idatim[3];
    *kmin = idatim[4];
}

fn ntodat(istep: Int, datch: &mut String, cal: &mut Calendar) {
    let mona = vec![
        "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
    ];

    let mut imin: Int = 1;
    let mut ihou: Int = 1;
    let mut iday: Int = 1;
    let mut imon: Int = 1;
    let mut iyea: Int = 1;

    ntomin(
        istep, &mut imin, &mut ihou, &mut iday, &mut imon, &mut iyea, cal,
    );
    *datch = format!(
        "{}-{}-{} {}:{}",
        iday,
        mona[imon as usize - 1],
        iyea,
        ihou,
        imin
    )
}

// !     =================
// !     SUBROUTINE NTODAT
// !     =================
//
//       subroutine ntodat(istep,datch)
//       character(len=18) datch
//       character(len=3) mona(12)
//       data mona /'Jan','Feb','Mar','Apr','May','Jun',                   &
//      &           'Jul','Aug','Sep','Oct','Nov','Dec'/
//       call ntomin(istep,imin,ihou,iday,imon,iyea)
//       write (datch,20030) iday,mona(imon),iyea,ihou,imin
// 20030 format(i2,'-',a3,'-',i4.4,2x,i2,':',i2.2)
//       end
//
//
// !     =================
// !     SUBROUTINE DTODAT
// !     =================
//
//       subroutine dtodat(id,datch)
//       integer :: id(6)
//       character(len=18) datch
//       character(len=3) mona(12)
//       data mona /'Jan','Feb','Mar','Apr','May','Jun',                   &
//      &           'Jul','Aug','Sep','Oct','Nov','Dec'/
//       write (datch,20030) id(3),mona(id(2)),id(1),id(4),id(5)
// 20030 format(i2,'-',a3,'-',i4.4,2x,i2,':',i2.2)
//       end
//
//
// !     =================
// !     SUBROUTINE MOMINT
// !     =================
//
// !     Compute month indices and weights for time interpolation from
// !     monthly mean data to current timestep
//
//       subroutine momint(kperp,kstep,kmona,kmonb,pweight)
//       use calmod
//       implicit none
//       integer, intent(in ) :: kperp   ! perpetual mode ?
//       integer, intent(in ) :: kstep   ! target step
//       integer, intent(out) :: kmona   ! current month (1-12)
//       integer, intent(out) :: kmonb   ! next or previous month (0-13)
//       real   , intent(out) :: pweight ! interpolation weight
//
//       integer :: idatim(7) ! date time array
//       integer :: iday
//       integer :: ihour
//       integer :: imin
//       integer :: jmonb  ! next or previous month (1-12)
//
//       real    :: zday   ! fractional day (including hour & minute)
//       real    :: zdpma  ! days per month a
//       real    :: zdpmb  ! days per month b
//       real    :: zmeda  ! median day of the month a
//       real    :: zmedb  ! median day of the month b
//
// !     convert time step to date / time
//
//       idatim(:) = 0
//       if (kperp > 0) then                     ! perpetual date
//          call yday2mmdd(kperp,idatim(2),idatim(3))
//       else if (n_days_per_year == 365) then   ! real calendar
//          call step2cal(kstep,ntspd,idatim)
//       else                                    ! simple calendar
//          call step2cal30(kstep,ntspd,idatim)
//       endif
//
//       kmona = idatim(2)
//       iday  = idatim(3)
//       ihour = idatim(4)
//       imin  = idatim(5)
//
// !     set fractional day
//
//       zday = iday + ((ihour * 60.0 + imin) * 60.0) / solar_day
//
// !     compute median of month a
//
//       zdpma = n_days_per_month
//       if (n_days_per_year == 365) then
//          zdpma = mondays(kmona)
//          if (kmona == 2) zdpma = zdpma + idatim(7) ! leap year
//       endif
//       zmeda = 0.5 * (zdpma + 1.0) ! median day a
//
// !     define neighbour month
//
//       if (zday > zmeda) then
//          kmonb = kmona + 1 !     next month (maybe 13)
//       else
//          kmonb = kmona - 1 ! previous month (maybe  0)
//       endif
//
// !     compute median of month b
//
//       zdpmb = n_days_per_month
//       if (n_days_per_year == 365) then
//          jmonb = mod(kmonb+11,12) + 1 ! convert month (0-13) -> (1-12)
//          zdpmb = mondays(jmonb)
//          if (jmonb == 2) zdpmb = zdpmb + idatim(7) ! leap year
//       endif
//       zmedb = 0.5 * (zdpmb + 1.0) ! median day b
//
// !     compute weight
//
//       pweight = abs(zday - zmeda) / (zmeda + zmedb - 1.0)
//
//       return
//       end subroutine momint
//
