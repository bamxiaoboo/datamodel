# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/control/infnan.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/control/infnan.F90"
module infnan
!-------------------------------------------------------------------------
!
! Purpose:
!
!	Set parameters for the floating point flags "inf" Infinity
!	and "nan" not-a-number. As well as "bigint" the point
!	at which integers start to overflow. These values are used
!	to initialize arrays with as a way to detect if arrays
!	are being used before being set.
!
! Author: CCM Core group
!
! $Id: infnan.F90,v 1.2.8.1 2002/06/15 13:50:08 erik Exp $
!
!-------------------------------------------------------------------------
    use shr_kind_mod, only: r8 => shr_kind_r8





    real(r8), parameter :: inf = O'0777600000000000000000'
    real(r8), parameter :: nan = O'0777700000000000000000'
    integer,  parameter :: bigint = O'17777777777' ! largest possible 32-bit integer

    real(r8), parameter :: uninit_r8 = inf
end module infnan

