# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/biogeophys/QSat.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/biogeophys/QSat.F90"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/biogeophys/QSat.F90" 2

# 1 "./preproc.h" 1






 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/biogeophys/QSat.F90" 2

subroutine QSat (T, p, es, esdT, qs, &
                 qsdT)

!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely
!  L                           M  available land surface process model.
!  M --COMMUNITY LAND MODEL--  C
!  C                           L
!  LMCLMCLMCLMCLMCLMCLMCLMCLMCLM
!
!-----------------------------------------------------------------------
! Purpose:
! Computes saturation mixing ratio and the change in saturation
! mixing ratio with respect to temperature.

! Method:
! Reference:  Polynomial approximations from:
!             Piotr J. Flatau, et al.,1992:  Polynomial fits to saturation 
!             vapor pressure.  Journal of Applied Meteorology, 31, 1507-1513.
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! April 2002: Vertenstein/Oleson/Levis; Final form
!
!-----------------------------------------------------------------------
! $Id: QSat.F90,v 1.1.10.5 2002/06/15 13:50:17 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_const_mod, ONLY: SHR_CONST_TKFRZ
  implicit none

!----Arguments----------------------------------------------------------

  real(r8), intent(in)  :: T        ! temperature (K)
  real(r8), intent(in)  :: p        ! surface atmospheric pressure (pa)

  real(r8), intent(out) :: es       ! vapor pressure (pa)
  real(r8), intent(out) :: esdT     ! d(es)/d(T)
  real(r8), intent(out) :: qs       ! humidity (kg/kg)
  real(r8), intent(out) :: qsdT     ! d(qs)/d(T)

!----Local Variables----------------------------------------------------

  real(r8) T_limit                  ! limitation on valid temperatures [C]

  real(r8) td,vp,vp1,vp2
  real(r8) a0,a1,a2,a3,a4,a5,a6,a7,a8
  real(r8) b0,b1,b2,b3,b4,b5,b6,b7,b8

  real(r8) c0,c1,c2,c3,c4,c5,c6,c7,c8
  real(r8) d0,d1,d2,d3,d4,d5,d6,d7,d8

!
! For water vapor (temperature range 0C-100C)
!

  data a0/6.11213476   /,a1/ 0.444007856    /,a2/0.143064234e-01/  &
      ,a3/0.264461437e-03/,a4/ 0.305903558e-05/,a5/0.196237241e-07/  &
      ,a6/0.892344772e-10/,a7/-0.373208410e-12/,a8/0.209339997e-15/

!
! For derivative:water vapor
!

  data b0/0.444017302  /,b1/ 0.286064092e-01/,b2/ 0.794683137e-03/ &
      ,b3/ 0.121211669e-04/,b4/ 0.103354611e-06/,b5/ 0.404125005e-09/ &
      ,b6/-0.788037859e-12/,b7/-0.114596802e-13/,b8/ 0.381294516e-16/

!
! For ice (temperature range -75C-0C) 
!

  data c0/6.11123516     /,c1/0.503109514    /,c2/0.188369801e-01/ &
      ,c3/0.420547422e-03/,c4/0.614396778e-05/,c5/0.602780717e-07/ &
      ,c6/0.387940929e-09/,c7/0.149436277e-11/,c8/0.262655803e-14/

!
! For derivative:ice  
!

  data d0/0.503277922    /,d1/0.377289173e-01/,d2/0.126801703e-02/ &
      ,d3/0.249468427e-04/,d4/0.313703411e-06/,d5/0.257180651e-08/ &
      ,d6/0.133268878e-10/,d7/0.394116744e-13/,d8/0.498070196e-16/

!----End Variable List--------------------------------------------------

  T_limit = T - SHR_CONST_TKFRZ
  if (T_limit > 100.0) T_limit=100.0
  if (T_limit < -75.0) T_limit=-75.0

  td       = T_limit
  if (td >= 0.0) then
     es   = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
          + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
     esdT = b0 + td*(b1 + td*(b2 + td*(b3 + td*(b4 &
          + td*(b5 + td*(b6 + td*(b7 + td*b8)))))))
  else
     es   = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
          + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
     esdT = d0 + td*(d1 + td*(d2 + td*(d3 + td*(d4 &
          + td*(d5 + td*(d6 + td*(d7 + td*d8)))))))
  endif

  es    = es    * 100.            ! pa
  esdT  = esdT  * 100.            ! pa/K

  vp    = 1.0   / (p - 0.378*es)
  vp1   = 0.622 * vp
  vp2   = vp1   * vp

  qs    = es    * vp1             ! kg/kg
  qsdT  = esdT  * vp2 * p         ! 1 / K

end subroutine QSat
