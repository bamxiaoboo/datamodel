# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/biogeophys/SnowCompaction.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/biogeophys/SnowCompaction.F90"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/biogeophys/SnowCompaction.F90" 2

# 1 "./preproc.h" 1






 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/biogeophys/SnowCompaction.F90" 2

subroutine SnowCompaction(clm)

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
! Determine the change in snow layer thickness due to compaction and 
! settling. 
!
! Method:
! Three metamorphisms of changing snow characteristics are implemented, 
! i.e., destructive, overburden, and melt. The treatments of the former 
! two are from SNTHERM.89 and SNTHERM.99 (1991, 1999). The contribution 
! due to melt metamorphism is simply taken as a ratio of snow ice 
! fraction after the melting versus before the melting. 
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! April 2002: Vertenstein/Oleson/Levis; Final form
!
!-----------------------------------------------------------------------
! $Id: SnowCompaction.F90,v 1.2.10.3 2002/06/15 13:50:18 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  use clm_varcon, only : denice, denh2o, tfrz
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm	 !CLM 1-D Module

!----Local Variables----------------------------------------------------

  integer i       ! do loop index
  real(r8) c2     ! = 23e-3 [m3/kg]
  real(r8) c3     ! = 2.777e-6 [1/s]
  real(r8) c4     ! = 0.04 [1/K]
  real(r8) c5     ! = 2.0 [-]
  real(r8) dm     ! Upper Limit on Destructive Metamorphism Compaction [kg/m3]
  real(r8) eta0   ! The Viscosity Coefficient Eta0 [kg-s/m2]
  real(r8) burden ! pressure of overlying snow [kg/m2]
  real(r8) ddz1   ! Rate of settling of snowpack due to destructive metamorphism [1/s]
  real(r8) ddz2   ! Rate of compaction of snowpack due to overburden [1/s]
  real(r8) ddz3   ! Rate of compaction of snowpack due to melt [1/s]
  real(r8) dexpf  ! expf=exp(-c4*(tfrz[K] -t_soisno)).
  real(r8) fi     ! Fraction of ice relative to the total water content at current time step
  real(r8) td     ! t_soisno - tfrz [K] 
  real(r8) pdzdtc ! Nodal rate of change in fractional-thickness due to compaction [fraction/s]
  real(r8) void   ! void (1 - vol_ice - vol_liq) [-]
  real(r8) wx     ! water mass (ice+liquid) [kg/m2]
  real(r8) bi     ! partial density of ice [kg/m3]

  data c2,c3,c4,c5/23.e-3, 2.777e-6, 0.04, 2.0/
  data dm/100./         
  data eta0/9.e5/      

!----End Variable List--------------------------------------------------

  burden = 0.0

  do i = clm%snl+1, 0

     wx = clm%h2osoi_ice(i) + clm%h2osoi_liq(i)
     void=1.-(clm%h2osoi_ice(i)/denice+clm%h2osoi_liq(i)/denh2o)/clm%dz(i)

!
! Disallow compaction for water saturated node and lower ice lens node.
!

     if(void <= 0.001 .or. clm%h2osoi_ice(i) <= .1)then
        burden = burden+wx
        cycle
     endif

     bi = clm%h2osoi_ice(i) / clm%dz(i)
     fi = clm%h2osoi_ice(i) / wx
     td = tfrz-clm%t_soisno(i)
     dexpf = exp(-c4*td)

!
! Settling as a result of destructive metamorphism
!

     ddz1 = -c3*dexpf
     if(bi > dm) ddz1 = ddz1*exp(-46.0e-3*(bi-dm))

!
! Liquid water term
! 

     if(clm%h2osoi_liq(i) > 0.01*clm%dz(i)) ddz1=ddz1*c5

!
! Compaction due to overburden
!

     ddz2 = -burden*exp(-0.08*td - c2*bi)/eta0

!
! Compaction occurring during melt
!

     if(clm%imelt(i) == 1)then
        ddz3 = - 1./clm%dtime * max(0._r8,(clm%frac_iceold(i) - fi)/clm%frac_iceold(i))
     else
        ddz3 = 0.
     endif

!
! Time rate of fractional change in dz (units of s-1)
!

     pdzdtc = ddz1+ddz2+ddz3

!
! The change in dz due to compaction
!

     clm%dz(i) = clm%dz(i)*(1.+pdzdtc*clm%dtime)

!
! Pressure of overlying snow
!

     burden = burden+wx

  enddo

end subroutine SnowCompaction
