# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/biogeophys/SoilTemperature.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/biogeophys/SoilTemperature.F90"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/biogeophys/SoilTemperature.F90" 2

# 1 "./preproc.h" 1






 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/biogeophys/SoilTemperature.F90" 2

subroutine SoilTemperature (clm,   tssbef, htvp, emg, cgrnd, &
                            dlrad, tg,     xmf,  fact )

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
! Snow and soil temperatures including phase change
!
! Method:
! o The volumetric heat capacity is calculated as a linear combination 
!   in terms of the volumetric fraction of the constituent phases. 
! o The thermal conductivity of soil is computed from 
!   the algorithm of Johansen (as reported by Farouki 1981), and the 
!   conductivity of snow is from the formulation used in
!   SNTHERM (Jordan 1991).
! o Boundary conditions:  
!   F = Rnet - Hg - LEg (top),  F = 0 (base of the soil column).
! o Soil / snow temperature is predicted from heat conduction 
!   in 10 soil layers and up to 5 snow layers. 
!   The thermal conductivities at the interfaces between two 
!   neighboring layers (j, j+1) are derived from an assumption that 
!   the flux across the interface is equal to that from the node j 
!   to the interface and the flux from the interface to the node j+1. 
!   The equation is solved using the Crank-Nicholson method and 
!   results in a tridiagonal system of equations.
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! April 2002: Vertenstein/Oleson/Levis; Final form
!
!-----------------------------------------------------------------------
! $Id: SoilTemperature.F90,v 1.1.10.3 2002/06/15 13:50:18 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  use clm_varcon, only : sb
  use clm_varpar, only : nlevsoi
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout)  :: clm	 !CLM 1-D Module

  real(r8), intent(in) :: tssbef(-nlevsno:nlevsoi) ! soil/snow temperature before update [K]
  real(r8), intent(in) :: htvp  ! latent heat of vapor of water (or sublimation) [J/kg]
  real(r8), intent(in) :: emg   ! ground emissivity [-]
  real(r8), intent(in) :: cgrnd ! deriv. of soil energy flux wrt to soil temp [W/m2/K]
  real(r8), intent(in) :: dlrad ! downward longwave radiation below the canopy [W/m2]

  real(r8), intent(out) :: xmf  ! latent heat of phase change of ground water [W/m2]
  real(r8), intent(out) :: fact(clm%snl+1 : nlevsoi)  ! used in computing tridiagonal matrix

  real(r8), intent(inout) :: tg ! ground surface temperature [K]

!----Local Variables----------------------------------------------------

  integer i,j                        ! do loop indices
  real(r8) at(clm%snl+1 : nlevsoi)   ! "a" vector for tridiagonal matrix
  real(r8) bt(clm%snl+1 : nlevsoi)   ! "b" vector for tridiagonal matrix
  real(r8) ct(clm%snl+1 : nlevsoi)   ! "c" vector for tridiagonal matrix
  real(r8) rt(clm%snl+1 : nlevsoi)   ! "r" vector for tridiagonal solution
  real(r8) cv(clm%snl+1 : nlevsoi)   ! heat capacity [J/(m2 K)]
  real(r8) tk(clm%snl+1 : nlevsoi)   ! thermal conductivity [W/(m K)]
  real(r8) fn  (clm%snl+1 : nlevsoi) ! heat diffusion through the layer interface [W/m2]
  real(r8) fn1 (clm%snl+1 : nlevsoi) ! heat diffusion through the layer interface [W/m2]
  real(r8) dzm                       ! used in computing tridiagonal matrix
  real(r8) dzp                       ! used in computing tridiagonal matrix
  real(r8) hs                        ! net energy flux into the surface [W/m2]
  real(r8) brr(clm%snl+1 : nlevsoi)  ! temporary variable
  real(r8) dhsdT                     ! d(hs)/dT

!----End Variable List--------------------------------------------------

!
! [1] Ground surface and soil temperatures 
!

!
! [1.1] Thermal conductivity and Heat capacity
!

  call SoilThermProp (tk, cv, clm)

!
! [1.2] Net ground heat flux into the surface and its temperature derivative
!
     
  hs    = clm%sabg + dlrad &
       + (1-clm%frac_veg_nosno)*emg*clm%forc_lwrad - emg*sb*tg**4 &
       - (clm%eflx_sh_grnd+clm%qflx_evap_soi*htvp) 

  dhsdT = - cgrnd - 4.*emg * sb * tg**3

  j       = clm%snl+1
  fact(j) = clm%dtime / cv(j) &
       * clm%dz(j) / (0.5*(clm%z(j)-clm%zi(j-1)+clm%capr*(clm%z(j+1)-clm%zi(j-1))))

  do j = clm%snl+1 + 1, nlevsoi
     fact(j) = clm%dtime/cv(j)
  enddo

  do j = clm%snl+1, nlevsoi - 1
     fn(j) = tk(j)*(clm%t_soisno(j+1)-clm%t_soisno(j))/(clm%z(j+1)-clm%z(j))
  enddo
  fn(nlevsoi) = 0.

!
! [1.3] Set up vector r and vectors a, b, c that define tridiagonal matrix
!       and solve system
!

  j     = clm%snl+1
  dzp   = clm%z(j+1)-clm%z(j)
  at(j) = 0.
  bt(j) = 1+(1.-clm%cnfac)*fact(j)*tk(j)/dzp-fact(j)*dhsdT
  ct(j) =  -(1.-clm%cnfac)*fact(j)*tk(j)/dzp
  rt(j) = clm%t_soisno(j) +  fact(j)*( hs - dhsdT*clm%t_soisno(j) + clm%cnfac*fn(j) )

  do j    = clm%snl+1 + 1, nlevsoi - 1
     dzm   = (clm%z(j)-clm%z(j-1))
     dzp   = (clm%z(j+1)-clm%z(j))

     at(j) =   - (1.-clm%cnfac)*fact(j)* tk(j-1)/dzm
     bt(j) = 1.+ (1.-clm%cnfac)*fact(j)*(tk(j)/dzp + tk(j-1)/dzm)
     ct(j) =   - (1.-clm%cnfac)*fact(j)* tk(j)/dzp

     rt(j) = clm%t_soisno(j) + clm%cnfac*fact(j)*( fn(j) - fn(j-1) )
  enddo

  j     =  nlevsoi
  dzm   = (clm%z(j)-clm%z(j-1))
  at(j) =   - (1.-clm%cnfac)*fact(j)*tk(j-1)/dzm
  bt(j) = 1.+ (1.-clm%cnfac)*fact(j)*tk(j-1)/dzm
  ct(j) = 0.
  rt(j) = clm%t_soisno(j) - clm%cnfac*fact(j)*fn(j-1)

  i = size(at)
  call Tridiagonal (i, at, bt, ct, rt, &
                    clm%t_soisno(clm%snl+1:nlevsoi))

!
! [2] Melting or Freezing 
!

  do j = clm%snl+1, nlevsoi - 1
     fn1(j) = tk(j)*(clm%t_soisno(j+1)-clm%t_soisno(j))/(clm%z(j+1)-clm%z(j))
  enddo
  fn1(nlevsoi) = 0.

  j = clm%snl+1
  brr(j) = clm%cnfac*fn(j) + (1.-clm%cnfac)*fn1(j)

  do j = clm%snl+1 + 1, nlevsoi
     brr(j) = clm%cnfac*(fn(j)-fn(j-1)) + (1.-clm%cnfac)*(fn1(j)-fn1(j-1))
  enddo

  call PhaseChange (fact(clm%snl+1),   brr(clm%snl+1), hs, dhsdT, &
                    tssbef(clm%snl+1), xmf,            clm )

  tg = clm%t_soisno(clm%snl+1)

end subroutine SoilTemperature
