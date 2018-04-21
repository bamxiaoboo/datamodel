# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/biogeophys/SoilThermProp.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/biogeophys/SoilThermProp.F90"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/biogeophys/SoilThermProp.F90" 2

# 1 "./preproc.h" 1






 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/biogeophys/SoilThermProp.F90" 2

subroutine SoilThermProp (tk, cv, clm)

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
! Calculation of thermal conductivities and heat capacities of 
! snow/soil layers
!
! Method:
! (1) The volumetric heat capacity is calculated as a linear combination 
!     in terms of the volumetric fraction of the constituent phases. 
!
! (2) The thermal conductivity of soil is computed from the algorithm of
!     Johansen (as reported by Farouki 1981), and of snow is from the
!     formulation used in SNTHERM (Jordan 1991).
!
! The thermal conductivities at the interfaces between two neighboring 
! layers (j, j+1) are derived from an assumption that the flux across 
! the interface is equal to that from the node j to the interface and the 
! flux from the interface to the node j+1. 
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! April 2002: Vertenstein/Oleson/Levis; Final form
!
!-----------------------------------------------------------------------
! $Id: SoilThermProp.F90,v 1.2.10.3 2002/06/15 13:50:19 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  use clm_varcon, only : denh2o, denice, tfrz,   tkwat, tkice, tkair, &
                         cpice,  cpliq,  istice, istwet
  use clm_varpar, only : nlevsoi
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm             ! CLM 1-D Module

  real(r8), intent(out) :: cv(clm%snl+1:nlevsoi) ! heat capacity [J/(m2 K)]
  real(r8), intent(out) :: tk(clm%snl+1:nlevsoi) ! thermal conductivity [W/(m K)]

!----Local Variables----------------------------------------------------

  integer i                       ! do loop index
  real(r8) bw                     ! partial density of water (ice + liquid) [kg/m3]
  real(r8) dksat                  ! thermal conductivity for saturated soil [J/(K s m)]
  real(r8) dke                    ! kersten number [-]
  real(r8) fl                     ! fraction of liquid or unfrozen water to total water [-]
  real(r8) satw                   ! relative total water content of soil [-]
  real(r8) thk(clm%snl+1:nlevsoi) ! thermal conductivity of layer [J/(K s m)]

!----End Variable List--------------------------------------------------

!
! Thermal conductivity of soil from Farouki (1981)
!

  do i = 1, nlevsoi
     if (clm%itypwat/=istwet .AND. clm%itypwat/=istice) then
        satw = (clm%h2osoi_liq(i)/denh2o+  &
             clm%h2osoi_ice(i)/denice)/(clm%dz(i)*clm%watsat(i))
        satw = min(1._r8, satw)
        if (satw > .1e-6) then          
           fl = clm%h2osoi_liq(i)/(clm%h2osoi_ice(i)+clm%h2osoi_liq(i))
           if(clm%t_soisno(i) >= tfrz) then       ! Unfrozen soil
              dke = max(0._r8, log10(satw) + 1.0)
              dksat = clm%tksatu(i)
           else                                   ! Frozen soil
              dke = satw
              dksat = clm%tkmg(i)*0.249**(fl*clm%watsat(i))*2.29**clm%watsat(i)
           endif
           thk(i) = dke*dksat + (1.-dke)*clm%tkdry(i)
        else    
           thk(i) = clm%tkdry(i)
        endif
     else
        thk(i) = tkwat
        if (clm%t_soisno(i) < tfrz) thk(i) = tkice
     endif
  enddo

!
! Thermal conductivity of snow, from Jordan (1991) pp. 18
!

  if(clm%snl+1 < 1)then
     do i = clm%snl+1, 0
        bw = (clm%h2osoi_ice(i)+clm%h2osoi_liq(i))/clm%dz(i)
        thk(i) = tkair + (7.75e-5 *bw + 1.105e-6*bw*bw) &
             *(tkice-tkair)
     enddo
  endif

!
! Thermal conductivity at the layer interface
!

  do i = clm%snl+1, nlevsoi-1
     tk(i) = thk(i)*thk(i+1)*(clm%z(i+1)-clm%z(i)) &
          /(thk(i)*(clm%z(i+1)-clm%zi(i))+thk(i+1)*(clm%zi(i)-clm%z(i)))
  enddo
  tk(nlevsoi) = 0.

!
! Soil heat capacity, from de Vires (1963)
!

  do i = 1, nlevsoi
     if (clm%itypwat/=istwet .AND. clm%itypwat/=istice) then
        cv(i) = clm%csol(i)*(1-clm%watsat(i))*clm%dz(i) +   &
             (clm%h2osoi_ice(i)*cpice + clm%h2osoi_liq(i)*cpliq)
     else
        cv(i) = (clm%h2osoi_ice(i)*cpice + clm%h2osoi_liq(i)*cpliq)
     endif
  enddo
  if (clm%snl+1 == 1 .AND. clm%h2osno > 0.) cv(1) = cv(1) + cpice*clm%h2osno

!
! Snow heat capacity
!

  if (clm%snl+1 < 1)then
     do i = clm%snl+1, 0
        cv(i) = cpliq*clm%h2osoi_liq(i) + cpice*clm%h2osoi_ice(i)
     enddo
  endif

end subroutine SoilThermProp