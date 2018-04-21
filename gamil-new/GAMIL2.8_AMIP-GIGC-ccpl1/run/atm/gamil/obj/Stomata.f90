# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/biogeophys/Stomata.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/biogeophys/Stomata.F90"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/biogeophys/Stomata.F90" 2

# 1 "./preproc.h" 1






 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/biogeophys/Stomata.F90" 2

subroutine Stomata(mpe,  apar,   ei,    ea,   tgcm,   &
                   o2,   co2,    btran, rb,   rs,     &
                   psn,  qe25,   vcmx25,mp,   c3psn,  &
                   clm   ) 

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
! Leaf stomatal resistance and leaf photosynthesis.
!
! Method:
!
! Author:
! author:            Gordon Bonan
! standardized:      J. Truesdale, Feb. 1996
! reviewed:          G. Bonan, Feb. 1996
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! April 2002: Vertenstein/Oleson/Levis; Final form
!
!-----------------------------------------------------------------------
! $Id: Stomata.F90,v 1.5.6.5.6.1 2002/10/03 20:07:28 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  use clm_varcon   , only : tfrz
  use shr_const_mod, only : SHR_CONST_TKFRZ,SHR_CONST_RGAS
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm	 !CLM 1-D Module

  real(r8), intent(in) :: mpe    ! prevents division by zero errors
  real(r8), intent(in) :: ei     ! vapor pressure inside leaf (sat vapor press at t_veg) [pa]
  real(r8), intent(in) :: ea     ! vapor pressure of canopy air [pa]
  real(r8), intent(in) :: apar   ! par absorbed per unit lai [W/m2]
  real(r8), intent(in) :: o2     ! atmospheric o2 concentration [pa]
  real(r8), intent(in) :: co2    ! atmospheric co2 concentration [pa]
  real(r8), intent(in) :: tgcm   ! air temperature at agcm reference height [K]
  real(r8), intent(in) :: btran  ! soil water transpiration factor (0 to 1)
  real(r8), intent(in) :: qe25   ! quantum efficiency at 25c [umol co2 / umol photon]
  real(r8), intent(in) :: vcmx25 ! maximum rate of carboxylation at 25c [umol co2/m2/s]
  real(r8), intent(in) :: mp     ! slope for conductance-to-photosynthesis relationship 
  real(r8), intent(in) :: c3psn  ! photosynthetic pathway: 0. = c4, 1. = c3

  real(r8), intent(inout) :: rb  ! boundary layer resistance [s/m]

  real(r8), intent(out)   :: rs  ! leaf stomatal resistance [s/m]
  real(r8), intent(out)   :: psn ! foliage photosynthesis [umol co2/m2/s] [always +]

!----Local Variables----------------------------------------------------

  integer, parameter :: niter = 3  ! number of iterations
  integer  iter                    ! iteration index

  real(r8) ab      ! used in statement functions
  real(r8) bc      ! used in statement functions
  real(r8) f1      ! generic temperature response (statement function)
  real(r8) f2      ! generic temperature inhibition (statement function)
  real(r8) tc      ! foliage temperature [C]
  real(r8) cs      ! co2 concentration at leaf surface [pa]
  real(r8) kc      ! co2 michaelis-menten constant [pa]
  real(r8) ko      ! o2 michaelis-menten constant [pa]
  real(r8) a,b,c,q ! intermediate calculations for rs
  real(r8) r1,r2   ! roots for rs
  real(r8) ppf     ! absorbed photosynthetic photon flux [umol photons/m2/s]
  real(r8) wc      ! rubisco limited photosynthesis [umol co2/m2/s]
  real(r8) wj      ! light limited photosynthesis [umol co2/m2/s]
  real(r8) we      ! export limited photosynthesis [umol co2/m2/s]
  real(r8) cp      ! co2 compensation point [pa]
  real(r8) ci      ! internal co2 [pa]
  real(r8) awc     ! intermediate calculation for wc
  real(r8) vcmx    ! maximum rate of carboxylation [umol co2/m2/s]
  real(r8) j       ! electron transport [umol co2/m2/s]
  real(r8) cea     ! constrain ea or else model blows up
  real(r8) cf      ! s m**2/umol -> s/m
  real(r8) rsmax0  ! maximum stomatal resistance [s/m]

  real(r8) kc25    ! co2 michaelis-menten constant at 25c [pa]
  real(r8) akc     ! q10 for kc25
  real(r8) ko25    ! o2 michaelis-menten constant at 25c [pa]
  real(r8) ako     ! q10 for ko25
  real(r8) avcmx   ! q10 for vcmx25
  real(r8) bp      ! minimum leaf conductance [umol/m2/s]

!----End Variable List--------------------------------------------------

  f1(ab,bc) = ab**((bc-25.)/10.)
  f2(ab) = 1. + exp((-2.2e05+710.*(ab+SHR_CONST_TKFRZ))/(SHR_CONST_RGAS*0.001*(ab+SHR_CONST_TKFRZ)))

  kc25 = 30.
  akc = 2.1
  ko25 = 30000.
  ako = 1.2
  avcmx = 2.4
  bp = 2000.

!
! Initialize rs=rsmax and psn=0 because calculations are performed only
! when apar > 0, in which case rs <= rsmax and psn >= 0
! Set constants
!

  rsmax0 = 2.e4
  cf = clm%forc_pbot/(SHR_CONST_RGAS*0.001*tgcm)*1.e06 

  if (apar <= 0.) then          ! night time
     rs = min(rsmax0, 1./bp * cf)
     psn = 0.
     return
  else                          ! day time
     tc = clm%t_veg-tfrz                            
     ppf = 4.6*apar                  
     j = ppf*qe25
     kc = kc25 * f1(akc,tc)       
     ko = ko25 * f1(ako,tc)
     awc = kc * (1.+o2/ko)
     cp = 0.5*kc/ko*o2*0.21
     vcmx = vcmx25 * f1(avcmx,tc) / f2(tc) * btran
!
! First guess ci
!

     ci = 0.7*co2*c3psn + 0.4*co2*(1.-c3psn)  

!
! rb: s/m -> s m2 / umol
!

     rb = rb/cf 

!
! Constrain ea
!

     cea = max(0.25*ei*c3psn+0.40*ei*(1.-c3psn), min(ea,ei) ) 

!
! ci iteration for 'actual' photosynthesis
!

     do iter = 1, niter
        wj = max(ci-cp,0._r8)*j/(ci+2.*cp)*c3psn + j*(1.-c3psn)
        wc = max(ci-cp,0._r8)*vcmx/(ci+awc)*c3psn + vcmx*(1.-c3psn)
        we = 0.5*vcmx*c3psn + 4000.*vcmx*ci/clm%forc_pbot*(1.-c3psn) 
        psn = min(wj,wc,we) 
        cs = max( co2-1.37*rb*clm%forc_pbot*psn, mpe )
        a = mp*psn*clm%forc_pbot*cea / (cs*ei) + bp
        b = ( mp*psn*clm%forc_pbot/cs + bp ) * rb - 1.
        c = -rb
        if (b >= 0.) then
           q = -0.5*( b + sqrt(b*b-4.*a*c) )
        else
           q = -0.5*( b - sqrt(b*b-4.*a*c) )
        endif
        r1 = q/a
        r2 = c/q
        rs = max(r1,r2)
        ci = max( cs-psn*clm%forc_pbot*1.65*rs, 0._r8 )
     enddo

     ! rs, rb:  s m2 / umol -> s/m 
     rs = min(rsmax0, rs*cf)
     rb = rb*cf 

  endif

end subroutine Stomata
