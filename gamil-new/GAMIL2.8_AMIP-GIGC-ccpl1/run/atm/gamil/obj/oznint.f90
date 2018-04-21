# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/control/oznint.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/control/oznint.F90"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/control/oznint.F90" 2

# 1 "./params.h" 1
# 15 "./params.h"
 
 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/control/oznint.F90" 2

subroutine oznint
!----------------------------------------------------------------------- 
! 
! Purpose: Interpolate ozone mixing ratios to current time, reading in new monthly
!          data if necessary, and spatially interpolating it.
! 
! Method: Find next month of ozone data to interpolate.  Linearly interpolate 
!         vertically and horizontally
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use comozp
   use pspect
   use rgrid
   use commap
   use c_coupler_interface_mod

   use mpishorthand


   implicit none


# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/control/comctl.h" 1
!----------------------------------------------------------------------- 
! 
! Purpose: Model control variables
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
!! (wh 2003.04.30)
!! (wh 2003.12.27)

      common /comctl/ itsst   ,nsrest  ,iradsw  ,iradlw  ,iradae
      common /comctl/ anncyc  ,nlend   ,nlres   ,nlhst   ,lbrnch
      common /comctl/ ozncyc  ,sstcyc  ,icecyc
      common /comctl/ adiabatic
      common /comctl_l1/ flxave
      common /comctl/ trace_gas, trace_test1,trace_test2, trace_test3
!!    common /comctl/ readtrace,ideal_phys, nsplit, iord, jord, kord, use_eta, aqua_planet
      common /comctl/ readtrace,ideal_phys,                                    aqua_planet
      common /comctl/ doRamp_ghg, doRamp_so4, doRamp_scon, fullgrid, doIPCC_so4, &
                      doCmip5_scon,doCmip5_ghg  !!(wh)
      common /comctl/ print_step_cost
      common /comctl/ doabsems, dosw, dolw, indirect

!!    common /comctl_r8/ divdampn, precc_thresh, precl_thresh
      common /comctl_r8/           precc_thresh, precl_thresh

      integer itsst             ! Sea surf. temp. update freq. (iters)
      integer nsrest            ! Restart flag
      integer iradsw            ! Iteration freq. for shortwave radiation
      integer iradlw            ! Iteration freq. for longwave radiation
      integer iradae            ! Iteration freq. for absorptivity/emissivity

! f-v dynamics specific
! _ord = 1: first order upwind
! _ord = 2: 2nd order van Leer (Lin et al 1994)
! _ord = 3: standard PPM 
! _ord = 4: enhanced PPM (default)
!!      integer nsplit            ! Lagrangian time splits (Lin-Rood only)
!!      integer iord              ! scheme to be used in E-W direction
!!      integer jord              ! scheme to be used in N-S direction
!!      integer kord              ! scheme to be used for vertical mapping
!!      logical use_eta           ! Flag to use a's and b's set by dynamics/lr/set_eta.F90

      logical aqua_planet       ! Flag to run model in "aqua planet" mode

      logical anncyc            ! true => do annual cycle (otherwise perpetual)
      logical nlend             ! true => end of run
      logical nlres             ! true => continuation run
      logical nlhst             ! true => regeneration run
      logical lbrnch            ! true => branch run
      logical ozncyc            ! true => cycle ozone dataset
      logical sstcyc            ! true => cycle sst dataset
      logical icecyc            ! true => cycle ice fraction dataset
      logical adiabatic         ! true => no physics
      logical ideal_phys        ! true => run "idealized" model configuration
      logical(1) flxave            ! true => send to coupler only on radiation time steps

      logical trace_gas         ! true => turn on greenhouse gas code
      logical trace_test1       ! true => turn on test tracer code with 1 tracer
      logical trace_test2       ! true => turn on test tracer code with 2 tracers
      logical trace_test3       ! true => turn on test tracer code with 3 tracers
      logical readtrace         ! true => obtain initial tracer data from IC file

      logical doRamp_ghg        ! true => turn on ramping for ghg
      logical doRamp_so4        ! true => turn on ramping for so4
      logical doRamp_scon       ! true => turn on ramping for scon
      logical doIPCC_so4        ! true => turn on IPCC scenario for so4  !!(wh) 
      logical doCmip5_scon      !
      logical doCmip5_ghg       ! ljli2010-08-12
      logical fullgrid          ! true => no grid reduction towards poles

      logical print_step_cost   ! true => print per-timestep cost info

      logical doabsems          ! True => abs/emiss calculation this timestep
      logical dosw              ! True => shortwave calculation this timestep
      logical dolw              ! True => longwave calculation this timestep
      logical indirect          ! True => include indirect radiative effects of sulfate aerosols

!!    real(r8) divdampn         ! Number of days to invoke divergence damper
      real(r8) precc_thresh     ! Precipitation threshold for PRECCINT and PRECCFRQ
      real(r8) precl_thresh     ! Precipitation threshold for PRECLINT and PRECLFRQ
# 30 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/control/oznint.F90" 2

# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/control/comlun.h" 1
!----------------------------------------------------------------------- 
! 
! Purpose: Logical unit numbers and related variables
!
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------

      common /comlun/ nsds    ,nrg     ,nrg2
      common /comlun/ ncid_ini,ncid_oz ,ncid_sst, ncid_trc
      common /comlun/ luhrest

      integer nsds       ! restart dataset unit
      integer nrg        ! master regeneration dataset unit
      integer nrg2       ! abs/ems regeneration dataset units
      integer ncid_ini   ! initial dataset
      integer ncid_oz    ! ozone dataset
      integer ncid_trc   ! greenhouse gas tracer dataset
      integer ncid_sst   ! sst dataset
      integer luhrest    ! history restart unit
# 31 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/control/oznint.F90" 2
!
! Local workspace
!
   integer cnt4(4)                ! array of counts for each dimension
   integer strt4(4)               ! array of starting indices
   integer i, k, lat              ! longitude, level, latitude indices
   integer ntmp                   ! temporary
   integer :: yr, mon, day        ! components of a date
   integer :: ncdate              ! current date in integer format [yyyymmdd]
   integer :: ncsec               ! current time of day [seconds]

   real(r8) fact1, fact2          ! time interpolation factors
   real(r8) :: calday             ! current calendar day
   real(r8) caldayloc             ! calendar day (includes yr if no cycling)
   real(r8) tmpozmix(levsiz,plat) ! temporary ozmix array
   real(r8) deltat                ! time (days) between interpolating ozone data

   real(r8), allocatable :: oznbdyp(:,:,:) ! ozone data from next time sample
!
! Use year information only if a multiyear dataset
!
   call c_coupler_get_current_calendar_time(calday)
   call c_coupler_get_current_time(yr, mon, day, ncsec)
   ncdate = yr*10000 + mon*100 + day
   if (ozncyc) then
      caldayloc = calday
   else
      caldayloc = calday + yr*365.
   end if

   strt4(1) = 1
   strt4(2) = 1
   strt4(3) = 1
   cnt4(1)  = lonsiz
   cnt4(2)  = levsiz
   cnt4(3)  = latsiz
   cnt4(4)  = 1
!
! If model time is past current forward ozone timeslice, read in the next
! timeslice for time interpolation.  Messy logic is for ozncyc = .true. 
! interpolation between December and January (np1 == 1).  Note that 
! np1 is never 1 when ozncyc is .false.
!
   if (caldayloc > cdayozp .and. .not. (np1 == 1 .and. caldayloc > cdayozm)) then
      if (ozncyc) then
         np1 = mod(np1,12) + 1
      else
         np1 = np1 + 1
      end if
      if (np1 > timesiz) then
         write(6,*)'OZNINT: Attempt to read past end of O3 dataset'
         call endrun
      end if
      cdayozm = cdayozp
      call bnddyi(date_oz(np1), sec_oz(np1), cdayozp)
      if (.not.ozncyc) then
         yr = date_oz(np1)/10000
         cdayozp = cdayozp + yr*365.
      end if
      if (np1 == 1 .or. caldayloc <= cdayozp) then
         ntmp = nm
         nm = np
         np = ntmp
         strt4(4) = np1
         if (masterproc) then
            allocate (oznbdyp(lonsiz,levsiz,latsiz))
            call wrap_get_vara_realx (ncid_oz,oznid,strt4,cnt4,oznbdyp)
            write(6,*)'OZNINT: Read ozone for date (yyyymmdd) ', date_oz(np1),' sec ',sec_oz(np1)
!
! Spatial interpolation.  If ozone dataset is only 2-d (i.e. lonsiz = 1) and 
! thus only latitude interpolation is necessary, expand to 3-d after 
! interpolation.
!
            if (lonsiz == 1) then
               call lininterp (oznbdyp, ozlat, levsiz, latsiz, tmpozmix, &
                               latdeg, plat)
               do lat=1,plat
                  do k=1,levsiz
                     do i=1,nlon(lat)
                        ozmixm(i,k,lat,np) = tmpozmix(k,lat)
                     end do
                  end do
               end do
            else
            call bilin (oznbdyp ,ozlon   ,ozlat   ,lonsiz  ,lonsiz  , &
                        levsiz  ,levsiz  ,latsiz  ,ozmixm(1,1,1,np) , londeg  , &
                        latdeg  ,plond   ,nlon    ,levsiz  ,plat    )
            end if
            deallocate (oznbdyp)
         end if
      else
         write(6,*)'OZNINT: Input ozone for date',date_oz(np1),' sec ',sec_oz(np1), &
              'does not exceed model date',ncdate,' sec ',ncsec,' Stopping.'
         call endrun
      end if

      call mpibcast(ozmixm(1,1,1,np), plond*levsiz*plat, mpir8, 0, mpicom)
   end if
!
! Determine time interpolation factor.  Account for December-January 
! interpolation if cycling ozone dataset.  Again note that np1 is never 1 
! when ozncyc is false
!
   if (np1 == 1) then                    ! Dec-Jan interpolation
      deltat = cdayozp + 365. - cdayozm
      if (caldayloc > cdayozp) then      ! We're in December
         fact1 = (cdayozp + 365. - caldayloc)/deltat
         fact2 = (caldayloc - cdayozm)/deltat
      else                                ! We're in January
         fact1 = (cdayozp - caldayloc)/deltat
         fact2 = (caldayloc + 365. - cdayozm)/deltat
      end if
   else
      deltat = cdayozp - cdayozm
      fact1 = (cdayozp - caldayloc)/deltat
      fact2 = (caldayloc - cdayozm)/deltat
   end if
!
! Check sanity of time interpolation calculation to within 32-bit roundoff
!
   if (abs(fact1+fact2-1.) > 1.e-6 .or. fact1 > 1.000001 .or. fact1 < -1.e-6 .or. &
       fact2 > 1.000001 .or. fact2 < -1.e-6) then
      write(6,*)'OZNINT: Bad fact1 and/or fact2=',fact1,fact2
      call endrun
   end if
!
! Time interpolation.
!
!$OMP PARALLEL DO PRIVATE (lat,k,i)
   do lat=1,plat
      do k=1,levsiz
         do i=1,nlon(lat)
            ozmix(i,k,lat) = ozmixm(i,k,lat,nm)*fact1 + ozmixm(i,k,lat,np)*fact2
         end do
      end do
   end do

   return
end subroutine oznint

