# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/surface_subroutines_interface.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/surface_subroutines_interface.F90"
!***************************************************************
!  This is a source file of GAMIL, containing all model 
!  subroutine interfaces which will be registerred into 
!  C-Coupler library. This file was initially finished by
!  Dr. Li Liu. If you have any problem, please contact Dr. Li 
!  Liu via liuli-cess@tsinghua.edu.cn
!***************************************************************



# 1 "./misc.h" 1
# 11 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/surface_subroutines_interface.F90" 2


subroutine initialize_online_lsm_ccpl_interface

   use surface_subroutines_from_datafile_mod
   call initialize_online_lsm

end subroutine initialize_online_lsm_ccpl_interface



subroutine calculate_ocn_albedo_ccpl_interface

   use surface_subroutines_from_datafile_mod
   call calculate_ocn_albedo
   
end subroutine calculate_ocn_albedo_ccpl_interface



subroutine calculate_sice_albedo_ccpl_interface

   use surface_subroutines_from_datafile_mod
   call calculate_sice_albedo
   
end subroutine calculate_sice_albedo_ccpl_interface



subroutine initialize_data_ocn_ccpl_interface

   use surface_subroutines_from_datafile_mod
   call initialize_data_ocn

end subroutine initialize_data_ocn_ccpl_interface



subroutine initialize_data_sice_ccpl_interface

   use surface_subroutines_from_datafile_mod
   call initialize_data_sice

end subroutine initialize_data_sice_ccpl_interface



subroutine merge_ts_ccpl_interface 

   use surface_subroutines_from_datafile_mod
   call merge_ts

end subroutine merge_ts_ccpl_interface 



subroutine get_initial_surface_data_from_coupler_ccpl_interface
   use surface_subroutines_with_coupler_mod

   call initialize_with_coupler

end subroutine get_initial_surface_data_from_coupler_ccpl_interface


subroutine get_forcing_data_sst_ccpl_interface


    use sst_data,            only: sstint
    !
    ! Time interpolate sst data
    !
    call sstint ()


end subroutine get_forcing_data_sst_ccpl_interface



subroutine get_forcing_data_sice_ccpl_interface


    use ice_data,            only: iceint

    call iceint ()


end subroutine get_forcing_data_sice_ccpl_interface



subroutine srfflx_state_reset_ccpl_interface

    use comsrf

    call srfflx_state_reset(srfflx_state2d)

end subroutine srfflx_state_reset_ccpl_interface




subroutine update_srf_fractions_ccpl_interface

    use comsrf

    call update_srf_fractions

end subroutine update_srf_fractions_ccpl_interface



subroutine process_surface_data_with_online_lsm_ccpl_interface


    use comsrf
    use atm_lndMod,     only: atmlnd_drv
    use c_coupler_interface_mod

    implicit none


# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/physics/cam1/comsol.h" 1
!
!	Common's to do with solar radiation
!
!	$Id: comsol.h,v 1.3 2000/06/02 16:20:40 jet Exp $
!
! Visible optical depth
!
      real(r8) tauvis     ! Visible optical depth

      common /comvis/ tauvis
!
! Solar constant
!
      real(r8) scon       ! Solar constant

      common /comsol/ scon
!
! Earth's orbital characteristics
!	
      real(r8) eccen       ! Earth's eccentricity factor (unitless) (typically 0 to 0.1)
      real(r8) obliq       ! Earth's obliquity angle (degree's) (-90 to +90) (typically 22-26)
      real(r8) mvelp       ! Earth's moving vernal equinox at perhelion (degree's) (0 to 360.0)
      integer iyear_AD ! Year (AD) to simulate above earth's orbital parameters for
!
! Orbital information after processed by orbit_params
!
      real(r8) obliqr      ! Earth's obliquity in radians
      real(r8) lambm0      ! Mean longitude of perihelion at the 
!                          ! vernal equinox (radians)
      real(r8) mvelpp      ! Earth's moving vernal equinox longitude
!                          ! of perihelion plus pi (radians)
!
      common /comorb/ eccen   , obliq   , mvelp   , obliqr  
      common /comorb/ lambm0  , mvelpp  , iyear_AD

# 132 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/surface_subroutines_interface.F90" 2

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
# 133 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/surface_subroutines_interface.F90" 2

    integer nstep                              ! current timestep number

    nstep = c_coupler_get_nstep()

    if (.not. aqua_planet) then
        !
        ! 3.1 Call land model driving routine
        !
        call t_startf('atmlnd_drv')
        call atmlnd_drv(nstep, iradsw, eccen, obliqr, lambm0, mvelpp, surface_state2d, srfflx_parm2d)
        call t_stopf ('atmlnd_drv')

        call update_srf_fluxes(srfflx_state2d, srfflx_parm2d, landfrac)
    end if


end subroutine process_surface_data_with_online_lsm_ccpl_interface



subroutine process_surface_data_with_forcing_data_sst_ccpl_interface


    use comsrf

    implicit none

    call t_startf('camoce')
    call camoce(surface_state2d, srfflx_parm2d)
    call t_stopf('camoce')
    call update_srf_fluxes(srfflx_state2d, srfflx_parm2d, ocnfrac)


end subroutine process_surface_data_with_forcing_data_sst_ccpl_interface



subroutine process_surface_data_with_forcing_data_sice_ccpl_interface


    use comsrf

    implicit none

    call t_startf('camice')
    call camice(surface_state2d, srfflx_parm2d)
    call t_stopf('camice')
    call update_srf_fluxes(srfflx_state2d, srfflx_parm2d, icefrac)


end subroutine process_surface_data_with_forcing_data_sice_ccpl_interface



subroutine process_surface_data_with_coupler_ccpl_interface

    use surface_subroutines_with_coupler_mod
    use c_coupler_interface_mod

    implicit none
    
    integer,parameter :: R8  = selected_real_kind(12) ! 8 byte real

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
# 197 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/surface_subroutines_interface.F90" 2
    integer nstep                              ! current timestep number

    nstep = c_coupler_get_nstep()

    if (flxave) then
        !
        ! Average the precipitation input to lsm between radiation calls.
        !
       call average_atm_flux_variables(iradsw, nstep, dosw)
    end if
    !
    ! Send/recv data to/from the csm flux coupler.
    !
    call send_surface_data_to_coupler
    call recv_surface_data_from_coupler

end subroutine process_surface_data_with_coupler_ccpl_interface



subroutine adjust_land_ocn_sice_fraction_ccpl_interface

    use surface_subroutines_with_coupler_mod

    call adjust_land_ocn_sice_fraction

end subroutine adjust_land_ocn_sice_fraction_ccpl_interface
