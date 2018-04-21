# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/restart_mod.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/restart_mod.F90"
!***************************************************************
!  This is a source file of GAMIL, which is reponsible for 
!  retart function of GAMIL with C-Coupler library. This file was 
!  initially finished by Dr. Li Liu. If you have any problem, 
!  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
!***************************************************************


module restart_mod

    use c_coupler_interface_mod
    use register_private_variables_mod

    private exchange_dyn_restart_fields_boundary
    public  do_restart_write
    public  do_restart_read

contains

    subroutine do_restart_read
       use mpi_gamil
       use prognostics
       use comfm1
       use pmgrid, only: beglatex,beglatexdyn,endlatexdyn
       implicit none

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
# 27 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/restart_mod.F90" 2

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
# 28 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/restart_mod.F90" 2
       integer j

       if (.not. c_coupler_is_first_restart_step()) return

       call c_coupler_register_model_data(u3(:,:,:,n3m2),"gamil_2D_decomp_prog","gamil_u3",.true.)
       call c_coupler_register_model_data(v3(:,:,:,n3m2),"gamil_2D_decomp_prog","gamil_v3",.true.)
       call c_coupler_register_model_data(t3(:,:,:,n3m2),"gamil_2D_decomp_prog","gamil_t3",.true.)
       call c_coupler_register_model_data(q3(:,:,:,:,n3m2),"gamil_2D_decomp_prog","gamil_q3",.true.)
       call c_coupler_register_model_data(ps(:,:,n3m2),"gamil_2D_decomp_prog","gamil_ps",.true.)
       call register_phys_dynamic_variables 

       call c_coupler_do_restart_read
       call exchange_dyn_restart_fields_boundary

       call copy_from_phys_io_arrays

       do j=jbeg0, jend0
          nigw(j) = nigw_2D(ibeg0,j)
       end do

       call c_coupler_withdraw_model_data("gamil_2D_decomp_prog","gamil_u3")
       call c_coupler_withdraw_model_data("gamil_2D_decomp_prog","gamil_v3")
       call c_coupler_withdraw_model_data("gamil_2D_decomp_prog","gamil_t3")
       call c_coupler_withdraw_model_data("gamil_2D_decomp_prog","gamil_q3")
       call c_coupler_withdraw_model_data("gamil_2D_decomp_prog","gamil_ps")
       call withdraw_phys_dynamic_variables 

       call c_coupler_advance_timer

    end subroutine do_restart_read



    subroutine do_restart_write
       use mpi_gamil
       use prognostics
       use comfm1
       use pmgrid, only: beglatex,beglatexdyn,endlatexdyn
       implicit none

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
# 68 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/restart_mod.F90" 2

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
# 69 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/restart_mod.F90" 2
       integer :: j


       call c_coupler_register_model_data(u3(:,:,:,n3m2),"gamil_2D_decomp_prog","gamil_u3",.true.)
       call c_coupler_register_model_data(v3(:,:,:,n3m2),"gamil_2D_decomp_prog","gamil_v3",.true.)
       call c_coupler_register_model_data(t3(:,:,:,n3m2),"gamil_2D_decomp_prog","gamil_t3",.true.)
       call c_coupler_register_model_data(q3(:,:,:,:,n3m2),"gamil_2D_decomp_prog","gamil_q3",.true.)
       call c_coupler_register_model_data(ps(:,:,n3m2),"gamil_2D_decomp_prog","gamil_ps",.true.)
       call register_phys_dynamic_variables 

       do j=jbeg0, jend0
          nigw_2D(:,j) = nigw(j)
       end do
       call copy_to_phys_io_arrays
       call c_coupler_do_restart_write()

       call c_coupler_withdraw_model_data("gamil_2D_decomp_prog","gamil_u3")
       call c_coupler_withdraw_model_data("gamil_2D_decomp_prog","gamil_v3")
       call c_coupler_withdraw_model_data("gamil_2D_decomp_prog","gamil_t3")
       call c_coupler_withdraw_model_data("gamil_2D_decomp_prog","gamil_q3")
       call c_coupler_withdraw_model_data("gamil_2D_decomp_prog","gamil_ps")
       call withdraw_phys_dynamic_variables 

    end subroutine do_restart_write


    subroutine exchange_dyn_restart_fields_boundary
       use mpi_gamil
       use comfm1
       use pmgrid, only: beglatexdyn,endlatexdyn
       implicit none

       call gamil_arrays_comm(COMM_TO_LEFT, 1, u(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, v(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, t(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, q(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, ws(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, wpa(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, ghi(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, pes(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, ghs(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, su(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, sv(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, st(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, du(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, dv(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, dtt(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, dps(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, du0(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, dv0(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, dtt0(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, dps0(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, du1(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, dv1(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, dtt1(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, dps1(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, uu(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, vv(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, tt(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, p(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, ply2(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, uk(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, vk(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, ttk(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, psk(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, tb2(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, cb(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, dcb(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, hps(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, c0(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, cb0(:,beglatexdyn,1))

       call gamil_arrays_comm(COMM_TO_RIGHT, 1, u(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, v(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, t(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, q(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, ws(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, wpa(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, ghi(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, pes(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, ghs(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, su(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, sv(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, st(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, du(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, dv(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, dtt(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, dps(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, du0(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, dv0(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, dtt0(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, dps0(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, du1(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, dv1(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, dtt1(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, dps1(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, uu(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, vv(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, tt(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, p(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, ply2(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, uk(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, vk(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, ttk(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, psk(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, tb2(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, cb(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, dcb(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, hps(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, c0(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, cb0(:,beglatexdyn,1))

       call gamil_arrays_comm(COMM_TO_TOP, 1, u(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, v(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, t(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, q(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, ws(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, wpa(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, ghi(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, pes(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_TOP, 1, ghs(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_TOP, 1, su(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, sv(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, st(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, du(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, dv(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, dtt(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, dps(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_TOP, 1, du0(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, dv0(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, dtt0(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, dps0(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_TOP, 1, du1(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, dv1(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, dtt1(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, dps1(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_TOP, 1, uu(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, vv(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, tt(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, p(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_TOP, 1, ply2(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, uk(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, vk(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, ttk(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, psk(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_TOP, 1, tb2(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, cb(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, dcb(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, hps(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_TOP, 1, c0(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_TOP, 1, cb0(:,beglatexdyn,1))

       call gamil_arrays_comm(COMM_TO_BOT, 1, u(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, v(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, t(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, q(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, ws(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, wpa(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, ghi(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, pes(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_BOT, 1, ghs(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_BOT, 1, su(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, sv(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, st(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, du(:,beglatexdyn,1))

       call gamil_arrays_comm(COMM_TO_BOT, 1, dv(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, dtt(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, dps(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_BOT, 1, du0(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, dv0(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, dtt0(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, dps0(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_BOT, 1, du1(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, dv1(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, dtt1(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, dps1(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_BOT, 1, uu(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, vv(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, tt(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, p(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_BOT, 1, ply2(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, uk(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, vk(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, ttk(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, psk(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_BOT, 1, tb2(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, cb(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, dcb(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, hps(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_BOT, 1, c0(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_BOT, 1, cb0(:,beglatexdyn,1))

    end subroutine exchange_dyn_restart_fields_boundary


end module restart_mod
