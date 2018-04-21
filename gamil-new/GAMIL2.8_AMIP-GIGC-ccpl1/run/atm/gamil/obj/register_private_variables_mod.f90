# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/register_private_variables_mod.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/register_private_variables_mod.F90"
!***************************************************************
!  This is a source file of GAMIL, which registers all variables
!  into C-Coupler library for I/O. This file was initially 
!  finished by Dr. Li Liu. If you have any problem, please 
!  contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
!***************************************************************



# 1 "./misc.h" 1
# 10 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/register_private_variables_mod.F90" 2

# 1 "./params.h" 1
# 15 "./params.h"
 
 
# 11 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/register_private_variables_mod.F90" 2


module register_private_variables_mod

    use shr_kind_mod, only: r8 => shr_kind_r8
    use ppgrid
    use phys_grid,    only: read_chunk_from_field, write_field_from_chunk, get_ncols_p
    use pmgrid,       only: masterproc
    use prognostics,  only: ptimelevels, n3, n3m2
    use buffer
    use radae,        only: abstot_3d, absnxt_3d, emstot_3d, initialize_radbuffer
    use comsrf
    use ioFileMod
    use phys_buffer
    use c_coupler_interface_mod

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
# 29 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/register_private_variables_mod.F90" 2
    
    real(r8), allocatable :: qrs_io(:,:,:) 
    real(r8), allocatable :: qrl_io(:,:,:) 
    real(r8), allocatable :: qpert_io(:,:,:) 
    real(r8), allocatable :: cld_n3_io(:,:,:) 
    real(r8), allocatable :: cld_n3m2_io(:,:,:) 
    real(r8), allocatable :: qcwat_n3_io(:,:,:) 
    real(r8), allocatable :: qcwat_n3m2_io(:,:,:) 
    real(r8), allocatable :: tcwat_n3_io(:,:,:) 
    real(r8), allocatable :: tcwat_n3m2_io(:,:,:) 
    real(r8), allocatable :: lcwat_n3_io(:,:,:) 
    real(r8), allocatable :: lcwat_n3m2_io(:,:,:) 

    real(r8), allocatable :: state_wsx_io(:,:)
    real(r8), allocatable :: state_wsy_io(:,:)
    real(r8), allocatable :: state_sst_io(:,:)
    real(r8), allocatable :: state_tref_io(:,:)
    real(r8), allocatable :: state_cflx_io(:,:)
    real(r8), allocatable :: state_lhf_io(:,:)
    real(r8), allocatable :: state_shf_io(:,:)

    real(r8), allocatable :: state_asdir_io(:,:) 
    real(r8), allocatable :: state_asdif_io(:,:) 
    real(r8), allocatable :: state_aldir_io(:,:) 
    real(r8), allocatable :: state_aldif_io(:,:) 

    real(r8), allocatable :: state_lwup_io(:,:) 
    real(r8), allocatable :: state_ts_io(:,:) 
    real(r8), allocatable :: state_tssub_io(:,:,:) 

    real(r8), allocatable :: state_flwds_io(:,:) 
    real(r8), allocatable :: state_sols_io(:,:) 
    real(r8), allocatable :: state_soll_io(:,:) 
    real(r8), allocatable :: state_solsd_io(:,:) 
    real(r8), allocatable :: state_solld_io(:,:) 
    real(r8), allocatable :: state_zbot_io(:,:) 
    real(r8), allocatable :: state_ubot_io(:,:) 
    real(r8), allocatable :: state_vbot_io(:,:) 
    real(r8), allocatable :: state_thbot_io(:,:) 
    real(r8), allocatable :: state_qbot_io(:,:) 
    real(r8), allocatable :: state_pbot_io(:,:) 
    real(r8), allocatable :: state_tbot_io(:,:) 

    real(r8), allocatable :: parm_ts_io(:,:) 
    real(r8), allocatable :: parm_asdir_io(:,:) 
    real(r8), allocatable :: parm_aldir_io(:,:) 
    real(r8), allocatable :: parm_asdif_io(:,:) 
    real(r8), allocatable :: parm_aldif_io(:,:) 
    real(r8), allocatable :: parm_wsx_io(:,:) 
    real(r8), allocatable :: parm_wsy_io(:,:) 
    real(r8), allocatable :: parm_lhf_io(:,:) 
    real(r8), allocatable :: parm_shf_io(:,:) 
    real(r8), allocatable :: parm_lwup_io(:,:) 
    real(r8), allocatable :: parm_cflx_io(:,:) 
    real(r8), allocatable :: parm_tref_io(:,:) 
    real(r8), allocatable :: kvh_io(:,:,:) 
    real(r8), allocatable :: abstot_3d_io(:,:,:,:) 
    real(r8), allocatable :: absnxt_3d_io(:,:,:,:) 
    real(r8), allocatable :: emstot_3d_io(:,:,:) 


# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/physics/cam1/RK_or_MG.h" 1
 
  character(len=2)  RK_or_MG
  parameter (RK_or_MG='MG')

# 90 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/register_private_variables_mod.F90" 2

    !
    ! Public interfaces
    !
    public  copy_to_phys_io_arrays
    public  copy_from_phys_io_arrays
    public  register_phys_dynamic_variables 
    public  withdraw_phys_dynamic_variables 
    public  register_static_variables
    private initialize_phys_io_arrays
    private register_phys_static_variables 
    private register_dyn_variables

    !
    ! Private data
    !
    character(len=256) :: pname  ! Full abs-ems restart filepath
    !
    ! Filename specifier for restart abs-ems file
    ! (%c = caseid, $y = year, $m = month, $d = day, $s = seconds in day, %t = tape number)
    !
    character(len=256) :: rafilename_spec = '%c.cam2.ra.%y-%m-%d-%s'   ! abs-ems restart


CONTAINS

    subroutine initialize_phys_io_arrays

    allocate(qrs_io(pcols,begchunk:endchunk,pver)) 
    allocate(qrl_io(pcols,begchunk:endchunk,pver)) 
    allocate(qpert_io(pcols,begchunk:endchunk,pcnst+pnats)) 
    allocate(cld_n3_io(pcols,begchunk:endchunk,pver)) 
    allocate(cld_n3m2_io(pcols,begchunk:endchunk,pver)) 
    allocate(qcwat_n3_io(pcols,begchunk:endchunk,pver)) 
    allocate(qcwat_n3m2_io(pcols,begchunk:endchunk,pver)) 
    allocate(tcwat_n3_io(pcols,begchunk:endchunk,pver)) 
    allocate(tcwat_n3m2_io(pcols,begchunk:endchunk,pver)) 
    allocate(lcwat_n3_io(pcols,begchunk:endchunk,pver)) 
    allocate(lcwat_n3m2_io(pcols,begchunk:endchunk,pver)) 

    allocate(state_wsx_io(pcols,begchunk:endchunk)) 
    allocate(state_wsy_io(pcols,begchunk:endchunk)) 
    allocate(state_sst_io(pcols,begchunk:endchunk)) 
    allocate(state_tref_io(pcols,begchunk:endchunk)) 
    allocate(state_cflx_io(pcols,begchunk:endchunk)) 
    allocate(state_lhf_io(pcols,begchunk:endchunk)) 
    allocate(state_shf_io(pcols,begchunk:endchunk)) 

    allocate(state_asdir_io(pcols,begchunk:endchunk)) 
    allocate(state_asdif_io(pcols,begchunk:endchunk)) 
    allocate(state_aldir_io(pcols,begchunk:endchunk)) 
    allocate(state_aldif_io(pcols,begchunk:endchunk)) 

    allocate(state_lwup_io(pcols,begchunk:endchunk)) 
    allocate(state_ts_io(pcols,begchunk:endchunk)) 
    allocate(state_tssub_io(pcols,begchunk:endchunk,plevmx)) 

    allocate(state_flwds_io(pcols,begchunk:endchunk)) 
    allocate(state_sols_io(pcols,begchunk:endchunk)) 
    allocate(state_soll_io(pcols,begchunk:endchunk)) 
    allocate(state_solsd_io(pcols,begchunk:endchunk)) 
    allocate(state_solld_io(pcols,begchunk:endchunk)) 
    allocate(state_zbot_io(pcols,begchunk:endchunk)) 
    allocate(state_ubot_io(pcols,begchunk:endchunk)) 
    allocate(state_vbot_io(pcols,begchunk:endchunk)) 
    allocate(state_thbot_io(pcols,begchunk:endchunk)) 
    allocate(state_qbot_io(pcols,begchunk:endchunk)) 
    allocate(state_pbot_io(pcols,begchunk:endchunk)) 
    allocate(state_tbot_io(pcols,begchunk:endchunk)) 
    allocate(parm_ts_io(pcols,begchunk:endchunk)) 
    allocate(parm_asdir_io(pcols,begchunk:endchunk)) 
    allocate(parm_aldir_io(pcols,begchunk:endchunk)) 
    allocate(parm_asdif_io(pcols,begchunk:endchunk)) 
    allocate(parm_aldif_io(pcols,begchunk:endchunk)) 
    allocate(parm_wsx_io(pcols,begchunk:endchunk)) 
    allocate(parm_wsy_io(pcols,begchunk:endchunk)) 
    allocate(parm_lhf_io(pcols,begchunk:endchunk)) 
    allocate(parm_shf_io(pcols,begchunk:endchunk)) 
    allocate(parm_lwup_io(pcols,begchunk:endchunk)) 
    allocate(parm_cflx_io(pcols,begchunk:endchunk)) 
    allocate(parm_tref_io(pcols,begchunk:endchunk)) 
    allocate(kvh_io(pcols,begchunk:endchunk,pverp)) 

    allocate(abstot_3d_io(pcols,begchunk:endchunk,pverp,pverp))
    allocate(absnxt_3d_io(pcols,begchunk:endchunk,pver,4))
    allocate(emstot_3d_io(pcols,begchunk:endchunk,pverp))

    end subroutine initialize_phys_io_arrays


    subroutine copy_to_phys_io_arrays
    implicit none
    integer ichunk,k,m,kvh_idx,ncol
    integer nstep

    nstep = c_coupler_get_nstep()
    
    do ichunk=begchunk,endchunk
       state_wsx_io(:,ichunk) = srfflx_state2d(ichunk)%wsx(:)
       state_wsy_io(:,ichunk) = srfflx_state2d(ichunk)%wsy(:)
       state_sst_io(:,ichunk) = srfflx_state2d(ichunk)%sst(:)
       state_tref_io(:,ichunk) = srfflx_state2d(ichunk)%tref(:)
       state_cflx_io(:,ichunk) = srfflx_state2d(ichunk)%cflx(:,1)
       state_lhf_io(:,ichunk) = srfflx_state2d(ichunk)%lhf(:)
       state_shf_io(:,ichunk) = srfflx_state2d(ichunk)%shf(:)

       state_asdir_io(:,ichunk) = srfflx_state2d(ichunk)%asdir(:)
       state_asdif_io(:,ichunk) = srfflx_state2d(ichunk)%asdif(:)
       state_aldir_io(:,ichunk) = srfflx_state2d(ichunk)%aldir(:)
       state_aldif_io(:,ichunk) = srfflx_state2d(ichunk)%aldif(:)
       state_lwup_io(:,ichunk) = srfflx_state2d(ichunk)%lwup(:)
       state_ts_io(:,ichunk) = srfflx_state2d(ichunk)%ts(:)
       state_flwds_io(:,ichunk) = surface_state2d(ichunk)%flwds(:)
       state_sols_io(:,ichunk) = surface_state2d(ichunk)%sols(:)
       state_soll_io(:,ichunk) = surface_state2d(ichunk)%soll(:)
       state_solsd_io(:,ichunk) = surface_state2d(ichunk)%solsd(:)
       state_solld_io(:,ichunk) = surface_state2d(ichunk)%solld(:)
       state_zbot_io(:,ichunk) = surface_state2d(ichunk)%zbot(:)
       state_ubot_io(:,ichunk) = surface_state2d(ichunk)%ubot(:)
       state_vbot_io(:,ichunk) = surface_state2d(ichunk)%vbot(:)
       state_thbot_io(:,ichunk) = surface_state2d(ichunk)%thbot(:)
       state_qbot_io(:,ichunk) = surface_state2d(ichunk)%qbot(:)
       state_pbot_io(:,ichunk) = surface_state2d(ichunk)%pbot(:)
       state_tbot_io(:,ichunk) = surface_state2d(ichunk)%tbot(:)
       parm_ts_io(:,ichunk) = srfflx_parm2d(ichunk)%ts(:)
       parm_asdir_io(:,ichunk) = srfflx_parm2d(ichunk)%asdir(:)
       parm_aldir_io(:,ichunk) = srfflx_parm2d(ichunk)%aldir(:)
       parm_asdif_io(:,ichunk) = srfflx_parm2d(ichunk)%asdif(:)
       parm_aldif_io(:,ichunk) = srfflx_parm2d(ichunk)%aldif(:)
       parm_wsx_io(:,ichunk) = srfflx_parm2d(ichunk)%wsx(:)
       parm_wsy_io(:,ichunk) = srfflx_parm2d(ichunk)%wsy(:)
       parm_lhf_io(:,ichunk) = srfflx_parm2d(ichunk)%lhf(:)
       parm_shf_io(:,ichunk) = srfflx_parm2d(ichunk)%shf(:)
       parm_lwup_io(:,ichunk) = srfflx_parm2d(ichunk)%lwup(:)
       parm_cflx_io(:,ichunk) = srfflx_parm2d(ichunk)%cflx(:,1)
       parm_tref_io(:,ichunk) = srfflx_parm2d(ichunk)%tref(:)
       do k=1,pver
          qrs_io(:,ichunk,k) = qrs(:,k,ichunk)
          qrl_io(:,ichunk,k) = qrl(:,k,ichunk)
          cld_n3_io(:,ichunk,k) = cld(:,k,ichunk,n3)
          cld_n3m2_io(:,ichunk,k) = cld(:,k,ichunk,n3m2)
          qcwat_n3_io(:,ichunk,k) = qcwat(:,k,ichunk,n3)
          qcwat_n3m2_io(:,ichunk,k) = qcwat(:,k,ichunk,n3m2)
          tcwat_n3_io(:,ichunk,k) = tcwat(:,k,ichunk,n3)
          tcwat_n3m2_io(:,ichunk,k) = tcwat(:,k,ichunk,n3m2)
          lcwat_n3_io(:,ichunk,k) = lcwat(:,k,ichunk,n3)
          lcwat_n3m2_io(:,ichunk,k) = lcwat(:,k,ichunk,n3m2)
       end do
       do k=1,pcnst+pnats
          qpert_io(:,ichunk,k) = qpert(:,k,ichunk)
       end do
       do k=1,plevmx
          state_tssub_io(:,ichunk,k) = surface_state2d(ichunk)%tssub(:,k)
       end do
       if (RK_or_MG=='MG') then
          kvh_idx = pbuf_get_fld_idx('KVH')
          ncol = get_ncols_p(ichunk)
          do k=1,pverp
             kvh_io(:ncol,ichunk,k) = pbuf(kvh_idx)%fld_ptr(1,1:ncol,k,ichunk,1)
          end do
       end if
       if (mod(nstep,iradae).ne.0) then
          do m=1,pverp
          do k=1,pverp
             abstot_3d_io(:,ichunk,k,m)=abstot_3d(:,k,m,ichunk)
          end do
          end do
          do m=1,4
          do k=1,pver
             absnxt_3d_io(:,ichunk,k,m)=absnxt_3d(:,k,m,ichunk)
          end do
          end do
          do k=1,pverp
             emstot_3d_io(:,ichunk,k)=emstot_3d(:,k,ichunk)
          end do
       end if
    end do

    end subroutine copy_to_phys_io_arrays
    
    subroutine copy_from_phys_io_arrays
    implicit none
    integer ichunk,k,m,kvh_idx,ncol
    integer nstep

    nstep = c_coupler_get_nstep()
    
    do ichunk=begchunk,endchunk
       srfflx_state2d(ichunk)%wsx(:) = state_wsx_io(:,ichunk)
       srfflx_state2d(ichunk)%wsy(:) = state_wsy_io(:,ichunk)
       srfflx_state2d(ichunk)%sst(:) = state_sst_io(:,ichunk)
       srfflx_state2d(ichunk)%tref(:) = state_tref_io(:,ichunk)
       srfflx_state2d(ichunk)%cflx(:,1) = state_cflx_io(:,ichunk)
       srfflx_state2d(ichunk)%lhf(:) = state_lhf_io(:,ichunk)
       srfflx_state2d(ichunk)%shf(:) = state_shf_io(:,ichunk)

       srfflx_state2d(ichunk)%asdir(:) = state_asdir_io(:,ichunk)
       srfflx_state2d(ichunk)%asdif(:) = state_asdif_io(:,ichunk)
       srfflx_state2d(ichunk)%aldir(:) = state_aldir_io(:,ichunk)
       srfflx_state2d(ichunk)%aldif(:) = state_aldif_io(:,ichunk)
       srfflx_state2d(ichunk)%lwup(:) = state_lwup_io(:,ichunk)
       srfflx_state2d(ichunk)%ts(:) = state_ts_io(:,ichunk)
       surface_state2d(ichunk)%flwds(:) = state_flwds_io(:,ichunk)
       surface_state2d(ichunk)%sols(:) = state_sols_io(:,ichunk)
       surface_state2d(ichunk)%soll(:) = state_soll_io(:,ichunk)
       surface_state2d(ichunk)%solsd(:) = state_solsd_io(:,ichunk)
       surface_state2d(ichunk)%solld(:) = state_solld_io(:,ichunk)
       surface_state2d(ichunk)%zbot(:) = state_zbot_io(:,ichunk)
       surface_state2d(ichunk)%ubot(:) = state_ubot_io(:,ichunk)
       surface_state2d(ichunk)%vbot(:) = state_vbot_io(:,ichunk)
       surface_state2d(ichunk)%thbot(:) = state_thbot_io(:,ichunk)
       surface_state2d(ichunk)%qbot(:) = state_qbot_io(:,ichunk)
       surface_state2d(ichunk)%pbot(:) = state_pbot_io(:,ichunk)
       surface_state2d(ichunk)%tbot(:) = state_tbot_io(:,ichunk)
       srfflx_parm2d(ichunk)%ts(:) = parm_ts_io(:,ichunk)
       srfflx_parm2d(ichunk)%asdir(:) = parm_asdir_io(:,ichunk)
       srfflx_parm2d(ichunk)%aldir(:) = parm_aldir_io(:,ichunk)
       srfflx_parm2d(ichunk)%asdif(:) = parm_asdif_io(:,ichunk)
       srfflx_parm2d(ichunk)%aldif(:) = parm_aldif_io(:,ichunk)
       srfflx_parm2d(ichunk)%wsx(:) = parm_wsx_io(:,ichunk)
       srfflx_parm2d(ichunk)%wsy(:) = parm_wsy_io(:,ichunk)
       srfflx_parm2d(ichunk)%lhf(:) = parm_lhf_io(:,ichunk)
       srfflx_parm2d(ichunk)%shf(:) = parm_shf_io(:,ichunk)
       srfflx_parm2d(ichunk)%lwup(:) = parm_lwup_io(:,ichunk)
       srfflx_parm2d(ichunk)%cflx(:,1) = parm_cflx_io(:,ichunk)
       srfflx_parm2d(ichunk)%tref(:) = parm_tref_io(:,ichunk)
   
       do k=1,pver
          qrs(:,k,ichunk) = qrs_io(:,ichunk,k)
          qrl(:,k,ichunk) = qrl_io(:,ichunk,k)
          cld(:,k,ichunk,n3) = cld_n3_io(:,ichunk,k)
          cld(:,k,ichunk,n3m2) = cld_n3m2_io(:,ichunk,k)
          qcwat(:,k,ichunk,n3) = qcwat_n3_io(:,ichunk,k)
          qcwat(:,k,ichunk,n3m2) = qcwat_n3m2_io(:,ichunk,k)
          tcwat(:,k,ichunk,n3) = tcwat_n3_io(:,ichunk,k)
          tcwat(:,k,ichunk,n3m2) = tcwat_n3m2_io(:,ichunk,k)
          lcwat(:,k,ichunk,n3) = lcwat_n3_io(:,ichunk,k)
          lcwat(:,k,ichunk,n3m2) = lcwat_n3m2_io(:,ichunk,k)
       end do
       do k=1,pcnst+pnats
          qpert(:,k,ichunk) = qpert_io(:,ichunk,k)
       end do
       do k=1,plevmx
          surface_state2d(ichunk)%tssub(:,k) = state_tssub_io(:,ichunk,k)
       end do
       if (RK_or_MG=='MG') then
          do k=1,pverp
             fld_kvh(:,k,ichunk) = kvh_io(:,ichunk,k)
          end do
       end if
       if (mod(nstep,iradae).ne.0) then
          do m=1,pverp
          do k=1,pverp
             abstot_3d(:,k,m,ichunk)=abstot_3d_io(:,ichunk,k,m)
          end do
          end do
          do m=1,4
          do k=1,pver
             absnxt_3d(:,k,m,ichunk)=absnxt_3d_io(:,ichunk,k,m)
          end do
          end do
          do k=1,pverp
             emstot_3d(:,k,ichunk)=emstot_3d_io(:,ichunk,k)
          end do
       end if
    end do

    end subroutine copy_from_phys_io_arrays


    subroutine register_phys_static_variables 

    call c_coupler_register_model_data(pblht,"gamil_2D_decomp_phys","gamil_pblht",.true.)
    call c_coupler_register_model_data(tpert,"gamil_2D_decomp_phys","gamil_tpert",.true.)
    call c_coupler_register_model_data(qrs_io,"gamil_2D_decomp_phys","gamil_qrs",.true.)
    call c_coupler_register_model_data(qrl_io,"gamil_2D_decomp_phys","gamil_qrl",.true.)
    call c_coupler_register_model_data(qpert_io,"gamil_2D_decomp_phys","gamil_qpert",.true.)
    call c_coupler_register_model_data(cld_n3_io,"gamil_2D_decomp_phys","gamil_cld_n3",.true.)
    call c_coupler_register_model_data(cld_n3m2_io,"gamil_2D_decomp_phys","gamil_cld_n3m2",.true.)
    call c_coupler_register_model_data(qcwat_n3_io,"gamil_2D_decomp_phys","gamil_qcwat_n3",.true.)
    call c_coupler_register_model_data(qcwat_n3m2_io,"gamil_2D_decomp_phys","gamil_qcwat_n3m2",.true.)
    call c_coupler_register_model_data(tcwat_n3_io,"gamil_2D_decomp_phys","gamil_tcwat_n3",.true.)
    call c_coupler_register_model_data(tcwat_n3m2_io,"gamil_2D_decomp_phys","gamil_tcwat_n3m2",.true.)
    call c_coupler_register_model_data(lcwat_n3_io,"gamil_2D_decomp_phys","gamil_lcwat_n3",.true.)
    call c_coupler_register_model_data(lcwat_n3m2_io,"gamil_2D_decomp_phys","gamil_lcwat_n3m2",.true.)
    call c_coupler_register_model_data(fsnt,"gamil_2D_decomp_phys","gamil_fsnt",.true.)
    call c_coupler_register_model_data(fsns,"gamil_2D_decomp_phys","gamil_fsns",.true.)
    call c_coupler_register_model_data(flnt,"gamil_2D_decomp_phys","gamil_flnt",.true.)
    call c_coupler_register_model_data(flns,"gamil_2D_decomp_phys","gamil_flns",.true.)

    call c_coupler_register_model_data(state_wsx_io,"gamil_2D_decomp_phys","gamil_state_wsx",.true.)
    call c_coupler_register_model_data(state_wsy_io,"gamil_2D_decomp_phys","gamil_state_wsy",.true.)
    call c_coupler_register_model_data(state_sst_io,"gamil_2D_decomp_phys","gamil_state_sst",.true.)
    call c_coupler_register_model_data(state_tref_io,"gamil_2D_decomp_phys","gamil_state_tref",.true.)
    call c_coupler_register_model_data(state_cflx_io,"gamil_2D_decomp_phys","gamil_state_cflx",.true.)
    call c_coupler_register_model_data(state_lhf_io,"gamil_2D_decomp_phys","gamil_state_lhf",.true.)
    call c_coupler_register_model_data(state_shf_io,"gamil_2D_decomp_phys","gamil_state_shf",.true.)

    call c_coupler_register_model_data(state_asdir_io,"gamil_2D_decomp_phys","gamil_state_asdir",.true.)
    call c_coupler_register_model_data(state_asdif_io,"gamil_2D_decomp_phys","gamil_state_asdif",.true.)
    call c_coupler_register_model_data(state_aldir_io,"gamil_2D_decomp_phys","gamil_state_aldir",.true.)
    call c_coupler_register_model_data(state_aldif_io,"gamil_2D_decomp_phys","gamil_state_aldif",.true.)
    call c_coupler_register_model_data(asdirice,"gamil_2D_decomp_phys","gamil_asdirice",.true.)
    call c_coupler_register_model_data(asdifice,"gamil_2D_decomp_phys","gamil_asdifice",.true.)
    call c_coupler_register_model_data(aldirice,"gamil_2D_decomp_phys","gamil_aldirice",.true.)
    call c_coupler_register_model_data(aldifice,"gamil_2D_decomp_phys","gamil_aldifice",.true.)
    call c_coupler_register_model_data(tsice,"gamil_2D_decomp_phys","gamil_tsice",.true.)
    call c_coupler_register_model_data(state_lwup_io,"gamil_2D_decomp_phys","gamil_state_lwup",.true.)
    call c_coupler_register_model_data(landfrac,"gamil_2D_decomp_phys","gamil_landfrac",.true.)
    call c_coupler_register_model_data(landm,"gamil_2D_decomp_phys","gamil_landm",.true.)
    call c_coupler_register_model_data(sgh,"gamil_2D_decomp_phys","gamil_sgh",.true.)
    call c_coupler_register_model_data(state_ts_io,"gamil_2D_decomp_phys","gamil_state_ts",.true.)
    call c_coupler_register_model_data(state_tssub_io,"gamil_2D_decomp_phys","gamil_state_tssub",.true.)
    call c_coupler_register_model_data(sicthk,"gamil_2D_decomp_phys","gamil_sicthk",.true.)
    call c_coupler_register_model_data(snowhland,"gamil_2D_decomp_phys","gamil_snowhland",.true.)
    call c_coupler_register_model_data(snowhice,"gamil_2D_decomp_phys","gamil_snowhice",.true.)
    call c_coupler_register_model_data(state_flwds_io,"gamil_2D_decomp_phys","gamil_state_flwds",.true.)
    call c_coupler_register_model_data(state_sols_io,"gamil_2D_decomp_phys","gamil_state_sols",.true.)
    call c_coupler_register_model_data(state_soll_io,"gamil_2D_decomp_phys","gamil_state_soll",.true.)
    call c_coupler_register_model_data(state_solsd_io,"gamil_2D_decomp_phys","gamil_state_solsd",.true.)
    call c_coupler_register_model_data(state_solld_io,"gamil_2D_decomp_phys","gamil_state_solld",.true.)
    call c_coupler_register_model_data(trefmxav,"gamil_2D_decomp_phys","gamil_trefmxav",.true.)
    call c_coupler_register_model_data(trefmnav,"gamil_2D_decomp_phys","gamil_trefmnav",.true.)
    call c_coupler_register_model_data(icefrac,"gamil_2D_decomp_phys","gamil_icefrac",.true.)
    call c_coupler_register_model_data(ocnfrac,"gamil_2D_decomp_phys","gamil_ocnfrac",.true.)
    call c_coupler_register_model_data(state_zbot_io,"gamil_2D_decomp_phys","gamil_state_zbot",.true.)
    call c_coupler_register_model_data(state_ubot_io,"gamil_2D_decomp_phys","gamil_state_ubot",.true.)
    call c_coupler_register_model_data(state_vbot_io,"gamil_2D_decomp_phys","gamil_state_vbot",.true.)
    call c_coupler_register_model_data(state_thbot_io,"gamil_2D_decomp_phys","gamil_state_thbot",.true.)
    call c_coupler_register_model_data(state_qbot_io,"gamil_2D_decomp_phys","gamil_state_qbot",.true.)
    call c_coupler_register_model_data(state_pbot_io,"gamil_2D_decomp_phys","gamil_state_pbot",.true.)
    call c_coupler_register_model_data(state_tbot_io,"gamil_2D_decomp_phys","gamil_state_tbot",.true.)
    call c_coupler_register_model_data(parm_ts_io,"gamil_2D_decomp_phys","gamil_parm_ts",.true.)
    call c_coupler_register_model_data(parm_asdir_io,"gamil_2D_decomp_phys","gamil_parm_asdir",.true.)
    call c_coupler_register_model_data(parm_aldir_io,"gamil_2D_decomp_phys","gamil_parm_aldir",.true.)
    call c_coupler_register_model_data(parm_asdif_io,"gamil_2D_decomp_phys","gamil_parm_asdif",.true.)
    call c_coupler_register_model_data(parm_aldif_io,"gamil_2D_decomp_phys","gamil_parm_aldif",.true.)
    call c_coupler_register_model_data(parm_wsx_io,"gamil_2D_decomp_phys","gamil_parm_wsx",.true.)
    call c_coupler_register_model_data(parm_wsy_io,"gamil_2D_decomp_phys","gamil_parm_wsy",.true.)
    call c_coupler_register_model_data(parm_lhf_io,"gamil_2D_decomp_phys","gamil_parm_lhf",.true.)
    call c_coupler_register_model_data(parm_shf_io,"gamil_2D_decomp_phys","gamil_parm_shf",.true.)
    call c_coupler_register_model_data(parm_lwup_io,"gamil_2D_decomp_phys","gamil_parm_lwup",.true.)
    call c_coupler_register_model_data(parm_cflx_io,"gamil_2D_decomp_phys","gamil_parm_cflx",.true.)
    call c_coupler_register_model_data(parm_tref_io,"gamil_2D_decomp_phys","gamil_parm_tref",.true.)
    call c_coupler_register_model_data(kvh_io,"gamil_2D_decomp_phys","gamil_kvh",.true.)

    call c_coupler_register_model_data(flxave, "NULL", "gamil_flxave", .true.)

    end subroutine register_phys_static_variables 


    subroutine register_phys_dynamic_variables 
    implicit none
    integer nstep

    nstep = c_coupler_get_nstep()

    if (mod(nstep,iradae).ne.0) then
       call c_coupler_register_model_data(abstot_3d_io,"gamil_2D_decomp_phys","gamil_abstot",.true.)
       call c_coupler_register_model_data(absnxt_3d_io,"gamil_2D_decomp_phys","gamil_absnxt",.true.)
       call c_coupler_register_model_data(emstot_3d_io,"gamil_2D_decomp_phys","gamil_emstot",.true.)
    end if

    end subroutine register_phys_dynamic_variables 


    subroutine withdraw_phys_dynamic_variables 
    implicit none
    integer nstep

    nstep = c_coupler_get_nstep()

    if (mod(nstep,iradae).ne.0) then
       call c_coupler_withdraw_model_data("gamil_2D_decomp_phys","gamil_abstot")
       call c_coupler_withdraw_model_data("gamil_2D_decomp_phys","gamil_absnxt")
       call c_coupler_withdraw_model_data("gamil_2D_decomp_phys","gamil_emstot")
    end if

    end subroutine withdraw_phys_dynamic_variables 



    subroutine register_dyn_variables
       use mpi_gamil
       use prognostics
       use comfm1
       use pmgrid, only: beglatex,beglatexdyn,endlatexdyn
       implicit none

# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/comfm2.h" 1

!! (wanhui 2003.07.07)
!! (wanhui 2003.11.04)
!! (b.wang 2004.02.15)




!       real*8  up  (nx,ny,nl)  !
!       real*8  vp  (nx,ny,nl)  !
!       real*8  ttp (nx,ny,nl)  ! variables at step n-1
!       real*8  pps (nx,ny)     !


       real*8  dlt1
       real*8  dlt2

!       real*8  cs0 (nx,ny)
!       real*8  cbs (NX,NY,NL)
!       integer nzad(ny)

       common/comfm2/ dlt1,dlt2
# 479 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/register_private_variables_mod.F90" 2

       call c_coupler_register_model_data(u,"gamil_2D_decomp_dyn","gamil_u",.true.)
       call c_coupler_register_model_data(v,"gamil_2D_decomp_dyn","gamil_v",.true.)
       call c_coupler_register_model_data(t,"gamil_2D_decomp_dyn","gamil_t",.true.)
       call c_coupler_register_model_data(q,"gamil_2D_decomp_dyn","gamil_q",.true.)
       call c_coupler_register_model_data(ws,"gamil_2D_decomp_dyn","gamil_ws",.true.)
       call c_coupler_register_model_data(wpa,"gamil_2D_decomp_dyn","gamil_wpa",.true.)
       call c_coupler_register_model_data(ghi,"gamil_2D_decomp_dyn","gamil_ghi",.true.)
       call c_coupler_register_model_data(pes,"gamil_2D_decomp_dyn","gamil_pes",.true.)
       call c_coupler_register_model_data(ghs,"gamil_2D_decomp_dyn","gamil_ghs",.true.)
       call c_coupler_register_model_data(su,"gamil_2D_decomp_dyn","gamil_su",.true.)
       call c_coupler_register_model_data(sv,"gamil_2D_decomp_dyn","gamil_sv",.true.)
       call c_coupler_register_model_data(st,"gamil_2D_decomp_dyn","gamil_st",.true.)
       call c_coupler_register_model_data(du,"gamil_2D_decomp_dyn","gamil_du",.true.)
       call c_coupler_register_model_data(dv,"gamil_2D_decomp_dyn","gamil_dv",.true.)
       call c_coupler_register_model_data(dtt,"gamil_2D_decomp_dyn","gamil_dtt",.true.)
       call c_coupler_register_model_data(dps,"gamil_2D_decomp_dyn","gamil_dps",.true.)
       call c_coupler_register_model_data(du0,"gamil_2D_decomp_dyn","gamil_du0",.true.)
       call c_coupler_register_model_data(dv0,"gamil_2D_decomp_dyn","gamil_dv0",.true.)
       call c_coupler_register_model_data(dtt0,"gamil_2D_decomp_dyn","gamil_dtt0",.true.)
       call c_coupler_register_model_data(dps0,"gamil_2D_decomp_dyn","gamil_dps0",.true.)
       call c_coupler_register_model_data(du1,"gamil_2D_decomp_dyn","gamil_du1",.true.)
       call c_coupler_register_model_data(dv1,"gamil_2D_decomp_dyn","gamil_dv1",.true.)
       call c_coupler_register_model_data(dtt1,"gamil_2D_decomp_dyn","gamil_dtt1",.true.)
       call c_coupler_register_model_data(dps1,"gamil_2D_decomp_dyn","gamil_dps1",.true.)
       call c_coupler_register_model_data(uu,"gamil_2D_decomp_dyn","gamil_uu",.true.)
       call c_coupler_register_model_data(vv,"gamil_2D_decomp_dyn","gamil_vv",.true.)
       call c_coupler_register_model_data(tt,"gamil_2D_decomp_dyn","gamil_tt",.true.)
       call c_coupler_register_model_data(p,"gamil_2D_decomp_dyn","gamil_p",.true.)
       call c_coupler_register_model_data(ply2,"gamil_2D_decomp_dyn","gaiml_ply2",.true.)
       call c_coupler_register_model_data(uk,"gamil_2D_decomp_dyn","gamil_uk",.true.)
       call c_coupler_register_model_data(vk,"gamil_2D_decomp_dyn","gamil_vk",.true.)
       call c_coupler_register_model_data(ttk,"gamil_2D_decomp_dyn","gamil_ttk",.true.)
       call c_coupler_register_model_data(psk,"gamil_2D_decomp_dyn","gamil_psk",.true.)
       call c_coupler_register_model_data(tb2,"gamil_2D_decomp_dyn","gamil_tb2",.true.)
       call c_coupler_register_model_data(cb,"gamil_2D_decomp_dyn","gamil_cb",.true.)
       call c_coupler_register_model_data(dcb,"gamil_2D_decomp_dyn","gamil_dcb",.true.)
       call c_coupler_register_model_data(hps,"gamil_2D_decomp_dyn","gamil_hps",.true.)
       call c_coupler_register_model_data(c0,"gamil_2D_decomp_dyn","gamil_c0",.true.)
       call c_coupler_register_model_data(cb0,"gamil_2D_decomp_dyn","gamil_cb0",.true.)
       call c_coupler_register_model_data(phis,"gamil_2D_decomp_prog","gamil_phis",.true.)
       call c_coupler_register_model_data(omga,"gamil_2D_decomp_prog","gamil_omga",.true.)
       call c_coupler_register_model_data(t31,"gamil_2D_decomp_prog","gamil_t31",.true.)
       call c_coupler_register_model_data(q31,"gamil_2D_decomp_prog","gamil_q31",.true.)
       call c_coupler_register_model_data(nigw_2D,"gamil_2D_decomp_dyn","gamil_nigw",.true.)
       call c_coupler_register_model_data(itime, "NULL", "gamil_itime", .true.)
       call c_coupler_register_model_data(dlt1, "NULL", "gamil_dlt1", .true.)
       call c_coupler_register_model_data(dlt2, "NULL", "gamil_dlt2", .true.)

    end subroutine register_dyn_variables



    subroutine register_static_variables

       implicit none

       call initialize_phys_io_arrays
       call register_phys_static_variables 
       call register_dyn_variables

    end subroutine register_static_variables


end module register_private_variables_mod
