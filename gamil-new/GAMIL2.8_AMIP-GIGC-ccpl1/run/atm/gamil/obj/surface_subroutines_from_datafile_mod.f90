# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/surface_subroutines_from_datafile_mod.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/surface_subroutines_from_datafile_mod.F90"
!***************************************************************
!  This is a source file of GAMIL, containing all subroutines 
!  which get surface data from data files. This file was 
!  initially finished by Dr. Li Liu, through reorganizing 
!  existing source files of NCAR. If you have any problem, 
!  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
!***************************************************************



# 1 "./misc.h" 1
# 11 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/surface_subroutines_from_datafile_mod.F90" 2

module surface_subroutines_from_datafile_mod

    use c_coupler_interface_mod
    use ppgrid,         only: begchunk, endchunk
    use pmgrid

    public calculate_ocn_albedo
    public calculate_sice_albedo
    public initialize_online_lsm
    public initialize_data_ocn
    public initialize_data_sice
    public merge_ts

    integer,parameter,private :: R8  = selected_real_kind(12) ! 8 byte real
    real(r8), private, allocatable :: tssav(:,:) ! cam surface temperatures

contains

subroutine calculate_ocn_albedo

    use c_coupler_interface_mod
    use phys_grid,      only: get_ncols_p, get_rlat_all_p, get_rlon_all_p
    use comsrf

    implicit none

    integer  lchnk
    integer  ncol               ! number of columns in current chunk
    real(r8) coszrs(pcols)      ! Cosine solar zenith angle
    real(r8) clat1(pcols)       ! Current latitude(radians)
    real(r8) clon1(pcols)       ! Current longitude(radians)
    real(r8) calday             ! current calendar day

    call c_coupler_get_current_calendar_time(calday)

    if (c_coupler_is_first_step()) then
        do lchnk = begchunk, endchunk
            ncol = get_ncols_p(lchnk)
            call get_rlat_all_p(lchnk, ncol, clat1)
            call get_rlon_all_p(lchnk, ncol, clon1)
            call zenith(calday, clat1, clon1, coszrs, ncol)
            call albocean(lchnk, ncol, coszrs, &
                          srfflx_parm2d_ocn(lchnk)%asdir, srfflx_parm2d_ocn(lchnk)%aldir, &
                          srfflx_parm2d_ocn(lchnk)%asdif, srfflx_parm2d_ocn(lchnk)%aldif)
        end do

    end if

end subroutine calculate_ocn_albedo



subroutine calculate_sice_albedo

    use c_coupler_interface_mod
    use phys_grid,      only: get_ncols_p, get_rlat_all_p, get_rlon_all_p
    use comsrf

    implicit none

    integer  lchnk
    integer  ncol               ! number of columns in current chunk
    real(r8) coszrs(pcols)      ! Cosine solar zenith angle
    real(r8) clat1(pcols)       ! Current latitude(radians)
    real(r8) clon1(pcols)       ! Current longitude(radians)
    real(r8) calday             ! current calendar day

    call c_coupler_get_current_calendar_time(calday)
    if (c_coupler_is_first_step()) then
        do lchnk = begchunk, endchunk
            ncol = get_ncols_p(lchnk)
            call get_rlat_all_p(lchnk, ncol, clat1)
            call get_rlon_all_p(lchnk, ncol, clon1)
            call zenith (calday, clat1, clon1, coszrs, ncol)
            call albice(lchnk, ncol, tsice(1,lchnk), snowhice(1,lchnk), coszrs, &
                        srfflx_parm2d_sice(lchnk)%asdir, srfflx_parm2d_sice(lchnk)%aldir, &
                        srfflx_parm2d_sice(lchnk)%asdif, srfflx_parm2d_sice(lchnk)%aldif)
            !
            ! fill in ice albedoes for therm ice model
            !
            asdirice(:ncol,lchnk)= srfflx_parm2d_sice(lchnk)%asdir(:ncol)
            aldirice(:ncol,lchnk)= srfflx_parm2d_sice(lchnk)%aldir(:ncol)
            asdifice(:ncol,lchnk)= srfflx_parm2d_sice(lchnk)%asdif(:ncol)
            aldifice(:ncol,lchnk)= srfflx_parm2d_sice(lchnk)%aldif(:ncol)
        end do
    end if

end subroutine calculate_sice_albedo


subroutine read_land_inidat
      use pmgrid
      use rgrid, only: nlon
      use comsrf
!!    use commap, only: w
!!    use physconst, only: gravit
      use phys_grid
      use buffer

      implicit none

      include 'netcdf.inc'

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
# 115 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/surface_subroutines_from_datafile_mod.F90" 2

# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/control/comhyb.h" 1
!----------------------------------------------------------------------- 
! 
! Purpose: Hybrid level definitions: p = a*p0 + b*ps
!          interfaces   p(k) = hyai(k)*ps0 + hybi(k)*ps
!          midpoints    p(k) = hyam(k)*ps0 + hybm(k)*ps
! 
!-----------------------------------------------------------------------
!!!  vertical level definitions in LASG dynamical core: p = pes*sigma + pt
!!!        interfaces   ply(k) = ps*sig (k) + pmtop
!!!        midpoints    ply(k) = ps*sigl(k) + pmtop
!!!---------------------------------------------------------------------
!!!(wanhui 2003.04.30)
!!!(wanhui 2003.10.23)  (std.atm. variables removed)

      real(r8) hyai(plevp)       ! ps0 component of hybrid coordinate - interfaces
      real(r8) hybi(plevp)       ! ps component of hybrid coordinate - interfaces
      real(r8) hyam(plev)        ! ps0 component of hybrid coordinate - midpoints
      real(r8) hybm(plev)        ! ps component of hybrid coordinate - midpoints

!!    real(r8) hybd(plev)        ! difference  in b (hybi) across layers
      real(r8) hypi(plevp)       ! reference pressures at interfaces
      real(r8) hypm(plev)        ! reference pressures at midpoints
!!    real(r8) hypd(plev)        ! reference pressure layer thickness

      real(r8) ps0         ! base state sfc pressure for level definitions
!!    real(r8) psr         ! reference surface pressure for linearization
!!    real(r8) prsfac      ! log pressure extrapolation factor (time, space independent)

!!    integer nprlev       ! number of pure pressure levels at top

      real(r8) :: pmtop              !
      real(r8) :: sig (plevp)        !
      real(r8) :: sigl(plev)         !  fm2003 VPAR variables
      real(r8) :: dsig(plev)         !

!!(wanhui 2003.10.23)
!!------------------------------------------------------------
!!    real(r8) :: tbb (plevstd)         !
!!    real(r8) :: hbb (plevstd)         !
!!    real(r8) :: cbb (plevstd)         !
!!    real(r8) :: dcbb(plevstd)         !  fm2003 std. atm.
!!    real(r8) :: p00, t00              !
!!    real(r8) :: psb (plond,plat)      !
!!    real(r8) :: tsb (plond,plat)      !
!!------------------------------------------------------------
!!(2003.10.23)(these variables are in module stdatm now)


!!      common /comhyb/ hyai ,hyam  ,hybi ,hybm
!!      common /comhyb/ hybd ,hypi ,hypm  ,hypd
!!      common /comhyb/ ps0         ,psr         ,prsfac      ,nprlev

      common /comhyb/ hyai ,hybi ,hyam ,hybm
      common /comhyb/ hypi ,hypm
      common /comhyb/ ps0  ,pmtop
      common /comhyb/ sig  ,sigl , dsig
!!    common /comhyb/ tbb  ,hbb  , cbb  ,dcbb  ,p00 ,t00 ,psb ,tsb    !!(wh 2003.10.23)
 
# 116 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/surface_subroutines_from_datafile_mod.F90" 2
!!#include <comqfl.h>

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
# 118 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/surface_subroutines_from_datafile_mod.F90" 2

# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/control/perturb.h" 1
!----------------------------------------------------------------------- 
! 
! Purpose: Perturbation information
! 
! Author: Jim Rosinski
! 
!-----------------------------------------------------------------------
      common /perturb/ pertlim
      real(r8) pertlim        ! Bound for rectangular distribution of perturbation
!
 
# 119 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/surface_subroutines_from_datafile_mod.F90" 2
!
      integer landfracid        ! Variable ID's
      integer tsid, ts1id, ts2id, ts3id, ts4id ,tsiceid! Variable ID's
      integer snowhiceid        ! Variable ID's
      integer i, j
      integer strt2d(3)         ! start lon, lat, time indices for netcdf 2-d
      data strt2d/3*1/          ! Only index 2 will ever change
      integer cnt2d(3)          ! lon, lat, time counts for netcdf 2-d
      data cnt2d/plon,1,1/      ! 2-d arrs: Always grab only a "plon" slice
      real(r8), allocatable :: tmpchunk3d(:,:,:)
      real(r8), allocatable :: tmpchunk(:,:)
      real(r8), allocatable :: landfrac_tmp(:,:)
      real(r8), allocatable :: ts_tmp(:,:)
      real(r8), allocatable :: tsice_tmp(:,:)
      real(r8), allocatable :: tssub_tmp(:,:,:)
      real(r8), allocatable :: snowhice_tmp(:,:)


      allocate ( landfrac_tmp(plond,plat) )
      allocate ( tmpchunk(pcols,begchunk:endchunk) )
      allocate ( tmpchunk3d(pcols,plevmx,begchunk:endchunk) )
      allocate ( ts_tmp(plond,plat) )
      allocate ( tsice_tmp(plond,plat) )
      allocate ( tssub_tmp(plond,plevmx,plat) )
      allocate ( snowhice_tmp(plond,plat) )

      if (masterproc) then
!
! For land-fraction check if the variable name LANDFRAC is on the dataset if not assume FLAND
!
         if ( nf_inq_varid(ncid_ini, 'LANDFRAC', landfracid ) == NF_NOERR ) then
            call wrap_inq_varid (ncid_ini, 'LANDFRAC', landfracid)
         else if ( nf_inq_varid(ncid_ini, 'FLAND', landfracid ) == NF_NOERR ) then
            call wrap_inq_varid (ncid_ini, 'FLAND', landfracid)
         else 
            write(6,*) "WARNING: land fraction is not given by the initial input file"
            landfracid = -999
         end if
         if (landfracid .ne. -999) then
            call wrap_inq_varid (ncid_ini, 'TS', tsid)
            call wrap_inq_varid (ncid_ini, 'TSICE', tsiceid)
            call wrap_inq_varid (ncid_ini, 'TS1', ts1id)
            call wrap_inq_varid (ncid_ini, 'TS2', ts2id)
            call wrap_inq_varid (ncid_ini, 'TS3', ts3id)
            call wrap_inq_varid (ncid_ini, 'TS4', ts4id)
            call wrap_inq_varid (ncid_ini, 'SNOWHICE', snowhiceid)
            do j=1,plat
               strt2d(2) = j
               if (aqua_planet) then
                  do i=1,nlon(j)
                     landfrac_tmp(i,j) = 0.
                  end do
               else
                  call wrap_get_vara_realx (ncid_ini, landfracid, strt2d, cnt2d, landfrac_tmp(1,j))
               endif
               call wrap_get_vara_realx (ncid_ini, tsid, strt2d, cnt2d, ts_tmp(1,j))
               call wrap_get_vara_realx (ncid_ini, tsiceid, strt2d, cnt2d, tsice_tmp(1,j))
               call wrap_get_vara_realx (ncid_ini, ts1id, strt2d, cnt2d, tssub_tmp(1,1,j))
               call wrap_get_vara_realx (ncid_ini, ts2id, strt2d, cnt2d, tssub_tmp(1,2,j))
               call wrap_get_vara_realx (ncid_ini, ts3id, strt2d, cnt2d, tssub_tmp(1,3,j))
               call wrap_get_vara_realx (ncid_ini, ts4id, strt2d, cnt2d, tssub_tmp(1,4,j))
               call wrap_get_vara_realx(ncid_ini, snowhiceid, strt2d, cnt2d, snowhice_tmp(1,j))
            end do
         end if
      end if

      call scatter_field_to_chunk(1,1,1,plond,landfrac_tmp,landfrac(1,begchunk))
      call scatter_field_to_chunk(1,1,1,plond,ts_tmp,tmpchunk)
      call scatter_field_to_chunk(1,1,1,plond,tsice_tmp,tsice(1,begchunk))
      call scatter_field_to_chunk(1,1,1,plond,snowhice_tmp,snowhice(1,begchunk))
      call scatter_field_to_chunk(1,plevmx,1,plond,tssub_tmp,tmpchunk3d)
      do i =begchunk,endchunk
         surface_state2d(i)%tssub(:,:) = tmpchunk3d(:,:,i)
      end do

      deallocate ( tmpchunk3d)
      deallocate ( tmpchunk )
      deallocate ( landfrac_tmp )
      deallocate ( ts_tmp )
      deallocate ( tsice_tmp )
      deallocate ( tssub_tmp )
      deallocate ( snowhice_tmp )

end subroutine read_land_inidat



subroutine initialize_online_lsm

    use atm_lndMod,     only: atmlnd_ini
    use comsrf
    use c_coupler_interface_mod

    implicit none

    integer lchnk


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
# 217 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/surface_subroutines_from_datafile_mod.F90" 2
    allocate(tssav(pcols,begchunk:endchunk))
    if (.not. adiabatic .and. .not. ideal_phys .and. .not. aqua_planet) then
        call srfflx_parm_reset(srfflx_parm2d_lsm)
        call atmlnd_ini(srfflx_parm2d_lsm)
        do lchnk = begchunk, endchunk
            tssav(:,lchnk) = srfflx_parm2d_lsm(lchnk)%ts(:)
        end do

    end if



end subroutine initialize_online_lsm



subroutine read_SOM_inidat
      use pmgrid
      use rgrid, only: nlon
      use comsrf
      use phys_grid
      use buffer

      implicit none

      include 'netcdf.inc'

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
# 244 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/surface_subroutines_from_datafile_mod.F90" 2

# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/control/comhyb.h" 1
!----------------------------------------------------------------------- 
! 
! Purpose: Hybrid level definitions: p = a*p0 + b*ps
!          interfaces   p(k) = hyai(k)*ps0 + hybi(k)*ps
!          midpoints    p(k) = hyam(k)*ps0 + hybm(k)*ps
! 
!-----------------------------------------------------------------------
!!!  vertical level definitions in LASG dynamical core: p = pes*sigma + pt
!!!        interfaces   ply(k) = ps*sig (k) + pmtop
!!!        midpoints    ply(k) = ps*sigl(k) + pmtop
!!!---------------------------------------------------------------------
!!!(wanhui 2003.04.30)
!!!(wanhui 2003.10.23)  (std.atm. variables removed)

      real(r8) hyai(plevp)       ! ps0 component of hybrid coordinate - interfaces
      real(r8) hybi(plevp)       ! ps component of hybrid coordinate - interfaces
      real(r8) hyam(plev)        ! ps0 component of hybrid coordinate - midpoints
      real(r8) hybm(plev)        ! ps component of hybrid coordinate - midpoints

!!    real(r8) hybd(plev)        ! difference  in b (hybi) across layers
      real(r8) hypi(plevp)       ! reference pressures at interfaces
      real(r8) hypm(plev)        ! reference pressures at midpoints
!!    real(r8) hypd(plev)        ! reference pressure layer thickness

      real(r8) ps0         ! base state sfc pressure for level definitions
!!    real(r8) psr         ! reference surface pressure for linearization
!!    real(r8) prsfac      ! log pressure extrapolation factor (time, space independent)

!!    integer nprlev       ! number of pure pressure levels at top

      real(r8) :: pmtop              !
      real(r8) :: sig (plevp)        !
      real(r8) :: sigl(plev)         !  fm2003 VPAR variables
      real(r8) :: dsig(plev)         !

!!(wanhui 2003.10.23)
!!------------------------------------------------------------
!!    real(r8) :: tbb (plevstd)         !
!!    real(r8) :: hbb (plevstd)         !
!!    real(r8) :: cbb (plevstd)         !
!!    real(r8) :: dcbb(plevstd)         !  fm2003 std. atm.
!!    real(r8) :: p00, t00              !
!!    real(r8) :: psb (plond,plat)      !
!!    real(r8) :: tsb (plond,plat)      !
!!------------------------------------------------------------
!!(2003.10.23)(these variables are in module stdatm now)


!!      common /comhyb/ hyai ,hyam  ,hybi ,hybm
!!      common /comhyb/ hybd ,hypi ,hypm  ,hypd
!!      common /comhyb/ ps0         ,psr         ,prsfac      ,nprlev

      common /comhyb/ hyai ,hybi ,hyam ,hybm
      common /comhyb/ hypi ,hypm
      common /comhyb/ ps0  ,pmtop
      common /comhyb/ sig  ,sigl , dsig
!!    common /comhyb/ tbb  ,hbb  , cbb  ,dcbb  ,p00 ,t00 ,psb ,tsb    !!(wh 2003.10.23)
 
# 245 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/surface_subroutines_from_datafile_mod.F90" 2

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
# 246 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/surface_subroutines_from_datafile_mod.F90" 2

# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/control/perturb.h" 1
!----------------------------------------------------------------------- 
! 
! Purpose: Perturbation information
! 
! Author: Jim Rosinski
! 
!-----------------------------------------------------------------------
      common /perturb/ pertlim
      real(r8) pertlim        ! Bound for rectangular distribution of perturbation
!
 
# 247 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/surface_subroutines_from_datafile_mod.F90" 2
!
      integer strt2d(3)         ! start lon, lat, time indices for netcdf 2-d
      data strt2d/3*1/          ! Only index 2 will ever change
      integer cnt2d(3)          ! lon, lat, time counts for netcdf 2-d
      data cnt2d/plon,1,1/      ! 2-d arrs: Always grab only a "plon" slice
      integer j
      integer sicid
      real(r8), allocatable :: sicthk_tmp(:,:)
!
! Set sea-ice thickness and snow cover:
!
      allocate ( sicthk_tmp(plond,plat) )

      if (masterproc) then
         call wrap_inq_varid (ncid_ini, 'SICTHK', sicid)
         do j=1,plat
            strt2d(2) = j
            call wrap_get_vara_realx(ncid_ini, sicid, strt2d, cnt2d, sicthk_tmp(1,j))
         end do
      end if
      call scatter_field_to_chunk(1,1,1,plond,sicthk_tmp,sicthk(1,begchunk))

      deallocate ( sicthk_tmp )

end subroutine read_SOM_inidat



subroutine initialize_SOM
    use c_coupler_interface_mod
    use pmgrid
    use mpishorthand

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
# 282 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/surface_subroutines_from_datafile_mod.F90" 2
    include 'netcdf.inc'

    integer  ret                ! NetCDF returned status
    integer  attlen             ! NetCDF attribute length
    character(256) text         ! NetCDF attribute
    integer  sghid              ! NetCDF sgh field id
    logical  oro_hires          ! true => ORO came from high res topo file


    call read_SOM_inidat
    !
    !----------------------------------------------------------------------
    ! 3. Determine if SGH field came from hi-res dataset
    !----------------------------------------------------------------------
    !
    if (c_coupler_is_first_step()) then
        if (masterproc) then
            call wrap_inq_varid(ncid_ini, 'SGH', sghid)
            ret = nf_inq_attlen(ncid_ini, sghid, 'from_hires', attlen)
            if (ret == nf_noerr .and. attlen > 256) then
                write(6, "('Error: initext: attribute length of ""from_hires"" is too long')")
                call endrun
            end if
            ret = nf_get_att_text(ncid_ini, sghid, 'from_hires', text)
            if (ret == nf_noerr .and. text(1:4) == 'true') then
                oro_hires = .true.
                write(6, "('Notice: initext: attribute ""from_hires"" is true')")
                write(6, "('Notice: initext: ""tssub"" will be used to guess sea ice')")
            else
                oro_hires = .false.
                write(6, "('Notice: initext: attribute ""from_hires"" is either false or not present')")
                write(6, "('Notice: initext: where sea ice exists, its initial temperature will be just below freezing')")
            end if
        end if
        call mpibcast(oro_hires, 1, mpilog, 0, mpicom)
    end if
    !
    ! Slab ocean model: set initial surf temps for initial run. Read in 2 time slices of
    ! mixed layer depths and q fluxes from boundary dataset whether initial or restart
    !
    !call somini(oro_hires)   ! We do not have slab ocean currently

    call calculate_ocn_albedo
    call calculate_ocn_albedo

end subroutine initialize_SOM



subroutine initialize_data_ocn

    use pmgrid
    use filenames,      only: bndtvo, bndtvs
    use comsrf
    use sst_data,       only: sstini, sstint, sstan
    use ice_data,       only: iceini, iceint
    use ioFileMod

    use ice_constants,  only: Tffresh


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
# 346 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/surface_subroutines_from_datafile_mod.F90" 2

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
# 347 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/surface_subroutines_from_datafile_mod.F90" 2
    character(256) locfn        ! netcdf local filename to open
    integer        lchnk

    if (.not. adiabatic .and. .not. ideal_phys) then
        call srfflx_parm_reset(srfflx_parm2d_ocn)
        if (ncid_sst .eq. -1 .and. masterproc) then
            call getfil(bndtvs, locfn)
            call wrap_open(locfn, 0, ncid_sst)
            write(6, "('Notice: initext: ')", advance="no")
            write(6, "('wrap_open returns ncid ', I5)", advance="no") ncid_sst
            write(6, "(' for file ', A)") trim(locfn)
        end if
        call sstini
        call sstint
        do lchnk = begchunk, endchunk
            srfflx_parm2d_ocn(lchnk)%ts(:) = srfflx_state2d(lchnk)%sst(:)
        end do
    end if

end subroutine initialize_data_ocn



subroutine initialize_data_sice

    use pmgrid
    use filenames,      only: bndtvo, bndtvs
    use comsrf
    use sst_data,       only: sstini, sstint, sstan
    use ice_data,       only: iceini, iceint
    use ioFileMod

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
# 382 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/surface_subroutines_from_datafile_mod.F90" 2

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
# 383 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/surface_subroutines_from_datafile_mod.F90" 2
    character(256) locfn        ! netcdf local filename to open
    integer        lchnk

    if (.not. adiabatic .and. .not. ideal_phys) then
        call srfflx_parm_reset(srfflx_parm2d_sice)
        if (ncid_sst .eq. -1 .and. masterproc) then
            call getfil(bndtvs, locfn)
            call wrap_open(locfn, 0, ncid_sst)
            write(6, "('Notice: initext: ')", advance="no")
            write(6, "('wrap_open returns ncid ', I5)", advance="no") ncid_sst
            write(6, "(' for file ', A)") trim(locfn)
        end if
        call iceini
        call iceint
        do lchnk = begchunk, endchunk
            srfflx_parm2d_sice(lchnk)%ts(:) = tsice(:,lchnk) 
        end do
    end if


end subroutine initialize_data_sice



subroutine merge_ts
    use comsrf
    use physconst,      only: stebol

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
# 414 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/surface_subroutines_from_datafile_mod.F90" 2

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
# 415 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/surface_subroutines_from_datafile_mod.F90" 2

    integer i, lchnk, ncol

    if (c_coupler_is_first_step()) then
        call update_srf_fluxes(srfflx_state2d, srfflx_parm2d_lsm, landfrac, .true.)
        call update_srf_fluxes(srfflx_state2d, srfflx_parm2d_ocn, ocnfrac, .true.)
        call update_srf_fluxes(srfflx_state2d,srfflx_parm2d_sice,icefrac,.true.)

        if (.not. adiabatic .and. .not. ideal_phys) then
            do lchnk = begchunk, endchunk
                ncol = get_ncols_p(lchnk)
                do i = 1, ncol
                    if (landfrac(i,lchnk).ne.1.) then
                        srfflx_state2d(lchnk)%lwup(i) = &
                            stebol*(srfflx_state2d(lchnk)%ts(i)**4)
                    end if
                end do
            end do
        end if
    end if
end subroutine merge_ts


end module surface_subroutines_from_datafile_mod



