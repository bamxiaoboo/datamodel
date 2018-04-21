# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/physics/cam1/prescribed_aerosols.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/physics/cam1/prescribed_aerosols.F90"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/physics/cam1/prescribed_aerosols.F90" 2

# 1 "./params.h" 1
# 15 "./params.h"
 
 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/physics/cam1/prescribed_aerosols.F90" 2

module prescribed_aerosols
!----------------------------------------------------------------------- 
! 
! Purposes: 
!       read, store, interpolate, and return fields
!         of aerosols to CAM.  The initialization
!         file (mass.nc) is assumed to be a monthly climatology
!         of aerosols from MATCH (on a sigma pressure
!         coordinate system).
!       also provide a "background" aerosol field to correct
!         for any deficiencies in the physical parameterizations
!         This fields is a "tuning" parameter.
!       Public methods:
!       (1) - initialization
!          read aerosol masses from external file
!             also pressure coordinates
!          convert from monthly average values to mid-month values
!       (2) - interpolation (time and vertical)
!          interpolate onto pressure levels of CAM
!          interpolate to time step of CAM
!          return mass of aerosols 
!          provide background aerosol based on "tauback"
!       (3) - diagnostics
!          write out various diagnostic fields
! 
! Calling Heirarchy
!   aerosol_initialize (public)
!   aerosol_mass_get (public)
!     vert_interpolate
!     background
!     scale_aerosols
!   aerosol_indirect (public)
!
!-----------------------------------------------------------------------

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use shr_const_mod,  only: pie => SHR_CONST_PI
  !use spmd_utils,     only: masterproc     ! in gamil masterproc belong to module pmgrid
  USE pmgrid,      ONLY: masterproc,iam     !! for debug-test 
  use ppgrid,         only: pcols, pver, pverp, begchunk, endchunk
  use phys_grid,      only: get_ncols_p, scatter_field_to_chunk
  use infnan,         only: inf, bigint
  !use abortutils,     only: endrun         ! in gamil , endrun can be called directly
  use aerosol_index,  only: naer, naer_all, aerosol_name, wgt_sscm, &
                            idxSUL,idxSSLTA, idxSSLTC, idxOCPHO, idxBCPHO, &
                            idxOCPHI, idxBCPHI, idxBG, idxVOLC, &
                            idxDUSTfirst, numDUST, idxCARBONfirst, numCARBON, &
                            idxSSLTfirst, numSSLT
  use scamMod,        only: scmlon,scmlat,single_column
  use error_messages, only: handle_ncerr
  use physics_types,  only: physics_state
  use boundarydata,   only: boundarydata_init, boundarydata_type
  use c_coupler_interface_mod
  !use perf_mod     ! such as t_start  sxj---

   use mpishorthand

!  use netcdf
  implicit none
   
  include 'netcdf.inc'  !!!!
!
! this boundarydata_type is used for datasets in the ncols format only.
!
  type(boundarydata_type) :: aerosol_datan

! "background" aerosol species mmr
  real(r8), public :: tauback = 0._r8

! 
! default values for namelist variables
!

! compute carbons scaling from population
  character(len=16), public :: scenario_carbon_scale   = 'FIXED' ! 'FIXED' or 'RAMPED'
  character(len=16), public :: scenario_prescribed_sulfur   = 'FIXED' ! 'FIXED' or 'RAMPED'
  character(len=16), public :: prescribed_sulfur   = 'direct' ! 'off', 'passive' or 'direct'
  integer, public :: rampyear_prescribed_sulfur   = bigint ! this is the only acceptable value

  real(r8), public :: density_aer(naer_all) ! density of aerosol (kg/m3)
  real(r8), public :: hygro_aer(naer_all) ! hygroscopicity of aerosol
  real(r8), public :: dryrad_aer(naer_all) ! number mode radius (m) of aerosol size distribution
  real(r8), public :: dispersion_aer(naer_all) ! geometric standard deviation of aerosol size distribution
  real(r8), public :: num_to_mass_aer(naer_all) ! ratio of number concentration to mass concentration (#/kg)
  character*8, public :: aername(naer_all)

  private
  save

  integer :: aernid = -1           ! netcdf id for aerosol file (init to invalid)
  integer :: species_id(naer) = -1 ! netcdf_id of each aerosol species (init to invalid)
  integer :: Mpsid                 ! netcdf id for MATCH PS
  integer :: nm = 1                ! index to prv month in array. init to 1 and toggle between 1 and 2
  integer :: np = 2                ! index to nxt month in array. init to 2 and toggle between 1 and 2
  integer :: mo_nxt = bigint       ! index to nxt month in file

  integer :: date_yr,date_mon,date_day,date_tod   !sxj 20100105
  integer :: ndecade ,timesize                             !sxj 20100105

  real(r8) :: cdaym = inf          ! calendar day of prv month
  real(r8) :: cdayp = inf          ! calendar day of next month
  

!
! This variable is dumb, the dates are in the dataset to be read in but they are
! slightly different than this so getting rid of it causes a change which 
! exceeds roundoff.
!
  real(r8) :: Mid(12)              ! Days into year for mid month date
  data Mid/16.5_r8, 46.0_r8, 75.5_r8, 106.0_r8, 136.5_r8, 167.0_r8, 197.5_r8, 228.5_r8, 259.0_r8, 289.5_r8, 320.0_r8, 350.5_r8 /
  
  public aerosol_initialize  ! read from file, interpolate onto horiz grid
  public aerosol_mass_get         ! interpolate onto pressure levels, and time
  !public aerosol_indirect    ! compute change in effective liqw radius from sulfate
  public aerint              ! read next month info
!
!  values read from file and temporary values used for interpolation
!
!  AEROSOLc is:
!  Cumulative Mass at midpoint of each month
!    on CAM's horizontal grid (col)
!    on MATCH's levels (lev)
!  AEROSOLc
!
  integer, parameter :: paerlev = 28           ! number of levels for aerosol fields (MUST = naerlev)
  integer naerlev                              ! size of level dimension in MATCH data
  integer :: naerlon
  integer :: naerlat
  real(r8), pointer :: M_hybi(:)                ! MATCH hybi
  real(r8), pointer :: M_ps(:,:)                  ! surface pressure from MATCH file
  real(r8), pointer :: AEROSOLc(:,:,:,:,:) ! Aerosol cumulative mass from MATCH
  real(r8), pointer :: M_ps_cam_col(:,:,:) ! PS from MATCH on Cam Columns

  integer :: mxaerl                            ! Maximum level of background aerosol



contains

subroutine aerosol_initialize(phys_state)
!------------------------------------------------------------------
!  Reads in:
!     file from which to read aerosol Masses on CAM grid. Currently
!        assumed to be MATCH ncep runs, averaged by month.
!     NOTE (Data have been externally interpolated onto CAM grid 
!        and backsolved to provide Mid-month values)
!     
!  Populates:
!     module variables:
!       AEROSOLc(pcols,paerlev+1,begchunk:endchunk,naer,2))
!       AEROSOLc(  column_index
!                , level_index (match levels)
!                , chunk_index 
!                , species_index
!                , month = 1:2 )
!       M_hybi(level_index = Lev_MATCH) = pressure at mid-level.
!       M_ps_cam_col(column,chunk,month) ! PS from MATCH on Cam Columns
!
!  Method:
!    read data from file
!    allocate memory for storage of aerosol data on CAM horizontal grid
!    distribute data to remote nodes
!    populates the module variables
!
!------------------------------------------------------------------
   use ioFileMod,    only: getfil
   use filenames,    only: bndtvaer
   !use volcanicmass, only: volcanic_initialize
   !use carbonscales, only: init_scale
   use rgrid,        only: nlon    !, fullgrid
   use physconst,    only: rhoh2o, mwh2o
   !use cam_history,  only: addfld, add_default,phys_decomp
   use history,     only: addfld, add_default, phys_decomp


   use mpishorthand



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
# 183 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/physics/cam1/prescribed_aerosols.F90" 2

   type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

! local variables

   integer, allocatable :: nlon_aer(:)            ! number of lons per lat on bdy dataset
   integer :: naerlev

   integer dateid                       ! netcdf id for date variable
   integer secid                        ! netcdf id for seconds variable
   integer londimid                     ! netcdf id for longitude dimension
   integer latdimid                     ! netcdf id for latitude dimension
   integer levdimid                     ! netcdf id for level dimension
   integer :: nlonid = -1               ! netcdf id for number of lons per lat

   integer timesiz                      ! number of time samples (=12) in netcdf file
   integer latid                        ! netcdf id for latitude variable
   integer Mhybiid                      ! netcdf id for MATCH hybi
   integer timeid                       ! netcdf id for time variable
   !integer dimids(nf90_max_var_dims)    ! variable shape
   integer dimids(nf_max_var_dims)    ! variable shape
   integer :: start(4)                  ! start vector for netcdf calls
   integer :: kount(4)                  ! count vector for netcdf calls
   integer mo                           ! month index
   integer mm                           ! month index
   integer m                            ! constituent index
   integer :: n                         ! loop index
   integer :: i,j,k                     ! spatial indices
   integer :: date_aer(12)              ! Date on aerosol dataset (YYYYMMDD)
   integer :: date_aer216(216)           ! Date on aerosol dataset (YYYYMMDD) for multi-years
   integer :: date_aer336(336)           ! Date on aerosol dataset (YYYYMMDD) for multi-years
   integer :: attnum                    ! attribute number
   integer :: ierr                      ! netcdf return code
   real(r8) ::  coldata(paerlev)        ! aerosol field read in from dataset
   integer :: ret
   integer mo_prv                       ! index to previous month
   integer latidx,lonidx

   character(len=8) :: aname                   ! temporary aerosol name
   character(len=8) :: tmp_aero_name(naer_all) ! name for input to boundary data

   character(len=256) :: locfn          ! netcdf local filename to open
!
! aerosol_data will be read in from the aerosol boundary dataset, then scattered to chunks
! after filling in the bottom level with zeros
! 
   real(r8), allocatable :: aerosol_data(:,:,:)    ! aerosol field read in from dataset
   real(r8), allocatable :: aerosol_field(:,:,:)   ! (plon,paerlev+1,plat)  aerosol field to be scattered
   real(r8) :: caldayloc                           ! calendar day of current timestep
   real(r8) :: closelat,closelon

   call t_startf ('aerosol_init')
   !call print_memusage ('Start aerosol_initialize')

   ! set new aerosol names because input file has 1 seasalt bin
   do m = 1, naer_all
      tmp_aero_name(m)=aerosol_name(m)
      if (aerosol_name(m)=='MSSLTA_V') tmp_aero_name(m) = 'MSSLT_V'
      if (aerosol_name(m)=='MSSLTC_V') tmp_aero_name(m) = 'MSSLT_V'
   end do

   if(scenario_carbon_scale.eq."RAMPED") then
     if (masterproc) &
        write(6,*)"aerosol_init: carbon ramping overrides carscl"
   else if(scenario_carbon_scale.eq."FIXED") then
     if (masterproc) &
        write(6,*) "aerosol_init: using namelist value for carscl, or default"
   else
     call endrun ('scenario_carbon_scale must be either RAMPED or FIXED')
   endif

   if(scenario_prescribed_sulfur.eq."RAMPED") then
     call endrun ('aerosol_init: using sulfur from external historical sequence is not implemented')
   else if(scenario_prescribed_sulfur.eq."FIXED") then
     if (masterproc) &
        write(6,*) "aerosol_init: Using sulfate from cyclic climatology"
     if ( rampyear_prescribed_sulfur .ne. bigint ) then
       call endrun ('aerosol_init: rampyear_prescribed_sulfur set but there is only one year in the climatology')
     endif
   else
     call endrun ('scenario_prescribed_sulfur must be either RAMPED or FIXED')
   endif

   if (prescribed_sulfur  == 'off') then
      call endrun ('aerosol_init: prescribed sulfur = off is not implemented, try passive')
   else if (.not.(prescribed_sulfur == 'direct' .or. prescribed_sulfur == 'passive')) then
      call endrun ('prescribed sulfur must be one of off, direct, or passive')
   endif

   ! compute mxaerl

   allocate (AEROSOLc(pcols,paerlev+1,begchunk:endchunk,naer,2))
   AEROSOLc(:,:,:,:,:) = 0._r8

 !  call setmxaerl ()
 !  set levels for background aerosol


   call c_coupler_get_current_calendar_time(caldayloc)
   if (caldayloc < Mid(1)) then
      mo_prv = 12
      mo_nxt =  1
   else if (caldayloc >= Mid(12)) then
      mo_prv = 12
      mo_nxt =  1
   else
      do i = 2 , 12
         if (caldayloc < Mid(i)) then
            mo_prv = i-1
            mo_nxt = i
            exit
         end if
      end do
   end if
   ! Set initial calendar day values
   cdaym = Mid(mo_prv)
   cdayp = Mid(mo_nxt)
   !write(6,*) "caldayloc:  ",caldayloc
   !write(6,*) "Mid : ",Mid
   write(6,*) "mo_prv: ",mo_prv
   write(6,*) "mo_nxt: ",mo_nxt
   !write(6,*) "debug prescribed_aerosols.F90",304

   call c_coupler_get_current_time(date_yr,date_mon,date_day,date_tod)
   write(6,*) "date_yr ",date_yr
   write(6,*) "date_mon",date_mon
   write(6,*) "date_day",date_day
   !write(6,*) "date_tod",date_tod
   !write(6,*) "debug prescribed_aerosols.F90",311

   aerosol_datan%isncol=.false.

   call getfil (bndtvaer, locfn, 0)

   write(6,*) "bndtvaer: ",bndtvaer
   !write(6,*) "debug prescribed_aerosols.F90",318
!
!   write(6,*) bndtvaer,locfn
!  write(6,*) "debug prescribed_aerosols.F90_line291"
!  call endrun !! sxj--
!

   !call handle_ncerr( nf90_open (locfn, 0, aernid),&
   !   'prescribed_aerosols.F90', 326)
   call handle_ncerr( nf_open (locfn, 0, aernid),&
      'prescribed_aerosols.F90', 328)
!
   if (single_column) &
      call setlatlonidx(aernid,scmlat,scmlon,closelat,closelon,latidx,lonidx)

   ! Check to see if this dataset is in ncol format. 
   !ierr = nf90_inq_dimid( aernid,  'ncol', londimid )
    ierr = nf_inq_dimid( aernid,  'ncol', londimid )
!
!  write(6,*) aernid,londimid,ierr,NF_NOERR
!  write(6,*) "debug prescribed_aerosols.F90_line310"
!  call endrun !! sxj--
!
   !if ( ierr==NF90_NOERR ) then
   if ( ierr==NF_NOERR ) then
      aerosol_datan%isncol=.true.
      !call handle_ncerr(nf90_close(aernid),'prescribed_aerosols.F90', 344)
      call handle_ncerr(nf_close(aernid),'prescribed_aerosols.F90', 345)
      write(6,*) "debug prescribed_aerosols.F90_line318"
      call endrun !! sxj-- 
      call boundarydata_init(bndtvaer, phys_state, tmp_aero_name, naer, &
                             aerosol_datan, 3)

      aerosolc(:,1:paerlev,:,:,:)=aerosol_datan%fields

      M_ps_cam_col=>aerosol_datan%ps
      M_hybi=>aerosol_datan%hybi

   else 

      ! Allocate memory for dynamic arrays local to this module
      allocate (M_ps_cam_col(pcols,begchunk:endchunk,2))
      allocate (M_hybi(paerlev+1))
      ! TBH:  HACK to avoid use of uninitialized values when ncols < pcols
      M_ps_cam_col(:,:,:) = 0._r8

      if (masterproc) then

         ! find and open file; abort if fail (getfil(,,0)).

         write(6,*)'aerosol_initialize: reading aerosol dataset....  If this seems to be taking too'
         write(6,*)'long, perhaps the dataset is being read from an nfs-mounted filesystem.'

         ! First ensure dataset is CAM-ready

        ! call handle_ncerr(nf90_inquire_attribute (aernid, nf90_global, 'cam-ready', attnum=attnum),&
         !     'interpaerosols needs to be run to create a cam-ready aerosol dataset')
              ! NF_INQ_ATTID (INTEGER NCID, INTEGER VARID,CHARACTER*(*) NAME, INTEGER attnum)

        call handle_ncerr(nf_inq_attid (aernid, nf_global, 'cam-ready', attnum),&  !! sxj-- 
              'interpaerosols needs to be run to create a cam-ready aerosol dataset')


         ! Get and check dimension info

         !call handle_ncerr( nf90_inq_dimid( aernid,  'lon', londimid ),&
         call handle_ncerr( nf_inq_dimid( aernid,  'lon', londimid ),&
              'prescribed_aerosols.F90', 385)
         !call handle_ncerr( nf90_inq_dimid( aernid,  'lev', levdimid ),&
         !call handle_ncerr( nf_inq_dimid( aernid,  'lev', levdimid ),&  !!
		 call handle_ncerr( nf_inq_dimid( aernid,  'lev', levdimid ),&  !! 
              'prescribed_aerosols.F90', 389)
		 !call handle_ncerr( nf90_inq_dimid( aernid, 'time',   timeid ),&
         call handle_ncerr( nf_inq_dimid( aernid, 'time',   timeid ),&
              'prescribed_aerosols.F90', 392)
         !call handle_ncerr( nf90_inq_dimid( aernid,  'lat', latdimid ),&
         call handle_ncerr( nf_inq_dimid( aernid,  'lat', latdimid ),&
              'prescribed_aerosols.F90', 395)

         !NF_INQ_DIMLEN (INTEGER NCID, INTEGER DIMID,INTEGER len)
         !call handle_ncerr( nf90_inquire_dimension( aernid, londimid, len=naerlon ),&
         call handle_ncerr( NF_INQ_DIMLEN( aernid, londimid, naerlon ),&
              'prescribed_aerosols.F90', 400)
		 ! call handle_ncerr( nf90_inquire_dimension( aernid, levdimid, len=naerlev ),&
         call handle_ncerr( NF_INQ_DIMLEN( aernid, levdimid, naerlev ),&
              'prescribed_aerosols.F90', 403)
		 ! call handle_ncerr( nf90_inquire_dimension( aernid, latdimid, len=naerlat ),&
         call handle_ncerr( NF_INQ_DIMLEN( aernid, latdimid, naerlat ),&
              'prescribed_aerosols.F90', 406)
		 ! call handle_ncerr( nf90_inquire_dimension( aernid,   timeid, len=timesiz ),&
         call handle_ncerr( NF_INQ_DIMLEN( aernid,   timeid, timesiz ),&
              'prescribed_aerosols.F90', 409)

         !INTEGER FUNCTION NF_INQ_VARID(INTEGER NCID, CHARACTER*(*) NAME,INTEGER varid)
         !call handle_ncerr( nf90_inq_varid( aernid, 'date',   dateid ),&
         call handle_ncerr( nf_inq_varid( aernid, 'date',   dateid ),&  !! sxj no 'date'
              'prescribed_aerosols.F90', 414)
		 ! call handle_ncerr( nf90_inq_varid( aernid, 'datesec', secid ),&  !!sxj no 'datesec'
         call handle_ncerr( nf_inq_varid( aernid, 'datesec', secid ),&
              'prescribed_aerosols.F90', 417)

	timesize=timesiz  ! timesize belong to module "prescribed_aerosols"

    !write(6,*) 'lon', londimid, naerlon
	!write(6,*) 'lev', levdimid, naerlev
    !write(6,*) 'time_sxj',  timeid, timesiz
	!write(6,*) 'lat', latdimid, naerlat
	!write(6,*) 'date', dateid
	!write(6,*) 'datesec', secid
    !write(6,*) "debug prescribed_aerosols.F90",427
    !call endrun !! sxj-- 




         do m = 1, naer
            aname=aerosol_name(m)
            ! rename because file has only one seasalt field
            if (aname=='MSSLTA_V') aname = 'MSSLT_V'
            if (aname=='MSSLTC_V') aname = 'MSSLT_V'
            !call handle_ncerr( nf90_inq_varid( aernid, TRIM(aname), species_id(m)), &
			call handle_ncerr( nf_inq_varid( aernid, TRIM(aname), species_id(m)), &
               'prescribed_aerosols.F90', 440)
         end do


    !write(6,*) "get species_id(m)  "
    !write(6,*) "debug prescribed_aerosols.F90",445
    !call endrun !! sxj-- 

         ! Check for reduced grid

            !nlon_aer(:)   ---- number of lons per lat on bdy dataset
         allocate (nlon_aer(naerlat))
         if (fullgrid) then
            !ret = nf90_inq_varid (aernid, 'nlon', nlonid)
			ret = nf_inq_varid (aernid, 'nlon', nlonid)
            !if (ret == NF90_NOERR) then
			if (ret == NF_NOERR) then
               !call handle_ncerr( nf90_get_var (aernid, nlonid, nlon_aer),&
			    call handle_ncerr( nf_get_var_int (aernid, nlonid, nlon_aer),&
                    'prescribed_aerosols.F90', 459)
               do j=1,naerlat
                  if (nlon_aer(j) /= naerlon) then
                     call endrun ('AEROSOL_INITIALIZE: model grid does not match dataset grid')
                  end if
               end do
            end if
    !write(6,*) "nlon_aer:",nlon_aer
    !write(6,*) "debug prescribed_aerosols.F90",467
    !call endrun !! sxj-- 
         else
            !call handle_ncerr( nf90_inq_varid (aernid, 'nlon', nlonid),&
			call handle_ncerr( nf_inq_varid (aernid, 'nlon', nlonid),&
                    'prescribed_aerosols.F90', 472)
            !call handle_ncerr( nf90_get_var (aernid, nlonid, nlon_aer),&
			call handle_ncerr( nf_get_var_int (aernid, nlonid, nlon_aer),&
                    'prescribed_aerosols.F90', 475)
            do j=1,naerlat
               if (nlon_aer(j) /= nlon(j)) then
                  call endrun ('AEROSOL_INITIALIZE: model grid does not match dataset grid')
               end if
            end do
         end if
         deallocate(nlon_aer)

         !call handle_ncerr( nf90_inq_varid( aernid, 'lat', latid   ),&
		 call handle_ncerr( nf_inq_varid( aernid, 'lat', latid   ),&
              'prescribed_aerosols.F90', 486)

         ! quick sanity check on one field
         !      NF_INQ_VARNDIMS (INTEGER NCID, INTEGER VARID,INTEGER ndims)
         !call handle_ncerr( nf90_inquire_variable (aernid, species_id(1), dimids=dimids),&
		 call handle_ncerr( NF_INQ_VARDIMID (aernid, species_id(1), dimids),&
              'prescribed_aerosols.F90', 492)

  

         if ( (dimids(4) /= timeid) .or. &
              (dimids(3) /= levdimid) .or. &
              (dimids(2) /= latdimid) .or. &
              (dimids(1) /= londimid) ) then
            write(6,*)'AEROSOL READ: Data must be ordered time, lev, lat, lon'
            write(6,*)'data are       ordered as', dimids(4), dimids(3), dimids(2), dimids(1)
            write(6,*)'data should be ordered as', timeid, levdimid, latdimid, londimid
            call endrun ()
         end if

    !write(6,*) 'lon', londimid, naerlon
	!write(6,*) 'lev', levdimid, naerlev
    !write(6,*) 'time',  timeid, timesiz
	!write(6,*) 'lat', latdimid, naerlat
	!write(6,*)  dimids(1:4)
    !write(6,*) "debug prescribed_aerosols.F90",511
    !call endrun !! sxj-- 

    

         ! use hybi,PS from MATCH

         !call handle_ncerr( nf90_inq_varid( aernid, 'hybi', Mhybiid   ),&
		 call handle_ncerr( nf_inq_varid( aernid, 'hybi', Mhybiid   ),&
              'prescribed_aerosols.F90', 520)
         !call handle_ncerr( nf90_inq_varid( aernid, 'PS', Mpsid   ),&
		 call handle_ncerr( nf_inq_varid( aernid, 'PS', Mpsid   ),&
              'prescribed_aerosols.F90', 523)

         ! check dimension order for MATCH's surface pressure

         !call handle_ncerr( nf90_inquire_variable (aernid, Mpsid, dimids=dimids),&
		 call handle_ncerr( NF_INQ_VARDIMID (aernid, Mpsid, dimids),&
              'prescribed_aerosols.F90', 529)
         if ( (dimids(3) /= timeid) .or. &
              (dimids(2) /= latdimid) .or. &
              (dimids(1) /= londimid) ) then
            write(6,*)'AEROSOL READ: Pressure must be ordered time, lat, lon'
            write(6,*)'data are       ordered as', dimids(3), dimids(2), dimids(1)
            write(6,*)'data should be ordered as', timeid, levdimid, latdimid, londimid
            call endrun ()
         end if

    !write(6,*) 'lon', londimid, naerlon
    !write(6,*) 'lev', levdimid, naerlev
    !write(6,*) 'time',  timeid, timesiz
	!write(6,*) 'lat', latdimid, naerlat
	!write(6,*)  dimids(1:4)
    !write(6,*) "debug prescribed_aerosols.F90",544
     !call endrun !! sxj-- 

    if (timesiz /= 12 .and. timesiz /= 216.and. timesiz /= 336) call endrun ('AEROSOL_timesiz: should be 12  216 336')

         ! read in hybi from MATCH
          !         NF_GET_VAR_DOUBLE (NCID, VARID, dval)
         !call handle_ncerr( nf90_get_var (aernid, Mhybiid, M_hybi),&
		  call handle_ncerr( NF_GET_VAR_DOUBLE(aernid, Mhybiid, M_hybi),&
              'prescribed_aerosols.F90', 553)

         ! Retrieve date and sec variables.
         !call handle_ncerr( nf90_get_var (aernid, dateid, date_aer),&
       if (timesiz==216) then
		  call handle_ncerr( nf_get_var_int (aernid, dateid, date_aer216),&
              'prescribed_aerosols.F90', 559)
		  write(6,*)  'date_aer216'
 	      write(6,*)  date_aer216 
          write(6,*) "debug prescribed_aerosols.F90",562
          !call endrun !! sxj-- 
	   else
	   
	   if (timesiz==336) then
              call handle_ncerr( nf_get_var_int (aernid, dateid, date_aer336),&
              'prescribed_aerosols.F90', 568)
		  write(6,*)  'date_aer336'
 	      write(6,*)  date_aer336 
	      write(6,*) "timesiz==336  prescribed_aerosols.F90",571
              !call endrun !! sxj--  
	   else
	   
               call handle_ncerr( nf_get_var_int (aernid, dateid, date_aer),&
              'prescribed_aerosols.F90', 576)
           endif

        endif
			 
         if (timesiz < 12) then
            write(6,*)'AEROSOL READ: When cycling aerosols, dataset must have 12 consecutive ', &
                 'months of data starting with Jan'
            write(6,*)'Current dataset has only ',timesiz,' months'
            call endrun ()
         end if

      
         if (timesiz==12 ) date_aer336(1:12) =date_aer(1:12)
         if (timesiz==216) date_aer336(1:216)=date_aer216(1:216)
         do mo = 1,timesiz
              mm=mod(mo,12)
              if (mm.eq.0) mm=12
              if (mod(date_aer336(mo),10000)/100 /= mm) then
               write(6,*)'AEROSOL READ: When cycling aerosols, dataset must have 12 consecutive ', &
                    'months of data starting with Jan'
               write(6,*)'Month ',mm,' of dataset says date=',date_aer(mo)
               call endrun ()
              end if
         end do
         !write(*,*) "momm"
         !call endrun !! sxj--  

         
         if (single_column) then
            naerlat=1
            naerlon=1
         endif
         kount(:) = (/naerlon,naerlat,paerlev,1/)
      end if          ! masterproc

      ! broadcast hybi to nodes


      call mpibcast (M_hybi, paerlev+1, mpir8, 0, mpicom)
     !call mpibcast (kount,3 , mpiint, 0, mpicom)   ! 3
	  call mpibcast (kount,4 , mpiint, 0, mpicom)   ! 3
      naerlon = kount(1)
      naerlat = kount(2)


      allocate(aerosol_field(kount(1),kount(3)+1,kount(2)))
      allocate(M_ps(kount(1),kount(2)))
      if (masterproc) allocate(aerosol_data(kount(1),kount(2),kount(3)))

      ! Retrieve Aerosol Masses (kg/m^2 in each layer), transpose to model order (lon,lev,lat),
      ! then scatter  to slaves.
	  write(6,*) "nm: ",nm
	  write(6,*) "np: ",np
      write(6,*) "debug prescribed_aerosols.F90",630
      !call endrun !! sxj-- 

      if (nm /= 1 .or. np /= 2) call endrun ('AEROSOL_INITIALIZE: bad nm or np value')
      do n=nm,np
         if (n == 1) then
            mo = mo_prv
         else
            mo = mo_nxt
         end if
         
         ndecade=0
         if (timesiz==216) then 
             if (date_yr<=1849) ndecade=1
	     if (date_yr>=2006) ndecade=18
	     if (date_yr>=2000.and.date_yr<=2005) ndecade=17
	     if (date_yr>=1850.and.date_yr<=1999) then
                ndecade=(date_yr-1830)/10
	     endif
	      mo=mo+12*(ndecade-1)
	endif
        if (timesiz==336) then 
             if (date_yr<=1849) ndecade=1
	     if (date_yr>=2106) ndecade=28
	     if (date_yr>=2100.and.date_yr<=2105) ndecade=27
	     if (date_yr>=1850.and.date_yr<=2099) then
                ndecade=(date_yr-1830)/10
	     endif
	     mo=mo+12*(ndecade-1)
	endif
	write(6,*) "ndecade: ",ndecade
        write(6,*) "mo: ",mo

         do m=1,naer
            if (masterproc) then
               if (single_column) then
                  start(:) = (/lonidx,latidx,1,mo/)
               else
                  start(:) = (/1,1,1,mo/)
               endif
               kount(:) = (/naerlon,naerlat,paerlev,1/)
               ! NF_GET_VARA_DOUBLE(NCID, VARID, START, COUNT, dvals)
               !call handle_ncerr( nf90_get_var (aernid, species_id(m),aerosol_data, start, kount),&
			   call handle_ncerr( NF_GET_VARA_DOUBLE (aernid, species_id(m), start, kount,aerosol_data),&
                    'prescribed_aerosols.F90', 674)
               do j=1,naerlat
                  do k=1,paerlev
                     aerosol_field(:,k,j) = aerosol_data(:,j,k)
                  end do
                  aerosol_field(:,paerlev+1,j) = 0._r8     !!!   value at bottom     
               end do
               
            end if
            call scatter_field_to_chunk (1, paerlev+1, 1, naerlon, aerosol_field, &
                 AEROSOLc(:,:,:,m,n)) !AEROSOLc(pcols,paerlev+1,begchunk:endchunk,naer,2
         end do !m loop

         ! Retrieve PS from Match

         if (masterproc) then
            if (single_column) then
               start(:) = (/lonidx,latidx,mo,-1/)
            else
               start(:) = (/1,1,mo,-1/)
            endif
            kount(:) = (/naerlon,naerlat,1,-1/)
            !call handle_ncerr( nf90_get_var(aernid, Mpsid, M_ps,start,kount),&
			call handle_ncerr( NF_GET_VARA_DOUBLE (aernid, Mpsid, start,kount,M_ps),&
                 'prescribed_aerosols.F90', 698)
         end if
	      write(6,*) "mo: ",mo
         call scatter_field_to_chunk (1, 1, 1, naerlon, M_ps(:,:), M_ps_cam_col(:,:,n))

     !write(6,*) "debug prescribed_aerosols.F90",703
     !call endrun !! sxj-- 
	     if (masterproc) write(6,*) "SXJ-----mo: ",mo

      end do     ! n=nm,np (=1,2)

      if(masterproc) deallocate(aerosol_data)
      deallocate(aerosol_field)

   end if   ! Check to see if this dataset is in ncol format. 

   ! read in volcanic grid data (structural only, no masses)

   !call volcanic_initialize(phys_state)        !!!!!!!!!sxj---

   if ( scenario_carbon_scale == "RAMPED") then
     !! call init_scale()   !!!! sxj----
	  write(*,*) "prescribed_aerosols.F90  init_scale can not be called",720 !!! sxj---
	  call endrun  !sxj
   endif

   ! set aerosol properties. From Ghan

   aername(idxSUL)='SULFATE '
   aername(idxSSLTA)='SEASALT1'
   aername(idxSSLTC)='SEASALT2'
   aername(idxOCPHO)='OCPHOB'
   aername(idxBCPHO)='BCPHOB'
   aername(idxOCPHI)='OCPHIL'
   aername(idxBCPHI)='BCPHIL'
   aername(idxBG)='BACKGRND'
   aername(idxVOLC)='VOLCANIC'
   aername(idxDUSTfirst)='DUST1'
   aername(idxDUSTfirst+1)='DUST2'
   aername(idxDUSTfirst+2)='DUST3'
   aername(idxDUSTfirst+3)='DUST4'
   do i=1,naer_all
      call addfld(aername(i), 'kg/kg', pver, 'A', 'aerosol number concentration', phys_decomp)
   end do

   ! specify aerosol microphysical properties
   density_aer(idxSUL) = 1.77e3
   do i=idxSSLTfirst,idxSSLTfirst+numSSLT-1
      density_aer(i) = 2.20e3
      hygro_aer(i) = 2. * 1.0 * mwh2o * density_aer(i) / (59. * rhoh2o)
   enddo
   density_aer(idxOCPHO) = 1.8e3
   density_aer(idxBCPHO) = 1.7e3
   density_aer(idxOCPHI) = 1.8e3
   density_aer(idxBCPHI) = 1.7e3
   density_aer(idxBG) = density_aer(idxSUL)
   density_aer(idxVOLC) = 1.7e3
   hygro_aer(idxSUL) = 3. * 0.7 * mwh2o * density_aer(idxSUL) / (132. * rhoh2o)
   hygro_aer(idxOCPHO) = 1.e-10 ! zero is a singularity in activation routine
   hygro_aer(idxBCPHO) = 1.e-10
   hygro_aer(idxOCPHI) = 0.14
   hygro_aer(idxBCPHI) = 1.e-10
   hygro_aer(idxBG) = hygro_aer(idxSUL)
   hygro_aer(idxVOLC) = 2. * 0.7 * mwh2o * density_aer(idxVOLC) / (98. * rhoh2o) ! sulfuric acid
   do i=idxDUSTfirst,idxDUSTfirst+numDUST-1
      hygro_aer(i) = 0.14
      density_aer(i) = 2.60e3
   enddo
   dryrad_aer(idxSUL) = 0.05e-6
   dryrad_aer(idxSSLTA) = 0.21e-6
   dryrad_aer(idxSSLTC) = 1.75e-6
   dryrad_aer(idxOCPHO) = 0.0212e-6
   dryrad_aer(idxBCPHO) = 0.0118e-6
   dryrad_aer(idxOCPHI) = 0.0212e-6
   dryrad_aer(idxBCPHI) = 0.0118e-6
   dryrad_aer(idxBG) = dryrad_aer(idxSUL)
   dryrad_aer(idxVOLC) = 0.375e-6
   ! these values for dust shouldn't matter, because num_to_aer is prescribed from truncated distribution
   dryrad_aer(idxDUSTfirst) = 0.2e-6
   dryrad_aer(idxDUSTfirst+1) = 3.e-6
   dryrad_aer(idxDUSTfirst+2) = 15.e-6
   dryrad_aer(idxDUSTfirst+3) = 30.e-6
   dispersion_aer(idxSUL) = 2.0
   dispersion_aer(idxSSLTA) = 2.
   dispersion_aer(idxSSLTC) = 2.
   dispersion_aer(idxOCPHO) = 2.24
   dispersion_aer(idxBCPHO) = 2.0
   dispersion_aer(idxOCPHI) = 2.24
   dispersion_aer(idxBCPHI) = 2.0
   dispersion_aer(idxBG) = dispersion_aer(idxSUL)
   dispersion_aer(idxVOLC) = 1.25
   ! these values are still needed in the aerosol activation scheme
   ! even though they aren't needed for num_to_mass_aer
   ! geometric dispersion for a uniform distribution between two sizes
   dispersion_aer(idxDUSTfirst) = (1.0/0.1)**(1./sqrt(12.))
   dispersion_aer(idxDUSTfirst+1) = (2.5/1.0)**(1./sqrt(12.))
   dispersion_aer(idxDUSTfirst+2) = (5.0/2.5)**(1./sqrt(12.))
   dispersion_aer(idxDUSTfirst+3) = (10./5.0)**(1./sqrt(12.))
   do i=1,naer_all
      num_to_mass_aer(i) = 3./(density_aer(i)*4.*pie*dryrad_aer(i)**3*exp(4.5*log(dispersion_aer(i))**2))
   enddo
   !  these values are from Zender et al JDR 2003 Table 2 assuming volume mode radius of 2.5 micron, sigma=2
   num_to_mass_aer(idxDUSTfirst) = 3.484e15
   num_to_mass_aer(idxDUSTfirst+1) = 2.138e14
   num_to_mass_aer(idxDUSTfirst+2) = 2.205e13
   num_to_mass_aer(idxDUSTfirst+3) = 3.165e12
   if (masterproc) then
      do i=1,naer_all
         write(6,'(a,i3,1pg12.3)')'mode, num_to_mass_aer=',i,num_to_mass_aer(i)
      enddo
   end if

   !call print_memusage('End aerosol_initialize')
   call t_stopf ('aerosol_init')

   !write(*,*) "aerosol_initialize___endrun"
   !call endrun

   return
end subroutine aerosol_initialize


subroutine aerosol_mass_get(c, pint, AEROSOLt, scale)
!------------------------------------------------------------------
!  Input:
!     time at which aerosol masses are needed (get_curr_calday())
!     chunk index
!     CAM's vertical grid (pint)
!
!  Output:
!     values for Aerosol Mass at time specified by get_curr_calday
!     on vertical grid specified by pint (AEROSOLt) :: aerosol at time t
!
!  Method:
!     first determine which indexs of aerosols are the bounding data sets
!     interpolate both onto vertical grid aerm(),aerp().
!     from those two, interpolate in time.
!------------------------------------------------------------------
   !use volcanicmass, only: get_volcanic_mass  !!! sxj--- do not take volcanic aerosol into account
   use interpolate_data, only: get_timeinterp_factors
!
! aerosol fields interpolated to current time step
!   on pressure levels of this time step.
! these should be made read-only for other modules
! Is allocation done correctly here?
!
   integer, intent(in) :: c                   ! Chunk Id.
   real(r8), intent(in) :: pint(pcols,pverp)  ! interface pres.
   real(r8), intent(in) :: scale(naer_all)    ! scale each aerosol by this amount

   real(r8), intent(out) :: AEROSOLt(pcols, pver, naer_all) ! aerosols
!
! Local workspace
!
   real(r8) caldayloc                     ! calendar day of current timestep
   real(r8) fact1, fact2                  ! time interpolation factors

   integer i, k, j                        ! spatial indices
   integer m                              ! constituent index
   integer lats(pcols),lons(pcols)        ! latitude and longitudes of column
   integer ncol                           ! number of columns
   
   real(r8) speciesmin(naer)              ! minimal value for each species
!
! values before current time step "the minus month"
! aerosolm(pcols,pver) is value of preceeding month's aerosol masses
! aerosolp(pcols,pver) is value of next month's aerosol masses
!  (think minus and plus or values to left and right of point to be interpolated)
!
   real(r8) AEROSOLm(pcols,pver,naer) ! aerosol mass from MATCH in column,level at previous (minus) month
!
! values beyond (or at) current time step "the plus month"
!
   real(r8) AEROSOLp(pcols,pver,naer) ! aerosol mass from MATCH in column,level at next (plus) month 

   logical error_found
!------------------------------------------------------------------

   call c_coupler_get_current_calendar_time(caldayloc)
!
! Determine time interpolation factors.  1st arg says we are cycling 1 year of data
!
   call get_timeinterp_factors (.true., mo_nxt, cdaym, cdayp, caldayloc, &
                    fact1, fact2, 'GET_AEROSOL:')

!! sxj--debug
     !write(6,*) caldayloc,fact1,fact2
     !write(6,*) "debug prescribed_aerosols.F90",885
     !call endrun !! sxj-- 
!
! interpolate (prv and nxt month) bounding datasets onto cam vertical grid.
! compute mass mixing ratios on CAMS's pressure coordinate
!  for both the "minus" and "plus" months
!
   ncol = get_ncols_p(c)
!
!#ifdef                !  for debug sxj--
!   CALL mpibarrier (mpicom)  
!   write(*,*) "prescribed_aerosols.F90_line808"                     
!   write(*,*) iam,c,"M_ps_cam_col"
!   write(*,*) M_ps_cam_col(:,c,nm)
!   write(*,*) iam,c,"pint(5,:)"
!   write(*,*) pint(5,:)
!#endif
!
   call vert_interpolate (M_ps_cam_col(:,c,nm), pint, nm, AEROSOLm, ncol, c)
   call vert_interpolate (M_ps_cam_col(:,c,np), pint, np, AEROSOLp, ncol, c)
  
   !write(*,*) "prescribed_aerosols.F90_line806"
   !write(*,*) nm,np,caldayloc,fact1,fact2
   !call endrun

!#ifdef                !  for debug sxj--
!   CALL mpibarrier (mpicom)  
!   write(*,*) "prescribed_aerosols.F90_line829"                     
!   write(*,*) iam,c,"AEROSOLm(5,1:pver,1)"
!   write(*,*) AEROSOLm(5,1:pver,1)
!   write(*,*) "prescribed_aerosols.F90 line 833"
!   call endrun
!#endif

!
! Time interpolate.
!
   do m=1,naer
      do k=1,pver
         do i=1,ncol
            AEROSOLt(i,k,m) = AEROSOLm(i,k,m)*fact1 + AEROSOLp(i,k,m)*fact2
         end do
      end do
      ! Partition seasalt aerosol mass
      if (m .eq. idxSSLTA) then
	      AEROSOLt(:ncol,:,m) = (1.-wgt_sscm)*AEROSOLt(:ncol,:,m) ! fraction of seasalt mass in accumulation mode
      elseif (m .eq. idxSSLTC) then
	      AEROSOLt(:ncol,:,m) = wgt_sscm*AEROSOLt(:ncol,:,m)      ! fraction of seasalt mass in coarse mode
      endif
   end do


!
! get background aerosol (tuning) field
!
   !call background (c, ncol, AEROSOLt(:, :, idxBG))   !!!--sxj--

! 
! find volcanic aerosol masses    
!  !!! sxj  no volcanic aerosol data
   !!!!call get_volcanic_mass(c, AEROSOLt(:,:,idxVOLC))   !!! --sxj---

!
! exit if mass is negative (we have previously set
!  cumulative mass to be a decreasing function.)
!
   speciesmin(:) = 0._r8 ! speciesmin(m) = 0 is minimum mass for each species
 
   error_found = .false.
   do m=1,naer
      do k=1,pver
         do i=1,ncol
            if (AEROSOLt(i, k, m) < speciesmin(m)) error_found = .true.
         end do
      end do
   end do
   if (error_found) then
      do m=1,naer
         do k=1,pver
            do i=1,ncol
               if (AEROSOLt(i, k, m) < speciesmin(m)) then
                  write(6,*) 'AEROSOL_INTERPOLATE: negative mass mixing ratio, exiting'
                  write(6,*) 'm, column, pver',m, i, k ,AEROSOLt(i, k, m)
                  call endrun ()
               end if
            end do
         end do
      end do
   end if
!
! scale any AEROSOLS as required
!
   call scale_aerosols (AEROSOLt, ncol, c, scale)

    ! write(6,*) "debug prescribed_aerosols.F90",979
    ! call endrun !! sxj-- 

   return
end subroutine aerosol_mass_get

subroutine scale_aerosols(AEROSOLt, ncol, lchnk, scale)
!-----------------------------------------------------------------
! scale each species as determined by scale factors
!-----------------------------------------------------------------
  integer, intent(in) :: ncol, lchnk ! number of columns and chunk index
  real(r8), intent(in) :: scale(naer_all) ! scale each aerosol by this amount
  real(r8), intent(inout) :: AEROSOLt(pcols, pver, naer_all) ! aerosols
  integer m

  do m = 1, naer_all
     AEROSOLt(:ncol, :, m) = scale(m)*AEROSOLt(:ncol, :, m)
  end do

  return
end subroutine scale_aerosols



subroutine vert_interpolate (Match_ps, pint, n, AEROSOL_mass, ncol, c)
!--------------------------------------------------------------------
! Input: match surface pressure, cam interface pressure, 
!        month index, number of columns, chunk index
! 
! Output: Aerosol mass mixing ratio (AEROSOL_mass)
!
! Method:
!         interpolate column mass (cumulative) from match onto
!           cam's vertical grid (pressure coordinate)
!         convert back to mass mixing ratio
!
!--------------------------------------------------------------------

   real(r8), intent(out) :: AEROSOL_mass(pcols,pver,naer)  ! aerosol mass from MATCH
   real(r8), intent(in) :: Match_ps(pcols)                ! surface pressure at a particular month
   real(r8), intent(in) :: pint(pcols,pverp)              ! interface pressure from CAM

   integer, intent(in) :: ncol,c                          ! chunk index and number of columns
   integer, intent(in) :: n                               ! prv or nxt month index
!
! Local workspace
!
   integer m                           ! index to aerosol species
   integer kupper(pcols)               ! last upper bound for interpolation
   integer i, k, kk, kkstart, kount    ! loop vars for interpolation
   integer isv, ksv, msv               ! loop indices to save

   logical bad                         ! indicates a bad point found
   logical lev_interp_comp             ! interpolation completed for a level 
   logical error_found

   real(r8) AEROSOL(pcols,pverp,naer)  ! cumulative mass of aerosol in column beneath upper 
                                       ! interface of level in column at particular month
   real(r8) dpl, dpu                   ! lower and upper intepolation factors
   real(r8) v_coord                    ! vertical coordinate
   real(r8) AER_diff                   ! temp var for difference between aerosol masses

   call t_startf ('vert_interpolate')
!
! Initialize index array 
!
   do i=1,ncol
      kupper(i) = 1
   end do
!
! assign total mass to topmost level
!
   AEROSOL(:,1,:) = AEROSOLc(:,1,c,:,n)
!
! At every pressure level, interpolate onto that pressure level
!
   do k=2,pver
!
! Top level we need to start looking is the top level for the previous k
! for all longitude points
!
      kkstart = paerlev+1
      do i=1,ncol
         kkstart = min0(kkstart,kupper(i))
      end do
      kount = 0
!
! Store level indices for interpolation
!
! for the pressure interpolation should be comparing
! pint(column,lev) with M_hybi(lev)*M_ps_cam_col(month,column,chunk)
!


      lev_interp_comp = .false.
      do kk=kkstart,paerlev
         if(.not.lev_interp_comp) then
           do i=1,ncol
            v_coord = pint(i,k)
            if (M_hybi(kk)*Match_ps(i) .lt. v_coord .and. v_coord .le. M_hybi(kk+1)*Match_ps(i)) then
               kupper(i) = kk
               kount = kount + 1
            end if
            end do
!
! If all indices for this level have been found, do the interpolation and
! go to the next level
!
! Interpolate in pressure.
!
           if (kount.eq.ncol) then
            do m=1,naer
               do i=1,ncol
                  dpu = pint(i,k) - M_hybi(kupper(i))*Match_ps(i)
                  dpl = M_hybi(kupper(i)+1)*Match_ps(i) - pint(i,k)
                  AEROSOL(i,k,m) = &
                     (AEROSOLc(i,kupper(i)  ,c,m,n)*dpl + &
                     AEROSOLc(i,kupper(i)+1,c,m,n)*dpu)/(dpl + dpu)
               enddo !i
            end do
            lev_interp_comp = .true.
            end if
         end if
      end do
!
! If we've fallen through the kk=1,levsiz-1 loop, we cannot interpolate and
! must extrapolate from the bottom or top pressure level for at least some
! of the longitude points.
!

      if(.not.lev_interp_comp) then
         do m=1,naer
            do i=1,ncol
               if (pint(i,k) .lt. M_hybi(1)*Match_ps(i)) then
                  AEROSOL(i,k,m) =  AEROSOLc(i,1,c,m,n)
               else if (pint(i,k) .gt. M_hybi(paerlev+1)*Match_ps(i)) then
                  AEROSOL(i,k,m) = 0.0_r8
               else
                  dpu = pint(i,k) - M_hybi(kupper(i))*Match_ps(i)
                  dpl = M_hybi(kupper(i)+1)*Match_ps(i) - pint(i,k)
                  AEROSOL(i,k,m) = &
                     (AEROSOLc(i,kupper(i)  ,c,m,n)*dpl + &
                     AEROSOLc(i,kupper(i)+1,c,m,n)*dpu)/(dpl + dpu)
               end if
            end do
         end do

         if (kount.gt.ncol) then
            call endrun ('VERT_INTERPOLATE: Bad data: non-monotonicity suspected in dependent variable')
         end if
      end if
   end do

!   call t_startf ('vi_checks')
!
! aerosol mass beneath lowest interface (pverp) must be 0
!
   AEROSOL(1:ncol,pverp,:) = 0._r8
!
! Set mass in layer to zero whenever it is less than 
!   1.e-40 kg/m^2 in the layer
!
   do m = 1, naer
      do k = 1, pver
         do i = 1, ncol
            if (AEROSOL(i,k,m) < 1.e-40_r8) AEROSOL(i,k,m) = 0._r8
         end do
      end do
   end do
!
! Set mass in layer to zero whenever it is less than 
!   10^-15 relative to column total mass
!
   error_found = .false.
   do m = 1, naer
      do k = 1, pver
         do i = 1, ncol
            AER_diff = AEROSOL(i,k,m)  - AEROSOL(i,k+1,m)       !! sxj -- aaa
            if( abs(AER_diff) < 1e-15_r8*AEROSOL(i,1,m)) then
               AER_diff = 0._r8
            end if
            AEROSOL_mass(i,k,m)= AER_diff 
            if (AEROSOL_mass(i,k,m) < 0) error_found = .true.
         end do
      end do
   end do
   if (error_found) then
      do m = 1, naer
         do k = 1, pver
            do i = 1, ncol
               if (AEROSOL_mass(i,k,m) < 0) then
                  write(6,*)'vert_interpolate: mass < 0, m, col, lev, mass',m, i, k, AEROSOL_mass(i,k,m)
                  write(6,*)'vert_interpolate: aerosol(k),(k+1)',AEROSOL(i,k,m),AEROSOL(i,k+1,m)
                  write(6,*)'vert_interpolate: pint(k+1),(k)',pint(i,k+1),pint(i,k)
                  write(6,*)'n,c',n,c
                  call endrun()
               end if
            end do
         end do
      end do
   end if

!   call t_stopf ('vi_checks')
   call t_stopf ('vert_interpolate')

   return
end subroutine vert_interpolate



subroutine aerint (phys_state)
   implicit none
   type(physics_state), intent(in) :: phys_state(begchunk:endchunk)
   ! out AEROSOLc(:,:,:,m,np),M_ps_cam_col(:,:,np)
   integer :: ntmp                                ! used in index swapping
   integer :: start(4)                            ! start vector for netcdf calls
   integer :: kount(4)                            ! count vector for netcdf calls
   integer :: i,j,k                               ! spatial indices
   integer :: m                                   ! constituent index
   integer :: cols, cole
   integer :: lchnk, ncol
   real(r8) :: caldayloc                          ! calendar day of current timestep
   real(r8) :: aerosol_data(naerlon,naerlat,paerlev)    ! aerosol field read in from dataset
   real(r8) :: aerosol_field(naerlon,paerlev+1,naerlat) ! aerosol field to be scattered
   integer latidx,lonidx
   real(r8) closelat,closelon

   integer :: mo
!
   if (single_column) &
       call setlatlonidx(aernid,scmlat,scmlon,closelat,closelon,latidx,lonidx)
!
! determine if need to read in next month data
! also determine time interpolation factors
!
   call c_coupler_get_current_calendar_time(caldayloc)
   !write(6,*) "caldayloc:",caldayloc
   !write(6,*) "mo_nxt:",mo_nxt
   !write(6,*) "debug prescribed_aerosols.F90",1217
   !call endrun !! sxj--
!
! If model time is past current forward ozone timeslice, then
! masterproc reads in the next timeslice for time interpolation.  Messy logic is 
! for interpolation between December and January (mo_nxt == 1).  Just like
! ozone_data_timestep_init, sstint.
!
   if (caldayloc > cdayp .and. .not. (mo_nxt == 1 .and. caldayloc >= cdaym)) then
      mo_nxt = mod(mo_nxt,12) + 1
      cdaym = cdayp
      cdayp = Mid(mo_nxt)
!
! Check for valid date info
!
      if (.not. (mo_nxt == 1 .or. caldayloc <= cdayp)) then
         call endrun ('AERINT: Non-monotonicity suspected in input aerosol data')
      end if

      ntmp = nm
      nm = np
      np = ntmp

!--
      !write(6,*) "mo_nxt:",mo_nxt
	  !write(6,*) "timesize:",timesize
	  !write(6,*) "date_yr:",date_yr
      !write(6,*) "aerosol_datan%isncol:",aerosol_datan%isncol


      mo=mo_nxt
      if (timesize==216) then 
	  call c_coupler_get_current_time(date_yr,date_mon,date_day,date_tod)
	  ndecade=0
          if (date_yr<=1849) ndecade=1
	  if (date_yr>=2006) ndecade=18
          if (date_yr>=2000.and.date_yr<=2005) ndecade=17
          if (date_yr>=1850.and.date_yr<=1999) then
             ndecade=(date_yr-1830)/10
	  endif
          mo=mo+12*(ndecade-1)
      endif
      if (timesize==336) then 
	  call c_coupler_get_current_time(date_yr,date_mon,date_day,date_tod)
	  ndecade=0
          if (date_yr<=1849) ndecade=1
	  if (date_yr>=2106) ndecade=28
          if (date_yr>=2100.and.date_yr<=2105) ndecade=27
          if (date_yr>=1850.and.date_yr<=2099) then
             ndecade=(date_yr-1830)/10
	  endif
          mo=mo+12*(ndecade-1)
      endif
      !write(6,*) "ndecade: ",ndecade
      !write(6,*) "mo: ",mo
      !write(6,*) "debug prescribed_aerosols.F90",1272

!--
      if(aerosol_datan%isncol) then
         do lchnk=begchunk,endchunk
            ncol=phys_state(lchnk)%ncol
            cols=1
            cole=cols+aerosol_datan%count(cols,lchnk)-1
            do while(cole<=ncol)
               !start=(/aerosol_datan%start(cols,lchnk),mo_nxt,1,-1/)
			   start=(/aerosol_datan%start(cols,lchnk),mo,1,-1/)
               kount=(/aerosol_datan%count(cols,lchnk),1,-1,-1/)
			   !NF_GET_VARA_DOUBLE(NCID, VARID, START, COUNT, dvals)
               !call handle_ncerr( nf90_get_var(aerosol_datan%ncid, aerosol_datan%psid , &
                !    aerosol_datan%ps(cols:cole,lchnk,np), start(1:2), &
                 !   kount(1:2)),&
                 !   'prescribed_aerosols.F90', 1288)
               call handle_ncerr( NF_GET_VARA_DOUBLE(aerosol_datan%ncid, aerosol_datan%psid , &
                     start(1:2), &
                    kount(1:2),aerosol_datan%ps(cols:cole,lchnk,np)), &
                    'prescribed_aerosols.F90', 1292)
               start(2)=1
               !start(3)=mo_nxt
               start(3)=mo
               kount(2)=paerlev
               kount(3)=1
               do m=1,naer
                  !call handle_ncerr( nf90_get_var(aerosol_datan%ncid, aerosol_datan%dataid(m) , &
                      ! aerosol_datan%fields(cols:cole,:,lchnk,m,np),  &
                      ! start(1:3), kount(1:3)),&
                      ! 'prescribed_aerosols.F90', 1302)
				  call handle_ncerr( NF_GET_VARA_DOUBLE(aerosol_datan%ncid, aerosol_datan%dataid(m) , &
                       start(1:3), kount(1:3),aerosol_datan%fields(cols:cole,:,lchnk,m,np)) , &
                       'prescribed_aerosols.F90', 1305)

               end do
               if(cols==ncol) exit
               cols=cols+aerosol_datan%count(cols,lchnk)
               cole=cols+aerosol_datan%count(cols,lchnk)-1
            end do
         end do
         aerosolc(:,1:paerlev,:,:,np)=aerosol_datan%fields(:,:,:,:,np)
      else
         do m=1,naer
            if (masterproc) then
               if (single_column) then
                  naerlon=1
                  naerlat=1
                 ! start(:) = (/lonidx,latidx,1,mo_nxt/)
				 start(:) = (/lonidx,latidx,1,mo/)
               else
                 !start(:) = (/1,1,1,mo_nxt/)
				 start(:) = (/1,1,1,mo/)
               endif
               kount(:) = (/naerlon,naerlat,paerlev,1/)
              ! call handle_ncerr( nf90_get_var (aernid, species_id(m), aerosol_data, start, kount),&
               !     'prescribed_aerosols.F90', 1328)
               call handle_ncerr(NF_GET_VARA_DOUBLE(aernid, species_id(m), start, kount,aerosol_data),&
                    'prescribed_aerosols.F90', 1330)
               do j=1,naerlat
                  do k=1,paerlev
                     aerosol_field(:,k,j) = aerosol_data(:,j,k)
                  end do
                  aerosol_field(:,paerlev+1,j) = 0._r8   ! value at bottom
               end do
            end if
            call scatter_field_to_chunk (1, paerlev+1, 1, naerlon, aerosol_field, &
                 AEROSOLc(:,:,:,m,np))
         end do
!
! Retrieve PS from Match
!
         if (masterproc) then
               if (single_column) then
                  naerlon=1
                  naerlat=1
                  !start(:) = (/lonidx,latidx,mo_nxt,-1/)
				  start(:) = (/lonidx,latidx,mo,-1/)
               else
                  !start(:) = (/1,1,mo_nxt,-1/)
				  start(:) = (/1,1,mo,-1/)
               endif
               kount(:) = (/naerlon,naerlat,1,-1/)
               !call handle_ncerr( nf90_get_var (aernid, Mpsid, M_ps, start, kount),&
               !     'prescribed_aerosols.F90', 1356)
			   call handle_ncerr(NF_GET_VARA_DOUBLE(aernid, Mpsid, start, kount,M_ps),&
                    'prescribed_aerosols.F90', 1358)
               write(6,*)'AERINT: Read aerosols data for julian day', Mid(mo_nxt)
         end if

         call scatter_field_to_chunk (1, 1, 1, naerlon, M_ps(:,:), M_ps_cam_col(:,:,np))
       
	   
	  end if

      if (masterproc) write(6,*) "SXJ-----mo: ",mo

      end if


      return
end subroutine aerint




subroutine setmxaerl()
!------------------------------------
!
! set levels for background aerosol
!
!-----------------------------------
    !use hycoef, only : hypm, hypi
    use pmgrid   !!for comhyb.h 
	implicit none

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
 
# 1388 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/physics/cam1/prescribed_aerosols.F90" 2
   integer k ! index through vertical levels
!
! mxaerl = max number of levels (from bottom) for background aerosol
! Limit background aerosol height to regions below 900 mb
!
   mxaerl = 0
   do k=pver,1,-1
      if (hypm(k) >= 9.e4_r8) mxaerl = mxaerl + 1
   end do
   mxaerl = max(mxaerl,1)
   if (masterproc) then
      write(6,*)'AEROSOLS:  Background aerosol will be limited to ', &
                'bottom ',mxaerl,' model interfaces. Top interface is ', &
                hypi(pverp-mxaerl),' pascals'
   end if

   return
end subroutine setmxaerl


end module prescribed_aerosols



