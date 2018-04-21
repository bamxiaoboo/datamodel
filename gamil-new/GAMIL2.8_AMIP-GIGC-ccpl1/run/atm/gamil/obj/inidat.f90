# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/inidat.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/inidat.F90"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/inidat.F90" 2

# 1 "./params.h" 1
# 15 "./params.h"
 
 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/inidat.F90" 2

module inidat

! (wanhui 2003.05.16)
!-----------------------------------------------------------------------
!
! Purpose:
!
! Method:
!
! Author:
!
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use chemistry, only: chem_init_mix

   real(r8), allocatable :: ps_tmp(:,:)
   real(r8), allocatable :: u3_tmp(:,:,:)
   real(r8), allocatable :: v3_tmp(:,:,:)
   real(r8), allocatable :: t3_tmp(:,:,:)
   real(r8), allocatable :: q3_tmp(:,:,:,:)
   real(r8), allocatable :: qcwat_tmp(:,:,:)
   real(r8), allocatable :: lcwat_tmp(:,:,:)
   real(r8), allocatable :: phis_tmp(:,:)
   real(r8), allocatable :: landm_tmp(:,:)
   real(r8), allocatable :: sgh_tmp(:,:)

!! real(r8) tmassf_tmp
!! real(r8) qmass1_tmp
!! real(r8) qmass2_tmp
!! real(r8) zgsint_tmp
!! real(r8) qmassf_tmp

contains




   subroutine read_inidat
!-----------------------------------------------------------------------
!
! Purpose:
! Read initial dataset and spectrally truncate as appropriate.
!
!-----------------------------------------------------------------------
!
! $Id: inidat.F90,v 1.20.4.4 2002/06/15 13:47:47 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

      use pmgrid
      use rgrid, only: nlon
      use comsrf
!!    use commap, only: w
!!    use physconst, only: gravit
      use phys_grid
      use constituents, only: pcnst, pnats, cnst_name, qmin
      use buffer
      use tracers, only: nusr_adv, nusr_nad, ixuadv, ixunad, ixcldw
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
# 69 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/inidat.F90" 2

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
 
# 70 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/inidat.F90" 2
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
# 72 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/inidat.F90" 2

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
 
# 73 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/inidat.F90" 2
!
! Local workspace
!
      integer i,j,k,m,lat       ! grid and constituent indices
!!    integer i,j,k,m,lat,irow  ! grid and constituent indices
!!    integer ihem              ! hemisphere index
!!    real(r8) pdelb(plond,plev) ! pressure diff between interfaces
!!! using "B" part of hybrid grid only
!!      real(r8) hyad (plev)      ! del (A)
!!      real(r8) pssum            ! surface pressure sum
!!      real(r8) dotproda         ! dot product
!!      real(r8) dotprodb         ! dot product
      real(r8) pertval          ! perturbation value
!!      real(r8) zgssum           ! partial sums of phis
!
! Netcdf related variables
!
      integer lonsiz, latsiz, levsiz       ! Dimension sizes
      integer londimid, levdimid, latdimid ! Dimension ID's
      integer uid, vid, tid, qid           ! Variable ID's
      integer tracid(pcnst+pnats)          ! Variable ID's
      integer phisid, sghid, psid          ! Variable ID's
      integer landmid

      integer strt2d(3)         ! start lon, lat, time indices for netcdf 2-d
      integer strt3d(4)         ! start lon, lev, lat, time for netcdf 3-d
      data strt2d/3*1/          ! Only index 2 will ever change
      data strt3d/4*1/          ! Only indices 2,3 will ever change

      integer cnt2d(3)          ! lon, lat, time counts for netcdf 2-d
      integer cnt3d(4)          ! lon, lat, lev, time counts for netcdf 2-d
      data cnt2d/plon,1,1/      ! 2-d arrs: Always grab only a "plon" slice
      data cnt3d/plon,plev,plat,1/ ! 3-d arrs: Always grab a full time slice

      integer ndims2d           ! number of dimensions
      integer dims2d(NF_MAX_VAR_DIMS) ! variable shape
      integer ndims3d           ! number of dimensions
      integer dims3d(NF_MAX_VAR_DIMS) ! variable shape
      integer tmptype
      integer natt, ret, attlen ! netcdf return values
      logical phis_hires        ! true => PHIS came from hi res topo
      real(r8) arr3d(plon,plev,plat)
      character*(NF_MAX_NAME) tmpname
      character*256 text
      character*80 trunits      ! tracer untis

      real(r8) ply, wk5,wk6
      integer  kpp

!-----------------------------------------------------------------------
! Allocate memory for temporary arrays
!-----------------------------------------------------------------------
!
! Note if not masterproc still might need to allocate array for spmd case
! since each processor calls MPI_scatter
!
      allocate ( ps_tmp(plond,plat) )
      allocate ( u3_tmp(plond,plev,plat) )
      allocate ( v3_tmp(plond,plev,plat) )
      allocate ( t3_tmp(plond,plev,plat) )
      allocate ( q3_tmp(plond,plev,pcnst+pnats,plat) )
      allocate ( qcwat_tmp(plond,plev,plat) )
      allocate ( lcwat_tmp(plond,plev,plat) )
      allocate ( phis_tmp(plond,plat) )
      allocate ( landm_tmp(plond,plat) )
      allocate ( sgh_tmp(plond,plat) )
!!    allocate ( dpsl_tmp(plond,plat) )
!!    allocate ( dpsm_tmp(plond,plat) )
!!    allocate ( vort_tmp(plond,plev,plat) )
!!    allocate ( div_tmp(plond,plev,plat) )
!
!-----------------------------------------------------------------------
! Read in input variables
!-----------------------------------------------------------------------

      if (masterproc) then
!
! Get dimension IDs and lengths
!
         call wrap_inq_dimid  (ncid_ini, 'lat', latdimid)
         call wrap_inq_dimlen (ncid_ini, latdimid, latsiz)
         call wrap_inq_dimid  (ncid_ini, 'lev', levdimid)
         call wrap_inq_dimlen (ncid_ini, levdimid, levsiz)
         call wrap_inq_dimid  (ncid_ini, 'lon', londimid)
         call wrap_inq_dimlen (ncid_ini, londimid, lonsiz)
!
! Get variable id's
! Check that all tracer units are in mass mixing ratios
!
         call wrap_inq_varid (ncid_ini, 'U'   , uid)
         call wrap_inq_varid (ncid_ini, 'V'   , vid)
         call wrap_inq_varid (ncid_ini, 'T'   , tid)
         call wrap_inq_varid (ncid_ini, 'Q'   , qid)
         call wrap_inq_varid (ncid_ini, 'PS'  , psid)
         call wrap_inq_varid (ncid_ini, 'PHIS', phisid)
         call wrap_inq_varid (ncid_ini, 'SGH' , sghid)
         call wrap_inq_varid (ncid_ini, 'LANDM', landmid)
         if (readtrace) then
 !           do m=2,pcnst+pnats
               do m=2,2   !sxj--only read cloud water    !  sxj
               write(6,*) 'read trace', cnst_name(m),"CWAT"
               call wrap_inq_varid (NCID_INI,"CWAT", tracid(m))
               call wrap_get_att_text (NCID_INI,tracid(m),'units', &
                  trunits)
			   write(6,*) "CWAT",tracid(m),"m=",m,"--sxj--"
               if (trunits(1:5) .ne. 'KG/KG' .and. &
                  trunits(1:5) .ne. 'kg/kg') then
                  write(6,*)'INIDAT: tracer units for tracer = ', &
                     cnst_name(m),' must be in KG/KG'
                  call endrun
               endif
            end do
         end if
!
! Check dimension ordering for one 2-d and one 3-d field.
! Assume other arrays of like rank will have dimensions ordered the same.
!
         call wrap_inq_var (ncid_ini, uid, tmpname, tmptype, &
            ndims3d, dims3d, natt)
         if (dims3d(1).ne.londimid .or. dims3d(2).ne.levdimid .or. &
            dims3d(3).ne.latdimid .or. ndims3d.gt.4) then
            write(6,*)'INIDAT: Bad number of dims or ordering on 3d fld'
            call endrun
         end if
         call wrap_inq_var (ncid_ini, psid, tmpname, tmptype, &
            ndims2d, dims2d, natt)
         if (dims2d(1).ne.londimid .or. dims2d(2).ne.latdimid .or. &
            ndims2d.gt.3) then
            write(6,*)'INIDAT: Bad number of dims or ordering on 2d fld'
            call endrun
         end if
!
! Check for presence of 'from_hires' attribute to decide whether to filter
!
         ret = nf_inq_attlen (ncid_ini, phisid, 'from_hires', attlen)
         if (ret.eq.NF_NOERR .and. attlen.gt.256) then
            write(6,*)'INIDAT: from_hires attribute length is too long'
            call endrun
         end if
         ret = nf_get_att_text (ncid_ini, phisid, 'from_hires', text)
         if (ret.eq.NF_NOERR .and. text(1:4).eq.'true') then
            phis_hires = .true.
            write(6,*)'INIDAT: Will filter input PHIS: attribute ',  &
               'from_hires is true'
         else
            phis_hires = .false.
            write(6,*)'INIDAT: Will not filter input PHIS: attribute ', &
               'from_hires is either false or not present'
         end if
!
! Read in 2d fields.
! For stand alone run: get surface temp and 4 (sub)surface temp fields
!
         do j=1,plat
            strt2d(2) = j
            if (ideal_phys .or. aqua_planet) then
             do i=1,nlon(j)
                phis_tmp(i,j) = 0.
                sgh_tmp (i,j) = 0.
             end do
            else
               call wrap_get_vara_realx (ncid_ini, phisid, strt2d, cnt2d, phis_tmp(1,j))
               call wrap_get_vara_realx (ncid_ini, sghid , strt2d, cnt2d, sgh_tmp(1,j))
            end if
            call wrap_get_vara_realx(ncid_ini, landmid, strt2d, cnt2d, landm_tmp(1,j))
            call wrap_get_vara_realx(ncid_ini, psid, strt2d, cnt2d, ps_tmp(1,j))
         end do

!
! Read in 3d fields.
! Copies are done instead of reading directly into
! prognostic arrays to address netcdf slowness on Cray.
! Array syntax would be really nice here.
! Initialize tracers if not read in from input data.
! Initialize all user tracers (advected and non-advectec to 0.)
!
         call wrap_get_vara_realx(ncid_ini, uid, strt3d, cnt3d, arr3d)
         u3_tmp(:plon,:plev,:plat) = arr3d(:plon,:plev,:plat)

         call wrap_get_vara_realx(ncid_ini, vid, strt3d, cnt3d, arr3d)
         v3_tmp(:plon,:plev,:plat) = arr3d(:plon,:plev,:plat)

         call wrap_get_vara_realx(ncid_ini, tid, strt3d, cnt3d, arr3d)
         t3_tmp(:plon,:plev,:plat) = arr3d(:plon,:plev,:plat)

!         if (ideal_phys) then
!            do k=1,plev
!               PLY = (PS0-pmtop)*SIGL(K) + PMTOP
!               WK5 = PLY*0.01/DPALIB
!               KPP = INT(WK5)
!               WK6 = WK5-DFLOAT(KPP)
!               t3_tmp(:plon,k,:plat)=(1.D0-WK6)*TBB(KPP)+WK6*TBB(KPP+1)
!            enddo
!
!             do j=1,plat
!               do k=1,plev
!                 do i=1,plon
!                    t3_tmp(i,k,j)= t3_tmp(i,k,j)+ 1.0-ABS(real(j-plat/2)/plat)
!                 enddo
!               enddo
!             enddo
!         endif

         call wrap_get_vara_realx(ncid_ini, qid, strt3d, cnt3d, arr3d)
         q3_tmp(:plon,:plev,1,:plat) = arr3d(:plon,:plev,:plat)

! sxj---if RK readtrace CWAT ( 3 testtrace not read)
! 2008-11-10 if MG only read CLDLIQ
         if (readtrace) then
            !do m=2,pcnst+pnats
            do m=2,2
			    write(6,*) "sxj--tracid(m)",m
               call wrap_get_vara_realx(ncid_ini, tracid(m), strt3d, cnt3d, arr3d)
               q3_tmp(:plon,:plev,m,:plat) = arr3d(:plon,:plev,:plat)
            end do
			do m=3,pcnst+pnats ! 3~5 CLDICE NUMLIQ NUMICE not read
			   q3_tmp(:plon,:plev,m,:plat) = 0.
			enddo
         else
            do m=2,pcnst+pnats
               q3_tmp(:plon,:plev,m,:plat) = 0.
            end do
         endif
!
! Add random perturbation to temperature if required
!
         if (pertlim.ne.0.0) then
            write(6,*)'INIDAT: Adding random perturbation bounded by +/-', &
                      pertlim,' to initial temperature field'
            do lat=1,plat
               do k=1,plev
                  do i=1,nlon(lat)
                     call random_number (pertval)
                     pertval = 2.*pertlim*(0.5 - pertval)
                     t3_tmp(i,k,lat) = t3_tmp(i,k,lat)*(1. + pertval)
                  end do
               end do
            end do
         endif
!
!
!
! Initialize tracers if not read in from input data.
! Initialize all user tracers (advected and non-advectec to 0.)
! Ensure sufficient constituent concentration at all gridpoints
! The following appears here for consistency with the SLD branch (mvertens).
!
         if (.not. readtrace) then
            do lat=1,plat
               q3_tmp(:plon,:plev,ixcldw,lat) = 0.
               if (nusr_adv .gt. 0) then
                  do m = ixuadv,ixuadv+nusr_adv-1
                     do k=1,plev
                        do i=1,nlon(lat)
                           q3_tmp(i,k,m,lat) = q3_tmp(i,k,1,lat)*10.**(m-ixuadv)
                        end do
                     end do
                  end do
               endif
               if (nusr_nad .gt. 0) then
                  do m = ixunad,ixunad+nusr_nad-1
                     do k=1,plev
                        do i=1,nlon(lat)
                           q3_tmp(i,k,m,lat) = q3_tmp(i,k,1,lat)*10.**(m-ixunad)
                        end do
                     end do
                  end do
               end if
               if (trace_gas) then
                  if (doRamp_ghg ) call ramp_ghg
                  call chem_init_mix(lat, ps_tmp(1,lat), q3_tmp(1,1,1,lat), nlon(lat))
               endif
               if (trace_test1 .or. trace_test2 .or. trace_test3) then
                  call initesttr( q3_tmp(1,1,1,lat),nlon(lat) )
               endif
            end do
         endif

         do lat=1,plat
            call qneg3('INIDAT  ',lat   ,nlon(lat),plond   ,plev    , &
                       pcnst+pnats,qmin ,q3_tmp(1,1,1,lat))
         end do


!
      endif                     ! end of if-masterproc

!!      write(6,*) 'all variables read'
!
!-----------------------------------------------------------------------
! Copy temporary arrays to model arrays
!-----------------------------------------------------------------------
!
!!      write(6,*) 'inidat copyint...'
      call copy_inidat
!
!-----------------------------------------------------------------------
! Deallocate memory for temporary arrays
!-----------------------------------------------------------------------
!
      deallocate ( ps_tmp )
      deallocate ( u3_tmp )
      deallocate ( v3_tmp )
      deallocate ( t3_tmp )
      deallocate ( q3_tmp )
      deallocate ( qcwat_tmp )
      deallocate ( lcwat_tmp )
      deallocate ( phis_tmp )
      deallocate ( landm_tmp )
      deallocate ( sgh_tmp )
!!    deallocate ( dpsl_tmp )
!!    deallocate ( dpsm_tmp )
!!    deallocate ( vort_tmp )
!!    deallocate ( div_tmp )
!
      return
   end subroutine read_inidat

!*********************************************************************C

   subroutine copy_inidat
!-----------------------------------------------------------------------
! Copy temporary arrays to model arrays
! note that the use statements below contain the definitions
! of the model arrays
!-----------------------------------------------------------------------
      use pmgrid
      use prognostics
      use buffer
      use comsrf
      use phys_grid
      use mpi_gamil
      use tracers, only: ixcldw

      use spmd_dyn, only: npes

!-----------------------------------------------------------------------
      implicit none
!------------------------------Commons----------------------------------
!!#include <comqfl.h>
!
! Local workspace
!
      integer, parameter :: iend = i1+plon-1
      integer, parameter :: jend = j1+plat-1
      integer begj,endj
      integer n,i,j,k


      integer :: numperlat         ! number of values per latitude band
      integer :: numsend(0:npes-1) ! number of items to be sent
      integer :: numrecv           ! number of items to be received
      integer :: displs(0:npes-1)  ! displacement array

!-----------------------------------------------------------------------







      begj = beglatex + numbnd
      endj = endlat

!PW Dynamics fields
      call gamil_scatter_2D_array_phys(ps_tmp, ilbnd, ihbnd, beglatex, endlat, ps(ilbnd,beglatex,1)) 
      call gamil_scatter_2D_array_phys(phis_tmp, ilbnd, ihbnd, beglatex, endlatex, phis(ilbnd,beglatex)) 
      call gamil_scatter_3D_array_phys(u3_tmp, ilbnd, ihbnd, beglatex, endlatex, plev, u3(ilbnd,beglatex,1,1))
      call gamil_scatter_3D_array_phys(v3_tmp, ilbnd, ihbnd, beglatex, endlatex, plev, v3(ilbnd,beglatex,1,1))
      call gamil_scatter_3D_array_phys(t3_tmp, ilbnd, ihbnd, beglatex, endlatex, plev, t3(ilbnd,beglatex,1,1))
      call gamil_scatter_3D_array_phys(q3_tmp, ilbnd, ihbnd, beglatex, endlatex, plev*(pcnst+pnats), q3(ilbnd,beglatex,1,1,1))
!!      write(6,*) 'prognostics ok'


!!      write(6,*) 'tmpchunks allocated'
!!      write(6,*) 'landfrac scattered'

      call scatter_field_to_chunk(1,1,1,plond,landm_tmp,landm(1,begchunk))
      call scatter_field_to_chunk(1,1,1,plond,sgh_tmp,sgh(1,begchunk))

!!      write(6,*) 'landfrac, landm, sgh, tsice, ts  : scattered'

!
!JR cloud and cloud water initialization.  Does this belong somewhere else?
!

      if (masterproc) then
         do k = 1, plev
         do j = 1, plat
         do i = 1,plond
            qcwat_tmp(i,k,j) = i * k*j
         enddo
         enddo
         enddo
         write(6,*) 'check scatter gather'
      endif

      call scatter_field_to_chunk(1,plev,1,plond,qcwat_tmp,qcwat(1,1,begchunk,1))
      call gather_chunk_to_field(1,plev,1,plond, qcwat(1,1,begchunk,1), qcwat_tmp)


      if (masterproc) then
         do k = 1, plev
         do j = 1, plat
         do i = 1,plond
            if (qcwat_tmp(i,k,j) .ne. i*k*j) then
               write (6,*) 'error of scatter gather chunk', i, j, k
            endif
         enddo
         enddo
         enddo
         write(6,*) 'check scatter gather'
      endif


      if (masterproc) then
         qcwat_tmp(:plon,:,:) = q3_tmp(:plon,:,1,:)
         lcwat_tmp(:plon,:,:) = q3_tmp(:plon,:,ixcldw,:)
      endif
      call scatter_field_to_chunk(1,plev,1,plond,qcwat_tmp,qcwat(1,1,begchunk,1))
      call scatter_field_to_chunk(1,plev,1,plond,lcwat_tmp,lcwat(1,1,begchunk,1))
      call scatter_field_to_chunk(1,plev,1,plond,t3_tmp,tcwat(1,1,begchunk,1))
      cld   (:,:,:,1) = 0.
      do n=2,3
         cld(:,:,:,n) = 0.
         qcwat(:,:,:,n) = qcwat(:,:,:,1)
         tcwat(:,:,:,n) = tcwat(:,:,:,1)
         lcwat(:,:,:,n) = lcwat(:,:,:,1)
      end do
!!!
!!! Global integerals
!!!
!!      if (masterproc) then
!!         tmassf = tmassf_tmp
!!         qmass1 = qmass1_tmp
!!         qmass2 = qmass2_tmp
!!         qmassf = qmassf_tmp
!!         zgsint = zgsint_tmp
!!      endif

!!#if ( defined  )
!!      call mpibcast (tmass0,1,mpir8,0,mpicom)
!!      call mpibcast (tmassf,1,mpir8,0,mpicom)
!!      call mpibcast (qmass1,1,mpir8,0,mpicom)
!!      call mpibcast (qmass2,1,mpir8,0,mpicom)
!!      call mpibcast (qmassf,1,mpir8,0,mpicom)
!!      call mpibcast (zgsint,1,mpir8,0,mpicom)
!!#endif

   end subroutine copy_inidat

# 624 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/inidat.F90"

end module inidat
