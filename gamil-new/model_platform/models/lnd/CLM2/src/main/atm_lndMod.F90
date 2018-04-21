#include <misc.h>
#include <preproc.h>

module atm_lndMod

!#define CHECKING
#if (defined COUP_CAM)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Atm - Land interface module
! 
! Method: 
! If running as part of cam, the land surface model must use the same 
! grid as the cam. The land surface model calculates its own net solar 
! radiation and net longwave radiation at the surface. The net longwave 
! radiation at the surface will differ somewhat from that calculated in 
! the atmospheric model because the atm model will use the upward 
! longwave flux (or radiative temperature) from the previous time
! step whereas the land surface model uses the flux for the current
! time step. The net solar radiation should equal that calculated
! in the atmospheric model. If not, there is a problem in how the models
! are coupled.
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id: atm_lndMod.F90,v 1.1.2.8 2002/08/28 15:41:10 erik Exp $
!-----------------------------------------------------------------------

use shr_kind_mod, only: r8 => shr_kind_r8
use pmgrid, only: plon, plond, plat
use tracers, only: pcnst, pnats
use rgrid, only: nlon
use ppgrid, only: pcols, begchunk, endchunk
use phys_grid
use comsrf, only :snowhland, srfflx_state2d, srfflx_parm2d,srfflx_parm,surface_state,landfrac
use history, only :  ctitle, inithist, nhtfrq, mfilt
use filenames, only: caseid
use shr_const_mod, only: SHR_CONST_PI
implicit none

private              ! By default make data private
integer :: landmask(plon,plat) !2d land mask
integer, allocatable, dimension(:,:) :: landmask_chunk

integer , private, parameter :: nsend_atm = 16
real(r8), private :: send2d(plon,nsend_atm, plat) !output to clm
real(r8), allocatable, dimension(:,:,:) :: send2d_chunk

integer , private, parameter :: nrecv_atm = 13
real(r8), private :: recv2d(plon,nrecv_atm, plat) !input from clm
real(r8), allocatable, dimension(:,:,:) :: recv2d_chunk

!Added by Li Ruizhe
real(r8), allocatable, private, dimension(:,:) :: atm_send_buffer
real(r8), allocatable, private, dimension(:,:) :: atm_recv_buffer
real(r8), allocatable, private, dimension(:,:) :: lnd_send_buffer
real(r8), allocatable, private, dimension(:,:) :: lnd_recv_buffer
type buffer_index
    integer :: lchunk
    integer :: column
end type buffer_index
type(buffer_index), allocatable, private, dimension(:) :: atm_send_index
integer, allocatable, private, dimension(:,:)  :: atm_recv_index
integer, allocatable, private, dimension(:)    :: lnd_recv_index
integer, allocatable, private, dimension(:)    :: lnd_send_index
integer, allocatable, private, dimension(:)    :: atm_buffer_pos
integer, allocatable, private, dimension(:)    :: lnd_buffer_pos
integer, allocatable, private, dimension(:)    :: atm_buffer_len
integer, allocatable, private, dimension(:)    :: lnd_buffer_len
integer, private :: atm_buffer_total_len
integer, private :: lnd_buffer_total_len

public atmlnd_ini, atmlnd_drv ! Public interfaces


!===============================================================================
CONTAINS
!===============================================================================

  subroutine atmlnd_ini(srfflx2d)

    use initializeMod, only : initialize           !initialization of clm 
    use lnd_atmMod, only : allocate_atmlnd_ini, lnd_to_atm_mapping_ini 
    use error_messages, only: alloc_err
#if ( defined SPMD )
    use mpishorthand
#endif
    use commap    
    use c_coupler_interface_mod
    use filenames,    only: mss_irt
    !Added by Li Ruizhe
#if ( defined SPMD )
    use spmdMod    , only : masterproc, iam, npes, proc_landi, proc_landf, proc_landpts, &
                          proc_patchi, proc_patchf, proc_patchpts, spmd_init_patch  
    use clm_varmap , only : numland, numpatch, begpatch, endpatch,  &
                               begland, endland, landvec, patchvec, &
                               mapvar_ini   
    use clm_varpar , only : maxpatch
#endif
    !End

#include <comsol.h>
#include <comctl.h>

!-----------------------------------------------------------------------
! Initialize land surface model and obtain relevant atmospheric model 
! arrays back from (i.e. albedos, surface temperature and snow cover over land)
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
    integer  :: i,lat,n,lchnk,ncols !indices
    integer  :: istat               !error return
    integer  :: nstep               !current timestep number
    integer  :: lats(pcols)         !chunk latitudes
    integer  :: lons(pcols)         !chunk longitudes
    real(r8) :: oro_glob(plon,plat)!global oro field
    real(r8) :: lsmlandfrac(plon,plat) !2d fractional land
    real(r8) :: latixy(plon,plat)   !2d latitude  grid (degrees)
    real(r8) :: longxy(plon,plat)   !2d longitude grid (degrees)
    real(r8) :: pi
    type(srfflx_parm), intent(inout), dimension(begchunk:endchunk) :: srfflx2d
    !Added by Li Ruizhe
    integer, dimension(1:plon,1:plat) :: grid_to_chunk, grid_to_column
    integer, dimension(1:plon,1:plat) :: grid_to_land
    integer, allocatable, dimension(:)   :: land_pid
    integer, allocatable, dimension(:)   :: buffer_count
    integer  :: j,k,nlchunks, gchnk

    !End

!-----------------------------------------------------------------------
! Time management variables.

    nstep = c_coupler_get_nstep()

! Allocate land model chunk data structures

   allocate( landmask_chunk(pcols,begchunk:endchunk), stat=istat )
   call alloc_err( istat, 'atmlnd_ini', 'landmask_chunk', &
                   pcols*(endchunk-begchunk+1) )
   allocate( send2d_chunk(pcols,nsend_atm,begchunk:endchunk), stat=istat )
   call alloc_err( istat, 'atmlnd_ini', 'send2d_chunk', &
                   pcols*nsend_atm*(endchunk-begchunk+1) )
   allocate( recv2d_chunk(pcols,nrecv_atm,begchunk:endchunk), stat=istat )
   call alloc_err( istat, 'atmlnd_ini', 'recv2d_chunk', &
                   pcols*nrecv_atm*(endchunk-begchunk+1) )

! Initialize land model

    call gather_chunk_to_field(1,1,1,plon,landfrac,oro_glob)
#if (defined SPMD)
    call mpibcast (oro_glob, size(oro_glob), mpir8, 0, mpicom)
#endif

    pi = SHR_CONST_PI
    longxy(:,:) = 1.e36
    do lat = 1,plat
       do i = 1,nlon(lat)
          longxy(i,lat) = (i-1)*360.0/nlon(lat)
          latixy(i,lat) = (180./pi)*clat(lat)
          if (oro_glob(i,lat) > 0.) then
             landmask(i,lat) = 1
             lsmlandfrac(i,lat) = oro_glob(i,lat)
          else
             landmask(i,lat) = 0
             lsmlandfrac(i,lat) = 0.
          endif
       end do
    end do

    do lchnk=begchunk,endchunk
       ncols = get_ncols_p(lchnk)
       call get_lat_all_p(lchnk,pcols,lats)
       call get_lon_all_p(lchnk,pcols,lons)
       do i=1,ncols
          landmask_chunk(i,lchnk) = landmask(lons(i),lats(i))
       enddo
    enddo

! Initialize albedos, surface temperature, upward longwave radiation,
! and snow depth for land points (used for initial run only)

    call  initialize(eccen    , obliqr   , lambm0  , mvelpp  , caseid  , &
                     ctitle   , nsrest   , nstep   , iradsw  , inithist, &
                     nhtfrq(1), mfilt(1) , longxy  , latixy  , nlon    , &
                     landmask , lsmlandfrac , mss_irt)

! Allocate dynamic memory for atm to/from land exchange

    call allocate_atmlnd_ini()
	
! For initial run only - get 2d data back from land model (Note that 
! in SPMD case, only masterproc contains valid recv2d data) and 
! split 2d data into appropriate arrays contained in module comsrf. 

    if (nstep == 0) then

       call lnd_to_atm_mapping_ini(recv2d)
       call scatter_field_to_chunk(1,nrecv_atm,1,plon,recv2d,recv2d_chunk)

       do lchnk=begchunk,endchunk
          ncols = get_ncols_p(lchnk)
          do i=1,ncols
             if (landmask_chunk(i,lchnk) == 1) then
                srfflx2d(lchnk)%ts(i)    = recv2d_chunk(i, 1,lchnk) 
                srfflx2d(lchnk)%asdir(i) = recv2d_chunk(i, 2,lchnk) 
                srfflx2d(lchnk)%aldir(i) = recv2d_chunk(i, 3,lchnk) 
                srfflx2d(lchnk)%asdif(i) = recv2d_chunk(i, 4,lchnk) 
                srfflx2d(lchnk)%aldif(i) = recv2d_chunk(i, 5,lchnk) 
                snowhland(i,lchnk)     = recv2d_chunk(i, 6,lchnk) 
                srfflx2d(lchnk)%lwup(i)  = recv2d_chunk(i,11,lchnk) 
             endif
          end do
       end do

    endif

#if ( defined SPMD )
    !Added by Li Ruizhe
    ! Build grid to chunks & columns mapping
    grid_to_chunk(:,:) = -1
    grid_to_column(:,:) = -1
    do i = 1, nchunks
       do j = 1, chunk_ncols(i)
          grid_to_chunk(chunk_lon(j,i), chunk_lat(j,i)) = i
          grid_to_column(chunk_lon(j,i), chunk_lat(j,i)) = j
       enddo
    enddo

    ! Build grid to land vector mapping
    grid_to_land(:,:) = -1
    do i = 1,numland
       grid_to_land(landvec%ixy(i), landvec%jxy(i)) = i
    enddo

    ! Build mapping from land to pid
    allocate(land_pid(1:numland))
    do i=0,npes-1
       do j=proc_landi(i),proc_landf(i)
          land_pid(j) = i
       enddo
    enddo

    ! == atm send&recv start ==
    allocate(atm_buffer_pos(0:(npes - 1)))
    allocate(atm_buffer_len(0:(npes - 1)))
    atm_buffer_len(:) = 0
    do lchnk = begchunk, endchunk
       gchnk = lchunk_to_chunk(lchnk)
       ncols = chunk_ncols(gchnk)
       do i = 1, ncols
          if (landmask_chunk(i, lchnk) == 1) then
              k=grid_to_land(chunk_lon(i,gchnk), &
                chunk_lat(i,gchnk))
              j=land_pid(k)
              atm_buffer_len(j) = atm_buffer_len(j) + 1
          endif
       enddo
    enddo

    allocate(buffer_count(0:npes-1))
    do i = 0, npes - 1
        if (i /= 0) then
            atm_buffer_pos(i) = atm_buffer_pos(i-1) + atm_buffer_len(i-1)
        else
            atm_buffer_pos(i) = 0
        endif
        buffer_count(i) = atm_buffer_pos(i)
    enddo
    atm_buffer_total_len = atm_buffer_pos(npes - 1) + atm_buffer_len(npes - 1)

    allocate(atm_send_index(0:(atm_buffer_total_len - 1)))
    allocate(atm_recv_index(1:pcols, begchunk : endchunk))
    do lchnk = begchunk, endchunk
       gchnk = lchunk_to_chunk(lchnk)
       ncols = chunk_ncols(gchnk)
       do i = 1, ncols
          if (landmask_chunk(i, lchnk) == 1) then
              k = grid_to_land(chunk_lon(i,gchnk), &
                chunk_lat(i,gchnk))
              j = land_pid(k)
              atm_send_index(buffer_count(j))%lchunk = lchnk
              atm_send_index(buffer_count(j))%column = i
              atm_recv_index(i, lchnk) = buffer_count(j)
              buffer_count(j) = buffer_count(j) + 1
          endif
       enddo
    enddo
    allocate(atm_send_buffer(1:nsend_atm, 0:(atm_buffer_total_len - 1)))
    allocate(atm_recv_buffer(1:nrecv_atm, 0:(atm_buffer_total_len - 1)))

    !== atm send&recv end ==

    !== lnd send&recv start ==

    nlchunks=sizeof(lchunk_to_chunk)/sizeof(lchunk_to_chunk(1))
    allocate(lnd_buffer_pos(0:(npes - 1)))
    allocate(lnd_buffer_len(0:(npes - 1)))
    lnd_buffer_len(:) = 0
    do lchnk=1,nlchunks
       gchnk = lchunk_to_chunk(lchnk)
       if (gchnk /= -1) then
           ncols = chunk_ncols(gchnk)
           do i = 1, ncols
              j = grid_to_land(chunk_lon(i, gchnk), chunk_lat(i, gchnk))
              if (j /= -1) then 
                  if (land_pid(j) == iam) then
                      lnd_buffer_len(chunk_pid(i, gchnk)) = &
                        lnd_buffer_len(chunk_pid(i, gchnk)) + 1
                  endif
              endif
           enddo
       endif
    enddo

    do i = 0, npes - 1
        if (i /= 0) then
            lnd_buffer_pos(i) = lnd_buffer_pos(i-1) + lnd_buffer_len(i-1)
        else
            lnd_buffer_pos(i) = 0
        endif
        buffer_count(i) = lnd_buffer_pos(i)
    enddo
    lnd_buffer_total_len = lnd_buffer_pos(npes - 1) + lnd_buffer_len(npes - 1)

    allocate(lnd_recv_index(begpatch:endpatch))
    allocate(lnd_send_index(0:(lnd_buffer_total_len - 1)))
    lnd_recv_index(:) = -1
    lnd_send_index(:) = -1
    do lchnk=1,nlchunks
       gchnk = lchunk_to_chunk(lchnk)
       if (gchnk /= -1) then
           ncols = chunk_ncols(gchnk)
           do i = 1, ncols
              j = grid_to_land(chunk_lon(i, gchnk), chunk_lat(i, gchnk))
              if (j /= -1) then
                  if (land_pid(j) == iam) then
                      n=chunk_pid(i, gchnk)
                      do k = 1, maxpatch
                        lnd_recv_index(landvec%patch(j,k)) = buffer_count(n)
                      enddo
                      lnd_send_index(buffer_count(n)) = j
                      buffer_count(n) = &
                        buffer_count(n) + 1
                  endif
              endif
           enddo
       endif
    enddo
    allocate(lnd_send_buffer(1:nrecv_atm, 0:(lnd_buffer_total_len - 1)))
    allocate(lnd_recv_buffer(1:nsend_atm, 0:(lnd_buffer_total_len - 1)), &
        stat=istat)
   call alloc_err( istat, 'atmlnd_ini', 'lnd_recv_buffer', &
                   pcols*(endchunk-begchunk+1) )

    !== lnd send & recv end ==
#endif 
    return
  end subroutine atmlnd_ini

!===============================================================================

  subroutine atmlnd_drv (nstep, iradsw, eccen, obliqr, lambm0, mvelpp,&
                         srf_state,srfflx2d)

!-----------------------------------------------------------------------
! Pack data to be sent to land model into a single array. 
! Send data to land model and call land model driver. 
! Receive data back from land model in a single array.
! Unpack this data into component arrays. 
! NOTE: component arrays are contained in module comsrf.
! When coupling to an atmospheric model: solar radiation depends on 
! surface albedos from the previous time step (based on current
! surface conditions and solar zenith angle for next time step).
! Longwave radiation depends on upward longwave flux from previous
! time step.
!-----------------------------------------------------------------------

#if ( defined SPMD )
    use mpishorthand
#endif
    use lnd_atmMod  !mapping from atm grid space <-> clm tile space
    use comsrf, only:surface_state
    !Added by Li Ruizhe
    use clm_varmap , only : numland, numpatch, begpatch, endpatch,  &
                               begland, endland, landvec, patchvec, &
                               mapvar_ini   
    use clm_varder, only : clm1d
    use spmdMod    , only : masterproc, iam, npes, proc_landi, proc_landf, proc_landpts, &
                          proc_patchi, proc_patchf, proc_patchpts, spmd_init_patch  
    use clm_varcon  , only : rair, cpair, po2, pco2
    use clm_varder  , only : clm
!---------------------------Arguments----------------------------------- 
    integer , intent(in) :: nstep    !Current time index
    integer , intent(in) :: iradsw   !Iteration frequency for shortwave radiation
    real(r8), intent(in) :: eccen    !Earth's orbital eccentricity
    real(r8), intent(in) :: obliqr   !Earth's obliquity in radians
    real(r8), intent(in) :: lambm0   !Mean longitude of perihelion at the vernal equinox (radians)
    real(r8), intent(in) :: mvelpp   !Earth's moving vernal equinox longitude of perihelion + pi (radians)
   type(srfflx_parm), intent(inout), dimension(begchunk:endchunk) :: srfflx2d
   type(surface_state), intent(inout), dimension(begchunk:endchunk) :: srf_state
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
    integer :: i,lat,m,n,lchnk,ncols !indices
    logical doalb          !true if surface albedo calculation time step
    !Added by Li Ruizhe
    integer :: j,k,ierr,stat(mpi_status_size) 
    integer, dimension(0:(npes-1)) :: isend, irecv
    integer :: status
    real(r8) :: forc_snowc, forc_snowl, forc_rainc, forc_rainl, wt
    real(r8), dimension(nrecv_atm) :: tmp
    integer  :: atm_send_type, atm_recv_type
!-----------------------------------------------------------------------

! -----------------------------------------------------------------
! Determine doalb
! [doalb] is a logical variable that is true when the next time
! step is a radiation time step. This allows for the fact that
! an atmospheric model may not do the radiative calculations 
! every time step. For example:
!      nstep dorad doalb
!        1     F     F
!        2     F     T
!        3     T     F
!        4     F     F
!        5     F     T
!        6     T     F
! The following expression for doalb is for example only (it is 
! specific to the NCAR CAM). This variable must be calculated
! appropriately for the host atmospheric model
! -----------------------------------------------------------------

    doalb = iradsw==1 .or. (mod(nstep,iradsw)==0 .and. nstep+1/=1)

! Condense the 2d atmospheric data needed by the land surface model into 
! one array. Note that precc and precl precipitation rates are in units 
! of m/sec. They are turned into fluxes by multiplying by 1000 kg/m^3.

#if (defined CHECKING)
!$OMP PARALLEL DO PRIVATE(lchnk,ncols,i)
    do lchnk=begchunk,endchunk
       ncols = get_ncols_p(lchnk)
       do i=1,ncols
          send2d_chunk(i, 1,lchnk)  =  srf_state(lchnk)%zbot(i)  ! Atmospheric state variable m
          send2d_chunk(i, 2,lchnk)  =  srf_state(lchnk)%ubot(i)  ! Atmospheric state variable m/s
          send2d_chunk(i, 3,lchnk)  =  srf_state(lchnk)%vbot(i)  ! Atmospheric state variable m/s
          send2d_chunk(i, 4,lchnk)  =  srf_state(lchnk)%thbot(i) ! Atmospheric state variable K
          send2d_chunk(i, 5,lchnk)  =  srf_state(lchnk)%qbot(i)  ! Atmospheric state variable kg/kg
          send2d_chunk(i, 6,lchnk)  =  srf_state(lchnk)%pbot(i)  ! Atmospheric state variable Pa
          send2d_chunk(i, 7,lchnk)  =  srf_state(lchnk)%tbot(i)  ! Atmospheric state variable K
          send2d_chunk(i, 8,lchnk)  =  srf_state(lchnk)%flwds(i) ! Atmospheric flux W/m^2
          send2d_chunk(i, 9,lchnk)  =  srf_state(lchnk)%precsc(i)*1000.                  !convert from m/sec to mm/sec
          send2d_chunk(i,10,lchnk)  =  srf_state(lchnk)%precsl(i)*1000.                  !convert from m/sec to mm/sec
          send2d_chunk(i,11,lchnk)  =  (srf_state(lchnk)%precc(i) - srf_state(lchnk)%precsc(i))*1000. !convert from m/sec to mm/sec
          send2d_chunk(i,12,lchnk)  =  (srf_state(lchnk)%precl(i) - srf_state(lchnk)%precsl(i))*1000. !convert from m/sec to mm/sec
          send2d_chunk(i,13,lchnk)  =  srf_state(lchnk)%soll(i)  ! Atmospheric flux W/m^2
          send2d_chunk(i,14,lchnk)  =  srf_state(lchnk)%sols(i)  ! Atmospheric flux W/m^2
          send2d_chunk(i,15,lchnk)  =  srf_state(lchnk)%solld(i) ! Atmospheric flux W/m^2
          send2d_chunk(i,16,lchnk)  =  srf_state(lchnk)%solsd(i) ! Atmospheric flux W/m^2
       end do
    end do

    call t_startf('atmlnd_gather_c2f');
    call gather_chunk_to_field(1,nsend_atm,1,plon,send2d_chunk,send2d)
    call t_stopf('atmlnd_gather_c2f');

! Convert two dimensional atm input data to one dimensional land model data
 
    call t_startf('atm_to_lnd_mapping');
    call atm_to_lnd_mapping(send2d)
    call t_stopf('atm_to_lnd_mapping');
#endif

    !Added by Li Ruizhe
    call t_startf('atm_to_lnd_pc');
!$OMP PARALLEL DO PRIVATE(lchnk,j,i)
    do i = 0, atm_buffer_total_len - 1
        lchnk = atm_send_index(i)%lchunk
        j = atm_send_index(i)%column
        atm_send_buffer(1, i)  =  srf_state(lchnk)%zbot(j)  ! Atmospheric state variable m
        atm_send_buffer(2, i)  =  srf_state(lchnk)%ubot(j)  ! Atmospheric state variable m/s
        atm_send_buffer(3, i)  =  srf_state(lchnk)%vbot(j)  ! Atmospheric state variable m/s
        atm_send_buffer(4, i)  =  srf_state(lchnk)%thbot(j) ! Atmospheric state variable K
        atm_send_buffer(5, i)  =  srf_state(lchnk)%qbot(j)  ! Atmospheric state variable kg/kg
        atm_send_buffer(6, i)  =  srf_state(lchnk)%pbot(j)  ! Atmospheric state variable Pa
        atm_send_buffer(7, i)  =  srf_state(lchnk)%tbot(j)  ! Atmospheric state variable K
        atm_send_buffer(8, i)  =  srf_state(lchnk)%flwds(j) ! Atmospheric flux W/m^2
        atm_send_buffer(9, i)  =  srf_state(lchnk)%precsc(j)*1000.                  !convert from m/sec to mm/sec
        atm_send_buffer(10, i)  =  srf_state(lchnk)%precsl(j)*1000.                  !convert from m/sec to mm/sec
        atm_send_buffer(11, i)  =  (srf_state(lchnk)%precc(j) - srf_state(lchnk)%precsc(j))*1000. !convert from m/sec to mm/sec
        atm_send_buffer(12, i)  =  (srf_state(lchnk)%precl(j) - srf_state(lchnk)%precsl(j))*1000. !convert from m/sec to mm/sec
        atm_send_buffer(13, i)  =  srf_state(lchnk)%soll(j)  ! Atmospheric flux W/m^2
        atm_send_buffer(14, i)  =  srf_state(lchnk)%sols(j)  ! Atmospheric flux W/m^2
        atm_send_buffer(15, i)  =  srf_state(lchnk)%solld(j) ! Atmospheric flux W/m^2
        atm_send_buffer(16, i)  =  srf_state(lchnk)%solsd(j) ! Atmospheric flux W/m^2
    enddo

    call t_startf('atm_to_lnd_comm');
#if (1)
    do i = 0, npes - 1
        if (atm_buffer_len(i) > 0) then
            call mpi_isend(atm_send_buffer(1,atm_buffer_pos(i)), &
                atm_buffer_len(i) * nsend_atm, &
                mpir8, i, iam * npes + i, mpicom, isend(i), ierr)
        endif
    enddo

    do i = 0, npes - 1
        if (lnd_buffer_len(i) > 0) then
            call mpi_irecv(lnd_recv_buffer(1,lnd_buffer_pos(i)), &
                lnd_buffer_len(i) * nsend_atm, &
                mpir8, i, i * npes + iam, mpicom, irecv(i), ierr)
        endif
    enddo

    do i = 0, npes - 1
        if (atm_buffer_len(i) > 0) then
            call MPI_WAIT(isend(i),status,ierr)                         
        endif
    enddo

    do i = 0, npes - 1
        if (lnd_buffer_len(i) > 0) then
            call MPI_WAIT(irecv(i),status,ierr)                         
        endif
    enddo

#else
    call mpi_type_contiguous(nsend_atm, mpir8, atm_send_type, ierr)
    call mpi_type_commit(atm_send_type, ierr)
    call mpi_alltoallv(atm_send_buffer, atm_buffer_len, atm_buffer_pos, &
        atm_send_type, lnd_recv_buffer, lnd_buffer_len, lnd_buffer_pos, &
        atm_send_type, mpicom, ierr)
    call mpi_type_free(atm_send_type, ierr)
#endif
    call t_stopf('atm_to_lnd_comm');

#if (!defined CHECKING)
!$OMP PARALLEL DO PRIVATE(k, forc_rainc, forc_rainl, forc_snowc, forc_snowl)
    do k = begpatch, endpatch
       clm(k)%forc_hgt      = lnd_recv_buffer( 1,lnd_recv_index(k))       !zgcmxy  Atm state m
       clm(k)%forc_u        = lnd_recv_buffer( 2,lnd_recv_index(k))       !forc_uxy  Atm state m/s
       clm(k)%forc_v        = lnd_recv_buffer( 3,lnd_recv_index(k))       !forc_vxy  Atm state m/s
       clm(k)%forc_th       = lnd_recv_buffer( 4,lnd_recv_index(k))       !forc_thxy Atm state K
       clm(k)%forc_q        = lnd_recv_buffer( 5,lnd_recv_index(k))       !forc_qxy  Atm state kg/kg
       clm(k)%forc_pbot     = lnd_recv_buffer( 6,lnd_recv_index(k))       !ptcmxy  Atm state Pa
       clm(k)%forc_t        = lnd_recv_buffer( 7,lnd_recv_index(k))       !forc_txy  Atm state K
       clm(k)%forc_lwrad    = lnd_recv_buffer( 8,lnd_recv_index(k))       !flwdsxy Atm flux  W/m^2
       forc_snowc           = lnd_recv_buffer( 9,lnd_recv_index(k))       !mm/s
       forc_snowl           = lnd_recv_buffer(10,lnd_recv_index(k))       !mm/s
       forc_rainc           = lnd_recv_buffer(11,lnd_recv_index(k))       !mm/s 
       forc_rainl           = lnd_recv_buffer(12,lnd_recv_index(k))       !mm/s 
#if defined(PERGRO)
       ! For error-growth only allow rain not snowfall

       forc_rainc           = forc_rainc + forc_snowc
       forc_rainl           = forc_rainl + forc_snowl
       forc_snowc           = 0.0_r8
       forc_snowl           = 0.0_r8
#endif
       clm(k)%forc_solad(2) = lnd_recv_buffer(13,lnd_recv_index(k))       !forc_sollxy  Atm flux  W/m^2
       clm(k)%forc_solad(1) = lnd_recv_buffer(14,lnd_recv_index(k))       !forc_solsxy  Atm flux  W/m^2 
       clm(k)%forc_solai(2) = lnd_recv_buffer(15,lnd_recv_index(k))       !forc_solldxy Atm flux  W/m^2
       clm(k)%forc_solai(1) = lnd_recv_buffer(16,lnd_recv_index(k))       !forc_solsdxy Atm flux  W/m^2

       ! determine derived quantities

       clm(k)%forc_hgt_u = clm(k)%forc_hgt          !observational height of wind [m] 
       clm(k)%forc_hgt_t = clm(k)%forc_hgt          !observational height of temperature [m]  
       clm(k)%forc_hgt_q = clm(k)%forc_hgt          !observational height of humidity [m]      
       clm(k)%forc_vp    = clm(k)%forc_q*clm(k)%forc_pbot / (0.622+0.378*clm(k)%forc_q)   
       clm(k)%forc_rho   = (clm(k)%forc_pbot-0.378*clm(k)%forc_vp) / (rair*clm(k)%forc_t) 
       clm(k)%forc_co2   = pco2*clm(k)%forc_pbot                                          
       clm(k)%forc_o2    = po2*clm(k)%forc_pbot                                           

       ! Determine precipitation needed by clm

       clm(k)%forc_rain = forc_rainc + forc_rainl
       clm(k)%forc_snow = forc_snowc + forc_snowl

       if ( clm(k)%forc_snow > 0.0_r8  .and. clm(k)%forc_rain > 0.0_r8 ) then
          write(6,*) 'kpatch= ',k,' snow= ',clm(k)%forc_snow,' rain= ',clm(k)%forc_rain, &
               ' CLM cannot currently handle both non-zero rain and snow'
          call endrun
       elseif (clm(k)%forc_rain > 0.) then
          clm(k)%itypprc = 1
       elseif (clm(k)%forc_snow > 0.) then
          clm(k)%itypprc = 2
       else
          clm(k)%itypprc = 0
       endif
    end do
#else
    call atm_to_lnd_checking(lnd_recv_buffer, lnd_recv_index, &
        lnd_buffer_total_len)
#endif
    call t_stopf('atm_to_lnd_pc');
    !Added by Li Ruizhe end

! Call land model driver

!!    write(6,*) 'calling land model driver...'
    call driver (doalb, eccen, obliqr, lambm0, mvelpp)

#if (defined CHECKING)
! Convert one dimensional land model output data to two dimensional atm data 

    call t_startf('lnd_to_mapping');
    call lnd_to_atm_mapping(recv2d) 
    call t_stopf('lnd_to_mapping');
    call t_startf('atmlnd_scatter_f2c');
    call scatter_field_to_chunk(1,nrecv_atm,1,plon,recv2d,recv2d_chunk)
    call t_stopf('atmlnd_scatter_f2c');

! Split 2d recv array into component arrays (in module comsrf)

!$OMP PARALLEL DO PRIVATE(lchnk,ncols,i)
    do lchnk=begchunk,endchunk
       ncols = get_ncols_p(lchnk)
       do i=1,ncols
          if (landmask_chunk(i,lchnk) == 1) then
             srfflx2d(lchnk)%ts(i)     =  recv2d_chunk(i, 1,lchnk) 
             srfflx2d(lchnk)%asdir(i)  =  recv2d_chunk(i, 2,lchnk) 
             srfflx2d(lchnk)%aldir(i)  =  recv2d_chunk(i, 3,lchnk) 
             srfflx2d(lchnk)%asdif(i)  =  recv2d_chunk(i, 4,lchnk) 
             srfflx2d(lchnk)%aldif(i)  =  recv2d_chunk(i, 5,lchnk) 
             snowhland(i,lchnk)        =  recv2d_chunk(i, 6,lchnk) 
             srfflx2d(lchnk)%wsx(i)    =  recv2d_chunk(i, 7,lchnk) 
             srfflx2d(lchnk)%wsy(i)    =  recv2d_chunk(i, 8,lchnk) 
             srfflx2d(lchnk)%lhf(i)    =  recv2d_chunk(i, 9,lchnk) 
             srfflx2d(lchnk)%shf(i)    =  recv2d_chunk(i,10,lchnk) 
             srfflx2d(lchnk)%lwup(i)   =  recv2d_chunk(i,11,lchnk) 
             srfflx2d(lchnk)%cflx(i,1) =  recv2d_chunk(i,12,lchnk) 
             srfflx2d(lchnk)%tref(i)   =  recv2d_chunk(i,13,lchnk) 
          endif
       end do
    end do
#endif

    !Added by Li Ruizhe
    call t_startf('lnd_to_atm_pc');
    lnd_send_buffer (:,:) = 0.0_r8
!$OMP PARALLEL DO PRIVATE(k,n,i,j,wt,tmp,m)
    do k= 0, lnd_buffer_total_len - 1
       n = lnd_send_index(k)
       do i=1,maxpatch
           j=landvec%patch(n,i)
           wt=landvec%wtxy(n,i)
           if (wt /= .0) then
               tmp( 1) = clm(j)%t_rad              !tsxy 
               tmp( 2) = clm(j)%albd(1)            !asdir
               tmp( 3) = clm(j)%albd(2)            !aldir
               tmp( 4) = clm(j)%albi(1)            !asdif
               tmp( 5) = clm(j)%albi(2)            !aldif
               tmp( 6) = clm(j)%h2osno/1000.       !snow (convert mm->m)
               tmp( 7) = clm(j)%taux               !taux 
               tmp( 8) = clm(j)%tauy               !tauy
               tmp( 9) = clm(j)%eflx_lh_tot        !lhflx 
               tmp(10) = clm(j)%eflx_sh_tot        !shflx 
               tmp(11) = clm(j)%eflx_lwrad_out     !lwup
               tmp(12) = clm(j)%qflx_evap_tot      !qflx 
               tmp(13) = clm(j)%t_ref2m            !tref
               do m=1,13
               lnd_send_buffer(m,k) = lnd_send_buffer(m,k) + tmp(m) * wt
               enddo
           endif
       enddo
    end do

    call t_startf('lnd_to_atm_comm');

#if (1)
    do i = 0, npes - 1
        if (lnd_buffer_len(i) > 0) then
            call mpi_isend(lnd_send_buffer(1,lnd_buffer_pos(i)), &
                lnd_buffer_len(i) * nrecv_atm, &
                mpir8, i, iam * npes + i, mpicom, isend(i), ierr)
        endif
    enddo

    do i = 0, npes - 1
        if (atm_buffer_len(i) > 0) then
            call mpi_irecv(atm_recv_buffer(1,atm_buffer_pos(i)), &
                atm_buffer_len(i) * nrecv_atm, &
                mpir8, i, i * npes + iam, mpicom, irecv(i), ierr)
        endif
    enddo

    do i = 0, npes - 1
        if (lnd_buffer_len(i) > 0) then
            call MPI_WAIT(isend(i),status,ierr)                         
        endif
    enddo

    do i = 0, npes - 1
        if (atm_buffer_len(i) > 0) then
            call MPI_WAIT(irecv(i),status,ierr)                         
        endif
    enddo

#else
    call mpi_type_contiguous(nrecv_atm, mpir8, atm_recv_type, ierr)
    call mpi_type_commit(atm_recv_type, ierr)
    call mpi_alltoallv(lnd_send_buffer, lnd_buffer_len, lnd_buffer_pos, &
        atm_recv_type, atm_recv_buffer, atm_buffer_len, atm_buffer_pos, &
        atm_recv_type, mpicom, ierr)
    call mpi_type_free(atm_recv_type, ierr)
#endif
    call t_stopf('lnd_to_atm_comm');

#if (!defined CHECKING)
!$OMP PARALLEL DO PRIVATE(lchnk,ncols,i)
    do lchnk=begchunk,endchunk
       ncols = get_ncols_p(lchnk)
       do i=1,ncols
          if (landmask_chunk(i,lchnk) == 1) then
             srfflx2d(lchnk)%ts(i)     =  atm_recv_buffer(1,atm_recv_index(i, lchnk))
             srfflx2d(lchnk)%asdir(i)  =  atm_recv_buffer(2,atm_recv_index(i, lchnk))
             srfflx2d(lchnk)%aldir(i)  =  atm_recv_buffer(3,atm_recv_index(i, lchnk))
             srfflx2d(lchnk)%asdif(i)  =  atm_recv_buffer(4,atm_recv_index(i, lchnk))
             srfflx2d(lchnk)%aldif(i)  =  atm_recv_buffer(5,atm_recv_index(i, lchnk))
             snowhland(i,lchnk)        =  atm_recv_buffer(6,atm_recv_index(i, lchnk))
             srfflx2d(lchnk)%wsx(i)    =  atm_recv_buffer(7,atm_recv_index(i, lchnk))
             srfflx2d(lchnk)%wsy(i)    =  atm_recv_buffer(8,atm_recv_index(i, lchnk))
             srfflx2d(lchnk)%lhf(i)    =  atm_recv_buffer(9,atm_recv_index(i, lchnk))
             srfflx2d(lchnk)%shf(i)    =  atm_recv_buffer(10,atm_recv_index(i, lchnk))
             srfflx2d(lchnk)%lwup(i)   =  atm_recv_buffer(11,atm_recv_index(i, lchnk))
             srfflx2d(lchnk)%cflx(i,1) =  atm_recv_buffer(12,atm_recv_index(i, lchnk))
             srfflx2d(lchnk)%tref(i)   =  atm_recv_buffer(13,atm_recv_index(i, lchnk))
          endif
       end do
    end do
#else
    call lnd_to_atm_checking();
#endif
    call t_stopf('lnd_to_atm_pc');

    !Added by Li Ruizhe end
    
! Reset all other consitutent surfaces fluxes to zero over land

    do lchnk=begchunk,endchunk
       ncols = get_ncols_p(lchnk)
       do i=1,ncols
          if (landmask_chunk(i,lchnk) == 1) then
             do m = 2,pcnst+pnats
                srfflx2d(lchnk)%cflx(i,m) = 0.
             end do
          endif
       end do
    end do
    
    return
  end subroutine atmlnd_drv

!===============================================================================

  subroutine lnd_to_atm_checking()
    integer :: lchnk, i, j, ncols
    write (*,*) "lnd_to_atm_checking!"
    do lchnk=begchunk,endchunk
        ncols = get_ncols_p(lchnk)
        do i=1,ncols
            if (landmask_chunk(i,lchnk) == 1) then
                do j=1,13
                    if (recv2d_chunk(i, j,lchnk) /= &
                        atm_recv_buffer(j,atm_recv_index(i, lchnk))) then
                        write (6,*) "[ERROR2] lchunk, column, index=", lchnk, i, &
                            j, "value=", recv2d_chunk(i, j,lchnk), &
                            atm_recv_buffer(j,atm_recv_index(i, lchnk))
                    endif
                enddo
            endif
        enddo
    enddo
  end subroutine
!===============================================================================


#endif        


end module atm_lndMod

