# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/dp_coupling.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/dp_coupling.F90"
!-------------------------------------------------------------------------------
! dynamics - physics coupling module
!-------------------------------------------------------------------------------
module dp_coupling
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,        only: pcols, pver
   use rgrid,         only: nlon
   use pmgrid
   use phys_grid
   use physics_types, only: physics_state, physics_tend
   use constituents,  only: pcnst, pnats

   implicit none

!===============================================================================
contains
!===============================================================================

!===============================================================================
    subroutine d_p_coupling(ps, t3, u3, v3, q3,q31, t31, q32, t32, &
                          omga, phis, phys_state, fully_coupling)
!------------------------------------------------------------------------------
! Coupler for converting dynamics output variables into physics input variables
! also writes dynamics variables (on physics grid) to history file
!------------------------------------------------------------------------------
    use physconst,     only: cappa
    use mpi_gamil
    use comfm1,        only: itime
!------------------------------Arguments--------------------------------
    real(r8), intent(in) :: ps  (ilbnd:ihbnd, beglatex:endlatex)        ! surface pressure
    real(r8), intent(in) :: t3  (ilbnd:ihbnd, beglatex:endlatex, plev)  ! temperature        !!
    real(r8), intent(in) :: t31 (ilbnd:ihbnd, beglatex:endlatex, plev)  ! temperature
    real(r8), intent(in) :: t32 (ilbnd:ihbnd, beglatex:endlatex, plev)  ! temperature

    real(r8), intent(in) :: u3  (ilbnd:ihbnd, beglatex:endlatex, plev)  ! u-wind component   !!
    real(r8), intent(in) :: v3  (ilbnd:ihbnd, beglatex:endlatex, plev)  ! v-wind component   !!(wh 03.10.28)
    real(r8), intent(in) :: q3  (ilbnd:ihbnd, beglatex:endlatex, plev, pcnst+pnats) ! constituents !!
    real(r8), intent(in) :: q31 (ilbnd:ihbnd, beglatex:endlatex, plev) ! constituents
    real(r8), intent(in) :: q32 (ilbnd:ihbnd, beglatex:endlatex, plev) ! constituents

    real(r8), intent(in) :: omga(ilbnd:ihbnd, beglatex:endlatex, plev)      ! vertical velocity
    real(r8), intent(in) :: phis(ilbnd:ihbnd, beglatex:endlatex)            ! Surface geopotential

    logical, intent(in)  :: fully_coupling

    type(physics_state), intent(out), dimension(begchunk:endchunk) :: phys_state
!
!---------------------------Local workspace-----------------------------
    real(r8), allocatable, dimension(:) :: &
       bbuffer, cbuffer              ! transpose buffers

    integer :: i,k,j,m,lchnk         ! indices
    integer :: ncol                  ! number of columns in current chunk
    integer :: lats(pcols)           ! array of latitude indices
    integer :: lons(pcols)           ! array of longitude indices
    integer :: tsize                 ! amount of data per grid point passed to physics
    integer :: bpter(plon,0:plev)    ! offsets into block buffer for packing data
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data
    integer :: m_size
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! copy data from dynamics data structure to physics data structure
!-----------------------------------------------------------------------

    m_size = pcnst
    if (fully_coupling) m_size = pcnst+pnats

    if (local_dp_map) then

!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, I, K, M, LONS, LATS)
       do lchnk = begchunk,endchunk
          ncol = get_ncols_p(lchnk)
          call get_lon_all_p(lchnk, ncol, lons)
          call get_lat_all_p(lchnk, ncol, lats)
          phys_state(lchnk)%ncol  = ncol
          phys_state(lchnk)%lchnk = lchnk

          do i=1,ncol
             phys_state(lchnk)%ps   (i)     = ps  (lons(i),lats(i))
             phys_state(lchnk)%phis (i)     = phis(lons(i),lats(i))
          end do

          do k=1,plev
             do i=1,ncol
                phys_state(lchnk)%t    (i,k)   = t3  (lons(i),lats(i),k)
                phys_state(lchnk)%u    (i,k)   = u3  (lons(i),lats(i),k)
                phys_state(lchnk)%v    (i,k)   = v3  (lons(i),lats(i),k)
                phys_state(lchnk)%omega(i,k)   = omga(lons(i),lats(i),k)
                if (fully_coupling) then
                   phys_state(lchnk)%t1   (i,k)   = t31 (lons(i),lats(i),k)
                   phys_state(lchnk)%q1(i,k)      = q31 (lons(i),lats(i),k)
                endif
             end do
          end do

          do m=1,m_size
             do k=1,plev
                do i=1,ncol
                   phys_state(lchnk)%q(i,k,m) = q3  (lons(i),lats(i),k,m)
                end do
             end do
          end do

       end do

   else

       tsize = 4 + m_size
       if (fully_coupling) tsize = 6 + m_size

       allocate(bbuffer(tsize*block_buf_nrecs))
       allocate(cbuffer(tsize*chunk_buf_nrecs))

!$OMP PARALLEL DO PRIVATE (I, J, K, M, BPTER)
       do j=beglat,endlat

          call block_to_chunk_send_pters(j,plon,plev+1,tsize,bpter)

          do i=1,nlon(j)

           if (i.ge.ibeg0 .and. i.le.iend0) then
             bbuffer(bpter(i,0))   = ps  (i,j)
             bbuffer(bpter(i,0)+1) = phis(i,j)

             do k=1,plev

                bbuffer(bpter(i,k))   = t3  (i,j,k)
                bbuffer(bpter(i,k)+1) = u3  (i,j,k)
                bbuffer(bpter(i,k)+2) = v3  (i,j,k)
                bbuffer(bpter(i,k)+3) = omga(i,j,k)

                do m=1,m_size
                   bbuffer(bpter(i,k)+3+m) = q3  (i,j,k,m)
                end do
                if (fully_coupling) then
                   bbuffer(bpter(i,k)+3+m_size+1) = q31 (i,j,k)
                   bbuffer(bpter(i,k)+3+m_size+2) = t31 (i,j,k)
                endif
             end do
           endif
          end do

       end do

       call transpose_block_to_chunk(tsize, bbuffer, cbuffer)

!$OMP PARALLEL DO PRIVATE (LCHNK, I, K, M, CPTER, NCOL)
       do lchnk = begchunk,endchunk
          ncol = get_ncols_p(lchnk)
          phys_state(lchnk)%ncol  = ncol
          phys_state(lchnk)%lchnk = lchnk

          call block_to_chunk_recv_pters(lchnk,pcols,plev+1,tsize,cpter)

          do i=1,ncol

             phys_state(lchnk)%ps   (i)     = cbuffer(cpter(i,0))
             phys_state(lchnk)%phis (i)     = cbuffer(cpter(i,0)+1)

             do k=1,plev

                phys_state(lchnk)%t    (i,k)   = cbuffer(cpter(i,k))
                phys_state(lchnk)%u    (i,k)   = cbuffer(cpter(i,k)+1)
                phys_state(lchnk)%v    (i,k)   = cbuffer(cpter(i,k)+2)
                phys_state(lchnk)%omega (i,k)   = cbuffer(cpter(i,k)+3)

                do m=1,m_size
                   phys_state(lchnk)%q (i,k,m) = cbuffer(cpter(i,k)+3+m)
                end do
                if (fully_coupling) then 
                   phys_state(lchnk)%q1(i,k)   = cbuffer(cpter(i,k)+m_size+4)
                   phys_state(lchnk)%t1(i,k)   = cbuffer(cpter(i,k)+m_size+5)
                endif

             end do

          end do

       end do

       deallocate(bbuffer)
       deallocate(cbuffer)
   endif

!-----------------------------------------------------------------------
! Fill auxilliary arrays in physics data structure
!-----------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, I, K, M, LONS, LATS)

    do lchnk = begchunk,endchunk
       ncol = get_ncols_p(lchnk)

! pressure arrays
       call plevs0(ncol, pcols, pver, &
                   phys_state(lchnk)%ps,   phys_state(lchnk)%pint,    &
                   phys_state(lchnk)%pmid, phys_state(lchnk)%pdel)

! log(pressure) arrays and Exner function
       do k=1,pver+1
          do i=1,ncol
             phys_state(lchnk)%lnpint(i,k) = log(phys_state(lchnk)%pint(i,k))
          end do
       end do
       do k=1,pver
          do i=1,ncol
             phys_state(lchnk)%lnpmid(i,k) = log(phys_state(lchnk)%pmid(i,k))
             phys_state(lchnk)%exner (i,k) = (phys_state(lchnk)%pint(i,pver+1) / phys_state(lchnk)%pmid(i,k))**cappa
          end do
       end do

    end do

    return
  end subroutine d_p_coupling

!===============================================================================
  subroutine p_d_coupling(phys_state, phys_tend, t2, fu, fv, qminus, q3, q31, t31,fully_coupling)

    use mpi_gamil

!------------------------------------------------------------------------------
! Coupler for converting physics output variables into dynamics input variables
!------------------------------Arguments--------------------------------
    type(physics_state),intent(in), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend), intent(in), dimension(begchunk:endchunk) :: phys_tend

    real(r8), intent(out) :: t2(ilbnd:ihbnd, beglat:endlat, plev)        ! temp tendency
    real(r8), intent(out) :: t31(ilbnd:ihbnd, beglatex:endlatex, plev)       ! temp (K)
    real(r8), intent(out) :: q31(ilbnd:ihbnd, beglatex:endlatex, plev)       ! moisture
    real(r8), intent(out) :: fu(ilbnd:ihbnd, beglat:endlat, plev)        ! u wind tendency
    real(r8), intent(out) :: fv(ilbnd:ihbnd, beglat:endlat, plev)        ! v wind tendency
    real(r8), intent(out) :: qminus(ilbnd:ihbnd, beglatex:endlatex, plev, pcnst) ! constituents
    real(r8), intent(out) :: q3(ilbnd:ihbnd, beglatex:endlatex, plev, pcnst+pnats) ! non-adv constituents

    logical, intent(in)  :: fully_coupling
!
!---------------------------Local workspace-----------------------------
    real(r8), allocatable, dimension(:) :: &
       bbuffer, cbuffer              ! transpose buffers

    integer :: i,k,j,m,lchnk         ! indices
    integer :: ncol                  ! number of columns in current chunk
    integer :: lats(pcols)           ! array of latitude indices
    integer :: lons(pcols)           ! array of longitude indices
    integer :: tsize                 ! amount of data per grid point passed to physics
    integer :: bpter(plon,0:plev)    ! offsets into block buffer for packing data
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data
    integer :: jdyn, jcam
    integer :: m_size
!-----------------------------------------------------------------------

    m_size = pcnst
    if (fully_coupling) m_size = pcnst+pnats

    if (local_dp_map) then

!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, I, K, M, LONS, LATS)

       do lchnk = begchunk,endchunk
          ncol = get_ncols_p(lchnk)
          call get_lon_all_p(lchnk, ncol, lons)
          call get_lat_all_p(lchnk, ncol, lats)

          do k=1,plev
             do i=1,ncol
                t2(lons(i),lats(i),k)   = phys_tend(lchnk)%dTdt (i,k)
                if (fully_coupling) then
                   t31(lons(i),lats(i),k)   = phys_state(lchnk)%t1 (i,k)
                   q31(lons(i),lats(i),k)   = phys_state(lchnk)%q1 (i,k)
                endif
                fu(lons(i),lats(i),k)   = phys_tend(lchnk)%dudt (i,k)
                fv(lons(i),lats(i),k)   = phys_tend(lchnk)%dvdt (i,k)
             end do
          end do

          do m=1,pcnst
             do k=1,plev
                do i=1,ncol
                   qminus(lons(i),lats(i),k,m) = phys_state(lchnk)%q(i,k,m)
                end do
             end do
          end do

          do m=pcnst+1,m_size
             do k=1,plev
                do i=1,ncol
                   q3(lons(i),lats(i),k,m) = phys_state(lchnk)%q(i,k,m)
                end do
             end do
          end do
       end do

    else

       tsize = 3 + m_size
       if (fully_coupling) tsize = 5 + m_size

       allocate(bbuffer(tsize*block_buf_nrecs))
       allocate(cbuffer(tsize*chunk_buf_nrecs))

!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, I, K, M, CPTER)
       do lchnk = begchunk,endchunk
          ncol = get_ncols_p(lchnk)

          call chunk_to_block_send_pters(lchnk,pcols,pver+1,tsize,cpter)

          do i=1,ncol

             do k=1,plev

                cbuffer(cpter(i,k))   = phys_tend(lchnk)%dTdt (i,k)
                cbuffer(cpter(i,k)+1) = phys_tend(lchnk)%dudt (i,k)
                cbuffer(cpter(i,k)+2) = phys_tend(lchnk)%dvdt (i,k)

                if (fully_coupling) then
                   cbuffer(cpter(i,k)+3+m_size)   = phys_state(lchnk)%t1 (i,k)
                   cbuffer(cpter(i,k)+4+m_size)   = phys_state(lchnk)%q1 (i,k)
                endif

                do m=1,m_size
                   cbuffer(cpter(i,k)+2+m) = phys_state(lchnk)%q(i,k,m)
                end do

             end do

          end do

       end do

       call transpose_chunk_to_block(tsize, cbuffer, bbuffer)

!$OMP PARALLEL DO PRIVATE (I, J, K, M, BPTER)
       do j=beglat,endlat

          call chunk_to_block_recv_pters(j,plon,plev+1,tsize,bpter)

          do i=1,nlon(j)

           if (i.ge.ibeg0 .and. i.le.iend0) then
             do k=1,plev

                t2(i,j,k) = bbuffer(bpter(i,k))
                fu(i,j,k) = bbuffer(bpter(i,k)+1)
                fv(i,j,k) = bbuffer(bpter(i,k)+2)

                if (fully_coupling) then
                   t31(i,j,k) = bbuffer(bpter(i,k)+3+m_size)
                   q31(i,j,k) = bbuffer(bpter(i,k)+4+m_size)
                endif

                do m=1,pcnst
                   qminus(i,j,k,m) = bbuffer(bpter(i,k)+2+m)
                end do

                do m=pcnst+1,m_size
                   q3(i,j,k,m)  = bbuffer(bpter(i,k)+2+m)
                end do
             end do
           endif

          end do

       end do

       deallocate(bbuffer)
       deallocate(cbuffer)

    endif

    return
  end subroutine p_d_coupling
end module dp_coupling
