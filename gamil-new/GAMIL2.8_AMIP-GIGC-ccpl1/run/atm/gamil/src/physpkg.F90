#include <misc.h>
#include <params.h>

!-----------------------------------------------------------------------
!
! Purpose:
! Loop over time, calling driving routines for physics
!
! Method:
! sequence for running the CSM model
!
! Author:
! Original version:  CCM3
!-----------------------------------------------------------------------

subroutine physpkg(phys_state, phys_state0, gw,     ztodt,  &
                   phys_tend,  cldo,  cldn, tcwato, tcwatn, &
                   qcwato,     qcwatn,      lcwato, lcwatn)

    use shr_kind_mod,   only: r8 => shr_kind_r8
    use pmgrid,         only: plon, plat, masterproc
    use ppgrid,         only: pcols, pver
    use buffer,         only: pblht, tpert, qpert, qrs, qrl
    use comsrf
    use mpishorthand
    use phys_grid,      only: get_ncols_p, get_lat_all_p, get_lon_all_p
    use physics_types,  only: physics_state, physics_tend
    use comsrf
    use diagnostics,    only: diag_surf
    use c_coupler_interface_mod
    use coupling_chemistry_model_mod, only:out_fld_for_coupling_chem,out_caculated_flds_for_coupling_chem
    use phys_buffer,    only: pbuf   ! added by SHI Xiangjun

    implicit none

#include <comctl.h>
#include <comsol.h>

    real(r8), intent(in) :: gw(plat) ! Gaussian weights
    real(r8), intent(in) :: ztodt    ! physics time step unless nstep=0

    type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
    type(physics_state), intent(inout) :: phys_state0(begchunk:endchunk)  ! WAN Hui
    type(physics_tend ), intent(out  ) :: phys_tend(begchunk:endchunk)

    real(r8), intent(inout) :: cldo(pcols, pver, begchunk:endchunk)   ! old cloud
    real(r8), intent(inout) :: cldn(pcols, pver, begchunk:endchunk)   ! new cloud
    real(r8), intent(inout) :: tcwato(pcols, pver, begchunk:endchunk) ! old temperature
    real(r8), intent(inout) :: tcwatn(pcols, pver, begchunk:endchunk) ! new temperature
    real(r8), intent(inout) :: qcwato(pcols, pver, begchunk:endchunk) ! old moisture
    real(r8), intent(inout) :: qcwatn(pcols, pver, begchunk:endchunk) ! new moisture
    real(r8), intent(inout) :: lcwato(pcols, pver, begchunk:endchunk) ! cloud liquid water
    real(r8), intent(inout) :: lcwatn(pcols, pver, begchunk:endchunk) ! cloud liquid water

    integer i, m, lat, c, lchnk                ! indices
    integer lats(pcols)                        ! array of latitude indices
    integer lons(pcols)                        ! array of longitude indices
    integer ncol                               ! number of columns
    integer nstep                              ! current timestep number
    integer ncdate                             ! current date in integer format [yyyymmdd]
    integer ncsec                              ! current time of day [seconds]
    integer yr, mon, day                       ! year, month, and day components of a date

    real(r8) tsgave                            ! TS global average
#ifdef SPMD
    real(r8) tsgridpt_glob(plon,plat)          ! TS global summed at each grid point
#endif
    real(r8), save :: tsgridpt(plon,plat)      ! TS summed at each grid point
    real(r8), save :: tszonal(plat)            ! TS summed along each latitude
    integer,  save :: numts                    ! number of samples for monthly ave

    !*** BAB's FV kludge
    ! DONG Li: check this out
    real(r8) tin(pcols, pver, begchunk:endchunk) ! input T, to compute FV output T

    !! (wanhui 2003.06.11)
    real(r8) dudtm, dvdtm, dtdtm
    integer  iic, kk

    call t_startf('physpkg_st')

    nstep = c_coupler_get_nstep()

    !-----------------------------------------------------------------------
    ! 1. Advance time information
    !-----------------------------------------------------------------------

    call advnce(phys_state)

    call c_coupler_execute_procedure("process_surface_data_before_tphysbc", "kernel")

    call t_stopf('physpkg_st')
    !
    !-----------------------------------------------------------------------
    ! 2. Calculate physical tendencies before flux coupler invokation
    !-----------------------------------------------------------------------
    !
    call t_startf('bc_physics')

!$OMP PARALLEL DO PRIVATE (C) schedule(dynamic, 1)
    do c = begchunk, endchunk

!        call t_startf('tphysbc')
        call tphysbc(ztodt, pblht(1,c), tpert(1,c), srfflx_state2d(c)%ts, &
                     qpert(1,1,c), surface_state2d(c)%precl,              &
                     surface_state2d(c)%precc, surface_state2d(c)%precsl, &
                     surface_state2d(c)%precsc,                           &
                     srfflx_state2d(c)%asdir, srfflx_state2d(c)%asdif,    &
                     srfflx_state2d(c)%aldir, srfflx_state2d(c)%aldif,    &
                     snowhland(1,c),                                      &
                     qrs(1,1,c), qrl(1,1,c), surface_state2d(c)%flwds,    &
                     fsns(1,c),  fsnt(1,c),  flns(1,c),    flnt(1,c),     &
                     srfflx_state2d(c)%lwup,   surface_state2d(c)%srfrad, &
                     surface_state2d(c)%sols,  surface_state2d(c)%soll,   &
                     surface_state2d(c)%solsd, surface_state2d(c)%solld,  &
                     cldo(1,1,c), cldn(1,1,c),                            &
                     tcwato(1,1,c), tcwatn(1,1,c), qcwato(1,1,c),         &
                     qcwatn(1,1,c), lcwato(1,1,c), lcwatn(1,1,c),         &
                     phys_state(c), phys_tend(c),                         &
                     icefrac(1,c), landfrac(1,c), ocnfrac(1,c),           &
                     tin(1,1,c),                                          &
                     srfflx_state2d(c)%cflx(1,1),                         & ! added by WAN Hui, according to P Liu 2003
                     prcsnw(1,c),                                         & ! WAN Hui, according to P Liu 2003
                     phys_state0(c),                                      & ! WAN Hui, according to P Liu 2003
                     srfflx_state2d(c)%sst,                               & ! for FGOALS2.0
                     pbuf)                                                  ! added by SHI Xiangjun

!        call t_stopf ('tphysbc')

    end do

!
! Output field information for coupling chemistry model
!
    call t_stopf ('bc_physics')

    call c_coupler_execute_procedure("process_surface_data_after_tphysbc", "kernel")
    !
    !-----------------------------------------------------------------------
    ! 4. Calculate physical tendencies after calling of flux coupler
    !    Not necessary at terminal timestep.
    !-----------------------------------------------------------------------
    !
    call t_startf ('ac_physics')
!$OMP PARALLEL DO PRIVATE (C, NCOL) schedule(dynamic, 1)
    do c = begchunk, endchunk
        ncol = get_ncols_p(c)
        !
        ! 4.1 Surface diagnostics for history files
        !
        call diag_surf(c, ncol,                                                   &
            srfflx_state2d(c)%shf, srfflx_state2d(c)%lhf, srfflx_state2d(c)%cflx, & ! surface flux
            srfflx_state2d(c)%tref, trefmxav(1,c), trefmnav(1,c),                 & ! temperature at 2m height
            srfflx_state2d(c)%qref,                                               & ! specific humidity at 2m height ! For FGOALS2.0
            srfflx_state2d(c)%rhref, phys_state(c)%ps,                            & ! relative humidity at 2m height ! For FGOALS2.0
            srfflx_state2d(c)%wsx, srfflx_state2d(c)%wsy,                         & ! wind stress
            icefrac(1,c), ocnfrac(1,c), landfrac(1,c),                            & ! fractions
            surface_state2d(c)%tssub, tsnam, srfflx_state2d(c)%ts,                & ! surface/subsurface temperature
            sicthk(1,c), snowhland(1,c), snowhice(1,c))                             ! surface snow

        ! 4.2
        !call t_startf('tphysac')
        call tphysac(ztodt, pblht(1,c), qpert(1,1,c), tpert(1,c), srfflx_state2d(c)%shf,    &
                     srfflx_state2d(c)%wsx, srfflx_state2d(c)%wsy, srfflx_state2d(c)%cflx, sgh(1,c), &
                     srfflx_state2d(c)%lhf, landfrac(1,c), snowhland(1,c), srfflx_state2d(c)%tref,   &
                     surface_state2d(c)%precc, surface_state2d(c)%precl, tin(1,1,c), phys_state(c),  &
                     phys_tend(c), ocnfrac(1,c))
        !call t_stopf('tphysac')
    end do
    call t_stopf('ac_physics')

    !-----------------------------------------------------------------------
    ! Calculate the monthly averaged TS
    ! NOTE: the following is only valid if restart on month boundary
    ! Initialize partial sums of global avg ts to 0.
    ! Sum TS pointwise for this timestep and save for monthly ave TS.
    !
    !-----------------------------------------------------------------------
    !
    call t_startf('global_ts')
    if (c_coupler_is_first_step() .or. c_coupler_is_first_restart_step()) then
        tsgridpt(:,:) = 0.
        numts = 0
    end if

!$OMP PARALLEL DO PRIVATE (c, ncol,i,lats,lons)
    do c = begchunk, endchunk
        ncol = get_ncols_p(c)
        call get_lat_all_p(c, ncol, lats)
        call get_lon_all_p(c, ncol, lons)
        do i = 1, ncol
            tsgridpt(lons(i),lats(i)) = tsgridpt(lons(i),lats(i))+srfflx_state2d(c)%ts(i)*gw(lats(i))
        end do
    end do

    numts = numts+1     ! Increment number of time samples

    if (c_coupler_is_end_current_month()) then
#ifdef SPMD
#ifdef TIMING_BARRIERS
        call t_startf('sync_tszonal')
        call mpibarrier(mpicom)
        call t_stopf('sync_tszonal')
#endif
        call mpisum(tsgridpt, tsgridpt_glob, plon*plat, mpir8, 0, mpicom)
        if (masterproc) tsgridpt(:,:) = tsgridpt_glob(:,:)
#endif
        if (masterproc) then
            call c_coupler_get_current_time(yr, mon, day, ncsec)
            ncdate = yr*10000+mon*100+day
            tsgave = 0.0
            
!$OMP PARALLEL DO PRIVATE (lat,i)
            do lat = 1, plat
                tszonal(lat) = 0.
                do i = 1, plon
                    tszonal(lat) = tszonal(lat)+tsgridpt(i,lat)
                end do
            end do
            
            do lat = 1, plat
                tsgave = tsgave+tszonal(lat)
            end do
            if (numts.gt.0) tsgave = tsgave/(2.*plon*numts)
            write(6, "('Notice: At the end of month ', I8.8, ' ', I5.5, ', ')", advance="no") ncdate, ncsec
            write(6, "('averaged Ts of ', I, ' samples is ', F)") numts, tsgave
        endif
        tszonal(:) = 0.
        numts = 0
    end if

    call t_stopf ('global_ts')

    call out_fld_for_coupling_chem('CLDF',cldn)
    call out_fld_for_coupling_chem('FROCEAN',ocnfrac)
    call out_fld_for_coupling_chem('SLP',psl)
    do c = begchunk, endchunk
        call out_fld_for_coupling_chem('EFLUX',srfflx_state2d(c)%lhf(:),c)
        call out_fld_for_coupling_chem('HFLUX',srfflx_state2d(c)%shf(:),c)
        call out_fld_for_coupling_chem('EVAP',srfflx_state2d(c)%cflx(:,1),c)
        call out_fld_for_coupling_chem('TS',srfflx_state2d(c)%tref(:),c)
        call out_fld_for_coupling_chem('TSKIN',srfflx_state2d(c)%ts(:),c)
        call out_fld_for_coupling_chem('PS',phys_state(c)%ps,c)
        call out_fld_for_coupling_chem('PARDR',surface_state2d(c)%sols,c)
        call out_fld_for_coupling_chem('PARDF',surface_state2d(c)%solsd,c)
    end do
    call out_caculated_flds_for_coupling_chem


    return
end subroutine physpkg
