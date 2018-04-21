#include <misc.h>
#include <params.h>

!-----------------------------------------------------------------------
!
! Purpose:
! Loop over time, calling driving routines for physics, dynamics,
! transport
!
! Method:
!
! Author:
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, D. Williamson, August 1992
! Reviewed:          B. Boville, D. Williamson, April 1996
! Restructured:      J. Truesdale, May 1999
!
!-----------------------------------------------------------------------

subroutine stepon

    use shr_kind_mod,   only: r8 => shr_kind_r8
    use history,        only: wshist, wrapup
    use pmgrid
    ! DONG Li: What is "rgrid" for?
    use rgrid
    use prognostics
    use comfm1
    use buffer
    use c_coupler_interface_mod
    use restart_mod,    only: do_restart_write

    use mpishorthand
    use ppgrid,         only: begchunk, endchunk
    use physics_types,  only: physics_state, physics_tend
    use dp_coupling,    only: d_p_coupling, p_d_coupling
    use commap
    use physconst,      only: gravit
    use time_manager,   only: dtdy ! added by WANG Hui
    use moistconvection,only: convection_scheme

    implicit none

#include <comctl.h>
#include <comhyb.h>
#include <comlun.h>

    type(physics_state), allocatable :: phys_state(:)
    type(physics_state), allocatable :: phys_state0(:) ! added by WAN Hui, according to P Liu 2003)
    type(physics_tend ), allocatable :: phys_tend(:)

    real(r8), allocatable :: t2(:,:,:) ! temp tendency
    real(r8), allocatable :: fu(:,:,:) ! u wind tendency
    real(r8), allocatable :: fv(:,:,:) ! v wind tendency

    real(r8) dtime               ! timestep size  (physics package)
    real(r8) ztodt               ! twice time step unless nstep=0
    real(r8) wcstart, wcend      ! wallclock timestamp at start, end of timestep
    real(r8) usrstart, usrend    ! user timestamp at start, end of timestep
    real(r8) sysstart, sysend    ! sys timestamp at start, end of timestep

    real(r8) calday              ! current calendar day
    real(r8) dsghl(plev)
    integer nseq

    integer i, k, lat, j, begj   ! longitude,level,latitude indices
    logical fully_dp_coupling, fully_pd_coupling
    !
    ! Externals
    !
    logical, external :: rstwr  ! whether or not to write restart files
    !
    !-----------------------------------------------------------------------
    call t_startf('stepon_startup'); if(masterproc) write(6,*) '+++++ stepon_startup +++++'
    dtime = c_coupler_get_step_size();         if(masterproc) write(6,*) 'dtime = ', dtime
                                     if(masterproc) write(6,*) 'dtdy  = ', dtdy
    nseq  = dtime/(dtdy-0.01);       if(masterproc) write(6,*) 'nseq  = ', nseq
    !!
    !! fm2003 : calculate dsghl for subroutine 'avnegq' at the end of 'qpdata'
    !!
    dsghl(1) = 0.0
    do k=2,plev
        dsghl(k) = dsig(k-1)/dsig(k)
    enddo

    if (masterproc) write(6, "('Notice: stepon: dsghl set')")

    pmtop = pmtop*0.01d0

    ! WAN Hui 2003.07.08)

    if (c_coupler_is_first_step()) then
        if(masterproc) write(6, "('Notice: stepon: first step start')")
        itime = 0
        !
        ! Calculate vertical motion field
        !
        omga(:,:,:) = 0.0
        if(masterproc) write(6, "('Notice: stepon: set omga to zero')")
        call init_ac_switching(pmtop) ! added by WAN Hui 2003.10.28
    else if (c_coupler_is_first_restart_step()) then
        !sq(:,:,:) = 0.0
        ! DONG Li: clean this out
        if(masterproc) write(6, "('Notice: stepon: sq set to 0.0 in comfm1 (Please clarify this)')")
    else
        if(masterproc) write(6, "('Error: neither first_step nor first_restart_step')")
        call endrun
    end if

    allocate(phys_state(begchunk:endchunk))
    allocate(phys_state0(begchunk:endchunk)) ! added by WAN Hui
    allocate(phys_tend(begchunk:endchunk))
    allocate(t2(ilbnd:ihbnd,beglat:endlat,plev))
    allocate(fu(ilbnd:ihbnd,beglat:endlat,plev))
    allocate(fv(ilbnd:ihbnd,beglat:endlat,plev))
    !
    ! Beginning of basic time step loop
    !
    call t_stopf ('stepon_startup')

    ! Begin time loop.

    do
        call t_startf('stepon_st')
        if (masterproc .and. print_step_cost) then
            call t_stampf(wcstart, usrstart, sysstart)
        end if

        ! DONG Li: clarify this
        !ztodt = 2.0*dtime
        ztodt = dtime

        call c_coupler_get_current_calendar_time(calday)

        if (masterproc) then
            write(6, *)
            write(6, *) 'date:', calday
        end if

        !----------------------------------------------------------
        ! PHYSPKG  Call the Physics package
        !----------------------------------------------------------
        if (masterproc) write(6, *) '------physpkg------'

        begj = beglatex+numbnd

        call t_stopf('stepon_st')
        call t_startf('d_p_coupling')
        fully_dp_coupling = c_coupler_is_first_step() .or. c_coupler_is_first_restart_step()
        call d_p_coupling(ps(ilbnd,beglatex,n3m2), t3(ilbnd,beglatex,1,n3m2), u3(ilbnd,beglatex,1,n3m2), &
                          v3(ilbnd,beglatex,1,n3m2), q3(ilbnd,beglatex,1,1,n3m2), &
                          q31(ilbnd,beglatex,1), t31(ilbnd,beglatex,1), q32(ilbnd,beglatex,1), t32(ilbnd,beglatex,1),&
                          omga, phis, phys_state, fully_dp_coupling)
        call t_stopf('d_p_coupling')

        !!(wh, according to P Liu 2003)
        if ( convection_scheme == 'Tiedtke' ) then
            !for Tiedtke scheme
            call t_startf('d_p_coupling')
            call d_p_coupling(ps(ilbnd,beglatex,n3), t3(ilbnd,beglatex,1,n3), u3(ilbnd,beglatex,1,n3), &
                              v3(ilbnd,beglatex,1,n3), q3(ilbnd,beglatex,1,1,n3), &
                              q31(ilbnd,beglatex,1), t31(ilbnd,beglatex,1), q32(ilbnd,beglatex,1), t32(ilbnd,beglatex,1),&!(ljli)
                              omga, phis, phys_state0, fully_dp_coupling)
            call t_stopf('d_p_coupling')
        end if
        !!(wh, according to P Liu 2003)

        call t_startf('phys_driver')

        if (masterproc) then
            write(6, *) 'ideal_phys =', ideal_phys
            write(6, *) 'adiabatic  =', adiabatic
        end if

        if (ideal_phys) then
            call phys_idealized(phys_state, phys_tend, ztodt, sigl)
        else if (adiabatic) then
            call phys_adiabatic(phys_state, phys_tend)
        else
            call physpkg( &
                phys_state, phys_state0, w, ztodt, phys_tend,     &
                cld(1,1,begchunk,n3m2),   cld(1,1,begchunk,n3),   &
                tcwat(1,1,begchunk,n3m2), tcwat(1,1,begchunk,n3), &
                qcwat(1,1,begchunk,n3m2), qcwat(1,1,begchunk,n3), &
                lcwat(1,1,begchunk,n3m2), lcwat(1,1,begchunk,n3))
        end if
        call t_stopf('phys_driver')

        fully_pd_coupling = .true.
        call t_startf('p_d_coupling')
        call p_d_coupling(phys_state, phys_tend, t2, fu, fv, &
            qminus(ilbnd,beglatex,1,1), q3(ilbnd,beglatex,1,1,n3), q31(ilbnd,beglatex,1), t31(ilbnd,beglatex,1),fully_pd_coupling)
        call t_stopf('p_d_coupling')

        !----------------------------------------------------------
        ! DYNPKG Call the Dynamics Package
        !----------------------------------------------------------
        if (masterproc) write(6, *) '------dynpkg------'

        call t_startf('dynpkg')

        ! accumulate su, sv, st and update q

        call a_c_switching(fu, fv, t2, beglat, endlat)   !!(wh 2003.10.28)

        call dynpkg(dtdy, nseq, dsghl)        !!(wh 2003.10.23)
        ! prepare data for physics

        call c_a_switching(pmtop)            !!(wh 2003.10.28)

        call t_stopf('dynpkg')
        !
        ! Shift time pointers
        !
        call shift_time_indices

        call t_startf('stepon_st')
        if (c_coupler_is_first_restart_step()) then
            call print_memusage
        end if

        ! Set end of run flag.

        !
        !----------------------------------------------------------
        ! History and restart logic: Write and/or dispose history tapes if required
        !----------------------------------------------------------
        !
        call t_startf ('wshist')
        call wshist ()
        call t_stopf ('wshist')

        call c_coupler_execute_procedure("output_field", "kernel")
        !
        ! Write restart file
        !
        call do_restart_write
        !
        ! Dispose necessary files
        !
        call t_startf ('wrapup')
        call wrapup
        call t_stopf ('wrapup')

        if (masterproc .and. print_step_cost) then
            call t_stampf (wcend, usrend, sysend)
            write(6,'(a,3f8.3,a)')'Prv timestep wallclock, usr, sys=', &
                wcend-wcstart, usrend-usrstart, sysend-sysstart, ' seconds'
        end if
        !
        ! Advance timestep before returning to top of loop
        !
        call c_coupler_advance_timer()
        call t_stopf('stepon_st')

        if (c_coupler_check_coupled_run_finished()) then
            if ( masterproc ) write(6,*)'atm: Stopping now'
            nlend = .true.
        end if
        !
        ! Check for end of run
        !
        if (nlend)  then
            deallocate(phys_state)
            deallocate(phys_state0)   !!(wh)
            deallocate(phys_tend)
            deallocate(t2)
            deallocate(fu)
            deallocate(fv)
            return
        end if

    end do  ! End of timestep loop

end subroutine stepon
