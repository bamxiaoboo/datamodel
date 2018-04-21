#include <misc.h>
#include <params.h>

!-----------------------------------------------------------------------
!
! Purpose: Entry point for GAMIL
!
!-----------------------------NOTICE------------------------------------
!
! Gridpoint Atmosphere Model of IAP/LASG, version 1.0
!
! Purpose:
!
!   Call initialization, time-stepping, and finalization routines.
!
!-----------------------------------------------------------------------

program gamil

    use pmgrid
    use dycore
    use history,      only: bldfld, intht
    use units
    use phys_buffer ! added by SHI Xiangjun and LIU Li
    use ppgrid,       only: pcols, pverp, begchunk, endchunk
    use comsrf,       only: fld_kvh ! added by LIU Li
    use phys_grid,    only: phys_grid_init
    use comfm1!,       only: initialize_comfm1 ! WAN Hui 2003/10/23
    use prognostics
    use mpishorthand, only: mpicom, nsend, nrecv, nwsend, nwrecv

    use mpi_gamil
    use c_coupler_interface_mod
    use restart_mod,    only: do_restart_read
    use register_decompositions_mod
    use register_all_variables_mod 
    use coupling_chemistry_model_mod

    implicit none

#include <comctl.h>
#include <comlun.h>
#include <gpt.inc>
! added by SHI Xiangjun
! DONG Li: try to merge it with the namelist variable
#include <RK_or_MG.h>

#ifdef SUNOS
!#include <floatingpoint.h>
#endif

#ifdef OSF1
#include <for_fpe_flags.f>
    integer(4) old_fpe_flags   ! old settings of floating point exception flags
    integer(4) new_fpe_flags   ! new settings of floating point exception flags
    integer(4) for_set_fpe     ! function to set the floating point exceptions
#endif
    character*8 cdate          ! System date
    character*8 ctime          ! System time
    character*13 filenam
    integer iu
    integer nstep           ! Current timestep number.
    integer kvh_idx ! added by LIU Li
    integer i,j,k,m
    !------------------------------Externals--------------------------------
#if ( defined SUNOS )
    !integer iexcept, ieee_handler
#endif

#ifdef OSF1
    !
    ! Compaq floating point exception handler
    ! Terminate if hit invalid, divide by zero, or overflow.
    !
    new_fpe_flags = FPE_M_TRAP_INV + FPE_M_TRAP_DIV0 + FPE_M_TRAP_OVF
    old_fpe_flags = for_set_fpe(new_fpe_flags)
#endif

#if ( defined SUNOS )
    !
    ! SUN: Trap ieee exceptions for debugging purposes
    !      iexcept = ieee_handler( 'set', 'common', SIGFPE_ABORT )
    !      if ( iexcept /= 0 ) write(6,*)'ieee trapping not supported here'
    !
#endif
    !
    ! Initialize timing library.  2nd arg 0 means disable, 1 means enable
    !
    call t_setoptionf(usrsys, 0)
    call t_initializef

    call t_startf('total')
    call t_startf('initialization')
    !
    ! Initialize internal/external MPI if appropriate
    !
    call c_coupler_initialize(mpicom)
    call register_gamil_surface_subroutines
    !
    ! Initialize SPMD environment if applicable
    !
    call spmdinit
    !
    if (masterproc) then
        write(6, *) "----------------------------------------------------------"
        write(6, *) "    Gridpoint Atmosphere Model of IAP/LASG  (GAMIL)       "
        write(6, *) "                    version 2.0                           "
        write(6, *) "----------------------------------------------------------"
    end if
    !
    ! Fetch and print current date and time
    !
    call datetime(cdate, ctime)
    if (is_rootproc) then
        write(6, *) "DATE ", cdate, " TIME ", ctime
        write(6, *) "----------------------------------------------------------"
        if (dycore_is('EUL')) then
            write(6, *) 'DYCORE is EUL'
        else if (dycore_is('SLD')) then
            write(6, *) 'DYCORE is SLD'
        else if (dycore_is('LR')) then
            write(6, *) 'DYCORE is LR'
        end if
    end if
    !
    ! Set defaults then override with user-specified input
    !
    call preset
    call parse_namelist
    call gamil_2D_decomp()
    !
    ! Define fortran unit numbers
    !
    nsds    = getunit()
    nrg     = getunit()
    nrg2    = getunit()
    luhrest = getunit()

    if (masterproc) then
        write(6, *)
        write(6, "('=========================================')")
        write(6, "('***Summary of Logical Unit assignments***')")
        write(6, *)
        write(6, "('   Restart pointer unit (nsds)     : ', I2)") nsds
        write(6, "('   Master restart unit (nrg)       : ', I2)") nrg
        write(6, "('   Abs/ems unit for restart (nrg2) : ', I2)") nrg2
        write(6, "('   History restart unit (luhrest)  : ', I2)") luhrest
        write(6, "('=========================================')")
        write(6, *)
    end if

    !
    ! Initialize index values for advected and non-advected tracers
    !

    call initindx
    call inital          ! dynamics (mostly) init
    call inti            ! physics init
    call bldfld          ! master field list
    call intht           ! set up history tape contents for this run

    call register_decompositions
    call register_all_variables
    call add_most_flds_for_coupling_chem

    !
    ! Initialize external models or datasets depending upon whether coupled
    !
    call initext

    call do_restart_read
    call t_stopf('initialization')
    !
    ! Invoke driving routine for time integration
    !
    if (RK_or_MG == 'MG') then
        call pbuf_allocate('global') ! added by SHI Xiangjun
        if (c_coupler_is_first_restart_step()) then ! added by LIU Li
            kvh_idx = pbuf_get_fld_idx('KVH')
            do i = begchunk, endchunk
                pbuf(kvh_idx)%fld_ptr(1,1:pcols,1:pverp,i,1) = fld_kvh(1:pcols,1:pverp,i)
            end do
        end if
    end if
    call t_startf('stepon')
    call stepon
    call t_stopf('stepon')
    if (RK_or_MG == 'MG') call pbuf_deallocate('global')  ! added by SHI Xiangjun
    !
    ! End the run cleanly
    !
    call t_stopf('total')
    call t_prf(iam)

    if (masterproc) then
        nstep = c_coupler_get_nstep()
        write (6,9300) nstep-1,nstep
9300    format (//'Number of completed timesteps:',i6,/,'Time step ',i6, &
            ' partially done to provide convectively adjusted and ', &
            'time filtered values for history tape.')
        write(6,*)'------------------------------------------------------------'
        write(6,*)'******* END OF MODEL RUN *******'
    end if

    write(6, *) 'gamil finish run'
    call c_coupler_finalize()

#if ( defined SPMD )
    iu = getunit ()
    write(filenam,'(a10,i3.3)') 'spmdstats.', iam
    open (unit=iu, file=filenam, form='formatted', status='replace')
    write (iu,*)'iam ',iam,' msgs  sent =',nsend
    write (iu,*)'iam ',iam,' msgs  recvd=',nrecv
    write (iu,*)'iam ',iam,' words sent =',nwsend
    write (iu,*)'iam ',iam,' words recvd=',nwrecv
#endif

    stop

end program gamil
