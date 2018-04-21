#include <misc.h>
#include <params.h>

!-----------------------------------------------------------------------
!
! Purpose:
!
!   Initialize external models and/or boundary dataset information
!
! Method:
!
! Author:
!
!   CCM Core Group
!
!-----------------------------------------------------------------------





subroutine initext
!!(wh 2003.12.27)

    use shr_kind_mod,   only: r8 => shr_kind_r8
    use pmgrid
    use ppgrid,         only: begchunk, endchunk
    use phys_grid,      only: get_ncols_p, get_rlat_all_p, get_rlon_all_p,get_lat_all_p, get_lon_all_p
    use comsrf
    use rgrid
    use shr_orb_mod
    use ioFileMod
    use so4bnd
    use so4bnd_IPCC ! added by WAN Hui
    use commap
    use filenames,      only: bndtvo, bndtvs
    use physconst,      only: stebol
    use c_coupler_interface_mod
    use mpishorthand
    use surface_subroutines_from_datafile_mod

    implicit none

#include <comlun.h>
#include <comctl.h>
#include <comsol.h>
    include 'netcdf.inc'

    character(256) locfn        ! netcdf local filename to open
    character(4) ncnam(5)
    integer  yr, mon, day, tod  ! components of a date
    real(r8) calday             ! current calendar day
    integer  lchnk
    real(r8) tssav(pcols,begchunk:endchunk) ! cam surface temperatures
    logical  log_print          ! Flag to print out log information or not

    call c_coupler_get_current_calendar_time(calday)
    !
    !----------------------------------------------------------------------
    ! 1. Obtain datasets
    !----------------------------------------------------------------------
    !
    ! Obtain time-variant ozone and sst datatsets and do initial read of
    ! ozone dataset
    !
    if (.not. ideal_phys) then
        if (masterproc) then
            call getfil(bndtvo, locfn)
            call wrap_open(locfn, 0, ncid_oz)
            write(6, "('Notice: initext: ')", advance="no")
            write(6, "('wrap_open returns ncid ', I5)", advance="no") ncid_oz
            write(6, "(' for file ', A)") trim(locfn)
        end if
        call oznini
    end if
    !
    !----------------------------------------------------------------------
    ! 2. Obtain sulfate aerosol datasets
    !----------------------------------------------------------------------
    !
    if (doRamp_so4) then
        call sulfini
    end if

    if (doIPCC_so4) then
        call sulfini_IPCC ! added by WAN Hui
    end if

    if (masterproc) then
        log_print = .true.
    else
        log_print = .false.
    end if
    call shr_orb_params(iyear_AD, eccen, obliq , mvelp, obliqr, lambm0, mvelpp, log_print)

    call srfflx_state_reset(srfflx_state2d)
    call read_land_inidat

    if (adiabatic .or. ideal_phys) then
        icefrac(:pcols,begchunk:endchunk) = 0.0
        call update_srf_fractions
    end if

    call c_coupler_execute_procedure("get_initial_surface_data", "initialize")

    return
end subroutine initext
