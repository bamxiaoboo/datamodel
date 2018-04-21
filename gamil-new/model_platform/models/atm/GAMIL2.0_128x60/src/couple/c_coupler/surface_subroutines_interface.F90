!***************************************************************
!  This is a source file of GAMIL, containing all model 
!  subroutine interfaces which will be registerred into 
!  C-Coupler library. This file was initially finished by
!  Dr. Li Liu. If you have any problem, please contact Dr. Li 
!  Liu via liuli-cess@tsinghua.edu.cn
!***************************************************************


#include <misc.h>


subroutine initialize_online_lsm_ccpl_interface

   use surface_subroutines_from_datafile_mod
   call initialize_online_lsm

end subroutine initialize_online_lsm_ccpl_interface



subroutine calculate_ocn_albedo_ccpl_interface

   use surface_subroutines_from_datafile_mod
   call calculate_ocn_albedo
   
end subroutine calculate_ocn_albedo_ccpl_interface



subroutine calculate_sice_albedo_ccpl_interface

   use surface_subroutines_from_datafile_mod
   call calculate_sice_albedo
   
end subroutine calculate_sice_albedo_ccpl_interface



subroutine initialize_data_ocn_ccpl_interface

   use surface_subroutines_from_datafile_mod
   call initialize_data_ocn

end subroutine initialize_data_ocn_ccpl_interface



subroutine initialize_data_sice_ccpl_interface

   use surface_subroutines_from_datafile_mod
   call initialize_data_sice

end subroutine initialize_data_sice_ccpl_interface



subroutine merge_ts_ccpl_interface 

   use surface_subroutines_from_datafile_mod
   call merge_ts

end subroutine merge_ts_ccpl_interface 



subroutine get_initial_surface_data_from_coupler_ccpl_interface
   use surface_subroutines_with_coupler_mod

   call initialize_with_coupler

end subroutine get_initial_surface_data_from_coupler_ccpl_interface


subroutine get_forcing_data_sst_ccpl_interface

#ifdef DATA_OCN 
    use sst_data,            only: sstint
    !
    ! Time interpolate sst data
    !
    call sstint ()
#endif

end subroutine get_forcing_data_sst_ccpl_interface



subroutine get_forcing_data_sice_ccpl_interface

#ifdef DATA_SICE 
    use ice_data,            only: iceint

    call iceint ()
#endif

end subroutine get_forcing_data_sice_ccpl_interface



subroutine srfflx_state_reset_ccpl_interface

    use comsrf

    call srfflx_state_reset(srfflx_state2d)

end subroutine srfflx_state_reset_ccpl_interface




subroutine update_srf_fractions_ccpl_interface

    use comsrf

    call update_srf_fractions

end subroutine update_srf_fractions_ccpl_interface



subroutine process_surface_data_with_online_lsm_ccpl_interface

#ifdef ONLINE_LSM 
    use comsrf
    use atm_lndMod,     only: atmlnd_drv
    use c_coupler_interface_mod

    implicit none

#include <comsol.h>
#include <comctl.h>

    integer nstep                              ! current timestep number

    nstep = c_coupler_get_nstep()

    if (.not. aqua_planet) then
        !
        ! 3.1 Call land model driving routine
        !
        call t_startf('atmlnd_drv')
        call atmlnd_drv(nstep, iradsw, eccen, obliqr, lambm0, mvelpp, surface_state2d, srfflx_parm2d)
        call t_stopf ('atmlnd_drv')

        call update_srf_fluxes(srfflx_state2d, srfflx_parm2d, landfrac)
    end if
#endif

end subroutine process_surface_data_with_online_lsm_ccpl_interface



subroutine process_surface_data_with_forcing_data_sst_ccpl_interface

#ifdef DATA_OCN 
    use comsrf

    implicit none

    call t_startf('camoce')
    call camoce(surface_state2d, srfflx_parm2d)
    call t_stopf('camoce')
    call update_srf_fluxes(srfflx_state2d, srfflx_parm2d, ocnfrac)
#endif

end subroutine process_surface_data_with_forcing_data_sst_ccpl_interface



subroutine process_surface_data_with_forcing_data_sice_ccpl_interface

#ifdef DATA_SICE 
    use comsrf

    implicit none

    call t_startf('camice')
    call camice(surface_state2d, srfflx_parm2d)
    call t_stopf('camice')
    call update_srf_fluxes(srfflx_state2d, srfflx_parm2d, icefrac)
#endif

end subroutine process_surface_data_with_forcing_data_sice_ccpl_interface



subroutine process_surface_data_with_coupler_ccpl_interface

    use surface_subroutines_with_coupler_mod
    use c_coupler_interface_mod

    implicit none
    
    integer,parameter :: R8  = selected_real_kind(12) ! 8 byte real
#include <comctl.h>
    integer nstep                              ! current timestep number

    nstep = c_coupler_get_nstep()

    if (flxave) then
        !
        ! Average the precipitation input to lsm between radiation calls.
        !
       call average_atm_flux_variables(iradsw, nstep, dosw)
    end if
    !
    ! Send/recv data to/from the csm flux coupler.
    !
    call send_surface_data_to_coupler
    call recv_surface_data_from_coupler

end subroutine process_surface_data_with_coupler_ccpl_interface



subroutine adjust_land_ocn_sice_fraction_ccpl_interface

    use surface_subroutines_with_coupler_mod

    call adjust_land_ocn_sice_fraction

end subroutine adjust_land_ocn_sice_fraction_ccpl_interface
