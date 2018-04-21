!***************************************************************
!  This is a source file of GAMIL, which is reponsible for 
!  coupling GAMIL and other component models with C-Coupler. 
!  This file was initially finished by Dr. Li Liu, based on the
!  ccsm_msg.F90 source file from NCAR. If you have any problem, 
!  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
!***************************************************************


#include <misc.h>
#include <params.h>


module surface_subroutines_with_coupler_mod

    use shr_kind_mod,   only: r8 => shr_kind_r8                                     ! atmospheric model precision
    use pmgrid,         only: plat, plon, beglat, endlat, plond, masterproc, iam    ! model grid
    use ppgrid,         only: pcols, pver, begchunk, endchunk                       ! physics grid
    use phys_grid,      only: get_ncols_p, nlcols, ngcols, &                        ! physics parallel decomposition
                              get_chunk_owner_p, &
                              read_chunk_from_field, write_field_from_chunk
    use shr_sys_mod,    only: shr_sys_flush                                         ! standardized system subroutines
    use shr_kind_mod,   only: shr_kind_in                                           ! defines real & integer kinds
    use c_coupler_interface_mod
    use infnan
    use comsrf

!--------------------------------------------------------------------------
! NOTE:  if nlon is not the same as nlon_p in phys_grid, this module will
!            need to be modified  -RLJ
!--------------------------------------------------------------------------
    use rgrid,          only: nlon                                                  ! reduced grid
!
    use mpishorthand,   only: mpicom, mpir8, mpiint, mpilog                         ! MPI interface
    use history,        only: outfld                                                ! history output

    implicit none

    public initialize_with_coupler                       ! Initialization
    public send_surface_data_to_coupler                  ! Send information to coupler
    public recv_surface_data_from_coupler                ! Receive information from coupler
    public average_atm_flux_variables                    ! Average flux data for coupling
    public copy_out_coupling_recv_fields
    public register_local_coupling_variables
    public adjust_land_ocn_sice_fraction

    private
    !
    ! Buffer information
    !
    integer, parameter :: atm_output_field_z_atm_bot = 1 
    integer, parameter :: atm_output_field_u_atm_bot = 2
    integer, parameter :: atm_output_field_v_atm_bot = 3
    integer, parameter :: atm_output_field_tbot = 4
    integer, parameter :: atm_output_field_ptem = 5
    integer, parameter :: atm_output_field_shum = 6
    integer, parameter :: atm_output_field_dens = 7
    integer, public, parameter :: atm_output_field_pbot = 8
    integer, parameter :: atm_output_field_pslv = 9
    integer, parameter :: atm_output_field_lwdn = 10
    integer, parameter :: atm_output_field_rainc = 11
    integer, parameter :: atm_output_field_rainl = 12
    integer, parameter :: atm_output_field_snowc = 13
    integer, parameter :: atm_output_field_snowl = 14
    integer, parameter :: atm_output_field_swndr = 15
    integer, parameter :: atm_output_field_swvdr = 16
    integer, parameter :: atm_output_field_swndf = 17
    integer, parameter :: atm_output_field_swvdf = 18
    integer, parameter :: atm_output_field_swnet = 19

    integer, parameter :: atm_input_field_tref = 1
    integer, parameter :: atm_input_field_qref = 2
    integer, parameter :: atm_input_field_avsdr = 3
    integer, parameter :: atm_input_field_anidr = 4
    integer, parameter :: atm_input_field_avsdf = 5
    integer, parameter :: atm_input_field_anidf = 6
    integer, parameter :: atm_input_field_t = 7
    integer, parameter :: atm_input_field_sst = 8
    integer, parameter :: atm_input_field_taux = 9
    integer, parameter :: atm_input_field_tauy = 10
    integer, parameter :: atm_input_field_heat_srf_flux_latent = 11
    integer, parameter :: atm_input_field_heat_srf_flux_sen = 12
    integer, parameter :: atm_input_field_lwup = 13
    integer, parameter :: atm_input_field_evap = 14

    integer, parameter :: grid_field_lon = 1
    integer, parameter :: grid_field_lat = 2
    integer, parameter :: grid_field_area = 3
    integer, parameter :: grid_field_total = 3

    !
    ! send/recv buffers
    !
    integer, private, parameter :: nsnd=19              ! number of send variables
    integer, private, parameter :: nrcv=14              ! number of recv variables
    real(r8), allocatable, private :: sbuf(:,:)         ! array for holding grid data to be sent to coupler

    real(r8), allocatable :: recv2d_chunk(:,:,:)        ! chunked recv array
    real(r8), public, allocatable :: send2d_chunk(:,:,:)        ! chunked send array
    !
    ! flux accumulator
    !
    integer, private :: countfa      ! counter for flux accumulators
    !
    ! Surface data that needs to be averaged
    !
    real(r8), allocatable:: precca(:,:)   ! average convective precipitation
    real(r8), allocatable:: precla(:,:)   ! average large-scale precipation
    real(r8), allocatable:: precsca(:,:)  ! average convective snow-fall
    real(r8), allocatable:: precsla(:,:)  ! average large-scale snow-fall
    real(r8), allocatable:: rainconv(:,:) ! convective rainfall
    real(r8), allocatable:: rainlrsc(:,:) ! large-scale rainfall
    real(r8), allocatable:: snowconv(:,:) ! convective snowfall
    real(r8), allocatable:: snowlrsc(:,:) ! larse-scale snowfall
    real(r8), allocatable:: prc_err(:,:)  ! error in precipitation sent to coupler

    integer albshift              ! albedo calculation time shift

contains


    subroutine allocate_local_coupling_variables

       implicit none
       integer sizebuf         ! size of buffer for sending grid data to coupler
       integer i,j


       sizebuf=0
       do j=1,plat
          do i=1,nlon(j)
            if(get_chunk_owner_p(i,j) .eq. iam) then
	      sizebuf=sizebuf+1
            end if
          enddo
       enddo

       allocate(sbuf(sizebuf,grid_field_total))
       allocate(send2d_chunk(pcols,begchunk:endchunk,nsnd))
       allocate(recv2d_chunk(pcols,begchunk:endchunk,nrcv))
       allocate(precca(pcols,begchunk:endchunk))
       allocate(precla(pcols,begchunk:endchunk))
       allocate(precsca(pcols,begchunk:endchunk))
       allocate(precsla(pcols,begchunk:endchunk))
       allocate(rainconv(pcols,begchunk:endchunk))
       allocate(rainlrsc(pcols,begchunk:endchunk))
       allocate(snowconv(pcols,begchunk:endchunk))
       allocate(snowlrsc(pcols,begchunk:endchunk))
       allocate(prc_err(pcols,begchunk:endchunk))
       
       precca  (:,:) = inf
       precla  (:,:) = inf
       precsca (:,:) = inf
       precsla (:,:) = inf
       snowconv(:,:) = inf
       snowlrsc(:,:) = inf
       rainconv(:,:) = inf
       rainlrsc(:,:) = inf
       prc_err (:,:) = inf

    end subroutine allocate_local_coupling_variables



    subroutine register_local_coupling_variables

       implicit none
#include <comlun.h>
#include <comctl.h>

       call allocate_local_coupling_variables

       call c_coupler_register_model_data(albshift, "NULL", "albedo_sec_shift", .true.)
       call c_coupler_register_model_data(countfa, "NULL", "gamil_countfa", .true.)
       call c_coupler_register_model_data(sbuf(:,grid_field_area),"gamil_gamil_grid_decomp", "areac", .true.)
       call c_coupler_register_model_data(sbuf(:,grid_field_lon),"gamil_gamil_grid_decomp", "lon", .true.)
       call c_coupler_register_model_data(sbuf(:,grid_field_lat),"gamil_gamil_grid_decomp", "lat", .true.)
       call c_coupler_register_model_data(send2d_chunk(:,:,atm_output_field_z_atm_bot), "gamil_2D_decomp_phys", "z_atm_bot", .true.)
       call c_coupler_register_model_data(send2d_chunk(:,:,atm_output_field_u_atm_bot), "gamil_2D_decomp_phys", "u_atm_bot", .true.)
       call c_coupler_register_model_data(send2d_chunk(:,:,atm_output_field_v_atm_bot), "gamil_2D_decomp_phys", "v_atm_bot", .true.)
       call c_coupler_register_model_data(send2d_chunk(:,:,atm_output_field_tbot), "gamil_2D_decomp_phys", "tbot", .true.)
       call c_coupler_register_model_data(send2d_chunk(:,:,atm_output_field_ptem), "gamil_2D_decomp_phys", "ptem", .true.)
       call c_coupler_register_model_data(send2d_chunk(:,:,atm_output_field_pbot), "gamil_2D_decomp_phys", "pbot", .true.)
       call c_coupler_register_model_data(send2d_chunk(:,:,atm_output_field_pslv), "gamil_2D_decomp_phys", "pslv", .true.)
       call c_coupler_register_model_data(send2d_chunk(:,:,atm_output_field_shum), "gamil_2D_decomp_phys", "shum", .true.)
       call c_coupler_register_model_data(send2d_chunk(:,:,atm_output_field_dens), "gamil_2D_decomp_phys", "dens", .true.)
       call c_coupler_register_model_data(send2d_chunk(:,:,atm_output_field_swnet), "gamil_2D_decomp_phys", "swnet", .true.)
       call c_coupler_register_model_data(send2d_chunk(:,:,atm_output_field_lwdn), "gamil_2D_decomp_phys", "lwdn", .true.)
       call c_coupler_register_model_data(send2d_chunk(:,:,atm_output_field_rainc), "gamil_2D_decomp_phys", "rainc", .true.)
       call c_coupler_register_model_data(send2d_chunk(:,:,atm_output_field_rainl), "gamil_2D_decomp_phys", "rainl", .true.)
       call c_coupler_register_model_data(send2d_chunk(:,:,atm_output_field_snowc), "gamil_2D_decomp_phys", "snowc", .true.)
       call c_coupler_register_model_data(send2d_chunk(:,:,atm_output_field_snowl), "gamil_2D_decomp_phys", "snowl", .true.)
       call c_coupler_register_model_data(send2d_chunk(:,:,atm_output_field_swndr), "gamil_2D_decomp_phys", "swndr", .true.)
       call c_coupler_register_model_data(send2d_chunk(:,:,atm_output_field_swvdr), "gamil_2D_decomp_phys", "swvdr", .true.)
       call c_coupler_register_model_data(send2d_chunk(:,:,atm_output_field_swndf), "gamil_2D_decomp_phys", "swndf", .true.)
       call c_coupler_register_model_data(send2d_chunk(:,:,atm_output_field_swvdf), "gamil_2D_decomp_phys", "swvdf", .true.)
       call c_coupler_register_model_data(recv2d_chunk(:,:,atm_input_field_taux), "gamil_2D_decomp_phys", "wtaux", .true.)
       call c_coupler_register_model_data(recv2d_chunk(:,:,atm_input_field_tauy), "gamil_2D_decomp_phys", "wtauy", .true.)
       call c_coupler_register_model_data(recv2d_chunk(:,:,atm_input_field_heat_srf_flux_latent), "gamil_2D_decomp_phys", "heat_srf_flux_latent", .true.)
       call c_coupler_register_model_data(recv2d_chunk(:,:,atm_input_field_heat_srf_flux_sen), "gamil_2D_decomp_phys", "heat_srf_flux_sen", .true.)
       call c_coupler_register_model_data(recv2d_chunk(:,:,atm_input_field_lwup), "gamil_2D_decomp_phys", "lwup", .true.)
       call c_coupler_register_model_data(recv2d_chunk(:,:,atm_input_field_evap), "gamil_2D_decomp_phys", "evap", .true.)
       call c_coupler_register_model_data(recv2d_chunk(:,:,atm_input_field_avsdr), "gamil_2D_decomp_phys", "avsdr", .true.)
       call c_coupler_register_model_data(recv2d_chunk(:,:,atm_input_field_anidr), "gamil_2D_decomp_phys", "anidr", .true.)
       call c_coupler_register_model_data(recv2d_chunk(:,:,atm_input_field_avsdf), "gamil_2D_decomp_phys", "avsdf", .true.)
       call c_coupler_register_model_data(recv2d_chunk(:,:,atm_input_field_anidf), "gamil_2D_decomp_phys", "anidf", .true.)
       call c_coupler_register_model_data(recv2d_chunk(:,:,atm_input_field_t), "gamil_2D_decomp_phys", "t", .true.)
       call c_coupler_register_model_data(recv2d_chunk(:,:,atm_input_field_sst), "gamil_2D_decomp_phys", "sst", .true.)
       call c_coupler_register_model_data(snowhland, "gamil_2D_decomp_phys", "snowh", .true.)
       call c_coupler_register_model_data(landfrac, "gamil_2D_decomp_phys", "lfrac", .true.)
       call c_coupler_register_model_data(icefrac, "gamil_2D_decomp_phys", "ifrac", .true.)
       call c_coupler_register_model_data(ocnfrac, "gamil_2D_decomp_phys", "ofrac", .true.)
       call c_coupler_register_model_data(recv2d_chunk(:,:,atm_input_field_tref), "gamil_2D_decomp_phys", "tref", .true.)
       call c_coupler_register_model_data(recv2d_chunk(:,:,atm_input_field_qref), "gamil_2D_decomp_phys", "qref", .true.)
       call c_coupler_register_model_data(precca, "gamil_2D_decomp_phys", "gamil_precca", .true.)
       call c_coupler_register_model_data(precla, "gamil_2D_decomp_phys", "gamil_precla", .true.)
       call c_coupler_register_model_data(precsca, "gamil_2D_decomp_phys", "gamil_precsca", .true.)
       call c_coupler_register_model_data(precsla, "gamil_2D_decomp_phys", "gamil_precsla", .true.)

    end subroutine register_local_coupling_variables



    subroutine initialize_with_coupler

        use physconst,      only: stebol
        use constituents,   only: pcnst, pnats

#include <comctl.h>

        integer i, m, lchnk, n  ! indices
        integer ncols           ! number of columns
        integer ierr            ! allocation error signal

        call t_startf('coupling initialize')
        !
        ! for now set all tracer fluxes to zero
        !
        do m = 2, pcnst+pnats
            do lchnk = begchunk, endchunk
                ncols = get_ncols_p(lchnk)
                do i=1,ncols
                    srfflx_state2d(lchnk)%cflx(i,m) = 0.
                end do
            end do
        end do
        !
        ! require the short and longwave radiation frequencies to match, since these
        ! fluxes will be sent as instantaneous fluxes to the coupler, valid over the
        ! next interval.
        !
        if (masterproc) then
            if (flxave) then
                if (iradsw == iradlw) then
                    write(6, *) 'coupling initialize: coupling will take place every ',iradsw, ' steps'
                else
                    write(6, *) 'coupling initialize: iradsw != iradlw ', iradsw, iradlw
                    call endrun("coupling initialize: bad irad")
                end if
            else
                write(6, *) 'coupling initialize: coupling will take place every time step'
            end if
            call shr_sys_flush(6)
        end if
        call send_atm_grid_to_coupler
        
        !
        ! initial run only: get albedos and ice fraction
        !
        if (c_coupler_is_first_step()) then
            call get_initial_albedo_from_coupler
            !
            ! Initial run only: determine longwave up flux from the surface temperature.
            !
            do lchnk = begchunk, endchunk
                ncols = get_ncols_p(lchnk)
                do i = 1, ncols
                    srfflx_state2d(lchnk)%lwup(i) = stebol*(srfflx_state2d(lchnk)%ts(i)**4)
                end do
            end do
        end if

        call t_stopf('coupling initialize')

        return
    end subroutine initialize_with_coupler



    subroutine recv_surface_data_from_coupler

        use comsrf, only: srfflx_state2d, icefrac, ocnfrac, landfrac, &
                          snowhice, snowhland, verify_fractions

#include <comctl.h>

        integer i, lat, n, lchnk    ! indices
        integer ncols               ! Number of columns
        integer len                 ! temporary variable length

        real(r8) newlandfrac        ! land fraction computed as residual of icefrac + ocnfrac
        real(r8) delta              ! land fraction difference across timesteps
                                    ! code needs to be rewritten to guarantee that this is zero

        call t_startf('receiv coupling data')
        !
        ! get data from flux coupler.
        !
        if (c_coupler_get_nstep() == 0) then
           call c_coupler_execute_procedure("init_recv", "initialize")
        else
           call c_coupler_execute_procedure("kernel_recv", "kernel")
        endif

        call copy_out_coupling_recv_fields

        call t_stopf('receiv coupling data')

        return
    end subroutine recv_surface_data_from_coupler



    subroutine send_surface_data_to_coupler

        use comsrf, only: surface_state2d, srfflx_state2d

#include <comctl.h>

        integer i, lchnk, n, lat    ! indices
        integer ncols               ! Number of columns
        integer len                 ! temporary length variable

        call t_startf('send coupling data')
        !
        ! Divide total precipitation and snowfall into rain and snowfall
        !
        if (flxave) then
            do lchnk = begchunk, endchunk
                ncols = get_ncols_p(lchnk)
                do i = 1, ncols
                    rainconv(i,lchnk) = ((precca(i,lchnk)-precsca(i,lchnk)))*1000.
                    rainlrsc(i,lchnk) = ((precla(i,lchnk)-precsla(i,lchnk)))*1000.
                    snowconv(i,lchnk) = precsca(i,lchnk)*1000.
                    snowlrsc(i,lchnk) = precsla(i,lchnk)*1000.
                end do
            end do
        else
            do lchnk = begchunk, endchunk
                ncols = get_ncols_p(lchnk)
                do i = 1, ncols
                    rainconv(i,lchnk) = ((surface_state2d(lchnk)%precc(i)-surface_state2d(lchnk)%precsc(i)))*1000.
                    rainlrsc(i,lchnk) = ((surface_state2d(lchnk)%precl(i)-surface_state2d(lchnk)%precsl(i)))*1000.
                    snowconv(i,lchnk) = surface_state2d(lchnk)%precsc(i)*1000.
                    snowlrsc(i,lchnk) = surface_state2d(lchnk)%precsl(i)*1000.
                end do
            end do
        end if

        do lchnk = begchunk, endchunk
            ncols = get_ncols_p(lchnk)
            prc_err(1:ncols,lchnk)  = 0.
        end do
        !
        ! Copy from component arrays into one chunk array.
        ! Note that coupler has convention that fluxes are positive downward.
        !
        do lchnk = begchunk, endchunk
            ncols = get_ncols_p(lchnk)
            do i = 1, ncols
                send2d_chunk(i,lchnk,atm_output_field_z_atm_bot)     = surface_state2d(lchnk)%zbot(i) ! Atmospheric state variable m
                send2d_chunk(i,lchnk,atm_output_field_u_atm_bot)     = surface_state2d(lchnk)%ubot(i) ! Atmospheric state variable m/s
                send2d_chunk(i,lchnk,atm_output_field_v_atm_bot)     = surface_state2d(lchnk)%vbot(i) ! Atmospheric state variable m/s
                send2d_chunk(i,lchnk,atm_output_field_tbot)  = surface_state2d(lchnk)%tbot(i) ! Atmospheric state variable K
                send2d_chunk(i,lchnk,atm_output_field_ptem)  = surface_state2d(lchnk)%thbot(i)! Atmospheric state variable K
                send2d_chunk(i,lchnk,atm_output_field_pbot)  = surface_state2d(lchnk)%pbot(i) ! Atmospheric state variable Pa
                send2d_chunk(i,lchnk,atm_output_field_pslv)  = psl(i,lchnk)                   ! Atmospheric state variable Pa
                send2d_chunk(i,lchnk,atm_output_field_shum)  = surface_state2d(lchnk)%qbot(i) ! Atmospheric state variable kg/kg
                send2d_chunk(i,lchnk,atm_output_field_dens)  = rho(i,lchnk)                   ! Atmospheric state variable kg/m^3
                send2d_chunk(i,lchnk,atm_output_field_swnet) = surface_state2d(lchnk)%srfrad(i) - surface_state2d(lchnk)%flwds(i)    ! Atmospheric flux W/m^2
                send2d_chunk(i,lchnk,atm_output_field_lwdn)  = surface_state2d(lchnk)%flwds(i)! Atmospheric flux W/m^2
                send2d_chunk(i,lchnk,atm_output_field_rainc) = rainconv(i,lchnk)              ! Atmospheric flux kg/s/m^2
                send2d_chunk(i,lchnk,atm_output_field_rainl) = rainlrsc(i,lchnk)              ! Atmospheric flux kg/s/m^2
                send2d_chunk(i,lchnk,atm_output_field_snowc) = snowconv(i,lchnk)              ! Atmospheric flux kg/s/m^2
                send2d_chunk(i,lchnk,atm_output_field_snowl) = snowlrsc(i,lchnk)              ! Atmospheric flux kg/s/m^2
                send2d_chunk(i,lchnk,atm_output_field_swndr) = surface_state2d(lchnk)%soll(i) ! Atmospheric flux W/m^2
                send2d_chunk(i,lchnk,atm_output_field_swvdr) = surface_state2d(lchnk)%sols(i) ! Atmospheric flux W/m^2
                send2d_chunk(i,lchnk,atm_output_field_swndf) = surface_state2d(lchnk)%solld(i)! Atmospheric flux W/m^2
                send2d_chunk(i,lchnk,atm_output_field_swvdf) = surface_state2d(lchnk)%solsd(i)! Atmospheric flux W/m^2
            end do
        end do
        !
        ! Output to history file the snow and rain actually sent to coupler as well as the
        ! error between what is sent and what is reported on history file in PRECT/PRECS
        !
        do lchnk = begchunk, endchunk
            call outfld('CPLRAINC', rainconv(1,lchnk), pcols, lchnk)
            call outfld('CPLRAINL', rainlrsc(1,lchnk), pcols, lchnk)
            call outfld('CPLSNOWC', snowconv(1,lchnk), pcols, lchnk)
            call outfld('CPLSNOWL', snowlrsc(1,lchnk), pcols, lchnk)
            call outfld('CPLPRCER', prc_err (1,lchnk), pcols, lchnk)
        end do
        !
        ! Send buffer to coupler
        !
        call msgsnd

        call t_stopf('send coupling data')

        return
    end subroutine send_surface_data_to_coupler



    subroutine copy_out_coupling_recv_fields

        use comsrf, only: srfflx_state2d, icefrac, ocnfrac, landfrac, &
                          snowhland, snowhice, verify_fractions

        implicit none
#include <comlun.h>
#include <comctl.h>

        integer i, n, lat, lchnk     ! indices
        integer ncols                ! Number of columns
        integer ierr            ! allocation error signal
        real    newlandfrac


        if (c_coupler_is_model_data_renewed_in_current_time_step(recv2d_chunk(:,:,atm_input_field_taux))) then
            do lchnk = begchunk, endchunk
                srfflx_state2d(lchnk)%wsx(:)    = -recv2d_chunk(:,lchnk,atm_input_field_taux)   ! Atmosphere-surface flux
            end do
        end if

        if (c_coupler_is_model_data_renewed_in_current_time_step(recv2d_chunk(:,:,atm_input_field_tauy))) then
            do lchnk = begchunk, endchunk
                srfflx_state2d(lchnk)%wsy(:)    = -recv2d_chunk(:,lchnk,atm_input_field_tauy)   ! Atmosphere-surface flux
            end do
        end if

        if (c_coupler_is_model_data_renewed_in_current_time_step(recv2d_chunk(:,:,atm_input_field_lwup))) then
            do lchnk = begchunk, endchunk
                srfflx_state2d(lchnk)%lwup(:)   = -recv2d_chunk(:,lchnk,atm_input_field_lwup)   ! Atmosphere-surface flux
            end do
        end if

        if (c_coupler_is_model_data_renewed_in_current_time_step(recv2d_chunk(:,:,atm_input_field_avsdr))) then
            do lchnk = begchunk, endchunk
                srfflx_state2d(lchnk)%asdir(:)  =  recv2d_chunk(:,lchnk,atm_input_field_avsdr)  ! Surface state variable
            end do
        end if

        if (c_coupler_is_model_data_renewed_in_current_time_step(recv2d_chunk(:,:,atm_input_field_anidr))) then
            do lchnk = begchunk, endchunk
                srfflx_state2d(lchnk)%aldir(:)  =  recv2d_chunk(:,lchnk,atm_input_field_anidr)  ! Surface state variable
            end do
        end if

        if (c_coupler_is_model_data_renewed_in_current_time_step(recv2d_chunk(:,:,atm_input_field_avsdf))) then
            do lchnk = begchunk, endchunk
                srfflx_state2d(lchnk)%asdif(:)  =  recv2d_chunk(:,lchnk,atm_input_field_avsdf)  ! Surface state variable
            end do
        end if

        if (c_coupler_is_model_data_renewed_in_current_time_step(recv2d_chunk(:,:,atm_input_field_anidf))) then
            do lchnk = begchunk, endchunk
                srfflx_state2d(lchnk)%aldif(:)  =  recv2d_chunk(:,lchnk,atm_input_field_anidf)  ! Surface state variable
            end do
        end if

        if (c_coupler_is_model_data_renewed_in_current_time_step(recv2d_chunk(:,:,atm_input_field_t))) then
            do lchnk = begchunk, endchunk
                srfflx_state2d(lchnk)%ts(:)     =  recv2d_chunk(:,lchnk,atm_input_field_t)       ! Surface state variable
            end do
        end if

        if (c_coupler_is_model_data_renewed_in_current_time_step(recv2d_chunk(:,:,atm_input_field_sst))) then
            do lchnk = begchunk, endchunk
                srfflx_state2d(lchnk)%sst(:)    =  recv2d_chunk(:,lchnk,atm_input_field_sst)     ! Surface state variable
            end do
        end if

        if (c_coupler_is_model_data_renewed_in_current_time_step(recv2d_chunk(:,:,atm_input_field_evap))) then
            do lchnk = begchunk, endchunk
                srfflx_state2d(lchnk)%cflx(:,1) = -recv2d_chunk(:,lchnk,atm_input_field_evap)   ! Atmosphere-surface flux
            end do
        end if

        if (c_coupler_is_model_data_renewed_in_current_time_step(recv2d_chunk(:,:,atm_input_field_heat_srf_flux_latent))) then
            do lchnk = begchunk, endchunk
                srfflx_state2d(lchnk)%lhf(:)    = -recv2d_chunk(:, lchnk,atm_input_field_heat_srf_flux_latent)   ! Atmosphere-surface flux
            end do
        end if

        if (c_coupler_is_model_data_renewed_in_current_time_step(recv2d_chunk(:,:,atm_input_field_heat_srf_flux_sen))) then
            do lchnk = begchunk, endchunk
                srfflx_state2d(lchnk)%shf(:)    = -recv2d_chunk(:, lchnk,atm_input_field_heat_srf_flux_sen)   ! Atmosphere-surface flux
            end do
        end if

        if (c_coupler_is_model_data_renewed_in_current_time_step(recv2d_chunk(:,:,atm_input_field_tref))) then
            do lchnk = begchunk, endchunk
                srfflx_state2d(lchnk)%tref(:)   =  recv2d_chunk(:,lchnk,atm_input_field_tref)    ! Surface state variable
            end do
        end if

        snowhice(:,:) = 0.0

        return
    end subroutine copy_out_coupling_recv_fields



    subroutine adjust_land_ocn_sice_fraction

        use comsrf, only: srfflx_state2d, icefrac, ocnfrac, landfrac, &
                          snowhland, snowhice, verify_fractions

        implicit none
#include <comlun.h>
#include <comctl.h>

        integer i, n, lat, lchnk     ! indices
        integer ncols                ! Number of columns
        integer ierr            ! allocation error signal
        real    newlandfrac

        if (c_coupler_is_model_data_renewed_in_current_time_step(icefrac) .or. &
            c_coupler_is_model_data_renewed_in_current_time_step(ocnfrac)) then
        do lchnk = begchunk, endchunk
            ncols = get_ncols_p(lchnk)
            do i = 1, ncols
                newlandfrac = 1.0-icefrac(i,lchnk)-ocnfrac(i,lchnk)
                landfrac(i,lchnk) = newlandfrac
                if (icefrac(i,lchnk)+landfrac(i,lchnk) > 1.0) then
                    icefrac(i,lchnk) = 1.0-landfrac(i,lchnk)
                end if
                ocnfrac(i,lchnk) = 1.0-landfrac(i,lchnk)-icefrac(i,lchnk)
            end do
            !
            ! Ensure that fractions are valid
            !
            call verify_fractions(lchnk, ncols)
        end do
        end if

    end subroutine adjust_land_ocn_sice_fraction 



    subroutine msgsnd

#include <comctl.h>

        integer n                      ! count indices
        integer nstep                  ! current time step
        integer nstepcsm               ! time step sent to flux coupler
        logical nextsw                 ! set to true for next sw calculation
        real(r8) dtime                 ! timestep size

        nstep = c_coupler_get_nstep()
        dtime = c_coupler_get_step_size()
        !
        ! Determine time step sent to flux coupler and corresponding date.
        !
        if (nstep==0) then
            nstepcsm = nstep
        else
            nstepcsm = nstep - 1
        end if
        !
        ! Determine albedo calculation time shift, which is the time interval
        ! from nstepcsm until the next short wave calculation.
        !
        if (nstep /= 0) then
            if (flxave) then
                albshift = nint((nstep+iradsw-nstepcsm)*dtime)
            else
                nextsw = .false.
                n = 1
                do while (.not. nextsw)
                    nextsw = (mod((nstep+n-1),iradsw)==0)
                    if (nextsw) albshift = nint((nstep+n-nstepcsm)*dtime)
                    n = n+1
                end do
            end if
        else
            albshift = nint(iradsw*dtime)+dtime
        end if
        !
        ! Send data to coupler.
        !
        if (nstep == 0) then
           call c_coupler_execute_procedure("init_send", "initialize")
        else
           call c_coupler_execute_procedure("kernel_send", "kernel")
        endif

        return
    end subroutine msgsnd



    subroutine average_atm_flux_variables(iradsw, nstep, dosw)

    use comsrf, only: surface_state2d

!------------------------------Arguments--------------------------------
    integer, intent(in) :: iradsw  ! solar radiation interval
    integer, intent(in) ::  nstep  ! time step number
    logical, intent(in) ::  dosw   ! time to compute averages (solar radiation time)
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
    integer i,lchnk  ! indices
    integer ncols    ! Number of columns
    real(r8) rcount  ! reciprocal of count
!-----------------------------------------------------------------------
!
! If iradsw == 1, then no averaging is required
!
    if (iradsw == 1) return
!
! Set the counter and normalizing factor
!
    if (nstep == 0) countfa = 0
    countfa = countfa + 1
    if (dosw) then
       rcount = 1./countfa
    end if

    do lchnk=begchunk,endchunk
       ncols = get_ncols_p(lchnk)
       if (countfa == 1) then
          do i = 1, ncols
             precca(i,lchnk)  = surface_state2d(lchnk)%precc(i)
             precla(i,lchnk)  = surface_state2d(lchnk)%precl(i)
             precsca(i,lchnk) = surface_state2d(lchnk)%precsc(i)
             precsla(i,lchnk) = surface_state2d(lchnk)%precsl(i)
          end do
!
! Final call of averaging interval, complete averaging and copy data back
!
       else if (dosw) then
          do i = 1, ncols
             precca(i,lchnk)  = rcount*(precca(i,lchnk) + surface_state2d(lchnk)%precc(i))
             precla(i,lchnk)  = rcount*(precla(i,lchnk) + surface_state2d(lchnk)%precl(i))
             precsca(i,lchnk) = rcount*(precsca(i,lchnk) + surface_state2d(lchnk)%precsc(i))
             precsla(i,lchnk) = rcount*(precsla(i,lchnk) + surface_state2d(lchnk)%precsl(i))
          end do
!
! Intermediate call, add data to accumulators
!
       else
          do i = 1, ncols
             precca(i,lchnk)  = precca(i,lchnk) + surface_state2d(lchnk)%precc(i)
             precla(i,lchnk)  = precla(i,lchnk) + surface_state2d(lchnk)%precl(i)
             precsca(i,lchnk) = precsca(i,lchnk) + surface_state2d(lchnk)%precsc(i)
             precsla(i,lchnk) = precsla(i,lchnk) + surface_state2d(lchnk)%precsl(i)
          end do
       end if
    end do
!
! Reset the counter if the average was just computed
!
    if (dosw) then
       countfa = 0
    end if

    return
  end subroutine average_atm_flux_variables



  subroutine send_atm_grid_to_coupler

    use infnan
    use commap, only: latdeg, londeg
    use dycore, only: dycore_is

#include <comctl.h>

!--------------------------Local Variables------------------------------
    integer lat, lon, i, j, n     ! loop indices
    integer nstep                  ! current time step
    integer(SHR_KIND_IN) ::  mask(plon,plat)       ! Mask of valid data
    real(r8) area(plon,plat)      ! Area in radians squared for each grid point
    real(r8) clondeg(plon,plat)   ! Longitude grid
    real(r8) clatdeg(plon,plat)   ! latitude grid as 2 dimensional array
    real(r8) ns_vert(4,plon,plat) ! latitude grid vertices
    real(r8) ew_vert(4,plon,plat) ! longitude grid vertices
    real(r8) del_theta            ! difference in latitude at a grid point
    real(r8) del_phi              ! difference in longitude at a grid point
    real(r8) pie                  ! mathmatical constant 3.1415...
    real(r8) degtorad             ! convert degrees to radians
!-----------------------------------------------------------------------


       nstep = c_coupler_get_nstep()
!
! Constants
!
       pie       = acos(-1.)
       degtorad  = pie / 180.0
!
! Mask for which cells are active and inactive and 2D latitude grid
!
       mask(:,:)    = 0        ! Initialize mask so that cells are inactive
       do lat = 1, plat
         mask(1:nlon(lat),lat)    = 1     ! Active cells
         clatdeg(1:nlon(lat),lat) = latdeg(lat) ! Put latitude in 2D array
         clondeg(1:nlon(lat),lat) = londeg(1:nlon(lat),lat)
       end do
!
! Send vertices of each grid point
! Verticies are ordered as follows:
! 1=lower left, 2 = upper left, 3 = upper right, 4 = lower right
!
! Longitude vertices
!
       do lat = 1, plat
         ew_vert(1,1,lat)             = (londeg(1,lat) - 360.0 + londeg(nlon(lat),lat))*0.5
         ew_vert(1,2:nlon(lat),lat)   = (londeg(1:nlon(lat)-1,lat) + &
                                         londeg(2:nlon(lat),lat))*0.5
         ew_vert(2,:nlon(lat),lat)    = ew_vert(1,:nlon(lat),lat)  ! Copy lowleft corner to upleft
         ew_vert(3,:nlon(lat)-1,lat)  = ew_vert(1,2:nlon(lat),lat)
         ew_vert(3,nlon(lat),lat)     = (londeg(nlon(lat),lat) + (360.0 + londeg(1,lat)))*0.5
         ew_vert(4,:nlon(lat),lat)    = ew_vert(3,:nlon(lat),lat)  ! Copy lowright corner to upright
       end do
!
! Latitude
!
       if ( dycore_is('LR') )then
         ns_vert(1,:nlon(1),1)         = -90.0 + (latdeg(1) - latdeg(2))*0.5
         ns_vert(2,:nlon(plat),plat)   =  90.0 + (latdeg(plat) - latdeg(plat-1))*0.5
       else
         ns_vert(1,:nlon(1),1)         = -90.0
         ns_vert(2,:nlon(plat),plat)   =  90.0
       end if
       ns_vert(4,:nlon(1),1)       = ns_vert(1,nlon(1),1)        ! Copy lower left to lower right
       ns_vert(3,:nlon(plat),plat) = ns_vert(2,nlon(plat),plat)  ! Copy up left to up right
       do lat = 2, plat
         ns_vert(1,:nlon(lat),lat) = (latdeg(lat) + latdeg(lat-1) )*0.5
         ns_vert(4,:nlon(lat),lat) = ns_vert(1,:nlon(lat),lat)
       end do
       do lat = 1, plat-1
         ns_vert(2,:nlon(lat),lat) = (latdeg(lat) + latdeg(lat+1) )*0.5
         ns_vert(3,:nlon(lat),lat) = ns_vert(2,:nlon(lat),lat)
       end do
!
! Get area of grid cells (as radians squared)
!
       area(:,:) = 0.0
       do lat = 1, plat
         do lon = 1, nlon(lat)
           del_phi = sin( ns_vert(2,lon,lat)*degtorad ) - sin( ns_vert(1,lon,lat)*degtorad )
           del_theta = ( ew_vert(4,lon,lat) - ew_vert(1,lon,lat) )*degtorad
           area(lon,lat) = del_theta*del_phi
         end do
       end do
!
! If grid has a pole point (as in Lin-Rood dynamics
!
      if ( dycore_is('LR') )then
         lat = 1
!         mask(2:nlon(lat),lat) = 0   ! Only active one point on pole
         do lon = 1, nlon(lat)
           del_phi = -sin( latdeg(lat)*degtorad ) + sin( ns_vert(2,lon,lat)*degtorad )
           del_theta = ( ew_vert(4,lon,lat) - ew_vert(1,lon,lat) )*degtorad
           area(lon,lat)  = del_theta*del_phi
         end do
         lat = plat
!         mask(2:nlon(lat),lat) = 0   ! Only active one point on pole
         do lon = 1, nlon(lat)
           del_phi =  sin( latdeg(lat)*degtorad ) - sin( ns_vert(1,lon,lat)*degtorad )
           del_theta = ( ew_vert(4,lon,lat) - ew_vert(1,lon,lat) )*degtorad
           area(lon,lat)  = del_theta*del_phi
         end do
       end if
       if ( abs(sum(area) - 4.0*pie) > 1.e-12 )then
         write (6,*) ' sum of areas on globe does not = 4*pi'
         write (6,*) ' sum of areas = ', sum(area)
         call endrun
       end if

! load in the lat, lon, area, mask, and compute gridpoint numbers for
! points on this processor
       n=0
       do j=1,plat
          do i=1,nlon(j)
            if(get_chunk_owner_p(i,j) .eq. iam) then
	      n=n+1
	      sbuf(n,grid_field_lon) = clondeg(i,j)
	      sbuf(n,grid_field_lat) = clatdeg(i,j)
	      sbuf(n,grid_field_area) = area(i,j)
            end if
          enddo
       enddo
!
       call shr_sys_flush(6)

!       call c_coupler_check_grid_values("gamil_gamil_grid_decomp", "gamil_grid", "lon", sbuf(:,grid_field_lon))
       call c_coupler_execute_procedure("connect_stage", "initialize")

    return
  end subroutine send_atm_grid_to_coupler



  subroutine get_initial_albedo_from_coupler

    use comsrf, only: srfflx_state2d, icefrac, ocnfrac, landfrac, &
                      snowhice, snowhland, verify_fractions

#include <comctl.h>

!--------------------------Local Variables------------------------------
    integer i,m,n,lat,lchnk             ! indices
    integer ncols                       ! Number of columns
!-----------------------------------------------------------------------
! send a dummy message to coupler which is expecting an initial data message.
! Coupler will proceed to the albedo init portion.  After it sends the albedo
! to us, it will wait for the true initial message from physpkg.
!
! Receive merged surface state from flux coupler.
!
    call c_coupler_execute_procedure("check_stage", "initialize")

    call copy_out_coupling_recv_fields

        return
    end subroutine get_initial_albedo_from_coupler

end module surface_subroutines_with_coupler_mod

