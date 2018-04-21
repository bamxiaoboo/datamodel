!***************************************************************
!  This is a source file of GAMIL, which registers all variables
!  with chemistry model into C-Coupler library for coupling. 
!  This file was initially finished by Dr. Li Liu. If you have 
!  any problem, please contact Dr. Li Liu via 
!  liuli-cess@tsinghua.edu.cn
!***************************************************************


#include <misc.h>
#include <params.h>


module coupling_chemistry_model_mod

    use shr_kind_mod, only: r8 => shr_kind_r8
    use ppgrid
    use phys_grid,    only: read_chunk_from_field, write_field_from_chunk, get_ncols_p
    use pmgrid,       only: masterproc
    use prognostics,  only: ptimelevels, n3, n3m2
    use buffer
    use radae,        only: abstot_3d, absnxt_3d, emstot_3d, initialize_radbuffer
    use comsrf
    use ioFileMod
    use phys_buffer
    use c_coupler_interface_mod
    use surface_subroutines_with_coupler_mod

    implicit none
    !
    ! Public interfaces
    !

    type, private :: fld_container_for_coupling_chem
        character(16)        name
        integer              num_lev
        real(r8), pointer :: fld_buf(:,:,:)
    end type fld_container_for_coupling_chem

    integer, private, parameter                    :: max_num_chem_flds = 128
    integer, private                               :: num_registered_flds_for_chem=0
    type(fld_container_for_coupling_chem), private :: registered_flds_for_chem(max_num_chem_flds)

    real, private,allocatable :: PRECCON_array(:,:) 
    real, private,allocatable :: PRECTOT_array(:,:) 
    real, private,allocatable :: PRECSNO_array(:,:) 
    real, private,allocatable :: RADLWG_array(:,:) 
    real, private,allocatable :: RADSWG_array(:,:) 
    real, private,allocatable :: FRLAKE_array(:,:) 
    real, private,allocatable :: FRLANDIC_array(:,:) 

    interface out_fld_for_coupling_chem ; module procedure &
        out_fld_for_coupling_chem_3D, &
        out_fld_for_coupling_chem_2D, &
        out_fld_for_coupling_chem_2D_lchnk, &
        out_fld_for_coupling_chem_1D_lchnk
    end interface




CONTAINS


    subroutine add_fld_for_coupling_chem(fld_name, units, long_name, num_lev)
        implicit none
        character(len=*), intent(in) :: fld_name      
        character(len=*), intent(in) :: units 
        character(len=*), intent(in) :: long_name
        integer         , intent(in) :: num_lev
        integer                      :: i

        num_registered_flds_for_chem = num_registered_flds_for_chem + 1
        if (num_registered_flds_for_chem .gt. max_num_chem_flds) then
            call c_coupler_abort("GAMIL register too many fields for coupling chemistry model")
        endif
        call c_coupler_register_field_info(fld_name, units, long_name) 
        registered_flds_for_chem(num_registered_flds_for_chem)%name      = fld_name
        registered_flds_for_chem(num_registered_flds_for_chem)%num_lev   = num_lev
        allocate(registered_flds_for_chem(num_registered_flds_for_chem)%fld_buf(pcols,begchunk:endchunk,num_lev))
        registered_flds_for_chem(num_registered_flds_for_chem)%fld_buf(:,:,:) = 0.0
        if (num_lev .eq. 1) then
            call c_coupler_register_model_data(registered_flds_for_chem(num_registered_flds_for_chem)%fld_buf, &
                                               "gamil_2D_decomp_phys",fld_name,.true., grid_name="gamil_grid")
        else if (num_lev .eq. pver) then
            call c_coupler_register_model_data(registered_flds_for_chem(num_registered_flds_for_chem)%fld_buf, &
                                               "gamil_2D_decomp_phys",fld_name,.true., grid_name="gamil_3D_grid_lev26")
            call c_coupler_register_model_data(registered_flds_for_chem(num_registered_flds_for_chem)%fld_buf, &
                                               "gamil_2D_decomp_phys",fld_name,.true., grid_name="gamil_3D_grid_lev26_sigma1")
        else if (num_lev .eq. pverp) then
            call c_coupler_register_model_data(registered_flds_for_chem(num_registered_flds_for_chem)%fld_buf, &
                                               "gamil_2D_decomp_phys",fld_name,.true., grid_name="gamil_3D_grid_lev27")
            call c_coupler_register_model_data(registered_flds_for_chem(num_registered_flds_for_chem)%fld_buf, &
                                               "gamil_2D_decomp_phys",fld_name,.true., grid_name="gamil_3D_grid_lev27_sigma1")
        else if (num_lev .eq. 1) then 
            call c_coupler_register_model_data(registered_flds_for_chem(num_registered_flds_for_chem)%fld_buf, &
                                               "gamil_2D_decomp_phys",fld_name,.true., grid_name="gamil_grid")
        else 
            call c_coupler_abort("number of levels of fields for coupling chemistry model is not supported")
        endif

    end subroutine add_fld_for_coupling_chem



    subroutine copy_fld_for_coupling_chem_3D(field_in, field_out, num_lev)
        implicit none
        real(r8), intent(in)         :: field_in(pcols,begchunk:endchunk,num_lev) 
        real(r8), intent(out)        :: field_out(pcols,begchunk:endchunk,num_lev) 
        integer , intent(in)         :: num_lev


        field_out(:,:,:) = field_in(:,:,:)
 
    end subroutine copy_fld_for_coupling_chem_3D



    subroutine copy_fld_for_coupling_chem_2D(field_in, field_out, num_lev)
        implicit none
        real(r8), intent(in)         :: field_in(pcols,begchunk:endchunk) 
        real(r8), intent(out)        :: field_out(pcols,begchunk:endchunk,num_lev) 
        integer , intent(in)         :: num_lev


        field_out(:,:,1) = field_in(:,:)
 
    end subroutine copy_fld_for_coupling_chem_2D



    subroutine copy_fld_for_coupling_chem_2D_lchnk(field_in, field_out, num_lev, lchnk)
        implicit none
        real(r8), intent(in)         :: field_in(pcols,num_lev) 
        real(r8), intent(out)        :: field_out(pcols,begchunk:endchunk,num_lev) 
        integer , intent(in)         :: num_lev
        integer , intent(in)         :: lchnk


        field_out(:,lchnk,:) = field_in(:,:)
 
    end subroutine copy_fld_for_coupling_chem_2D_lchnk



    subroutine copy_fld_for_coupling_chem_1D_lchnk(field_in, field_out, num_lev, lchnk)
        implicit none
        real(r8), intent(in)         :: field_in(pcols) 
        real(r8), intent(out)        :: field_out(pcols,begchunk:endchunk,num_lev) 
        integer , intent(in)         :: num_lev
        integer , intent(in)         :: lchnk


        field_out(:,lchnk,1) = field_in(:)
 
    end subroutine copy_fld_for_coupling_chem_1D_lchnk



    subroutine search_fld_index(fld_name, indx)
        implicit none
        character(len=*), intent(in) :: fld_name      
        integer,          intent(out) :: indx

        
        do indx = 1, num_registered_flds_for_chem
            if (registered_flds_for_chem(indx)%name == fld_name) then
                goto 200
            endif
        enddo 

200     if (indx .gt. num_registered_flds_for_chem) then
            call c_coupler_abort("field has not been registerred when output it as a for coupling chemistry model")
        endif 

    end subroutine search_fld_index



    subroutine out_fld_for_coupling_chem_3D(fld_name, field_buf)
        implicit none
        character(len=*), intent(in) :: fld_name      
        real(r8), intent(in)         :: field_buf(:,:,:) ! Array containing field values
        integer                      :: indx

        
        call search_fld_index(fld_name, indx)
        call copy_fld_for_coupling_chem_3D(field_buf, registered_flds_for_chem(indx)%fld_buf, &
                                        registered_flds_for_chem(indx)%num_lev)

    end subroutine out_fld_for_coupling_chem_3D



    subroutine out_fld_for_coupling_chem_2D(fld_name, field_buf)
        implicit none
        character(len=*), intent(in) :: fld_name      
        real(r8), intent(in)         :: field_buf(:,:) ! Array containing field values
        integer                      :: indx

        
        call search_fld_index(fld_name, indx)
        if (registered_flds_for_chem(indx)%num_lev .ne. 1) then
            call c_coupler_abort("number of levels of for 2D field has not been registerred correctly")
        endif
        call copy_fld_for_coupling_chem_2D(field_buf, registered_flds_for_chem(indx)%fld_buf, &
                                        registered_flds_for_chem(indx)%num_lev)

    end subroutine out_fld_for_coupling_chem_2D



    subroutine out_fld_for_coupling_chem_1D_lchnk(fld_name, field_buf, lchnk)
        implicit none
        character(len=*), intent(in) :: fld_name      
        real(r8), intent(in)         :: field_buf(:) ! Array containing field values
        integer, intent(in)          :: lchnk
        integer                      :: indx

        
        call search_fld_index(fld_name, indx)
        if (registered_flds_for_chem(indx)%num_lev .ne. 1) then
            call c_coupler_abort("number of levels of for 2D field has not been registerred correctly")
        endif
        call copy_fld_for_coupling_chem_1D_lchnk(field_buf, registered_flds_for_chem(indx)%fld_buf, &
                                        registered_flds_for_chem(indx)%num_lev, lchnk)

    end subroutine out_fld_for_coupling_chem_1D_lchnk



    subroutine out_fld_for_coupling_chem_2D_lchnk(fld_name, field_buf, lchnk)
        implicit none
        character(len=*), intent(in) :: fld_name      
        real(r8), intent(in)         :: field_buf(:,:) ! Array containing field values
        integer, intent(in)          :: lchnk
        integer                      :: indx

        
        call search_fld_index(fld_name, indx)
        if (registered_flds_for_chem(indx)%num_lev .eq. 1) then
            call c_coupler_abort("number of levels of for 3D field has not been registerred correctly")
        endif
        call copy_fld_for_coupling_chem_2D_lchnk(field_buf, registered_flds_for_chem(indx)%fld_buf, &
                                        registered_flds_for_chem(indx)%num_lev, lchnk)

    end subroutine out_fld_for_coupling_chem_2D_lchnk



    subroutine add_most_flds_for_coupling_chem
    implicit none

        allocate(PRECCON_array(pcols,begchunk:endchunk))
        allocate(PRECTOT_array(pcols,begchunk:endchunk))
        allocate(PRECSNO_array(pcols,begchunk:endchunk))
        allocate(RADLWG_array(pcols,begchunk:endchunk))
        allocate(RADSWG_array(pcols,begchunk:endchunk))
        allocate(FRLAKE_array(pcols,begchunk:endchunk))
        allocate(FRLANDIC_array(pcols,begchunk:endchunk))

        call add_fld_for_coupling_chem('CLDF','fraction','Cloud fraction',pver)
        call add_fld_for_coupling_chem('CMFMC','kg m-2 s-1','Moist convection mass flux',pverp)
        call add_fld_for_coupling_chem('DQRCU','kg m-2 s-1','conv precip prod rate',pver)
        call add_fld_for_coupling_chem('DQRLSAN','kg m-2 s-1','LS precip prod rate',pver)
        call add_fld_for_coupling_chem('DQIDTMST','kg kg-1 s-1','ice tendency, mst proc',pver)
        call add_fld_for_coupling_chem('DQLDTMST','kg kg-1 s-1','H2O tendency, mst proc',pver)
        call add_fld_for_coupling_chem('DQVDTMST','kg kg-1 s-1','vapor tendency, mst proc',pver)
        call add_fld_for_coupling_chem('DTRAIN','kg m-2 s-1','detrainment flux',pver)
        call add_fld_for_coupling_chem('MOISTQ','g kg-1 day-1','tendency in sp. C17',pver)
        call add_fld_for_coupling_chem('OPTDEP','1','visible optical depth',pver)
        call add_fld_for_coupling_chem('PV','m2 kg-1 s-1','potential vort',pver)
        call add_fld_for_coupling_chem('QI','kg kg-1','cloud ice mixing ratio',pver)
        call add_fld_for_coupling_chem('QL','kg kg-1','cloud water mixing ratio',pver)
        call add_fld_for_coupling_chem('RH','fraction','relative humidity',pver)
        call add_fld_for_coupling_chem('SPHU','g kg-1','specific humidity',pver)
        call add_fld_for_coupling_chem('T','K','temperature',pver)
        call add_fld_for_coupling_chem('TAUCLI','dimensionless','opt depth of ice clouds',pver)
        call add_fld_for_coupling_chem('TAUCLW','dimensionless','opt depth of H2O cloud',pver)
        call add_fld_for_coupling_chem('U','m -s','E/W component of wind',pver)
        call add_fld_for_coupling_chem('V','m -s','N/S component of wind',pver)
        call add_fld_for_coupling_chem('ALBD','fraction','visible surface albedo',1)
        call add_fld_for_coupling_chem('CLDFRC','fraction','column cloud fraction',1)
        call add_fld_for_coupling_chem('CLDTOPS','fraction','column cloud fraction',1)
        call add_fld_for_coupling_chem('EFLUX','W m-2','latent heat flux',1)
        call add_fld_for_coupling_chem('EVAP','kg m-2 s-1','surface evaporation',1)
        call add_fld_for_coupling_chem("FRLAKE",   "fraction", "fraction of lake", 1) 
        call add_fld_for_coupling_chem('FRLAND','fraction','fraction of land', 1) 
        call add_fld_for_coupling_chem("FRLANDIC", "fraction", "fraction of land ice", 1) 
        call add_fld_for_coupling_chem('FROCEAN','fraction','fraction of ocean',1)
        call add_fld_for_coupling_chem('GRN','fraction','greenness fraction',1)
        call add_fld_for_coupling_chem('GWETROOT','fraction','root zone soil wetness',1)
        call add_fld_for_coupling_chem('GWETTOP','fraction','top soil moisture',1)
        call add_fld_for_coupling_chem('HFLUX','W m-2','sensible heat flux',1)
        call add_fld_for_coupling_chem('LAI','m2 m-2','leaf area index',1)
        call add_fld_for_coupling_chem('PARDR','W m-2','direct photsyn active rad',1)
        call add_fld_for_coupling_chem('PARDF','W m-2','diffuse photsyn active rad',1)
        call add_fld_for_coupling_chem('PBLH','m','PBL height',1)
        call add_fld_for_coupling_chem('PRECCON','kg m-2 s-1','conv precip @ ground',1)
        call add_fld_for_coupling_chem('PRECTOT','kg m-2 s-1','total precip @ ground',1)
        call add_fld_for_coupling_chem('PRECSNO','kg m-2 s-1','snow precip',1)
        call add_fld_for_coupling_chem('PS','Pa','sfc press at timestep',1)
        call add_fld_for_coupling_chem('RADLWG','W m-2','net LW radiation @ ground',1)
        call add_fld_for_coupling_chem('RADSWG','W m-2','solar radiation @ ground',1)
        call add_fld_for_coupling_chem('SLP','Pa','sea level pressure',1)
        call add_fld_for_coupling_chem('SNODP','m','snow depth',1)
        call add_fld_for_coupling_chem('SNOMAS','m','snow mass(total snow storage on the land)',1)
        call add_fld_for_coupling_chem('TROPP','Pa','tropopause pressure',1)
        call add_fld_for_coupling_chem('TS','K','surface temperature',1)
        call add_fld_for_coupling_chem('TSKIN','K','surface skin temperature',1)
        call add_fld_for_coupling_chem('U10M','m s-1','E/W wind speed @ 10m height',1)
        call add_fld_for_coupling_chem('V10M','m s-1','N/S wind speed @ 10m height',1)
        call add_fld_for_coupling_chem('USTAR','m s-1','friction velocity',1)
        call add_fld_for_coupling_chem('Z0','m','surface roughness height',1)

        call c_coupler_register_sigma_grid_bottom_field(send2d_chunk(:,:,atm_output_field_pbot),"gamil_3D_grid_lev26_sigma1")

    end subroutine add_most_flds_for_coupling_chem



    subroutine out_caculated_flds_for_coupling_chem()
        implicit none
        integer :: lchnk, ncols, i
        
        do lchnk = begchunk, endchunk
            ncols = get_ncols_p(lchnk)
            do i = 1, ncols
                PRECCON_array(i,lchnk) = surface_state2d(lchnk)%precc(i)*1000.
                PRECTOT_array(i,lchnk) = (surface_state2d(lchnk)%precl(i)+surface_state2d(lchnk)%precc(i))*1000.
                PRECSNO_array(i,lchnk) = (surface_state2d(lchnk)%precsc(i)+surface_state2d(lchnk)%precsl(i))*1000.
                RADLWG_array(i,lchnk)  = srfflx_state2d(lchnk)%lwup(i)-surface_state2d(lchnk)%flwds(i) 
                RADSWG_array(i,lchnk)  = surface_state2d(lchnk)%srfrad(i)-surface_state2d(lchnk)%flwds(i)
                send2d_chunk(i,lchnk,atm_output_field_pbot)  = surface_state2d(lchnk)%pbot(i) ! Atmospheric state variable Pa
            end do
        end do
        call out_fld_for_coupling_chem('PRECCON',PRECCON_array)
        call out_fld_for_coupling_chem('PRECTOT',PRECTOT_array)
        call out_fld_for_coupling_chem('PRECSNO',PRECSNO_array)
        call out_fld_for_coupling_chem('RADLWG',RADLWG_array)
        call out_fld_for_coupling_chem('RADSWG',RADSWG_array)
        call out_fld_for_coupling_chem('FRLAND',landfrac)

    end subroutine out_caculated_flds_for_coupling_chem


end module coupling_chemistry_model_mod
