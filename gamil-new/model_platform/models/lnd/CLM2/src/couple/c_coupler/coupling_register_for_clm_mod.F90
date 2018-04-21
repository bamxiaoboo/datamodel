!***************************************************************
!  This file was initially finished by Dr. Li Liu. If you have 
!  any problem, please contact Dr. Li Liu via 
!  liuli-cess@tsinghua.edu.cn
!***************************************************************

#include <misc.h>
#include <preproc.h>


module coupling_register_for_clm_mod 

#if (defined COUP_CAM)
    use clm_varmap , only : numland, landvec
    use pmgrid, only: plon, plond, plat
    use c_coupler_interface_mod
    use clm_varsur, only : pctlak, pctgla 

#if ( defined SPMD )
    use spmdMod    , only : masterproc, iam, npes, proc_landi, proc_landf
#endif



CONTAINS


    subroutine register_lnd_decomp_data_buffers
       use clm_varsur, only: numlon 
       use clm_varpar, only: lsmlat, lsmlon
       implicit none 
       integer, allocatable :: local_cell_indexes(:)
       integer              :: i, j, k

       allocate(local_cell_indexes(lsmlat*lsmlon))

       local_cell_indexes(:) = -1
       if (masterproc) then
          do j = 1,lsmlat
                write(6,*) 'okok1', numlon(j)
             do i = 1,numlon(j)
                local_cell_indexes((j-1)*plon+i) = (j-1)*plon+i
                write(6,*) 'okok', (j-1)*plon+i, numlon(j)
             enddo
          enddo
          call c_coupler_register_decomposition("clm_gamil_grid_scatter_decomp", "gamil_grid", lsmlat*lsmlon, local_cell_indexes)
       else
          call c_coupler_register_decomposition("clm_gamil_grid_scatter_decomp", "gamil_grid", 0, local_cell_indexes)
       endif

       call c_coupler_register_field_info("FRLAKE",   "fraction", "fraction of lake") 
       call c_coupler_register_field_info("FRLANDIC", "fraction", "fraction of land ice") 

       call c_coupler_register_model_data(pctlak,   "clm_gamil_grid_scatter_decomp", "FRLAKE", .true., grid_name="gamil_grid")
       call c_coupler_register_model_data(pctgla,   "clm_gamil_grid_scatter_decomp", "FRLANDIC",  .true., grid_name="gamil_grid")

       deallocate(local_cell_indexes)


    end subroutine register_lnd_decomp_data_buffers


#endif

end module coupling_register_for_clm_mod
