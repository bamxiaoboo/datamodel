!***************************************************************
!  This is a source file of GAMIL, which is reponsible for 
!  retart function of GAMIL with C-Coupler library. This file was 
!  initially finished by Dr. Li Liu. If you have any problem, 
!  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
!***************************************************************


module restart_mod

    use c_coupler_interface_mod
    use register_private_variables_mod

    private exchange_dyn_restart_fields_boundary
    public  do_restart_write
    public  do_restart_read

contains

    subroutine do_restart_read
       use mpi_gamil
       use prognostics
       use comfm1
       use pmgrid, only: beglatex,beglatexdyn,endlatexdyn
       implicit none
#include <comlun.h>
#include <comctl.h>
       integer j

       if (.not. c_coupler_is_first_restart_step()) return

       call c_coupler_register_model_data(u3(:,:,:,n3m2),"gamil_2D_decomp_prog","gamil_u3",.true.)
       call c_coupler_register_model_data(v3(:,:,:,n3m2),"gamil_2D_decomp_prog","gamil_v3",.true.)
       call c_coupler_register_model_data(t3(:,:,:,n3m2),"gamil_2D_decomp_prog","gamil_t3",.true.)
       call c_coupler_register_model_data(q3(:,:,:,:,n3m2),"gamil_2D_decomp_prog","gamil_q3",.true.)
       call c_coupler_register_model_data(ps(:,:,n3m2),"gamil_2D_decomp_prog","gamil_ps",.true.)
       call register_phys_dynamic_variables 

       call c_coupler_do_restart_read
       call exchange_dyn_restart_fields_boundary

       call copy_from_phys_io_arrays

       do j=jbeg0, jend0
          nigw(j) = nigw_2D(ibeg0,j)
       end do

       call c_coupler_withdraw_model_data("gamil_2D_decomp_prog","gamil_u3")
       call c_coupler_withdraw_model_data("gamil_2D_decomp_prog","gamil_v3")
       call c_coupler_withdraw_model_data("gamil_2D_decomp_prog","gamil_t3")
       call c_coupler_withdraw_model_data("gamil_2D_decomp_prog","gamil_q3")
       call c_coupler_withdraw_model_data("gamil_2D_decomp_prog","gamil_ps")
       call withdraw_phys_dynamic_variables 

       call c_coupler_advance_timer

    end subroutine do_restart_read



    subroutine do_restart_write
       use mpi_gamil
       use prognostics
       use comfm1
       use pmgrid, only: beglatex,beglatexdyn,endlatexdyn
       implicit none
#include <comlun.h>
#include <comctl.h>
       integer :: j


       call c_coupler_register_model_data(u3(:,:,:,n3m2),"gamil_2D_decomp_prog","gamil_u3",.true.)
       call c_coupler_register_model_data(v3(:,:,:,n3m2),"gamil_2D_decomp_prog","gamil_v3",.true.)
       call c_coupler_register_model_data(t3(:,:,:,n3m2),"gamil_2D_decomp_prog","gamil_t3",.true.)
       call c_coupler_register_model_data(q3(:,:,:,:,n3m2),"gamil_2D_decomp_prog","gamil_q3",.true.)
       call c_coupler_register_model_data(ps(:,:,n3m2),"gamil_2D_decomp_prog","gamil_ps",.true.)
       call register_phys_dynamic_variables 

       do j=jbeg0, jend0
          nigw_2D(:,j) = nigw(j)
       end do
       call copy_to_phys_io_arrays
       call c_coupler_do_restart_write()

       call c_coupler_withdraw_model_data("gamil_2D_decomp_prog","gamil_u3")
       call c_coupler_withdraw_model_data("gamil_2D_decomp_prog","gamil_v3")
       call c_coupler_withdraw_model_data("gamil_2D_decomp_prog","gamil_t3")
       call c_coupler_withdraw_model_data("gamil_2D_decomp_prog","gamil_q3")
       call c_coupler_withdraw_model_data("gamil_2D_decomp_prog","gamil_ps")
       call withdraw_phys_dynamic_variables 

    end subroutine do_restart_write


    subroutine exchange_dyn_restart_fields_boundary
       use mpi_gamil
       use comfm1
       use pmgrid, only: beglatexdyn,endlatexdyn
       implicit none

       call gamil_arrays_comm(COMM_TO_LEFT, 1, u(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, v(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, t(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, q(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, ws(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, wpa(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, ghi(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, pes(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, ghs(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, su(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, sv(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, st(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, du(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, dv(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, dtt(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, dps(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, du0(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, dv0(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, dtt0(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, dps0(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, du1(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, dv1(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, dtt1(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, dps1(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, uu(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, vv(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, tt(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, p(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, ply2(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, uk(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, vk(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, ttk(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, psk(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, tb2(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, cb(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, dcb(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, hps(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, c0(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_LEFT, 1, cb0(:,beglatexdyn,1))

       call gamil_arrays_comm(COMM_TO_RIGHT, 1, u(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, v(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, t(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, q(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, ws(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, wpa(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, ghi(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, pes(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, ghs(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, su(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, sv(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, st(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, du(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, dv(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, dtt(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, dps(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, du0(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, dv0(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, dtt0(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, dps0(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, du1(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, dv1(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, dtt1(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, dps1(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, uu(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, vv(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, tt(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, p(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, ply2(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, uk(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, vk(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, ttk(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, psk(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, tb2(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, cb(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, dcb(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, hps(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, c0(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_RIGHT, 1, cb0(:,beglatexdyn,1))

       call gamil_arrays_comm(COMM_TO_TOP, 1, u(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, v(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, t(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, q(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, ws(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, wpa(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, ghi(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, pes(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_TOP, 1, ghs(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_TOP, 1, su(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, sv(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, st(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, du(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, dv(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, dtt(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, dps(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_TOP, 1, du0(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, dv0(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, dtt0(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, dps0(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_TOP, 1, du1(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, dv1(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, dtt1(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, dps1(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_TOP, 1, uu(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, vv(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, tt(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, p(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_TOP, 1, ply2(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, uk(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, vk(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, ttk(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, psk(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_TOP, 1, tb2(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, cb(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, dcb(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_TOP, 1, hps(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_TOP, 1, c0(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_TOP, 1, cb0(:,beglatexdyn,1))

       call gamil_arrays_comm(COMM_TO_BOT, 1, u(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, v(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, t(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, q(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, ws(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, wpa(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, ghi(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, pes(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_BOT, 1, ghs(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_BOT, 1, su(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, sv(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, st(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, du(:,beglatexdyn,1))

       call gamil_arrays_comm(COMM_TO_BOT, 1, dv(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, dtt(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, dps(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_BOT, 1, du0(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, dv0(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, dtt0(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, dps0(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_BOT, 1, du1(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, dv1(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, dtt1(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, dps1(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_BOT, 1, uu(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, vv(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, tt(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, p(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_BOT, 1, ply2(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, uk(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, vk(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, ttk(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, psk(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_BOT, 1, tb2(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, cb(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, dcb(:,beglatexdyn,1))
       call gamil_arrays_comm(COMM_TO_BOT, 1, hps(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_BOT, 1, c0(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_BOT, 1, cb0(:,beglatexdyn,1))

    end subroutine exchange_dyn_restart_fields_boundary


end module restart_mod
