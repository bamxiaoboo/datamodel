# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/register_all_variables_mod.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/couple/c_coupler/register_all_variables_mod.F90"
!***************************************************************
!  This is a source file of GAMIL, which registers all variables 
!  into C-Coupler library. This file was initially finished by
!  Dr. Li Liu. If you have any problem, please contact Dr. Li 
!  Liu via liuli-cess@tsinghua.edu.cn
!***************************************************************


module register_all_variables_mod


    public register_all_variables


contains


    subroutine register_all_variables
       use register_private_variables_mod
       use surface_subroutines_with_coupler_mod

       call register_local_coupling_variables
       call register_static_variables

    end subroutine register_all_variables


end module register_all_variables_mod

