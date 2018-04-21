# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/physics/cam1/virtem.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/physics/cam1/virtem.F90"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/physics/cam1/virtem.F90" 2

# 1 "./params.h" 1
# 15 "./params.h"
 
 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/physics/cam1/virtem.F90" 2
subroutine virtem(ncol    ,ncold   ,nver    ,t       ,q       ,zvir    , &
                  tv      )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute the virtual temperature.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: B. Boville
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: ncol  ! number of atmospheric columns
   integer, intent(in) :: ncold ! atmospheric column index dimension
   integer, intent(in) :: nver  ! number of vertical levels in a column

   real(r8) t(ncold,nver)       ! temperature
   real(r8) q(ncold,nver)       ! specific humidity
   real(r8) zvir                ! virtual temperature constant
!
! Output arguments
!
   real(r8), intent(out) :: tv(ncold,nver)      ! virtual temperature
!
!---------------------------Local storage-------------------------------
!
   integer i,k              ! column and level indexes
!
   do k=1,nver
      do i=1,ncol
         tv(i,k) = t(i,k)*(1.0 + zvir*q(i,k))
      end do
   end do
!
   return
end subroutine virtem

