# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/control/time_manager.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/control/time_manager.F90"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/control/time_manager.F90" 2

module time_manager

   use shr_kind_mod, only: r8 => shr_kind_r8
   use dycore, only: dycore_is

   implicit none
   private
   save

! Public methods

   public ::&
      timemgr_preset          ! time manager initialization before namelist input

! Public data for namelist input

   character(len=32), public ::&
      calendar   = 'NO_LEAP'     ! Calendar to use in date calculations.
                                 ! 'NO_LEAP' or 'GREGORIAN'
   integer, parameter :: uninit_int = -999999999
!!
   real(r8),public  ::  dtdy = -9999.9     ! timestep of the dycore in seconds !!(wh 2004.04.14)
!!
! Public data for communicating with modules that don't have 'get' methods.

!=========================================================================================
contains
!=========================================================================================

subroutine timemgr_preset()

! Initialize variables before namelist input.

   implicit none

! Local variables
   character(len=*), parameter :: sub = 'timemgr_preset'
!-----------------------------------------------------------------------------------------

   if (dtdy .eq. -9999.9) then
      if (dycore_is ('EUL')) then
         dtdy   = 240.0                 !!(wh 2004.04.14)
      end if
   end if

end subroutine timemgr_preset
!=========================================================================================
end module time_manager
