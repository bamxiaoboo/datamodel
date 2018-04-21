# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/stdatm.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/stdatm.F90"
module stdatm

!! (wanhui 2003.10.23)
!-----------------------------------------------------------------------
! Purpose:  variables of the LASG std.atm.
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: plond, beglatexdyn, endlatexdyn, plevstd
   use mpi_gamil, only: ilbnd, ihbnd
   use infnan
   use comfm1

   implicit none

      real(r8) :: tbb (plevstd)         
      real(r8) :: hbb (plevstd)         
      real(r8) :: cbb (plevstd)         
      real(r8) :: dcbb(plevstd)         
      real(r8) :: p00, t00              


CONTAINS

    subroutine initialize_stdatm
!
! Purpose:  Allocate and initialize the arrays of the standard atmosphere.
!

      tbb (:) = inf
      hbb (:) = inf
      cbb (:) = inf
      dcbb(:) = inf
      
      p00     = inf
      t00     = inf

      psb (:,:) = inf
      tsb (:,:) = inf
 
    end subroutine initialize_stdatm

end module stdatm


