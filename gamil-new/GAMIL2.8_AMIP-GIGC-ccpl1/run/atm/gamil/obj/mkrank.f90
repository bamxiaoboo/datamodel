# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/mksrfdata/mkrank.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/mksrfdata/mkrank.F90"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/mksrfdata/mkrank.F90" 2

# 1 "./preproc.h" 1






 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/mksrfdata/mkrank.F90" 2

subroutine mkrank (n, a, miss, iv, num)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! return indices of largest [num] values in array [a] 
! 
! Method: 
! 
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------
!
! $Id: mkrank.F90,v 1.2.6.1 2002/06/15 13:50:40 erik Exp $ 
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none

! ------------------------ input variables ------------------------
  integer , intent(in) :: n        !array length
  real(r8), intent(in) :: a(0:n)   !array to be ranked
  integer , intent(in) :: miss     !missing data value
  integer , intent(in) :: num      !number of largest values requested
! -----------------------------------------------------------------

! ------------------------ output variables -----------------------
  integer iv(num)      !index to [num] largest values in array [a]
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------
  real(r8) a_max       !maximum value in array
  integer i            !array index
  real(r8) delmax      !tolerance for finding if larger value
  integer m            !do loop index
  integer k            !do loop index
  logical exclude      !true if data value has already been chosen
! -----------------------------------------------------------------

  delmax = 1.e-06

! -----------------------------------------------------------------
! Find index of largest non-zero number
! -----------------------------------------------------------------

  iv(1) = miss
  a_max = -9999.

  do i = 0, n
     if (a(i)>0. .and. (a(i)-a_max)>delmax) then
        a_max = a(i)
        iv(1)  = i
     end if
  end do

! iv(1) = miss indicates no values > 0. this is an error

  if (iv(1) == miss) then
     write (6,*) 'MKRANK error: iv(1) = missing'
     call endrun
  end if

! -----------------------------------------------------------------
! Find indices of the next [num]-1 largest non-zero number.
! iv(m) = miss if there are no more values > 0
! -----------------------------------------------------------------

  do m = 2, num
     iv(m) = miss
     a_max = -9999.
     do i = 0, n

! exclude if data value has already been chosen

        exclude = .false.
        do k = 1, m-1
           if (i == iv(k)) exclude = .true.
        end do

! if not already chosen, see if it is the largest of the remaining values

        if (.not. exclude) then
           if (a(i)>0. .and. (a(i)-a_max)>delmax) then
              a_max = a(i)
              iv(m)  = i
           end if
        end if
     end do
  end do

  return
end subroutine mkrank


