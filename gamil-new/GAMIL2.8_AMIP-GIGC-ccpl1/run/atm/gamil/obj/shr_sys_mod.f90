# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/libs/shr/shr_sys_mod.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/libs/shr/shr_sys_mod.F90"
!===============================================================================
! CVS $Id: shr_sys_mod.F90,v 1.2 2003/11/22 00:27:09 tcraig Exp $
! CVS $Source: /fs/cgd/csm/models/CVS.REPOS/shared/csm_share/shr/shr_sys_mod.F90,v $
! CVS $Name: ccsm3_0_1_beta14 $
!===============================================================================

MODULE shr_sys_mod

   use shr_kind_mod, only: SHR_KIND_IN, SHR_KIND_R8, SHR_KIND_I8  ! defines real & integer kinds
   use shr_mpi_mod          ! MPI layer

   implicit none

! PUBLIC: Public interfaces

   private
   public :: shr_sys_system, shr_sys_chdir, &
             shr_sys_getenv, shr_sys_flush, &
             shr_sys_abort , shr_sys_irtc,  &
             shr_sys_sleep 

CONTAINS

!===============================================================================

SUBROUTINE shr_sys_system(str,rcode)

   IMPLICIT none

   !----- arguments ---
   character(len=*)    ,intent(in)  :: str    ! system/shell command string
   integer(SHR_KIND_IN),intent(out) :: rcode  ! function return error code






   !----- functions -----
   integer(SHR_KIND_IN),external    :: system ! function to envoke shell command


!-------------------------------------------------------------------------------
! PURPOSE: an architecture independant system call
!-------------------------------------------------------------------------------















   rcode = system(str)







END SUBROUTINE shr_sys_system

!===============================================================================

SUBROUTINE shr_sys_chdir(path, rcode)

   IMPLICIT none

   !----- arguments -----
   character(len=*)    ,intent(in)  :: path    ! chdir to this dir
   integer(SHR_KIND_IN),intent(out) :: rcode   ! return code

   !----- local -----
   integer(SHR_KIND_IN)             :: lenpath ! length of path

   integer(SHR_KIND_IN),external    :: chdir   ! AIX system call


!-------------------------------------------------------------------------------
! PURPOSE: an architecture independant system call
!-------------------------------------------------------------------------------

   lenpath=len_trim(path)










   rcode=chdir(path(1:lenpath))







END SUBROUTINE shr_sys_chdir

!===============================================================================

SUBROUTINE shr_sys_getenv(name, val, rcode)

   IMPLICIT none

   !----- arguments -----
   character(len=*)    ,intent(in)  :: name    ! env var name
   character(len=*)    ,intent(out) :: val     ! env var value
   integer(SHR_KIND_IN),intent(out) :: rcode   ! return code

   !----- local -----
   integer(SHR_KIND_IN)             :: lenname ! length of env var name
   integer(SHR_KIND_IN)             :: lenval  ! length of env var value
   character(len=512)               :: tmpval  ! temporary env var value

!-------------------------------------------------------------------------------
! PURPOSE: an architecture independant system call
!-------------------------------------------------------------------------------

   lenname=len_trim(name)






   call getenv(trim(name),tmpval)
   val=trim(tmpval)
   rcode = 0
   if (len_trim(val) ==  0 ) rcode = 1
   if (len_trim(val) > 512 ) rcode = 2







END SUBROUTINE shr_sys_getenv

!===============================================================================

SUBROUTINE shr_sys_flush(unit)

   IMPLICIT none

   !----- arguments -----
   integer(SHR_KIND_IN) :: unit  ! flush output buffer for this unit

!-------------------------------------------------------------------------------
! PURPOSE: an architecture independant system call
!-------------------------------------------------------------------------------


   call flush(unit)









END SUBROUTINE shr_sys_flush

!===============================================================================

SUBROUTINE shr_sys_abort(string,rcode)

   IMPLICIT none

   character*(*),optional :: string    ! error message string
   integer(SHR_KIND_IN),optional :: rcode   ! error code

   !----- local -----
   integer(SHR_KIND_IN) :: ierr
   logical              :: flag

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(shr_sys_abort) ',a)"

!-------------------------------------------------------------------------------
! PURPOSE: consistent stopping mechanism
!-------------------------------------------------------------------------------

   call shr_sys_flush(6)
   if (len_trim(string) > 0) write(6,F00) 'ERROR: '//trim(string)
   write(6,F00) 'WARNING: calling shr_mpi_abort() and stopping'
   call shr_sys_flush(6)
   call shr_mpi_initialized(flag)
   if (flag) then
     if (present(string).and.present(rcode)) then
       call shr_mpi_abort(trim(string),rcode)
     elseif (present(string)) then
       call shr_mpi_abort(trim(string))
     elseif (present(rcode)) then
       call shr_mpi_abort(rcode=rcode)
     else
       call shr_mpi_abort()
     endif
   endif
   call shr_sys_flush(6)
   call abort()
   stop

END SUBROUTINE shr_sys_abort

!===============================================================================

integer(SHR_KIND_I8) FUNCTION shr_sys_irtc( rate )

   IMPLICIT none
   !----- optional output argument -----
   integer(SHR_KIND_I8), optional :: rate

   !----- local -----
   integer(SHR_KIND_IN)          :: count
   integer(SHR_KIND_IN)          :: count_rate
   integer(SHR_KIND_IN)          :: count_max

   integer(SHR_KIND_IN),save :: last_count = -1
   integer(SHR_KIND_I8),save :: count_offset = 0

!-------------------------------------------------------------------------------
! emulates Cray/SGI irtc function (returns clock tick since last reboot)
!-------------------------------------------------------------------------------
   call system_clock(count=count,count_rate=count_rate, count_max=count_max)
   if ( present(rate) ) rate = count_rate
   shr_sys_irtc = count
!
! System clock is a 24-hour clock, so must check each time pass over midnight
!
   if ( last_count /= -1 )then
     if ( count < last_count ) count_offset = count_offset + count_max
   end if
   shr_sys_irtc = shr_sys_irtc + count_offset
   last_count = count

END FUNCTION shr_sys_irtc

!===============================================================================

SUBROUTINE shr_sys_sleep(sec)

   IMPLICIT none

   !----- input -----
   real   (shr_kind_r8),intent(in) :: sec  ! number of seconds to sleep

   !----- local -----
   integer(shr_kind_in) :: isec   ! integer number of seconds
   integer(shr_kind_in) :: rcode  ! return code
   character(len=90) :: sleep_var

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('sleep ',i8 )"

   save

!-------------------------------------------------------------------------------
! PURPOSE: Sleep for approximately sec seconds
!-------------------------------------------------------------------------------

   isec = nint(sec)

   if (isec < 0) then
      write(6,*) 'ERROR: seconds must be > 0, sec=',sec
   else if (isec == 0) then
!     Don't consider this an error and don't call system sleep
   else
      write(sleep_var,FMT=F00) isec
      call shr_sys_system( sleep_var, rcode )
   endif

END SUBROUTINE shr_sys_sleep

!===============================================================================

END MODULE shr_sys_mod
