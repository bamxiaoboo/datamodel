# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/libs/esmf/src/Infrastructure/TimeMgmt/ESMF_CalendarMod.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/libs/esmf/src/Infrastructure/TimeMgmt/ESMF_CalendarMod.F"
! $Id: ESMF_CalendarMod.F,v 1.1.2.5 2002/04/25 17:14:52 dneckels Exp $
	module ESMF_CalendarMod
!===============================================================================
!BOP
! !MODULE: ESMF_CalendarMod
!
! !USES:
!
! !PUBLIC TYPES:
	implicit none

	type ESMF_Calendar
          private
          sequence
          integer(8) type
          integer(8) dim(13) 
          integer(8) dimRunningSum(13)
          integer(8) diy
	end type ESMF_Calendar
!
! !PUBLIC MEMBER FUNCTIONS:
!
! !PUBLIC DATA MEMBERS:
      integer, parameter :: ESMF_CALENDAR_TYPE_UNDEFINED=0,
     &                      ESMF_NO_LEAP=1,
     &                      ESMF_GREGORIAN=2,
     &                      ESMF_360_DAY=3
!
! !DESCRIPTION:
! The calendar class contains implements GREGORIAN and no-leap-year calendars.
!EOP
!===============================================================================
	
	end module
