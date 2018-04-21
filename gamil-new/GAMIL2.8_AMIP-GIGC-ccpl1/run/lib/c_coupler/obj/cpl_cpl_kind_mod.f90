# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/libs/c_coupler/External_Algorithms/cpl_cpl_kind_mod.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/libs/c_coupler/External_Algorithms/cpl_cpl_kind_mod.F90"
!===============================================================================
! CVS $Id: cpl_kind_mod.F90,v 1.2 2003/11/22 00:27:06 tcraig Exp $
! CVS $Source: /fs/cgd/csm/models/CVS.REPOS/shared/csm_share/cpl/cpl_kind_mod.F90,v $
! CVS $Name: ccsm3_0_1_beta14 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: cpl_kind_mod -- F90 kind declarations
!
! !DESCRIPTION:
!   F90 kind declarations.
!
! !REVISION HISTORY:
!     2002-Nov-04 - B. Kauffman - created initial version
!
! !REMARKS:
!   This module does not use the standard cpl6 module variable naming convention
!   because this would results in excessively long variable declarations.
!   ie. we want to see real(R8) and not real(cpl_kind_r8)
!
! !INTERFACE: ------------------------------------------------------------------

module cpl_cpl_kind_mod

! !USES:

   use shr_kind_mod  !  shared kind declaration

   implicit none

   private ! except

! !PUBLIC TYPES: 
 
  ! none

! !PUBLIC MEMBER FUNCTIONS:

  ! none

! !PUBLIC DATA MEMBERS:

  integer,parameter,public :: R16= SHR_KIND_R16 ! 16 byte real
  integer,parameter,public :: R8 = SHR_KIND_R8  ! 8 byte real
  integer,parameter,public :: R4 = SHR_KIND_R4  ! 4 byte real
  integer,parameter,public :: RN = SHR_KIND_RN  ! native/default real
  integer,parameter,public :: I8 = SHR_KIND_I8  ! 8 byte integer
  integer,parameter,public :: I4 = SHR_KIND_I4  ! 4 byte integer
  integer,parameter,public :: IN = SHR_KIND_IN  ! native/default integer

  integer,parameter,public :: CL = SHR_KIND_CL  ! generic "long" char string

!EOP

end module cpl_cpl_kind_mod
