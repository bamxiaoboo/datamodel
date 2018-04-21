# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/mksrfdata/mksrfdatMod.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/mksrfdata/mksrfdatMod.F90"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/mksrfdata/mksrfdatMod.F90" 2

# 1 "./preproc.h" 1






 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/mksrfdata/mksrfdatMod.F90" 2

module mksrfdatMod

!=======================================================================
CONTAINS
!=======================================================================

  subroutine mksrfdat(cam_longxy, cam_latixy, cam_numlon, cam_landfrac, cam_landmask)

!-----------------------------------------------------------------------
!
! Purpose:
! make land model surface dataset from original "raw" data files
!
! Method:
!
! Author: Gordon Bonan, Sam Levis and Mariana Vertenstein
!
!-----------------------------------------------------------------------
! $Id: mksrfdatMod.F90,v 1.12.2.4 2002/08/28 15:41:16 erik Exp $
!-----------------------------------------------------------------------

    use shr_kind_mod, only: r8 => shr_kind_r8
    use infnan
    use clm_varpar                      !parameters
    use clm_varctl                      !run control variables
    use clm_varsur                      !land model grid and fractional land
    use pft_varcon                      !vegetation type (PFT) constants
    use areaMod                         !area averaging routines
    use spmdMod                         !spmd routines and variables
    use surfFileMod                     !write and read surface file
    use histFileMod                     !history file variables
    use mkgridMod                       !land model grid
    use c_coupler_interface_mod
    use fileutils, only : putfil, opnfil, getavu, relavu
    implicit none

! ------------------------ arguments -----------------------------------
    real(r8), optional, intent(in) :: cam_longxy(:,:)   !cam lon values
    real(r8), optional, intent(in) :: cam_latixy(:,:)   !cam lat values
    integer , optional, intent(in) :: cam_numlon(:)     !cam number of longitudes
    real(r8), optional, intent(in) :: cam_landfrac(:,:) !cam fractional land
    integer , optional, intent(in) :: cam_landmask(:,:) !cam land mask
! ----------------------------------------------------------------------

! ------------------------ local variables -----------------------------
    integer :: i,j,k,m                               !indices
    integer :: ier                                   !error status
    integer :: ndiag                                 !unit number for surface data summary
    character(len=256) :: loc_fn                      !local file name
    character(len=256) :: rem_dir                    !mass store file name
    character(len=256) :: rem_fn                     !mass store full path name
    character(len=  8) :: resol                      !resolution for file name
    real(r8) :: sum                                  !sum for error check
    real(r8) :: rmax                                 !maximum patch cover
    logical  :: lremov = .false.                     !true => remove file after dispose
    integer , allocatable :: pft(:,:,:)              !PFT data: PFT values
    real(r8), allocatable :: pctpft(:,:,:)           !PFT data: % of vegetated area for PFTs
    real(r8), allocatable :: pctlnd_pft(:,:)         !PFT data: % land per gridcell
    real(r8), allocatable :: mlai (:,:,:,:)          !monthly lai
    real(r8), allocatable :: msai (:,:,:,:)          !monthly sai
    real(r8), allocatable :: mhgtt(:,:,:,:)          !monthly hgt at top
    real(r8), allocatable :: mhgtb(:,:,:,:)          !monthly hgt at bottom
! ----------------------------------------------------------------------

    if (masterproc) then

       write (6,*) 'Attempting to create surface boundary data .....'
       write (6,'(72a1)') ("-",i=1,60)

! allocate dynamic memory

       allocate(pft(lsmlon,lsmlat,maxpatch_pft), stat=ier)
       if (ier/=0) call endrun
       allocate(pctpft(lsmlon,lsmlat,maxpatch_pft), stat=ier)
       if (ier/=0) call endrun
       allocate(pctlnd_pft(lsmlon,lsmlat), stat=ier)
       if (ier/=0) call endrun

! ----------------------------------------------------------------------
! Open diagnostic output log file
! ----------------------------------------------------------------------

       loc_fn = './surface-data.log'
       ndiag = getavu()
       call opnfil (loc_fn, ndiag, 'f')

# 99 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/mksrfdata/mksrfdatMod.F90"
       write (ndiag,*)'using fractional land data from cam', &
            ' model to create the surface dataset'




       write (ndiag,*) 'PFTs from:         ',trim(mksrf_fvegtyp)
       write (ndiag,*) 'glaciers from:     ',trim(mksrf_fglacier)
       write (ndiag,*) 'urban from:        ',trim(mksrf_furban)
       write (ndiag,*) 'inland water from: ',trim(mksrf_flanwat)
       write (ndiag,*) 'soil texture from: ',trim(mksrf_fsoitex)
       write (ndiag,*) 'soil color from:   ',trim(mksrf_fsoicol)

! ----------------------------------------------------------------------
! Initialize surface variables with unusable values
! ----------------------------------------------------------------------

       soic2d(:,:)   = -999
       sand3d(:,:,:) = 1.e36
       clay3d(:,:,:) = 1.e36
       pctlak(:,:)   = 1.e36
       pctwet(:,:)   = 1.e36
       pcturb(:,:)   = 1.e36
       pctgla(:,:)   = 1.e36
       pft(:,:,:)    = 0
       pctpft(:,:,:) = 0.

! ----------------------------------------------------------------------
! Determine land model grid, fractional land and land mask
! ----------------------------------------------------------------------

! Initialize grid variables with unusable values

       numlon(:)     = 0
       latixy(:,:)   = 1.e36
       longxy(:,:)   = 1.e36
       landmask(:,:) = -999
       landfrac(:,:) = 1.e36




       call mkgrid_cam(cam_longxy, cam_latixy, cam_numlon, cam_landfrac, cam_landmask)



! ----------------------------------------------------------------------
! Make PFTs [pft, pctpft] from dataset [fvegtyp] (1/2 degree PFT data)
! ----------------------------------------------------------------------

       call mkpft (mksrf_fvegtyp, ndiag,  noveg, pctlnd_pft, pft, pctpft)

! ----------------------------------------------------------------------
! Make inland water [pctlak, pctwet] from Cogley's one degree data [flanwat]
! ----------------------------------------------------------------------

       call mklanwat (mksrf_flanwat, ndiag,  pctlak, pctwet)

! ----------------------------------------------------------------------
! Make glacier fraction [pctgla] from [fglacier] dataset
! ----------------------------------------------------------------------

       call mkglacier (mksrf_fglacier, ndiag, pctgla)

! ----------------------------------------------------------------------
! Make soil texture [sand3d, clay3d] from IGBP 5 minute data [fsoitex]
! ----------------------------------------------------------------------

       call mksoitex (mksrf_fsoitex, ndiag, pctgla, sand3d, clay3d)

! ----------------------------------------------------------------------
! Make soil color classes [soic2d] from BATS T42 data [fsoicol]
! ----------------------------------------------------------------------

       call mksoicol (mksrf_fsoicol, ndiag, pctgla, soic2d)

! ----------------------------------------------------------------------
! Make LAI and SAI from 1/2 degree data
! ----------------------------------------------------------------------

       allocate(mlai (lsmlon,lsmlat,maxpatch_pft,12), stat=ier)
       if (ier/=0) call endrun
       allocate(msai (lsmlon,lsmlat,maxpatch_pft,12), stat=ier)
       if (ier/=0) call endrun
       allocate(mhgtt(lsmlon,lsmlat,maxpatch_pft,12), stat=ier)
       if (ier/=0) call endrun
       allocate(mhgtb(lsmlon,lsmlat,maxpatch_pft,12), stat=ier)
       if (ier/=0) call endrun

       call mklai (mksrf_flai, ndiag, pft, mlai, msai, mhgtt, mhgtb)

! ----------------------------------------------------------------------
! Make urban fraction [pcturb] from [furban] dataset
! ----------------------------------------------------------------------

       call mkurban (mksrf_furban, ndiag, pcturb)

! ----------------------------------------------------------------------
! Set LAND values on Ross ice shelf to glacier
! ----------------------------------------------------------------------

       do j = 1,lsmlat
          do i = 1,numlon(j)
             if (latixy(i,j) < -79. .and. landmask(i,j) == 1) then
                soic2d(i,j)   = 0
                pctlak(i,j)   = 0.
                pctwet(i,j)   = 0.
                pcturb(i,j)   = 0.
                pctgla(i,j)   = 100.
                pft(i,j,1)    = noveg
                pctpft(i,j,1) = 100.
                do k = 1,nlevsoi
                   sand3d(i,j,k) = 0.
                   clay3d(i,j,k) = 0.
                end do
                do m = 2,maxpatch_pft
                   pft(i,j,m) = noveg
                   pctpft(i,j,m) = 0.
                end do
             end if
          end do
       end do

! ----------------------------------------------------------------------
! Assume 100% wetland where there is a significant missmatch between the
! land mask and the pft dataset land mask. Also assume medium soil
! color (4) and loamy texture.
! ----------------------------------------------------------------------

       do j = 1,lsmlat
          do i = 1,numlon(j)
             if (landmask(i,j)==1 .and. nint(pctlnd_pft(i,j))==0) then
                soic2d(i,j)   = 4
                pctlak(i,j)   = 0.
                pctwet(i,j)   = 100.
                pcturb(i,j)   = 0.
                pctgla(i,j)   = 0.
                pctpft(i,j,1) = 100.
                pft(i,j,1)    = noveg
                do k = 1, nlevsoi
                   sand3d(i,j,k) = 43.
                   clay3d(i,j,k) = 18.
                end do
                do m = 2,maxpatch_pft
                   pctpft(i,j,m) = 0.
                   pft(i,j,m)    = noveg
                end do
             end if
          end do
       end do

! ----------------------------------------------------------------------
! If have pole points on grid - set south pole to glacier
! north pole is as assumed as non-land
! ----------------------------------------------------------------------


       if (pole_points) then
          do i = 1,numlon(1)
             soic2d(i,1)   = 0
             pctlak(i,1)   = 0.
             pctwet(i,1)   = 0.
             pcturb(i,1)   = 0;
             sand3d(i,1,:) = 0.
             clay3d(i,1,:) = 0.
             pctgla(i,1)   = 100.
             pft   (i,1,:             ) = noveg
             pctpft(i,1,1             ) = 100.
             pctpft(i,1,2:maxpatch_pft) = 0.
          end do
       endif


! ----------------------------------------------------------------------
! Truncate all percentage fields on output grid. This is needed to
! insure that wt is not nonzero (i.e. a very small number such as
! 1e-16) wehre it really should be zero
! ----------------------------------------------------------------------

       do j = 1,lsmlat
          do i = 1,numlon(j)
             do k = 1,nlevsoi
                sand3d(i,j,k) = float(nint(sand3d(i,j,k)))
                clay3d(i,j,k) = float(nint(clay3d(i,j,k)))
             end do
             pctlak(i,j) = float(nint(pctlak(i,j)))
             pctwet(i,j) = float(nint(pctwet(i,j)))
             pcturb(i,j) = float(nint(pcturb(i,j)))
             pctgla(i,j) = float(nint(pctgla(i,j)))
             do m = 1,maxpatch_pft
                pctpft(i,j,m) = float(nint(pctpft(i,j,m)))
             end do
          end do
       end do

! ----------------------------------------------------------------------
! Make sure sum of land cover types does not exceed 100. If it does,
! subtract excess from most dominant land cover.
! ----------------------------------------------------------------------

       do j = 1, lsmlat
          do i = 1, numlon(j)
             rmax = -9999.
             k    = -9999
             if (pctlak(i,j) > rmax) then
                k = 1
                rmax = pctlak(i,j)
             end if
             if (pctwet(i,j) > rmax) then
                k = 2
                rmax = pctwet(i,j)
             end if
             if (pcturb(i,j) > rmax) then
                k = 3
                rmax = pcturb(i,j)
             end if
             if (pctgla(i,j) > rmax) then
                k = 4
                rmax = pctgla(i,j)
             end if
             sum = pctlak(i,j)+pctwet(i,j)+pcturb(i,j)+pctgla(i,j)
             if (k == -9999) then
                write (6,*) 'MKSRFDAT error: largest patch not found'
                call endrun
             else if (sum > 105.) then
                write (6,*) 'MKSRFDAT error: sum of pctlak, pctwet,', &
                     'pcturb and pctgla is greater than 105%'
                write (6,*)'i,j,pctlak,pctwet,pcturb,pctgla= ', &
                     i,j,pctlak(i,j),pctwet(i,j),pcturb(i,j),pctgla(i,j)
                 pctlak(i,j) = pctlak(i,j)/sum*100.
                 pctwet(i,j) = pctwet(i,j)/sum*100.
                 pcturb(i,j) = pcturb(i,j)/sum*100.
                 pctgla(i,j) = pctgla(i,j)/sum*100.
!ljli2012                call endrun     
             else if (sum > 100.) then
                if (k==1) pctlak(i,j) = pctlak(i,j) - (sum-100.)
                if (k==2) pctwet(i,j) = pctwet(i,j) - (sum-100.)
                if (k==3) pcturb(i,j) = pcturb(i,j) - (sum-100.)
                if (k==4) pctgla(i,j) = pctgla(i,j) - (sum-100.)

             end if
          end do
       end do

! ----------------------------------------------------------------------
! Make sure sum of PFT cover equals 100 for land points. If it does not,
! subtract excess from most dominant PFT.
! ----------------------------------------------------------------------

       do j = 1, lsmlat
          do i = 1, numlon(j)
             rmax = -9999.
             k = -9999
             sum = 0.
             do m = 1, maxpatch_pft
                sum = sum + pctpft(i,j,m)
                if (pctpft(i,j,m) > rmax) then
                   k = m
                   rmax = pctpft(i,j,m)
                end if
             end do
             if (k == -9999) then
                write (6,*) 'MKSRFDAT error: largest PFT patch not found'
                call endrun
             else if (landmask(i,j) == 1) then
                if (sum < 95 .or. sum > 105.) then
                   write (6,*) 'MKSRFDAT error: sum of PFT cover is ',sum
                   call endrun
                else if (sum /= 100.) then
                   pctpft(i,j,k) = pctpft(i,j,k) - (sum-100.)
                endif
             endif
          end do
       end do

! ----------------------------------------------------------------------
! Write and dispose surface data file
! ----------------------------------------------------------------------

       write (resol,'(i4.4,"x",i3.3)') lsmlon,lsmlat
!        resol = '1440x616'
       fsurdat = './surface-data.'//trim(resol)//'.nc'

       call surfwrt(fsurdat, pft, pctpft, mlai, msai, mhgtt, mhgtb)

       write (6,'(72a1)') ("-",i=1,60)
       write (6,'(a46,f5.1,a4,f5.1,a5)') 'land model surface data set successfully created for ', &
            360./lsmlon,' by ',180./lsmlat,' grid'

       if (mss_irt > 0) then
          rem_dir = trim(archive_dir) // '/surf/'
          rem_fn = trim(rem_dir)//'surface-data.'//trim(resol)//'.nc'
          call putfil (fsurdat, rem_fn, mss_wpass, mss_irt, lremov)
       endif

! ----------------------------------------------------------------------
! Close and dispose diagnostic log file
! ----------------------------------------------------------------------

       write (6,*)
       write (6,*) 'Surface data output file = ',trim(fsurdat)
       write (6,*) '   This file contains the land model surface data'
       write (6,*) 'Diagnostic log file      = ',trim(loc_fn)
       write (6,*) '   See this file for a summary of the dataset'
       write (6,*)
       close (ndiag)
       call relavu(ndiag)
       if (mss_irt > 0) then
          rem_dir = trim(archive_dir) // '/surf/'
          rem_fn = trim(rem_dir) // 'surface-data.log'
          call putfil (loc_fn, rem_fn, mss_wpass, mss_irt, lremov)
       endif

! deallocate dynamic memory

       deallocate(pft)
       deallocate(pctpft)
       deallocate(pctlnd_pft)
       deallocate(mlai)
       deallocate(msai)
       deallocate(mhgtt)
       deallocate(mhgtb)

    endif    !end of if-masterproc block



! ----------------------------------------------------------------------
! End run if only making surface dataset
! ----------------------------------------------------------------------

! Note that nestep is determined by the flux coupler and not by
! the namelist for a coupled model run

    if (c_coupler_check_coupled_run_finished()) then
       write (6,*)
       write (6,*)'model stopped because run length is zero'
       call endrun
    end if


! ----------------------------------------------------------------------
! Reset real arrays to 1.e36 and integer arrays to -999 since all
! these arrays will be read back in to insure that bit for bit results
! are obtained for a run where a surface dataset file is generated and
! a run where a surface dataset is read in
! ----------------------------------------------------------------------

    lsmedge(:) = inf
    lats(:) = inf
    lonw(:,:) = inf
    numlon(:) =  -999
    latixy(:,:) = 1.e36
    longxy(:,:) = 1.e36
    landmask(:,:) =  -999
    landfrac(:,:) = 1.e36
    soic2d(:,:) = -999
    sand3d(:,:,:) = 1.e36
    clay3d(:,:,:) = 1.e36
    pctwet(:,:) = 1.e36
    pctlak(:,:) = 1.e36
    pctgla(:,:) = 1.e36
    pcturb(:,:) = 1.e36

    return
  end subroutine mksrfdat

!=======================================================================

end module mksrfdatMod


