# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/main/controlMod.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/main/controlMod.F90"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/main/controlMod.F90" 2

# 1 "./preproc.h" 1






 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/main/controlMod.F90" 2

module controlMod

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varpar        !parameter statements
  use clm_varctl        !run control variables
  use histFileMod       !history file variables
  use spmdMod           !spmd routines and variables

  implicit none
  save

! Namelist variables only used locally

  integer :: hist_ndens                         !output density of netcdf history files
  logical :: mkfsurdat                          !true => make surface data from raw data
  character(len=256) :: rpntpath                !full UNIX pathname of restart pointer file
  character(len=  7) :: runtyp(4)               !run type
  character(len  =8) :: hist_fldaux1(maxalflds) !fields for first  auxillary history file
  character(len  =8) :: hist_fldaux2(maxalflds) !fields for second auxillary history file





!=======================================================================
CONTAINS
!=======================================================================

  subroutine control_init (cam_caseid , cam_ctitle, cam_irad , cam_nsrest, &
                           cam_crtinic, cam_nhtfrq, cam_mfilt, cam_irt )

!-----------------------------------------------------------------------
!
! Purpose:
! initialize run control variables
!
! Method:
! When running in cam mode, the base calendar info, nstep, nestep,
! nsrest, and time step are input to the land model from CAM.
! The values in the clmexp namelist are not used. The minimum
! namelist parameters are:
!    o fsurdat
!    o finidat
!    o fpftcon
! When running in offline or csm mode, the minimum namelist parameters are:
!    o nsrest
!    o nestep or nelapse
!    o fsurdat
!    o finidat
!    o dtime
!
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------
! $Id: controlMod.F90,v 1.9.2.9.6.1 2002/10/03 20:07:36 erik Exp $
!-----------------------------------------------------------------------





  use c_coupler_interface_mod


! ------------------------ includes ------------------------------------
    include 'netcdf.inc'
! ----------------------------------------------------------------------

! ------------------------ arguments -----------------------------------
    character(len=*), optional, intent(in) :: cam_caseid    !cam caseid
    character(len=*), optional, intent(in) :: cam_ctitle    !cam title
    integer         , optional, intent(in) :: cam_irad      !cam radiation frequency
    integer         , optional, intent(in) :: cam_nsrest    !cam run type
    character(len=*), optional, intent(in) :: cam_crtinic   !cam initial dataset frequency
    integer         , optional, intent(in) :: cam_nhtfrq    !cam history write freq for tape 1
    integer         , optional, intent(in) :: cam_mfilt     !cam number of files per tape for tape 1
    integer         , optional, intent(in) :: cam_irt       !cam mss retention time
! ----------------------------------------------------------------------

! ------------------------ local variables -----------------------------
    character(len=256) :: homedir   !full UNIX filepath name of home directory
    character(len=256) :: logid     !logid part of file path name
    character(len=256) :: cap       !upper case logid
    character(len=  1) :: ctmp      !character temporary
    integer :: i,j,n                !loop indices
    integer :: iundef               !integer undefined value
    real(r8):: rundef               !real undefined value
    integer :: ierr                 !error code
! -----------------------------------------------------------------



! The following can be read in but are overwritten with values from the
! cam time_manager module - consequently they are only declared as local
! variables here

    integer :: dtime      ! timestep in seconds
    integer :: nestep     ! final timestep (or day if negative) number
    integer :: nelapse    ! number of timesteps (or days if negative) to extend a run
    integer :: start_ymd  ! starting date for run in yearmmdd format
    integer :: start_tod  ! starting time of day for run in seconds
    integer :: stop_ymd   ! stopping date for run in yearmmdd format
    integer :: stop_tod   ! stopping time of day for run in seconds
    integer :: ref_ymd    ! reference date for time coordinate in yearmmdd format
    integer :: ref_tod    ! reference time of day for time coordinate in seconds
    character(len=32) :: calendar ! Calendar in date calculations ('NO_LEAP' or 'GREGORIAN')

! ----------------------------------------------------------------------

! ------------------------ namelist variables --------------------------
    namelist /clmexp/  &
         ctitle, caseid, nsrest,  &
         calendar, dtime, nelapse, nestep, start_ymd, start_tod,  &
         stop_ymd, stop_tod, ref_ymd, ref_tod, &
         nrevsn, rpntpath, hist_ndens, hist_dov2xy, &
         hist_nhtfrq, hist_mfilt, hist_fldadd, &
         hist_chntyp, hist_fldaux1, hist_fldaux2, hist_crtinic, &
         archive_dir, mss_wpass, mss_irt, &
         finidat, fsurdat, fpftcon, frivinp_rtm, offline_atmdir, &
         mksrf_fvegtyp, mksrf_fsoitex, mksrf_fsoicol, mksrf_flanwat, &
         mksrf_fglacier, mksrf_furban, mksrf_flai, mksrf_offline_fgrid, &
         mksrf_offline_edgen, mksrf_offline_edgee, mksrf_offline_edges, &
         mksrf_offline_edgew, mksrf_offline_fnavyoro, &
         conchk, irad, wrtdia, csm_doflxave, rtm_nsteps


! === define run =======================
!
!    o caseid     = 32 character case name
!    o ctitle     = 80 character case title
!    o nsrest     = integer flag. 0: initial run. 1: restart: 3: branch
!
! === model time =======================
!
!    o dtime      = real model time step (s)
!    o calendar   = Calendar to use in date calculations.
!                  'no_leap' (default) or 'gregorian'
!    o start_ymd  = Starting date for run encoded in yearmmdd format.
!                   Default value is read from initial conditions file.
!    o start_tod  = Starting time of day for run in seconds since 0Z.
!                   Default value is read from initial conditions file.
!    o stop_ymd   = Stopping date for run encoded in yearmmdd format.
!                   No default.
!    o stop_tod   = Stopping time of day for run in seconds since 0Z.
!                   Default: 0.
!    o nelapse    = nnn, Specify the ending time for the run as an interval
!                   starting at the current time in either timesteps
!                   (if positive) or days (if negative).
!                   Either nestep or (stop_ymd,stop_tod) take precedence.
!    o nestep     = nnnn, Specify the ending time for the run as an interval
!                   starting at (start_ymd,start_tod) in either timesteps
!                   (if positive) or days (if negative).
!                   (stop_ymd,stop_tod) takes precedence if set.
!    o ref_ymd    = Reference date for time coordinate encoded in yearmmdd format.
!                   Default value is start_ymd.
!    o ref_tod    = Reference time of day for time coordinate in seconds since 0Z.
!                   Default value is start_tod.
!
! === input data ===
!
!    o finidat         = 256 character initial conditions file name
!    o fsurdat         = 256 character surface data file name
!    o fpftcon         = 256 character data file with PFT physiological constants
!    o frivinp_rtm     = 256 character input data file for rtm
!    o nrevsn          = 256 character restart file name for use with branch run
!
! === offline forcing data ===
!
!    o offline_atmdir  = 256 character directory for input atm data files (can be Mass Store)
!
! === input data when making surface data [fsurdat] ===
!
!    o mksrf_offline_fgrid   = offline - land grid dataset to use instead of generating grid
!    o mksrf_offline_fnavyoro= offline - 20 min navy orography dataset
!    o mksrf_offline_edgen   = offline - northern edge of grid (degrees): >  -90 and <= 90
!    o mksrf_offline_edgee   = offline - eastern edge of grid (degrees) : see following notes
!    o mksrf_offline_edges   = offline - southern edge of grid (degrees): >= -90 and <  90
!    o mksrf_offline_edgew   = offline - western edge of grid (degrees) : see following notes
!    o mksrf_fvegtyp         = 256 character vegetation type data file name
!    o mksrf_fsoitex         = 256 character soil texture data file name
!    o mksrf_fsoicol         = 256 character soil color data file name
!    o mksrf_flanwat         = 256 character inland water data file name
!    o mksrf_furban          = 256 character urban data file name
!    o mksrf_fglacier        = 256 character glacier data file name
!    o mksrf_flai            = 256 character lai data file file name
!
! === history and restart files ===
!
!    o hist_ndens    = integer, can have value of 1 (nc_double) or 2 (nf_float)
!    o hist_dov2xy   = true if want grid-average history field (false = vector)
!    o hist_nhtfrq   = integer history interval (+ = iterations,  - = hours, 0=monthly ave)
!    o hist_mfilt    = integer number of time samples per history file
!    o hist_fldadd   = 8 character name of fields to change to active
!    o hist_chntyp   = paired 8 character field name and field type.
!                      OVERRIDES default settings in routine histlst(e.g., 'TV','maximum')
!    o hist_fldaux1  = 8 character name of fields for first  auxillary history file
!    o hist_fldaux2  = 8 character name of fields for second auxillary history file
!    o hist_crtinic  = 8 character frequency to generate initial dataset
!                      can be set to 'MONTHLY', 'YEARLY' or 'NONE'.
!    o rpntpath      = 256 character full UNIX pathname of the local restart pointer
!                      file. This file must exist when the model is restarted. This
!                      file is overwritten and updated every time new restart data
!                      files are output.
!
! === long term archiving =====
!
!    o archive_dir = 256 character long term archive directory (can be MSS directory)
!    o mss_irt     = integer mass store retention period (days)
!    o mss_wpass   = 8 character mass store write password for output data sets
!
! === model physics ===
!
!    o conchk     = true if want error energy and water conservation checks
!    o irad       = integer solar radiation frequency (+ = iteration. - = hour)
!    o wrtdia     = true if want output written
!    o csm_doflxave = true => flux averaging is to be performed (only used for csm mode)
!
! === rtm control variables ===
!
!    o rtm_nsteps  = if > 1, average rtm over rtm_nsteps time steps
!
! ----------------------------------------------------------------------

    if (masterproc) then
       write(6,*) 'Attempting to initialize run control settings .....'
    endif

    runtyp(0 + 1) = 'initial'
    runtyp(1 + 1) = 'restart'
    runtyp(3 + 1) = 'branch '

    iundef = -9999999
    rundef = -9999999.

! ----------------------------------------------------------------------
! Default values
! ----------------------------------------------------------------------

! control variables

    caseid  = ' '
    ctitle  = ' '
    nsrest  = iundef

! initial data

    fsurdat = ' '
    finidat = ' '
    fpftcon = ' '
    frivinp_rtm = ' '
    nrevsn  = ' '

! offline mode

    offline_atmdir   = ' '

! surface generation

    mksrf_offline_fgrid    = ' '
    mksrf_offline_fnavyoro = ' '
    mksrf_offline_edgen    =   90.
    mksrf_offline_edgee    =  180.
    mksrf_offline_edges    =  -90.
    mksrf_offline_edgew    = -180.

    mksrf_fvegtyp  = ' '
    mksrf_fsoitex  = ' '
    mksrf_fsoicol  = ' '
    mksrf_flanwat  = ' '
    mksrf_furban   = ' '
    mksrf_fglacier = ' '
    mksrf_flai     = ' '

! long term archive settings

    archive_dir = ' '
    mss_irt = 0
    mss_wpass = ' '

! history file variables

    hist_ndens  = 1
    hist_dov2xy(1) = .true.
    hist_dov2xy(2:maxhist) = .true.
    hist_nhtfrq(1) = -24
    hist_nhtfrq(2:maxhist) = iundef
    hist_mfilt(1) = 1
    hist_mfilt(2:maxhist)  = iundef
    hist_fldadd(:) = ' '
    hist_chntyp(:,:) = ' '
    hist_fldaux1(:) = ' '
    hist_fldaux2(:) = ' '
    hist_crtinic = 'YEARLY'
    rpntpath = 'not_specified'

! other namelist variables

    irad = -1
    conchk = .true.
    wrtdia = .false.
    csm_doflxave = .true.

! RTM control variables





! ----------------------------------------------------------------------
! Read namelist from standard input. Override if coupled to CAM
! ----------------------------------------------------------------------

    if (masterproc) then
       read(5, clmexp, iostat=ierr)
       if (ierr /= 0) then
          if (masterproc) then
             write(6,*)'error: namelist input resulted in error code ',ierr
          endif
          call endrun
       endif
    endif


    caseid  = cam_caseid
    ctitle  = cam_ctitle
    irad    = cam_irad
    nsrest  = cam_nsrest
    hist_crtinic  = cam_crtinic



    hist_nhtfrq(1) = cam_nhtfrq

    hist_mfilt(1) = cam_mfilt
    mss_irt = cam_irt


# 357 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/main/controlMod.F90"

!if archive directory not input in namelist - set default from caseid

    if (archive_dir == ' ') then
       logid  = ' '
!!     call getenv ('LOGNAME', logid)    !!(2003.10.22)
       call getenv ('HOME'   , logid)
       if (logid(1:1) == ' ') then
          write (6,*) 'error: logname not defined'
          call endrun
       end if
       cap = ' '
       do i = 1, len_trim(logid)
          cap(i:i) = logid(i:i)
          ctmp = cap(i:i)
          if (ichar(logid(i:i))>=97 .and. ichar(logid(i:i))<=122) then
             cap(i:i) = char(ichar(ctmp) - 32)
          endif
       end do
       archive_dir = '/' // trim(cap) // '/csm/' // trim(caseid) // '/lnd'
    end if


    call control_spmd()


! ----------------------------------------------------------------------
! Define run
! ----------------------------------------------------------------------

! Determine run type

    if (nsrest == iundef) then
       if (masterproc) write(6,*) 'error: must set nsrest'
       call endrun
    end if

! ----------------------------------------------------------------------
! Surface data
! ----------------------------------------------------------------------

    if (fsurdat == ' ') then
       mkfsurdat = .true.
    else
       mkfsurdat = .false.
    endif

    if (mkfsurdat) then

       if (mksrf_fvegtyp  == ' ' .or. &
           mksrf_fsoitex  == ' ' .or. &
           mksrf_fsoicol  == ' ' .or. &
           mksrf_flanwat  == ' ' .or. &
           mksrf_furban   == ' ' .or. &
           mksrf_fglacier == ' ' .or. &
           mksrf_flai     == ' ') then
          if (masterproc) then
             write(6,*) 'error: need to set data file name'
             write(6,*)'mksrf_fvegtyp = ',mksrf_fvegtyp
             write(6,*)'mksrf_fsoitex = ',mksrf_fsoitex
             write(6,*)'mksrf_fsoicol = ',mksrf_fsoicol
             write(6,*)'mksrf_flanwat = ',mksrf_flanwat
             write(6,*)'mksrf_furban  = ',mksrf_furban
             write(6,*)'mksrf_fglacier= ',mksrf_fglacier
             write(6,*)'mksrf_flai    = ',mksrf_flai
          endif
          call endrun
       end if

       if (nsrest > 0) then
          if (masterproc) then
             write(6,*) 'error: can not make surface data ', &
                  'during a continuation run'
          endif
          call endrun
       end if

       if (finidat /= ' ') then
          if (masterproc) then
             write(6,*) 'error: can not make surface data ', &
                  'when finidat is already specified'
             write(6,*) 'set finidat to empty string in namelist'
          endif
          call endrun
       end if

    end if

# 458 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/main/controlMod.F90"

! ----------------------------------------------------------------------
! Model physics
! ----------------------------------------------------------------------


!time manager initialization done in cam code
    dtime = c_coupler_get_step_size()







    if (irad < 0) irad = nint(-irad*3600./dtime)

# 484 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/main/controlMod.F90"

! ----------------------------------------------------------------------
! History and restart files
! ----------------------------------------------------------------------

    mss_irt = min(mss_irt,1825)

    fldaux(:,:) = ' '
    if (maxhist-1 > 3) then
       if (masterproc) write(6,*) 'CLM_CTLI error: must create additional fldaux'
       call endrun
    end if
    do j = 1, maxhist-1
       do i = 1, maxalflds
          if (j == 1) fldaux(i,j) = hist_fldaux1(i)
          if (j == 2) fldaux(i,j) = hist_fldaux2(i)
       end do
    end do

    nhist = 0
    do j = 1, maxhist-1
       do i = 1, maxalflds
          if (fldaux(i,j) /= ' ') nhist = j
       end do
    end do
    nhist = nhist + 1

    do i = 1, nhist
       if (hist_mfilt(i) == iundef ) then
          if (masterproc) then
             write(6,*)'error: must set hist_mfilt for file ',i
          endif
          call endrun
       end if
       if (hist_nhtfrq(i) == iundef) then
          if (masterproc) then
             write(6,*)'error: must set hist_nhtfrq for file ',i
          endif
          call endrun
       else if (hist_nhtfrq(i) < 0) then
          hist_nhtfrq(i) = nint(-hist_nhtfrq(i)*SHR_CONST_CDAY/(24.*dtime))
       endif
    end do

    if (hist_ndens == 1) then
       ncprec = nf_double
    else if (hist_ndens == 2) then
       ncprec = nf_float
    else
       if (masterproc) then
          write(6,*)'error: history tape hist_ndens must be 1 or 2'
       endif
       call endrun
    end if

    if (rpntpath == 'not_specified') then
       call getenv ('HOME', homedir)
       rpntpath = './lnd.'//trim(caseid)//'.rpointer'
    endif

    do i = 1, nhist
       if (hist_nhtfrq(i)==0) then
          hist_mfilt(i) = 1
       endif
    end do

    if (nsrest == 0) nrevsn = ' '
    if (nsrest == 1) nrevsn = 'set by restart pointer file file'
    if (nsrest == 3 .and. nrevsn == ' ') then
       if (masterproc) write(6,*) 'error: need to set restart data file name'
       call endrun
    end if

    if (trim(hist_crtinic) /= 'MONTHLY' .and. trim(hist_crtinic) /= 'YEARLY') then
       hist_crtinic = 'NONE'
    endif

! ----------------------------------------------------------------------
! Restart pointer file
! ----------------------------------------------------------------------

! split the full pathname of the restart pointer file into a
! directory name and a file name
! check if the directory exists and if not, make it

    rpntdir = ' '
    rpntfil = ' '
    do n = len_trim(rpntpath),1,-1
       if (rpntpath(n:n) ==  '/') then
          rpntdir = rpntpath(1:n-1)
          rpntfil = rpntpath(n+1:len_trim(rpntpath))
          go to 100
       endif
    enddo
    rpntdir = '.'        ! no "/" found, set path = "."
    rpntfil = rpntpath   ! no "/" found, use whole input string.
100 continue

    if (masterproc) then
       write(6,*) 'Successfully initialized run control settings'
       write(6,*)
    endif

    return
  end subroutine control_init

!=======================================================================



  subroutine control_spmd()

!-----------------------------------------------------------------------
!
! Purpose:
! Distribute namelist data all processors. The cpp  definition
! provides for the funnelling of all program i/o through the master
! processor. Processor 0 either reads restart/history data from the
! disk and distributes it to all processors, or collects data from
! all processors and writes it to disk.
!
! Method:
!
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------





    use mpishorthand

! ------------------------ local variables -----------------------------
    integer ier       !error code
!-----------------------------------------------------------------------

! run control variables

    call mpi_bcast (caseid, len(caseid), mpichar, 0, mpicom, ier)
    call mpi_bcast (ctitle, len(ctitle), mpichar, 0, mpicom, ier)
    call mpi_bcast (nsrest, 1, mpiint, 0, mpicom, ier)

# 639 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/main/controlMod.F90"

! initial file variables

    call mpi_bcast (nrevsn, len(nrevsn), mpichar, 0, mpicom, ier)
    call mpi_bcast (finidat, len(finidat), mpichar, 0, mpicom, ier)
    call mpi_bcast (fsurdat, len(fsurdat), mpichar, 0, mpicom, ier)




! surface dataset generation variables

    if (fsurdat == ' ') then
       call mpi_bcast (mksrf_fvegtyp, len(mksrf_fvegtyp), mpichar, 0, mpicom, ier)
       call mpi_bcast (mksrf_fsoitex, len(mksrf_fsoitex), mpichar, 0, mpicom, ier)
       call mpi_bcast (mksrf_fsoicol, len(mksrf_fsoicol), mpichar, 0, mpicom, ier)
       call mpi_bcast (mksrf_flanwat, len(mksrf_flanwat), mpichar, 0, mpicom, ier)
       call mpi_bcast (mksrf_furban, len(mksrf_furban), mpichar, 0, mpicom, ier)
       call mpi_bcast (mksrf_fglacier, len(mksrf_fglacier), mpichar, 0, mpicom, ier)
       call mpi_bcast (mksrf_flai, len(mksrf_flai), mpichar, 0, mpicom, ier)
    endif

! physics variables

    call mpi_bcast (conchk, 1, mpilog, 0, mpicom, ier)
    call mpi_bcast (irad, 1, mpiint, 0, mpicom, ier)
    call mpi_bcast (csm_doflxave, 1, mpilog, 0, mpicom, ier)
    call mpi_bcast (rtm_nsteps, 1, mpiint, 0, mpicom, ier)
    call mpi_bcast (wrtdia, 1, mpilog, 0, mpicom, ier)

! history and restart file variables

    call mpi_bcast (hist_dov2xy, size(hist_dov2xy), mpilog, 0, mpicom, ier)
    call mpi_bcast (hist_nhtfrq, size(hist_nhtfrq), mpiint, 0, mpicom, ier)
    call mpi_bcast (hist_mfilt, size(hist_mfilt), mpiint, 0, mpicom, ier)
    call mpi_bcast (hist_ndens, 1, mpiint, 0, mpicom, ier)
    call mpi_bcast (hist_chntyp, len(hist_chntyp(1,1))*size(hist_chntyp), mpichar, 0, mpicom, ier)
    call mpi_bcast (hist_fldadd, len(hist_fldadd(1))*size(hist_fldadd), mpichar, 0, mpicom, ier)
    call mpi_bcast (hist_fldaux1, len(hist_fldaux1(1))*size(hist_fldaux1), mpichar, 0, mpicom, ier)
    call mpi_bcast (hist_fldaux2, len(hist_fldaux2(1))*size(hist_fldaux2), mpichar, 0, mpicom, ier)
    call mpi_bcast (hist_crtinic, len(hist_crtinic), mpichar, 0, mpicom, ier)
    call mpi_bcast (rpntpath, len(rpntpath), mpichar, 0, mpicom, ier)

! long term archiving variables

    call mpi_bcast (mss_irt, 1, mpiint, 0, mpicom, ier)
    call mpi_bcast (mss_wpass, len(mss_wpass), mpichar, 0, mpicom, ier)
    call mpi_bcast (archive_dir, len(archive_dir), mpichar, 0, mpicom, ier)

    return
  end subroutine control_spmd


!=======================================================================

  subroutine control_print

!-----------------------------------------------------------------------
!
! Purpose:
! Write out run control variables
!
! Method:
!
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------

! ------------------------ local variables -----------------------------
    integer i  !loop index
!-----------------------------------------------------------------------

    write(6,*) 'define run:'
    write(6,*) '   run type              = ',runtyp(nsrest+1)
    write(6,*) '   case title            = ',trim(ctitle)
    write(6,*) 'input data files:'
    write(6,*) '   PFT physiology = ',trim(fpftcon)
    if (mkfsurdat) then
       write(6,*) '   generated surface dataset using raw data'
       write(6,*) '     plant types  = ',trim(mksrf_fvegtyp)
       write(6,*) '     inland water = ',trim(mksrf_flanwat)
       write(6,*) '     glacier      = ',trim(mksrf_fglacier)
       write(6,*) '     urban        = ',trim(mksrf_furban)
       write(6,*) '     soil texture = ',trim(mksrf_fsoitex)
       write(6,*) '     soil color   = ',trim(mksrf_fsoicol)
       write(6,*) '     lai and sai  = ',trim(mksrf_flai)
# 738 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/main/controlMod.F90"
    else
       write(6,*) '   surface data   = ',trim(fsurdat)
    end if
    if (nsrest == 0 .and. finidat == ' ') write(6,*) '   initial data created by model'
    if (nsrest == 0 .and. finidat /= ' ') write(6,*) '   initial data   = ',trim(finidat)
    if (nsrest /= 0) write(6,*) '   restart data   = ',trim(nrevsn)





    write(6,*) '   atmosperhic forcing data is from cam model'






    write(6,*) 'history and restart parameters:'
    if (hist_ndens == 1) then
       write(6,*)'   history tape data will be double precision'
    else if (hist_ndens == 2) then
       write(6,*)'   history tape data will be single precision'
    end if
    write(6,101) (i,hist_dov2xy(i),i=1,nhist)
    write(6,*) '   there will be ',nhist,' history files'
    do i = 1, nhist
       if (hist_nhtfrq(i)==0) then
          write(6,*) '   history file ',i,' is monthly averaged'
       else
          write(6,*) '   history file ',i,' time interval (iterations)= ', hist_nhtfrq(i)
       endif
    end do
    write(6,104) (i,hist_mfilt(i),i=1,nhist)
    if (mss_irt /= 0) then
       write(6,*)'   mass store path                    = ',trim(archive_dir)
       write(6,*)'   mass store retention (days)        = ',mss_irt
       write(6,*)'   mass store write password          = ',mss_wpass
    endif
    write(6,*)'   restart pointer file directory     = ',trim(rpntdir)
    write(6,*)'   restart pointer file name          = ',trim(rpntfil)
    if (hist_crtinic == 'MONTHLY') then
       write(6,*)'initial datasets will be written monthly'
    else if (hist_crtinic == 'YEARLY') then
       write(6,*)'initial datasets will be written yearly'
    else
       write(6,*)'initial datasets will not be produced'
    endif
    write(6,*) 'model physics parameters:'



    write(6,*) '   flag for random perturbation test is not set'

    write(6,*) '   energy and water conservation checks   = ',conchk
    write(6,*) '   solar radiation frequency (iterations) = ',irad
# 813 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/main/controlMod.F90"
    if (nsrest == 1) then
       write(6,*) 'restart warning:'
       write(6,*) '   Namelist not checked for agreement with initial run.'
       write(6,*) '   Namelist should not differ except for ending time step and run type'
    end if
    if (nsrest == 3) then
       write(6,*) 'branch warning:'
       write(6,*) '   Namelist not checked for agreement with initial run.'
       write(6,*) '   Surface data set and reference date should not differ from initial run'
    end if




101 format (1x,'   history fields are grid-average    = ',4(i1,':',l1,' '))
104 format (1x,'   time samples per history file      = ',4(i1,':',i2,' '))

  end subroutine control_print

!=======================================================================

end module controlMod

