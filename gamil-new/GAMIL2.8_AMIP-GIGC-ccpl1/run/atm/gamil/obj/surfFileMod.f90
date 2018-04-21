# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/main/surfFileMod.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/main/surfFileMod.F90"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/main/surfFileMod.F90" 2

# 1 "./preproc.h" 1






 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/main/surfFileMod.F90" 2

module surfFileMod

!=======================================================================
contains
!=======================================================================

  subroutine surfrd(veg, wt,  &
       cam_longxy, cam_latixy, cam_numlon, cam_landfrac, cam_landmask)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Read surface data file and make subgrid patches
! 
! Method: 
! The model's surface dataset recognizes 5 basic land cover types within 
! a grid cell: lake, wetland, urban, glacier, and vegetated. The vegetated 
! portion of the grid cell is comprised of up to [maxpatch_pft] PFTs 
! (patches). These subgrid patches are read in explicitly from the surface
! dataset at the model resolution for each grid cell. This is in contrast 
! to LSMv1, where the PFTs were built implicitly from biome types.
!
! The monthly surface dataset contains the following (on the 
! [lsmlon] x [lsmlat] grid). 
!    o real edges of grid
!    o real latitude  of grid cell (degrees)
!    o real longitude of grid cell (degrees)
!    o integer number of longitudes per latitude
!    o integer surface type: 0 = ocean or 1 = land
!    o integer soil color (1 to 9) for use with soil albedos
!    o real soil texture, %sand, for thermal and hydraulic properties
!    o real soil texture, %clay, for thermal and hydraulic properties
!    o real % of cell covered by lake    for use as subgrid patch
!    o real % of cell covered by wetland for use as subgrid patch
!    o real % of cell that is urban      for use as subgrid patch
!    o real % of cell that is glacier    for use as subgrid patch
!    o integer PFTs
!    o real % abundance PFTs (as a percent of vegetated area)
!
! OFFLINE MODE ONLY:
! Surface grid edges -- Grids do not have to be global. 
! If the model grid is read in from the surface dataset, the grid is 
! assumed to be global although it does not have to be regular (it
! can be a gaussian grid, for example). If the model grid is generated 
! by clm2, the grid does not have to be global but it but must be 
! regular and the user must define the north, east, south, and west edges
! as noted below.
!
!    o lsmedge(1) = northern edge of grid (degrees): >  -90 and <= 90
!    o lsmedge(2) = eastern edge of grid (degrees) : see following notes
!    o lsmedge(3) = southern edge of grid (degrees): >= -90 and <  90
!    o lsmedge(4) = western edge of grid (degrees) : see following notes
!
!    For partial grids, the northern and southern edges are any latitude
!    between 90 (North Pole) and -90 (South Pole). Western and eastern
!    edges are any longitude between -180 and 180, with longitudes
!    west of Greenwich negative. That is, western edge >= -180 and < 180;
!    eastern edge > western edge and <= 180.
!
!    For global grids, the northern and southern edges are 90 (North Pole)
!    and -90 (South Pole). The western and eastern edges depend on
!    whether the grid starts at Dateline or Greenwich. Regardless,
!    these edges must span 360 degrees. 
!
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id: surfFileMod.F90,v 1.13.2.5.6.1 2002/10/03 20:07:43 erik Exp $ 
!-----------------------------------------------------------------------

    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varpar                      !parameters
    use clm_varctl                      !control variables 
    use clm_varsur                      !surface data 
    use pft_varcon                      !vegetation type (PFT) 
    use fileutils, only : getfil
    use spmdMod                         !spmd routines and variables       
    use areaMod                         !area avg routines
    implicit none
    include 'netcdf.inc'

! ------------------------ arguments------------------------------------
    integer , intent(out) :: veg(lsmlon,lsmlat,maxpatch) !PFT 
    real(r8), intent(out) :: wt(lsmlon,lsmlat,maxpatch)  !subgrid weights
    real(r8), optional, intent(in) :: cam_longxy(:,:)    !cam lon values
    real(r8), optional, intent(in) :: cam_latixy(:,:)    !cam lat values 
    integer , optional, intent(in) :: cam_numlon(:)      !cam number of longitudes 
    real(r8), optional, intent(in) :: cam_landfrac(:,:)  !cam fractional land
    integer , optional, intent(in) :: cam_landmask(:,:)  !cam land mask
! ----------------------------------------------------------------------

! ------------------------ local variables -----------------------------
    character(len=256) :: locfn                    !local file name
    integer  :: i,j,k,m                            !indices
    integer  :: ier                                !error status 
    integer  :: pft(lsmlon,lsmlat,maxpatch_pft)    !PFT 
    integer  :: ncid,dimid,varid                   !input netCDF id's
    integer  :: beg4d(4),len4d(4)                  !netCDF variable edges
    integer  :: nlon_i                             !number of input data longitudes
    integer  :: nlat_i                             !number of input data latitudes
    integer  :: nlev_i                             !number of input data levels
    integer  :: maxpatch_pft_i                     !number of input data pft types
    real(r8) :: pctpft(lsmlon,lsmlat,maxpatch_pft) !percent of vegetated area for PFTs
    real(r8) :: sumscl                             !temporory scalar sum 
    real(r8) :: sumvec(lsmlon,lsmlat)              !temporary vector sum 
! ----------------------------------------------------------------------

    if (masterproc) then

! Initialize surface data to unused value

       landmask(:,:) = -999
       landfrac(:,:) = -999.
       latixy(:,:)   = -999.
       longxy(:,:)   = -999.
       soic2d(:,:)   = -999
       sand3d(:,:,:) = -999.
       clay3d(:,:,:) = -999.
       pctlak(:,:)   = -999.
       pctwet(:,:)   = -999.
       pcturb(:,:)   = -999.
       pctgla(:,:)   = -999.

! Obtain netcdf file and read surface data

       write (6,*) 'Attempting to read surface boundary data .....'

       call getfil (fsurdat, locfn, 0)
       call wrap_open(locfn, 0, ncid)

       call wrap_inq_dimid  (ncid, 'lsmlon', dimid)
       call wrap_inq_dimlen (ncid, dimid, nlon_i)
       if (nlon_i /= lsmlon) then
          write(6,*)'SURFRD: parameter lsmlon= ',lsmlon, &
               'does not equal input nlon_i= ',nlon_i
          call endrun
       endif

       call wrap_inq_dimid  (ncid, 'lsmlat', dimid)
       call wrap_inq_dimlen (ncid, dimid, nlat_i)
       if (nlat_i /= lsmlat) then
          write(6,*)'SURFRD: parameter lsmlat= ',lsmlat, &
               'does not equal input nlat_i= ',nlat_i
          call endrun
       endif

       call wrap_inq_dimid  (ncid, 'nlevsoi', dimid)
       call wrap_inq_dimlen (ncid, dimid, nlev_i)
       if (nlev_i /= nlevsoi) then
          write(6,*)'SURFRD: parameter nlevsoi= ',nlevsoi, &
               'does not equal input nlev_i= ',nlev_i
          call endrun
       endif

       call wrap_inq_dimid  (ncid, 'lsmpft', dimid)
       call wrap_inq_dimlen (ncid, dimid, maxpatch_pft_i)
       if (maxpatch_pft_i /= maxpatch_pft) then
          write(6,*)'SURFRD: parameter maxpatch_pft',maxpatch_pft, &
               'does not equal input maxpatch_pft_i= ',maxpatch_pft_i
          call endrun
       endif

# 188 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/main/surfFileMod.F90"

       call wrap_inq_varid (ncid, 'NUMLON', varid)
       call wrap_get_var_int   (ncid, varid  , numlon)

       fullgrid = .true.
       do j = 1,lsmlat
          if (numlon(j) < lsmlon) fullgrid = .false.
       end do

       if (fullgrid) then
          call wrap_inq_varid (ncid, 'LONGXY' , varid)
       else
          call wrap_inq_varid (ncid, 'RLONGXY', varid)
       endif
       call wrap_get_var_realx (ncid, varid, longxy)

       call wrap_inq_varid (ncid, 'LATIXY', varid)
       call wrap_get_var_realx (ncid, varid, latixy)

       call wrap_inq_varid (ncid, 'LANDMASK', varid) 
       call wrap_get_var_int (ncid, varid, landmask)

       call wrap_inq_varid (ncid, 'LANDFRAC', varid) 
       call wrap_get_var_realx (ncid, varid, landfrac)

       call wrap_inq_varid (ncid, 'SOIL_COLOR', varid)
       call wrap_get_var_int   (ncid, varid, soic2d)

       call wrap_inq_varid (ncid, 'PCT_SAND', varid)
       call wrap_get_var_realx (ncid, varid, sand3d)

       call wrap_inq_varid (ncid, 'PCT_CLAY', varid)
       call wrap_get_var_realx (ncid, varid, clay3d) 

       call wrap_inq_varid (ncid, 'PCT_WETLAND', varid)
       call wrap_get_var_realx (ncid, varid, pctwet)

       call wrap_inq_varid (ncid, 'PCT_LAKE', varid)
       call wrap_get_var_realx (ncid, varid, pctlak)

       call wrap_inq_varid (ncid, 'PCT_GLACIER', varid)
       call wrap_get_var_realx (ncid, varid, pctgla)

       call wrap_inq_varid (ncid, 'PCT_URBAN', varid)
       call wrap_get_var_realx (ncid, varid, pcturb)

       call wrap_inq_varid (ncid, 'PFT', varid)
       call wrap_get_var_int (ncid, varid, pft)

       call wrap_inq_varid (ncid, 'PCT_PFT', varid)
       call wrap_get_var_realx (ncid, varid, pctpft)

       call wrap_close(ncid)

! Define edges and area of grid cells

# 263 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/main/surfFileMod.F90"

      ! CLM2 is run as part of CAM model or through the CCSM flux coupler
      ! The input grid, land mask and fractional land must match the corresponding
      ! CAM values. Consequently, the following check must be done on input 
      ! latitudes and longitudes since these are computed each time by the CAM  
      ! model and might be different to roundoff depending on the platform

       do j = 1,lsmlat
          do i = 1,numlon(j)
             if ( abs(cam_latixy(i,j)-latixy(i,j)) > 1.e-6 ) then
                write(6,*)'CAM latitude ',cam_latixy(i,j),' and clm2 input latitude ', &
                     latixy(i,j),' has difference too large at i,j= ',i,j
!                call endrun
             else
                latixy(i,j) = cam_latixy(i,j)
             endif
             if ( abs(cam_longxy(i,j)-longxy(i,j)) > 1.e-6 ) then
                write(6,*)'CAM longitude ',cam_longxy(i,j),' and clm2 input longitude ', &
                     longxy(i,j),' has difference too large at i,j= ',i,j
                call endrun
             else
                longxy(i,j) = cam_longxy(i,j)
             endif
             if (cam_numlon(j) /= numlon(j)) then                 
                write(6,*)'CAM numlon array different from input CLM2 value' 
                write(6,*)'lat index= ',j,' cam numlon= ',cam_numlon(j), &
                     ' clm2 numlon= ',numlon(j)
                call endrun
             else if (cam_landmask(i,j) /= landmask(i,j)) then         
                write(6,*)'CAM land mask different from input CLM2 value' 
                write(6,*)'lat index= ',j,' lon index= ',i,&
                     ' cam landmask= ',cam_landmask(i,j),' clm2 landmask= ',landmask(i,j)
                call endrun
             elseif (cam_landfrac(i,j) /= landfrac(i,j)) then
                write(6,*)'CAM fractional land different from CLM2 value' 
                write(6,*)'lat index= ',j,' lon index= ',i,&
                     ' cam landmask= ',cam_landfrac(i,j),' clm2 landfrac= ',landfrac(i,j)
                call endrun
             endif
          end do
       end do

       call celledge (lsmlat, lsmlon, numlon, longxy, latixy, &
                      lats  , lonw  )

       call cellarea (lsmlat, lsmlon, numlon, lats, lonw, &
                      area   )
  


! Error check: valid PFTs and sum of cover must equal 100

       sumvec(:,:) = abs(sum(pctpft,dim=3)-100.)
       do j=1,lsmlat
          do i=1,numlon(j)
             do m = 1, maxpatch_pft
                if (pft(i,j,m)<0 .or. pft(i,j,m)>numpft) then
                   write(6,*)'SURFRD error: invalid PFT for i,j,m=',i,j,m,pft(i,j,m)
                   call endrun
                end if
             end do
             if (sumvec(i,j)>1.e-04 .and. landmask(i,j)==1) then
                write(6,*)'SURFRD error: PFT cover ne 100 for i,j=',i,j
                do m=1,maxpatch_pft
                   write(6,*)'m= ',m,' pft= ',pft(i,j,m)
                end do
                write(6,*)'sumvec= ',sumvec(i,j)
                call endrun
             end if
          end do
       end do

! Error check: percent glacier, lake, wetland, urban sum must be less than 100

       do j=1,lsmlat
          do i=1,numlon(j)
             sumscl = pctlak(i,j)+pctwet(i,j)+pcturb(i,j)+pctgla(i,j)
             if (sumscl > 100.+1.e-04) then
                write(6,*)'SURFRD error: PFT cover>100 for i,j=',i,j
                call endrun
             end if
          end do
       end do

! Check that urban parameterization is not yet implemented

       do j=1,lsmlat
          do i=1,numlon(j)
             if (pcturb(i,j) /= 0.) then
                write (6,*) 'urban parameterization not yet implemented'
                call endrun
             end if
          end do
       end do

    endif                     !end of if-masterproc block








    call mpi_bcast (numlon  , size(numlon)  , mpiint, 0, mpicom, ier)
    call mpi_bcast (latixy  , size(latixy)  , mpir8 , 0, mpicom, ier)
    call mpi_bcast (longxy  , size(longxy)  , mpir8 , 0, mpicom, ier)
    call mpi_bcast (landmask, size(landmask), mpiint, 0, mpicom, ier)
    call mpi_bcast (landfrac, size(landfrac), mpir8 , 0, mpicom, ier)
    call mpi_bcast (soic2d  , size(soic2d)  , mpiint, 0, mpicom, ier)
    call mpi_bcast (sand3d  , size(sand3d)  , mpir8 , 0, mpicom, ier)
    call mpi_bcast (clay3d  , size(clay3d)  , mpir8 , 0, mpicom, ier)
    call mpi_bcast (pctwet  , size(pctwet)  , mpir8 , 0, mpicom, ier)
    call mpi_bcast (pctlak  , size(pctlak)  , mpir8 , 0, mpicom, ier)
    call mpi_bcast (pctgla  , size(pctgla)  , mpir8 , 0, mpicom, ier)
    call mpi_bcast (pcturb  , size(pcturb)  , mpir8 , 0, mpicom, ier)
    call mpi_bcast (pft     , size(pft)     , mpiint, 0, mpicom, ier)
    call mpi_bcast (pctpft  , size(pctpft)  , mpir8 , 0, mpicom, ier)


! Make patch arrays, [veg] and [wt]:
! [veg] sets the PFT for each of the [maxpatch] patches on the 2d model grid.
! [wt]  sets the relative abundance of the PFT on the 2d model grid.
! Fill in PFTs for vegetated portion of grid cell. Fractional areas for
! these points [pctpft] pertain to "vegetated" area not to total grid area.
! So need to adjust them for fraction of grid that is vegetated.
! Next, fill in urban, lake, wetland, and glacier patches.

    veg(:,:,:) = 0
    wt(:,:,:)  = 0.
    do j=1,lsmlat
       do i=1,numlon(j)
          if (landmask(i,j) == 1) then
             sumscl = pcturb(i,j)+pctlak(i,j)+pctwet(i,j)+pctgla(i,j)
             do m = 1, maxpatch_pft
                veg(i,j,m) = pft(i,j,m)
                wt(i,j,m) = pctpft(i,j,m) * (100.-sumscl)/10000.
             end do
             veg(i,j,npatch_urban) = noveg
             wt(i,j,npatch_urban) = pcturb(i,j)/100.
             veg(i,j,npatch_lake) = noveg
             wt(i,j,npatch_lake) = pctlak(i,j)/100.
             veg(i,j,npatch_wet) = noveg
             wt(i,j,npatch_wet) = pctwet(i,j)/100.
             veg(i,j,npatch_gla) = noveg
             wt(i,j,npatch_gla) = pctgla(i,j)/100.
          end if
       end do
    end do

    sumvec(:,:) = abs(sum(wt,dim=3)-1.)
    do j=1,lsmlat
       do i=1,numlon(j)
          if (sumvec(i,j) > 1.e-06 .and. landmask(i,j)==1) then
             write (6,*) 'SURFRD error: WT > 1 occurs at i,j= ',i,j 
             call endrun
          endif
       end do
    end do

    if ( masterproc )then
       write (6,*) 'Successfully read surface boundary data'
       write (6,*)
    end if

    return
  end subroutine surfrd

!=======================================================================

  subroutine surfwrt(fname, pft, pctpft, mlai, msai, mhgtt, mhgtb)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Write surface data file
!
! Method: 
! 
! Author: Mariana Vertenstein
!
! -----------------------------------------------------------------

    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varpar           
    use clm_varsur           
    use clm_varctl
    use fileutils, only : get_filename
    implicit none
    include 'netcdf.inc'

! ------------------------ arguments ------------------------------
    character(len=*), intent(in) :: fname                        !filename to create
    integer , intent(in) :: pft(lsmlon,lsmlat,maxpatch_pft)      !vegetation type
    real(r8), intent(in) :: pctpft(lsmlon,lsmlat,maxpatch_pft)   !vegetation type subgrid weights
    real(r8), intent(in) :: mlai (lsmlon,lsmlat,maxpatch_pft,12) !monthly lai
    real(r8), intent(in) :: msai (lsmlon,lsmlat,maxpatch_pft,12) !monthly sai
    real(r8), intent(in) :: mhgtt(lsmlon,lsmlat,maxpatch_pft,12) !monthly hgt at top
    real(r8), intent(in) :: mhgtb(lsmlon,lsmlat,maxpatch_pft,12) !monthly hgt at bottom
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------
    integer i,m                  !indices
    integer ncid                 !netcdf id
    integer omode                !netcdf output mode
    integer ret                  !netcdf return status

    integer dimtim_id            !id for time dimension
    integer dimlon_id            !id for grid longitude 
    integer dimlat_id            !id for grid latitude 
    integer dimlev_id            !id for soil layer dimension
    integer dimpft_id            !id for plft
    integer dimstr_id            !id for character string variables







    integer longxy_id            !variable id
    integer latixy_id            !variable id
    integer numlon_id            !variable id
    integer landmask_id          !variable id
    integer landfrac_id          !variable id
    integer soic2d_id            !variable id
    integer sand3d_id            !variable id
    integer clay3d_id            !variable id
    integer pctlak_id            !variable id
    integer pctwet_id            !variable id
    integer pctgla_id            !variable id
    integer pcturb_id            !variable id
    integer pft_id               !variable id
    integer pctpft_id            !variable id
    integer mlai_id              !variable id
    integer msai_id              !variable id
    integer mhgtt_id             !variable id 
    integer mhgtb_id             !variable id

    integer dim1_id(1)           !dim id for 1-d variables
    integer dim2_id(2)           !dim id for 2-d variables
    integer dim3_id(3)           !dim id for 3-d variables
    integer dim4_id(4)           !dim id for 3-d variables
    integer beg4d(4),len4d(4)    !netCDF variable edges

    character(len=256) str       !global attribute string 
    character(len=256) name      !name of attribute
    character(len=256) unit      !units of attribute

    integer values(8)
    character(len=18) datetime
    character(len= 8) date
    character(len=10) time
    character(len= 5) zone
! -----------------------------------------------------------------

! Create new netCDF file. File will be in define mode
! Set fill mode to "no fill" to optimize performance

    call wrap_create (trim(fname), nf_clobber, ncid)

    ret = nf_set_fill (ncid, nf_nofill, omode)
    if (ret .ne. nf_noerr) then
       write (6,*) ' netCDF error = ',nf_strerror(ret)
       call endrun
    end if

! Create global attributes. Attributes are used to store information
! about the data set. Global attributes are information about the
! data set as a whole, as opposed to a single variable

    str = 'NCAR-CSM'
    call wrap_put_att_text (ncid, NF_GLOBAL, 'Conventions', trim(str))

    call date_and_time (date, time, zone, values)
    datetime(1:8) =        date(5:6) // '/' // date(7:8) // '/' // date(3:4)
    datetime(9:)  = ' ' // time(1:2) // ':' // time(3:4) // ':' // time(5:6) // ' '
    str = 'created on: ' // datetime 
    call wrap_put_att_text (ncid, NF_GLOBAL, 'History', trim(str))

!!  call getenv ('LOGNAME', str)  !!(2003.10.22)
    call getenv ('HOME'   , str)
    call wrap_put_att_text (ncid, NF_GLOBAL, 'Logname',trim(str))

    call getenv ('HOST', str)
    call wrap_put_att_text (ncid, NF_GLOBAL, 'Host', trim(str))

    str = 'Community Land Model: CLM2'
    call wrap_put_att_text (ncid, NF_GLOBAL, 'Source', trim(str))

    str = '$Name: GAMIL1.0 $' 
    call wrap_put_att_text (ncid, NF_GLOBAL, 'Version', trim(str))

    str = '$Id: surfFileMod.F90,v 1.13.2.5.6.1 2002/10/03 20:07:43 erik Exp $'
    call wrap_put_att_text (ncid, NF_GLOBAL, 'Revision_Id', trim(str))

    if (offline_rdgrid) then
       str = mksrf_offline_fgrid
       call wrap_put_att_text(ncid, NF_GLOBAL, 'Input_grid_dataset', trim(str))
    else
       str = mksrf_offline_fnavyoro
       call wrap_put_att_text(ncid, NF_GLOBAL, 'Input_navy_oro_dataset', trim(str))
    endif

    str = get_filename(mksrf_fvegtyp)
    call wrap_put_att_text(ncid, NF_GLOBAL, 'Vegetation_type_raw_data_filename', trim(str))
         
    str = get_filename(mksrf_fsoitex)
    call wrap_put_att_text(ncid, NF_GLOBAL, 'Soil_texture_raw_data_file_name', trim(str))

    str = get_filename(mksrf_fsoicol)
    call wrap_put_att_text(ncid, NF_GLOBAL, 'Soil_color_raw_data_file_name', trim(str))

    str = get_filename(mksrf_flanwat)
    call wrap_put_att_text(ncid, NF_GLOBAL, 'Inland_water_raw_data_file_name', trim(str))

    str = get_filename(mksrf_fglacier)
    call wrap_put_att_text(ncid, NF_GLOBAL, 'Glacier_raw_data_file_name', trim(str))

    str = get_filename(mksrf_furban)
    call wrap_put_att_text(ncid, NF_GLOBAL, 'Urban_raw_data_file_name', trim(str))

    str = get_filename(mksrf_flai)
    call wrap_put_att_text(ncid, NF_GLOBAL, 'Lai_raw_data_file_name', trim(str))




    str = 'run through cam'



    call wrap_put_att_text(ncid, NF_GLOBAL, 'Run_mode', trim(str))

! Define dimensions. Array dimensions are referenced by an
! associated dimenision id: e.g., lsmlon_id -> lsmlon.
! o Time is an unlimited dimension.
! o Character string is treated as an array of characters. 

    call wrap_def_dim (ncid, 'lsmlon' , lsmlon      , dimlon_id)
    call wrap_def_dim (ncid, 'lsmlat' , lsmlat      , dimlat_id)
    call wrap_def_dim (ncid, 'nlevsoi', nlevsoi     , dimlev_id)
    call wrap_def_dim (ncid, 'lsmpft' , maxpatch_pft, dimpft_id)
    call wrap_def_dim (ncid, 'time'   , nf_unlimited, dimtim_id)
    call wrap_def_dim (ncid, 'nchar'  , 128         , dimstr_id) 

! Define time-independent variables and their attributes. 
! Variables are referenced by an associated variable id, e.g., 
! the netCDF variable 'LONGXY' is referenced by the id longxy_id

# 640 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/main/surfFileMod.F90"

    dim1_id(1) = dimlat_id

    name = 'number of longitudes for each latitude'
    unit = 'unitless'
    call wrap_def_var (ncid, 'NUMLON', nf_int, 1, dim1_id, numlon_id)
    call wrap_put_att_text (ncid, numlon_id, 'long_name', name)
    call wrap_put_att_text (ncid, numlon_id, 'units'    , unit)

! Model grids MUST BE double precision so that when running through 
! the cam mode the latitude and longitude grid are bfb the same 
! as the cam grid

    dim2_id(1) = dimlon_id
    dim2_id(2) = dimlat_id

    if (fullgrid) then
       name = 'longitude'
    else
       name = 'rlongitude'
    endif
    unit = 'degrees east'
    if (fullgrid) then
       call wrap_def_var (ncid, 'LONGXY' , nf_double, 2, dim2_id, longxy_id)
    else
       call wrap_def_var (ncid, 'RLONGXY', nf_double, 2, dim2_id, longxy_id)
    endif
    call wrap_put_att_text (ncid, longxy_id, 'long_name', name)
    call wrap_put_att_text (ncid, longxy_id, 'units'    , unit)

    name = 'latitude'
    unit = 'degrees north'
    call wrap_def_var (ncid, 'LATIXY', nf_double, 2, dim2_id, latixy_id)
    call wrap_put_att_text (ncid, latixy_id, 'long_name', name)
    call wrap_put_att_text (ncid, latixy_id, 'units'    , unit)

! Land mask and Land fraction

    name = 'land/ocean mask'
    unit = '0=ocean and 1=land'
    call wrap_def_var (ncid, 'LANDMASK', nf_int, 2, dim2_id, landmask_id)
    call wrap_put_att_text (ncid, landmask_id, 'long_name', name)
    call wrap_put_att_text (ncid, landmask_id, 'units'    , unit)

    name = 'land fraction'
    unit = 'unitless'
    call wrap_def_var (ncid, 'LANDFRAC', nf_double, 2, dim2_id, landfrac_id)
    call wrap_put_att_text (ncid, landfrac_id, 'long_name', name)
    call wrap_put_att_text (ncid, landfrac_id, 'units'    , unit)

! Surface variables

    name = 'soil color'
    unit = 'unitless'
    call wrap_def_var (ncid, 'SOIL_COLOR', nf_int, 2, dim2_id, soic2d_id)
    call wrap_put_att_text (ncid, soic2d_id, 'long_name', name)
    call wrap_put_att_text (ncid, soic2d_id, 'units'    , unit)

    dim3_id(1) = dimlon_id
    dim3_id(2) = dimlat_id
    dim3_id(3) = dimlev_id

    name = 'percent sand'
    unit = 'unitless'
    call wrap_def_var (ncid ,'PCT_SAND' ,nf_float, 3, dim3_id, sand3d_id)
    call wrap_put_att_text (ncid, sand3d_id, 'long_name', name)
    call wrap_put_att_text (ncid, sand3d_id, 'units'    , unit)

    name = 'percent clay'
    unit = 'unitless'
    call wrap_def_var (ncid ,'PCT_CLAY' ,nf_float, 3, dim3_id, clay3d_id)
    call wrap_put_att_text (ncid, clay3d_id, 'long_name', name)
    call wrap_put_att_text (ncid, clay3d_id, 'units'    , unit)

    name = 'percent wetland'
    unit = 'unitless'
    call wrap_def_var (ncid ,'PCT_WETLAND' ,nf_float, 2, dim2_id, pctwet_id)
    call wrap_put_att_text (ncid, pctwet_id, 'long_name', name)
    call wrap_put_att_text (ncid, pctwet_id, 'units'    , unit)

    name = 'percent lake'
    unit = 'unitless'
    call wrap_def_var (ncid ,'PCT_LAKE' ,nf_float, 2, dim2_id, pctlak_id)
    call wrap_put_att_text (ncid, pctlak_id, 'long_name', name)
    call wrap_put_att_text (ncid, pctlak_id, 'units'    , unit)

    name = 'percent glacier'
    unit = 'unitless'
    call wrap_def_var (ncid ,'PCT_GLACIER' ,nf_float, 2, dim2_id, pctgla_id)
    call wrap_put_att_text (ncid, pctgla_id, 'long_name', name)
    call wrap_put_att_text (ncid, pctgla_id, 'units'    , unit)

    name = 'percent urban'
    unit = 'unitless'
    call wrap_def_var (ncid ,'PCT_URBAN' ,nf_float, 2, dim2_id, pcturb_id)
    call wrap_put_att_text (ncid, pcturb_id, 'long_name', name)
    call wrap_put_att_text (ncid, pcturb_id, 'units'    , unit)

    dim3_id(1) = dimlon_id
    dim3_id(2) = dimlat_id
    dim3_id(3) = dimpft_id

    name = 'plant functional type' 
    unit = 'unitless'
    call wrap_def_var (ncid ,'PFT' ,nf_int, 3, dim3_id, pft_id)
    call wrap_put_att_text (ncid, pft_id, 'long_name', name)
    call wrap_put_att_text (ncid, pft_id, 'units'    , unit)

    name = 'percent plant functional type'
    unit = 'unitless'
    call wrap_def_var (ncid ,'PCT_PFT' ,nf_float, 3, dim3_id, pctpft_id)
    call wrap_put_att_text (ncid, pctpft_id, 'long_name', name)
    call wrap_put_att_text (ncid, pctpft_id, 'units'    , unit)

! LAI/SAI/HEIGHT data

    dim4_id(1) = dimlon_id
    dim4_id(2) = dimlat_id
    dim4_id(3) = dimpft_id
    dim4_id(4) = dimtim_id

    name = 'monthly leaf area index'
    unit = 'unitless'
    call wrap_def_var (ncid ,'MONTHLY_LAI', nf_float, 4, dim4_id, mlai_id)
    call wrap_put_att_text (ncid, mlai_id, 'long_name', name)
    call wrap_put_att_text (ncid, mlai_id, 'units'    , unit)

    name = 'monthly stem area index'
    unit = 'unitless'
    call wrap_def_var (ncid ,'MONTHLY_SAI', nf_float, 4, dim4_id, msai_id)
    call wrap_put_att_text (ncid, msai_id, 'long_name', name)
    call wrap_put_att_text (ncid, msai_id, 'units'    , unit)

    name = 'monthly height top'
    unit = 'meters'
    call wrap_def_var (ncid ,'MONTHLY_HEIGHT_TOP', nf_float, 4, dim4_id, mhgtt_id)
    call wrap_put_att_text (ncid, mhgtt_id, 'long_name', name)
    call wrap_put_att_text (ncid, mhgtt_id, 'units'    , unit)

    name = 'monthly height bottom'
    unit = 'meters'
    call wrap_def_var (ncid ,'MONTHLY_HEIGHT_BOT', nf_float, 4, dim4_id, mhgtb_id)
    call wrap_put_att_text (ncid, mhgtb_id, 'long_name', name)
    call wrap_put_att_text (ncid, mhgtb_id, 'units'    , unit)

! Finish creating netcdf file

    ret = nf_enddef(ncid)
    if (ret /= 0) then
       write (6,*)'failed to end define mode'
       call endrun
    end if

! Write out data

# 803 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/main/surfFileMod.F90"
    call wrap_put_var_int   (ncid, numlon_id  , numlon    )
    call wrap_put_var_realx (ncid, longxy_id  , longxy    )
    call wrap_put_var_realx (ncid, latixy_id  , latixy    )
    call wrap_put_var_int   (ncid, landmask_id, landmask  )
    call wrap_put_var_realx (ncid, landfrac_id, landfrac  )
    call wrap_put_var_int   (ncid, soic2d_id  , soic2d    )
    call wrap_put_var_realx (ncid, sand3d_id  , sand3d    )
    call wrap_put_var_realx (ncid, clay3d_id  , clay3d    ) 
    call wrap_put_var_realx (ncid, pctwet_id  , pctwet    )
    call wrap_put_var_realx (ncid, pctlak_id  , pctlak    )
    call wrap_put_var_realx (ncid, pctgla_id  , pctgla    )
    call wrap_put_var_realx (ncid, pcturb_id  , pcturb    )
    call wrap_put_var_int   (ncid, pft_id     , pft       )
    call wrap_put_var_realx (ncid, pctpft_id  , pctpft    )

    do m=1,12
       beg4d(1) = 1     
       len4d(1) = lsmlon
       beg4d(2) = 1     
       len4d(2) = lsmlat
       beg4d(3) = 1
       len4d(3) = maxpatch_pft
       beg4d(4) = m
       len4d(4) = 1
       call wrap_put_vara_realx (ncid, mlai_id , beg4d, len4d, mlai(1,1,1,m))
       call wrap_put_vara_realx (ncid, msai_id , beg4d, len4d, msai(1,1,1,m))
       call wrap_put_vara_realx (ncid, mhgtt_id, beg4d, len4d, mhgtt(1,1,1,m))
       call wrap_put_vara_realx (ncid, mhgtb_id, beg4d, len4d, mhgtb(1,1,1,m))
    end do

! Close output file

    call wrap_close(ncid)

    return
  end subroutine surfwrt

!=======================================================================

end module surfFileMod
