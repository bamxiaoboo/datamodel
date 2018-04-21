# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/main/clm_map.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/main/clm_map.F90"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/main/clm_map.F90" 2

# 1 "./preproc.h" 1






 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/main/clm_map.F90" 2

subroutine clm_map (wtxy) 

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Build 1d subgrid patch <-> 2d grid mapping indices and weights
! 
! Method: 
! Build mapping indices and weights: [lsmlon]x[lsmlat] 2d grid <->
! [numland] vector of land points <-> [numpatch] vector of subgrid patches. 
! Allow for variable longitudinal resolution: [numlon] <= [lsmlon]
! The land surface model works by gathering all the land points on a
! [lon]x[lat] grid into a vector of [numland] land points. This
! is then expanded into a vector of [numpatch] subgrid patches, allowing
! for up to [maxpatch] subgrid patches per land point. 
! [ixy], [jxy], [patch], and [land] are indices for the mapping: 
! [lsmlon]x[lsmlat] grid <-> [numland] vector of land points <-> 
! [numpatch] vector of subgrid points. [landvec%wtxy] are the weights 
! to obtain the grid average from the subgrid patches.
!
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id: clm_map.F90,v 1.2.6.6.6.1 2002/10/03 20:07:34 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varpar , only : lsmlon, lsmlat, maxpatch
  use clm_varsur , only : numlon, landmask
  use clm_varmap , only : numland, numpatch, begpatch, endpatch, &
                          begland, endland, landvec, patchvec, mapvar_ini
  use clm_varder , only : clm_varder_ini  
  use histFileMod, only : histvar_ini
  use mvegFileMod, only : monthveg_ini   

  use spmdMod    , only : masterproc, iam, npes, proc_landi, proc_landf, proc_landpts, &
                          proc_patchi, proc_patchf, proc_patchpts, spmd_init_patch  



  implicit none

! ------------------------ arguments----------------------------------
  real(r8), intent(in) :: wtxy(lsmlon,lsmlat,maxpatch) !subgrid patch weights
! --------------------------------------------------------------------

! ------------------------ local variables ---------------------------
  integer  :: i,j,l,m,n,p                 !indices
  integer  :: li,lf                       !land indices
  integer  :: nland                       !number of land points
  integer  :: npatch                      !number of patch points 
  real(r8) :: sumwts                      !sum of wtxy over patches
  real(r8) :: procwt                      !weight of processor patches
  real(r8), allocatable :: land_wt(:)     !weight of patches for land point
  integer , allocatable :: land_patchs(:) !numger of patches for land point
! --------------------------------------------------------------------



! ----------------------------------------------------------------------
! Initialize spmd arrays for number of land/patch points per processor
! ----------------------------------------------------------------------

  call spmd_init_patch



! --------------------------------------------------------------------
! Find total number of land grid cells [numland] and total number of 
! patches [numpatch] allowing for multiple subgrid patches in a grid cell.
! --------------------------------------------------------------------

  npatch = 0
  nland = 0
  do j = 1, lsmlat
     do i = 1, numlon(j)
        if (landmask(i,j) == 1) then         
           nland = nland+1                               !land point number
           do m = 1, maxpatch                
              if (wtxy(i,j,m) > 0.)  npatch = npatch+1   !subgrid patch number
           end do
        end if
     end do
  end do
  numland = nland
  numpatch = npatch

  if (masterproc) then
     write (6,*)' Surface Grid Characteristics'
     write (6,*)'   longitude points         = ',lsmlon
     write (6,*)'   latitude points          = ',lsmlat
     write (6,*)'   minimum longitude points = ',minval(numlon)
     write (6,*)'   maximum longitude points = ',maxval(numlon)
     write (6,*)'   land points on grid                    = ',numland
     write (6,*)'   total points including subgrid patches = ',numpatch
     write (6,*)
  endif
       
! --------------------------------------------------------------------
! Find beginning and ending patch index for each process
! --------------------------------------------------------------------

# 114 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/main/clm_map.F90"

  if (npes > numpatch) then
     write(6,*)'CLM_MAP: asking for more processors than patches'
     call endrun
  endif

  allocate (land_wt(numland)) 
  allocate (land_patchs(numland))

  procwt = 1./npes 

  nland  = 0
  sumwts = 0.
  do j = 1, lsmlat
     do i = 1, numlon(j)
        if (landmask(i,j) == 1) then                 
           nland  = nland+1
           npatch = 0
           do m = 1, maxpatch                           
              if (wtxy(i,j,m) > 0.) npatch = npatch+1
           end do
           land_wt(nland) = real(npatch)/real(numpatch)
           sumwts = sumwts + land_wt(nland)
           land_patchs(nland) = npatch
        end if
     end do
  end do
  write(6,*)'SUM: iam= ',iam,' sumwts= ',sumwts

  nland = 0
  npatch = 0
  do p = 0, npes-1
     proc_landi(p) = nland + 1
     proc_landf(p) = nland
     proc_patchi(p) = npatch + 1 
     proc_patchf(p) = npatch
     sumwts = 0.
     do l = proc_landi(p), numland
        sumwts = sumwts + land_wt(l)
        npatch = npatch + land_patchs(l)
        if ((sumwts>=procwt .and. (p/=npes-1)) .or. (npatch==numpatch)) then
           nland = l
           proc_landf(p) = nland
           proc_patchf(p) = npatch
           EXIT
        elseif (npatch > numpatch) then
           write(6,*)'CLM_MAP: iam= ',iam,' p= ',p,' max patches exceeded'
           write(6,*)'npatch = ',npatch,' numpatch= ',numpatch
           call endrun
        endif
     end do
     proc_patchpts(p) = proc_patchf(p) - proc_patchi(p) + 1 
     proc_landpts(p) = proc_landf(p) - proc_landi(p) + 1
  end do
  begpatch = proc_patchi(iam)
  endpatch = proc_patchf(iam) 
  begland  = proc_landi(iam)
  endland  = proc_landf(iam)

  write(6,*)'iam= ',iam,' begland = ',proc_landi(iam) ,' endland = ',proc_landf(iam),&
       ' total land points = ',proc_landpts(iam)
  write(6,*)'iam= ',iam,' begpatch= ',proc_patchi(iam),' endpatch= ',endpatch, &
       ' total patch pionts= ',proc_patchpts(iam)

  deallocate (land_wt)
  deallocate (land_patchs)



! --------------------------------------------------------------------
! Allocate dynamic memory
! --------------------------------------------------------------------

  call mapvar_ini
  call clm_varder_ini
  call histvar_ini
  call monthveg_ini

! --------------------------------------------------------------------
! Build 1d land vector and 1d patch vector mapping components
! --------------------------------------------------------------------

! Determine land vector and patch vector mapping components

  landvec%wtxy(:,:) = 0._r8
  patchvec%wtxy(:)  = 0._r8
  npatch = 0
  nland  = 0
  do j = 1, lsmlat
     do i = 1, numlon(j)
        if (landmask(i,j) == 1) then                 
           nland = nland+1                           
           landvec%ixy(nland) = i                       !land longitude index
           landvec%jxy(nland) = j                       !land latitude index
           do m = 1, maxpatch                           
              if (wtxy(i,j,m) > 0.) then                
                 npatch = npatch+1                      
                 landvec%patch(nland,m) = npatch        !land subgrid patch number
                 landvec%wtxy(nland,m)  = wtxy(i,j,m)   !land subgrid weights
                 patchvec%ixy(npatch)   = i             !patch longitude index
                 patchvec%jxy(npatch)   = j             !patch latitude index
                 patchvec%mxy(npatch)   = m             !patch subgrid index of lnd point
                 patchvec%wtxy(npatch)  = wtxy(i,j,m)   !patch weight
                 patchvec%land(npatch)  = nland         !patch land index
              end if
           end do
        end if
     end do
  end do

! Initialize land vector patches with zero weights to first patch with 
! non-zero weight

  do l = 1, numland
     n = 0
     do m = maxpatch, 1, -1
        if (landvec%wtxy(l,m) > 0.) n = m
     end do
     if (n == 0) then
        write (6,*) 'CLM_MAP error: n = 0'
        call endrun
     end if
     do m = 1, maxpatch
        if (landvec%wtxy(l,m) == 0.) then
           landvec%patch(l,m) = landvec%patch(l,n)
        endif
     end do
  end do

! Error check: make sure weights sum to one for each land cell

  do l = 1, numland
     i = landvec%ixy(l)
     j = landvec%jxy(l)
     sumwts = 0
     do m = 1, maxpatch
        sumwts = sumwts + landvec%wtxy(l,m)
     end do
     if (abs(sumwts-1.) > 1.e-05) then
        write (6,*) 'CLM_MAP error: weights do not sum to 1'
        write (6,*) 'lon = ',i,' lat = ',j,' :sum = ',sumwts
        call endrun
     end if
  end do

  return
end subroutine clm_map

subroutine clm_map_rebuild_from_atm (wtxy) 

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varpar , only : lsmlon, lsmlat, maxpatch
  use clm_varsur , only : numlon, landmask
  use clm_varmap , only : numland, numpatch, begpatch, endpatch, &
                          begland, endland, landvec, patchvec, mapvar_ini
  use clm_varder , only : clm_varder_ini
  use histFileMod, only : histvar_ini
  use mvegFileMod, only : monthveg_ini   

  use phys_grid
  use spmdMod    , only : masterproc, iam, npes, proc_landi, proc_landf, proc_landpts, &
                          proc_patchi, proc_patchf, proc_patchpts, spmd_init_patch  



  implicit none

! ------------------------ arguments----------------------------------
  real(r8), intent(in) :: wtxy(lsmlon,lsmlat,maxpatch) !subgrid patch weights
! --------------------------------------------------------------------

! ------------------------ local variables ---------------------------
  integer  :: i,j,k,l,m,n,p               !indices
  integer  :: gchnk, ncols,nlchunks
  integer  :: nland                       !number of land points
  integer  :: npatch                      !number of patch points 
  real(r8) :: sumwts                      !sum of wtxy over patches
  real(r8) :: procwt                      !weight of processor patches
  integer  :: proc_land_count(0:npes-1)
  integer  :: proc_patch_count(0:npes-1)
  integer, dimension(lsmlon, lsmlat) :: grid_to_chunk, grid_to_column
! --------------------------------------------------------------------



! ----------------------------------------------------------------------
! Initialize spmd arrays for number of land/patch points per processor
! ----------------------------------------------------------------------

  call spmd_init_patch



    grid_to_chunk(:,:) = -1
    grid_to_column(:,:) = -1
    do i = 1, nchunks
       do j = 1, chunk_ncols(i)
          grid_to_chunk(chunk_lon(j,i), chunk_lat(j,i)) = i
          grid_to_column(chunk_lon(j,i), chunk_lat(j,i)) = j
       enddo
    enddo

  npatch = 0
  nland = 0
  proc_landpts(:) = 0
  proc_patchpts(:) = 0
  do j = 1, lsmlat
     do i = 1, numlon(j)
        if (landmask(i,j) == 1) then         
           nland = nland+1                               !land point number
           k = chunk_pid(grid_to_column(i,j), &
                grid_to_chunk(i,j))
           proc_landpts(k) = proc_landpts(k) + 1
           do m = 1, maxpatch                
              if (wtxy(i,j,m) > 0.) then
                  npatch = npatch+1   !subgrid patch number
                  proc_patchpts(k) = proc_patchpts(k) + 1
              endif
           end do
        end if
     end do
  end do
  numland = nland
  numpatch = npatch

  if (masterproc) then
     write (6,*)' Surface Grid Characteristics'
     write (6,*)'   longitude points         = ',lsmlon
     write (6,*)'   latitude points          = ',lsmlat
     write (6,*)'   minimum longitude points = ',minval(numlon)
     write (6,*)'   maximum longitude points = ',maxval(numlon)
     write (6,*)'   land points on grid                    = ',numland
     write (6,*)'   total points including subgrid patches = ',numpatch
     write (6,*)
  endif
 
  proc_landi(0) = 1
  proc_landf(0) = proc_landpts(0)
  proc_patchi(0) = 1
  proc_patchf(0) = proc_patchpts(0) 
  do i = 1, npes - 1
    proc_landi(i) = proc_landf(i-1) + 1
    proc_landf(i) = proc_landi(i) - 1 + proc_landpts(i)
    proc_patchi(i) = proc_patchf(i-1) + 1
    proc_patchf(i) = proc_patchi(i) - 1 + proc_patchpts(i)
  enddo

  begpatch = proc_patchi(iam)
  endpatch = proc_patchf(iam)
  begland = proc_landi(iam)
  endland = proc_landf(iam)

  write(6,*)'iam= ',iam,' begland = ',proc_landi(iam) ,' endland = ',proc_landf(iam),&
       ' total land points = ',proc_landpts(iam)
  write(6,*)'iam= ',iam,' begpatch= ',proc_patchi(iam),' endpatch= ',endpatch, &
       ' total patch pionts= ',proc_patchpts(iam)

  call mapvar_ini
  call clm_varder_ini
  call histvar_ini
  call monthveg_ini

  nlchunks=sizeof(lchunk_to_chunk)/sizeof(lchunk_to_chunk(1))
  do i = 0, npes - 1
    proc_land_count(i) = proc_landi(i)
    proc_patch_count(i) = proc_patchi(i)
  enddo
  landvec%wtxy(:,:) = 0._r8
  patchvec%wtxy(:)  = 0._r8
# 413 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/main/clm_map.F90"
  do j = 1, lsmlat
     do i = 1, numlon(j)
        if (landmask(i,j) == 1) then                 
           k = chunk_pid(grid_to_column(i,j), &
                grid_to_chunk(i,j))
           landvec%ixy(proc_land_count(k)) = i
           landvec%jxy(proc_land_count(k)) = j
           do m = 1, maxpatch                           
              if (wtxy(i,j,m) > 0.) then                
                  landvec%patch(proc_land_count(k),m) = proc_patch_count(k)  !land subgrid patch number
                  landvec%wtxy(proc_land_count(k), m) = wtxy(i,j,m)   !land subgrid weights
                  patchvec%ixy(proc_patch_count(k)) = i             !patch longitude index
                  patchvec%jxy(proc_patch_count(k)) = j             !patch latitude index
                  patchvec%mxy(proc_patch_count(k)) = m             !patch subgrid index of lnd point
                  patchvec%wtxy(proc_patch_count(k)) = wtxy(i,j,m)  !patch weight
                  patchvec%land(proc_patch_count(k)) = proc_land_count(k) !patch land index
                  proc_patch_count(k) = proc_patch_count(k) + 1
              endif
           enddo
           proc_land_count(k) = proc_land_count(k) + 1
        endif
    enddo
 enddo


  !do i = 0, npes - 1
  !  write(*,*) "[TAT]", iam, i, proc_land_count(i), proc_patch_count(i) 
  !enddo

  do l = 1, numland
     n = 0
     do m = maxpatch, 1, -1
        if (landvec%wtxy(l,m) > 0.) n = m
     end do
     if (n == 0) then
        write (6,*) 'CLM_MAP error: n = 0'
        call endrun
     end if
     do m = 1, maxpatch
        if (landvec%wtxy(l,m) == 0.) then
           landvec%patch(l,m) = landvec%patch(l,n)
        endif
     end do
  end do

! Error check: make sure weights sum to one for each land cell

  do l = 1, numland
     i = landvec%ixy(l)
     j = landvec%jxy(l)
     sumwts = 0
     do m = 1, maxpatch
        sumwts = sumwts + landvec%wtxy(l,m)
     end do
     if (abs(sumwts-1.) > 1.e-05) then
        write (6,*) 'CLM_MAP error: weights do not sum to 1'
        write (6,*) 'lon = ',i,' lat = ',j,' :sum = ',sumwts
        call endrun
     end if
  end do

# 482 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/main/clm_map.F90"

  return
end subroutine
