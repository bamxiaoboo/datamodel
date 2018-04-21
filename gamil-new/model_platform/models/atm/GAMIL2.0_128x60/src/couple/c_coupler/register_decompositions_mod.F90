!***************************************************************
!  This is a source file of GAMIL, which registers all parallel 
!  decompositions into C-Coupler library. This file was initially 
!  finished by Dr. Li Liu. If you have any problem, please 
!  contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
!***************************************************************


module register_decompositions_mod

    public register_decompositions

contains

    subroutine register_decompositions
       use mpi_gamil
       use pmgrid
       use phys_grid
       use rgrid,          only: nlon                                                  ! reduced grid
       use c_coupler_interface_mod
       implicit none
       integer,allocatable :: decomp_cell_indexes(:)
       integer             :: n, i, j, startpoint, bufsize
       integer :: lchnk         ! indices
       integer :: ncol                  ! number of columns in current chunk
       integer :: lats(pcols)           ! array of latitude indices
       integer :: lons(pcols)           ! array of longitude indices

       bufsize=0
       do j=1,plat
          bufsize = bufsize + nlon(j)
       end do
       allocate(decomp_cell_indexes(bufsize))
       n = 0
       startpoint = 0
       do j=1,plat
          do i=1,nlon(j)
             if(get_chunk_owner_p(i,j) .eq. iam) then
                n=n+1
                decomp_cell_indexes(n) = startpoint + i
             end if
          enddo
          startpoint = startpoint + nlon(j)
       enddo
       call c_coupler_register_decomposition("gamil_gamil_grid_decomp", "gamil_grid", n, decomp_cell_indexes)
       deallocate(decomp_cell_indexes)

       call c_coupler_register_decomposition("gamil_2D_decomp_dyn", "gamil_grid", &
                                  (ihbnd-ilbnd+1)*(endlatexdyn-beglatexdyn+1), dyn_cell_global_index)

       bufsize = (endlatex-beglatex+1)*(ihbnd-ilbnd+1)
       allocate(decomp_cell_indexes(bufsize))
       decomp_cell_indexes=0
       do j=beglat,endlat
          do i=ibeg0,iend0
             if (i.ge.1 .and. i.le.NX_LON-2) then
                 decomp_cell_indexes((i-ilbnd+1)+(ihbnd-ilbnd+1)*(j-beglatex)) = i+(j-1)*(NX_LON-2)
             end if
          end do
       end do
       call c_coupler_register_decomposition("gamil_2D_decomp_prog", "gamil_grid", &
                                  (ihbnd-ilbnd+1)*(endlatex-beglatex+1), decomp_cell_indexes)
       deallocate(decomp_cell_indexes)

       allocate(decomp_cell_indexes(pcols*(endchunk-begchunk+1)))
       decomp_cell_indexes=0
       do lchnk = begchunk,endchunk
          ncol = get_ncols_p(lchnk)
          call get_lon_all_p(lchnk, ncol, lons)
          call get_lat_all_p(lchnk, ncol, lats)
          do i=1,ncol
              decomp_cell_indexes((lchnk-begchunk)*pcols+i)=(lats(i)-1)*(NX_LON-2)+lons(i)
          end do
       end do
       call c_coupler_register_decomposition("gamil_2D_decomp_phys", "gamil_grid", &
                                  pcols*(endchunk-begchunk+1), decomp_cell_indexes)

       deallocate(decomp_cell_indexes)
    end subroutine register_decompositions


end module register_decompositions_mod

