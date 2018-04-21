# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/control/rgrid.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/control/rgrid.F90"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/control/rgrid.F90" 2

module rgrid

    use pmgrid, only: plat, platd
    use pspect, only: pmmax, pmax
    use infnan, only: bigint

    implicit none

    integer :: nlon(plat)        = bigint ! num longitudes per latitude
    integer :: nlonex(platd)     = bigint ! num longitudes per lat (extended grid)
    integer :: beglatpair(pmmax) = bigint
    integer :: nmmax(plat/2)     = bigint
    integer :: wnummax(plat)     = bigint ! cutoff Fourier wavenumber




end module rgrid
