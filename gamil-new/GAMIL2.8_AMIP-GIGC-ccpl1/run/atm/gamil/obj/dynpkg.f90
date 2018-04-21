# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/dynpkg.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/dynpkg.F90"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/dynpkg.F90" 2

# 1 "./params.h" 1
# 15 "./params.h"
 
 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/dynpkg.F90" 2

subroutine dynpkg(dtdy, nseq, dsghl)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Driving routines for dynamics and transport.
! 
!-----------------------------------------------------------------------

    use shr_kind_mod, only: r8 => shr_kind_r8
    use pmgrid
    use prognostics
    use qadv       
    use commap
    use stdatm  
    use comfm1
    use comhd
    use fspan      !!(wh 2003.11.04)
    use mpi_gamil
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------

# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/control/comhyb.h" 1
!----------------------------------------------------------------------- 
! 
! Purpose: Hybrid level definitions: p = a*p0 + b*ps
!          interfaces   p(k) = hyai(k)*ps0 + hybi(k)*ps
!          midpoints    p(k) = hyam(k)*ps0 + hybm(k)*ps
! 
!-----------------------------------------------------------------------
!!!  vertical level definitions in LASG dynamical core: p = pes*sigma + pt
!!!        interfaces   ply(k) = ps*sig (k) + pmtop
!!!        midpoints    ply(k) = ps*sigl(k) + pmtop
!!!---------------------------------------------------------------------
!!!(wanhui 2003.04.30)
!!!(wanhui 2003.10.23)  (std.atm. variables removed)

      real(r8) hyai(plevp)       ! ps0 component of hybrid coordinate - interfaces
      real(r8) hybi(plevp)       ! ps component of hybrid coordinate - interfaces
      real(r8) hyam(plev)        ! ps0 component of hybrid coordinate - midpoints
      real(r8) hybm(plev)        ! ps component of hybrid coordinate - midpoints

!!    real(r8) hybd(plev)        ! difference  in b (hybi) across layers
      real(r8) hypi(plevp)       ! reference pressures at interfaces
      real(r8) hypm(plev)        ! reference pressures at midpoints
!!    real(r8) hypd(plev)        ! reference pressure layer thickness

      real(r8) ps0         ! base state sfc pressure for level definitions
!!    real(r8) psr         ! reference surface pressure for linearization
!!    real(r8) prsfac      ! log pressure extrapolation factor (time, space independent)

!!    integer nprlev       ! number of pure pressure levels at top

      real(r8) :: pmtop              !
      real(r8) :: sig (plevp)        !
      real(r8) :: sigl(plev)         !  fm2003 VPAR variables
      real(r8) :: dsig(plev)         !

!!(wanhui 2003.10.23)
!!------------------------------------------------------------
!!    real(r8) :: tbb (plevstd)         !
!!    real(r8) :: hbb (plevstd)         !
!!    real(r8) :: cbb (plevstd)         !
!!    real(r8) :: dcbb(plevstd)         !  fm2003 std. atm.
!!    real(r8) :: p00, t00              !
!!    real(r8) :: psb (plond,plat)      !
!!    real(r8) :: tsb (plond,plat)      !
!!------------------------------------------------------------
!!(2003.10.23)(these variables are in module stdatm now)


!!      common /comhyb/ hyai ,hyam  ,hybi ,hybm
!!      common /comhyb/ hybd ,hypi ,hypm  ,hypd
!!      common /comhyb/ ps0         ,psr         ,prsfac      ,nprlev

      common /comhyb/ hyai ,hybi ,hyam ,hybm
      common /comhyb/ hypi ,hypm
      common /comhyb/ ps0  ,pmtop
      common /comhyb/ sig  ,sigl , dsig
!!    common /comhyb/ tbb  ,hbb  , cbb  ,dcbb  ,p00 ,t00 ,psb ,tsb    !!(wh 2003.10.23)
 
# 27 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/dynpkg.F90" 2
!----------------------------------------------------------------------

# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/control/comctl.h" 1
!----------------------------------------------------------------------- 
! 
! Purpose: Model control variables
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
!! (wh 2003.04.30)
!! (wh 2003.12.27)

      common /comctl/ itsst   ,nsrest  ,iradsw  ,iradlw  ,iradae
      common /comctl/ anncyc  ,nlend   ,nlres   ,nlhst   ,lbrnch
      common /comctl/ ozncyc  ,sstcyc  ,icecyc
      common /comctl/ adiabatic
      common /comctl_l1/ flxave
      common /comctl/ trace_gas, trace_test1,trace_test2, trace_test3
!!    common /comctl/ readtrace,ideal_phys, nsplit, iord, jord, kord, use_eta, aqua_planet
      common /comctl/ readtrace,ideal_phys,                                    aqua_planet
      common /comctl/ doRamp_ghg, doRamp_so4, doRamp_scon, fullgrid, doIPCC_so4, &
                      doCmip5_scon,doCmip5_ghg  !!(wh)
      common /comctl/ print_step_cost
      common /comctl/ doabsems, dosw, dolw, indirect

!!    common /comctl_r8/ divdampn, precc_thresh, precl_thresh
      common /comctl_r8/           precc_thresh, precl_thresh

      integer itsst             ! Sea surf. temp. update freq. (iters)
      integer nsrest            ! Restart flag
      integer iradsw            ! Iteration freq. for shortwave radiation
      integer iradlw            ! Iteration freq. for longwave radiation
      integer iradae            ! Iteration freq. for absorptivity/emissivity

! f-v dynamics specific
! _ord = 1: first order upwind
! _ord = 2: 2nd order van Leer (Lin et al 1994)
! _ord = 3: standard PPM 
! _ord = 4: enhanced PPM (default)
!!      integer nsplit            ! Lagrangian time splits (Lin-Rood only)
!!      integer iord              ! scheme to be used in E-W direction
!!      integer jord              ! scheme to be used in N-S direction
!!      integer kord              ! scheme to be used for vertical mapping
!!      logical use_eta           ! Flag to use a's and b's set by dynamics/lr/set_eta.F90

      logical aqua_planet       ! Flag to run model in "aqua planet" mode

      logical anncyc            ! true => do annual cycle (otherwise perpetual)
      logical nlend             ! true => end of run
      logical nlres             ! true => continuation run
      logical nlhst             ! true => regeneration run
      logical lbrnch            ! true => branch run
      logical ozncyc            ! true => cycle ozone dataset
      logical sstcyc            ! true => cycle sst dataset
      logical icecyc            ! true => cycle ice fraction dataset
      logical adiabatic         ! true => no physics
      logical ideal_phys        ! true => run "idealized" model configuration
      logical(1) flxave            ! true => send to coupler only on radiation time steps

      logical trace_gas         ! true => turn on greenhouse gas code
      logical trace_test1       ! true => turn on test tracer code with 1 tracer
      logical trace_test2       ! true => turn on test tracer code with 2 tracers
      logical trace_test3       ! true => turn on test tracer code with 3 tracers
      logical readtrace         ! true => obtain initial tracer data from IC file

      logical doRamp_ghg        ! true => turn on ramping for ghg
      logical doRamp_so4        ! true => turn on ramping for so4
      logical doRamp_scon       ! true => turn on ramping for scon
      logical doIPCC_so4        ! true => turn on IPCC scenario for so4  !!(wh) 
      logical doCmip5_scon      !
      logical doCmip5_ghg       ! ljli2010-08-12
      logical fullgrid          ! true => no grid reduction towards poles

      logical print_step_cost   ! true => print per-timestep cost info

      logical doabsems          ! True => abs/emiss calculation this timestep
      logical dosw              ! True => shortwave calculation this timestep
      logical dolw              ! True => longwave calculation this timestep
      logical indirect          ! True => include indirect radiative effects of sulfate aerosols

!!    real(r8) divdampn         ! Number of days to invoke divergence damper
      real(r8) precc_thresh     ! Precipitation threshold for PRECCINT and PRECCFRQ
      real(r8) precl_thresh     ! Precipitation threshold for PRECLINT and PRECLFRQ
# 29 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/dynpkg.F90" 2
!------------------------------Arguments--------------------------------
    real(r8), intent(in) :: dtdy               ! timestep size ( dyn core )
    integer, intent(in) :: nseq
    real(r8), intent(in) :: dsghl(plev)

!---------------------------Local workspace-----------------------------

    integer i, j, k
    integer begj

	! TODE: Clarify this!
    REAL*8  tmp_3d(plond,plat,plev*20)

!----------------------------------------------------------
!....  PERFORM THE DYNAMIC INTEGRATON CYCLE
!----------------------------------------------------------

    begj = beglatexdyn

    call t_startf('DYFRAM')

    CALL DYFRAM2(NSEQ, DTDY, ITIME,                             &
        U, V, T, Q, WS, PES, WPA, GHS, GHI, PLY, TB,            &
        SU, SV, ST, SUT, SVT, STT,                              &
        NONOS, IORD, ISOR, EP, IPQ, DTDLN, DTDLT, DTDSG, DSGHL, &
        PMTOP, SIG, SIGL, DSIG,                                 &
        TBB, HBB, CBB, DCBB, PSB, TSB,                          &
        DY, WTGU(BEGJ), WTGV(BEGJ),                             &
        DX, SINU, SINV, OUX, OUY, OVX, OVY, FF, CUR,            &
        MM1, MP1, MM2, MP2, MM3, MP3, MDJ,                      &
        U0, V0, WS0, QT, DP, FAC, FBC, PP,                      &
        UUK, HHK, DUS, DPS2, PLY2, TB2, CB, DCB, CB0,           &
        P, C0, NIGW, UU, VV, TT, DPS0, DPS1, HPS, HH, TTZ,      &
        UZ, VZ, TTV, DPS, DU, DV, DTT, DU1, DV1, DTT1,          &
        UK, VK, TTK, PSK)

    call t_stopf('DYFRAM')
!
!----------------------------------------------------------
!....  DO FIRST HALF-STEP HORIZONTAL DIFFUSION
!----------------------------------------------------------

    if (.not. aqua_planet)  then
!$OMP PARALLEL DO PRIVATE (I, J, K)
        DO K = 1, plev
            DO J = jbeg0, jend0
                DO I = beglonex, endlonex
                    UK(I,J,K) = U(I,J,K)
                    VK(I,J,K) = V(I,J,K)
                    ttk(I,J,K) = T(I,J,K)
                    QK(I,J,K) = Q(I,J,K)
                END DO
            END DO
        END DO

        call t_startf('HDIFUS')

        CALL HDIFUS(U, V, T, Q, FRDT, FRDS, FRDU, &
            FRDV, FRDP, TB, PLY, DXVPN, DXVPS)
    
        call t_stopf('HDIFUS')

!$OMP PARALLEL DO PRIVATE (I, J, K)
        DO K = 1, plev
            DO J = jbeg0, jend0
                DO I = beglonex, endlonex
                    SU(I,J,K)=(U(I,J,K)-UK(I,J,K))/DTHDFS
                    SV(I,J,K)=(V(I,J,K)-VK(I,J,K))/DTHDFS
                    ST(I,J,K)=(T(I,J,K)-ttk(I,J,K))/DTHDFS
                    U(I,J,K)=UK(I,J,K)
                    V(I,J,K)=VK(I,J,K)
                    T(I,J,K)=ttk(I,J,K)
                    Q(I,J,K)=QK(I,J,K)
                END DO
            END DO
        END DO
    else
        write(6, *) "[Notice]: dynpkg: No horizontal diffusion."
        SU(:,:,:) = 0.0
        SV(:,:,:) = 0.0
        ST(:,:,:) = 0.0
    end if

    return
end subroutine dynpkg

