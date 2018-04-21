#include <misc.h>
#include <params.h>

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
#include <comhyb.h>
!----------------------------------------------------------------------
#include <comctl.h>
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

