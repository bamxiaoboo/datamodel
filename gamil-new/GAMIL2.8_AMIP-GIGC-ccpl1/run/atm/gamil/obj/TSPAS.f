# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/TSPAS.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/TSPAS.F"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/TSPAS.F" 2

# 1 "./params.h" 1
# 15 "./params.h"
 
 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/TSPAS.F" 2

       ! -----------------------------------------------------------------------
       ! Description:
       !
       !    This is the main subroutine for Two-tep Shape-Preserving Advection
       !    Scheme (TSPAS) developed by Rucong Yu. The 2D advection is performed
       !    in the global spherical geometry with the even-area weighting
       !    C-grid mesh.
       !
       ! Revisions:
       !
       !   2012-05~06: Correct the problem of inconservation by Li Dong with the
       !     following changes:
       !       - the position of weights (WTGU, WTGV)
       !       - the usage of SINU and SINV
       !       - the area at Poles
       ! -----------------------------------------------------------------------

      SUBROUTINE TSPAS(Q, U, V, SINU, SINV, WTGU, WTGV, DTDLT, DTDLN)

      use pmgrid, only: beglatexdyn, endlatexdyn, plon, plat
      use mpi_gamil

      use comfm1, only: cu, cv, tuu, tvv, A, QHSTAR, VSTAR, BETA, FX, FY

      IMPLICIT NONE


# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/PARADYN" 1



!     Define the parameters related to the model resolution


      INTEGER
     _        IM   ! the grid number along the longitude
     _       ,NL   ! the vertical layers
     _       ,NZ   ! NZ = NL + 1, considering the adding boundary in the top atmosphere
     _       ,NA

      PARAMETER(IM=128,NL=26,NZ=NL+1)

!     Define the paramters about the earth and the atmosphere, required by
!     the model atmosphere
!
      REAL*8
     _       RAD    ! the earth radius
     _      ,OMGA   ! the angular velocity of the earth	rotation
     _      ,GRAVIT ! the gravity
     _      ,RD     ! the dry air specific gas constant
     _      ,CP     ! specific heat at constant pressure
     _      ,CPD    ! specific heat at constant pressure
     _      ,CAPA   ! CAPA=RD/CP
!     _      ,P0    ! The sea level pressure of the standard atmosphere
!     _      ,T0    ! The sea level temperature of the standard atmosphere
     _      ,PI     ! the ratio of the circumference of a circle to its diameter
     _      ,PEALIB ! the maxium pressure of the standard atmoshere
     _      ,DPALIB ! the interval of two adjoining levels
!
      PARAMETER(RAD=6371000.0D0, OMGA=0.7292D-4, GRAVIT=9.806D0
!     _         ,RD =287.0D0,CP=1004.6D0,CAPA=RD/CP,T0=288.15D0
!     _         ,P0 =1013.25D0, PI=3.141592653589793D0)
     _         ,RD =287.0D0,CP=1004.6D0,CAPA=RD/CP,CPD=CP
     _         ,PI=3.141592653589793D0)
!      PARAMETER ( PEALIB=1160.0D0,DPALIB=2.5D0,NA=PEALIB/DPALIB )
*     PARAMETER ( PEALIB=1160.0D0,DPALIB=5.0D0,NA=PEALIB/DPALIB )
      PARAMETER ( PEALIB=1160.0D0,DPALIB=0.5D0,NA=PEALIB/DPALIB )
!
# 31 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/TSPAS.F" 2
       ! -----------------------------------------------------------------------
       ! INPUT
      REAL*8 U(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
      REAL*8 V(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
      REAL*8 SINU(beglatexdyn:endlatexdyn), SINV(beglatexdyn:endlatexdyn)
      REAL*8 WTGU(beglatexdyn:endlatexdyn), WTGV(beglatexdyn:endlatexdyn)
      REAL*8 DTDLN(beglatexdyn:endlatexdyn), DTDLT(beglatexdyn:endlatexdyn)
       ! -----------------------------------------------------------------------
       ! OUTPUT
      REAL*8 Q(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
       ! -----------------------------------------------------------------------
      REAL*8 QMIN, QMAX
      REAL*8 C(ilbnd:ihbnd,2,NL)
      REAL*8 GAMA, CXSTAR, CYSTAR, CX, CY, TEMP1, TEMP2, TEMP3, TEMP4
      REAL*16 XNP(NL), XNP_TMP(8,NL), XSP(NL), XSP_TMP(8,NL), TMP_SUM
      REAL*8  QH_pole_1(NL), QH_pole_1_south(NL), QH_pole_1_north(NL)
      REAL*8  ZERO, HALF, FOURTH, EPSM
      DATA ZERO, HALF, FOURTH, EPSM /0.0D0,0.5D0,0.25D0,1.0D-80/
      INTEGER I, J, K
      INTEGER comm_request1, comm_request2, comm_request3
       ! -----------------------------------------------------------------------

      call gamil_arrays_comm(COMM_TO_RIGHT, 1, Q(:,beglatexdyn,1),
     $                       request_id=comm_request1)
      call gamil_arrays_comm(COMM_TO_TOP,   1, Q(:,beglatexdyn,1),
     $                       request_id=comm_request2)

!$OMP PARALLEL
!$OMP DO PRIVATE (J,I)
      DO J = jbeg0, jend0
         IF (J .GE. 2 .AND. J .LE. 60-1) THEN
            DO I = beglonex, endlonex
               CU(I,J) = DTDLN(J)*SINU(J)
            END DO
         ELSE
            DO I = beglonex, endlonex
               CU(I,J) = ZERO
            END DO
         END IF
      END DO
!$OMP END DO NOWAIT
!$OMP DO PRIVATE (I,J)      
      DO J = jbeg0, jend1
         DO I = beglonex, endlonex
            CV(I,J) = DTDLT(J)*SINV(J)
         END DO
      END DO
!$OMP END PARALLEL
       ! -----------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (I,J,K)
      DO K = 1, NL
         DO J = jbeg0, jend0
            DO I = beglonex, endlonex
               TUU(I,J,K) = CU(I,J)*U(I,J,K)
               TVV(I,J,K) = CV(I,J)*V(I,J,K)
            END DO
         END DO
      END DO

      call wait_icomm_request(comm_request1)
      call wait_icomm_request(comm_request2)
      call gamil_arrays_comm(COMM_TO_LEFT, 1, TUU(:,beglatexdyn,1),
     $    Q(:,beglatexdyn,1), request_id=comm_request1)
      call gamil_arrays_comm(COMM_TO_BOT,  1, TVV(:,beglatexdyn,1),
     $    Q(:,beglatexdyn,1), request_id=comm_request2)
       ! -----------------------------------------------------------------------
       ! LAX-WENDROFF SCHEME
       ! =======================================================================
       ! CALCULATE FLUXES
!$OMP PARALLEL DO PRIVATE (I,J,K)
      DO K = 1, NL
         DO J = jbeg1, jend1
            DO I = ibeg1, iend1
               FX(I,J,K) = HALF*TUU(I,J,K)*(Q(I,J,K)+Q(I-1,J,K))
     $         -HALF*TUU(I,J,K)*TUU(I,J,K)*(Q(I,J,K)-Q(I-1,J,K))/SINU(J)
            END DO
         END DO
         DO J = jbeg0, jend1
            DO I = ibeg1, iend1
               FY(I,J,K) = HALF*TVV(I,J,K)*(Q(I,J+1,K)+Q(I,J,K))
     $         -WTGV(J)*HALF*TVV(I,J,K)*TVV(I,J,K)*(Q(I,J+1,K)-Q(I,J,K))/SINV(J)
            END DO
         END DO
      END DO
      call wait_icomm_request(comm_request1)
      call wait_icomm_request(comm_request2)
      call gamil_arrays_comm(COMM_TO_LEFT, 1, FX(:,beglatexdyn,1),
     $                       request_id=comm_request1)
      call gamil_arrays_comm(COMM_TO_BOT,  1, FY(:,beglatexdyn,1),
     $                       request_id=comm_request2)
       ! =======================================================================
       ! CALCULATE INTERMEDIATE QHSTAR
!$OMP PARALLEL DO PRIVATE (I,J,K,TEMP1,TEMP2,TEMP3,TEMP4,GAMA)
      DO K = 1, NL
         DO J = jbeg1, jend1
            DO I = ibeg1, iend1
               TEMP1 = ABS(TUU(I,J,K)/SINU(J))*(1-ABS(TUU(I,J,K)/SINU(J)))
               TEMP2 = ABS(TUU(I+1,J,K)/SINU(J))*(1-ABS(TUU(I+1,J,K)/SINU(J)))
               TEMP3 = ABS(TVV(I,J-1,K)/SINV(J-1))*(1-ABS(TVV(I,J-1,K)/SINV(J-1)))
               TEMP4 = ABS(TVV(I,J,K)/SINV(J))*(1-ABS(TVV(I,J,K)/SINV(J)))
               GAMA = MAX(TEMP1, TEMP2, TEMP3, TEMP4)
               BETA(I,J,K) = 2.0D0/(2.0D0-2.0D0*GAMA)
            END DO
         END DO
      END DO
      call wait_icomm_request(comm_request1)
      call wait_icomm_request(comm_request2)
!$OMP PARALLEL DO PRIVATE (I,J,K)
      DO K = 1, NL
         DO J = jbeg1, jend1
            DO I = ibeg1, iend1
               QHSTAR(I,J,K) = Q(I,J,K)-BETA(I,J,K)*(FX(I+1,J,K)-FX(I,J,K)
     $              +WTGU(J)*(FY(I,J,K)-FY(I,J-1,K)))/SINU(J)
            END DO
         END DO
      END DO
       ! =======================================================================
       ! HANDLE POLES
      IF (jbeg0 .EQ. 1) THEN
!$OMP PARALLEL DO PRIVATE (I,J,K,TMP_SUM)
         DO K = 1, NL
            TMP_SUM = ZERO
            DO I = ibeg1, iend1
               TMP_SUM = TMP_SUM+BETA(I,2,K)*FY(I,1,K)
            END DO
            XNP_TMP(1,K) = TMP_SUM
         END DO
         DO K = 1, NL
            XNP(K) = XNP_TMP(1,K)
            QH_pole_1(K) = Q(ibeg1,jbeg0,K)
         END DO
         call gamil_sum_pole_data_phys(1, XNP, NL)
         call broadcast_lon_data(2, 1, QH_pole_1, NL) 
!$OMP PARALLEL DO PRIVATE (I,J,K)
         DO K = 1, NL
            XNP_TMP(1,K) = QH_pole_1(K)-XNP(K)/SINV(1)*4*WTGV(1)/plon
            DO I = ibeg1, iend1
               QHSTAR(I,1,K) = XNP_TMP(1,K)
            END DO
         END DO
      END IF

      IF (jend0 .eq. 60) then
!$OMP PARALLEL DO PRIVATE (I,J,K,TMP_SUM)
         DO K = 1, NL
            TMP_SUM = ZERO
            DO I = ibeg1, iend1
               TMP_SUM = TMP_SUM+BETA(I,60-1,K)*FY(I,60-1,K)
            END DO
            XSP_TMP(1,K) = TMP_SUM
         END DO
         DO K = 1, NL
            XSP(K) = XSP_TMP(1,K)
            QH_pole_1(K) = Q(ibeg1,jend0,K)
         END DO
         call gamil_sum_pole_data_phys(60, XSP, NL)
         call broadcast_lon_data(2, 60, QH_pole_1, NL) 
!$OMP PARALLEL DO PRIVATE (I,J,K)
         DO K = 1, NL
            XSP_TMP(1,K) = QH_pole_1(K)+XSP(K)/SINV(60-1)*4*WTGV(60-1)/plon
            DO I = ibeg1, iend1
               QHSTAR(I,60,K) = XSP_TMP(1,K)
            END DO
         END DO
      END IF
       ! -----------------------------------------------------------------------
       ! TSPR JUDGEMENT (A <= 0 IS GOOD)
!$OMP PARALLEL DO PRIVATE (I,J,K,QMIN,QMAX,TEMP1,TEMP2,TEMP3,TEMP4,CXSTAR,CX)
      DO K = 1, NL
         DO J = jbeg0, jend0
            DO I = ibeg1, iend1
               QMIN = 1.0E15
               QMAX = -1.0E15
               IF (J .EQ. 1) THEN
                  QMIN = MIN(Q(I,J,K), Q(I,J+1,K), QMIN)
                  QMAX = MAX(Q(I,J,K), Q(I,J+1,K), QMAX)
               ELSE IF (J .EQ. 60) THEN
                  QMIN = MIN(Q(I,J,K), Q(I,J-1,K), QMIN)
                  QMAX = MAX(Q(I,J,K), Q(I,J-1,K), QMAX)
               ELSE
                  QMIN = MIN(Q(I,  J,K), Q(I+1,  J,K), Q(I-1,J,K),
     $                 Q(I,J+1,K), Q(  I,J-1,K), QMIN)
                  QMAX = MAX(Q(I,  J,K), Q(I+1,  J,K), Q(I-1,J,K),
     $                 Q(I,J+1,K), Q(  I,J-1,K), QMAX)
               END IF
               A(I,J,K) = (QHSTAR(I,J,K)-QMAX)*(QHSTAR(I,J,K)-QMIN)
            END DO
         END DO
      END DO
      call gamil_arrays_comm(COMM_TO_RIGHT, 1, A(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_TOP, 1, A(:,beglatexdyn,1))
       ! ------------------------------------------------------------------------
       ! UPWIND SCHEME
       ! =======================================================================
       ! CALCULATE FLUXES
!$OMP PARALLEL DO PRIVATE (I,J,K,TEMP1,TEMP2,TEMP3,TEMP4,CXSTAR,CX,CYSTAR,CY)
      DO K = 1, NL
         DO J = jbeg1, jend1
            DO I = ibeg1, iend1
               TEMP1 = (ABS(A(I-1,J,K))+A(I-1,J,K))/(ABS(A(I-1,J,K))+EPSM)
               TEMP2 = (ABS(A(I,J,K))+A(I,J,K))/(ABS(A(I,J,K))+EPSM)
               TEMP3 = (ABS(A(I-1,J,K))+A(I-1,J,K))*(ABS(A(I,J,K))+A(I,J,K))
               TEMP4 = ABS(A(I-1,J,K))*ABS(A(I,J,K))+EPSM
               CXSTAR = HALF*(TEMP1+TEMP2)-FOURTH*TEMP3/TEMP4
               CX = CXSTAR+(1-CXSTAR)*ABS(TUU(I,J,K)/SINU(J))
               FX(I,J,K) = HALF*TUU(I,J,K)*(Q(I,J,K)+Q(I-1,J,K))
     $              -HALF*ABS(CX*TUU(I,J,K))*(Q(I,J,K)-Q(I-1,J,K))
            END DO
         END DO
         DO J = jbeg0, jend1
            DO I = ibeg1, iend1
               TEMP1 = (ABS(A(I,J,K))+A(I,J,K))/(ABS(A(I,J,K))+EPSM)
               TEMP2 = (ABS(A(I,J+1,K))+A(I,J+1,K))/(ABS(A(I,J+1,K))+EPSM)
               TEMP3 = (ABS(A(I,J,K))+A(I,J,K))*(ABS(A(I,J+1,K))+A(I,J+1,K))
               TEMP4 = ABS(A(I,J,K))*ABS(A(I,J+1,K))+EPSM
               CYSTAR = HALF*(TEMP1+TEMP2)-FOURTH*TEMP3/TEMP4
               CY = CYSTAR+(1-CYSTAR)*ABS(TVV(I,J,K)/SINV(J))
               FY(I,J,K) = HALF*TVV(I,J,K)*(Q(I,J+1,K)+Q(I,J,K))
     $              -HALF*ABS(CY*TVV(I,J,K))*(Q(I,J+1,K)-Q(I,J,K))
            END DO
         END DO
      END DO
      call gamil_arrays_comm(COMM_TO_LEFT, 1, FX(:,beglatexdyn, 1))
      call gamil_arrays_comm(COMM_TO_BOT,  1, FY(:,beglatexdyn, 1))
       ! -----------------------------------------------------------------------
       ! HANDLE POLES
!$OMP PARALLEL DO PRIVATE (I,J,K,TMP_SUM)
      DO K = 1, NL
         IF (jbeg0 .EQ. 1) THEN
            TMP_SUM = ZERO
            DO I = ibeg1, iend1
               TMP_SUM = TMP_SUM+FY(I,1,K)
            END DO
            XNP_TMP(1,K) = TMP_SUM
         END IF
         IF (jend0 .EQ. 60) THEN
            TMP_SUM = ZERO
            DO I = ibeg1, iend1
               TMP_SUM = TMP_SUM+FY(I,60-1,K)
            ENDDO
            XSP_TMP(1,K) = TMP_SUM
         END IF
      END DO
      DO K = 1, NL
         XNP(K) = XNP_TMP(1,K)
         XSP(K) = XSP_TMP(1,K)
      END DO
      call gamil_sum_pole_data_phys(   1, XNP, NL)
      call gamil_sum_pole_data_phys(60, XSP, NL)

      IF (jbeg0 .eq. 1) THEN
         DO K = 1, NL
            QH_pole_1_north(K) = Q(ibeg1,jbeg0,K)
         END DO
         call broadcast_lon_data(2, 1, QH_pole_1_north, NL) 
      END IF
      IF (jend0 .eq. 60) then
         DO K = 1, NL
            QH_pole_1_south(K) = Q(ibeg1,jend0,K)
         END DO
         call broadcast_lon_data(2, 60, QH_pole_1_south, NL) 
      END IF


!$OMP PARALLEL DO PRIVATE (I,J,K,TMP_SUM)
      DO K = 1, NL
         DO J = jbeg1, jend1
            DO I = ibeg1, iend1
               Q(I,J,K) = Q(I,J,K)-(FX(I+1,J,K)-FX(I,J,K)
     $              +WTGU(J)*(FY(I,J,K)-FY(I,J-1,K)))/SINU(J)
            END DO
         END DO
      END DO

!$OMP PARALLEL DO PRIVATE (I,J,K,TEMP1,TEMP2,TEMP3)
      DO K = 1, NL
         IF (jbeg0 .EQ. 1) THEN
            XNP_TMP(1,K) = QH_pole_1_north(K)-XNP(K)/SINV(1)*4*WTGV(1)/plon
            DO I = beglonex, endlonex
               Q(I,1,K) = XNP_TMP(1,K)
            END DO
         END IF

         IF (jend0 .EQ. 60) THEN
            XSP_TMP(1,K) = QH_pole_1_south(K)+XSP(K)/SINV(60-1)*4*WTGV(60-1)/plon
            DO I = beglonex, endlonex
               Q(I,60,K) = XSP_TMP(1,K)
            END DO
         END IF
      END DO

      call gamil_arrays_comm(COMM_TO_LEFT,  1, Q(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_RIGHT, 1, Q(:,beglatexdyn,1))

      RETURN
      END

