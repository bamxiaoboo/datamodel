# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DYFRAM.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DYFRAM.F"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DYFRAM.F" 2

# 1 "./params.h" 1
# 15 "./params.h"
 
 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DYFRAM.F" 2

      SUBROUTINE DYFRAM2(NSEQ, DTDY, ITIME,
     _     U, V, T, Q, WS, PES, WPA, GHS, GHI, PLY, TB,
     _     SU, SV, ST, SUT, SVT, STT,
     _     NONOS, IORD, ISOR, EP, IPQ, DTDLN, DTDLT, DTDSG, DSGHL,
     _     PMTOP, SIG, SIGL, DSIG,
     _     TBB, HBB, CBB, DCBB, PSB, TSB,
     _     DY, WTGU, WTGV,
     _     DX, SINU, SINV, OUX, OUY, OVX, OVY, FF, CUR,
     _     MM1, MP1, MM2, MP2, MM3, MP3, MDJ,
     _     U0, V0, WS0, QT, DP, FAC, FBC, PP,
     _     UUK, HHK, DUS, DPS2, PLY2, TB2, CB, DCB, CB0,
     _     P, C0, NIGW, UU, VV, TT, DPS0, DPS1, HPS, HH, TTZ,
     _     UZ, VZ, TTV, DPS, DU, DV, DTT, DU1, DV1, DTT1,
     _     UK, VK, TTK, PSK)

!
!     MAIN ROUTINE  OF THE DYNAMIC FRAME CALCULATION

      use pmgrid, only: beglatexdyn, endlatexdyn, plat
      use mpi_gamil

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
# 28 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DYFRAM.F" 2

!
      REAL*8 DTDY
      INTEGER*4 ITIME

      REAL*8 PMTOP, SIG(NZ), DSIG(NL), SIGL(NL)
      REAL*8 DX, DY
      REAL*8 WTGU(BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 WTGV(BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 SINU(BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 SINV(BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 OUX(BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 OUY(BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 OVX(BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 OVY(BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 FF(BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 CUR(BEGLATEXDYN:ENDLATEXDYN)

      INTEGER MM1(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN)     !    mm1=i-mdj      , if mm1<2, mm1=NX_LON-2+mm1
      INTEGER MP1(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN)     !    mp1=i+mdj      , if mp1>NX_LON-1, mp1=mp1-NX_LON+2
      INTEGER MM2(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN)     !    mm2=i-(mdj-1)/2, if mm2<2, mm2=NX_LON-2+mm2
      INTEGER MP2(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN)     !    mp2=i+(mdj+1)/2, if mp2>NX_LON-1, mp2=mp2-NX_LON+2
      INTEGER MM3(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN)     !    mm3=i-(mdj+1)/2, if mm3<2, mm3=NX_LON-2+mm3
      INTEGER MP3(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN)     !    mp3=i+(mdj-1)/2, if mp3>NX_LON-1, mp3=mp3-NX_LON+2
      INTEGER MDJ(BEGLATEXDYN:ENDLATEXDYN)        !    leaping span of the difference

      REAL*8 CBB(NA), TBB(NA), HBB(NA), DCBB(NA)
      REAL*8 PSB(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 TSB(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN)

      INTEGER NONOS, IORD, ISOR, IPQ(NX_LON)

      REAL*8 EP, DTDSG(NL)
      REAL*8 DTDLN(BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 DTDLT(BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 DSGHL(NL)

      REAL*8 U0(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 V0(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 WS0(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 CBT
      REAL*8 TB(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 WS(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NZ)
      REAL*8 PLY(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NZ)
      REAL*8 GHI(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NZ)

      REAL*8 GHS(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 DP(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN)

      REAL*8 PES(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 U(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 V(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 T(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 Q(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 QT(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)

C     SOURCE & SINK  OF THE MODEL PREDICTED VARIABLES   DYNAMICS/PHYS
      REAL*8 SUT(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 SVT(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 STT(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 SU(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 SV(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 ST(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)

      REAL*8 WPA(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 FAC(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NZ)
      REAL*8 FBC(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 PP(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN)

      REAL*8 UUK(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 HHK(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 DUS(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 DPS2(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN)

      REAL*8 PLY2(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 TB2(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 CB(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 DCB(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 CB0(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 P(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 C0(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN)
      INTEGER NIGW(BEGLATEXDYN:ENDLATEXDYN)

      REAL*8 UU(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 VV(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 TT(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 DPS0(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 DPS1(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 HPS(ILBND:IHBND, BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 HH(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NZ)
      REAL*8 TTZ(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NZ+1)
      REAL*8 UZ(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NZ+1)
      REAL*8 VZ(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NZ+1)
      REAL*8 TTV(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 DPS(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN)
      
      REAL*8 DU(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 DV(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 DTT(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 DU1(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 DV1(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 DTT1(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)

      REAL*8 UK(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 VK(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 TTK(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 PSK(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN)

      REAL*8 DQ(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 UQ(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 VQ(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 WQ(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 PQ(ILBND:IHBND,BEGLATEXDYN:ENDLATEXDYN,NL)

      REAL*8 WK5, WK6, PX, PY, QM
      INTEGER NSEQ, I, J, K, NCYC, KPP
      INTEGER BEGJ, ENDJ
!
!     START   THE DYNAMICAL INTEGRATING CYCLE
!

!$OMP PARALLEL DO PRIVATE (I,J)
      DO J = JBEG0, JEND0
         DO I = IBEG1, IEND1
            PP(I,J) = SQRT(PES(I,J))
         END DO
      END DO

      call gamil_arrays_comm(COMM_TO_RIGHT,1,pp(:,beglatexdyn))
      call gamil_arrays_comm(COMM_TO_TOP,1,pp(:,beglatexdyn))

!$OMP PARALLEL DO PRIVATE (I,J,K)
      DO K=1,NL
      DO J=jbeg0,jend0
      DO I=beglonex,endlonex
      IF (I .eq. 1 .or. I .eq. NX_LON) THEN
         WPA(I,J,K)=0.0
      ELSE
         QT(I,J,K)=Q(I,J,K)*pes(I,J)
         WPA(I,J,K)=0.0
      ENDIF
      END DO
      END DO
      END DO

!
!     TRANSFORM THE FORCING TERMS OF WIND, TEMPERATURE AND MOISTURE
!     TO THE CONTRIBUTIONS TO THE TENDENCIES OF THE BAISC MODEL
!     PREDICTION VARIABLES
!

!$OMP PARALLEL DO PRIVATE (I,J,K,PY,PX,WK5,KPP,WK6,cbt)
      DO K=1,NL
       DO J=jbeg0,jend0

         DO I=ibeg1,iend1
!
! 	 CALCULATING THE AVERAGE VALUE OF pp AT V-GRID.
!
         IF(J.LT.plat) THEN
 	    PY=0.5D0*(pp(I,J)+pp(I,J+1))
         ELSE
 	    PY=0.0D0
 	 ENDIF
!
!	 CALCULATING THE AVERAGE VALUE OF pp AT U-GRID.
!
         PX=0.5D0*(pp(I,J)+pp(I-1,J))
!
 	    sut(I,J,K)=PX*su(I,J,K)
            svt(I,J,K)=PY*sv(I,J,K)
!
!           CALCULATING STT=stt*pp*RD/cbt.
!
            WK5=(pes(I,J)*SIGL(K)+PMTOP)/DPALIB
            KPP=INT(WK5)
            WK6=WK5-DFLOAT(KPP)
            cbt=(1.D0-WK6)*CBB(KPP)+WK6*CBB(KPP+1)
            stt(I,J,K)=st(I,J,K)*pp(I,J)*RD/cbt
	 END DO
      END DO
!
      END DO

      call gamil_arrays_comm(COMM_TO_LEFT,1,QT(:,beglatexdyn,1),sut(:,beglatexdyn,1),svt(:,beglatexdyn,1),stt(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_RIGHT,1,QT(:,beglatexdyn,1),sut(:,beglatexdyn,1),svt(:,beglatexdyn,1),stt(:,beglatexdyn,1))

!
!     THE END OF THE TRANSFORMATION
!
      DO NCYC = 1, NSEQ
!
!     PREDICT DRY-ADIABATIC SYSTEM
!     ___________________________
!
!$OMP PARALLEL DO PRIVATE (I,J,K)
         DO K = 1, NL
            DO J = JBEG0, JEND0
               DO I = BEGLONEX, ENDLONEX
                  U0(I,J,K) = U(I,J,K)
                  V0(I,J,K) = V(I,J,K)
                  WS0(I,J,K) = WS(I,J,K)
                  FAC(I,J,K) = 1.0
                  FBC(I,J,K) = 1.0
               END DO
            END DO
         END DO

!$OMP PARALLEL DO PRIVATE (I,J)
         DO J = JBEG0, JEND0
           DO I = BEGLONEX, ENDLONEX
              FAC(I,J,NZ) = 1.0
           END DO
         END DO

         CALL T_STARTF("DYNAMICS")
         CALL DYNAMICS( DTDY,ITIME,NCYC,
     _               U,V,WS,WS0,PES,T,GHI,GHS,DP,SUT,SVT,STT,FAC,FBC,
     _               PMTOP,SIGL,DSIG,
     _               TBB,HBB,CBB,DCBB,PSB,TSB,
     _               DY,WTGU,WTGV,
     _               DX,SINU,SINV,OUX,OUY,OVX,OVY,FF,CUR,
     _               MM1,MP1,MM2,MP2,MM3,MP3,MDJ,
     _               UUK,HHK,DUS,DPS2,PLY2,TB2,CB,DCB,CB0,
     _               P,C0,NIGW,UU,VV,TT,DPS0,DPS1,HPS,HH,TTZ,
     _               UZ,VZ,TTV,DPS,DU,DV,DTT,DU1,DV1,DTT1,
     _               UK,VK,TTK,PSK,NCYC)
         CALL T_STOPF("DYNAMICS")

         ITIME = ITIME+INT(DTDY+0.1)
!
!     PREDICT WATER VAPOR MIXING RATIO
!     BY POSITIVE DEFINITE ADVECTION TRANSPORT ALGORITHM
!     ____________________________
!

         CALL T_STARTF("QPDATA1")
         CALL QPDATA1(QT, U0, V0, WS0, DSGHL, U, V, WS, DTDLN, DTDLT, 
     _                SINU, SINV, WTGU, WTGV, DTDSG, DQ, UQ, VQ, WQ)
         CALL T_STOPF("QPDATA1")

# 277 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DYFRAM.F"

!
!     GET THE TIME AVERAGED pp-SURFACE VERTICAL VELOCITY
!
!$OMP PARALLEL DO PRIVATE (I,J,K)
         DO K = 1, NL
            DO J = JBEG0, JEND0
               DO I = IBEG1, IEND1
                  WPA(I,J,K) = WPA(I,J,K)+(0.5*(WS(I,J,K+1)+WS(I,J,K))
     _                 *PES(I,J)+DP(I,J)*SIGL(K))/REAL(NSEQ)
               END DO
            END DO
         END DO

      END DO

      call gamil_arrays_comm(COMM_TO_LEFT,1,WPA(:,beglatexdyn,1)) 
      call gamil_arrays_comm(COMM_TO_RIGHT,1,WPA(:,beglatexdyn,1)) 

!     --------------------
!
!$OMP PARALLEL DO PRIVATE (I,J,K,WK5,KPP,WK6)
      DO J = jbeg0,jend0
      DO K=1,NL
      DO I=beglonex,endlonex
         Q(I,J,K)=QT(I,J,K)/pes(I,J)
         PLY(I,J,K)=pes(I,J)*SIGL(K)+PMTOP
         WK5=PLY(I,J,K)/DPALIB
         KPP=INT(WK5)
         WK6=WK5-DFLOAT(KPP)
         TB(I,J,K)=(1.D0-WK6)*TBB(KPP)+WK6*TBB(KPP+1)
      END DO
      END DO
      DO I=beglonex,endlonex
      PLY(I,J,NZ)=pes(I,J)+PMTOP
      END DO
      END DO

      RETURN
      END
