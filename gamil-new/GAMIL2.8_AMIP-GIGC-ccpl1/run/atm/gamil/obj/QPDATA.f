# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/QPDATA.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/QPDATA.F"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/QPDATA.F" 2

# 1 "./params.h" 1
# 15 "./params.h"
 
 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/QPDATA.F" 2

!##################################################################################
!! (2003.11.30)

      SUBROUTINE QPDATA1(QT, U0, V0, W0, DSGHL, U, V, WS,
     _                   DTDLN, DTDLT, SINU, SINV, WTGU, WTGV, DTDSG,
     _                   DQ, UQ, VQ, WQ)

      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil
!
!     PREDICT POSITIVE DEFINITE FIELD Q        DUE TO 3-D ADVECTION
!           1)  BY USING THE SCHEME     PROPOSED BY R.C.Yu
!           2)  BY USING THE SCHEME     PROPOSED BY P.K.Smolarkiewicz
!
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
# 21 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/QPDATA.F" 2

C
      REAL*8 U0(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
      REAL*8 V0(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
      REAL*8 W0(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
      REAL*8 U(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
      REAL*8 V(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
      REAL*8 WS(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NZ)
      REAL*8 QT(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
      REAL*8 DQ(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
      REAL*8 DSGHL(NL)
      REAL*8 SINU(beglatexdyn:endlatexdyn)
      REAL*8 SINV(beglatexdyn:endlatexdyn)
      REAL*8 WTGU(beglatexdyn:endlatexdyn)
      REAL*8 WTGV(beglatexdyn:endlatexdyn)
!
      INTEGER NONOS,IORD,ISOR
!
      REAL*8 DSNP,DSSP,DTDLN(beglatexdyn:endlatexdyn)
      REAL*8 DTDLT(beglatexdyn:endlatexdyn)
      REAL*8 GC(beglatexdyn:endlatexdyn),DTDSG(NL)
      REAL*8 UQ(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
      REAL*8 VQ(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
      REAL*8 WQ(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
      REAL*8 HALF,ONE
      DATA    HALF,ONE / 0.5E0,1.0E0 /
      INTEGER I,J,K
!
!
!     GET THE ADVECTION VELOCITY
!

!$OMP PARALLEL DO PRIVATE (I,J,K)
      DO K   = 1 ,NL
         DO J   = jbeg0 , jend0
            DO I   = beglonex , endlonex
               UQ(I,J,K)  = HALF*(U (I,J,K)+U0(I,J,K))
               VQ(I,J,K)  = HALF*(V (I,J,K)+V0(I,J,K))
               WQ(I,J,K)  = HALF*(WS(I,J,K)+W0(I,J,K))
               DQ(I,J,K)    = QT(I,J,K)
            ENDDO
         ENDDO
      ENDDO
!
!     PERFORM HORIZONTAL ADVECTION IN SPHERICAL GEOMETRY
!
      CALL TSPAS(QT, UQ, VQ, SINU, SINV, WTGU, WTGV, DTDLT, DTDLN)
!
!     PERFORM THE VERTICAL ADVECTION
!       BY  R.C.Yu
      CALL TSPASW(QT,WQ,DTDSG)
!
!     PERFORM VERTICAL REDISTRIBUTION TO AVOID NEGATIVE Q-H2O
!
      CALL AVNEGQ(QT,DSGHL)

      RETURN
      END
