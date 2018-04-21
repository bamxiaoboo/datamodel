# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/CONPDA.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/CONPDA.F"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/CONPDA.F" 2

# 1 "./params.h" 1
# 15 "./params.h"
 
 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/CONPDA.F" 2


      SUBROUTINE CONPDA(DTIMEQ, DX, DY, SINU, WTGV, DSIG, FIRST, ISOR,
     $                  IORD, IP, DTDLT, DTDLN, DTDSG)

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
# 14 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/CONPDA.F" 2

      REAL*8 DTIMEQ, DX, DY, FIRST
      REAL*8 SINU(beglatexdyn:endlatexdyn)
      REAL*8 WTGV(beglatexdyn:endlatexdyn)
      REAL*8 DSIG(NL)
      REAL*8 DTDLT(beglatexdyn:endlatexdyn)
      REAL*8 DTDLN(beglatexdyn:endlatexdyn)
      REAL*8 DTDSG(NL)
      INTEGER ISOR, IORD, IP(NX_LON)
      INTEGER I, J, K

      IF (FIRST .LE. 0.0D0) THEN
         IF (ISOR .EQ. 3) IORD = MAX0(IORD, 3)
         DO I = 1, NX_LON
           IP(I) = MOD(I+IM/2-1, IM)+1
         END DO
      END IF

      DO J = jbeg0, jend1
         DTDLT(J) = DTIMEQ/(RAD*DY)
      END DO
      DO J = jbeg1, jend1
         DTDLN(J) = DTIMEQ/(RAD*SINU(J)*DX)
      END DO
      DO K  = 1 ,NL
         DTDSG(K)  = DTIMEQ / DSIG(K)
      END DO

      RETURN
      END