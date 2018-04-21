# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/TSPASW.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/TSPASW.F"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/TSPASW.F" 2

# 1 "./params.h" 1
# 15 "./params.h"
 
 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/TSPASW.F" 2

C     =======================
      SUBROUTINE TSPASW(Q,W,DTDSG)
C     =======================
C
      use pmgrid, only: beglatexdyn,endlatexdyn, plat
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
# 14 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/TSPASW.F" 2

C
      REAL*8 Q(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
      REAL*8 W(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
      REAL*8 QWMIN,QWMAX
!
      REAL*8 WW(NL),FW(NZ),BETAW(NL),QW(NL),QWSTAR(NL)
     _      ,WSTAR(NL),AW(NL),HS(NL),HW(NL),DTDSG(NL)
      REAL*8 north_pole_q_nx(NL), south_pole_q_1(NL)
!
      REAL*8 GAMA,CWSTAR,CW,TEMP1,TEMP2,TEMP3,TEMP4
!
      REAL*8  ZERO,HALF,FOURTH,EPSM
      DATA ZERO,HALF,FOURTH,EPSM/ 0.0D0,0.5D0,0.25D0,
     $     1.0D-80/
*    $     1.0E-6/
      INTEGER I,J,K,IS,IT
      integer begj,endj
!
      DO K=1,NL
	HS(K)=1.0D0/DTDSG(K)
      ENDDO
C
      DO K=2,NL
	HW(K)=HALF*(HS(K)+HS(K-1))
      ENDDO
C

!$OMP PARALLEL DO PRIVATE (I,J,K,IS,IT,QW,WW,FW,GAMA,BETAW,QWSTAR,
!$   &                                       QWMIN,QWMAX,AW,TEMP1,TEMP2,TEMP3,TEMP4,
!$   &                                       CWSTAR,CW,WSTAR)
      DO J=jbeg0,jend0

      IS = NX_LON
      IT = -1
      IF(J.EQ.1) THEN
        IF (endlonex .eq. NX_LON) THEN
	      IS=NX_LON
	      IT=NX_LON
        ENDIF
      ELSE IF(J.EQ.60) THEN
        IF (beglonex .eq. 1) THEN
	      IS=1
	      IT=1
        ENDIF
      ELSE
	    IS=beglonex
	    IT=endlonex
      ENDIF

C
      DO I=IS,IT
C
      DO K=1,NL
	QW(K)=Q(I,J,K)
	WW(K)=W(I,J,K)
      ENDDO
C
      DO K=2,NL
	FW(K)=HALF*WW(K)*(QW(K)+QW(K-1))
     $ -HALF*WW(K)*WW(K)*(QW(K)-QW(K-1))/HW(K)
      ENDDO
C
      DO K=2,NL-1
 	TEMP1=ABS(WW(K)/HW(K))*(1-ABS(WW(K)/HW(K)))
	TEMP2=ABS(WW(K+1)/HW(K+1))*(1-ABS(WW(K+1)/HW(K+1)))
        GAMA=MAX(TEMP1,TEMP2)
	BETAW(K)=2.0D0/(2.0D0-GAMA)
	QWSTAR(K)=QW(K)-BETAW(K)*(FW(K+1)-FW(K))*DTDSG(K)
      ENDDO
C
      QWSTAR(1)=QW(1)-BETAW(2)*FW(2)*DTDSG(1)
      QWSTAR(NL)=QW(NL)+BETAW(NL-1)*FW(NL)*DTDSG(NL)
C
      DO K=1,NL
       QWMIN=1.0E15
       QWMAX=-1.0E15
	IF(K.EQ.1) THEN
	  QWMIN=MIN(QW(K),QW(K+1),QWMIN)
	  QWMAX=MAX(QW(K),QW(K+1),QWMAX)
CCCC	ELSE IF(J.EQ.NL) THEN
	ELSE IF(K.EQ.NL) THEN
	  QWMIN=MIN(QW(K),QW(K-1),QWMIN)
	  QWMAX=MAX(QW(K),QW(K-1),QWMAX)
	ELSE
          QWMIN=MIN(QW(K),QW(K-1),QW(K+1),QWMIN)
	  QWMAX=MAX(QW(K),QW(K-1),QW(K+1),QWMAX)
	ENDIF
	  AW(K)=(QWSTAR(K)-QWMIN)*(QWSTAR(K)-QWMAX)
      ENDDO
C
      DO K=2,NL
	TEMP1=(ABS(AW(K))+AW(K))/(AW(K)+EPSM)
	TEMP2=(ABS(AW(K-1))+AW(K-1))/(AW(K-1)+EPSM)
	TEMP3=(ABS(AW(K))+AW(K))*(ABS(AW(K-1))+AW(K-1))
	TEMP4=ABS(AW(K))*ABS(AW(K-1))+EPSM
	CWSTAR=HALF*(TEMP1+TEMP2)-FOURTH*TEMP3/TEMP4
	CW=CWSTAR+(1-CWSTAR)*ABS(WW(K)/HW(K))
	WSTAR(K)=CW*WW(K)
      ENDDO
C
      DO K=2,NL
	FW(K)=HALF*WW(K)*(QW(K)+QW(K-1))
     $ -HALF*ABS(WSTAR(K))*(QW(K)-QW(K-1))
      ENDDO
C
      FW(1)=ZERO
      FW(NZ)=ZERO
C
      DO K=1,NL
	QW(K)=QW(K)-(FW(K+1)-FW(K))*DTDSG(K)
      ENDDO
C
      DO K=1,NL
	Q(I,J,K)=QW(K)
      ENDDO
C
      ENDDO


      IF(J.EQ.1 .and. endlonex .EQ. NX_LON) THEN
        DO K = 1 ,NL
          north_pole_q_nx(K)= Q(NX_LON,J,K)
        ENDDO
      ENDIF

      IF(J.EQ.60 .and. beglonex .EQ. 1) THEN
        DO K = 1 ,NL
          south_pole_q_1(K)= Q(1,J,K)
        ENDDO
      ENDIF
C
      ENDDO

      call broadcast_lon_data(NX_LON,1,north_pole_q_nx,nl)
      call broadcast_lon_data(1,plat,south_pole_q_1,nl)

      IF(jbeg0.EQ.1) THEN
        J = 1
!$OMP PARALLEL DO PRIVATE (I,K)
        DO K = 1 ,NL
          DO I = beglonex ,endlonex
            Q(I,J,K)= north_pole_q_nx(K)
          ENDDO
        ENDDO
      ENDIF

      IF(jend0.EQ.60) THEN
        J = jend0
!$OMP PARALLEL DO PRIVATE (I,K)
        DO K = 1 ,NL
          DO I = beglonex ,endlonex
            Q(I,J,K)= south_pole_q_1(K)
          ENDDO
        ENDDO
      ENDIF


      RETURN
      END

