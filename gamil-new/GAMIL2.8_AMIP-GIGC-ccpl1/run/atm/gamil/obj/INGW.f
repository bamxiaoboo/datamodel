# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/INGW.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/INGW.F"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/INGW.F" 2

# 1 "./params.h" 1
# 15 "./params.h"
 
 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/INGW.F" 2

      SUBROUTINE INGW(UK,P,TTK,CB0,DSIG,OUX,DT,SINU,WTGU,DU0,DTT0,ITN,EE)
      use mpi_gamil
!
      IMPLICIT NONE
!
!	This subroutine is to calculate the tendency of the wind and the temperature:
!     DU, DV, DTT
!


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
# 14 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/INGW.F" 2

!
!	The file PARA is to define the parameters related to the model resolution:
!     NX is the grid number in longitude
!     NY is the grid number in latitude
!     NL is the number of vertical layers
!
	  REAL*8
     _       U (ilbnd:ihbnd,NL )  !  input variable, zonal wind velocity
     _      ,UP(ilbnd:ihbnd,NL )  !
     _      ,UK(ilbnd:ihbnd,NL )  !
     _      ,P  (ilbnd:ihbnd   )  !  input variable, P  = sqrt(Ps)
     _      ,TT (ilbnd:ihbnd,NL)  !  input variable, TT=R*T'*Sqrt(Ps)/CB,
!                              where T'=T-TB, T is the temperatur,
!                              TBB	is Temperature of the standard atmosphere
     _	  ,TTK(ilbnd:ihbnd,NL )  !  input variable,
     _	  ,TTP(ilbnd:ihbnd,NL )  !  input variable,
!							 H=gz-HBB, gz is the geopotential height,
!                              HBB is the geopotential height of the standard atmopshere
     _      ,CB0(ilbnd:ihbnd,NL)     !  input variable, CB=Sqrt(R*(KK*TBB-dTB/dlnp)),
!                           !  CB0=CB*P/PLY
!                              where, KK=R/Cp, R is a constant
     _      ,DSIG(NL     )  !  input constant, the vertical stepsizes
     _      ,OUX            !  input constant, OUX=1/(RAD*SINU*DX*MDJ)
!                              where, DX is the horizontal stepsize in zonal direction,
!                              MDJ is the leaping length of the central difference
!                              SINU is sin(theta) at intger grid j
!
      INTEGER I,K,KWB,JWB,ITN
!
      REAL*8
     _       DU0 (ilbnd:ihbnd,NL)	!  input variables
     _      ,DTT0(ilbnd:ihbnd,NL)  !  input variables
     _      ,DU1 (ilbnd:ihbnd,NL)	!  output variables
     _      ,DTT1(ilbnd:ihbnd,NL)  !  output variables
     _      ,DU2 (ilbnd:ihbnd,NL)	!  output variables
     _      ,DTT2(ilbnd:ihbnd,NL)  !  output variables
     _      ,DU  (ilbnd:ihbnd,NL)	!  output variables
     _      ,DTT (ilbnd:ihbnd,NL)  !  output variables
     _      ,DUX (ilbnd:ihbnd,NL)	!  output variables
     _      ,DTTX(ilbnd:ihbnd,NL)  !  output variables
!
	REAL*8 DTGW,DT2,DT,DJ,DS,SINU,WTGU
        REAL*16 Y1,Y2,Y3,EE,Y1_ARR(NL),Y2_ARR(NL),Y3_ARR(NL)
!
	DTGW=DT/FLOAT(ITN)
	DT2=0.5*DTGW
!
!$OMP PARALLEL DO PRIVATE (I,K)
	DO K=1,NL
	DO I=beglonex,endlonex
	   U (I,K)=UK (I,K)
	   TT(I,K)=TTK(I,K)
!
	   DU (I,K)=0.0
	   DTT(I,K)=0.0
!
	   DU1 (I,K)=0.0
	   DTT1(I,K)=0.0
!
	   DU2 (I,K)=0.0
	   DTT2(I,K)=0.0
	ENDDO
	ENDDO
!
	DO KWB=1,ITN
!
!$OMP PARALLEL DO PRIVATE (I,K)
	DO K=1,NL
	DO I=beglonex,endlonex
	   UP (I,K)=U (I,K)
	   TTP(I,K)=TT(I,K)
	ENDDO
	ENDDO
!
	DO JWB=1,2
!
	CALL DIFUTX(U,P,TT,CB0,DSIG,OUX,DUX,DTTX)
!
!$OMP PARALLEL DO PRIVATE (I,K)
	DO K=1,NL
	DO I=beglonex,endlonex
	   U (I,K)=UP (I,K)+DT2*(DUX (I,K)+DU0 (I,K))
	   TT(I,K)=TTP(I,K)+DT2*(DTTX(I,K)+DTT0(I,K))
	ENDDO
	ENDDO
!
	END DO
!
!$OMP PARALLEL DO PRIVATE (I,K)
	DO K=1,NL
	DO I=beglonex,endlonex
	   DU1 (I,K)=DU1 (I,K)+DUX (I,K)
	   DTT1(I,K)=DTT1(I,K)+DTTX(I,K)
	ENDDO
	ENDDO
!
	CALL DIFUTX(U,P,TT,CB0,DSIG,OUX,DUX,DTTX)
!
!$OMP PARALLEL DO PRIVATE (I,K)
	DO K=1,NL
	DO I=beglonex,endlonex
	   U (I,K)=UP (I,K)+DTGW*(DUX (I,K)+DU0 (I,K))
	   TT(I,K)=TTP(I,K)+DTGW*(DTTX(I,K)+DTT0(I,K))
!
	   DU2 (I,K)=DU2 (I,K)+DU (I,K)
	   DTT2(I,K)=DTT2(I,K)+DTT(I,K)
!
	   DU  (I,K)=DU  (I,K)+DUX (I,K)
	   DTT (I,K)=DTT (I,K)+DTTX(I,K)
	ENDDO
	ENDDO
!
	END DO
!
        DJ=SINU/WTGU
!$OMP PARALLEL DO PRIVATE (I,K,Y1,Y2,Y3,DS)
	DO K=1,NL
	   Y1=0.0
	   Y2=0.0
	   Y3=0.0
	   DS=DJ*DSIG(K)
	   DO I=ibeg1,iend1
	      DU  (I,K)=DU  (I,K)/FLOAT(ITN)
	      DTT (I,K)=DTT (I,K)/FLOAT(ITN)
!
		  Y1=Y1+(DU (I,K)*DU2 (I,K)+DTT(I,K)*DTT2(I,K))*DS
		  Y2=Y2+(DU (I,K)*DU1 (I,K)+DTT(I,K)*DTT1(I,K))*DS
		  Y3=Y3+(DU (I,K)*DU0 (I,K)+DTT(I,K)*DTT0(I,K))*DS
!
	      DU0 (I,K)=DU0 (I,K)+DU (I,K)
	      DTT0(I,K)=DTT0(I,K)+DTT(I,K)
	   ENDDO
           Y1_ARR(K)=Y1
           Y2_ARR(K)=Y2
           Y3_ARR(K)=Y3
	ENDDO

       call register_comm_array(ilbnd,ihbnd,1,1,1,NL,1,1,DU0(:,1)) 
       call register_comm_array(ilbnd,ihbnd,1,1,1,NL,1,1,DTT0(:,1)) 
       call gamil_arrays_comm(COMM_TO_LEFT,1,DU0(:,1),DTT0(:,1))
       call gamil_arrays_comm(COMM_TO_RIGHT,1,DU0(:,1),DTT0(:,1))
       call remove_comm_array(DTT0(:,1)) 
       call remove_comm_array(DU0(:,1)) 

	Y1=0.0
	Y2=0.0
	Y3=0.0
	DO K=1,NL
           Y1 = Y1+Y1_ARR(K)
           Y2 = Y2+Y2_ARR(K)
           Y3 = Y3+Y3_ARR(K)
	ENDDO

	EE=(Y1+Y1+Y2)/FLOAT(ITN*ITN)+Y3
!
	RETURN
	END
!
!!!#include <misc.h>
!!!#include <params.h>

!!(2004.02.15)
!***********************************************************************
!
      SUBROUTINE MINUS_INGW(UU,P,TT,CB0,DSIG,OUX,DU,DTT,NIGW)
! 
      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil
!
      IMPLICIT NONE
!
!	This subroutine is to deduct the inner gravity waves from the tendency
!     of the zonal wind and the temperature: DU, DTT
!


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
# 191 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/INGW.F" 2

!
!	The file PARA is to define the parameters related to the model resolution:
!     NX is the grid number in longitude
!     NY is the grid number in latitude
!     NL is the number of vertical layers
!
	  REAL*8
     _       UU (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL )  !  input array, zonal wind velocity
     _      ,UJ (ilbnd:ihbnd,   NL) !  working array,
     _      ,P  (ilbnd:ihbnd,beglatexdyn:endlatexdyn)  !  input array, P  = sqrt(Ps)
     _      ,TT (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  input array, TT=R*T'*Sqrt(Ps)/CB,
     _      ,TJ (ilbnd:ihbnd,    NL)  !  input array, TT=R*T'*Sqrt(Ps)/CB,
!                              where T'=T-TB, T is the temperatur,
!                              TBB	is Temperature of the standard atmosphere
     _      ,CBJ(ilbnd:ihbnd   ,NL)  !  input array, CB=Sqrt(R*(KK*TBB-dTB/dlnp)),
     _      ,CB0(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  input array, CB=Sqrt(R*(KK*TBB-dTB/dlnp)),
!                           !  CB0=CB*P/PLY
!                              where, KK=R/Cp, R is a constant
     _      ,DSIG(NL     )  !  input constant array, the vertical stepsizes
     _      ,OUX(beglatexdyn:endlatexdyn)        !  input constant array, OUX=1/(RAD*SINU*DX*MDJ)
!                              where, DX is the horizontal stepsize in zonal direction,
!                              MDJ is the leaping length of the central difference
!                              SINU is sin(theta) at intger grid j
!
      REAL*8
     _       DU (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)   !  input & output array
     _      ,DTT (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  input & output array
     _      ,DUX (ilbnd:ihbnd   ,NL)  !  working array: inner gravity waves in DU
     _      ,DTTX(ilbnd:ihbnd   ,NL)  !  working array: inner gravity waves in DTT
!
      INTEGER
     _       I               !  working variable
     _      ,J			   !  working variable
     _      ,K			   !  working variable
     _      ,NIGW(beglatexdyn:endlatexdyn)		   !  times of moving length of inner gravity waves to
!						   !  the zonal gridsize
!

       call gamil_arrays_comm(COMM_TO_LEFT,1,p(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_RIGHT,1,p(:,beglatexdyn))

       DO j=jbeg1, jend1
	IF (NIGW(J).GT.1) THEN
!
!$OMP PARALLEL DO PRIVATE (I,K)
	   DO K=1,NL
	   DO I=beglonex,endlonex
	      UJ (I,K)=UU  (I,J,K)
	      TJ (I,K)=TT (I,J,K)
	      CBJ(I,K)=CB0(I,J,K)
	   END DO
	   END DO
!
	   CALL DIFUTX(UJ,P(ilbnd,J),TJ,CBJ,DSIG,OUX(J),DUX,DTTX)
!
!$OMP PARALLEL DO PRIVATE (I,K)
	   DO K=1,NL
	   DO I=beglonex,endlonex
	      DU (I,J,K)=DU (I,J,K)-DUX (I,K)
	      DTT(I,J,K)=DTT(I,J,K)-DTTX(I,K)
	   END DO
	   END DO
!
	ENDIF
!
	END DO
!
	RETURN
	END
!
!!#include <misc.h>
!!#include <params.h>

!!(wb 2004.02.15)
!***********************************************************************
!
      SUBROUTINE PLUS_INGW(UK,P,TTK,CB0,DSIG,OUX,DT,SINU,WTGU,DU0,DTT0,EE,NIGW)
!
      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil
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
# 276 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/INGW.F" 2
!
!	This subroutine is to calculate the tendency of the wind and the temperature:
!     DU, DV, DTT
!
!	The file PARA is to define the parameters related to the model resolution:
!     NX is the grid number in longitude
!     NY is the grid number in latitude
!     NL is the number of vertical layers
!
	  REAL*8
     _       UK (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL )  !
     _      ,UKJ(ilbnd:ihbnd   ,NL )  !
     _      ,P  (ilbnd:ihbnd,beglatexdyn:endlatexdyn   )  !  input variable, P  = sqrt(Ps)
     _	    ,TTK(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL )  !  input variable,
     _	    ,TKJ(ilbnd:ihbnd   ,NL )  !  input variable,
     _      ,CBJ(ilbnd:ihbnd   ,NL)  !  input variable, CB=Sqrt(R*(KK*TBB-dTB/dlnp)),
     _      ,CB0(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  input variable, CB=Sqrt(R*(KK*TBB-dTB/dlnp)),
!                           !  CB0=CB*P/PLY
!                              where, KK=R/Cp, R is a constant
     _      ,DSIG(NL     )  !  input constant, the vertical stepsizes
     _      ,OUX(beglatexdyn:endlatexdyn)        !  input constant, OUX=1/(RAD*SINU*DX*MDJ)
!                              where, DX is the horizontal stepsize in zonal direction,
!                              MDJ is the leaping length of the central difference
!                              SINU is sin(theta) at intger grid j
!
      REAL*8
     _       DU0 (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)	  !  input variables
     _      ,DTT0(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)     !  input variables
     _      ,DUJ (ilbnd:ihbnd   ,NL)	  !  working variables
     _      ,DTJ (ilbnd:ihbnd   ,NL)     !  working variables
      real*16 EK(60)
!
      REAL*8 DT,EE,SINU(beglatexdyn:endlatexdyn),WTGU(beglatexdyn:endlatexdyn)
      REAL*16 EJ
!
      INTEGER I,J,K,MM,NIGW(beglatexdyn:endlatexdyn)
!
!
       call gamil_arrays_comm(COMM_TO_LEFT,1,p(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_RIGHT,1,p(:,beglatexdyn))

      DO J=jbeg1,jend1
!
      MM=NIGW(J)
      IF (MM.GT.1) THEN
!
!$OMP PARALLEL DO PRIVATE (I,K)
      DO K=1,NL
        DO I=beglonex,endlonex
          UKJ(I,K)=UK (I,J,K)
          TKJ(I,K)=TTK(I,J,K)
          CBJ(I,K)=CB0(I,J,K)
          DUJ(I,K)=DU0 (I,J,K)
          DTJ(I,K)=DTT0(I,J,K)
        ENDDO
      ENDDO


      CALL INGW(UKJ,P(ilbnd,J),TKJ,CBJ,DSIG,OUX(J),DT
     _         ,SINU(J),WTGU(J),DUJ,DTJ,MM,EK(J))
!
!$OMP PARALLEL DO PRIVATE (I,K)
      DO K=1,NL
        DO I=beglonex,endlonex
          DU0 (I,J,K)=DUJ(I,K)
          DTT0(I,J,K)=DTJ(I,K)
        ENDDO
      ENDDO
!
      END IF
!
      ENDDO
!
      EJ=0.0
      DO J=jbeg1,jend1
      IF (NIGW(J) .GT. 1) THEN
         EJ=EJ+EK(J)
      ENDIF
      ENDDO
      call gamil_all_reduce(EJ)
    
      EE=EE+EJ
!
      RETURN
      END
