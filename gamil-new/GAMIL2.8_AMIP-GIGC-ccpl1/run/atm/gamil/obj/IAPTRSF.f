# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/IAPTRSF.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/IAPTRSF.F"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/IAPTRSF.F" 2

# 1 "./params.h" 1
# 15 "./params.h"
 
 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/IAPTRSF.F" 2

      SUBROUTINE IAPTRSF(US,VS,WS,PES,T,UU,VV,P,TT,TB2,CB,KP)
!
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
# 12 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/IAPTRSF.F" 2

!
!     The file PARA is to define the parameters related to the model resolution:
!     NX is the grid number in longitude
!     NY is the grid number in latitude
!     NL is the number of vertical layers
!
	REAL*8
     _	     US(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL )  !  US = UU, zonal wind,    input variable
     _      ,VS(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL )  !  VS = VV, meridional wind, input variable
     _      ,WS(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NZ )  !  vertical velocity, input variable
     _      ,PES (ilbnd:ihbnd,beglatexdyn:endlatexdyn)  !  Surface pressure, input variable
     _      ,T (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  temperature, input variable
     _      ,TB2 (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  the temperature of the standard atmosphere a
!                           !  at the sigma layers, input constant
     _      ,CB (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  CB=Sqrt(R*(KK*TBB-dTB/dlnp)), input variable,
!                           !  where, KK=R/Cp, R is a constant
     _      ,UU  (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  UU=UU*sqrt(PES), output variable
     _      ,VV  (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !	 VV=VV*sqrt(PES), output variable
     _      ,P  (ilbnd:ihbnd,beglatexdyn:endlatexdyn   )  !  P=sqrt(PES)  , output variable
     _      ,TT (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  TT=R*T'*Sqrt(PES)/CB, output variable
!                           !  where T'=T-TB2, T is the temperatur,
     _      ,PX,PY          !  working variables
!
	INTEGER I,J,K,KP


!$OMP PARALLEL DO PRIVATE (I, J)
        DO J=jbeg0,jend0
	  DO I=ibeg1,iend1
	     P(I,J)=SQRT(PES(I,J))
	  END DO
        END DO

      call gamil_arrays_comm(COMM_TO_LEFT,1,p(:,beglatexdyn),pes(:,beglatexdyn))
      call gamil_arrays_comm(COMM_TO_RIGHT,1,p(:,beglatexdyn),pes(:,beglatexdyn))
      call gamil_arrays_comm(COMM_TO_TOP,1,p(:,beglatexdyn))


      IF (KP.EQ.1) THEN

!$OMP PARALLEL DO PRIVATE (I, J)
        DO J=jbeg0,jend0
	  DO I=beglonex,endlonex
             WS(I,J,1)=0.0D0
             WS(I,J,NZ)=0.0D0
	  END DO
	END DO
      ENDIF
        

!$OMP PARALLEL DO PRIVATE (I, J, K, PY, PX)
      DO J=jbeg0,jend0
        DO I=ibeg1,iend1
!
!         CALCULATING THE AVERAGE VALUE OF P AT VV-GRID.
!
          IF(J.LT.60) THEN
            PY=0.5D0*(P(I,J)+P(I,J+1))
          ELSE
            PY=0.0D0
          ENDIF
!
!         CALCULATING THE AVERAGE VALUE OF P AT UU-GRID.
!
          PX=0.5D0*(P(I,J)+P(I-1,J))
!
          DO K=1,NL
!
             UU (I,J,K)=PX*US(I,J,K)
             VV (I,J,K)=PY*VS(I,J,K)
!
!         CALCULATING TT=(T-TB2)*P*RD/CB.
!
	     TT (I,J,K)=(T(I,J,K)-TB2(I,J,K))*P(I,J)*RD/CB(I,J,K)
!
          ENDDO
        ENDDO
      ENDDO

      call gamil_arrays_comm(COMM_TO_LEFT,1,US(:,beglatexdyn,1),VS(:,beglatexdyn,1),
     &                       T(:,beglatexdyn,1),UU(:,beglatexdyn,1),VV(:,beglatexdyn,1),TT(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_RIGHT,1,US(:,beglatexdyn,1),VS(:,beglatexdyn,1),
     &                       T(:,beglatexdyn,1),UU(:,beglatexdyn,1),VV(:,beglatexdyn,1),TT(:,beglatexdyn,1))

!
	RETURN
	END
