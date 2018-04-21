# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DIFPS1.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DIFPS1.F"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DIFPS1.F" 2

# 1 "./params.h" 1
# 15 "./params.h"
 
 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DIFPS1.F" 2

!!(2003.11.10)
!***********************************************************************
!
      SUBROUTINE DIFPS1(U,V,P,PS,WS,DPS,DPSU
     _                ,DSIG,DY,OUX,OUY,SINV,MP1,MP2,MM1,MM2,WTGV)
     _


      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil
      
	IMPLICIT NONE
!
!	This subroutine is to calculate the tendency of the surface prressure DPS,
!     the vertical velocity WS, the zonal wind US, meridional wind VS and the
!	departure of the geopotential height from the standard atmopshere


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
# 22 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DIFPS1.F" 2

!
!	The file PARA is to define the parameters related to the model resolution:
!     NX_LON is the grid number in longitude
!     NY is the grid number in latitude
!     NL is the number of vertical layers
!
	REAL*8
     _       U  (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  U=u*sqrt(Ps), input variable
     _      ,V  (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !	 V=v*sqrt(Ps), input variable
     _      ,P  (ilbnd:ihbnd,beglatexdyn:endlatexdyn   )  !  P=sqrt(Ps)  , input variable
     _      ,PS (ilbnd:ihbnd,beglatexdyn:endlatexdyn   )  !	 Surface pressure, input variable
     _      ,DSIG(NL     )  !  The vertical stepsizes, input constant
     _      ,DY             !  The horizontal stepsize in meridional direction
!                              input constant
     _      ,SINV(beglatexdyn:endlatexdyn)       !  sin(theta) at half grid j+1/2, input constant
     _      ,OUX(beglatexdyn:endlatexdyn)        !  OUX=1/(RAD*SINU*DX*MDJ), input constant
!                              where, DX is the horizontal stepsize in zonal direction,
!                              MDJ is the leaping length of the central difference
!                              SINU is sin(theta) at intger grid j
     _      ,OUY(beglatexdyn:endlatexdyn)        !  OUY=1/(RAD*SINU*DY*WTGU), input constant,
!                              where WTGU is the weighting at the integer grid j
     _      ,WTGV(beglatexdyn:endlatexdyn)       !	 the weighting at the half grid j+1/2,
!                              input constant
     _      ,WS  (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NZ) !  WS = w, vertical velocity, output variable
     _      ,DPS (ilbnd:ihbnd,beglatexdyn:endlatexdyn   ) !  the tendency of the surface pressure,
     _      ,DPSU(ilbnd:ihbnd,beglatexdyn:endlatexdyn   ) !  the tendency of the surface pressure,
!                              output variable
      REAL*8 WK0,WK1,WK2,WK3 !  working variables
     _      ,PXP,PXM,PYP	  !  working variables
     _      ,PYM,DPSP,WKQ   !  working variables
      REAL*8 W1(NZ),D1_8(NZ), TMP_DATA
      REAL*16 D1_16(NZ),D1_16_TMP(8,NZ),TMP_SUM  !  working variables
      REAL*8 pole_ps_1(1)
!
	INTEGER
     _       MM1(ilbnd:ihbnd,beglatexdyn:endlatexdyn)     !
     _      ,MP1(ilbnd:ihbnd,beglatexdyn:endlatexdyn)     !
     _      ,MM2(ilbnd:ihbnd,beglatexdyn:endlatexdyn)     !
     _      ,MP2(ilbnd:ihbnd,beglatexdyn:endlatexdyn)     !
!
      INTEGER I,J,K,IP1,IM1,IP2,IM2 !  working variables
      INTEGER J1,J2
      INTEGER icomm_request1, icomm_request2

      call gamil_arrays_comm(COMM_TO_TOP,1,p(:,beglatexdyn),request_id=icomm_request1)
      call gamil_arrays_comm(COMM_TO_BOT,1,p(:,beglatexdyn),v(:,beglatexdyn,1),request_id=icomm_request2)
      call wait_icomm_request(icomm_request1)
      call wait_icomm_request(icomm_request2)

C     DEALING WITH INTERNAL GRID POINTS.
C     CALCULATING FACTORS RELATIVE TO J IN DIVERGENCE FORMULA.
C

      call t_startf("DIFPS1 COMP")
!$OMP PARALLEL DO PRIVATE (I,J,K,WK1,WK2,WK3,IP1,IM1,IP2,IM2,PXP,PXM,
!$   &           PYP,PYM,DPSP,D1_8,W1,WK0,WKQ,J1)
      DO J=jbeg1,jend1
        DO I=ibeg1,iend1
          WS(I,J,1)=0.0D0
          WS(I,J,NZ)=0.0D0
        ENDDO
        
          WK1=0.5D0*OUX(J)
          WK2=0.5D0*OUY(J)*SINV(J)
          WK3=0.5D0*OUY(J)*SINV(J-1)
C     CALCULATING FACTORS RELATIVE TO I AND J IN DIVERGENCE FORMULA.
C
          DO I=ibeg1,iend1
            IP1=MP1(I,J)
            IM1=MM1(I,J)
            IP2=MP2(I,J)
            IM2=MM2(I,J)
            PXP=WK1*(P(IP1,J)+P(I,J))
            PXM=WK1*(P(I,J)+P(IM1,J))
            PYP=WK2*(P(I,J+1)+P(I,J))
            PYM=WK3*(P(I,J)+P(I,J-1))
C     CALCULATING DIVERGENCES AS WELL AS SUM OF THEM.
            DPSP=0.0D0
            DO K=1,NL
              D1_8(K)=PYP*V(I,J,K)-PYM*V(I,J-1,K)
              DPSP=DPSP-DSIG(K)*D1_8(K)
              D1_8(K)=D1_8(K)+PXP*U(IP2,J,K)-PXM*U(IM2,J,K)
            ENDDO
C
C     CALCULATING DPS/DT, DPS/DT AND D(SIGMA)/DT.
C
            DPS(I,J)=DPSP
	    DPSP=DPSP+DPSU(I,J)
            WKQ=1.D0/PS(I,J)
            DO K=2,NL
              WS(I,J,K)=WS(I,J,K-1)-DSIG(K-1)*WKQ*(DPSP+D1_8(K-1))
            ENDDO
          ENDDO
        ENDDO
      call t_stopf("DIFPS1 COMP")

      DO J=jbeg0,jend0 
      IF (J.EQ.1 .OR. J.EQ.60) THEN
          call t_startf("DIFPS1 COMP")
          DO I=ibeg1,iend1
             WS(I,J,1)=0.0D0
             WS(I,J,NZ)=0.0D0
          ENDDO
 
          IF(J.EQ.1) THEN
             J1 = J+1
             J2 = J
             WK0=2.0D0/(DFLOAT(NX_LON-2)*RAD*DY)*WTGV(J2)
          ELSE
             J1 = J-1
             J2 = J1
             WK0=-2.0D0/(DFLOAT(NX_LON-2)*RAD*DY)*WTGV(J2)
          ENDIF
C
C     FINDING DP/DT AND D(SIGMA)/DT AT POLES.
C     IN BELOW, SUBSCRIPTS 1 AND 2 REPRESENT J=1 AND 60 RESPECTIVELY.

          W1(1)=0.0D0
C     CALCULATING DIVERGENCE AT POLES.
!$OMP PARALLEL DO PRIVATE (I,K,TMP_SUM,TMP_DATA)
          DO K=1,NL
            TMP_SUM=0.0D0
            DO I=ibeg1,iend1
              TMP_DATA=(P(I,J)+P(I,J1))*V(I,J2,K)
              TMP_SUM=TMP_SUM+TMP_DATA
            ENDDO
            D1_16_TMP(1,K)=WK0*TMP_SUM
          ENDDO
          call t_stopf("DIFPS1 COMP")

          DO K=1,NL
            D1_16(K)=D1_16_TMP(1,K)
          ENDDO

          call gamil_sum_pole_data_phys(J,D1_16,NL)

          call t_startf("DIFPS1 COMP")
          DPSP=0.0D0
          DO K=1,NL
            D1_8(K) = D1_16(K)
            DPSP=DPSP-DSIG(K)*D1_8(K)
          ENDDO
C
C     CALCULATING DPS/DT AND D(SIGMA)/DT AT POLES.
C
          IF (ibeg0 .eq. 1) pole_ps_1 = PS(1,J)
          call broadcast_lon_data(1,j,pole_ps_1,1)
          DO K=2,NL
            W1(K)=W1(K-1)-DSIG(K-1)/pole_ps_1(1)*(DPSP+D1_8(K-1))
          ENDDO
C
          DO I=ibeg1,iend1
            DPS(I,J)=DPSP
          ENDDO

!$OMP PARALLEL DO PRIVATE (I,K)
          DO K=2,NL
            DO I=ibeg1,iend1
              WS(I,J,K)=W1(K)
            ENDDO
          ENDDO

        call t_stopf("DIFPS1 COMP")
        ENDIF
      ENDDO

      call gamil_arrays_comm(COMM_TO_LEFT,1,DPS(:,beglatexdyn),WS(:,beglatexdyn,1),request_id=icomm_request1)
      call gamil_arrays_comm(COMM_TO_RIGHT,1,DPS(:,beglatexdyn),WS(:,beglatexdyn,1),request_id=icomm_request2)
      call wait_icomm_request(icomm_request1)
      call wait_icomm_request(icomm_request2)

      RETURN
      END
