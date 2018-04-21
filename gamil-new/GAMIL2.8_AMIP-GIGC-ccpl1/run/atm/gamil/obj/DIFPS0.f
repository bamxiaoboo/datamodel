# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DIFPS0.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DIFPS0.F"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DIFPS0.F" 2

# 1 "./params.h" 1
# 15 "./params.h"
 
 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DIFPS0.F" 2

!!(2003.11.29)
!***********************************************************************
!
        SUBROUTINE DIFPS0(UU,P,DPS0,DSIG,OUX,MP1,MP2,MM1,MM2)
        
      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil 
      
	IMPLICIT NONE
!
!	This subroutine is to calculate the tendency of the surface prressure DPS0,
!     the vertical velocity WS, the zonal wind US, meridional wind VS and the
!	departure of the geopotential height from the standard atmopshere
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
# 19 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DIFPS0.F" 2

!
!	The file PARA is to define the parameters related to the model resolution:
!     NX is the grid number in longitude
!     NY is the grid number in latitude
!     NL is the number of vertical layers
!
	REAL*8
     _       UU  (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  UU=UU*sqrt(Ps), input variable
     _      ,P  (ilbnd:ihbnd,beglatexdyn:endlatexdyn)  !  P=sqrt(Ps)  , input variable
     _      ,DSIG(NL)  !  The vertical stepsizes, input constant
     _      ,OUX(beglatexdyn:endlatexdyn)        !  OUX=1/(RAD*SINU*DX*MDJ), input constant
!                              where, DX is the horizontal stepsize in zonal direction,
!                              MDJ is the leaping length of the central difference
!                              SINU is sin(theta) at intger grid j
     _      ,DPS0(ilbnd:ihbnd,beglatexdyn:endlatexdyn)  !  the tendency of the surface pressure,
!                              output variable
      REAL*8   WK1,WK2,WK3    !  working variables
     _      ,PXP,PXM,DPSP   !  working variables
!
	INTEGER
     _       MM1(ilbnd:ihbnd,beglatexdyn:endlatexdyn)     !
     _      ,MP1(ilbnd:ihbnd,beglatexdyn:endlatexdyn)     !
     _      ,MM2(ilbnd:ihbnd,beglatexdyn:endlatexdyn)     !
     _      ,MP2(ilbnd:ihbnd,beglatexdyn:endlatexdyn)     !
!
      INTEGER I,J,K,IP1,IM1,IP2,IM2 !  working variables


C     DEALING WITH INTERNAL GRID POINTS.
C     CALCULATING FACTORS RELATIVE TO J IN DIVERGENCE FORMULA.
C

!$OMP PARALLEL DO PRIVATE (I,J,K,WK1,IP1,IM1,IP2,IM2,PXP,PXM,DPSP)
      do j=jbeg0,jend0
        IF(J.GE.2.AND.J.LE.60-1) THEN
          WK1=0.5D0*OUX(J)
C     CALCULATING FACTORS RELATIVE TO I AND J IN DIVERGENCE FORMULA.
C
          DO I=ibeg1,iend1
            IP1=MP1(I,J)
            IM1=MM1(I,J)
            IP2=MP2(I,J)
            IM2=MM2(I,J)
            PXP=WK1*(P(IP1,J)+P(I,J))
            PXM=WK1*(P(I,J)+P(IM1,J))
C     CALCULATING DIVERGENCES AS WELL AS SUM OF THEM.
            DPSP=0.0D0
            DO K=1,NL
              DPSP=DPSP-DSIG(K)*(PXP*UU(IP2,J,K)-PXM*UU(IM2,J,K))
            ENDDO
C
C     CALCULATING DPS0/DT, DPS0/DT AND D(SIGMA)/DT.
C
            DPS0(I,J)=DPSP
          ENDDO
        ELSE IF(J.EQ.1) THEN
C
C     FINDING DP/DT AND D(SIGMA)/DT AT POLES.
C     IN BELOW, SUBSCRIPTS 1 AND 2 REPRESENT J=1 AND 60 RESPECTIVELY.

          DO I=ibeg1,iend1
            DPS0(I,J)=0.0
          ENDDO
          
       ELSE IF(J.EQ.60) THEN
C
C     FINDING DP/DT AND D(SIGMA)/DT AT POLES.
C     IN BELOW, SUBSCRIPTS 1 AND 2 REPRESENT J=1 AND 60 RESPECTIVELY.
          DO I=ibeg1,iend1
            DPS0(I,J)=0.0
          ENDDO
       ENDIF
       END DO

      call gamil_arrays_comm(COMM_TO_LEFT,1,DPS0(:,beglatexdyn))
      call gamil_arrays_comm(COMM_TO_RIGHT,1,DPS0(:,beglatexdyn))

      RETURN
      END
