#include <misc.h>
#include <params.h>

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
#include <PARADYN>

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
        IF(J.GE.2.AND.J.LE.PLAT-1) THEN
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
C     IN BELOW, SUBSCRIPTS 1 AND 2 REPRESENT J=1 AND PLAT RESPECTIVELY.

          DO I=ibeg1,iend1
            DPS0(I,J)=0.0
          ENDDO
          
       ELSE IF(J.EQ.PLAT) THEN
C
C     FINDING DP/DT AND D(SIGMA)/DT AT POLES.
C     IN BELOW, SUBSCRIPTS 1 AND 2 REPRESENT J=1 AND PLAT RESPECTIVELY.
          DO I=ibeg1,iend1
            DPS0(I,J)=0.0
          ENDDO
       ENDIF
       END DO

      call gamil_arrays_comm(COMM_TO_LEFT,1,DPS0(:,beglatexdyn))
      call gamil_arrays_comm(COMM_TO_RIGHT,1,DPS0(:,beglatexdyn))

      RETURN
      END
