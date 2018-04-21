      SUBROUTINE DIFUVTX(US,US1,U,V,TT,OUX,OVX
     _                 ,MP1,MP2,MM1,MM2,DU,DV,DTT,J)
!
      use mpi_gamil
      IMPLICIT NONE
!	This subroutine is to calculate the tendency of the wind and the temperature:
!     DU, DV, DTT
!
#include <PARADYN>
!
!	The file PARA is to define the parameters related to the model resolution:
!     NX_LON is the grid number in longitude
!     NY is the grid number in latitude
!     NL is the number of vertical layers
!
	  REAL*8
     _       U  (NX_LON,NL)  !  input variable, U  = u*sqrt(Ps)
     _	    ,US (NX_LON,NL)  !  input variable, US = u, zonal wind
     _	    ,US1(NX_LON,NL)  !  input variable, US = u, zonal wind
     _      ,V  (NX_LON,NL)  !  input variable, V  = v*sqrt(Ps)
     _      ,TT (NX_LON,NL)  !  input variable, TT=R*T'*Sqrt(Ps)/CB,
!                           where T'=T-TB, T is the temperatur,
!                           TBB	is Temperature of the standard atmosphere
     _      ,OUX         !  input constant, OUX=1/(RAD*SINU*DX*MDJ),where,
!                           DX is the horizontal stepsize in zonal direction,
!                           MDJ is the leaping length of the central difference
!                           SINU is sin(theta) at intger grid j
     _      ,OVX         !  input constant, OUX=1/(RAD*SINV*DX*MDJ)
!
	INTEGER
     _       MM1(NX_LON)     !	 input constant
     _      ,MP1(NX_LON)     !	 input constant
     _      ,MM2(NX_LON)     !	 input constant
     _      ,MP2(NX_LON)     !	 input constant
!
      REAL*8
     _       DU (NX_LON,NL)	  !  output variables
     _      ,DV (NX_LON,NL)   !  output variables
     _      ,DTT(NX_LON,NL)  !  output variables
!
      INTEGER I,J,K
      INTEGER IM1,IP2,IP1,IM2,I1,J0,IM3,IP3
      REAL*8  OUX2,OUX4,OVX4
C
C     AS FOLLOWS, THE IMPROVED LEAP-FROG AND THE REGENERATION OF
C     VARIABLES WILL BE FINISHED IN THE SAME CYCLES.
C
      I1=NX_LON-1
      OUX2=OUX*0.5
      OUX4=OUX*0.25
      OVX4=OVX*0.25
C
C     CALCULATING DU/DT, DV/DT AND DTT/DT.
C
      IF(J.EQ.1) THEN
	 DO K=1,NL
            DO I=2,I1
              IP1=MP1(I)
              IM1=MM1(I)
              IP2=MP2(I)
              IM2=MM2(I)
C
              DTT(I,K)=0.0D0
              DU (I,K)=0.0D0
              DV (I,K)=-OVX4*(US1(IP2,K)*V(IP1,K)
     &                       -US1(IM2,K)*V(IM1,K))
            ENDDO
         ENDDO
      ELSE
         DO K=1,NL
            DO I=2,I1
              IP1=MP1(I)
              IM1=MM1(I)
              IP2=MP2(I)
              IM2=MM2(I)
              DTT(I,K)=-OUX2*(US(IP2,K)*TT(IP1,K)-US(IM2,K)*TT(IM1,K))
              DU (I,K)=-OUX4*((US(I,K)+US(IP1,K))*U(IP1,K)
     &                         -(US(I,K)+US(IM1,K))*U(IM1,K))
              DV (I,K)=-OVX4*((US1(IP2,K)+US(IP2,K))*V(IP1,K)
     &                       -(US1(IM2,K)+US(IM2,K))*V(IM1,K))
            ENDDO
         ENDDO
      ENDIF
!
      DO K=1,NL
         DTT(1,K)=DTT(I1,K)
         DTT(NX_LON,K)=DTT(2,K)
         DU (1,K)=DU (I1,K)
         DU (NX_LON,K)=DU (2,K)
         DV (1,K)=DV (I1,K)
         DV (NX_LON,K)=DV (2,K)
      ENDDO
!
      RETURN
      END
