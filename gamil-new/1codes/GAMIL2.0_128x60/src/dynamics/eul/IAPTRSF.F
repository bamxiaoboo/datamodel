#include <misc.h>
#include <params.h>

      SUBROUTINE IAPTRSF(US,VS,WS,PES,T,UU,VV,P,TT,TB2,CB,KP)
!
      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil

      IMPLICIT NONE

#include <PARADYN>

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
          IF(J.LT.PLAT) THEN
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
