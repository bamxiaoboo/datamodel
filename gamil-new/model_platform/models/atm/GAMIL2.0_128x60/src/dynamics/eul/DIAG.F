#include <misc.h>
#include <params.h>

!!(2003.11.29)
!!-------------------

	SUBROUTINE DIAG(UU,VV,P,PES,PLY2,TT,US,VS,TS,GHI,HPS
     _               ,PMTOP,PSB,TSB,TB2,CB,FAC,DSIG)
!
      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil
	IMPLICIT NONE

#include <PARADYN>

!
!	The file PARA is to define the parameters related to the model resolution:
!     NX is the grid number in longitude
!     NY is the grid number in latitude
!     NL is the number of vertical layers
!
	REAL*8
     _     UU  (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  UU=UU*sqrt(PES), input variable
     _    ,VV  (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !	 VV=VV*sqrt(PES), input variable
     _    ,P  (ilbnd:ihbnd,beglatexdyn:endlatexdyn   )  !  P=sqrt(PES)  , input variable
     _    ,PES (ilbnd:ihbnd,beglatexdyn:endlatexdyn   )  !	 Surface pressure, input variable
     _    ,PLY2(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  PLY2=p, Pressure in Sigma Layer, input variable
     _    ,TS (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  input variable, TS=T, TEMPERATURE
     _    ,TT (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  TT=R*T'*Sqrt(PES)/CB, input variable
!                            where T'=T-TB2, T is the temperatur,
!                            TBB is Temperature of the standard atmosphere
     _    ,TB2 (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  input variable,
     _    ,CB (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  CB=Sqrt(R*(KK*TBB-dTB/dlnp)), input variable,
!                            where, KK=R/Cp, R is a constant
     _    ,PMTOP          !  PMTOP=10hPa
     _    ,TSB(ilbnd:ihbnd,beglatexdyn:endlatexdyn   )  !  TBB at the surface, input constant
     _    ,PSB(ilbnd:ihbnd,beglatexdyn:endlatexdyn   )  !  PSB is the surface pressure of the standard
!					 atmosphere, input constant
     _    ,DSIG(NL     )  !  The vertical stepsizes, input constant
     _	  ,US(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL )  !  US = UU, zonal wind,    output variable
     _    ,VS(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL )  !  VS = VV, meridional wind,output variable
     _    ,HPS(ilbnd:ihbnd,beglatexdyn:endlatexdyn   )  !  the surface geopotential height deviation
     _	  ,GHI (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NZ )  !  GHI=gz-HBB, gz is the geopotential height,
!                              HBB is the geopotential height of the standard atmopshere
     _    ,FAC(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NZ)
     _    ,WK1,WK2(ilbnd:ihbnd,beglatexdyn:endlatexdyn)        !  working variables
!
	INTEGER I,J,K

      call t_startf("DIAG COMP")

!$OMP PARALLEL DO PRIVATE(I,J,K,WK1)
      DO J=jbeg0,jend0
        DO I=ibeg1,iend1
	      HPS(I,J   )=RD*TSB(I,J)/PSB(I,J)*(PES(I,J)+PMTOP-PSB(I,J))
          GHI  (I,J,NZ)=HPS(I,J)
          IF (J .eq. PLAT) THEN
            WK2(I,J)=0.0D0
          ELSE
            WK2(I,J)=2.0D0/(P(I,J)+P(I,J+1))
          ENDIF
        ENDDO

        DO K=NL,1,-1
          DO I=ibeg1,iend1
            WK1=2.0D0/(P(I,J)+P(I-1,J))
            US(I,J,K)=WK1*UU(I,J,K)
            VS(I,J,K)=WK2(I,J)*VV(I,J,K)
            TS(I,J,K)=TT(I,J,K)*CB(I,J,K)/(P(I,J)*RD)+TB2(I,J,K)
            GHI(I,J,K)=GHI(I,J,K+1)+DSIG(K)*P(I,J)*CB(I,J,K)
     &       /PLY2(I,J,K)*TT(I,J,K)*.5*(FAC(I,J,K+1)+FAC(I,J,K))
          ENDDO
        ENDDO
      ENDDO

      call t_stopf("DIAG COMP")

      call gamil_arrays_comm(COMM_TO_LEFT,1,US(:,beglatexdyn,1),TS(:,beglatexdyn,1),
     &                       VS(:,beglatexdyn,1),HPS(:,beglatexdyn),GHI(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_RIGHT,1,US(:,beglatexdyn,1),TS(:,beglatexdyn,1),
     &                       VS(:,beglatexdyn,1),HPS(:,beglatexdyn),GHI(:,beglatexdyn,1))


	RETURN
	END
