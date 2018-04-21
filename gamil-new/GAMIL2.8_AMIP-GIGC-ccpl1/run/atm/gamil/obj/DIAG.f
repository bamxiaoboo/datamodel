# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DIAG.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DIAG.F"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DIAG.F" 2

# 1 "./params.h" 1
# 15 "./params.h"
 
 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DIAG.F" 2

!!(2003.11.29)
!!-------------------

	SUBROUTINE DIAG(UU,VV,P,PES,PLY2,TT,US,VS,TS,GHI,HPS
     _               ,PMTOP,PSB,TSB,TB2,CB,FAC,DSIG)
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
# 15 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DIAG.F" 2

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
          IF (J .eq. 60) THEN
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
