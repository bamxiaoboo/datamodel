# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/INNER.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/INNER.F"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/INNER.F" 2

# 1 "./params.h" 1
# 15 "./params.h"
 
 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/INNER.F" 2

!!(wh 2003.11.10)
!!-------------------


	FUNCTION INNER(U1,V1,T1,P1,L1,U2,V2,T2,P2,L2
     _                ,DSIG,TBS,PSB,PMTOP,SINU,SINV,WTGU,WTGV)

      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil
      

        use mpishorthand, only: mpicom


	IMPLICIT NONE
!
!	This subroutine is to define the inner product of the dynamical system
!     of the atmosphere. Here, INNER is the inner product of the vector
!     (U1,V1,T1,P1) and the vector (U2,V2,T2,P2)
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
# 25 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/INNER.F" 2

!
!	The file PARA is to define the parameters related to the model resolution:
!     NX is the grid number in longitude
!     NY is the grid number in latitude
!     NL is the number of vertical layers
!
	REAL*8
     _       U1(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !	input variable
     _      ,V1(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !	input variable
     _      ,T1(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !	input variable
     _      ,P1(ilbnd:ihbnd,beglatexdyn:endlatexdyn)     !	input variable
     _      ,U2(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !	input variable
     _      ,V2(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !	input variable
     _      ,T2(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !	input variable
     _      ,P2(ilbnd:ihbnd,beglatexdyn:endlatexdyn)     !	input variable
     _      ,TBS(ilbnd:ihbnd,beglatexdyn:endlatexdyn)    !	input variable
     _      ,PSB(ilbnd:ihbnd,beglatexdyn:endlatexdyn)    !	input variable
     _      ,SINU(beglatexdyn:endlatexdyn)      !	input variable
     _      ,SINV(beglatexdyn:endlatexdyn)      !	input variable
     _      ,WTGU(beglatexdyn:endlatexdyn)		 !	input variable
     _      ,WTGV(beglatexdyn:endlatexdyn)      !	input variable
     _      ,DSIG(NL)		 !	input variable
     _      ,PMTOP		 !	input variable
     _      ,INNER         !  output variable
!
	INTEGER
     _        L1			 !	input variable
     _       ,L2			 !	input variable
     _       ,I,J,K     !  working variables
!
      REAL*8  TMPU(ilbnd:ihbnd,jbeg0:jend0)	 !  working variables
     _       ,TMPV(ilbnd:ihbnd,jbeg0:jend0)	 !  working variables
     _       ,TMPT(ilbnd:ihbnd,jbeg0:jend0)	 !  working variables
     _       ,TMP_DATA

      REAL*16 EUJ(8,jbeg0:jend0)      !  working variables
     _       ,EVJ(8,jbeg0:jend0)      !  working variables
     _       ,ETJ(8,jbeg0:jend0)      !  working variables
     _       ,EPJ(8,jbeg0:jend0)      !  working variables
     _       ,DS           !  working variables
     _       ,DSU          !  working variables
     _       ,DSV          !  working variables
     _       ,DPS1         !  working variables
     _       ,DPS2         !  working variables

        REAL*16 TMP_INNER

!
!
!$OMP PARALLEL DO PRIVATE (I,J,K,DPS1,DPS2,TMP_DATA)
       DO J=jbeg0,jend0
        DO I = ibeg1,iend1
          TMPU(I,J) = 0.0D0
          TMPV(I,J) = 0.0D0
          TMPT(I,J) = 0.0D0
        ENDDO
        DO K = 1,NL
          DO I = ibeg1,iend1
            TMPU(I,J) = TMPU(I,J)+U1(I,J,K)*U2(I,J,K)*DSIG(K)
            TMPV(I,J) = TMPV(I,J)+V1(I,J,K)*V2(I,J,K)*DSIG(K)
            TMPT(I,J) = TMPT(I,J)+T1(I,J,K)*T2(I,J,K)*DSIG(K)
          ENDDO
        ENDDO
        EUJ(1,J) = 0.0D0
        EVJ(1,J) = 0.0D0
        ETJ(1,J) = 0.0D0
        EPJ(1,J) = 0.0D0
        DO I = ibeg1,iend1
!
	    IF (L1.EQ.1) THEN
	       DPS1 = P1(I,J)+PMTOP-PSB(I,J)
	    ELSE
	       DPS1 = P1(I,J)
	    END IF
!
	    IF (L2.EQ.1) THEN
	       DPS2 = P2(I,J)+PMTOP-PSB(I,J)
	    ELSE
	       DPS2 = P2(I,J)
	    END IF
!
          TMP_DATA = RD*TBS(I,J)/PSB(I,J)*DPS1*DPS2
          EPJ(1,J) = EPJ(1,J) + TMP_DATA
          EUJ(1,J) = EUJ(1,J) + TMPU(I,J)
          EVJ(1,J) = EVJ(1,J) + TMPV(I,J)
          ETJ(1,J) = ETJ(1,J) + TMPT(I,J)
        ENDDO
      ENDDO

      
!
      TMP_INNER=0.0D0
!

      DO J=jbeg0,jend0
        IF(J.EQ.1) THEN
          DS=0.25D0*SINV(1)/WTGV(1)
          DSU=DS
          DSV=4.0D0*DSU
        ELSE IF(J.EQ.60) THEN
          DS=0.25D0*SINV(J-1)/WTGV(J-1)
          DSU=DS
          DSV=0.0D0
        ELSE
          DS=SINU(J)/WTGU(J)
          DSU=DS
          DSV=SINV(J)/WTGV(J)
        ENDIF
        TMP_INNER=TMP_INNER+DSU*EUJ(1,J)+DSV*EVJ(1,J)+DSU*ETJ(1,J)+DS*EPJ(1,J)
      ENDDO

      call gamil_all_reduce(TMP_INNER)
      INNER = TMP_INNER


!
	RETURN
	END
