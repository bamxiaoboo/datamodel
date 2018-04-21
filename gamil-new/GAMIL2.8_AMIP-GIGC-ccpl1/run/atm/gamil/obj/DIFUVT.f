# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DIFUVT.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DIFUVT.F"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DIFUVT.F" 2

# 1 "./params.h" 1
# 15 "./params.h"
 
 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DIFUVT.F" 2

!!(2003.11.10)
!!---------------------


      SUBROUTINE DIFUVT(UU,US,VV,VS,WS,P,PES,PLY2,DPS,TT,GHI,HPS,CB,DCB
     _                ,SIGL,DSIG,DY,OUX,OVX,OUY,OVY,SINV,FF
     _                ,CUR,DLT1,DLT2,MP1,MP2,MP3,MM1,MM2,MM3,WTGV
     _                ,DU,DV,DTT,SUT,SVT,STT,FBC,HH,TTZ,UZ,VZ,TTV)

      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil

      IMPLICIT NONE
!
!	This subroutine is to calculate the tendency of the wind and the temperature:
!     DU, DV, DTT


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
# 22 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DIFUVT.F" 2

!
!	The file PARA is to define the parameters related to the model resolution:
!     NX_LON is the grid number in longitude
!     NY is the grid number in latitude
!     NL is the number of vertical layers
!
	  REAL*8
     _       UU  (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  input variable, UU  = UU*sqrt(PES)
     _	  ,US(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL )  !  input variable, US = UU, zonal wind
     _      ,VV  (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !	 input variable, VV  = VV*sqrt(PES)
     _      ,VS(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL )  !  input variable, VS = VV, meridional wind
     _      ,WS(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NZ )  !  input variable, WS = w, vertical velocity
     _      ,P  (ilbnd:ihbnd,beglatexdyn:endlatexdyn   )  !  input variable, P  = sqrt(PES)
     _      ,PES (ilbnd:ihbnd,beglatexdyn:endlatexdyn   )  !	 input variable, Surface pressure
     _      ,PLY2(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  input variable, PLY2=p, Pressure in Sigma Layer
     _      ,DPS(ilbnd:ihbnd,beglatexdyn:endlatexdyn   )  !  input variable,
!                              the tendency of the surface pressure
     _      ,TT (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  input variable, TT=R*T'*Sqrt(PES)/CB,
!                              where T'=T-TB, T is the temperatur,
!                              TBB	is Temperature of the standard atmosphere
     _	  ,GHI (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NZ )  !  input variable,
     _	  ,HH(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NZ )  !  input variable,
!  			     GHI=gz-HBB, gz is the geopotential height,
!                              HBB is the geopotential height of the standard atmopshere
     _      ,CB (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  input variable, CB=Sqrt(R*(KK*TBB-dTB/dlnp)),
!                              where, KK=R/Cp, R is a constant
     _      ,DCB(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  input variable,
     _      ,SIGL(NL     )  !  input constant, the vertical layers
     _      ,DSIG(NL     )  !  input constant, the vertical stepsizes
     _      ,DY             !  input constant,
!                              the horizontal stepsize in meridional direction
     _      ,SINV(beglatexdyn:endlatexdyn)       !  input constant, sin(theta) at half grid j+1/2
     _      ,OUX(beglatexdyn:endlatexdyn)        !  input constant, OUX=1/(RAD*SINU*DX*MDJ)
!                              where, DX is the horizontal stepsize in zonal direction,
!                              MDJ is the leaping length of the central difference
!                              SINU is sin(theta) at intger grid j
     _      ,OVX(beglatexdyn:endlatexdyn)        !  input constant, OUX=1/(RAD*SINV*DX*MDJ)
     _      ,OUY(beglatexdyn:endlatexdyn)        !  input constant, OUY=1/(RAD*SINU*DY*WTGU)
!                              where WTGU is the weighting at the integer grid j
     _      ,OVY(beglatexdyn:endlatexdyn)        !  input constant, OUY=1/(RAD*SINV*DY*WTGU)
     _      ,WTGV(beglatexdyn:endlatexdyn)       !	 input constant,
!                              the weighting at the half grid j+1/2,
     _	  ,DLT1           !  input constant
     _      ,DLT2           !	 input constant
     _      ,FF(beglatexdyn:endlatexdyn)         !	 input constant
     _      ,CUR(beglatexdyn:endlatexdyn)		  !	 input constant
!
	INTEGER
     _       MM1(ilbnd:ihbnd,beglatexdyn:endlatexdyn)     !	 input constant
     _      ,MP1(ilbnd:ihbnd,beglatexdyn:endlatexdyn)     !	 input constant
     _      ,MM2(ilbnd:ihbnd,beglatexdyn:endlatexdyn)     !	 input constant
     _      ,MP2(ilbnd:ihbnd,beglatexdyn:endlatexdyn)     !	 input constant
     _      ,MM3(ilbnd:ihbnd,beglatexdyn:endlatexdyn)     !	 input constant
     _      ,MP3(ilbnd:ihbnd,beglatexdyn:endlatexdyn)     !	 input constant
!
      REAL*8
     _       DU(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)	  !  output variables
     _      ,DV(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)   !  output variables
     _      ,DTT(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  output variables
     _      ,SUT(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)   !  input variables
     _      ,SVT(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)   !  input variables
     _      ,STT(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)   !  input variables
!
      REAL*8  TTZ(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NZ+1)
      REAL*8  UZ(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NZ+1)
      REAL*8  VZ(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NZ+1)
      INTEGER I,J,K
      INTEGER IM1,IP2,IP1,IM2,IM3,IP3
      REAL*8  OY1,R22,R21,R12,DYP0,DXP1,DXP0,DTTP,TO3P
      REAL*8  R11,OUX2,DYP2,O2,O1,TL3,TO3,TO2,TO1
      REAL*8  TL2,TL1,OPK0,OZ2
      REAL*8  OY4,RI0,TO1P,TL3P,TO2P,O2P_TMP,TLP_TMP
      REAL*16 O2P_16(NL),TLP_16(NL),O2P_SUM,TLP_SUM
      REAL*16 O2P_BUF(8,NL),TLP_BUF(8,NL)
      REAL*8  O1P,WKP
      REAL*8  RV1,RV2,OVX4,OVY4,PDXP1,PDXP2,R14
      REAL*8  OUX4,R10,R20,R24,UL2,UL3
      REAL*8  H0,UL1,PX2,FSV,PX1,AXP0
      REAL*8  OPK1,PDYP2,AYP0,PDYP1
      REAL*8  VL2,VL3,OPK3,VL1,PY1
      REAL*8  FSU,PY2,FFYY,OZ4,OPK2
      REAL*8  FS0,FFXX,H1,OYY4
      REAL*8  HPS(ilbnd:ihbnd,beglatexdyn:endlatexdyn)
      REAL*8  FBC(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
      REAL*8  TTV(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
      REAL*8  PLY2_pole_2(NL), PES_pole_2(1), WS_pole_2(NZ), DPS_pole_2(1)
      REAL*8  TTZ_pole_2(NZ+1), P_pole_2(1), CB_pole_2(NL), DCB_pole_2(NL), TTV_pole_2(NL)
      integer :: icomm_request1, icomm_request2, icomm_request3, icomm_request4


      call t_startf("DIFUVT COMM")
      call gamil_arrays_comm(COMM_TO_BOT,1,p(:,beglatexdyn),uu(:,beglatexdyn,1),
     &                       vv(:,beglatexdyn,1),tt(:,beglatexdyn,1),vs(:,beglatexdyn,1),
     &                       request_id=icomm_request1)
      call gamil_arrays_comm(COMM_TO_TOP,1,ghi(:,beglatexdyn,1),ws(:,beglatexdyn,1),
     &                        cb(:,beglatexdyn,1),ply2(:,beglatexdyn,1),uu(:,beglatexdyn,1),
     &                        vv(:,beglatexdyn,1),tt(:,beglatexdyn,1), us(:,beglatexdyn,1),
     &                        vs(:,beglatexdyn,1),request_id=icomm_request2)
      call t_stopf("DIFUVT COMM")

C     AS FOLLOWS, THE IMPROVED LEAP-FROG AND THE REGENERATION OF
C     VARIABLES WILL BE FINISHED IN THE SAME CYCLES.
C
      RI0=DFLOAT(NX_LON-2)
      OY1=1.0D0/(RAD*DY)
      OY4=0.25D0*OY1

C
C     CALCULATING DU/DT, DV/DT AND DTT/DT.
C
C     FORMING THE FACTORS OUT OF ALL THE CYCLES AND CLEANING SOME
C     SPECIAL ARRAIES.
C

      call t_startf("DIFUVT COMP")
!$OMP PARALLEL DO PRIVATE (I,J,K)
      DO K=1,NL
      DO J=jbeg0,jend0
      DO I=beglonex,endlonex
	    TTV(I,J,K)=TT(I,J,K)*FBC(I,J,K)
      END DO
      END DO
      END DO

      call t_startf("DIFUVT COMM")
      call gamil_arrays_comm(COMM_TO_TOP,1,ttv(:,beglatexdyn,1),request_id=icomm_request3)
      call t_stopf("DIFUVT COMM")

!$OMP PARALLEL DO PRIVATE (I,J,K)
      DO K=1,NZ
      DO J=jbeg0,jend0
      DO I=beglonex,endlonex
        HH(I,J,K)=GHI(I,J,K)-HPS(I,J)
      END DO
      END DO
      END DO

      call t_stopf("DIFUVT COMP")


      call t_startf("DIFUVT COMP")
!$OMP PARALLEL DO PRIVATE (I,J,K)
        DO K = 1,NZ+1
          IF(K.EQ.1.OR.K.EQ.NZ+1) THEN
          DO j=jbeg0,jend0
            DO I = beglonex,endlonex
              TTZ(I,J,K) = 0.0D0
              UZ(I,J,K) = 0.0D0
              VZ(I,J,K) = 0.0D0
            ENDDO
          ENDDO
          ELSE
          DO j=jbeg0,jend0
            DO I = beglonex,endlonex
              TTZ(I,J,K) = TT(I,J,K-1)
              UZ(I,J,K) = UU(I,J,K-1)
              VZ(I,J,K) = VV(I,J,K-1)
            ENDDO
          ENDDO
          ENDIF
        ENDDO
      call t_stopf("DIFUVT COMP")

      call t_startf("DIFUVT COMM")
      call wait_icomm_request(icomm_request1)
      call wait_icomm_request(icomm_request2)
      call wait_icomm_request(icomm_request3)
      call t_stopf("DIFUVT COMM")

      DO j=jbeg0,jend0
C
C     CALCULATING DTT/DT AT POLES.
C
       IF(J.EQ.1 .OR. J.EQ.60) THEN
          call t_startf("DIFUVT COMP")
          IF(ibeg1.EQ.2) THEN
             DO K=1,NL
                PLY2_pole_2(K) = PLY2(ibeg1,J,K)
                CB_pole_2(K) = CB(ibeg1,J,K)
                DCB_pole_2(K) = DCB(ibeg1,J,K)
                TTV_pole_2(K) = TTV(ibeg1,J,K)
                WS_pole_2(K) = WS(ibeg1,J,K)
                TTZ_pole_2(K) = TTZ(ibeg1,J,K)
             ENDDO    
             WS_pole_2(NZ) = WS(ibeg1,J,NZ)
             TTZ_pole_2(NZ) = TTZ(ibeg1,J,NZ)
             TTZ_pole_2(NZ+1) = TTZ(ibeg1,J,NZ+1)
             PES_pole_2(1) = PES(ibeg1,J)
             P_pole_2(1) = P(ibeg1,J)
             DPS_pole_2(1) = DPS(ibeg1,J)
          ENDIF
          call t_stopf("DIFUVT COMP")
          call broadcast_lon_data(2,J,PLY2_pole_2,NL,CB_pole_2,NL,DCB_pole_2,NL,TTV_pole_2,NL,
     &                       WS_pole_2,NZ,TTZ_pole_2,NZ+1,PES_pole_2,1,P_pole_2,1,DPS_pole_2,1)
       ENDIF

       IF(J.EQ.1) THEN   ! for the north pole

          call t_startf("DIFUVT COMP")
C
C       CALCULATING ADVECTION OF TT AND VERTICAL MOTION OVER PRESSURE
C       AT POLES.
!$OMP PARALLEL DO PRIVATE (K,I,TLP_SUM,O2P_SUM,TLP_TMP,O2P_TMP)
          DO K=1,NL
            TLP_SUM=0.0D0
            O2P_SUM=0.0D0
            DO I=ibeg1,iend1
              TLP_TMP=VS(I,J,K)*TT(I,J+1,K)
              TLP_SUM=TLP_SUM+TLP_TMP
              O2P_TMP=(P(I,J+1)-P(I,J))*VV(I,J,K)
              O2P_SUM=O2P_SUM+O2P_TMP
            ENDDO
            TLP_BUF(1,K)=TLP_SUM
            O2P_BUF(1,K)=O2P_SUM
          ENDDO
          
          DO K=1,NL
            TLP_16(K)=TLP_BUF(1,K)
            O2P_16(K)=O2P_BUF(1,K)
          ENDDO
          call t_stopf("DIFUVT COMP")

          call t_startf("DIFUVT COMM")
          call gamil_sum_pole_data_phys(J,TLP_16,NL,O2P_16,NL)
          call t_stopf("DIFUVT COMM")

          call t_startf("DIFUVT COMP")
!$OMP PARALLEL DO PRIVATE (K,I,OZ2,OZ4,WKP,O1P,TL3P,TO1P,TO2P,TO3P,DTTP,IP1,IM1,IP2,IM2,
!$   &                     DYP0,PDYP1,PDYP2,AYP0,OPK0,OPK3,VL1,VL2,VL3,FSU,PY1,PY2,FFYY,
!$   &                     TLP_TMP,O2P_TMP)
          DO K=1,NL
            OZ2=0.5D0/DSIG(K)
            OZ4=0.5D0*OZ2
            WKP=2.0D0/RI0*OY1*WTGV(J)
            TLP_TMP=TLP_16(K)
            TLP_TMP=WKP*TLP_TMP
            O2P_TMP=O2P_16(K)
            O2P_TMP=2.0D0*WKP*O2P_TMP*SIGL(K)/PLY2_pole_2(K)
            O1P=(0.5D0*PES_pole_2(1)*(WS_pole_2(K+1)+WS_pole_2(K))
     &        +DPS_pole_2(1)*SIGL(K))/PLY2_pole_2(K)
            TL3P=OZ2*(WS_pole_2(K+1)*TTZ_pole_2(K+2)
     &               -WS_pole_2(K)*TTZ_pole_2(K))
C
            TO1P=CB_pole_2(K)*P_pole_2(1)*O1P
            TO2P=CB_pole_2(K)*P_pole_2(1)*O2P_TMP
            TO3P=(DLT1*CAPA-DLT2*DCB_pole_2(K))*TTV_pole_2(K)*(O1P+O2P_TMP)
C
C     CALCULATING DTT/DT AND FORMING BOTH POLAR BOUNDARIES.
            DTTP=-TLP_TMP-TL3P+TO1P+TO2P+TO3P
            DO I=beglonex,endlonex
              DTT(I,J,K)=DTTP+STT(I,J,K)
              DU(I,J,K)=SUT(I,J,K)
            ENDDO
            DO I=ibeg1,iend1
              IP1=MP1(I,J)
              IM1=MM1(I,J)
              IP2=MP2(I,J)
              IM2=MM2(I,J)
C
              DYP0=P(I,J+1)-P(I,J)
              PDYP1=P(I,J+1)*DYP0
              PDYP2=P(I,J)*DYP0
              AYP0=0.25D0*(P(I,J)+P(I,J+1))
C
              OPK0=1.0D0/PLY2(I,J,K)
              OPK3=1.0D0/PLY2(I,J+1,K)
C
              VL1=0.25D0*OVX(J)*(US(IP2,J+1,K)*VV(IP1,J,K)
     &                          -US(IM2,J+1,K)*VV(IM1,J,K))
              VL2=0.25D0*OVY(J)*(SINV(J+1)*VS(I,J+1,K)
     &                          +SINV(J)*VS(I,J,K))*VV(I,J+1,K)
              VL3=OZ4*((WS(I,J+1,K+1)+WS(I,J,K+1))*VZ(I,J,K+2)
     &                -(WS(I,J+1,K)+WS(I,J,K))*VZ(I,J,K))
C
              FSU=0.25D0*((FF(J+1)+CUR(J+1)*US(I,J+1,K))*UU(I,J+1,K)
     &                   +(FF(J+1)+CUR(J+1)*US(I+1,J+1,K))*UU(I+1,J+1,K))
              PY1=OY1*AYP0*(GHI(I,J+1,K+1)+GHI(I,J+1,K)-GHI(I,J,K+1)-GHI(I,J,K))
              PY2=OY1*(PDYP1*TTV(I,J+1,K)*OPK3*CB(I,J+1,K)
     &                +PDYP2*TTV(I,J,  K)*OPK0*CB(I,J,  K))*SIGL(K)
C
              FFYY=PY1+PY2
              DV(I,J,K)=-VL1-VL2-VL3+FSU-FFYY*WTGV(J)+SVT(I,J,K)

C
            ENDDO
          ENDDO
          call t_stopf("DIFUVT COMP")
        ELSE IF(J.EQ.60) THEN  ! for the south pole

          call t_startf("DIFUVT COMP")
C
C       CALCULATING ADVECTION OF TT AND VERTICAL MOTION OVER PRESSURE
C       AT POLES.
!$OMP PARALLEL DO PRIVATE (K,I,TLP_SUM,O2P_SUM,TLP_TMP,O2P_TMP)
          DO K=1,NL
            TLP_SUM=0.0D0
            O2P_SUM=0.0D0
            DO I=ibeg1,iend1
              TLP_TMP=VS(I,J-1,K)*TT(I,J-1,K)
              TLP_SUM=TLP_SUM+TLP_TMP
              O2P_TMP=(P(I,J)-P(I,J-1))*VV(I,J-1,K)
              O2P_SUM=O2P_SUM+O2P_TMP
            ENDDO
            TLP_BUF(1,K)=TLP_SUM
            O2P_BUF(1,K)=O2P_SUM
          ENDDO
          DO K=1,NL
            TLP_16(K)=TLP_BUF(1,K)
            O2P_16(K)=O2P_BUF(1,K)
          ENDDO

          call t_stopf("DIFUVT COMP")

          call t_startf("DIFUVT COMM")
          call gamil_sum_pole_data_phys(J,TLP_16,NL,O2P_16,NL)
          call t_stopf("DIFUVT COMM")

          call t_startf("DIFUVT COMP")
!$OMP PARALLEL DO PRIVATE (K,I,OZ2,OZ4,WKP,O1P,TL3P,TO1P,TO2P,TO3P,DTTP,TLP_TMP,O2P_TMP)
          DO K=1,NL
            OZ2=0.5D0/DSIG(K)
            OZ4=0.5D0*OZ2
            WKP=2.0D0/RI0*OY1*WTGV(J-1)
            TLP_TMP=TLP_16(K)
            TLP_TMP=-WKP*TLP_TMP
            O2P_TMP=O2P_16(K)
            O2P_TMP=2.0D0*WKP*O2P_TMP*SIGL(K)/PLY2_pole_2(K)
            O1P=(0.5D0*PES_pole_2(1)*(WS_pole_2(K+1)+WS_pole_2(K))
     &        +DPS_pole_2(1)*SIGL(K))/PLY2_pole_2(K)
            TL3P=OZ2*(WS_pole_2(K+1)*TTZ_pole_2(K+2)
     &             -WS_pole_2(K)*TTZ_pole_2(K))
C
            TO1P=CB_pole_2(K)*P_pole_2(1)*O1P
            TO2P=CB_pole_2(K)*P_pole_2(1)*O2P_TMP
            TO3P=(DLT1*CAPA-DLT2*DCB_pole_2(K))*TTV_pole_2(K)*(O1P+O2P_TMP)
C
C     CALCULATING DTT/DT AND FORMING BOTH POLAR BOUNDARIES.
            DTTP=-TLP_TMP-TL3P+TO1P+TO2P+TO3P
            DO I=beglonex,endlonex
              DTT(I,J,K)=DTTP+STT(I,J,K)
              DU(I,J,K)=SUT(I,J,K)
              DV(I,J,K)=SVT(I,J,K)
            ENDDO
          ENDDO
          call t_stopf("DIFUVT COMP")
        ENDIF
      ENDDO


      call t_startf("DIFUVT COMP")
!$OMP PARALLEL DO PRIVATE (I,J,K,OZ2,OZ4,WKP,O1P,TL3P,TO1P,TO2P,
!$   &                           TO3P,DTTP,IP1,IM1,IP2,IM2,DYP0,PDYP1,PDYP2,AYP0,OPK0,OPK3,
!$   &                           VL1,VL2,VL3,FSU,PY1,PY2,FFYY,OUX2,OUX4,R11,R12,R21,R22,R14,
!$   &                           R24,R10,R20,OVX4,OVY4,RV1,RV2,OYY4,DXP0,DXP1,DYP2,TL1,TL2,
!$   &                           TL3,O1,O2,TO1,TO2,TO3,IP3,IM3,PDXP1,PDXP2,AXP0,OPK1,OPK2,
!$   &                           H0,H1,UL1,UL2,UL3,FS0,FSV,PX1,PX2,FFXX)
      DO j=jbeg1,jend1
C
C     CALCULATING DTT/DT FROM J=2 TO J=60-1.
C
C     FORMING THE FACTORS INDEPENDENT ON I, K
C
          OUX2=0.5D0*OUX(J)
          R11=SINV(J)*OUY(J)
          R12=0.5D0*R11
          R21=SINV(J-1)*OUY(J)
          R22=0.5D0*R21
          OUX4=0.5D0*OUX2
          R14=0.5D0*R12
          R24=0.5D0*R22
          R10=R14/(OY1*WTGV(J))
          R20=R24/(OY1*WTGV(J-1))
          OVX4=0.25D0*OVX(J)
          OVY4=0.25D0*OVY(J)
          RV1=SINV(J+1)*OVY4
          RV2=SINV(J-1)*OVY4
          OYY4=OY4*WTGV(J)
C
          DO K=1,NL
            OZ2=0.5D0/DSIG(K)
            OZ4=0.5D0*OZ2
C
            DO I=ibeg1,iend1
C
              IP1=MP1(I,J)
              IM1=MM1(I,J)
              IP2=MP2(I,J)
              IM2=MM2(I,J)
C
C     FORMING THE FACTORS INDEPENDENT ON K.
C
              DXP0=P(I,J)-P(IM1,J)
              DXP1=P(IP1,J)-P(I,J)
              DYP0=P(I,J+1)-P(I,J)
              DYP2=P(I,J)-P(I,J-1)
C
C     TAKING THE ARRAY ELEMENTS APPEARING IN FOLLOWING FORMULAS
C     REPEATEDLY AND PLACING THEM INTO WORKING UNITS.
C
              OPK0=1.0D0/PLY2(I,J,K)
C
C     CALCULATING DTT/DT FOR J=2--NY-1, I=2--NX_LON-1 ANDK=1--NL.
C
              TL1=OUX2*(US(IP2,J,K)*TT(IP1,J,K)-US(IM2,J,K)*TT(IM1,J,K))
              TL2=R12*VS(I,J,K)*TT(I,J+1,K)-R22*VS(I,J-1,K)*TT(I,J-1,K)
              TL3=OZ2*(WS(I,J,K+1)*TTZ(I,J,K+2)-WS(I,J,K)*TTZ(I,J,K))
              O1=OPK0*(0.5D0*PES(I,J)*(WS(I,J,K+1)+WS(I,J,K))
     &          +DPS(I,J)*SIGL(K))
              O2=OPK0*SIGL(K)*(OUX(J)*(DXP1*UU(IP2,J,K)+DXP0*UU(IM2,J,K))
     &                  +R11*DYP0*VV(I,J,K)+R21*DYP2*VV(I,J-1,K))
              TO1=CB(I,J,K)*P(I,J)*O1
              TO2=CB(I,J,K)*P(I,J)*O2

              TO3=(DLT1*CAPA-DLT2*DCB(I,J,K))*TTV(I,J,K)*(O1+O2)
              DTT(I,J,K)=-TL1-TL2-TL3+TO1+TO2+TO3+STT(I,J,K)
            ENDDO
C
            DO I=ibeg1,iend1
C
              IP1=MP1(I,J)
              IM1=MM1(I,J)
              IP2=MP2(I,J)
              IM2=MM2(I,J)
              IP3=MP3(I,J)
              IM3=MM3(I,J)
C
C     FORMING THE FACTORS INDEPENDENT ON K.
C
              DXP0=P(IP3,J)-P(IM3,J)
              PDXP1=P(IP3,J)*DXP0
              PDXP2=P(IM3,J)*DXP0
              DXP1=P(IP2,J)-P(IM2,J)
              AXP0=0.25D0*(P(IP3,J)+P(IM3,J))
              DYP0=P(I,J+1)-P(I,J)
              PDYP1=P(I,J+1)*DYP0
              PDYP2=P(I,J)*DYP0
              DYP2=P(I,J)-P(I,J-1)
              AYP0=0.25D0*(P(I,J+1)+P(I,J))
C
C     TAKING THE ARRAY ELEMENTS APPEARING IN FOLLOWING FORMULAS
C     REPEATEDLY AND PLACING THEM INTO WORKING UNITS.
C
              OPK0=1.0D0/PLY2(I,J,K)
              OPK1=1.0D0/PLY2(IP3,J,K)
              OPK2=1.0D0/PLY2(IM3,J,K)
              OPK3=1.0D0/PLY2(I,J+1,K)
              H1=GHI(I,J,K+1)+GHI(I,J,K)
!             H0=GHI(IP3,J,K+1)+GHI(IP3,J,K)
              H0=HH(IP3,J,K+1)+HH(IP3,J,K)
C
C     CALCULATING DU/DT, DV/DT J=2--NY-1, I=2--NX_LON-1 AND K=1--NL.
C
              UL1=OUX4*((US(I,J,K)+US(IP1,J,K))*UU(IP1,J,K)
     &                 -(US(I,J,K)+US(IM1,J,K))*UU(IM1,J,K))
              UL2=R14*(VS(I,J,K)+VS(I-1,J,K))*UU(I,J+1,K)
     &           -R24*(VS(I,J-1,K)+VS(I-1,J-1,K))*UU(I,J-1,K)
              UL3=OZ4*((WS(I,J,K+1)+WS(I-1,J,K+1))*UZ(I,J,K+2)
     &                -(WS(I,J,K)+WS(I-1,J,K))*UZ(I,J,K))
              FS0=FF(J)+CUR(J)*US(I,J,K)
              FSV=FS0*(R10*(VV(I,J,K)+VV(I-1,J,K))
     &                +R20*(VV(I,J-1,K)+VV(I-1,J-1,K)))
!             PX1=OUX(J)*(AXP0*(H0-GHI(IM3,J,K+1)-GHI(IM3,J,K)))
              PX1=OUX(J)*(AXP0*(H0-HH(IM3,J,K+1)-HH(IM3,J,K)))
              PX2=OUX(J)*(PDXP1*TTV(IP3,J,K)*OPK1*CB(IP3,J,K)
     &                   +PDXP2*TTV(IM3,J,K)*OPK2*CB(IM3,J,K))*SIGL(K)
              FFXX=PX1+PX2
C
              VL1=OVX4*((US(IP2,J+1,K)+US(IP2,J,K))*VV(IP1,J,K)
     &                 -(US(IM2,J+1,K)+US(IM2,J,K))*VV(IM1,J,K))
              VL2=(RV1*VS(I,J+1,K)+OYY4*VS(I,J,K))*VV(I,J+1,K)
     &           -(RV2*VS(I,J-1,K)+OYY4*VS(I,J,K))*VV(I,J-1,K)
              VL3=OZ4*((WS(I,J+1,K+1)+WS(I,J,K+1))*VZ(I,J,K+2)
     &                -(WS(I,J+1,K)+WS(I,J,K))*VZ(I,J,K))
              FSU=0.25D0*(FS0*UU(I,J,K)
     &                   +(FF(J)+CUR(J)*US(I+1,J,K))*UU(I+1,J,K)
     &                   +(FF(J+1)+CUR(J+1)*US(I,J+1,K))*UU(I,J+1,K)
     &                   +(FF(J+1)+CUR(J+1)*US(I+1,J+1,K))*UU(I+1,J+1,K))
              PY1=OY1*(AYP0*(GHI(I,J+1,K+1)+GHI(I,J+1,K)-H1))
              PY2=OY1*(PDYP1*TTV(I,J+1,K)*OPK3*CB(I,J+1,K)
     &                +PDYP2*TTV(I,J,K)*OPK0*CB(I,J,K))*SIGL(K)
              FFYY=PY1+PY2
C
              DU(I,J,K)=-UL1-UL2-UL3-FSV-FFXX+SUT(I,J,K)
              DV(I,J,K)=-VL1-VL2-VL3+FSU-FFYY*WTGV(J)+SVT(I,J,K)
            ENDDO
          ENDDO
	END DO

      call t_stopf("DIFUVT COMP")

      call t_startf("DIFUVT COMM")
      call gamil_arrays_comm(COMM_TO_LEFT,1,DU(:,beglatexdyn,1),DV(:,beglatexdyn,1),DTT(:,beglatexdyn,1),request_id=icomm_request1)
      call gamil_arrays_comm(COMM_TO_RIGHT,1,DU(:,beglatexdyn,1),DV(:,beglatexdyn,1),DTT(:,beglatexdyn,1),request_id=icomm_request2)
      call wait_icomm_request(icomm_request1)
      call wait_icomm_request(icomm_request2)
      call t_stopf("DIFUVT COMM")

!
	RETURN
	END
