# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DYNAMICS.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DYNAMICS.F"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DYNAMICS.F" 2

# 1 "./params.h" 1
# 15 "./params.h"
 
 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DYNAMICS.F" 2

!! (2003.07)
!! (2003.10.23-24)
!! (2003.11.03-04)
!! (2003.10.09)
!! (wb 2004.02.15)
!------------------------

!!      SUBROUTINE DYNAMICS(US,VS,WS,ws0,pes,t,ghi,DP,sut,svt,stt,FAC,FBC,
!!     _                   PMTOP,ghs,SIGL,DSIG,DTDY,ITIME,IHDIFUS,contn)


!!      subroutine dynamics (us,vs,ws,ws0,pes,t,ghi,dp,sut,svt,stt,fac,fbc,  &
!!                             ghs,          dtdy,itime,ihdifus,contn)

        subroutine dynamics( dtdy,itime,ihdifus,
     _                       us,vs,ws,ws0,pes,t,ghi,ghs,dp,sut,svt,stt,fac,fbc,
     _                       pmtop,sigl,dsig,
     _                       tbb,hbb,cbb,dcbb,psb,tsb,
     _                       dy,wtgu,wtgv,
     _                       dx,sinu,sinv,oux,ouy,ovx,ovy,ff,cur,
     _                       mm1,mp1,mm2,mp2,mm3,mp3,mdj,
     _                       uuk,hhk,dus,dps2,ply2,tb2,cb,dcb,cb0,
     _                       p,c0,nigw,uu,vv,tt,dps0,dps1,hps,hh,ttz,
     _                       uz,vz,ttv,dps,du,dv,dtt,du1,dv1,dtt1,
     _                       uk,vk,ttk,psk, counter)
     

!---------------------------------------------------------------------------------
!
!	    This subroutine is the main part of the model dynamical
!	framework developed by Dr. Bin Wang from LASG/IAP in September
!	of 2001 recoded in May of 2002.
!         In this framework, some new numerical methods are used,
!     including:
!         1) the explicit difference scheme with exact linear and square
!            conservations, developed by Bin Wang, Zhongzhen Ji and
!            Qingcun Zeng, used to solve the atmopsheric equations. By
!            using this scheme, the model can conserve the total mass
!            exactly all the time and keep the conservation of the total
!            available energy when ignoring the outer forcing and friction,
!            and DLT1=DLT2=0.0
!         2) the weighted even-area coordinate along the latitude,
!            developed by Bin Wang, used to reduce the instability
!            of the model at the poles. The coordinate is produced by
!            the subroutine LATMESH
!         3) the flexible leaping-grid method, developed by Bin Wang to
!            further reduce the instability at the ploes. By using the
!            weighted even-area coordinate and the flexible leap-grid method,
!            no filter or smoother is needed in the model.
!         4) the reduction of the standard atmosphere, developed by
!            Qingcun Zeng, to improve the prediction accuracy.
!         This dynamica framework is normalized by Dr. Bin Wang and Rucong Yu,
!     Parallelized by Xin Zhang, translated to F90 by Pu Ye.
!
!---------------------------------------------------------------------------------
      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil
        implicit none


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
# 64 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DYNAMICS.F" 2


# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/comfm2.h" 1

!! (wanhui 2003.07.07)
!! (wanhui 2003.11.04)
!! (b.wang 2004.02.15)




!       real*8  up  (nx,ny,nl)  !
!       real*8  vp  (nx,ny,nl)  !
!       real*8  ttp (nx,ny,nl)  ! variables at step n-1
!       real*8  pps (nx,ny)     !


       real*8  dlt1
       real*8  dlt2

!       real*8  cs0 (nx,ny)
!       real*8  cbs (NX,NY,NL)
!       integer nzad(ny)

       common/comfm2/ dlt1,dlt2
# 66 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/DYNAMICS.F" 2
!
!     1) INPUT CONSTANTS:	PMTOP,ghs,SIGL,DSIG,DTDY,ITIME,IHDIFUS
!
        real*8   PMTOP          !  PMTOP=10hPa, the pressure of the model top
        real*8   SIGL(NL )      !  half-move vertical layers (sigma(k+1/2)) (unit: 1)
        real*8  DSIG(NL )      !  vertical stepsizes (unit: 1)

        real*8   ghs(ilbnd:ihbnd,beglatexdyn:endlatexdyn)      !  geopotential of the ground elevation  (unit: m^2/s^2)
        real*8   DTDY           !  time stepsize      (unit: s)

        INTEGER*4  ITIME    !  integrated time
        INTEGER   IHDIFUS   !  control parameter related to the horizontal


! 2) INPUT ARRAYS:sut,svt,stt

        real*8  sut(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  the forcing on u
        real*8  svt(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  the forcing on v
        real*8  stt(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  the forcing on T


! 3) BOTH INPUT AND OUTPUT ARRAYS: US, VS, pes, t

        real*8  US (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL )  !  zonal wind        (unit: m/s)
        real*8  VS (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL )  !  meridional wind   (unit: m/s)
        real*8  t (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)   !  Temperature       (unit: K  )
        real*8  pes (ilbnd:ihbnd,beglatexdyn:endlatexdyn)   !  Surface pressure  (unit: hPa)

! 4) OUTPUT ARRAYS ONLY: WS, ws0, ghi, DPS

        real*8  WS (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NZ)  !  diagnosed vertical velocity at the present step (unit: 1/s)
        real*8  ws0 (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)  !  diagnosed vertical velocity at the previous step (unit: 1/s)
        real*8  ghi  (ilbnd:ihbnd,beglatexdyn:endlatexdyn,NZ)  !  diagnosed geopotential height (unit: m^2/s^2)
        real*8  DP (ilbnd:ihbnd,beglatexdyn:endlatexdyn)  !  the tendency of the surface pressure (unit: hPa/s)
        real*8  DPS (ilbnd:ihbnd,beglatexdyn:endlatexdyn)


! 5) WORKING VARIABLES AND ARRAYS:

!   (5-1) FOR IAP TRANSFORMATION:                      ! in module comfm2

!   (5-2) FOR STANDARD ATMOSPHERE AT THE P-LAYERS
        real*8  TBB (NA)       !  temperature         (unit: K)
        real*8  HBB (NA)       !  geopotential height (unit: m^2/s^2)
        real*8  CBB (NA)       !
        real*8  DCBB(NA)       !

!   (5-3) FOR STANDARD ATMOSPHERE AT THE SIGMA-LAYERS  ! in module comfm2

!   (5-4) SURFACE VARIBALES
        real*8  PSB(ilbnd:ihbnd,beglatexdyn:endlatexdyn)   !  Surface pressure of the standard atmopshere  (unit: hPa)
        real*8  TSB(ilbnd:ihbnd,beglatexdyn:endlatexdyn)   !  Surface temperature of the standard atmopshere (unit: K  )

!   (5-5) FOR FELXIBLE LEAPING-POINT ZONAL DIFFERENCE  ! in module commap
        integer  MM1(ilbnd:ihbnd,beglatexdyn:endlatexdyn)     !    MM1=i-MDJ      , if MM1<2, MM1=NX_LON-2+MM1
        integer  MP1(ilbnd:ihbnd,beglatexdyn:endlatexdyn)     !    MP1=i+MDJ      , if MP1>NX_LON-1, MP1=MP1-NX_LON+2
        integer  MM2(ilbnd:ihbnd,beglatexdyn:endlatexdyn)     !    MM2=i-(MDJ-1)/2, if MM2<2, MM2=NX_LON-2+MM2
        integer  MP2(ilbnd:ihbnd,beglatexdyn:endlatexdyn)     !    MP2=i+(MDJ+1)/2, if MP2>NX_LON-1, MP2=MP2-NX_LON+2
        integer  MM3(ilbnd:ihbnd,beglatexdyn:endlatexdyn)     !    MM3=i-(MDJ+1)/2, if MM3<2, MM3=NX_LON-2+MM3
        integer  MP3(ilbnd:ihbnd,beglatexdyn:endlatexdyn)     !    MP3=i+(MDJ-1)/2, if MP3>NX_LON-1, MP3=MP3-NX_LON+2
        integer  MDJ(beglatexdyn:endlatexdyn)        !    leaping span of the difference

!   (5-6) FOR EVEN-AREA MERIDIONAL PARTITION
        real*8  WTGU(beglatexdyn:endlatexdyn)       !    area weighting at the normal grid    (unit: 1)
        real*8  WTGV(beglatexdyn:endlatexdyn)       !    area weighting at the half-move grid (unit: 1)
        real*8  DY             !    meridional stepsize in computing mesh

!   (5-7) TENDENCIES                                   ! in module comfm2

!   (5-8) OTHERS
        real*8  DX             !  input constant (unit: s)
        real*8  SINU(beglatexdyn:endlatexdyn)       !  input constant, sin(theta) at INTEGER grid j
        real*8  SINV(beglatexdyn:endlatexdyn)       !  input constant, sin(theta) at half grid j+1/2
        real*8  OUX (beglatexdyn:endlatexdyn)       !  input constant, OUX=1/(RAD*SINU*DX*MDJ)
        real*8  OVX (beglatexdyn:endlatexdyn)       !  input constant, OUX=1/(RAD*SINV*DX*MDJ)
        real*8  OUY (beglatexdyn:endlatexdyn)       !  input constant, OUY=1/(RAD*SINU*DY*WTGU)
        real*8  OVY (beglatexdyn:endlatexdyn)       !  input constant, OUY=1/(RAD*SINV*DY*WTGU)
        real*8  FF  (beglatexdyn:endlatexdyn)       !  input constant
        real*8  CUR (beglatexdyn:endlatexdyn)       !  input constant

        real*8  WK0,WK5,WK6,BYY1,BYY2,BYY3,DT2,TE,TM,TMS(beglatexdyn:endlatexdyn)
        real*8  TMS_BUF(16,beglatexdyn:endlatexdyn)
        real*8  UUK(ilbnd:ihbnd,beglatexdyn:endlatexdyn)
        real*8  HHK(ilbnd:ihbnd,beglatexdyn:endlatexdyn)
        real*8  DUS(ilbnd:ihbnd,beglatexdyn:endlatexdyn)
        real*8  DPS2(ilbnd:ihbnd,beglatexdyn:endlatexdyn)
        real*8  INNER
        real*8  FAC(ilbnd:ihbnd,beglatexdyn:endlatexdyn,nz)
        real*8  FBC(ilbnd:ihbnd,beglatexdyn:endlatexdyn,nl)
        integer  I,J,K,KPP,KWB

        real*8  ply2(ilbnd:ihbnd,beglatexdyn:endlatexdyn,nl)
        real*8  tb2(ilbnd:ihbnd,beglatexdyn:endlatexdyn,nl)
        real*8  cb(ilbnd:ihbnd,beglatexdyn:endlatexdyn,nl)
        real*8  dcb(ilbnd:ihbnd,beglatexdyn:endlatexdyn,nl)
        real*8  cb0(ilbnd:ihbnd,beglatexdyn:endlatexdyn,nl)
        real*8  p(ilbnd:ihbnd,beglatexdyn:endlatexdyn)
        real*8  c0(ilbnd:ihbnd,beglatexdyn:endlatexdyn)
        integer nigw(beglatexdyn:endlatexdyn)

        real*8  uu(ilbnd:ihbnd,beglatexdyn:endlatexdyn,nl)
        real*8  vv(ilbnd:ihbnd,beglatexdyn:endlatexdyn,nl)
        real*8  tt(ilbnd:ihbnd,beglatexdyn:endlatexdyn,nl)
        real*8  dps0(ilbnd:ihbnd,beglatexdyn:endlatexdyn)
        real*8  dps1(ilbnd:ihbnd,beglatexdyn:endlatexdyn)
        real*8  hps(ilbnd:ihbnd,beglatexdyn:endlatexdyn)
        real*8  hh(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NZ)
        real*8  ttz(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NZ+1)
        real*8  uz(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NZ+1)
        real*8  vz(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NZ+1)
        real*8  ttv(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)

        real*8 du(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
        real*8 dv(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
        real*8 dtt(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
        real*8 du1(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
        real*8 dv1(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
        real*8 dtt1(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
        
        real*8 uk(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
        real*8 vk(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
        real*8 ttk(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
        real*8 psk(ilbnd:ihbnd,beglatexdyn:endlatexdyn)
        integer counter

! added by liuli begin
        real*8  tmp_3d(NX_LON,plat,NL*20)
! added by liuli end

        EXTERNAL INNER

!--------------------------------------------------------------------------
! AS FOLLOWS, THE IMPROVED LEAP-FROG AND THE REGENERATION OF
! VARIABLES WILL BE FINISHED IN THE SAME CYCLES.
!--------------------------------------------------------------------------


        DT2=0.5D0*DTDY
!
	IF (ITIME.EQ.0) THEN
!
!	FOR THE NORMAL RUN, DLT1 AND DTL2 MUST BE 1.0D0
!
 	DLT1=1.0d0
 	DLT2=1.0d0
!
!       DLT1=0.OD0 AND DLT2=0.0D0 ONLY FOR TESTING
!       THE CONSERVATION OF AVAILABLE ENERGY
!
!	DLT1=0.0d0
!	DLT2=0.0d0
!

      call gamil_arrays_comm(COMM_TO_LEFT,1,pes(:,beglatexdyn)) 
      call gamil_arrays_comm(COMM_TO_RIGHT,1,pes(:,beglatexdyn)) 

!$OMP PARALLEL DO PRIVATE (I,J,K,WK5,KPP,WK6)
        DO J=jbeg0,jend0
	      TMS_BUF(1,J)=0.0
          DO I=beglonex,endlonex
	        P (I,J)=SQRT(pes(I,J))
	        C0(I,J)=SQRT(RD*TSB(I,J)/PSB(I,J))
          ENDDO
          DO K=1,nl
            DO I=beglonex,endlonex
              ply2(I,J,K) = pes(I,J)*SIGL(K) + PMTOP
              WK5=ply2(I,J,K)/DPALIB
              KPP=INT(WK5)
              WK6=WK5-DFLOAT(KPP)
              CB(I,J,K)=(1.D0-WK6)*CBB(KPP)+WK6*CBB(KPP+1)
              DCB(I,J,K)=(1.D0-WK6)*DCBB(KPP)+WK6*DCBB(KPP+1)
              tb2(I,J,K)=(1.D0-WK6)*TBB(KPP)+WK6*TBB(KPP+1)
              CB0(I,J,K)=CB(I,J,K)*P(I,J)/ply2(I,J,K)
              WK5=CB0(I,J,K)*P(I,J)*DSIG(K)
              IF (TMS_BUF(1,J).LT.WK5) TMS_BUF(1,J)=WK5
            ENDDO
          ENDDO
        ENDDO

        DO J=jbeg0,jend0
          TMS(J)=TMS_BUF(1,J)
        ENDDO
        call gamil_min_lat_row_data(TMS(jbeg0:),jend0-jbeg0+1)

        DO J=jbeg0,jend0
           KPP=DTDY*TMS(J)*OUX(J)*2.0+0.0001
           NIGW(J)=KPP+1
         ENDDO
!
        call t_startf("IAPTRSF")
        CALL IAPTRSF(US,VS,WS,pes,t,UU,VV,P,TT,tb2,CB,1)
        call t_stopf("IAPTRSF")

        call t_startf("DIFPS0")
        CALL DIFPS0(UU,P,DPS0,DSIG,OUX,MP1,MP2,MM1,MM2)
        call t_stopf("DIFPS0")
!        
        call t_startf("DIFPS1")
        CALL DIFPS1(UU,VV,P,pes,WS,DPS1,DPS0
     _             ,DSIG,DY,OUX,OUY,SINV,MP1,MP2,MM1,MM2,WTGV)
        call t_stopf("DIFPS1")
!       
        call t_startf("DIAG")
        CALL DIAG(UU,VV,P,pes,ply2,TT,US,VS,t,ghi
     _           ,HPS,PMTOP,PSB,TSB,tb2,CB,FAC,DSIG)
        call t_stopf("DIAG")
!

!$OMP PARALLEL DO PRIVATE(I,J,K)
	DO K=1,nl
	DO J=jbeg0,jend0
	DO I=beglonex,endlonex
	   ws0(I,J,K)=WS(I,J,K)
	END DO
	END DO
	END DO
!
	END IF


!
!$OMP PARALLEL DO PRIVATE(I,J,K)
        DO J=jbeg0,jend0
        DO I=beglonex,endlonex
           PSK(I,J)=pes(I,J)
           DPS(I,J)=DPS0(I,J)+DPS1(I,J)
        ENDDO
        DO K=1,NL
        DO I=beglonex,endlonex
           UK (I,J,K)=UU(I,J,K)
           VK (I,J,K)=VV(I,J,K)
           TTK(I,J,K)=TT(I,J,K)
        ENDDO
        ENDDO
        DO I=beglonex,endlonex
           HHK(I,J)=C0(I,J)*(pes(I,J)+PMTOP-PSB(I,J))
           UUK(I,J)=0.0d0
           DO K=1,NL
             UUK(I,J)=UUK(I,J)+UU(I,J,K)*DSIG(K)
           ENDDO
        ENDDO
        ENDDO
!
        DO KWB=1,2
!
!-------------------------------------------------------------------------------
! CALCULATING DU/DT, DV/DT AND DTT/DT.
!  FORMING THE FACTORS OUT OF ALL THE CYCLES AND CLEANING SOME SPECIAL ARRAIES.
!-------------------------------------------------------------------------------
        call t_startf("DIFUVT")
        CALL DIFUVT(UU,US,VV,VS,WS,P,pes,ply2,DPS,TT,ghi,HPS,CB,DCB
     _           ,SIGL,DSIG,DY,OUX,OVX,OUY,OVY,SINV,FF
     _           ,CUR,DLT1,DLT2,MP1,MP2,MP3,MM1,MM2,MM3,WTGV
     _           ,DU,DV,DTT,sut,svt,stt,FBC,hh,ttz,uz,vz,ttv)
!!
        call t_stopf("DIFUVT")
        call t_startf("SEMIU")
        CALL SEMIU(UUK,HHK,P,DU,DPS1,DUS,DPS0,DTDY,OUX,DSIG,C0)
        call t_stopf("SEMIU")
!

!$OMP PARALLEL 
!$OMP DO PRIVATE (I,J,K,WK5,KPP,WK6)
        DO J=jbeg0,jend0
        DO I=beglonex,endlonex
           DPS2(I,J)=DPS1(I,J)+DPS0(I,J)
           pes (I,J)=PSK (I,J)+DT2*DPS2(I,J)
        ENDDO
!
        DO K=1,nl
        DO I=beglonex,endlonex
           DU(I,J,K)=DU (I,J,K)+    DUS(I,J  )
           UU(I,J,K)=UK (I,J,K)+DT2*DU (I,J,K)
           TT(I,J,K)=TTK(I,J,K)+DT2*DTT(I,J,K)
           VV(I,J,K)=VK (I,J,K)+DT2*DV (I,J,K)
           ply2(I,J,K) = pes(I,J)*SIGL(K) + PMTOP
           WK5=ply2(I,J,K)/DPALIB
           KPP=INT(WK5)
           WK6=WK5-DFLOAT(KPP)
           CB(I,J,K)=(1.D0-WK6)*CBB(KPP)+WK6*CBB(KPP+1)
           DCB(I,J,K)=(1.D0-WK6)*DCBB(KPP)+WK6*DCBB(KPP+1)
	   CB0(I,J,K)=CB(I,J,K)*P(I,J)/ply2(I,J,K)
        ENDDO
        ENDDO
        ENDDO
!$OMP END PARALLEL


!
!	BYY2=INNER(DU,DV,DTT,DPS0,0,UU,VV,TT,pes,1
!    _             ,DSIG,TSB,PSB,PMTOP,SINU,SINV,WTGU,WTGV)
! 	print *,BYY2
! 	STOP
!
        IF (KWB.LT.2) THEN
!
!$OMP PARALLEL DO PRIVATE (I,J)
        DO J=jbeg0,jend0
        DO I=beglonex,endlonex
           P  (I,J)=SQRT(pes(I,J))
        ENDDO
        ENDDO
!
        ENDIF
!
        call t_startf("DIFPS1")
        CALL DIFPS1(UU,VV,P,pes,WS,DPS1,DPS0
     _             ,DSIG,DY,OUX,OUY,SINV,MP1,MP2,MM1,MM2,WTGV)
        call t_stopf("DIFPS1")
!
        call t_startf("DIAG")
        CALL DIAG(UU,VV,P,pes,ply2,TT,US,VS,t,ghi
     _           ,HPS,PMTOP,PSB,TSB,tb2,CB,FAC,DSIG)
        call t_stopf("DIAG")

!
!$OMP PARALLEL DO PRIVATE (I,J)
        DO J=jbeg0,jend0
        DO I=beglonex,endlonex
           DPS (I,J)=DPS0(I,J)+DPS1(I,J)
        ENDDO
        ENDDO

        ENDDO
!


!     CALCULATING DU/DT, DV/DT AND DTT/DT.
!
!     FORMING THE FACTORS OUT OF ALL THE CYCLES AND CLEANING SOME
!     SPECIAL ARRAIES.
!
!
      call t_startf("DIFUVT")
      CALL DIFUVT(UU,US,VV,VS,WS,P,pes,ply2,DPS,TT,ghi,HPS,CB,DCB
     _           ,SIGL,DSIG,DY,OUX,OVX,OUY,OVY,SINV,FF
     _           ,CUR,DLT1,DLT2,MP1,MP2,MP3,MM1,MM2,MM3,WTGV
     _           ,DU1,DV1,DTT1,sut,svt,stt,FBC,hh,ttz,uz,vz,ttv)
      call t_stopf("DIFUVT")
!
!$OMP PARALLEL 
!$OMP DO PRIVATE (I,J,K)
        DO K=1,nl
        DO J=jbeg0,jend0
          DO I=beglonex,endlonex
             DU1(I,J,K)=DU1(I,J,K)+DUS(I,J)
          ENDDO
        ENDDO
        END DO
!$OMP END DO NOWAIT
!
!$OMP DO PRIVATE (I,J)
        DO J=jbeg0,jend0
        DO I=beglonex,endlonex
           DPS1(I,J)=DPS (I,J)
           DPS (I,J)=DPS2(I,J)
        ENDDO
        ENDDO
!$OMP END PARALLEL
!
!     To deduct the inner gravity waves from the tendences of the zonal wind
!     and the temperature during the long-time-stepsize integrations
!
       call t_startf("MINUS_INGW")
       CALL MINUS_INGW(UU,P,TT,CB0,DSIG,OUX,DU1,DTT1,NIGW)
       call t_stopf("MINUS_INGW")
!
       call t_startf("INNER")
        BYY1=INNER(DU1,DV1,DTT1,DPS1,0,DU ,DV ,DTT ,DPS ,0
     _            ,DSIG,TSB,PSB,PMTOP,SINU,SINV,WTGU,WTGV)
       call t_stopf("INNER")
!
!$OMP PARALLEL DO PRIVATE (I,J)
        DO J=jbeg0,jend0
        DO I=beglonex,endlonex
           pes(I,J)=PSK(I,J)+DT2*DPS1(I,J)
           P(I,J)=SQRT(pes(I,J))
        ENDDO
        ENDDO
!
!       To add the inner gravity waves to the tendences of the zonal wind
!       and the temperature by the method of short-time-stepsize integrations
!
        call t_startf("PLUS_INGW")
        CALL PLUS_INGW(UK,P,TTK,CB0,DSIG,OUX,DTDY,SINU,WTGU,DU1,DTT1,BYY1,NIGW)
        call t_stopf("PLUS_INGW")
!
        call t_startf("INNER")
	BYY3=INNER(DU1,DV1,DTT1,DPS1,0,DU1,DV1,DTT1,DPS1,0
     _            ,DSIG,TSB,PSB,PMTOP,SINU,SINV,WTGU,WTGV)
        call t_stopf("INNER")

        DT2=DTDY*BYY1/BYY3

!$OMP PARALLEL DO PRIVATE (I,J,K,WK5,KPP,WK6)
        DO J=jbeg0,jend0
        DO I=beglonex,endlonex
           pes(I,J)=PSK(I,J)+DT2*DPS1(I,J)
           P(I,J)=SQRT(pes(I,J))
        ENDDO
!
	    TMS_BUF(1,J)=0.0
        DO K=1,nl
        DO I=beglonex,endlonex
           ply2(I,J,K) = pes(I,J)*SIGL(K) + PMTOP
           WK5=ply2(I,J,K)/DPALIB
           KPP=INT(WK5)
           WK6=WK5-DFLOAT(KPP)
           tb2(I,J,K)=(1.D0-WK6)*TBB(KPP)+WK6*TBB(KPP+1)
           CB(I,J,K)=(1.D0-WK6)*CBB(KPP)+WK6*CBB(KPP+1)
           DCB(I,J,K)=(1.D0-WK6)*DCBB(KPP)+WK6*DCBB(KPP+1)
           CB0(I,J,K)=CB(I,J,K)*P(I,J)/ply2(I,J,K)
           WK5=CB0(I,J,K)*P(I,J)*DSIG(K)
           IF (TMS_BUF(1,J).LT.WK5) TMS_BUF(1,J)=WK5
           TT(I,J,K)=TTK(I,J,K)+DT2*DTT1(I,J,K)
           UU(I,J,K)=UK (I,J,K)+DT2*DU1 (I,J,K)
           VV(I,J,K)=VK (I,J,K)+DT2*DV1 (I,J,K)
        ENDDO
        ENDDO
        ENDDO
!
        DO J=jbeg0,jend0
          TMS(J)=TMS_BUF(1,J)
        ENDDO
        call gamil_min_lat_row_data(TMS(jbeg0:),jend0-jbeg0+1)

        DO J=jbeg0,jend0
        KPP=DTDY*TMS(J)*OUX(J)*2.0+0.0001
        NIGW(J)=KPP+1
        ENDDO

 	I=ITIME
        k=i/86400
!!**  	IF (K*86400.EQ.I) THEN
	   IF (DLT1*DLT2.EQ.0) THEN
          call t_startf("INNER")
 	      BYY2=INNER(UK,VK,TTK,PSK,1,UK,VK,TTK,PSK,1
     _                ,DSIG,TSB,PSB,PMTOP,SINU,SINV,WTGU,WTGV)
	      BYY1=INNER(UU,VV,TT,pes,1,UU,VV,TT,pes,1
     _                ,DSIG,TSB,PSB,PMTOP,SINU,SINV,WTGU,WTGV)
          call t_stopf("INNER")
     
          call t_startf("WRITE")
          if (is_rootproc) then
            WRITE(6,*) float(ITIME)+DTDY,BYY1,DT2
!           WRITE(6,*) BYY2,BYY1,DT2
!          WRITE(112,*) ITIME,BYY1,DT2
          endif
          call t_stopf("WRITE")
	   ELSE
            call t_startf("TEM")
            CALL TEM(UU,VV,pes,P,TT,ghs,CB,tb2,DX,DY,DSIG
     _              ,SINU,SINV,WTGU,WTGV,TE,TM)
            call t_stopf("TEM")
          call t_startf("WRITE")
          if (is_rootproc) then
            WRITE(6,'(1x,3e25.18)') TE,TM,DT2  !---sxj--
           ! WRITE(112,'(1x,3e25.18)') TE,TM,DT2
          endif
          call t_stopf("WRITE")
	   END IF
!!** 	END IF
!
        call t_startf("DIFPS0")
        CALL DIFPS0(UU,P,DPS0,DSIG,OUX,MP1,MP2,MM1,MM2)
        call t_stopf("DIFPS0")

        call t_startf("DIFPS1")
        CALL DIFPS1(UU,VV,P,pes,WS,DPS1,DPS0
     _                ,DSIG,DY,OUX,OUY,SINV,MP1,MP2,MM1,MM2,WTGV)
        call t_stopf("DIFPS1")
!
        call t_startf("DIAG")
        CALL DIAG(UU,VV,P,pes,ply2,TT,US,VS,t,ghi
     _           ,HPS,PMTOP,PSB,TSB,tb2,CB,FAC,DSIG)
        call t_stopf("DIAG")
!
!$OMP PARALLEL DO PRIVATE (I,J)
        DO J=jbeg0,jend0
        DO I=beglonex,endlonex
 	   DP(I,J)=DPS0(I,J)+DPS1(I,J)
 	END DO
 	END DO
        


      RETURN
      END
