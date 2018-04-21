# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/STDFSC.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/STDFSC.F"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/STDFSC.F" 2

# 1 "./params.h" 1
# 15 "./params.h"
 
 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/STDFSC.F" 2


C     +++++++++++++++++++++++++++
C     :::::::::::::::::::::::::::
      SUBROUTINE STDFSC(DFS0,DTHDFS,SINU,SINV,WTGU,WTGV,DX,DY
     _                 ,FRDT,FRDS,FRDU,FRDV,FRDP,DXVPN,DXVPS)
C     :::::::::::::::::::::::::::
C     +++++++++++++++++++++++++++


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
# 19 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/STDFSC.F" 2

!
       REAL*8  DFS0,DTHDFS,DX,DY
     _            ,SINU(beglatexdyn:endlatexdyn)
     _            ,SINV(beglatexdyn:endlatexdyn)
     _            ,WTGU(beglatexdyn:endlatexdyn)
     _            ,WTGV(beglatexdyn:endlatexdyn)
     _            ,FRDT(beglatexdyn:endlatexdyn,3)
     _            ,FRDS(beglatexdyn:endlatexdyn,3)
     _            ,FRDU(beglatexdyn:endlatexdyn,3)
     _            ,FRDV(beglatexdyn:endlatexdyn,3)
     _            ,FRDP(beglatexdyn:endlatexdyn,3)
     _            ,DXVPN(beglatexdyn:endlatexdyn)
     _            ,DXVPS(beglatexdyn:endlatexdyn)
     _            ,DXR,DYR,ALN,ALT,ONE,ANT,AVT,AUT,AUO,AVO,ANTS,DNSI,FRTP
     _            ,SNU,SNV,VST,SUI,SUJ,SVI,SVJ,FPO,CDFS,ZERO,HALF,FOUR,FIM
     _            ,DXYV(beglatexdyn:endlatexdyn),DAN,DAS,DAP
	DATA ZERO,HALF,ONE,FOUR/0.0D0,0.5D0,1.0D0,4.0D0/
	INTEGER J,K
C
C     SET CONSTANTS     FOR COMPUTATION OF THE HORIZONTAL DIFFUSION
C

      
      DO J  = jbeg0 ,jend1
         DXYV(J)   = HALF*(RAD*DX*(SINU(J+1)+SINU(J)))*(RAD*DY)
!         DXYV(J)   = (RAD*DX)*(RAD*DY)
      END DO

      if (jend0 .eq. 60) THEN
         DXYV(60)=0.D0
      ENDIF

      call register_comm_array(1,1,beglatexdyn,endlatexdyn,1,1,1,1,DXYV)
      call gamil_arrays_comm(COMM_TO_TOP,1,DXYV)
      call gamil_arrays_comm(COMM_TO_BOT,1,DXYV)
      call remove_comm_array(DXYV)
      
      DO J  = jbeg0,jend0
      IF (J .eq. 1 .or. J .eq. 60) THEN
        DXVPN(J)  = ZERO
        DXVPS(J)  = ZERO      
      ELSE
        DAN       = DXYV(J-1)
        DAS       = DXYV(J  )
        DAP       = DAN + DAS
        DXVPN(J)  = DAN / DAP
        DXVPS(J)  = DAS / DAP
      ENDIF
      ENDDO
      

      FIM=DBLE(IM)
      DO J  = jbeg0 ,jend0
      IF (J .eq. 1 .or. J .eq. 60) THEN
      DO K  = 1 ,3
        FRDT(J,K) = ZERO
        FRDS(J,K) = ZERO
        FRDU(J,K) = ZERO
        FRDV(J,K) = ZERO
        FRDP(J,K) = ZERO
      ENDDO
      ENDIF
      ENDDO
!
      DXR     = DX  * RAD
      DYR     = DY  * RAD
      ALN       = ONE   / DXR
      ALT       = ONE   / DYR
      ANT       = DX  / DY
      AVT       = DXR * ANT
      AUT       = AVT   / FOUR
      AVO       = DXR / FOUR
      AUO       = DXR
      ANTS      = ANT   * ANT
      DNSI      = DX*DX / FIM
      FRTP      = FOUR*ALT  / FIM
C
      DO J  = jbeg1,jend1
        SNU       = SINU(J)
        FRDT(J,1) = ALN   / SNU
        VST       = ALT   / SNU * WTGU(J)
        FRDT(J,2) = VST   * SINV(J)
        FRDT(J,3) = VST   * SINV(J-1)
      ENDDO

      IF (jbeg0 .eq. 1) FRDT(1 ,1)=-FRTP
      IF (jend0 .eq. 60) FRDT(60,1)= FRTP
C
      DO J  = jbeg0 ,jend1
        SNV       = SINV(J)
        VST       = ALT   / SNV * WTGV(J)
        FRDS(J,1) = ALN   / SNV
        FRDS(J,2) = VST   * SINU(J+1)
        FRDS(J,3) = VST   * SINU(J)
      ENDDO
C
      DO J  = jbeg1,jend1
        SVI       = SINV(J)
        SVJ       = SINV(J-1)
        FRDU(J,1) = AUO   * SINU(J)
        FRDU(J,2) = AUT   * SVI*SVI*WTGU(J)
        FRDU(J,3) = AUT   * SVJ*SVJ*WTGU(J)
      ENDDO
C
      DO J  = jbeg0 , jend1
        SUI       = SINU(J)
        SUJ       = SINU(J+1)
        FRDV(J,1) = AVO   * SINV(J)
        FRDV(J,2) = AVT   * SUJ*SUJ*WTGV(J)
        FRDV(J,3) = AVT   * SUI*SUI*WTGV(J)
      ENDDO
C
      DO J  = jbeg1,jend1
        SVI       = SINV(J)
        SVJ       = SINV(J-1)
        FPO       = ANTS  / SINU(J)
        FRDP(J,1) = ONE
        FRDP(J,2) = FPO   * SVI*SVI*SVI*WTGU(J)*WTGV(J)
        FRDP(J,3) = FPO   * SVJ*SVJ*SVJ*WTGU(J)*WTGV(J-1)
      ENDDO

      IF (jbeg0 .eq. 1) FRDP(1 ,1)= + DNSI
      IF (jend0 .eq. 60) FRDP(60,1) = - DNSI
C
      CDFS      = DTHDFS * DFS0*DFS0
      DO K  = 1 ,3
      DO J  = jbeg0 ,jend0
         FRDU(J,K) = FRDU(J,K) * CDFS
         FRDV(J,K) = FRDV(J,K) * CDFS
         FRDP(J,K) = FRDP(J,K) * CDFS
      ENDDO
      ENDDO

      DO K  = 1 ,3
      DO J  = jbeg0 ,jend0
      FRDP(J,K) = 3.0E0 * FRDP(J,K)
      ENDDO
      ENDDO

      RETURN
      END

