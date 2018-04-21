# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/HPAR.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/HPAR.F"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/HPAR.F" 2

# 1 "./params.h" 1
# 15 "./params.h"
 
 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/HPAR.F" 2


!-------------------------------------------------------
!
	SUBROUTINE HPAR(DX,DY,YTHU,YTHV,WTGU,WTGV,MDJ
     _               ,SINU,SINV,OUX,OUY,OVX,OVY,FF,CUR)
!
      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil
      IMPLICIT NONE
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
# 15 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/HPAR.F" 2
!
	INTEGER J
!
!	INPUT VARIBALES:
!          1) HORIZONTAL RESOLUTION PARAMETERS (SPATIAL STEPSIZE): DX,DY
	REAL DX,DY
!
!		 2) VARIBLES RELATED TO ISO-AREA COORDINATE:
	REAL YTHU(beglatexdyn:endlatexdyn)
	REAL YTHV(beglatexdyn:endlatexdyn)
	REAL WTGU(beglatexdyn:endlatexdyn)
	REAL WTGV(beglatexdyn:endlatexdyn)
!
!          3) THE PARAMETER RELATED TO THE FLEXIBLE LEAPING-GRID METHOD
	INTEGER MDJ(beglatexdyn:endlatexdyn)
!
!     OUTPUT VARIABLES:
!
	REAL*8 OUX(beglatexdyn:endlatexdyn)
	REAL*8 OUY(beglatexdyn:endlatexdyn)
	REAL*8 OVX(beglatexdyn:endlatexdyn)
	REAL*8 OVY(beglatexdyn:endlatexdyn)
	REAL*8 SINU(beglatexdyn:endlatexdyn)
	REAL*8 SINV(beglatexdyn:endlatexdyn)
	REAL*8 FF(beglatexdyn:endlatexdyn)
	REAL*8 CUR(beglatexdyn:endlatexdyn)
!
!	WORKING VARIBALES:
!
	REAL YU,YV
!
      DO J=jbeg0,jend0
        IF(J.GE.2.AND.J.LE.60-1) THEN
 	    YU=YTHU(J)
 	    YV=YTHV(J)
          FF(J)=2.0D0*OMGA*COS(YU)
          CUR(J)=COS(YU)/SIN(YU)/RAD
*         FF(J)=-2.0D0*OMGA*COS(YU)
*         CUR(J)=-COS(YU)/SIN(YU)/RAD
*         FF(NY-J+1)=2.0D0*OMGA*COS(YU)
*         CUR(NY-J+1)=COS(YU)/SIN(YU)/RAD
          SINU(J)=SIN(YU)
          SINV(J)=SIN(YV)
          OUX(J)=1.0D0/(RAD*SINU(J)*DX*MDJ(J))
          OUY(J)=1.0D0/(RAD*SINU(J)*DY)*WTGU(J)
          OVX(J)=1.0D0/(RAD*SINV(J)*DX*MDJ(J))
          OVY(J)=1.0D0/(RAD*SINV(J)*DY)*WTGV(J)
        ELSE IF(J.EQ.1) THEN
 		YV=YTHV(J)
 		SINV(1)=SIN(YV)
          OVX(1)=1.0D0/(RAD*SINV(1)*DX*MDJ(J))
          OVY(1)=1.0D0/(RAD*SINV(1)*DY)*WTGV(1)
          FF(J)=2.0D0*OMGA
          CUR(J)=0.0D0
*         FF(J)=-2.0D0*OMGA
*         CUR(J)=0.0D0
*         FF(NY)=2.0D0*OMGA
*         CUR(NY)=0.0D0
          SINU(J)=0.0D0
          OUX(J)=0.0D0
          OUY(J)=0.0D0
        ELSE
          FF(J)=-2.0D0*OMGA
          CUR(J)=0.0D0
*         FF(J)=2.0D0*OMGA
*         CUR(J)=0.0D0
*         FF(1)=-2.0D0*OMGA
*         CUR(1)=0.0D0
          SINU(J)=0.0D0
          SINV(J)=0.0D0
          OUX(J)=0.0D0
          OUY(J)=0.0D0
          OVX(J)=0.0D0
          OVY(J)=0.0D0
        ENDIF
      ENDDO

      call register_comm_array(1,1,beglatexdyn,endlatexdyn,1,1,1,1,SINV)
      call register_comm_array(1,1,beglatexdyn,endlatexdyn,1,1,1,1,SINU)
      call register_comm_array(1,1,beglatexdyn,endlatexdyn,1,1,1,1,FF)
      call register_comm_array(1,1,beglatexdyn,endlatexdyn,1,1,1,1,CUR)
      call register_comm_array(1,1,beglatexdyn,endlatexdyn,1,1,1,1,WTGV)
      call gamil_arrays_comm(COMM_TO_TOP,1,SINV,SINU,FF,CUR,WTGV)
      call gamil_arrays_comm(COMM_TO_BOT,1,SINV,SINU,FF,CUR,WTGV)
      call remove_comm_array(WTGV)
      call remove_comm_array(CUR)
      call remove_comm_array(FF)
      call remove_comm_array(SINU)
      call remove_comm_array(SINV)

	RETURN
	END

