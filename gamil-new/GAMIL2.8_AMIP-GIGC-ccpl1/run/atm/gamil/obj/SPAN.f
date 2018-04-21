# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/SPAN.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/SPAN.F"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/SPAN.F" 2

# 1 "./params.h" 1
# 15 "./params.h"
 
 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/SPAN.F" 2

!
!----------------------------------------------------------------
!
      SUBROUTINE SPAN(MM1,MM2,MM3,MP1,MP2,MP3,MDJ)
      use pmgrid, only: beglatexdyn, endlatexdyn
      use mpi_gamil
!
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
# 14 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/SPAN.F" 2
!
      INTEGER MM1(ilbnd:ihbnd,beglatexdyn:endlatexdyn)
      INTEGER MM2(ilbnd:ihbnd,beglatexdyn:endlatexdyn)
      INTEGER MM3(ilbnd:ihbnd,beglatexdyn:endlatexdyn)
      INTEGER MP1(ilbnd:ihbnd,beglatexdyn:endlatexdyn)
      INTEGER MP2(ilbnd:ihbnd,beglatexdyn:endlatexdyn)
      INTEGER MP3(ilbnd:ihbnd,beglatexdyn:endlatexdyn)
      INTEGER MDJ(beglatexdyn:endlatexdyn)
!
      INTEGER I,J
!

!************************************************************************
!               Set the values to MM1,MP1;MM2,MP2;MM3,MP3               *
!                    MM1=i-mdj,if MM1<2, MM1=NX-2+MM1                   *
!                  MP1=i+mdj,if MP1>NX-1, MP1=MP1-NX+2                  *
!                 MM2=i-(mdj-1)/2,if MM2<2, MM2=NX-2+MM2                *
!               MP2=i+(mdj+1)/2,if MP2>NX-1, MP2=MP2-NX+2               *
!                 MM3=i-(mdj+1)/2,if MM3<2, MM3=NX-2+MM3                *
!               MP3=i+(mdj-1)/2,if MP3>NX-1, MP3=MP3-NX+2               *
!************************************************************************
      DO J=jbeg0,jend0
         MDJ(J)=1
      END DO
!
      DO J=jbeg0,jend0
      DO I=beglonex,endlonex
	    MM1(I,J)=I-MDJ(J)
	    MP1(I,J)=I+MDJ(J)
	    MM2(I,J)=I-(MDJ(J)-1)/2
	    MP2(I,J)=I+(MDJ(J)+1)/2
	    MM3(I,J)=I-(MDJ(J)+1)/2
	    MP3(I,J)=I+(MDJ(J)-1)/2
      END DO
      END DO
!******************           Made by Wang Bin           *****************
!*************************************************************************

!- check -----------------------------------------------------------------
!      open (10,file='span-s.out')
!      write(10,11) mdj
!      write(10,*) '-------------------mp2-----------------'
!      do j=1,ny
!        write(10,*) 'j=',j,'----------------------'
!        write(10,11) (mp2(i,j),i=1,nx)
!11    format(1x,10i9)
!      enddo
!!      close (10)
!!      stop
!---------------------------------------------------------------(wh)------

	RETURN
	END

