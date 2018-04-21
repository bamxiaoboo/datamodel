# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/AVNEGQ.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/AVNEGQ.F"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/AVNEGQ.F" 2

# 1 "./params.h" 1
# 15 "./params.h"
 
 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/AVNEGQ.F" 2

!!(2003.11.12)
!!--------------------

      SUBROUTINE AVNEGQ(Q,DSGHL)

      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil
      
C     **********************
C     **********************
C
C     AVOIDE   NEGATIVE MIXING RATIO
C
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
# 20 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/AVNEGQ.F" 2


      REAL*8  DSGHL(NL)
      REAL*8  Q(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL),QR(NL),QI
      REAL*8 ZERO
      DATA    ZERO / 0.0E0 /
      INTEGER I,J,K


!$OMP PARALLEL DO PRIVATE (I,J,K,QR,QI)
      DO j = jbeg0, jend0
      if ( j .ge. 2 .and. j .le. 60-1) then
      DO I  = ibeg1, iend1
      DO K  = 1 ,NL
        QR(K)     = Q(I,J,K)
      ENDDO
      DO K  = 2 ,NL
      QI        = QR(K-1)
      IF( QI.LT.ZERO ) THEN
        QR(K-1)  = ZERO
        QR(K )  = QR(K ) + QI*DSGHL(K)
      ENDIF
      ENDDO
      IF( QR(NL).LT.ZERO ) QR(NL) = ZERO
      DO K  = 1 ,NL
      Q(I,J,K)  = QR(K)
      ENDDO
      ENDDO
      ELSE
      DO K  = 1 ,NL
         QR(K)     = Q(ibeg1,J,K)
      ENDDO
      DO K  = 2 ,NL
      QI        = QR(K-1)
      IF( QI.LT.ZERO ) THEN
        QR(K-1)  = ZERO
        QR(K )  = QR(K ) + QI*DSGHL(K)
      ENDIF
      ENDDO
      IF( QR(NL).LT.ZERO ) QR(NL) = ZERO
      DO K  = 1 ,NL
      DO I  = ibeg1,iend1
      Q(I,J,K)  = QR(K)
      ENDDO
      ENDDO
      ENDIF
      ENDDO

      call gamil_arrays_comm(COMM_TO_LEFT,1,Q(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_RIGHT,1,Q(:,beglatexdyn,1))

      RETURN
      END

