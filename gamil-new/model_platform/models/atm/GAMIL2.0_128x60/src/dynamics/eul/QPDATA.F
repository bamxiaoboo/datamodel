#include <misc.h>
#include <params.h>

!##################################################################################
!! (2003.11.30)

      SUBROUTINE QPDATA1(QT, U0, V0, W0, DSGHL, U, V, WS,
     _                   DTDLN, DTDLT, SINU, SINV, WTGU, WTGV, DTDSG,
     _                   DQ, UQ, VQ, WQ)

      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil
!
!     PREDICT POSITIVE DEFINITE FIELD Q        DUE TO 3-D ADVECTION
!           1)  BY USING THE SCHEME     PROPOSED BY R.C.Yu
!           2)  BY USING THE SCHEME     PROPOSED BY P.K.Smolarkiewicz
!
      IMPLICIT NONE

#include <PARADYN>

C
      REAL*8 U0(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
      REAL*8 V0(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
      REAL*8 W0(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
      REAL*8 U(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
      REAL*8 V(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
      REAL*8 WS(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NZ)
      REAL*8 QT(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
      REAL*8 DQ(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
      REAL*8 DSGHL(NL)
      REAL*8 SINU(beglatexdyn:endlatexdyn)
      REAL*8 SINV(beglatexdyn:endlatexdyn)
      REAL*8 WTGU(beglatexdyn:endlatexdyn)
      REAL*8 WTGV(beglatexdyn:endlatexdyn)
!
      INTEGER NONOS,IORD,ISOR
!
      REAL*8 DSNP,DSSP,DTDLN(beglatexdyn:endlatexdyn)
      REAL*8 DTDLT(beglatexdyn:endlatexdyn)
      REAL*8 GC(beglatexdyn:endlatexdyn),DTDSG(NL)
      REAL*8 UQ(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
      REAL*8 VQ(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
      REAL*8 WQ(ilbnd:ihbnd,beglatexdyn:endlatexdyn,NL)
      REAL*8 HALF,ONE
      DATA    HALF,ONE / 0.5E0,1.0E0 /
      INTEGER I,J,K
!
!
!     GET THE ADVECTION VELOCITY
!

!$OMP PARALLEL DO PRIVATE (I,J,K)
      DO K   = 1 ,NL
         DO J   = jbeg0 , jend0
            DO I   = beglonex , endlonex
               UQ(I,J,K)  = HALF*(U (I,J,K)+U0(I,J,K))
               VQ(I,J,K)  = HALF*(V (I,J,K)+V0(I,J,K))
               WQ(I,J,K)  = HALF*(WS(I,J,K)+W0(I,J,K))
               DQ(I,J,K)    = QT(I,J,K)
            ENDDO
         ENDDO
      ENDDO
!
!     PERFORM HORIZONTAL ADVECTION IN SPHERICAL GEOMETRY
!
      CALL TSPAS(QT, UQ, VQ, SINU, SINV, WTGU, WTGV, DTDLT, DTDLN)
!
!     PERFORM THE VERTICAL ADVECTION
!       BY  R.C.Yu
      CALL TSPASW(QT,WQ,DTDSG)
!
!     PERFORM VERTICAL REDISTRIBUTION TO AVOID NEGATIVE Q-H2O
!
      CALL AVNEGQ(QT,DSGHL)

      RETURN
      END
