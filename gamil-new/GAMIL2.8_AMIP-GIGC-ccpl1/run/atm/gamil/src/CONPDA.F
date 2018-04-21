#include <misc.h>
#include <params.h>


      SUBROUTINE CONPDA(DTIMEQ, DX, DY, SINU, WTGV, DSIG, FIRST, ISOR,
     $                  IORD, IP, DTDLT, DTDLN, DTDSG)

      use pmgrid, only: beglatexdyn, endlatexdyn, plat
      use mpi_gamil

      IMPLICIT NONE

#include <PARADYN>

      REAL*8 DTIMEQ, DX, DY, FIRST
      REAL*8 SINU(beglatexdyn:endlatexdyn)
      REAL*8 WTGV(beglatexdyn:endlatexdyn)
      REAL*8 DSIG(NL)
      REAL*8 DTDLT(beglatexdyn:endlatexdyn)
      REAL*8 DTDLN(beglatexdyn:endlatexdyn)
      REAL*8 DTDSG(NL)
      INTEGER ISOR, IORD, IP(NX_LON)
      INTEGER I, J, K

      IF (FIRST .LE. 0.0D0) THEN
         IF (ISOR .EQ. 3) IORD = MAX0(IORD, 3)
         DO I = 1, NX_LON
           IP(I) = MOD(I+IM/2-1, IM)+1
         END DO
      END IF

      DO J = jbeg0, jend1
         DTDLT(J) = DTIMEQ/(RAD*DY)
      END DO
      DO J = jbeg1, jend1
         DTDLN(J) = DTIMEQ/(RAD*SINU(J)*DX)
      END DO
      DO K  = 1 ,NL
         DTDSG(K)  = DTIMEQ / DSIG(K)
      END DO

      RETURN
      END
