# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/SEMIIMP.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/SEMIIMP.F"

# 1 "./misc.h" 1
# 2 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/SEMIIMP.F" 2

# 1 "./params.h" 1
# 15 "./params.h"
 
 
# 3 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/SEMIIMP.F" 2
!
        subroutine semiu(uuk,hhk,p,du,dps1
     _                  ,dus,dps0,dtdy,oux,dsig,c0)
c
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
# 12 "/data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/dynamics/eul/SEMIIMP.F" 2

        integer i,j,i1,i2,k
        real*8 uuk(ilbnd:ihbnd,beglatexdyn:endlatexdyn)
        real*8 hhk(ilbnd:ihbnd,beglatexdyn:endlatexdyn)
        real*8 p(ilbnd:ihbnd,beglatexdyn:endlatexdyn)
        real*8 du(ilbnd:ihbnd,beglatexdyn:endlatexdyn,nl)
        real*8 dps1(ilbnd:ihbnd,beglatexdyn:endlatexdyn)
        real*8 dus(ilbnd:ihbnd,beglatexdyn:endlatexdyn)
        real*8 dps0(ilbnd:ihbnd,beglatexdyn:endlatexdyn)
        real*8 c0(ilbnd:ihbnd,beglatexdyn:endlatexdyn)
        real*8 oux(beglatexdyn:endlatexdyn)
        real*8 dsig(nl)
        real*8 dtdy
        real*8 ai(ilbnd:ihbnd,jbeg0:jend0),
     _         bi(ilbnd:ihbnd,jbeg0:jend0),
     _         di(ilbnd:ihbnd,jbeg0:jend0)
        real*8 dx,us(ilbnd:ihbnd,jbeg0:jend0),duss(ilbnd:ihbnd,jbeg0:jend0)
        real*8 al(NX_LON,jbeg0:jend0,5),ak(NX_LON,jbeg0:jend0),wh(ilbnd:ihbnd,jbeg0:jend0)
        integer comm_request_id1, comm_request_id2
c
c       Semi-implicit scheme
c
c

        call register_comm_array(ilbnd,ihbnd,jbeg0,jend0,1,1,1,1,ai(:,jbeg0))
        call register_comm_array(ilbnd,ihbnd,jbeg0,jend0,1,1,1,1,bi(:,jbeg0))
        call register_comm_array(ilbnd,ihbnd,jbeg0,jend0,1,1,1,1,wh(:,jbeg0))
        call gamil_arrays_comm(COMM_TO_LEFT,1,p(:,beglatexdyn),uuk(:,beglatexdyn),request_id=comm_request_id1)
        call gamil_arrays_comm(COMM_TO_RIGHT,1,p(:,beglatexdyn),c0(:,beglatexdyn),uuk(:,beglatexdyn),request_id=comm_request_id2)

!$OMP PARALLEL DO PRIVATE (I,J,K)
        do j=jbeg1,jend1
          do i=ibeg1,iend1
            duss(i,j)=0.0d0
            do k=1,nl
	          duss(i,j)=duss(i,j)+du(i,j,k)*dsig(k)
            end do
          end do
        end do

        call wait_icomm_request(comm_request_id1)
        call wait_icomm_request(comm_request_id2)

!$OMP PARALLEL DO PRIVATE (I,J,dx)
        do j=jbeg1,jend1
          dx=0.25d0*oux(j)*dtdy
          do i=ibeg1,iend1
	      ai(i,j)=dx*c0(i,j)*(p(i+1,j)+p(i,j))
	      bi(i,j)=dx*c0(i,j)*(p(i,j)+p(i-1,j))
	      di(i,j)=dx*c0(i-1,j)*(p(i,j)+p(i-1,j))
              us(i,j)=uuk(i,j)+dtdy*duss(i,j)
     _                -bi(i,j)*hhk(i,j)+di(i,j)*hhk(i-1,j)
              wh(i,j)=hhk(i,j)+dtdy*dps1(i,j)*c0(i,j)
     _                -ai(i,j)*uuk(i+1,j)+bi(i,j)*uuk(i,j)
          end do
       end do
c
        call gamil_arrays_comm(COMM_TO_RIGHT,1,ai(:,jbeg0),bi(:,jbeg0),wh(:,jbeg0))

!$OMP PARALLEL DO PRIVATE (I,J,K,dx,i1,i2)
        do j=jbeg1,jend1
           do i=ibeg1,iend1
              ak(i,j) = 0.0
              al(i,j,:) = 0.0
              i1=i-1
              i2=i+1
              if (i1.eq.1) i1=NX_LON-1
              if (i2.eq.NX_LON) i2=2
              if (i.lt.NX_LON-1) then
                 al(i,j,1)=-bi(i-1,j)*di(i,j)
                 al(i,j,2)=1.0d0+bi(i,j)*bi(i,j)+ai(i-1,j)*di(i,j)
                 al(i,j,3)=-ai(i,j)*bi(i,j)
                 al(i,j,4)=0.0
                 al(i,j,5)=us(i,j)-bi(i,j)*wh(i,j)+di(i,j)*wh(i-1,j)
              else
                 ak(i,j)= 1.0d0+bi(i,j)*bi(i,j)+ai(i-1,j)*di(i,j)
                 ak(i1,j)=-bi(i-1,j)*di(i,j)
                 ak(i2,j)=-ai(i,j)*bi(i,j)
                 ak(NX_LON-3,j) = ak(2,j)
                 ak(NX_LON,j)=us(i,j)-bi(i,j)*wh(i,j)+di(i,j)*wh(i-1,j)
              end if
           end do
!
           al(2,j,4)=al(2,j,1)
           al(2,j,1)=0.0
           al(NX_LON-2,j,4)=al(NX_LON-2,j,3)
           al(NX_LON-2,j,3)=0.0
        end do

        call t_startf("GAUSS")
        call gauss(al,ak,NX_LON)

        do j=jbeg1,jend1  
           al(NX_LON-1,j,5)=ak(NX_LON,j)
        enddo
        call t_stopf("GAUSS")
        
 !$OMP PARALLEL DO PRIVATE (I,J)
        do j=jbeg1,jend1  
	    do i=ibeg1,iend1
	     dus(i,j)=al(i,j,5)
	    end do
        end do

        call gamil_arrays_comm(COMM_TO_LEFT,1,dus(:,beglatexdyn))

 !$OMP PARALLEL DO PRIVATE (I,J)
        do j=jbeg1,jend1  
	   do i=ibeg1,iend1
	     dps0(i,j)=wh(i,j)-ai(i,j)*dus(i+1,j)+bi(i,j)*dus(i,j)
	     dps0(i,j)=(dps0(i,j)-hhk(i,j))/(dtdy*c0(i,j))-dps1(i,j)
	   end do
c
	   do i=ibeg1,iend1
	     dus(i,j)=(dus(i,j)-uuk(i,j))/dtdy-duss(i,j)
	   end do
        enddo

      do j=jbeg0,jend0
       if(j .eq. 1 .or. j .eq. plat) then
        do i=ibeg0,iend0
           dus(i,j)=0.0
           dps0(i,j)=0.0
        end do
       endif
      end do        

      call gamil_arrays_comm(COMM_TO_LEFT,1,dus(:,beglatexdyn),dps0(:,beglatexdyn))
      call gamil_arrays_comm(COMM_TO_RIGHT,1,dus(:,beglatexdyn),dps0(:,beglatexdyn))

      call remove_comm_array(wh(:,jbeg0))
      call remove_comm_array(bi(:,jbeg0))
      call remove_comm_array(ai(:,jbeg0))

c
c       End of semi-implicit quadratic-conservation scheme
c
        return
        end
