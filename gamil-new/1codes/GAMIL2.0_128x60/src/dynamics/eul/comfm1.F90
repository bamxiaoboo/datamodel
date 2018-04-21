module comfm1

    !! (wh 2003.07.09)
    !! (wh 2003.10.23)  change the arrays into allocatable ones
    !! (wh 2003.12.01)  ghi added

    use shr_kind_mod, only: r8 => shr_kind_r8
    use pmgrid, only: beglatexdyn, endlatexdyn, plev, plevp
    use mpi_gamil
    use infnan

    implicit none

    integer  :: itime

    real(r8), allocatable :: su (:,:,:)
    real(r8), allocatable :: sv (:,:,:)
    real(r8), allocatable :: st (:,:,:)
    real(r8), allocatable :: sq (:,:,:)
    real(r8), allocatable :: sut (:,:,:)
    real(r8), allocatable :: svt (:,:,:)
    real(r8), allocatable :: stt (:,:,:)

    real(r8), allocatable :: u  (:,:,:)
    real(r8), allocatable :: v  (:,:,:)
    real(r8), allocatable :: t  (:,:,:)
    real(r8), allocatable :: q  (:,:,:)
    real(r8), allocatable :: ws (:,:,:)
    real(r8), allocatable :: wpa(:,:,:)
    real(r8), allocatable :: ghi(:,:,:)
    real(r8), allocatable :: pes(:,:)       ! ps - pmtop (for fm2003)
    real(r8), allocatable :: ghs(:,:)    
    real(r8), allocatable :: u0  (:,:,:)
    real(r8), allocatable :: v0  (:,:,:)
    real(r8), allocatable :: ws0  (:,:,:)
    real(r8), allocatable :: qt  (:,:,:)
    real(r8), allocatable :: dp  (:,:)
    real(r8), allocatable :: pp  (:,:)
    real(r8), allocatable :: fac  (:,:,:)
    real(r8), allocatable :: fbc  (:,:,:)
    real(r8), allocatable :: uuk  (:,:)
    real(r8), allocatable :: hhk  (:,:)
    real(r8), allocatable :: dus  (:,:)
    real(r8), allocatable :: dps2  (:,:)

    real(r8), allocatable :: ply  (:,:,:)
    real(r8), allocatable :: ply2 (:,:,:)
    real(r8), allocatable :: tb (:,:,:)
    real(r8), allocatable :: tb2 (:,:,:)
    real(r8), allocatable :: qk (:,:,:)
    real(r8), allocatable :: cb (:,:,:)
    real(r8), allocatable :: dcb (:,:,:)
    real(r8), allocatable :: cb0 (:,:,:)
    real(r8), allocatable :: p (:,:)
    real(r8), allocatable :: c0 (:,:)
    integer, allocatable :: nigw(:)
    integer, allocatable :: nigw_2D(:,:)
    
    real(r8), allocatable :: uu (:,:,:)
    real(r8), allocatable :: vv (:,:,:)
    real(r8), allocatable :: tt (:,:,:)
    real(r8), allocatable :: dps0 (:,:)
    real(r8), allocatable :: dps1 (:,:)
    real(r8), allocatable :: hps (:,:)
    real(r8), allocatable :: hh (:,:,:)
    real(r8), allocatable :: ttz (:,:,:)
    real(r8), allocatable :: uz (:,:,:)
    real(r8), allocatable :: vz (:,:,:)
    real(r8), allocatable :: ttv (:,:,:)
    real(r8), allocatable :: dps (:,:)


    real(r8), allocatable :: du (:,:,:)
    real(r8), allocatable :: dv (:,:,:)
    real(r8), allocatable :: dtt (:,:,:)
    real(r8), allocatable :: du0 (:,:,:)
    real(r8), allocatable :: dv0 (:,:,:)
    real(r8), allocatable :: dtt0 (:,:,:)
    real(r8), allocatable :: du1 (:,:,:)
    real(r8), allocatable :: dv1 (:,:,:)
    real(r8), allocatable :: dtt1 (:,:,:)

    real(r8), allocatable :: uk (:,:,:)
    real(r8), allocatable :: vk (:,:,:)
    real(r8), allocatable :: ttk (:,:,:)
    real(r8), allocatable :: psk (:,:)

    real(r8), allocatable :: dq (:,:,:)
    real(r8), allocatable :: uq (:,:,:)
    real(r8), allocatable :: vq (:,:,:)
    real(r8), allocatable :: wq (:,:,:)
    real(r8), allocatable :: pq (:,:)

    real(r8), allocatable :: hu (:,:)
    real(r8), allocatable :: hv (:,:)
    real(r8), allocatable :: cu (:,:)
    real(r8), allocatable :: cv (:,:)
    real(r8), allocatable :: h (:,:)

    real(r8), allocatable :: tuu (:,:,:)
    real(r8), allocatable :: tvv (:,:,:)
    real(r8), allocatable :: A (:,:,:)
    real(r8), allocatable :: QH (:,:,:)
    real(r8), allocatable :: QHSTAR (:,:,:)
    real(r8), allocatable :: USTAR (:,:,:)
    real(r8), allocatable :: VSTAR (:,:,:)
    real(r8), allocatable :: BETA (:,:,:)
    real(r8), allocatable :: FX (:,:,:)
    real(r8), allocatable :: FY (:,:,:)

    real(r8), allocatable :: QC (:,:,:)
    real(r8), allocatable :: QR (:,:,:)
    real(r8), allocatable :: QI (:,:,:)
    real(r8), allocatable :: QS (:,:,:)
    real(r8), allocatable :: QG (:,:,:)
    real(r8), allocatable :: D (:,:,:)
    real(r8), allocatable :: DT (:,:,:)
    real(r8), allocatable :: DS (:,:,:)
    real(r8), allocatable :: DA (:,:,:)
    real(r8), allocatable :: DB (:,:,:)
    real(r8), allocatable :: VR (:,:,:)
    real(r8), allocatable :: HQK (:,:,:)
    real(r8), allocatable :: TK (:,:,:)
    real(r8), allocatable :: HVK (:,:,:)
    real(r8), allocatable :: HUK (:,:,:)
    real(r8), allocatable :: ROT (:,:,:)
    real(r8), allocatable :: RLNT (:,:,:)
    real(r8), allocatable :: RDLN (:,:,:)
    real(r8), allocatable :: RDLT (:,:,:)
    real(r8), allocatable :: TW (:,:,:)
    
    real(r8), allocatable :: psb (:,:)
    real(r8), allocatable :: tsb (:,:)


contains

    subroutine initialize_comfm1
        !
        ! Purpose:  Allocate and initialize the comfm1 arrays.
        !
        allocate (su  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (sv  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (st  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (sq  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (sut  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (svt  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (stt  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))

        allocate (u   (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (v   (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (t   (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (q   (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (ws  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plevp))
        allocate (wpa (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (ghi (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plevp))
        allocate (pes (ilbnd:ihbnd, beglatexdyn:endlatexdyn))
        allocate (ghs (ilbnd:ihbnd, beglatexdyn:endlatexdyn))
        allocate (u0   (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (v0   (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (ws0   (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (qt   (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (dp   (ilbnd:ihbnd, beglatexdyn:endlatexdyn))
        allocate (pp   (ilbnd:ihbnd, beglatexdyn:endlatexdyn))
        allocate (fac   (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plevp))
        allocate (fbc   (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plevp))
        allocate (uuk   (ilbnd:ihbnd, beglatexdyn:endlatexdyn))
        allocate (hhk   (ilbnd:ihbnd, beglatexdyn:endlatexdyn))
        allocate (dus   (ilbnd:ihbnd, beglatexdyn:endlatexdyn))
        allocate (dps2   (ilbnd:ihbnd, beglatexdyn:endlatexdyn))

        allocate (ply   (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plevp))
        allocate (ply2  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (tb   (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (tb2  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (qk   (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (cb  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (dcb  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (cb0  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (p  (ilbnd:ihbnd, beglatexdyn:endlatexdyn))
        allocate (c0  (ilbnd:ihbnd, beglatexdyn:endlatexdyn))
        allocate (nigw  (beglatexdyn:endlatexdyn))
        allocate (nigw_2D  (ilbnd:ihbnd, beglatexdyn:endlatexdyn))

        allocate (uu  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (vv  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (tt  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))

        
        allocate (dps0  (ilbnd:ihbnd, beglatexdyn:endlatexdyn))
        allocate (dps1  (ilbnd:ihbnd, beglatexdyn:endlatexdyn))
        allocate (hps  (ilbnd:ihbnd, beglatexdyn:endlatexdyn))
        allocate (hh  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plevp))
        allocate (ttz  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plevp+1))
        allocate (uz  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plevp+1))
        allocate (vz  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plevp+1))
        allocate (ttv  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (dps  (ilbnd:ihbnd, beglatexdyn:endlatexdyn))


        allocate (du  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (dv  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (dtt  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (du0  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (dv0  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (dtt0  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (du1  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (dv1  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (dtt1  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))

        allocate (uk  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (vk  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (ttk  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (psk  (ilbnd:ihbnd, beglatexdyn:endlatexdyn))

        allocate (dq  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (uq  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (vq  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (wq  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (pq  (ilbnd:ihbnd, beglatexdyn:endlatexdyn))

        allocate (hu  (ilbnd:ihbnd, beglatexdyn:endlatexdyn))
        allocate (hv  (ilbnd:ihbnd, beglatexdyn:endlatexdyn))
        allocate (cu  (ilbnd:ihbnd, beglatexdyn:endlatexdyn))
        allocate (cv  (ilbnd:ihbnd, beglatexdyn:endlatexdyn))
        allocate (h  (ilbnd:ihbnd, beglatexdyn:endlatexdyn))

        allocate (tuu  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (tvv  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (A  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (QH  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (QHSTAR  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (USTAR  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (VSTAR  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (BETA  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (FX  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (FY  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))


        allocate (qc  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (qr  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (qi  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (qs  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (qg  (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (D (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (DT (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (DS (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (DA (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (DB (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (VR (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (HQK (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (TK (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (HVK (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (HUK (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (ROT (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (RLNT (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (RDLN (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (RDLT (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))
        allocate (TW (ilbnd:ihbnd, beglatexdyn:endlatexdyn, plev))

        allocate (psb (ilbnd:ihbnd, beglatexdyn:endlatexdyn))
        allocate (tsb (ilbnd:ihbnd, beglatexdyn:endlatexdyn))

        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,su(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,sv(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,st(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,sq(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,sut(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,svt(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,stt(:,beglatexdyn,1))

        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,u(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,v(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,t(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,q(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plevp,1,1,ws(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,wpa(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plevp,1,1,ghi(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,1,1,1,pes(:,beglatexdyn))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,1,1,1,ghs(:,beglatexdyn))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,u0(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,v0(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,ws0(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,qt(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,1,1,1,dp(:,beglatexdyn))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,1,1,1,pp(:,beglatexdyn))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plevp,1,1,fac(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plevp,1,1,fbc(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,1,1,1,uuk(:,beglatexdyn))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,1,1,1,hhk(:,beglatexdyn))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,1,1,1,dus(:,beglatexdyn))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,1,1,1,dps2(:,beglatexdyn))

        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plevp,1,1,ply(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,ply2(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,tb(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,tb2(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,qk(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,cb(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,dcb(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,cb0(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,1,1,1,p(:,beglatexdyn))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,1,1,1,c0(:,beglatexdyn))

        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,uu(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,vv(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,tt(:,beglatexdyn,1))

        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,1,1,1,dps0(:,beglatexdyn))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,1,1,1,dps1(:,beglatexdyn))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,1,1,1,hps(:,beglatexdyn))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plevp,1,1,hh(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plevp+1,1,1,ttz(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plevp+1,1,1,uz(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plevp+1,1,1,vz(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,ttv(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,1,1,1,dps(:,beglatexdyn))

        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,du(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,dv(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,dtt(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,du0(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,dv0(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,dtt0(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,du1(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,dv1(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,dtt1(:,beglatexdyn,1))

        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,uk(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,vk(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,ttk(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,1,1,1,psk(:,beglatexdyn))

        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,dq(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,uq(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,vq(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,wq(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,1,1,1,pq(:,beglatexdyn))

        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,1,1,1,hu(:,beglatexdyn))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,1,1,1,hv(:,beglatexdyn))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,1,1,1,cu(:,beglatexdyn))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,1,1,1,cv(:,beglatexdyn))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,1,1,1,h(:,beglatexdyn))

        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,tuu(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,tvv(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,A(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,QH(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,QHSTAR(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,USTAR(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,VSTAR(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,BETA(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,FX(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,FY(:,beglatexdyn,1))

        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,qc(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,qr(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,qi(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,qs(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,qg(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,D(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,DT(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,DS(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,DA(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,DB(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,VR(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,HQK(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,TK(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,HVK(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,HUK(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,ROT(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,RLNT(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,RDLN(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,RDLT(:,beglatexdyn,1))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,plev,1,1,TW(:,beglatexdyn,1))

        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,1,1,1,tsb(:,beglatexdyn))
        call register_comm_array(ilbnd,ihbnd,beglatexdyn,endlatexdyn,1,1,1,1,psb(:,beglatexdyn))

        su (:,:,:) = 0.0
        sv (:,:,:) = 0.0
        st (:,:,:) = 0.0
        sq (:,:,:) = 0.0

        u  (:,:,:) = inf
        v  (:,:,:) = inf 
        t  (:,:,:) = inf 
        q  (:,:,:) = inf
        ws (:,:,:) = inf
        wpa(:,:,:) = inf
        ghi(:,:,:) = inf
        pes(:,:)   = inf
        ghs(:,:)   = inf

        psb (:,:) = inf
        tsb (:,:) = inf

        return
    end subroutine initialize_comfm1

end module comfm1
