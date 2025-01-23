      MODULE technical

       include 'typedef.inc'

        PUBLIC :: kin_p, kin_n, levp, levn, lhfp, lhfn,
     &            Up_HFB, Vp_HFB, Un_HFB, Vn_HFB,
     &            Ap_HFB, Bp_HFB, An_HFB, Bn_HFB,
     &            rhop_HFB, rhon_HFB, kapp_HFB, kapn_HFB,
     &            Vpp,Vnn,Vpn,Fpp,Fnn,Fpn,tran_p,tran_n,
     &            V3B,V3BNO2, 
     &            trE0_p,trE0_n,trE1_p,trE1_n,trE2_p,trE2_n,
     &            trE3_p,trE3_n,trEN_p,trEN_n,trS1_p,trS1_n,
     &            trM1s_p,trM1s_n,trM1l_p,trM1l_n,
     &            trE1_p_dens,trE1_n_dens,
     &            rad_den1,H11p,H11n,zcross,weight,
     &            lnl, lev3, lpoint, idim3j,idim3jt,lev3ord,
     &            lev1pn,lp1,lp2,cg3,larrow1,larrow2,larrowm,iphase,
     &            lev1pnm,levtz,sixj1, itpoi1, klpoi1
     &            energy_HF,energy_PT2
     &            r4Z,r2Z,r4N,r2N,r4tot,r2tot
     &            idim_v3b,ipozjt1,ipozt1

        double precision :: energy_HF,energy_PT2
        double precision :: r4Z,r2Z,r4N,r2N,r4tot,r2tot

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: kin_p(:,:), kin_n(:,:)

        type(level_type),dimension(:),allocatable, save :: levp,levn
        type(level_type),dimension(:),allocatable, save :: lhfp,lhfn

        type(level_type),dimension(:),allocatable, save :: lnl

        type(level_type),dimension(:),allocatable, save :: lev1pn
        type(levelm_type),dimension(:),allocatable, save :: lev1pnm

        type(level3_type),dimension(:),allocatable,save :: lev3

        type(leveltz_type),dimension(:),allocatable,save :: levtz

        INTEGER, ALLOCATABLE, SAVE :: lpoint(:,:,:,:,:,:)
        integer, allocatable, save :: lev3ord(:)
        integer (kind=8), allocatable, save :: idim3j(:)
        integer (kind=8), allocatable, save :: idim3jt(:,:)
        INTEGER, ALLOCATABLE, SAVE :: lp1(:),lp2(:)

        INTEGER, ALLOCATABLE, SAVE :: itpoi1(:,:,:)
        INTEGER, ALLOCATABLE, SAVE :: klpoi1(:,:)

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Up_HFB(:,:), Vp_HFB(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Un_HFB(:,:), Vn_HFB(:,:)

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Ap_HFB(:,:), Bp_HFB(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: An_HFB(:,:), Bn_HFB(:,:)

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: rhop_HFB(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: rhon_HFB(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: kapp_HFB(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: kapn_HFB(:,:)

        real(kind=8), ALLOCATABLE, SAVE :: Vpp(:,:,:,:,:)
        real(kind=8), ALLOCATABLE, SAVE :: Vnn(:,:,:,:,:)
        real(kind=8), ALLOCATABLE, SAVE :: Vpn(:,:,:,:,:)

        real(kind=8), ALLOCATABLE, SAVE :: Vpp_2b(:,:,:,:,:)
        real(kind=8), ALLOCATABLE, SAVE :: Vnn_2b(:,:,:,:,:)
        real(kind=8), ALLOCATABLE, SAVE :: Vpn_2b(:,:,:,:,:)

!!!I'm working on it
        real(kind=4), ALLOCATABLE, SAVE :: V3B_ar(:), V3B_ar2(:)
!        INTEGER, ALLOCATABLE, SAVE :: V3B_bra(:),V3B_ket(:)
        INTEGER(kind=8), ALLOCATABLE, SAVE :: V3B_pair(:)        

        REAL, ALLOCATABLE, SAVE :: V3B(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: V3BNO2(:,:,:,:,:,:,:)

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Fpp(:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Fnn(:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Fpn(:,:,:,:)

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: tran_p(:,:),tran_n(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trE0_p(:,:),trE0_n(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trE1_p(:,:),trE1_n(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trE2_p(:,:),trE2_n(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trE3_p(:,:),trE3_n(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trEN_p(:,:),trEN_n(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trS1_p(:,:),trS1_n(:,:)

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trM1s_p(:,:),trM1s_n(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trM1l_p(:,:),trM1l_n(:,:)

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trE1_p_dens(:,:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trE1_n_dens(:,:,:)

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: rad_den1(:)

        double precision, allocatable, save :: H11p(:,:),H11n(:,:)

        double precision, allocatable, save :: zcross(:),weight(:)

        double precision, allocatable, save :: cg3(:,:,:,:,:)
        double precision, allocatable, save :: sixj1(:,:,:,:,:,:)

        type(elem_type) :: larrow1,larrow2,larrowm
        INTEGER :: iphase
        integer(kind=8) :: idim_v3b
        integer(kind=8) :: ipozjt1,ipozt1,ipoz,ipozjt

      END MODULE technical
