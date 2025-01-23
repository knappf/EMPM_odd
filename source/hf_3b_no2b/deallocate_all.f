      subroutine deallocate_all
       
       USE technical

       deallocate(levp,levn,lhfp,lhfn,lnl)
       deallocate(lev1pn,lev1pnm,levtz)
       deallocate(kin_p,kin_n)
       deallocate(Vpp,Vnn,Vpn)
       deallocate(Up_HFB,Un_HFB,Vp_HFB,Vn_HFB)
       deallocate(Ap_HFB,An_HFB,Bp_HFB,Bn_HFB)
       deallocate(rhop_HFB,rhon_HFB,kapp_HFB,kapn_HFB)
       deallocate(tran_p,tran_n)
       deallocate(trE0_p,trE0_n,trE1_p,trE1_n)
       deallocate(trE2_p,trE2_n,trE3_p,trE3_n)
       deallocate(trEN_p,trEN_n,trS1_p,trS1_n)
       deallocate(trM1s_p,trM1s_n,trM1l_p,trM1l_n)
       deallocate(trE1_p_dens,trE1_n_dens)
       deallocate(rad_den1)
       deallocate(zcross,weight)
       deallocate(H11p,H11n)
       deallocate(lp1,lp2,cg3)
       deallocate(itpoi1,klpoi1)
!       deallocate(V3BNO2)

       return
      end
