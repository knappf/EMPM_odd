      subroutine HFB_energy

       USE technical
       use v3body_no2b

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'
       type(twoquas_type), allocatable, save :: prho(:)
       integer :: iprho
       integer :: mpp(6),mtp(6)

       DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vpp_gen(:,:,:,:)
       DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vnn_gen(:,:,:,:)

       E_kin=0.d0
       E_prot=0.d0
       E_neut=0.d0
       E_pn=0.d0
       E_pair=0.d0
       E_HFB=0.d0

       do i=1,id
        do j=1,id
         E_kin=E_kin+kin_p(i,j)*rhop_HFB(j,i)*dble(levp(j)%j2+1)
     &              +kin_n(i,j)*rhon_HFB(j,i)*dble(levn(j)%j2+1)
        enddo
       enddo

       do i=1,id
        do j=1,id
         if(levp(i)%j2.eq.levp(j)%j2) then
         do k=1,id
          do l=1,id
           if(levp(k)%j2.eq.levp(l)%j2) then
            do Jp=0,jmax
            E_prot=E_prot+0.5d0*Vpp(i,k,j,l,Jp)
     &                       *rhop_HFB(l,k)*rhop_HFB(j,i)*dble(2*Jp+1)
            E_neut=E_neut+0.5d0*Vnn(i,k,j,l,Jp)
     &                       *rhon_HFB(l,k)*rhon_HFB(j,i)*dble(2*Jp+1)
            E_pn=E_pn+1.d0*Vpn(i,k,j,l,Jp)
     &                       *rhop_HFB(j,i)*rhon_HFB(l,k)*dble(2*Jp+1)
            enddo
           endif
          enddo
         enddo
         endif
        enddo
       enddo


      E_prot_2b=E_prot
      E_neut_2b=E_neut
      E_pn_2b=E_pn

  
!    Radek's version

!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP& PRIVATE(i,j,l,m,n,Jp)
!$OMP DO REDUCTION(+:E_prot,E_neut,E_pn)
       do k=1,id
        do n=1,id
         if(klpoi1(k,n).ne.0) then
          do i=1,id
          do l=1,id
            if(klpoi1(i,l).ne.0) then 

!           if(lev1pn(i)%l.eq.lev1pn(l)%l.and.
!     &                         lev1pn(i)%j2.eq.lev1pn(l)%j2) then
            do j=1,id
             do Jp=abs(lev1pn(i)%j2-lev1pn(j)%j2)/2,
     &                             (lev1pn(i)%j2+lev1pn(j)%j2)/2
             do m=1,id

               if(klpoi1(j,m).ne.0) then

          E_prot=E_prot
     &    +(1.d0/6.d0)*rhop_HFB(lp1(l),lp1(i))
     &     *rhop_HFB(lp1(m),lp1(j))*rhop_HFB(lp1(n),lp1(k))
     &            *V3BNO2_me(i,j,k,l,m,n,Jp,1,1,3)

          E_neut=E_neut
     &    +(1.d0/6.d0)*rhon_HFB(lp1(l),lp1(i))
     &     *rhon_HFB(lp1(m),lp1(j))*rhon_HFB(lp1(n),lp1(k))
     &            *V3BNO2_me(i,j,k,l,m,n,Jp,1,1,3)

            E_pn=E_pn
     &    +(1.d0/3.d0)*rhop_HFB(lp1(l),lp1(i))
     &     *rhop_HFB(lp1(m),lp1(j))*rhon_HFB(lp1(n),lp1(k))
     &            *V3BNO2_me(i,j,k,l,m,n,Jp,1,1,1)
     &    +(1.d0/6.d0)*rhop_HFB(lp1(l),lp1(i))
     &     *rhop_HFB(lp1(m),lp1(j))*rhon_HFB(lp1(n),lp1(k))
     &            *V3BNO2_me(i,j,k,l,m,n,Jp,1,1,3)
     

     &    +(1.d0/6.d0)*(rhop_HFB(lp1(l),lp1(i))
     &     *rhon_HFB(lp1(m),lp1(j))*rhon_HFB(lp1(n),lp1(k)))*
     &     ((3.d0/2.d0)*V3BNO2_me(i,j,k,l,m,n,Jp,0,0,1)
     &    + (dsqrt(3.d0)/2.d0)
     &     *V3BNO2_me(i,j,k,l,m,n,Jp,1,0,1)
     &    + (dsqrt(3.d0)/2.d0)
     &     *V3BNO2_me(i,j,k,l,m,n,Jp,0,1,1)
     &    + (1.d0/2.d0)*V3BNO2_me(i,j,k,l,m,n,Jp,1,1,1)
     &    + V3BNO2_me(i,j,k,l,m,n,Jp,1,1,3)
     & )

              endif
             enddo
             enddo
            enddo
            endif
          enddo
          enddo
         endif
        enddo
       enddo
!$OMP END DO
!$OMP END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      Calculation of the pairing energy

       if(ifp_hfb.or.ifn_hfb) then

       allocate(Vpp_gen(id,id,id,id),Vnn_gen(id,id,id,id))
       Vpp_gen=0.d0
       Vnn_gen=0.d0

       do i=1,id
        do j=1,id
         do k=1,id
          do l=1,id
           valp=0.d0
           valn=0.d0

           do m=1,id
            do n=1,id

             valp=valp
     &+rhop_HFB(lp1(n),lp1(m))*V3BNO2_me(i,j,m,k,l,n,0,1,1,3)
     &+rhon_HFB(lp1(n),lp1(m))*V3BNO2_me(i,j,m,k,l,n,0,1,1,3)/3.d0
     &+rhon_HFB(lp1(n),lp1(m))*V3BNO2_me(i,j,m,k,l,n,0,1,1,1)*2.d0/3.d0

             valn=valn
     &+rhon_HFB(lp1(n),lp1(m))*V3BNO2_me(i,j,m,k,l,n,0,1,1,3)
     &+rhop_HFB(lp1(n),lp1(m))*V3BNO2_me(i,j,m,k,l,n,0,1,1,3)/3.d0
     &+rhop_HFB(lp1(n),lp1(m))*V3BNO2_me(i,j,m,k,l,n,0,1,1,1)*2.d0/3.d0

            enddo
           enddo

           Vpp_gen(lp1(i),lp1(j),lp1(k),lp1(l))=
     &                         Vpp(lp1(i),lp1(j),lp1(k),lp1(l),0)+valp
           Vnn_gen(lp1(i),lp1(j),lp1(k),lp1(l))=
     &                         Vnn(lp1(i),lp1(j),lp1(k),lp1(l),0)+valn
          enddo
         enddo
        enddo
       enddo

!       do i=1,id
!        do j=1,id
!         do k=1,id
!          do l=1,id
!           if(ifp_hfb) E_pair=E_pair+0.25d0
!     &                       *Vpp(i,j,k,l,0)
!     &                                 *kapp_HFB(i,j)*kapp_HFB(l,k)
!     &        *dsqrt(dble((levp(i)%j2+1)*(levp(k)%j2+1)))
!           if(ifn_hfb) E_pair=E_pair+0.25d0
!     &                       *Vnn(i,j,k,l,0)
!     &                                 *kapn_HFB(i,j)*kapn_HFB(l,k)
!     &        *dsqrt(dble((levn(i)%j2+1)*(levn(k)%j2+1)))
!          enddo
!         enddo
!        enddo
!       enddo
       do i=1,id
        do j=1,id
         do k=1,id
          do l=1,id
           if(ifp_hfb) E_pair=E_pair+0.25d0
     &                       *Vpp_gen(i,j,k,l)
     &                                 *kapp_HFB(i,j)*kapp_HFB(l,k)
     &        *dsqrt(dble((levp(i)%j2+1)*(levp(k)%j2+1)))
           if(ifn_hfb) E_pair=E_pair+0.25d0
     &                       *Vnn_gen(i,j,k,l)
     &                                 *kapn_HFB(i,j)*kapn_HFB(l,k)
     &        *dsqrt(dble((levn(i)%j2+1)*(levn(k)%j2+1)))
          enddo
         enddo
        enddo
       enddo

       deallocate(Vpp_gen,Vnn_gen)

       endif

       E_HFB = E_kin + E_prot + E_neut + E_pn + E_pair
       
!       write(*,*) 'HFB energy                 =',E_HFB,'  MeV'
!       write(*,*) 'Kinetic energy             =',E_kin,'  MeV'
!       write(*,*) 'Proton interaction energy  =',E_prot,'  MeV'
!       write(*,*) 'Neutron interaction energy =',E_neut,'  MeV'
!       write(*,*) 'P-N interaction energy     =',E_pn,'  MeV'
!       write(*,*) 'Pairing energy             =',E_pair,'  MeV'

       return
      end
