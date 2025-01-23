      subroutine interaction

       USE technical
       USE geom
       use Tcm_2body
       USE v3body_no2b
       use NNN_no2b_bin

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       double precision :: unitar_matrix(id,id)

       integer(kind=1) :: ija,ita
       integer(kind=1) :: j2_abFF,j2_deFF,J2_3FF,it_abFF,it_deFF,IT2_3FF
       integer(kind=2) :: i1a,i2a,i3a,i4a,i5a,i6a
       real(kind=4) :: vint1,vint4

       character*2 :: ch

       id2=2*id

       allocate(tran_p(id,id),tran_n(id,id))
       tran_p=0.d0
       tran_n=0.d0

!       goto 1111
       allocate(Vpp(id,id,id,id,0:jmax))
       allocate(Vnn(id,id,id,id,0:jmax))
       allocate(Vpn(id,id,id,id,0:jmax))
       Vpp=0.d0
       Vnn=0.d0
       Vpn=0.d0

       ACM = dble(AZ+AN)

       open(unit=1,file='vlk.dat',status='old',form='unformatted')
        do while(.not.eof(1))
!         read(1,*) i1,ipar,Jc,i,j,k,l,vint1,vint2,vint3,vint4
         read(1) ita,ija,i1a,i2a,i3a,i4a,vint1,vint4
         i1=int(ita)
         Jc=int(ija)
         i=int(i1a)
         j=int(i2a)
         k=int(i3a)
         l=int(i4a)

!         vint=vint1-vint4*hbarom/ACM  !(dble(AZ+AN)+1.19d0)

          vint=vint1

!       2-body part of CM Hamiltonian

         if (ita.eq.0) 
     &   vint=vint1+hbarom/ACM*T_nn(i/2+1,j/2,k/2+1,l/2,Jc/2)
         if (ita.eq.1) 
     &   vint=vint1+hbarom/ACM*T_nn_asym(i/2,j/2,k/2,l/2,Jc/2)  
         if (ita.eq.-1) 
     &    vint=vint1+hbarom/ACM*T_nn_asym(i/2+1,j/2+1,k/2+1,l/2+1,Jc/2)  
      

         if(i1.eq.-1) vint=vint-vint1*(1.d0-quenp)
         if(i1.eq.1) vint=vint-vint1*(1.d0-quenn)
         if(i1.eq.0) vint=vint-vint1*(1.d0-quenpn)
         if(i1.eq.-1) then
          if((i.le.id2.and.j.le.id2).and.(k.le.id2.and.l.le.id2)) then
           if((if_2b.eq.1.and.(levp(i/2+1)%N+levp(j/2+1)%N.le.noscmax12
     &          .and.levp(k/2+1)%N+levp(l/2+1)%N.le.noscmax12)).or.
     &        (if_2b.eq.0.and.(levp(i/2+1)%N+levp(j/2+1)%N.le.2*noscmax
     &          .and.levp(k/2+1)%N+levp(l/2+1)%N.le.2*noscmax))) then
           factab=1.d0
           factcd=1.d0
           if(i.eq.j) factab=dsqrt(2.d0)
           if(k.eq.l) factcd=dsqrt(2.d0)
           xnorm=factab*factcd
           ij_phase=(-1)**((levp(i/2+1)%j2+levp(j/2+1)%j2)/2+Jc/2)
           kl_phase=(-1)**((levp(k/2+1)%j2+levp(l/2+1)%j2)/2+Jc/2)
           Vpp(i/2+1,j/2+1,k/2+1,l/2+1,Jc/2)=vint*xnorm
           Vpp(j/2+1,i/2+1,k/2+1,l/2+1,Jc/2)=-vint*ij_phase*xnorm
           Vpp(i/2+1,j/2+1,l/2+1,k/2+1,Jc/2)=-vint*kl_phase*xnorm
           Vpp(j/2+1,i/2+1,l/2+1,k/2+1,Jc/2)=vint*xnorm
     &                                         *ij_phase*kl_phase
           if(i.ne.k.or.j.ne.l) then
            Vpp(k/2+1,l/2+1,i/2+1,j/2+1,Jc/2)=vint*xnorm
            Vpp(l/2+1,k/2+1,i/2+1,j/2+1,Jc/2)=-vint*kl_phase*xnorm
            Vpp(k/2+1,l/2+1,j/2+1,i/2+1,Jc/2)=-vint*ij_phase*xnorm
            Vpp(l/2+1,k/2+1,j/2+1,i/2+1,Jc/2)=vint*xnorm
     &                                         *ij_phase*kl_phase
           endif
           endif
          endif
         endif
         if(i1.eq.1) then
          if((i.le.id2.and.j.le.id2).and.(k.le.id2.and.l.le.id2)) then
           if((if_2b.eq.1.and.(levp(i/2)%N+levp(j/2)%N.le.noscmax12
     &               .and.levp(k/2)%N+levp(l/2)%N.le.noscmax12)).or.
     &       (if_2b.eq.0.and.(levp(i/2)%N+levp(j/2)%N.le.2*noscmax
     &               .and.levp(k/2)%N+levp(l/2)%N.le.2*noscmax))) then
           factab=1.d0
           factcd=1.d0
           if(i.eq.j) factab=dsqrt(2.d0)
           if(k.eq.l) factcd=dsqrt(2.d0)
           xnorm=factab*factcd
           ij_phase=(-1)**((levn(i/2)%j2+levn(j/2)%j2)/2+Jc/2)
           kl_phase=(-1)**((levn(k/2)%j2+levn(l/2)%j2)/2+Jc/2)
           Vnn(i/2,j/2,k/2,l/2,Jc/2)=vint*xnorm
           Vnn(j/2,i/2,k/2,l/2,Jc/2)=-vint*ij_phase*xnorm
           Vnn(i/2,j/2,l/2,k/2,Jc/2)=-vint*kl_phase*xnorm
           Vnn(j/2,i/2,l/2,k/2,Jc/2)=vint*xnorm*ij_phase*kl_phase
           if(i.ne.k.or.j.ne.l) then
            Vnn(k/2,l/2,i/2,j/2,Jc/2)=vint*xnorm
            Vnn(l/2,k/2,i/2,j/2,Jc/2)=-vint*kl_phase*xnorm
            Vnn(k/2,l/2,j/2,i/2,Jc/2)=-vint*ij_phase*xnorm
            Vnn(l/2,k/2,j/2,i/2,Jc/2)=vint*xnorm*ij_phase*kl_phase
           endif
           endif
          endif
         endif
         if(i1.eq.0) then
          if((i.le.id2.and.j.le.id2).and.(k.le.id2.and.l.le.id2)) then
           if((if_2b.eq.1.and.(levp(i/2+1)%N+levp(j/2)%N.le.noscmax12
     &             .and.levp(k/2+1)%N+levp(l/2)%N.le.noscmax12)).or.
     &      (if_2b.eq.0.and.(levp(i/2+1)%N+levp(j/2)%N.le.2*noscmax
     &             .and.levp(k/2+1)%N+levp(l/2)%N.le.2*noscmax))) then
           Vpn(i/2+1,j/2,k/2+1,l/2,Jc/2)=vint
            if(i.ne.k.or.j.ne.l) then
             Vpn(k/2+1,l/2,i/2+1,j/2,Jc/2)=vint
            endif
           endif
          endif
         endif
        enddo
       close(1)

!***********************************************************************
!    Here the DD interaction elements are calculated                   *

       unitar_matrix=0.d0
       do i1=1,id
        unitar_matrix(i1,i1)=1.d0
       enddo

       igrid=20      ! number of the node points
       igrid2=250     ! number of the grid points

!       sizebox=10.d0    ! the interval in which the radial integral will be numerically summed
       sizebox=max(6.d0,2.0d0*dble(AZ+AN)**0.33333333d0)

       bos1=dsqrt(0.5d0*(zmp+zmn)*hbarom/(hbarc**2.d0))
       bos2=dsqrt(zmY*hbarom/(hbarc**2.d0))

       dx=sizebox/dble(igrid2)

       allocate(zcross(igrid))
       zcross(1)=0.070539889692d0 
       zcross(2)=0.372126818002d0
       zcross(3)=0.916582102483d0
       zcross(4)=1.70730653103d0
       zcross(5)=2.74919925531d0
       zcross(6)=4.04892531384d0
       zcross(7)=5.61517497087d0
       zcross(8)=7.45901745389d0
       zcross(9)=9.59439286749d0
       zcross(10)=12.0388025566d0
       zcross(11)=14.8142934155d0
       zcross(12)=17.9488955686d0
       zcross(13)=21.4787881904d0
       zcross(14)=25.4517028094d0
       zcross(15)=29.9325546634d0
       zcross(16)=35.0134341868d0
       zcross(17)=40.8330570974d0
       zcross(18)=47.6199940299d0
       zcross(19)=55.8107957541d0
       zcross(20)=66.5244165252d0
      
       do i=1,igrid
        zcross(i)=dsqrt(zcross(i))/bos1
       enddo

       allocate(weight(igrid))
       weight(1)=0.181080062419d0
       weight(2)=0.422556767879d0
       weight(3)=0.666909546702d0
       weight(4)=0.9153523727d0
       weight(5)=1.1695397071d0
       weight(6)=1.43135498624d0
       weight(7)=1.7029811359d0
       weight(8)=1.98701589585d0
       weight(9)=2.28663576323d0
       weight(10)=2.60583465152d0
       weight(11)=2.94978381794d0
       weight(12)=3.32539569477d0
       weight(13)=3.74225636246d0
       weight(14)=4.21424053477d0
       weight(15)=4.76252016007d0
       weight(16)=5.42172779036d0
       weight(17)=6.25401146407d0
       weight(18)=7.38731523837d0
       weight(19)=9.15132879607d0
       weight(20)=12.8933886244d0
  
       do i=1,igrid
        weight(i)=weight(i)/(2.d0*zcross(i)*bos1**2.d0)
       enddo

       allocate(rad_den1(igrid))
       rad_den1=0.d0

       do i=1,igrid
        val=0.d0
        radi=zcross(i)   !dble(i)*sizebox/dble(igrid)
        do j=1,id
         val_p=0.d0
         val_n=0.d0
         do k=1,id
          val_p=val_p+unitar_matrix(k,j)
     &                *R_val(levp(k)%nn,levp(k)%l,levp(k)%j2,bos1,radi)
          val_n=val_n+unitar_matrix(k,j)
     &                *R_val(levn(k)%nn,levn(k)%l,levn(k)%j2,bos1,radi)
         enddo
         val=val+lhfp(j)%vi**2.d0*val_p**2.d0*dble(lhfp(j)%j2+1)
         val=val+lhfn(j)%vi**2.d0*val_n**2.d0*dble(lhfn(j)%j2+1)
        enddo
        rad_den1(i)=val/(4.d0*pi)
       enddo

!   !***********************************************************************
!      Reading of 3-body NO2B elements from binary file generated by NuHamil code
       call read_3b_no2b_bin(noscmax,noscmax12,noscmax123)      




!***************************************************************************
!      Reading of 3-body elements

1111   continue 

       idim_v3b=0
       do it=1,3,2
       do i=0,jmax
        idim_v3b=idim_v3b+idim3jt(i,it)*(idim3jt(i,it)+1)/2
       enddo
       enddo

         ipozt1=0
         ipozjt1=0
       do jj=0,jmax
          ipozt1=ipozt1+(idim3jt(jj,1)*(idim3jt(jj,1)+1)/2)
          ipozjt1=ipozjt1+idim3jt(jj,1)
       enddo

!       idim_v3b=maxval(idim3j)
!       idim_v3b=int(idim_v3b,8)*(int(idim_v3b,8)+1)/2*(3*jmax+1)/2
       write(*,*)'Loading 3B m.e.'
!       idim_v3b=int(id3,8)*(int(id3,8)+1)/2

       write(*,*)'V3B array size (GB)',
     &  idim_v3b*4.d0/(1024.d0*1024.d0*1024.d0)
       allocate(V3B_ar(idim_v3b))
       V3B_ar=0.0                                      

        open(111,file='3belem_NO2.out',form='formatted',status='old')
        read(111,*)
        do while(.not.eof(111))
        read(111,'(10i4,f12.6)') ia,ib,ic,it_abFF,idd,ie,iff,it_deFF,
     &             j2_ab,IT2_3FF,xv3

!       write(234,'(10i4,f12.6)')ia,ib,ic,it_abFF,idd,ie,iff,it_deFF,
!     &             j2_ab,IT2_3FF,xv3


         Ni=lev1pn(ia)%N
         Nj=lev1pn(ib)%N
         Nk=lev1pn(ic)%N
         Nl=lev1pn(idd)%N
         Nm=lev1pn(ie)%N
         Nn=lev1pn(iff)%N

         if((Ni.le.noscmax.and.(Ni+Nj.le.noscmax12.and.Ni+Nj+Nk.le.
     &    noscmax123)).and.(Nl.le.noscmax.and.(Nl+Nm.le.noscmax12.and.
     &                                  Nl+Nm+Nn.le.noscmax123))) then

         ii1=lpoint(ia,ib,ic,j2_ab,it_abFF,IT2_3FF)
         ii2=lpoint(idd,ie,iff,j2_ab,it_deFF,IT2_3FF)

        if (ii1 == 0.or.ii2 == 0) then
          write(*,*) 'Error in 3-body matrix elements pointer'
          stop
        endif

         ipoz=0
         ipozjt=0

         if (IT2_3FF==1) then  
         do jj=0,j2_ab-1
          ipoz=ipoz+(idim3jt(jj,1)*(idim3jt(jj,1)+1)/2)
          ipozjt=ipozjt+idim3jt(jj,1)
         enddo
         endif 



         if (IT2_3FF==3) then  

!         do jj=0,jmax
!          ipoz=ipoz+(idim3jt(jj,1)*(idim3jt(jj,1)+1)/2)
!          ipozjt=ipozjt+idim3jt(jj,1)
!         enddo
          ipoz=ipozt1
          ipozjt=ipozjt1
         do jj=0,j2_ab-1
          ipoz=ipoz+(idim3jt(jj,3)*(idim3jt(jj,3)+1)/2)
          ipozjt=ipozjt+idim3jt(jj,3)
         enddo
         endif 

 
         ii1j=lev3(ii1)%ord-ipozjt
         ii2j=lev3(ii2)%ord-ipozjt

         if(ii1.ne.0.and.ii2.ne.0) then
         if (ii1j >= ii2j ) then
           V3B_ar((int(ii1j,8)*(int(ii1j,8)-1))/2+int(ii2j,8)+ipoz)=
     &real(xv3)
            else
           V3B_ar((int(ii2j,8)*(int(ii2j,8)-1))/2+int(ii1j,8)+ipoz)=
     &real(xv3)
         endif
!          V3B(ii1,ii2)=xv3*V3b_NNN
!          V3B(ii2,ii1)=xv3*V3b_NNN
         endif

!         if(j2_ab.le.jmax) then

!           if(klpoi1(ic,iff).ne.0.and.klpoi1(ia,idd).ne.0
!     &      .and.klpoi1(ib,ie).ne.0) then
!            if(itpoi1(it_abFF,it_deFF,IT2_3FF).ne.0) then

!            if ((ia>=ib).and.(idd>=ie)) then

!              V3BNO2(ia,ib,idd,ie,j2_ab,itpoi1(it_abFF,it_deFF,IT2_3FF),
!     &                                klpoi1(ic,iff))=xv3

!              V3BNO2(idd,ie,ia,ib,j2_ab,itpoi1(it_deFF,it_abFF,IT2_3FF),
!     &                             klpoi1(iff,ic))=xv3
!            endif 

!            if(idd.ne.ie.and.ia.eq.ib) then
!            if ((ia==ib).and.(idd/=ie)) then
             ph1=dble((-1)**(j2_ab-(lev1pn(idd)%j2+lev1pn(ie)%j2)/2))
             ph2=dble((-1)**(it_deFF))
!             V3BNO2(ia,ib,ie,idd,j2_ab,itpoi1(it_abFF,it_deFF,IT2_3FF),
!     &                      klpoi1(ic,iff))=xv3*ph2*ph1
!             V3BNO2(ie,idd,ia,ib,j2_ab,itpoi1(it_deFF,it_abFF,IT2_3FF),
!     &                      klpoi1(iff,ic))=xv3*ph2*ph1

!             endif

!            if ((ia/=ib).and.(idd==ie)) then
!             if(ia.ne.ib.and.idd.eq.ie) then
             ph1=dble((-1)**(j2_ab-(lev1pn(ia)%j2+lev1pn(ib)%j2)/2))
             ph2=dble((-1)**(it_abFF))
!             V3BNO2(ib,ia,idd,ie,j2_ab,itpoi1(it_abFF,it_deFF,IT2_3FF),
!     &                    klpoi1(ic,iff))=xv3*ph2*ph1

!             V3BNO2(idd,ie,ib,ia,j2_ab,itpoi1(it_deFF,it_abFF,IT2_3FF),
!     &                    klpoi1(iff,ic))=xv3*ph2*ph1

!            endif

!            if ((ia/=ib).and.(idd/=ie)) then
!             if(idd.ne.ie.and.ia.ne.ib) then
             ph1=dble((-1)**((lev1pn(ia)%j2+lev1pn(ib)%j2)/2
     &                    +(lev1pn(idd)%j2+lev1pn(ie)%j2)/2))
             ph2=dble((-1)**(it_abFF+it_deFF))
!             V3BNO2(ib,ia,ie,idd,j2_ab,itpoi1(it_abFF,it_deFF,IT2_3FF),
!     &                 klpoi1(ic,iff))=xv3*ph2*ph1

!             V3BNO2(ie,idd,ib,ia,j2_ab,itpoi1(it_deFF,it_abFF,IT2_3FF),
!     &                  klpoi1(iff,ic))=xv3*ph2*ph1

!            endif

!            endif
!            endif 


!         endif
!         endif
         endif
        enddo

!       allocate(cg3(0:2*jmax,-2*jmax:2*jmax,0:jmax,-jmax:3*jmax,
!     &                                                    0:3*jmax))
!       cg3=0.d0

!       do j1=0,2*jmax
!        do m1=-2*jmax,2*jmax
!         do j2=0,jmax
!          do m2=-jmax,jmax
!           do j3=0,3*jmax
!            cg3(j1,m1,j2,m2,j3)=cleb(j1,m1,j2,m2,j3,m1+m2)
!           enddo
!          enddo
!         enddo
!        enddo
!       enddo

!                                                                      *
!***********************************************************************

       return
      end
