      subroutine transf_interaction
       
       USE technical
       use v3body_no2b

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       real(kind=4), ALLOCATABLE, SAVE :: Vpp_HFB(:,:,:,:,:)
       real(kind=4), ALLOCATABLE, SAVE :: Vnn_HFB(:,:,:,:,:)
       real(kind=4), ALLOCATABLE, SAVE :: Vpn_HFB(:,:,:,:,:)

       double precision :: Tp(id,id),Tn(id,id)
       double precision :: bp(id,id),bn(id,id)
       double precision :: timef,timein

       type(twoquas_type), allocatable, save :: prho(:)

       integer :: point_p(id),point_n(id)
       integer :: pinv_p(id),pinv_n(id)

       integer :: mpp(6),mtp(6)

!       Vpp=Vpp+3.d0*Vpp_DD
!       Vnn=Vnn+3.d0*Vnn_DD
!       Vpn=Vpn+3.d0*Vpn_DD

!    transformation to VNO2B elements

       allocate(Vpp_HFB(id,id,id,id,0:jmax))
       allocate(Vnn_HFB(id,id,id,id,0:jmax))
       allocate(Vpn_HFB(id,id,id,id,0:jmax))
       Vpp_HFB=0.0
       Vnn_HFB=0.0
       Vpn_HFB=0.0


!  based on Radek's formulae

!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP& PRIVATE(j,l,m,k,n,Jp,valp,valn,
!$OMP& valpn,v3b1,v3b2,v3b3,v3b4,v3b5,fact) 
!$OMP DO
       do i=1,id
       write(*,*) 'No2b part of residual interaction calc. i = '
     &  ,i,' of ',id
        do l=1,id
         do j=1,id
          do m=1,id

           do Jp=abs(lev1pn(i)%j2-lev1pn(j)%j2)/2,
     &                             (lev1pn(i)%j2+lev1pn(j)%j2)/2

            fact=float(2*Jp+1)

            valp=0.d0
            valn=0.d0
            valpn=0.d0

            do k=1,id
             do n=1,id

             if(klpoi1(k,n).ne.0) then
                   
          v3b1=V3BNO2_me(i,j,k,l,m,n,Jp,0,0,1)
          v3b2=V3BNO2_me(i,j,k,l,m,n,Jp,1,0,1)
          v3b3=V3BNO2_me(i,j,k,l,m,n,Jp,0,1,1)
          v3b4=V3BNO2_me(i,j,k,l,m,n,Jp,1,1,1)
          v3b5=V3BNO2_me(i,j,k,l,m,n,Jp,1,1,3)


          valp=valp+1.d0/dble(fact)*(
     &    (v3b5*rhop_HFB(lp1(n),lp1(k))+
     &     (2.d0/3.d0*v3b4+1.d0/3.d0*v3b5)
     &    *rhon_HFB(lp1(n),lp1(k)))) 

          valn=valn+1.d0/dble(fact)*(
     &    (v3b5*rhon_HFB(lp1(n),lp1(k))+
     &     (2.d0/3.d0*v3b4+1.d0/3.d0*v3b5)
     &    *rhop_HFB(lp1(n),lp1(k)))) 

          valpn=valpn+1.d0/dble(fact)*
     &    (1.d0/2.d0*v3b1-1.d0/(2.d0*dsqrt(3.d0))*v3b2 
     &    -1.d0/(2.d0*dsqrt(3.d0))*v3b3+1.d0/6.d0*v3b4
     &     +1.d0/3.d0*v3b5)
     &    *rhop_HFB(lp1(n),lp1(k))

           valpn=valpn+1.d0/dble(fact)*
     &    (1.d0/2.d0*v3b1+1.d0/(2.d0*dsqrt(3.d0))*v3b2 
     &    +1.d0/(2.d0*dsqrt(3.d0))*v3b3+1.d0/6.d0*v3b4
     &     +1.d0/3.d0*v3b5)
     &    *rhon_HFB(lp1(n),lp1(k))


             endif
             enddo  ! n 
           enddo   ! k

            Vpp_HFB(i,j,l,m,Jp)=real(valp,4)
            Vnn_HFB(i,j,l,m,Jp)=real(valn,4)
            Vpn_HFB(i,j,l,m,Jp)=real(valpn,4)

           enddo  ! Jp 
          enddo    ! m
         enddo   ! j
        enddo  !  l
       enddo  !  i

!$OMP END DO
!$OMP END PARALLEL

       do i=1,id
        do j=1,id
         do k=1,id
          do l=1,id
           do Jp=0,jmax
            Vpp(i,j,k,l,Jp)=Vpp(i,j,k,l,Jp)
     &           +Vpp_HFB(lp2(i),lp2(j),lp2(k),lp2(l),Jp)
            Vnn(i,j,k,l,Jp)=Vnn(i,j,k,l,Jp)
     &           +Vnn_HFB(lp2(i),lp2(j),lp2(k),lp2(l),Jp)
            Vpn(i,j,k,l,Jp)=Vpn(i,j,k,l,Jp)
     &           +Vpn_HFB(lp2(i),lp2(j),lp2(k),lp2(l),Jp)
           enddo
          enddo
         enddo
        enddo
       enddo

       deallocate(Vpp_HFB,Vnn_HFB,Vpn_HFB)
       deallocate(V3B_ar,V3B_pair)

       Tp=tran_p
       Tn=tran_n

       bp=0.d0
       bn=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         do k=1,id
          do l=1,id
           val1=val1+trE0_p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+trE0_n(k,l)*Tn(k,i)*Tn(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
        enddo
       enddo
       trE0_p=bp
       trE0_n=bn

       open(61,file='r2Y0_n.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
         if (dabs(bn(i,j)).gt.1.d-10) then
          write(61,'(2i5,5x,2f30.15)')i,j,bn(i,j)
         endif
        enddo
       enddo
       close(61)

       open(61,file='r2Y0_p.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
        if (dabs(bp(i,j)).gt.1.d-10) then
         write(61,'(2i5,5x,2f30.15)')i,j,bp(i,j)
        endif
        enddo
       enddo
       close(61)

       bp=0.d0
       bn=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         do k=1,id
          do l=1,id
           val1=val1+trE1_p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+trE1_n(k,l)*Tn(k,i)*Tn(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
        enddo
       enddo
       trE1_p=bp
       trE1_n=bn

       open(61,file='r1Y1_n.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
         if (dabs(bn(i,j)).gt.1.d-10) then
          write(61,'(2i5,5x,2f30.15)')i,j,bn(i,j)
         endif
        enddo
       enddo
       close(61)

       open(61,file='r1Y1_p.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
        if (dabs(bp(i,j)).gt.1.d-10) then
         write(61,'(2i5,5x,2f30.15)')i,j,bp(i,j)
        endif
        enddo
       enddo
       close(61)

!************************
!       open(61,file='E1p.out',status='unknown',form='formatted')
!        do i=1,id
!         write(61,*) (trE1_p(i,j),j=1,id)
!        enddo
!       close(61)
!       open(62,file='E1n.out',status='unknown',form='formatted')
!        do i=1,id
!         write(62,*) (trE1_n(i,j),j=1,id)
!        enddo
!       close(62)
!************************

       bp=0.d0
       bn=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         do k=1,id
          do l=1,id
           val1=val1+trE2_p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+trE2_n(k,l)*Tn(k,i)*Tn(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
        enddo
       enddo
       trE2_p=bp
       trE2_n=bn

       open(61,file='r2Y2_n.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
         if (dabs(bn(i,j)).gt.1.d-10) then
          write(61,'(2i5,5x,2f30.15)')i,j,bn(i,j)
         endif
        enddo
       enddo
       close(61)

       open(61,file='r2Y2_p.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
        if (dabs(bp(i,j)).gt.1.d-10) then
         write(61,'(2i5,5x,2f30.15)')i,j,bp(i,j)
        endif
        enddo
       enddo
       close(61)


       bp=0.d0
       bn=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         do k=1,id
          do l=1,id
           val1=val1+trE3_p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+trE3_n(k,l)*Tn(k,i)*Tn(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
        enddo
       enddo
       trE3_p=bp
       trE3_n=bn

       open(61,file='r3Y3_n.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
         if (dabs(bn(i,j)).gt.1.d-10) then
          write(61,'(2i5,5x,2f30.15)')i,j,bn(i,j)
         endif
        enddo
       enddo
       close(61)

       open(61,file='r3Y3_p.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
        if (dabs(bp(i,j)).gt.1.d-10) then
         write(61,'(2i5,5x,2f30.15)')i,j,bp(i,j)
        endif
        enddo
       enddo
       close(61)


       bp=0.d0
       bn=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         do k=1,id
          do l=1,id
           val1=val1+trEN_p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+trEN_n(k,l)*Tn(k,i)*Tn(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
        enddo
       enddo
       trEN_p=bp
       trEN_n=bn

       open(61,file='EN_n.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
         if (dabs(bn(i,j)).gt.1.d-10) then
          write(61,'(2i5,5x,2f30.15)')i,j,bn(i,j)
         endif
        enddo
       enddo
       close(61)

       open(61,file='EN_p.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
        if (dabs(bp(i,j)).gt.1.d-10) then
         write(61,'(2i5,5x,2f30.15)')i,j,bp(i,j)
        endif
        enddo
       enddo
       close(61)


       bp=0.d0
       bn=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         do k=1,id
          do l=1,id
           val1=val1+trS1_p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+trS1_n(k,l)*Tn(k,i)*Tn(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
        enddo
       enddo
       trS1_p=bp
       trS1_n=bn

       open(61,file='r3Y1_n.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
         if (dabs(bn(i,j)).gt.1.d-10) then
          write(61,'(2i5,5x,2f30.15)')i,j,bn(i,j)
         endif
        enddo
       enddo
       close(61)

       open(61,file='r3Y1_p.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
        if (dabs(bp(i,j)).gt.1.d-10) then
         write(61,'(2i5,5x,2f30.15)')i,j,bp(i,j)
        endif
        enddo
       enddo
       close(61)


       bp=0.d0
       bn=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         do k=1,id
          do l=1,id
           val1=val1+trM1s_p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+trM1s_n(k,l)*Tn(k,i)*Tn(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
        enddo
       enddo
       trM1s_p=bp
       trM1s_n=bn


       open(61,file='M1s_n.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
         if (dabs(bn(i,j)).gt.1.d-10) then
          write(61,'(2i5,5x,2f30.15)')i,j,bn(i,j)
         endif
        enddo
       enddo
       close(61)

       open(61,file='M1s_p.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
        if (dabs(bp(i,j)).gt.1.d-10) then
         write(61,'(2i5,5x,2f30.15)')i,j,bp(i,j)
        endif
        enddo
       enddo
       close(61)

       bp=0.d0
       bn=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         do k=1,id
          do l=1,id
           val1=val1+trM1l_p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+trM1l_n(k,l)*Tn(k,i)*Tn(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
        enddo
       enddo
       trM1l_p=bp
       trM1l_n=bn

       open(61,file='M1l_n.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
         if (dabs(bn(i,j)).gt.1.d-10) then
          write(61,'(2i5,5x,2f30.15)')i,j,bn(i,j)
         endif
        enddo
       enddo
       close(61)

       open(61,file='M1l_p.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
        if (dabs(bp(i,j)).gt.1.d-10) then
         write(61,'(2i5,5x,2f30.15)')i,j,bp(i,j)
        endif
        enddo
       enddo
       close(61)


       do m=0,igrid2
        bp=0.d0
        bn=0.d0
        do i=1,id
         do j=1,id
          val1=0.d0
          val2=0.d0
          do k=1,id
           do l=1,id
            val1=val1+trE1_p_dens(k,l,m)*Tp(k,i)*Tp(l,j)
            val2=val2+trE1_n_dens(k,l,m)*Tn(k,i)*Tn(l,j)
           enddo
          enddo
          bp(i,j)=val1
          bn(i,j)=val2
         enddo
        enddo
        do i=1,id
         do j=1,id
          trE1_p_dens(i,j,m)=bp(i,j)
          trE1_n_dens(i,j,m)=bn(i,j)
         enddo
        enddo
       enddo

       bp=0.d0
       bn=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         do k=1,id
          do l=1,id
           val1=val1+H11p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+H11n(k,l)*Tn(k,i)*Tn(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
        enddo
       enddo
       H11p=bp
       H11n=bn

       do i=1,id
        do j=1,id
         H11p(i,j)=-(lhfp(i)%ui*lhfp(j)%vi+lhfp(i)%vi*lhfp(j)%ui)
     &                                                  *H11p(i,j)
         H11n(i,j)=-(lhfn(i)%ui*lhfn(j)%vi+lhfn(i)%vi*lhfn(j)%ui)
     &                                                  *H11n(i,j)
        enddo
       enddo
       do i=1,id
        H11p(i,i)=H11p(i,i)+(lhfp(i)%ui**2.d0-lhfp(i)%vi**2.d0)
     &                                        *(lhfp(i)%ei-ferp)
        H11n(i,i)=H11n(i,i)+(lhfn(i)%ui**2.d0-lhfn(i)%vi**2.d0)
     &                                        *(lhfn(i)%ei-fern)
       enddo

       open(4,file='H11p.dat',status='unknown',form='formatted')
       write(4,*) ferp
        do i1=1,id
         write(4,*) (H11p(i1,j1),j1=1,id) 
        enddo
       close(4)
       open(4,file='H11n.dat',status='unknown',form='formatted')
       write(4,*) fern
        do i1=1,id
         write(4,*) (H11n(i1,j1),j1=1,id) 
        enddo
       close(4)

       allocate(Vpp_HFB(id,id,id,id,0:jmax))
       allocate(Vnn_HFB(id,id,id,id,0:jmax))
       allocate(Vpn_HFB(id,id,id,id,0:jmax))
       Vpp_HFB=0.0
       Vnn_HFB=0.0
       Vpn_HFB=0.0


  
      write(*,*) 'Transf. interaction in 1st index'
!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP& PRIVATE(Jp,i1,i2,j1,k1,l1,valp,valn,val0)
!$OMP DO
       do Jp=0,jmax
        do i1=1,id
         do j1=1,id
          do k1=1,id
           do l1=1,id
            valp=0.d0
            valn=0.d0
            val0=0.d0
            do i2=1,id
             valp=valp+Vpp(i2,j1,k1,l1,Jp)*Tp(i2,i1)
             valn=valn+Vnn(i2,j1,k1,l1,Jp)*Tn(i2,i1)
             val0=val0+Vpn(i2,j1,k1,l1,Jp)*Tp(i2,i1)
            enddo
            Vpp_HFB(i1,j1,k1,l1,Jp)=real(valp,4)
            Vnn_HFB(i1,j1,k1,l1,Jp)=real(valn,4)
            Vpn_HFB(i1,j1,k1,l1,Jp)=real(val0,4)
           enddo
          enddo
         enddo
        enddo
       enddo
!$OMP END DO
!$OMP END PARALLEL

       Vpp=Vpp_HFB
       Vnn=Vnn_HFB
       Vpn=Vpn_HFB


      write(*,*) 'Transf. interaction in 2nd index'
!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP& PRIVATE(Jp,i1,j1,j2,k1,l1,valp,valn,val0)
!$OMP DO
       do Jp=0,jmax
        do i1=1,id
         do j1=1,id
          do k1=1,id
           do l1=1,id
            valp=0.d0
            valn=0.d0
            val0=0.d0
            do j2=1,id
             valp=valp+Vpp(i1,j2,k1,l1,Jp)*Tp(j2,j1)
             valn=valn+Vnn(i1,j2,k1,l1,Jp)*Tn(j2,j1)
             val0=val0+Vpn(i1,j2,k1,l1,Jp)*Tn(j2,j1)
            enddo
            Vpp_HFB(i1,j1,k1,l1,Jp)=real(valp,4)
            Vnn_HFB(i1,j1,k1,l1,Jp)=real(valn,4)
            Vpn_HFB(i1,j1,k1,l1,Jp)=real(val0,4)
           enddo
          enddo
         enddo
        enddo
       enddo
!$OMP END DO
!$OMP END PARALLEL

       Vpp=Vpp_HFB
       Vnn=Vnn_HFB
       Vpn=Vpn_HFB


       write(*,*) 'Transf. interaction in 3rd index'

!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP& PRIVATE(Jp,i1,j1,k1,k2,l1,valp,valn,val0)
!$OMP DO       
       do Jp=0,jmax
        do i1=1,id
         do j1=1,id
          do k1=1,id
           do l1=1,id
            valp=0.d0
            valn=0.d0
            val0=0.d0
            do k2=1,id
             valp=valp+Vpp(i1,j1,k2,l1,Jp)*Tp(k2,k1)
             valn=valn+Vnn(i1,j1,k2,l1,Jp)*Tn(k2,k1)
             val0=val0+Vpn(i1,j1,k2,l1,Jp)*Tp(k2,k1)
            enddo
            Vpp_HFB(i1,j1,k1,l1,Jp)=real(valp,4)
            Vnn_HFB(i1,j1,k1,l1,Jp)=real(valn,4)
            Vpn_HFB(i1,j1,k1,l1,Jp)=real(val0,4)
           enddo
          enddo
         enddo
        enddo
       enddo
!$OMP END DO
!$OMP END PARALLEL


       Vpp=Vpp_HFB
       Vnn=Vnn_HFB
       Vpn=Vpn_HFB

       write(*,*) 'Transf. interaction in 4th index'

!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP& PRIVATE(Jp,i1,j1,k1,l1,l2,valp,valn,val0)
!$OMP DO   

       do Jp=0,jmax
        do i1=1,id
         do j1=1,id
          do k1=1,id
           do l1=1,id
            valp=0.d0
            valn=0.d0
            val0=0.d0
            do l2=1,id
             valp=valp+Vpp(i1,j1,k1,l2,Jp)*Tp(l2,l1)
             valn=valn+Vnn(i1,j1,k1,l2,Jp)*Tn(l2,l1)
             val0=val0+Vpn(i1,j1,k1,l2,Jp)*Tn(l2,l1)
            enddo
            Vpp_HFB(i1,j1,k1,l1,Jp)=real(valp,4)
            Vnn_HFB(i1,j1,k1,l1,Jp)=real(valn,4)
            Vpn_HFB(i1,j1,k1,l1,Jp)=real(val0,4)
           enddo
          enddo
         enddo
        enddo
       enddo

!$OMP END DO
!$OMP END PARALLEL


       Vpp=Vpp_HFB
       Vnn=Vnn_HFB
       Vpn=Vpn_HFB

       open(1,file='vlk_hfb.dat',status='unknown',form='unformatted')
       open(2,file='int_HF_form.dat',status='unknown',form='formatted')
       write(2,*) '  Tz    J    a    b    c    d     <ab,J|V|cd,J> '
       do Jp=0,jmax
        do i1=1,id
         do i2=1,id
          do i3=1,id
           do i4=1,id
            i=2*i1-1
            j=2*i2-1
            k=2*i3-1
            l=2*i4-1
            if(i.le.j.and.k.le.l) then
             if(1000*i+j.le.1000*k+l) then
              factab=1.d0
              factcd=1.d0
              if(i.eq.j) factab=dsqrt(2.d0)
              if(k.eq.l) factcd=dsqrt(2.d0)
              xnorm=factab*factcd
              if(abs(Vpp_HFB(i1,i2,i3,i4,Jp)).gt.precis) then 
              write(1) 
     &   int(-1,1),int(2*Jp,1), int(i,2),int(j,2),int(k,2),int(l,2)
     &   ,real(Vpp_HFB(i1,i2,i3,i4,Jp),8)/xnorm
!              write(2,'(6i5,f10.5)')
!     &   -1,Jp,i1,i2,i3,i4,Vpp_HFB(i1,i2,i3,i4,Jp)
              endif
             endif
            endif
           enddo
          enddo
         enddo
        enddo
       enddo
       do Jp=0,jmax
        do i1=1,id
         do i2=1,id
          do i3=1,id
           do i4=1,id
            i=2*i1-1
            j=2*i2
            k=2*i3-1
            l=2*i4
!            if(i.le.j.and.k.le.l) then
             if(1000*i+j.le.1000*k+l) then
              if(abs(Vpn_HFB(i1,i2,i3,i4,Jp)).gt.precis) then 
              write(1) 
     &   int(0,1),int(2*Jp,1), int(i,2),int(j,2),int(k,2),int(l,2)
     &   ,real(Vpn_HFB(i1,i2,i3,i4,Jp),8)
!              write(2,'(6i5,f10.5)')
!     &   0,Jp,i1,i2,i3,i4,Vpn_HFB(i1,i2,i3,i4,Jp)

     
             endif
             endif 
!            endif
           enddo
          enddo
         enddo
        enddo
       enddo
       do Jp=0,jmax
        do i1=1,id
         do i2=1,id
          do i3=1,id
           do i4=1,id
            i=2*i1
            j=2*i2
            k=2*i3
            l=2*i4
            if(i.le.j.and.k.le.l) then
             if(1000*i+j.le.1000*k+l) then
              factab=1.d0
              factcd=1.d0
              if(i.eq.j) factab=dsqrt(2.d0)
              if(k.eq.l) factcd=dsqrt(2.d0)
              xnorm=factab*factcd
              if(abs(Vnn_HFB(i1,i2,i3,i4,Jp)).gt.precis) then
              write(1)
     &   int(1,1),int(2*Jp,1), int(i,2),int(j,2),int(k,2),int(l,2)
     &   ,real(Vnn_HFB(i1,i2,i3,i4,Jp),8)/xnorm
!              write(2,'(6i5,f10.5)')
!     &   1,Jp,i1,i2,i3,i4,Vnn_HFB(i1,i2,i3,i4,Jp)

             endif
             endif
            endif
           enddo
          enddo
         enddo
        enddo
       enddo
       close(1)

!****************************************************************************
!      The many-body perturbation theory E(2) energy is calculated here     *
       ener2pt=0.d0
       ener2ptY=0.d0
       if((.not.ifp_hfb).and.(.not.ifn_hfb)) then
        Vpp=Vpp_HFB
        Vnn=Vnn_HFB
        Vpn=Vpn_HFB
        call MBPT_energy(ener2pt)

        open(1,file='Energy_2.out',status='unknown',form='formatted')
         write(1,*) 'E_HF=',Energy_HF,' MeV'
         write(1,*) 'E^(2)=',ener2pt,' MeV'
         write(1,*) 'E_tot=',Energy_HF+ener2pt,' MeV'
        close(1)
       endif
!****************************************************************************


!****************************************************************************

      open(1,file='Summary.out',status='unknown',form='formatted')
        write(1,*) ' Single particle spectrum'
        write(1,*) ' ----------------------------------------------'
        
        write(1,*)'      neutrons             protons '
        write(1,*) '  i   l  2j  e [MeV]        l  2j  e [MeV] '
        do ii = 1, id
        write(1,'(3i4,f10.5,3x,2i4,f10.5)') 
     &    lhfn(ii)%index,lhfn(ii)%l,lhfn(ii)%j2,
     &    lhfn(ii)%ei,
     &    lhfp(ii)%l,lhfp(ii)%j2,
     &    lhfp(ii)%ei
        end do
        write(1,*) ' ----------------------------------------------'
        write(*,*)
        write(1,*) 'E_HF=',Energy_HF,' MeV'
        write(1,*) 'E^(2)=',ener2pt,' MeV'
        write(1,*) 'E_tot=',Energy_HF+ener2pt,' MeV'
        write(1,*) ' ----------------------------------------------'
        write(1,*) 'r =',dsqrt(r4tot/r2tot),'fm'
        write(1,*) 'r_p=',dsqrt(r4Z/r2Z),'fm'
        write(1,*) 'r_n=',dsqrt(r4N/r2N),'fm'

       close(1)

!****************************************************************************

       do ll=0,(jmax+1)/2
        do jj=1,jmax,2
         nn=-1
         do i=1,id
          if(lhfp(i)%j2.eq.jj.and.lhfp(i)%l.eq.ll) then
           nn=nn+1
           lhfp(i)%nn=nn
          endif
         enddo
        enddo
       enddo

       do ll=0,(jmax+1)/2
        do jj=1,jmax,2
         nn=-1
         do i=1,id
          if(lhfn(i)%j2.eq.jj.and.lhfn(i)%l.eq.ll) then
           nn=nn+1
           lhfn(i)%nn=nn
          endif
         enddo
        enddo
       enddo

       i1p=0
       i1n=0
       do i=1,id
        if(lhfp(i)%vi**2.d0.gt.0.98) then
         i1p=i1p+1
         point_p(i)=id+1-i1p
         pinv_p(id+1-i1p)=i
        else
         point_p(i)=i-i1p
         pinv_p(i-i1p)=i
        endif
        if(lhfn(i)%vi**2.d0.gt.0.98) then
         i1n=i1n+1
         point_n(i)=id+1-i1n
         pinv_n(id+1-i1n)=i
        else
         point_n(i)=i-i1n
         pinv_n(i-i1n)=i
        endif
       enddo

       open(2,file='proton_HF.dat',status='unknown',form='formatted')
        write(2,'(1x,a61)')'n,     l,    2*j,      Tz,     qei,
     &  ui,      vi'
        itz=-1
        do ii = 1, id
         jj=ii !pinv_p(ii)
         write(2,'(1x,4(i4,1x),3(f12.5,1x))') lhfp(jj)%nn,lhfp(jj)%l,lhf
     &p(jj)%j2,itz,lhfp(jj)%qei,DSQRT(1.D0-lhfp(jj)%vi**2.d0),
     &lhfp(jj)%vi
        end do
       close(2)
       open(2,file='neutron_HF.dat',status='unknown',form='formatted')
        write(2,'(1x,a61)')'n,     l,    2*j,      Tz,     qei,
     &  ui,      vi'
        itz=1
        do ii = 1, id
         jj=ii !pinv_n(ii)
         write(2,'(1x,4(i4,1x),3(f12.5,1x))') lhfn(jj)%nn,lhfn(jj)%l,lhf
     &n(jj)%j2,itz,lhfn(jj)%qei,DSQRT(1.D0-lhfn(jj)%vi**2.d0),
     &lhfn(jj)%vi
        end do
       close(2)

!       open(3,file='gmat_HF_pp.dat',status='unknown',form='formatted')
!       itz=-2
!        do Jp=0,jmax
!         do i1=1,id
!          do i2=1,id
!           do i3=1,id
!            do i4=1,id
!             i=i1 !pinv_p(i1)
!             j=i2 !pinv_p(i2)
!             k=i3 !pinv_p(i3)
!             l=i4 !pinv_p(i4)
!             if(i.le.j.and.k.le.l) then
!              if(1000*i+j.le.1000*k+l) then
!               if(dabs(Vpp_HFB(i,j,k,l,Jp)).gt.precis) then
!                factab=1.d0
!                factcd=1.d0
!                if(i.eq.j) factab=dsqrt(2.d0)
!                if(k.eq.l) factcd=dsqrt(2.d0)
!                xnorm=factab*factcd
!                write(3,'(1x,14(i4,1x),1(f12.5,1x))') 
!     &            lhfp(i)%nn,lhfp(i)%l,lhfp(i)%j2,
!     &            lhfp(j)%nn,lhfp(j)%l,lhfp(j)%j2,
!     &            lhfp(k)%nn,lhfp(k)%l,lhfp(k)%j2,
!     &            lhfp(l)%nn,lhfp(l)%l,lhfp(l)%j2,
!     &            itz,2*Jp,Vpp_HFB(i,j,k,l,Jp)/xnorm
!               endif
!              endif
!             endif
!            enddo
!           enddo
!          enddo
!         enddo
!        enddo
!       close(3)

!       open(3,file='gmat_HF_nn.dat',status='unknown',form='formatted')
!       itz=2
!        do Jp=0,jmax
!         do i1=1,id
!          do i2=1,id
!           do i3=1,id
!            do i4=1,id
!             i=i1 !pinv_n(i1)
!             j=i2 !pinv_n(i2)
!             k=i3 !pinv_n(i3)
!             l=i4 !pinv_n(i4)
!             if(i.le.j.and.k.le.l) then
!              if(1000*i+j.le.1000*k+l) then
!               if(dabs(Vnn_HFB(i,j,k,l,Jp)).gt.precis) then
!                factab=1.d0
!                factcd=1.d0
!                if(i.eq.j) factab=dsqrt(2.d0)
!                if(k.eq.l) factcd=dsqrt(2.d0)
!                xnorm=factab*factcd
!                write(3,'(1x,14(i4,1x),1(f12.5,1x))')
!     &            lhfn(i)%nn,lhfn(i)%l,lhfn(i)%j2,
!     &            lhfn(j)%nn,lhfn(j)%l,lhfn(j)%j2,
!     &            lhfn(k)%nn,lhfn(k)%l,lhfn(k)%j2,
!     &            lhfn(l)%nn,lhfn(l)%l,lhfn(l)%j2,
!     &            itz,2*Jp,Vnn_HFB(i,j,k,l,Jp)/xnorm
!               endif
!              endif
!             endif
!            enddo
!           enddo
!          enddo
!         enddo
!        enddo
!       close(3)

!       open(3,file='gmat_HF_pn.dat',status='unknown',form='formatted')
!       itz=0
!        do Jp=0,jmax
!         do i1=1,id
!          do i2=1,id
!           do i3=1,id
!            do i4=1,id
!             i=i1 !pinv_p(i1)
!             j=i2 !pinv_n(i2)
!             k=i3 !pinv_p(i3)
!             l=i4 !pinv_n(i4)
!             if(1000*i+j.le.1000*k+l) then
!              if(dabs(Vpn_HFB(i,j,k,l,Jp)).gt.precis) then
!               write(3,'(1x,14(i4,1x),1(f12.5,1x))')
!     &           lhfp(i)%nn,lhfp(i)%l,lhfp(i)%j2,
!     &           lhfn(j)%nn,lhfn(j)%l,lhfn(j)%j2,
!     &           lhfp(k)%nn,lhfp(k)%l,lhfp(k)%j2,
!     &           lhfn(l)%nn,lhfn(l)%l,lhfn(l)%j2,
!     &           itz,2*Jp,Vpn_HFB(i,j,k,l,Jp)
!              endif
!             endif
!            enddo
!           enddo
!          enddo
!         enddo
!        enddo
!       close(3)

       open(4,file='unitar_p.dat',status='unknown',form='formatted')
        do i1=1,id
         write(4,*) (Tp(i1,j1),j1=1,id) !(Tp(i1,pinv_p(j1)),j1=1,id)
        enddo
       close(4)

       open(4,file='unitar_n.dat',status='unknown',form='formatted')
        do i1=1,id
         write(4,*) (Tn(i1,j1),j1=1,id) !(Tn(i1,pinv_n(j1)),j1=1,id)
        enddo
       close(4)

       deallocate(Vpp_HFB,Vnn_HFB,Vpn_HFB)

       return
      end
