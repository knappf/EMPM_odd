      Subroutine hn_field(h1)

       USE technical
       use v3body_no2b

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'
   
       double precision :: h1(id,id)
       double precision :: h2(id,id), h3(id,id)
       type(twoquas_type), allocatable, save :: prho(:)
       integer :: iprho
       integer :: mpp(6),mtp(6)

       h1=0.d0
       h2=0.d0
       h3=0.d0

       do i=1,id
        do j=1,id
         val=0.d0
         do k=1,id
          do l=1,id
           if(levn(i)%j2.eq.levn(j)%j2
     &                .and.levn(k)%j2.eq.levn(l)%j2) then
            do Jp=0,jmax
            val=val+Vnn(i,k,j,l,Jp)
     &              *rhon_HFB(l,k)*dble(2*Jp+1)/dble(levn(i)%j2+1)
     &             +Vpn(k,i,l,j,Jp)
     &              *rhop_HFB(l,k)*dble(2*Jp+1)/dble(levn(i)%j2+1)
            enddo
           endif
          enddo
         enddo
         h1(i,j)=kin_n(i,j)+val
        enddo
       enddo

!    Petr's version
!       do i=1,id
!        do j=1,id
!         val2=0.d0
!         if(klpoi1(i,j).ne.0) then
!          do k=1,id
!          do m=1,id
!           if(lev1pn(k)%l.eq.lev1pn(m)%l.and.
!     &                         lev1pn(k)%j2.eq.lev1pn(m)%j2) then
!            do l=1,id
!             do Jp=abs(lev1pn(k)%j2-lev1pn(l)%j2)/2,
!     &                             (lev1pn(k)%j2+lev1pn(l)%j2)/2
!             do n=1,id
!              if(lev1pn(l)%l.eq.lev1pn(n)%l.and.
!     &                          lev1pn(l)%j2.eq.lev1pn(n)%j2) then

!               val2=val2
!     &          +0.5d0*rhon_HFB(lp1(m),lp1(k))*rhon_HFB(lp1(n),lp1(l))
!     &            *V3BNO2(k,l,m,n,Jp,itpoi1(1,1,3),klpoi1(i,j))
!     &            /dble(lev1pn(i)%j2+1)
!     &     +(1.d0/3.d0)*rhop_HFB(lp1(m),lp1(k))*rhop_HFB(lp1(n),lp1(l))
!     &            *V3BNO2(k,l,m,n,Jp,itpoi1(1,1,1),klpoi1(i,j))
!     &            /dble(lev1pn(i)%j2+1)
!     &     +(1.d0/6.d0)*rhop_HFB(lp1(m),lp1(k))*rhop_HFB(lp1(n),lp1(l))
!     &            *V3BNO2(k,l,m,n,Jp,itpoi1(1,1,3),klpoi1(i,j))
!     &            /dble(lev1pn(i)%j2+1)
!     &     +(1.d0/2.d0)*rhop_HFB(lp1(m),lp1(k))*rhon_HFB(lp1(n),lp1(l))
!     &            *V3BNO2(k,l,m,n,Jp,itpoi1(0,0,1),klpoi1(i,j))
!     &            /dble(lev1pn(i)%j2+1)
!     &     +(1.d0/6.d0)*rhop_HFB(lp1(m),lp1(k))*rhon_HFB(lp1(n),lp1(l))
!     &            *V3BNO2(k,l,m,n,Jp,itpoi1(1,1,1),klpoi1(i,j))
!     &            /dble(lev1pn(i)%j2+1)
!     &     +(1.d0/3.d0)*rhop_HFB(lp1(m),lp1(k))*rhon_HFB(lp1(n),lp1(l))
!     &            *V3BNO2(k,l,m,n,Jp,itpoi1(1,1,3),klpoi1(i,j))
!     &            /dble(lev1pn(i)%j2+1)
!     &     +(1.d0/(2.d0*dsqrt(3.d0)))

!     &             *rhop_HFB(lp1(m),lp1(k))*rhon_HFB(lp1(n),lp1(l))
!     &            *V3BNO2(k,l,m,n,Jp,itpoi1(0,1,1),klpoi1(i,j))
!     &            /dble(lev1pn(i)%j2+1)
!     &     +(1.d0/(2.d0*dsqrt(3.d0)))
!     &             *rhop_HFB(lp1(m),lp1(k))*rhon_HFB(lp1(n),lp1(l))
!     &            *V3BNO2(k,l,m,n,Jp,itpoi1(1,0,1),klpoi1(i,j))
!     &            /dble(lev1pn(i)%j2+1)
!              endif
!             enddo
!             enddo
!            enddo
!           endif
!          enddo
!          enddo
!         endif
!         h2(i,j)=val2
!        enddo
!       enddo


!      formulae according to Radek Folprecht
      do i=1,id
        do l=1,id
         if(klpoi1(i,l).ne.0) then
         val2=0.d0

!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP& PRIVATE(j,m,k,n,Jp,v3b1,v3b2,v3b3,v3b4,v3b5)
!$OMP DO REDUCTION(+:val2)

         do j=1,id
          do m=1,id
            if(klpoi1(j,m).ne.0) then

            do k=1,id
             do n=1,id

             if(klpoi1(k,n).ne.0) then

           do Jp=abs(lev1pn(i)%j2-lev1pn(j)%j2)/2,
     &                             (lev1pn(i)%j2+lev1pn(j)%j2)/2


          v3b1=V3BNO2_me(i,j,k,l,m,n,Jp,0,0,1)
          v3b2=V3BNO2_me(i,j,k,l,m,n,Jp,1,0,1)
          v3b3=V3BNO2_me(i,j,k,l,m,n,Jp,0,1,1)
          v3b4=V3BNO2_me(i,j,k,l,m,n,Jp,1,1,1)
          v3b5=V3BNO2_me(i,j,k,l,m,n,Jp,1,1,3)

          val2=val2+1.d0/dble(lev1pn(i)%j2+1)*(
     &    (0.25d0*(v3b1
     &     +1.d0/dsqrt(3.d0)*v3b2+1.d0/dsqrt(3.d0)*v3b3
     &     +1.d0/3.d0*v3b4+2.d0/3.d0*v3b5)*         
     &    rhop_HFB(lp1(m),lp1(j))*rhop_HFB(lp1(n),lp1(k))
     &    )

     &    +(1.d0/3.d0*(2.d0*v3b4
     &     +v3b5)*         
     &    rhon_HFB(lp1(m),lp1(j))*rhop_HFB(lp1(n),lp1(k))
     &    )
     &    +0.5d0*v3b5*
     &    rhon_HFB(lp1(m),lp1(j))*rhon_HFB(lp1(n),lp1(k))
     &    )

             enddo  
             endif 
           enddo   ! k
           enddo  ! m 
           endif          
          enddo    ! m
         enddo   ! j

!$OMP END DO
!$OMP END PARALLEL
!         endif
         h3(i,l)=val2
         endif
        enddo  !  l
       enddo  !  i 

 
       do i=1,id
        do j=1,id
         h1(i,j)=h1(i,j)+h3(lp2(i),lp2(j))
        enddo
       enddo

       return
      end
