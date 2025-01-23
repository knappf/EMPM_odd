       subroutine Dp_field(D1)

       USE technical
       use v3body_no2b

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       double precision :: D1(id,id)

       type(twoquas_type), allocatable, save :: prho(:)
       integer :: iprho
       integer :: mpp(6),mtp(6)

       DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vpp_gen(:,:,:,:)
       DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vnn_gen(:,:,:,:)

       D1=0.d0

       if(ifp_hfb) then

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

       do i=1,id
        do j=1,id
         if(levp(i)%j2.eq.levp(j)%j2.and.levp(i)%l.eq.levp(j)%l) then
          val=0.d0
          do k=1,id
           do l=1,id
            val=val+Vpp_gen(k,l,i,j)
     &                *kapp_HFB(l,k)*dsqrt(dble(levp(k)%j2+1))
           enddo
          enddo
          D1(i,j)=0.5d0*val/dsqrt(dble(levp(i)%j2+1))
         endif
        enddo
       enddo

       deallocate(Vpp_gen,Vnn_gen)

       endif 

       return
      end
