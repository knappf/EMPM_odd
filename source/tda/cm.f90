!  last update 28.7. 2011

module cmconst
!
 contains 
!
!
 subroutine cmcreate(idph,iph,tbase)

  implicit none

  include 'types_tda_cp.inc'
 
  double precision, dimension (:,:,:), allocatable :: muel
  double precision, dimension (:,:), allocatable :: tbase
  type(ph_typ),dimension(:), allocatable :: iph
  integer :: i, j, ip, ih, it, idph
 
  call read_mul(muel)

  allocate(tbase(idph,idph))
  tbase=0.0d0

  do i=1,idph
      ip=iph(i)%par
      ih=iph(i)%hol
      it=iph(i)%tz
      if (it == 0 ) then
       if ((mod(ip,2)==0).and.(mod(ih,2)==0)) tbase(i,1)=muel(ip/2,ih/2,2)
       if ((mod(ip,2)==1).and.(mod(ih,2)==1)) tbase(i,1)=muel((ip+1)/2,(ih+1)/2,1)    
      endif 
  enddo


  do j=2,idph
   tbase(j,j)=1.0d0
  enddo

 end subroutine cmcreate
!
 subroutine read_mul(muel)
 
  implicit double precision (a-h,o-z)
   
  include 'formats_tda_cp.inc'


  double precision, dimension (:,:,:), allocatable :: muel
  character*30 :: fname 

  imax=200
  allocate(muel(imax,imax,2))
  muel=0.d0

  fname='r1Y1_p.dat'
  open(21,file=fname,status='unknown',form='formatted')

    do while (.not.eof(21))
      read(21,1002)i,j,rs
      muel(i,j,1)=rs
    enddo

    close(21)

  fname='r1Y1_n.dat'
  open(21,file=fname,status='unknown',form='formatted')

    do while (.not.eof(21))
      read(21,1002)i,j,rs
      muel(i,j,2)=rs
    enddo
     
    close(21)
 

 end subroutine read_mul
!
 subroutine ortog(ndim,tbase)
 
  implicit double precision (a-h,o-z)

  double precision, dimension (:,:), allocatable :: tbase 
  double precision, dimension (:), allocatable :: a
 
 
  allocate(a(ndim))
  a=0.0d0
! normalization of first vector      

      xnorm=0.d0

      do i=1,ndim
        xnorm=xnorm+tbase(i,1)**2.0d0
      enddo
      
      do i=1,ndim
        tbase(i,1)=tbase(i,1)/dsqrt(xnorm)
      end do

      do i=2,ndim
      
      do j=1,i-1
          xovrl=0.d0
        do k=1,ndim
          xovrl=xovrl+tbase(k,j)*tbase(k,i)
        end do
          a(j)=-1.d0*xovrl      
      end do
      
      do jj=1,i-1
         do ii=1,ndim
           tbase(ii,i)=tbase(ii,i)+a(jj)*tbase(ii,jj)
         end do
      end do
       xnorm=0.d0
      do ii=1,ndim
           xnorm=xnorm+tbase(ii,i)*tbase(ii,i)
      end do
      do ii=1,ndim
           tbase(ii,i)=tbase(ii,i)/dsqrt(xnorm)
      end do
      
      end do    



!      if (nf.ne.1) then
!      xnorm=0.d0
!      do i=1,no
!      do j=1,no
!      xnorm=xnorm+dd(i,1)*d1(i,j)*dd(j,1)
!      enddo
!      enddo
      
!      do i=1,no
!      dd(i,1)=dd(i,1)/dsqrt(xnorm)
!      end do

!      do i=2,no
      
!      do j=1,i-1
!      xovrl=0.d0
!      do k=1,no
!      do l=1,no
!      xovrl=xovrl+dd(k,j)*d1(k,l)*dd(l,i)
!      end do
!      end do
!      a(j)=-1.d0*xovrl      
!      end do
      
!      do jj=1,i-1
!      do ii=1,no
!      dd(ii,i)=dd(ii,i)+a(jj)*dd(ii,jj)
!      end do
!      end do
!      xnorm=0.d0
!      do ii=1,no
!      do jj=1,no
!      xnorm=xnorm+dd(ii,i)*d1(ii,jj)*dd(jj,i)
!      end do
!      end do
!      do ii=1,no
!      dd(ii,i)=dd(ii,i)/dsqrt(xnorm)
!      end do
      
!      end do    
!      endif 
           
      return 
      end subroutine ortog
!
!
 
     
end module cmconst







