!c     last modification 11.5.2018

module hami

contains

subroutine ham_geev(ndim,ndimr,no,nor,ns,irow,wr,xr,hami,ipar,jcal,itzcal,mxtr,phonbs,nx)

use choleski

implicit double precision (a-h,o-z)
   
include 'types_eqm.inc'   
include 'formats_eqm.inc'

type(phonbase_typ), dimension (:), allocatable :: phonbs
double precision, dimension(:,:), allocatable :: d1,amatr,hami,hamid,dmatr,cq,vr,d1r,hamir,hamidr,hamd,xr
double precision, dimension(:,:), allocatable :: d_orig,h_orig,h_corr
 
double precision, dimension(:), allocatable ::  work,wr,wi,wro

integer, dimension(:), allocatable :: mxt,nxt,mxtr,ipoz,irow,nx,ifail

double precision, dimension (:), allocatable :: s
double precision, dimension (:,:), allocatable :: vt
integer, dimension (:), allocatable :: ind_red
logical :: decoup_cm



allocate(d1(ndimr,ndim))
d1=0.d0

allocate(amatr(ndimr,ndim))
amatr=0.d0

ndimtotal=ndim


open(6,file='d_mat.dat',status='old',form='unformatted')
read(6)ndimrt,ndimt
if (ndimrt.ne.ndimr.or.ndimt.ne.ndim) then
  write(*,*)' Dimensions does not match in D_m file'
  stop
endif

do iii=1,ndimrt
  read(6)(d1(iii,jjj),jjj=1,ndimt)
enddo
      
close(6)

open(6,file='a_mat.dat',status='old',form='unformatted')

read(6)ndimrt,ndimt
if (ndimrt.ne.ndimr.or.ndimt.ne.ndim) then
  write(*,*)' Dimensions does not match in A_m file'
  stop
endif

do iii=1,ndimrt
 read(6)(amatr(iii,jjj),jjj=1,ndimt)
enddo

close(6)

write(*,*)'Dimensions :',ndim,ndimr

allocate(hami(ndimr,ndimr))
hami=0.d0

!  AD matrix 
call dgemm('N','T',ndimr,ndimr,ndim,1.d0,amatr,ndimr,d1,ndimr,0.d0,hami,ndimr)

iout=0
if (iout.eq.1) then 
  write(998,*)    
  write(998,*)'******** matrix  AD *************'       
  write(998,*)
do i=1,ndimr
  write(998,102)(hami(i,j),j=1,ndimr)
enddo
endif

do i=1,ndimr
  do j=1,ndimr
    if (dabs(hami(i,j)-hami(j,i)).gt.0.001d0) then
       write(*,*)'Non-symmetric AD'
       write(*,'(2i5,3f10.5)')i,j,hami(i,j),hami(j,i),hami(i,j)-hami(j,i)
     endif    
   enddo
enddo

allocate(dmatr(ndimr,ndimr))
dmatr=0.0d0 

do i=1,ndimr
   do j=1,ndimr
     dmatr(i,j)=d1(i,mxtr(j))
   enddo
enddo
  
call reduce_mat(nx,hami,ndimr,no)
call reduce_mat(nx,dmatr,ndimr,no)

iout=0
if (iout.eq.1) then 
  write(998,*)    
  write(998,*)'******** matrix  D after Choleski *************'       
  write(998,*)
  do i=1,no
    write(998,102)(dmatr(i,j),j=1,no)
  enddo

  write(998,*)    
  write(998,*)'******** matrix  AD after Choleski *************'       
  write(998,*)
  do i=1,no
   write(998,102)(hami(i,j),j=1,no) 
  enddo
endif

!  Generalized EGV problem
   
write(*,*)' Diagonalisation '
lwork=20*no
allocate(work(lwork),ifail(no),wr(no),wro(no))
wr=0.d0
wro=0.d0
ifail=0
work=0.d0
vr=0.d0

!    call dsygvx(1,'V', 'V', 'U', no, hami, no, dmatr, no, vl, vu, il, iu, abstol, m, wr, vr, no, work, lwork, lwork, ifail, info)
!    write(*,*)' Number of eigevalues is interval (',vl,vu,')  :',m


call dsygv(1, 'V','U', no , hami, no, dmatr, no, wr, work, lwork, info)


write(99,*)'Tz = ',itzcal,' Parity = ',ipar, ' J = ',jcal
write(99,*)
!   write(99,*)(wro(i),i=1,no-ns) 
write(99,*)(wr(i),i=1,no) 
write(99,*)


!call dgemm('T','N',no,no,no,1.d0,vt,no,hami,no,0.d0,amatr,no)   !? spravne
!vr=amatr
!!!!!!!!!!!!!
! reduce rows of D
!deallocate(hamd)
allocate(hamd(no,ndim))

ii=0
do i=1,ndimr
if (nx(i).ne.0) then
ii=ii+1
!  do j=1,dim_base
  hamd(ii,:)=d1(i,:)
!  enddo
endif
enddo

deallocate(d1)

!write(99,*)' '
!write(99,*)' D matrix reduced rows'
!do i=1,dim_ind
! write(99,'(1000f15.10)')(hamd(i,j),j=1,dim_base)
!enddo

write(*,*)' check dim_ind =',ii
write(*,*)'X calculation'

allocate(xr(ndim,no))
!  X=DC

call dgemm('T','N',ndim,no,no,1.0d0,hamd,no,hami,no,0.d0,xr,ndim)

!write(998,*)' '
!write(998,*)' X matrix'
!do i=1,ndim
! write(998,'(1000f15.10)')(xr(i,j),j=1,no)
!enddo

xfact=1.d0*(dfloat(2*jcal+1))**0.5d0

xr=xr*xfact

!call normalize_c(no,ndim,ndimr,hami,xr,mxtr,nx,jcal)


return            
end subroutine ham_geev


!
      subroutine permutuj(ndim,d1,mxt)

      implicit double precision (a-h,o-z)

      double precision, dimension(:,:), allocatable :: d1
      double precision, dimension(:), allocatable :: work
      integer, dimension(:), allocatable :: mxt,nxt

      allocate (work(ndim))
      work=0.d0
      allocate(nxt(ndim))

      do i=1,ndim
       nxt(i)=i
      enddo

      do i=1,ndim
       work(:)=d1(i,:)
       do j=i,ndim
        if (nxt(j).eq.mxt(i)) then 
        d1(i,:)=d1(j,:)
        d1(j,:)=work(:)
        nxt(j)=nxt(i)
        nxt(i)=mxt(i)
       endif
       enddo
      enddo

      do i=1,ndim
       nxt(i)=i
      enddo

      do i=1,ndim
       work(:)=d1(:,i)
       do j=i,ndim
        if (nxt(j).eq.mxt(i)) then 
        d1(:,i)=d1(:,j)
        d1(:,j)=work(:)
        nxt(j)=nxt(i)
        nxt(i)=mxt(i)
        endif
       enddo
      enddo


 
      deallocate(work)

      return
      end subroutine permutuj


!***************************************************************************

      subroutine read_dmatch(dmatr)

      implicit double precision (a-h,o-z)


      double precision, dimension(:,:), allocatable :: dmatr
      double precision, dimension(:), allocatable :: dmm
      integer, dimension(:), allocatable :: nxtr

      
      dmatr=0.d0

      open(2,file='nxt.dat',status='unknown',form='unformatted')
      read(2)ndim,no
      allocate(nxtr(ndim))
      read(2)(nxtr(i),i=1,ndim)
      close(2)

      write(911,*)' idphon =', ndim
      write(911,*)(nxtr(i),i=1,ndim)


      open(66,file='d_mat.dat',status='old',form='unformatted')

      read(66)ndimrt,ndimt
      allocate(dmm(ndimt))
      dmm=0.0d0


!      do while (.not.eof(66))
!      read(66)i,j,dd
!      if (nxtr(i).ne.0) then 
!      dmatr(nxtr(i),j)=dd
!      endif

!      enddo

!      close(66)

      do iii=1,ndimrt
       read(66)(dmm(jjj),jjj=1,ndimt)
       do jjj=1,ndimt
         if (nxtr(iii).ne.0) then
          dmatr(nxtr(iii),jjj)=dmm(jjj)
         endif
       enddo
      enddo

      close(66)
      deallocate(dmm)


      deallocate(nxtr)

      end subroutine read_dmatch

!**********************************************************************
subroutine reduce_mat(nx,mat,ndim,ndim_red)
  implicit none
  double precision, dimension (:,:),allocatable :: mat,matc
  integer, dimension(:), allocatable :: nx
  integer i,j,ii,jj,ndim,ndim_red
  
  allocate(matc(ndim,ndim))
  matc=mat
  deallocate(mat)
  allocate(mat(ndim_red,ndim_red))
  
  
  ii=0
  do i=1,ndim
   if (nx(i).ne.0) then
     ii=ii+1
     jj=0
     do j=1,ndim
      if (nx(j).ne.0) then
         jj=jj+1
         mat(ii,jj)=matc(i,j)
       endif
      enddo
  
    endif
  enddo
  deallocate(matc)
  
  end subroutine reduce_mat
  


end module hami
