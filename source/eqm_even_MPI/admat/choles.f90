!     Cholesky    procedure      
!     LAPACK implementation 
!     last modification 25.3.2020

module choleski

! use read_admat      

 contains
       
 subroutine cholesk(i_spur,ndim,no,noo,dd,nx,mxt)

      implicit double precision (a-h,o-z)

      include 'formats_eqm.inc'
     
      
      double precision, dimension(:,:), allocatable :: dd
    
      double precision, dimension(:), allocatable :: r,work 

      integer, dimension(:), allocatable :: nx,mxt
 
      if (i_spur.eq.0) then
       allocate(dd(ndim,ndim))
       dd=0.d0
       call read_dmat(dd)
      endif  
 
      allocate(work(10*ndim))
      allocate(mxt(ndim),nx(ndim))

      mxt=0
      nx=0 

      
      tol=0.0001d0

      write(*,*)' dimension = ',ndim
            
      call dpstrf('U', ndim, dd, ndim, mxt, no, tol, work, info )

      write(*,*)'Cholesky info =', info


      do i=1,no
        nx(mxt(i))=1
      enddo

      open(6,file='mxt.dat',status='unknown',form='unformatted')

      write(6)ndim,no
      do i=1,ndim
        write(6)i,mxt(i)
      enddo
      close(6)

      
      write(*,*)' Number of linearly independent states ',no


 
      return 
 end subroutine cholesk

     subroutine read_dmat(dmatr)

      implicit double precision (a-h,o-z)

!c      include 'formats_eqm.inc'

      double precision, dimension(:,:), allocatable :: dmatr
       double precision, dimension(:), allocatable :: dmm

      integer, dimension(:), allocatable :: nxtr
!c      double precision, dimension(:), allocatable :: work,e

      
      dmatr=0.d0

      open(2,file='nxtr.dat',status='unknown',form='unformatted')
      read(2)idphon
      allocate(nxtr(idphon))
      read(2)(nxtr(i),i=1,idphon)
      close(2)

      write(991,*)' idphon =', idphon
      write(991,*)(nxtr(i),i=1,idphon)

!c      iold=-100


      open(66,file='d_mat.dat',status='old',form='unformatted')

      read(66)ndimrt,ndimt
      allocate(dmm(ndimt))
      dmm=0.0d0

!      if (ndimrt.ne.ndimr.or.ndimt.ne.ndim) then
!        write(*,*)' Dimensions does not match in d_m file'
!        stop
!      endif
      do iii=1,ndimrt
       read(66)(dmm(jjj),jjj=1,ndimt)
       do jjj=1,ndimt
         if (nxtr(jjj).ne.0) then
          dmatr(iii,nxtr(jjj))=dmm(jjj)
         endif
       enddo
      enddo

      close(66)
      deallocate(dmm)




!      do while (.not.eof(66))
!      read(66)i,j,dd
!      if (nxtr(j).ne.0) then 
!      dmatr(i,nxtr(j))=dd
!      endif

!      enddo

!      close(66)

      deallocate(nxtr)
     

      end subroutine read_dmat

 
end module choleski
     
