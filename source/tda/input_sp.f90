
 module input_sp


  contains 

  subroutine inp_sp(lev,nlev,jmax)

      implicit none

      include 'types_tda_cp.inc'

      integer :: i,ii,jmax,nlev
      integer, dimension(:), allocatable:: il_n,il_p   
      type(level_typ),dimension(:), allocatable :: lev     

      allocate(il_n(0:1000))
      il_n=0

      allocate(il_p(0:1000))
      il_p=0

      jmax=0
   
      open(1,file='singpart_coup.dat',status='old',form='formatted')

      read(1,*)
!      do i=1,ipnmx
!       read(1,16)nt,lt,jt,erspnt,ersppt
       i=0
       do while (.not.eof(1))
        i=i+1
        read(1,'(4i5,7x,f15.5)')lev(i)%tz,lev(i)%n,lev(i)%l,lev(i)%j,lev(i)%en

        if (lev(i)%j.gt.jmax) jmax=lev(i)%j



      enddo
      write(*,*)' Number of levels loaded ',i
      write(*,*)' Maximal value 2J',jmax
      nlev=i
 
      close(1)

   end subroutine inp_sp


   subroutine read_fmat(f_int,nlev,jmax)
   implicit none 


   double precision, dimension(:,:,:,:,:), allocatable :: f_int
   integer(kind=2) :: i_i,i_j,i_k,i_l
   integer(kind=1) :: j_f
   integer :: i,j,k,l,ijt,nlev,jmax
   

   allocate (f_int(nlev,nlev,nlev,nlev,0:jmax))
   f_int=0.d0 
  
   open(1,file='fmat.dat',status='old',form='unformatted')

     do while (.not.eof(1))
      read(1)j_f,i_i,i_j,i_k,i_l,f_int(i_i,i_j,i_k,i_l,j_f)  
     enddo 

   close(1)

   end subroutine read_fmat

!**************************************************************************
!*     phbase generates particle-hole combinations with parity
!*     ipar and J and stores indices ip and ih to arrays
!*     iphp, iphn for protons, neutrons respectively.
!**************************************************************************
      subroutine phbase(ipar,ijj,i_ph_int,lev,ph_base,idim_ph)

      implicit double precision (a-h,o-z)

      include 'types_tda_cp.inc'
      include 'input_tda_cp.inc'

      type(level_typ),dimension(:),allocatable  :: lev
      type(ph_typ),dimension(:),allocatable :: ph_base
      integer, dimension(:,:), allocatable :: i_ph_int

      i=0

!   proton - proton
      do ip=i_ph_int(3,-1),i_ph_int(4,-1)
        iparp=(-1)**lev(2*ip-1)%l

        do ih=i_ph_int(1,-1),i_ph_int(2,-1)
         iparh=(-1)**lev(2*ih-1)%l
         ipart=iparh*iparp
          jjmn=iabs(lev(2*ip-1)%j-lev(2*ih-1)%j)/2
          jjmx=(lev(2*ip-1)%j+lev(2*ih-1)%j)/2

         if (ipart.eq.ipar.and.ijj.ge.jjmn.and.ijj.le.jjmx) then
          i=i+1
          ph_base(i)%par=2*ip-1
          ph_base(i)%hol=2*ih-1 
          ph_base(i)%tz=0
         endif

       enddo
      enddo
! neutron neutron

      do ip=i_ph_int(3,1),i_ph_int(4,1)
        iparp=(-1)**lev(2*ip)%l

        do ih=i_ph_int(1,1),i_ph_int(2,1)
         iparh=(-1)**lev(2*ih)%l
         ipart=iparh*iparp
          jjmn=iabs(lev(2*ip)%j-lev(2*ih)%j)/2
          jjmx=(lev(2*ip)%j+lev(2*ih)%j)/2

         if (ipart.eq.ipar.and.ijj.ge.jjmn.and.ijj.le.jjmx) then
          i=i+1
          ph_base(i)%par=2*ip
          ph_base(i)%hol=2*ih
          ph_base(i)%tz=0
         endif

       enddo
      enddo
!     proton- neutron
      do ip=i_ph_int(3,-1),i_ph_int(4,-1)
        iparp=(-1)**lev(2*ip-1)%l

        do ih=i_ph_int(1,1),i_ph_int(2,1)
         iparh=(-1)**lev(2*ih)%l
         ipart=iparh*iparp
          jjmn=iabs(lev(2*ip-1)%j-lev(2*ih)%j)/2
          jjmx=(lev(2*ip-1)%j+lev(2*ih)%j)/2

         if (ipart.eq.ipar.and.ijj.ge.jjmn.and.ijj.le.jjmx) then
          i=i+1
          ph_base(i)%par=2*ip-1
          ph_base(i)%hol=2*ih
          ph_base(i)%tz=-1
         endif

       enddo
      enddo

    do ip=i_ph_int(3,1),i_ph_int(4,1)
        iparp=(-1)**lev(2*ip)%l

        do ih=i_ph_int(1,-1),i_ph_int(2,-1)
         iparh=(-1)**lev(2*ih-1)%l
         ipart=iparh*iparp
          jjmn=iabs(lev(2*ip)%j-lev(2*ih-1)%j)/2
          jjmx=(lev(2*ip)%j+lev(2*ih-1)%j)/2

         if (ipart.eq.ipar.and.ijj.ge.jjmn.and.ijj.le.jjmx) then
          i=i+1
          ph_base(i)%par=2*ip
          ph_base(i)%hol=2*ih-1
          ph_base(i)%tz=1
         endif

       enddo
      enddo


      idim_ph=i

      return
    end subroutine phbase

!*********************************************************************


end module input_sp



