
 module input_sp


  contains 

  subroutine inp_sp(lev,nlev,jmax)

      implicit none

      include 'types_eqm.inc'

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


!*********************************************************************


end module input_sp



