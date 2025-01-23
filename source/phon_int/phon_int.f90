!*     
!*     Program phon_dens1 computes 1phonon densities in proton-neutron 
!*     J-coupled formalism.

!*     last update 16.4.2013

program phon_int

use input_sp
use phoninterac

implicit double precision (a-h,o-z)

include 'input_phon_int.inc'
include 'formats_phon_int.inc'
include 'types_phon_int.inc'

type(level_typ),dimension(:), allocatable :: lev
integer, dimension(:,:), allocatable :: i_ph_int
integer, dimension(:), allocatable :: par_lev, hol_lev
integer :: dim_sp
character*3 :: type_phon
      
!     loading of input data 

xrotrunc=1.d-10
     
write(*,*)'Loading of input '

allocate(i_ph_int(4,-1:1))

  open(1,file='input_tda_coup.dat',status='old',form='formatted')
    read(1,'(30i6)')ia,iz
    read(1,'(30i6)')i_ph_int(1,1),i_ph_int(2,1)
    read(1,'(30i6)')i_ph_int(1,-1),i_ph_int(2,-1)
    read(1,'(30i6)')i_ph_int(3,1),i_ph_int(4,1)
    read(1,'(30i6)')i_ph_int(3,-1),i_ph_int(4,-1)
  close(1)


 dim_sp=1000 ! dimension of single particle space

 allocate(lev(dim_sp))

 call inp_sp(lev,nlev,jmax)

     allocate(par_lev(nlev),hol_lev(nlev))

      ii=0
      do i=i_ph_int(3,-1),i_ph_int(4,-1)
       ii=ii+1
       par_lev(ii)=2*i-1
      enddo

      do i=i_ph_int(3,1),i_ph_int(4,1)
       ii=ii+1
       par_lev(ii)=2*i
      enddo
      ip_lev=ii

      ii=0
      do i=i_ph_int(1,-1),i_ph_int(2,-1)
       ii=ii+1
       hol_lev(ii)=2*i-1
      enddo

      do i=i_ph_int(1,1),i_ph_int(2,1)
       ii=ii+1
       hol_lev(ii)=2*i
      enddo
      ih_lev=ii


!      jmaxn=0
!      do i=1,ippmx
!       if (lev(i)%j.gt.jmaxn) jmaxn=lev(i)%j
!      enddo


write(*,*)'jmax =',jmax

!      jmax=jmaxn
      
!c      call checkf
!c      call checkfmat
!c      call checkro
      
!c      stop      
      
ifmx=1000

!write(*,*)'energy threshold for 1 phonon states'
!read(*,*)xthrun_min,xthrun_max
!write(*,*)xthrun_min,xthrun_max
xthtun_min=-100000000.d0 
xthrun_max=100000000.d0

write(*,*)'Type of phonon-phonon intetaction: par/hol'
read(*,*)type_phon
write(*,*)type_phon

if (type_phon.eq.'par') call vintp(ifmx,jmax,lev,nlev,xthrun_min,xthrun_max,ih_lev,ip_lev,par_lev,hol_lev)  ! particle 
if (type_phon.eq.'hol') call vinth(ifmx,jmax,lev,nlev,xthrun_min,xthrun_max,ih_lev,ip_lev,par_lev,hol_lev) !  hole


      
end 
!*     
!*     END of the main program 

!* 

