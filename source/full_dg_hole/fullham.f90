program fullham

use full_mat_hole

implicit none

include 'types_fullham.inc'


integer :: ia,iz, nlev,ii,i,jmax,j
integer :: iparc,jcalc,itzc,nfon,iff,ih_lev,ip_lev
integer, dimension(:), allocatable :: ndmpho
integer, dimension(:,:), allocatable :: nbs
type (phon_typ), dimension (:,:), allocatable :: phonbs
type(level_typ),dimension(:), allocatable :: lev
integer, dimension(:,:), allocatable :: i_ph_int
integer, dimension(:), allocatable :: par_lev, hol_lev
 
type (amp_typ), dimension(:,:), allocatable :: xamp
integer, dimension (:), allocatable :: ndx
double precision, dimension(:,:), allocatable :: vlam
double precision, dimension(:), allocatable :: diag_h
character(len=30)fnamex
double precision :: fact


!write(*,*)'Loading of input '

!open(1,file='input_tda_coup.dat',status='old',form='formatted')

!read(1,'(30i6)')ia,iz
!read(1,'(30i6)')i_ph_int(1,1),i_ph_int(2,1)
!read(1,'(30i6)')i_ph_int(1,-1),i_ph_int(2,-1)
!read(1,'(30i6)')i_ph_int(3,1),i_ph_int(4,1)
!read(1,'(30i6)')i_ph_int(3,-1),i_ph_int(4,-1)

!close(1)



! nlev=300
! allocate(lev(nlev))

!      call inp_sp(lev,nlev,jmax)

!      allocate(par_lev(nlev),hol_lev(nlev))

!      ii=0
!      do i=i_ph_int(3,-1),i_ph_int(4,-1)
!       ii=ii+1
!       par_lev(ii)=2*i-1
!      enddo

!      do i=i_ph_int(3,1),i_ph_int(4,1)
!       ii=ii+1
!       par_lev(ii)=2*i
!      enddo
!      ip_lev=ii

!      ii=0
!      do i=i_ph_int(1,-1),i_ph_int(2,-1)
!       ii=ii+1
!       hol_lev(ii)=2*i-1
!      enddo

!      do i=i_ph_int(1,1),i_ph_int(2,1)
!       ii=ii+1
!       hol_lev(ii)=2*i
!      enddo
!      ih_lev=ii



write(*,*)' Parity ? '
read(*,*)iparc
write(*,*)' J ? '
read (*,*)jcalc
write(*,*)' Tz ?'
read(*,*)itzc
write(*,*)' Nf of full diagonalization?'
read(*,*)nfon



write(*,*)' Construction of total H matrix'
write(*,*)' Parity = ',iparc
write(*,*)' Angular momentum = ',jcalc

allocate(diag_h(50000))
ii=0
do iff=0,nfon
      call bstot(iff,iparc,jcalc,itzc,nbs,ndmpho,phonbs,par_lev,hol_lev)
! diagonal part of H matrix
      fact=1.0d0
      if (iff.eq.0) fact=-1.d0  ! for holes 
      do i=1,ndmpho(iff)
       ii=ii+1
       j=nbs(iff,i)
       diag_h(ii)=fact*phonbs(iff,j)%enf
      enddo
 enddo

      write(*,*)'Dimension of 0 phonon space =  ',ndmpho(0)
      write(*,*)'Dimension of 1 phonon space =  ',ndmpho(1)


      if (ndmpho(1).gt.0.and.ndmpho(0).gt.0) then
       write(*,*)'Calculating nondiagonal part 01'

        fnamex='1phon_hole/1f_x.dat'
        call readx(fnamex,xamp,ndx,iparc,jcalc,itzc)
        call inter_hol_phon1(ndmpho,nbs,xamp,ndx,vlam,jcalc,par_lev,hol_lev)
!       call nondiag2(2,nbs,ndmpho,phonbs,iparc,jcalc)
      endif


      call fulhamt(ndmpho,nbs,diag_h,vlam)



end





