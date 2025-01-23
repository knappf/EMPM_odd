! last modification 24.2.2023

program eqm_hole_phon

use input_sp
use base_p_phon
use dmatr

implicit none 

include 'types_eqm_hole.inc'

type(level_typ),dimension(:), allocatable :: lev
type(phon_typ), dimension (:), allocatable :: phon
type(phonbase_typ), dimension (:), allocatable :: phonbs
integer :: iout,i,j,nlev,jmax,ih_lev,ip_lev,ii
integer :: nf,dim_phon,dim_sp,ia,iz,ihnmn,ihnmx,ihpmn,ihpmx,ipnmn,ipnmx,ippmn,ippmx,isp_min,isp_max
integer :: ipar,ijj,dim_base,dim_baser,it_bs
integer, dimension(:,:), allocatable :: i_ph_int
integer, dimension(:), allocatable :: par_lev, hol_lev
double precision, dimension(:,:), allocatable :: amat,dmat,ham
integer, dimension (:), allocatable :: ind_red

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

 nf=1
 dim_phon=1000 ! dimension of phonon space
 

 isp_min=1
 isp_max=ip_lev


 write(*,*)' Parity?'
 read(*,*)ipar
 write(*,*)ipar
 write(*,*)' 2*J?'
 read(*,*)ijj
 write(*,*)ijj
 write(*,*) ' Tz =?   1 for neutron particle,  -1 for proton particle'
 read(*,*)it_bs
 write(*,*)it_bs

 !it_bs=-1

 call odd_phonbase(nf,ipar,ijj,it_bs,dim_phon,lev,phon,isp_min,isp_max,par_lev,dim_base,dim_baser,phonbs,ind_red)

 write(*,*)'Dimension of space =',dim_base
 write(*,*)'Dimension of Tz= ',it_bs,'   =',dim_baser


 call dmatrix(ijj,dim_base,dim_baser,phonbs,lev,phon,ih_lev,ip_lev,hol_lev,par_lev,i_ph_int,nlev,ind_red,dmat)


 call amatrix(ijj,dim_base,dim_baser,phonbs,lev,phon,ih_lev,ip_lev,hol_lev,par_lev,i_ph_int,nlev,ind_red,amat)

allocate(ham(dim_baser,dim_baser))
ham=0.d0

call dgemm('N','T',dim_baser,dim_baser,dim_base,1.d0,amat,dim_baser,dmat,dim_baser,0.d0,ham,dim_baser)

iout=1
     if (iout.eq.1) then
      write(99,*)
      write(99,*)'******** matrix  AD *************'
      write(99,*)
      do i=1,dim_baser
        write(99,'(1000f15.10)')(ham(i,j),j=1,dim_baser)
      enddo
     endif

do i=1,dim_baser
 do j=i,dim_baser
   if (dabs(ham(i,j)-ham(j,i)).gt.0.00001d0) write(*,'(2i5,3f10.5)'),i,j,dabs(ham(i,j)-ham(j,i)),dmat(i,i),dmat(j,j)
 enddo
enddo



end 
