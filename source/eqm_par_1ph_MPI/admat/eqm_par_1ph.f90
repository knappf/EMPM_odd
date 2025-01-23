! last modification 25.2.2023

program eqm_par_1phon

use input_sp
use base_h_phon
use admatr
use dens_list
use types_eqm

implicit none 

include 'mpif.h'
!include 'types_eqm_hole.inc'

type(level_typ),dimension(:), allocatable :: lev
type(phon_typ), dimension (:), allocatable :: phon
type(phonbase_typ), dimension (:), allocatable :: phonbs
integer :: iout,i,j,nlev,jmax,ih_lev,ip_lev,ii,ipozz
integer :: nf,dim_phon,dim_sp,ia,iz,ihnmn,ihnmx,ihpmn,ihpmx,ipnmn,ipnmx,ippmn,ippmx,isp_min,isp_max,n_spur
integer :: ipar,ijj,dim_base,dim_baser,dim_ind,it_bs,nlam
integer, dimension(:,:), allocatable :: i_ph_int
integer, dimension(:), allocatable :: par_lev, hol_lev
double precision, dimension(:,:), allocatable :: amat,dmat,ham,hamd,dinv,dmatc,vr,vl,xamp
integer, dimension (:), allocatable :: ind_red,ind_red_ind
integer, dimension(:), allocatable :: nx,ipoz
integer :: info,lwork
double precision :: xe 
double precision, dimension(:), allocatable ::  work,wr,wi,wro
character*30 namex,names,namec
integer :: myid,ierr,numprocs

call MPI_INIT( ierr )
call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

nlam=0

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
 dim_phon=500000 ! dimension of phonon space
 
 isp_min=1
 isp_max=ip_lev

 

open(23,file='AD_J_Pi_Tz.dat',status='old',form='formatted')
read(23,*)ipar,ijj,it_bs
close(23) 

write(*,*)' Parity?'
write(*,*)ipar
write(*,*)' 2*J?'
write(*,*)ijj
write(*,*) ' Tz =?   1 for neutron particle ,  -1 for proton particle '
write(*,*)it_bs

call odd_phonbase(nf,ipar,ijj,it_bs,dim_phon,lev,phon,isp_min,isp_max,par_lev,dim_base,dim_baser,phonbs,ind_red,n_spur)

write(*,*)'Dimension of space =',dim_base
write(*,*)'Dimension of Tz= ',it_bs,' subspace  =',dim_baser

call dmatrix(dim_phon,ijj,dim_base,dim_baser,phonbs,lev,phon,ih_lev,ip_lev,hol_lev,par_lev,i_ph_int,nlev,ind_red,myid,numprocs)
call amatrix(dim_phon,ijj,dim_base,dim_baser,phonbs,lev,phon,ih_lev,ip_lev,hol_lev,par_lev,i_ph_int,nlev,ind_red,myid,numprocs)

call MPI_FINALIZE(ierr)

end 
