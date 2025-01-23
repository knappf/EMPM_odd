! last modification 25.2.2023

program eqm_par_2phon

use input_sp
use base_h_phon
use dens_list
use choleski
use hami
use types_eqm

implicit none 

!include 'types_eqm_hole.inc'

type(level_typ),dimension(:), allocatable :: lev
type(phon_typ), dimension (:), allocatable :: phon
type(phonbase_typ), dimension (:), allocatable :: phonbs
integer :: iout,i,j,nlev,jmax,ih_lev,ip_lev,ii,ipozz
integer :: nf,dim_phon,dim_sp,ia,iz,ihnmn,ihnmx,ihpmn,ihpmx,ipnmn,ipnmx,ippmn,ippmx,isp_min,isp_max
integer :: ipar,ijj,dim_base,dim_baser,dim_ind,it_bs,nlam
integer, dimension(:,:), allocatable :: i_ph_int
integer, dimension(:), allocatable :: par_lev, hol_lev
double precision, dimension(:,:), allocatable :: amat,dmat,ham,hamd,dinv,dmatc,vr,vl,xamp,xr,h_corr
integer, dimension (:), allocatable :: ind_red,ind_red_ind
integer, dimension(:), allocatable :: nx,mxt,ipoz,mxtr,irow
integer :: ns, no, nor, n_spur
integer :: info,lwork
double precision :: xe 
double precision, dimension(:), allocatable ::  work,wr,wi,wro
character*30 namex,names,namec



nlam=0
namex='2phon_hole/2f_x.dat'
namec='2phon_hole/2f_c.dat'
names='2phon_hole/2f_states.dat'


open(9,file='egv_h2phon.dat',status='unknown',form='formatted')
open(12,file=namex,status='unknown',form='unformatted')
open(22,file=namec,status='unknown',form='unformatted')
open(13,file=names,status='unknown',form='unformatted')




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

 nf=2
 dim_phon=500000 ! dimension of phonon space
 

 isp_min=1
 isp_max=ih_lev


 write(*,*)' Parity?'
 read(*,*)ipar
 write(*,*)' 2*J?'
 read(*,*)ijj
 write(*,*) ' Tz =?   -1 for neutron hole,  1 for proton hole'
 read(*,*)it_bs
 write(*,*)it_bs

open(23,file='AD_J_Pi_Tz.dat',status='unknown',form='formatted')
write(23,*)ipar,ijj,it_bs
close(23) 

 call odd_phonbase(nf,ipar,ijj,it_bs,dim_phon,lev,phon,isp_min,isp_max,hol_lev,dim_base,dim_baser,phonbs,ind_red)

 !open(62,file='2_phon_dens_calc.dat',status='unknown',form='formatted')
 !close(62)
 call create_dens_list

 call execute_command_line('cp 2_phon_dens_calc_new.dat 2_phon_dens_list.dat')
 call execute_command_line('./run_dens2.sh' )

 call execute_command_line('cat 2_phon_dens_calc_myid* >> 2_phon_dens_calc.dat')
 call execute_command_line('rm 2_phon_dens_calc_myid*')

 write(*,*)'2-phon dens. calculated'

 write(*,*)'Dimension of space =',dim_base
 write(*,*)'Dimension of Tz= ',it_bs,' subspace  =',dim_baser

CALL execute_command_line('./run_admat_h2ph.sh' )

n_spur=0 ! no CM ort. for the moment 
open(3,file='chol_spur_subspace.dat',form='formatted',status='unknown')
write(3,*)n_spur
do i=1,n_spur
  write(3,*)nx(i)
enddo
close(3)

call cholesk(0,dim_baser,dim_ind,dim_baser,nx,mxt)

!stop

call ham_geev(dim_base,dim_baser,dim_ind,ns,irow,wr,xr,vr,ipar,ijj,it_bs,ind_red,phonbs,nx,h_corr)


!do i=1,dim_ind
!  if (dabs(wi(i)).gt.1.d-10) write(*,*)' Imaginary part',i,wi(i)
!enddo
  
      
!write(99,*)' '
!write(99,*)' X matrix'
!do i=1,dim_base
! write(99,'(1000f15.10)')(xr(i,j),j=1,dim_ind)
!enddo


!write(99,*)' '
!write(99,*)' C matrix normalized'
!do i=1,dim_ind
! write(99,'(1000f15.10)')(vr(i,j),j=1,dim_ind)
!enddo

allocate(ind_red_ind(dim_ind))

ii=0
do i=1,dim_baser
if (nx(i).eq.1) then
ii=ii+1
ind_red_ind(ii)=ind_red(i)
endif
enddo

open(1,file=' phonon_base_ind.dat',status='unknown',form='formatted')
   write(1,*)'i   2*Tz    i_sp      en_sp      lam      en_phon'
    do i=1,dim_ind
     write(1,'(3i5,f10.5,i5,f10.5)')i,phonbs(ind_red_ind(i))%tz,phonbs(ind_red_ind(i))%isp,lev(phonbs(ind_red_ind(i))%isp)%en,phonbs(ind_red_ind(i))%ila,phon(phonbs(ind_red_ind(i))%ila)%enf
    enddo

close(1)

do i=1,dim_ind
 write(999,'(3i5,f10.5,i5,3f10.5)')i,phonbs(ind_red_ind(i))%tz,phonbs(ind_red_ind(i))%isp,lev(phonbs(ind_red_ind(i))%isp)%en,phonbs(ind_red_ind(i))%ila,phon(phonbs(ind_red_ind(i))%ila)%enf,vr(i,1),vr(i,2)

enddo



do i=1,dim_ind
  write(13)nlam+i,ipar,ijj,it_bs,wr(i)
enddo


write(12)ipar,ijj,it_bs,dim_ind,dim_base
write(22)ipar,ijj,it_bs,dim_ind,dim_ind

do j=1,dim_ind
write(12)(phonbs(i)%isp,phonbs(i)%ila,xr(i,j),i=1,dim_base)


write(22)(phonbs(ind_red_ind(i))%isp,phonbs(ind_red_ind(i))%ila,vr(i,j),i=1,dim_ind)

enddo

nlam=nlam+dim_ind


close(9)
close(12)
close(22)
close(13)

end 
