! last modification 25.4.2019

program eqm_hole_phon

use input_sp
use base_h_phon
use dmatr

implicit none 

include 'types_eqm_hole.inc'

type(level_typ),dimension(:), allocatable :: lev
type(phon_typ), dimension (:), allocatable :: phon
type(phonbase_typ), dimension (:), allocatable :: phonbs
integer :: iout,i,j,nlev,jmax,ih_lev,ip_lev,ii,ipozz
integer :: nf,dim_phon,dim_sp,ia,iz,ihnmn,ihnmx,ihpmn,ihpmx,ipnmn,ipnmx,ippmn,ippmx,isp_min,isp_max
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



nlam=0
namex='1phon_hole/1f_x.dat'
namec='1phon_hole/1f_c.dat'
names='1phon_hole/1f_states.dat'

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

 nf=1
 dim_phon=1000 ! dimension of phonon space
 

 isp_min=1
 isp_max=ih_lev


 write(*,*)' Parity?'
 read(*,*)ipar
 write(*,*)' 2*J?'
 read(*,*)ijj
 write(*,*) ' Tz =?   -1 for neutron hole,  1 for proton hole'
 read(*,*)it_bs
 write(*,*)it_bs

 call odd_phonbase(nf,ipar,ijj,it_bs,dim_phon,lev,phon,isp_min,isp_max,hol_lev,dim_base,dim_baser,phonbs,ind_red)

 write(*,*)'Dimension of space =',dim_base
 write(*,*)'Dimension of Tz= ',it_bs,'   =',dim_baser


 call dmatrix(ijj,dim_base,dim_baser,phonbs,lev,phon,ih_lev,ip_lev,hol_lev,par_lev,i_ph_int,nlev,ind_red,dmat)
 call dmat_ind_set(dmat,dmatc,dim_base,dim_baser,dim_ind,nx,ind_red)


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
   if (dabs(ham(i,j)-ham(j,i)).gt.0.01d0) write(*,'(2i5,3f10.5)'),i,j,dabs(ham(i,j)-ham(j,i)),dmat(i,i),dmat(j,j)
 enddo
enddo

call reduce_mat(nx,ham,dim_baser,dim_ind)
call reduce_mat(nx,dmatc,dim_baser,dim_ind)

iout=1
     if (iout.eq.1) then
      write(99,*)
      write(99,*)'******** reduced D *************'
      write(99,*)
      do i=1,dim_ind
        write(99,'(1000f15.10)')(dmatc(i,j),j=1,dim_ind)
      enddo

      write(99,*)
      write(99,*)'******** reduced AD *************'
      write(99,*)
      do i=1,dim_ind
        write(99,'(1000f15.10)')(ham(i,j),j=1,dim_ind)
      enddo
     endif
! inverse 
      call dpotrf('U',dim_ind,dmatc,dim_ind,info)
      write(*,*)' Factorization info ',info

      call dpotri('U',dim_ind,dmatc,dim_ind,info)

      write(*,*)' Inverse info ',info
      do i=1,dim_ind
        do j=i+1,dim_ind
          dmatc(j,i)=dmatc(i,j)
        enddo
       enddo

! D^-1 AD 
!      deallocate(dmat)
      allocate(hamd(dim_ind,dim_ind))
      call dgemm('N','N',dim_ind,dim_ind,dim_ind,1.d0,dmatc,dim_ind,ham,dim_ind,0.d0,hamd,dim_ind)

!  diag
      write(*,*)' Diagonalisation '
      lwork=20*dim_ind 
      allocate(work(lwork),wi(dim_ind),wr(dim_ind),vr(dim_ind,dim_ind),wro(dim_ind))
      wr=0.d0
      wro=0.d0
      wi=0.d0
      work=0.d0
      vr=0.d0

      if (dim_ind.gt.0) then

      CALL DGEEV('N','V',dim_ind,hamd,dim_ind,wr,wi,vl,dim_ind,vr,dim_ind,work,lwork,info)
      write(*,*)' info=  ',info

write(99,*)' '
write(99,*)' C matrix non-normalized'
do i=1,dim_ind
 write(99,'(1000f15.10)')(vr(i,j),j=1,dim_ind)
enddo




      do i=1,dim_ind
        if (dabs(wi(i)).gt.1.d-10) write(*,*)' Imaginary part',i,wi(i)
      enddo

      allocate(ipoz(dim_ind))

      do i=1,dim_ind
        ipoz(i)=i
        wro(i)=wr(i)
      enddo

      do i=1,dim_ind
        do j=1,dim_ind-i
          if (wro(j).gt.wro(j+1)) then
            xe=wro(j)
            wro(j)=wro(j+1)
            wro(j+1)=xe
            ipozz=ipoz(j)
            ipoz(j)=ipoz(j+1)
            ipoz(j+1)=ipozz
          endif
        enddo
      enddo

       open(9,file='egv.dat',status='unknown',form='formatted')
        write(9,*)' Parity = ',ipar, ' J = ',ijj, ' Tz =',it_bs
        write(9,*)
        write(9,*)(wro(i),i=1,dim_ind)
        write(9,*)
       close(9)
      endif

! reduce rows of D
deallocate(hamd)
allocate(hamd(dim_ind,dim_base))

ii=0
do i=1,dim_baser
 if (nx(i).ne.0) then 
   ii=ii+1
!  do j=1,dim_base
      hamd(ii,:)=dmat(i,:)
!  enddo
 endif 
enddo

!write(99,*)' '
!write(99,*)' D matrix reduced rows'
!do i=1,dim_ind
! write(99,'(1000f15.10)')(hamd(i,j),j=1,dim_base)
!enddo


write(*,*)' check dim_ind =',ii
write(*,*)'X calculation'

allocate(xamp(dim_base,dim_ind))
!  X=DC

call dgemm('T','N',dim_base,dim_ind,dim_ind,1.0d0,hamd,dim_ind,vr,dim_ind,0.d0,xamp,dim_base) 

write(99,*)' '
write(99,*)' X matrix'
do i=1,dim_base
 write(99,'(1000f15.10)')(xamp(i,j),j=1,dim_ind)
enddo


!xamp=-1.0d0*(ijj+1)**0.5d0*xamp

call normal(dim_ind,dim_base,dim_baser,vr,xamp,ind_red,nx,ijj)

write(99,*)' '
write(99,*)' C matrix normalized'
do i=1,dim_ind
 write(99,'(1000f15.10)')(vr(i,j),j=1,dim_ind)
enddo




write(99,*)' '
write(99,*)' X matrix normalized'
do i=1,dim_base
 write(99,'(1000f15.10)')(xamp(i,j),j=1,dim_ind)
enddo


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
write(12)(phonbs(i)%isp,phonbs(i)%ila,xamp(i,j),i=1,dim_base)


write(22)(phonbs(ind_red_ind(i))%isp,phonbs(ind_red_ind(i))%ila,vr(i,j),i=1,dim_ind)

enddo

nlam=nlam+dim_ind


close(12)
close(22)
close(13)

end 
