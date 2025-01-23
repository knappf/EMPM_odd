module dmatr

use anglib 
use read_densfmat

contains  

subroutine dmatrix(j_tot,dim_base,dim_baser,phonbs,lev,phon,ih_lev,ip_lev,hol_lev,par_lev,i_ph_int,nlev,ind_red,d_mat)


implicit none 

include 'types_eqm_hole.inc'

type(level_typ),dimension(:), allocatable :: lev
type(phon_typ), dimension (:), allocatable :: phon
type(phonbase_typ), dimension (:), allocatable :: phonbs
integer :: i,j,dim_base,dim_baser,dim_int,nronh,ii,isi,j_tot,iout,ifaz,nlev,itz
integer :: i_sp,j_sp,i_lam,j_lam,i_lam_old,ip_lev,ih_lev
integer :: ihp_min,ihp_max,ihn_min,ihn_max,ipp_min,ipp_max,ipn_min,ipn_max,isi_max,ndla_max
integer :: lwork,info
integer, dimension(:,:), allocatable :: i_ph_int
integer, dimension(:), allocatable :: par_lev, hol_lev
integer, dimension (:), allocatable :: ind_red

double precision, dimension (:,:),allocatable :: d_mat,d_matc
double precision, dimension(:,:,:,:),allocatable :: ronh1
double precision :: rac
double precision, dimension(:), allocatable :: work,e


allocate(d_mat(dim_baser,dim_base))
d_mat=0.0d0
i_lam_old=-10
isi_max=15

call read_input_par(ihp_min,ihp_max,ihn_min,ihn_max,ipp_min,ipp_max,ipn_min,ipn_max,ndla_max)

itz=-1

  do ii=1,dim_baser
     i=ind_red(ii)
     i_sp=phonbs(i)%isp
     i_lam=phonbs(i)%ila 

    
     if (i_lam.ne.i_lam_old) then 
          call readro11('1f_rh.dat',i_lam,ronh1,ndla_max,1,nlev,1,nlev,isi_max)
          i_lam_old=i_lam
     endif 

   do j=1,dim_base 
     j_sp=phonbs(j)%isp
     j_lam=phonbs(j)%ila

!     if (levn(i_sp)%j.eq.levn(j_sp)%j.and.phon(i_lam)%j.eq.phon(j_lam)%j) d_mat(i,j)=1.0d0
     
     if (i_sp.eq.j_sp.and.i_lam.eq.j_lam) d_mat(ii,j)=1.0d0

 !    ifaz=(i_lam+j_lam)+(levn(i_sp)%j-levn(j_sp)%j)/2
     ifaz=(phon(j_lam)%j+phon(i_lam)%j)+(lev(i_sp)%j-lev(j_sp)%j)/2
     ifaz=(-1)**(ifaz)

     do isi=0,isi_max
       d_mat(ii,j)=d_mat(ii,j)+dfloat(ifaz)*ronh1(j_lam,isi,j_sp,i_sp)*(2*isi+1)**0.5d0*racah(2*isi,lev(j_sp)%j,2*phon(i_lam)%j,j_tot,lev(i_sp)%j,2*phon(j_lam)%j)
     enddo


   enddo
  enddo
 

!  allocate(d_matc(dim_baser,dim_baser))
!  do i=1,dim_baser
!   do j=1,dim_baser
!      d_matc(i,j)=d_mat(i,ind_red(j))
!   enddo
!  enddo


!    iout=0
!     if (iout.eq.1) then
!      write(99,*)
!      write(99,*)'******** matrix  D *************'
!      write(99,*)
!      do i=1,dim_baser
!        write(99,'(1000f15.10)')(d_matc(i,j),j=1,dim_baser)
!      enddo
!     endif


!  lwork=26*dim_baser
!  allocate(work(26*dim_baser))
!  allocate(e(dim_baser))

! call DSYEV('V','L',dim_baser,d_matc,dim_baser,e,WORK,LWORK,INFO )
! write(98,'(1000f10.5)')(e(i),i=1,dim_baser)

! deallocate(work,e)
! deallocate(d_matc)


end subroutine dmatrix
!**********************************************************************
subroutine dmat_sp_subset(dim_base,dim_baser,dim_ind,phonbs,d_mat,ind_red,nx,s,vt,ns)

implicit none 

include 'types_eqm_hole.inc'

integer :: dim_base,dim_baser,dim_spur,dim_ind,i,j,k,ii,jj,ns
integer, dimension (:), allocatable :: ind_spur 
type(phonbase_typ), dimension (:), allocatable :: phonbs
double precision, dimension (:,:),allocatable :: d_mat,d_mat_sp
double precision, dimension (:), allocatable :: s,work
double precision, dimension (:,:), allocatable :: u,vt,a
integer, dimension (:), allocatable :: ind_red,nx
integer :: lwork,info

double precision :: xpom

allocate(ind_spur(dim_base))

ns=0

ii=0
do i=1,dim_base
  if (phonbs(i)%spur.eq.1) then 
    ii=ii+1
    ind_spur(ii)=i
  endif 
enddo

dim_spur=ii
write(*,*)' Dimension of spurious subspace ',dim_spur

if (dim_spur ==0) then 
allocate(vt(dim_ind,dim_ind))
vt=0.d0
 do i=1,dim_ind
   vt(i,i)=1.0d0
 enddo
endif 

if (dim_spur >0 ) then 

allocate(d_mat_sp(dim_spur,dim_base))
d_mat_sp=0.d0

do i=1,dim_spur
 jj=0
 do j=1,dim_baser
  if (nx(j).eq.1) then 
   jj=jj+1
   d_mat_sp(i,jj)=d_mat(j,ind_spur(i))
  endif
  enddo
enddo

write(*,*)' Independent rows  ', jj

!    SVD
allocate(s(dim_spur))
allocate(u(dim_spur,dim_spur),vt(dim_ind,dim_ind))
lwork=10*dim_baser
allocate(work(lwork))

call dgesvd('A','A', dim_spur, dim_ind, d_mat_sp, dim_spur, s, u, dim_spur, vt, dim_ind, work, lwork, info)


write(*,*) 'SVD  info =',info

write(999,*)
write(999,*)' SVD '
write(999,*)(s(i),i=1,dim_spur)
do i=1,dim_ind
 write(999,'(1000f10.5)')(vt(i,j),j=1,dim_ind)
enddo


write(*,*)'Test SVD'

allocate(a(dim_spur,dim_ind))

call dgemm('N','T',dim_spur,dim_ind,dim_ind,1.d0,d_mat_sp,dim_spur,vt,dim_ind,0.d0,a,dim_spur)

write(999,*)'Test SVD'

do i=1,dim_spur
 write(999,'(1000f10.5)')(a(i,j),j=1,dim_ind)
enddo

ii=0
do i=1,dim_spur
if (dabs(s(i)).gt.0.000001d0) then 
 ii=ii+1
endif
enddo
ns=ii
write(*,*)' Number of nonzero dg. values ',ns

! nulify 

do i=1,ns
  vt(i,:)=0.d0
enddo

!do i=1,ns
!vt(i,i)=1.0
!enddo


endif 


return

end subroutine dmat_sp_subset
!***********************************************************************
subroutine dmat_ind_set(d_mat,d_matc,dim_base,dim_baser,dim_ind,nx,ind_red)
use choleski

implicit none 

integer :: dim_base,dim_baser,i,j,lwork,info,iout,dim_ind,ns
double precision, dimension (:,:),allocatable :: d_mat,d_matc,d_inv
integer, dimension(:), allocatable :: nx,mxt
double precision, dimension(:), allocatable :: work,e
integer, dimension (:), allocatable :: ind_red
double precision, dimension (:), allocatable :: s
double precision, dimension (:,:), allocatable :: u,vt,a



  allocate(d_matc(dim_baser,dim_baser))
  do i=1,dim_baser
   do j=1,dim_baser
      d_matc(i,j)=d_mat(i,ind_red(j))
   enddo
  enddo




    iout=1
     if (iout.eq.1) then
      write(99,*)
      write(99,*)'******** matrix  D *************'
      write(99,*)
      do i=1,dim_baser
        write(99,'(1000f11.6)')(d_matc(i,j),j=1,dim_baser)
      enddo
     endif

!  transformation to new basis
!  ns=2

!  vt(1,:)=0.d0
!  vt(2,:)=0.d0
  
!  allocate(a(dim_baser,dim_baser))
!  a=d_matc
!  d_matc=0.d0
!  call dgemm('N','T',dim_baser,dim_baser,dim_baser,1.d0,a,dim_baser,vt,dim_baser,0.d0,d_matc,dim_baser)
!  allocate(u(dim_baser,dim_baser))
!  a=0.d0
!  call dgemm('N','N',dim_baser,dim_baser,dim_baser,1.d0,vt,dim_baser,d_matc,dim_baser,0.d0,a,dim_baser)
 
!   deallocate(d_matc)
!   allocate(d_matc(dim_baser-ns,dim_baser-ns))
!   d_matc=0.d0

!   ns=0  
 


!   do i=1,dim_baser-ns
!    do j=1,dim_baser-ns
!        d_matc(i,j)=a(i+ns,j+ns)
!     enddo
!    enddo

!    dim_baser=dim_baser-ns

!     if (iout.eq.1) then
!      write(99,*)
!      write(99,*)'******** matrix  D transformed *************'
!      write(99,*)
!      do i=1,dim_baser
!        write(99,'(1000f11.6)')(d_matc(i,j),j=1,dim_baser)
!      enddo
!     endif

!  a=d_matc 


  lwork=26*dim_baser
  allocate(work(26*dim_baser))
  allocate(e(dim_baser))

 call DSYEV('V','L',dim_baser,d_matc,dim_baser,e,WORK,LWORK,INFO )

 open(98,file='dmat_egv.log',status='unknown',form='formatted')
 write(98,'(1000f10.5)')(e(i),i=1,dim_baser)

 j=0
 do i=1,dim_baser
  if (dabs(e(i)).gt.1.0d-10) j=j+1
 enddo

 dim_ind=j
 write(*,*)'Number of independent spates =',dim_ind

!  d_matc=a

  do i=1,dim_baser
   do j=1,dim_baser
      d_matc(i,j)=d_mat(i,ind_red(j))
   enddo
  enddo

  call cholesk(dim_baser,dim_ind,dim_baser,d_matc,nx,mxt,d_inv)

  do i=1,dim_baser
   do j=1,dim_baser
      d_matc(i,j)=d_mat(i,ind_red(j))
   enddo
  enddo



 deallocate(work,e)


return
end subroutine dmat_ind_set

!*************************************************************************
subroutine amatrix(j_tot,dim_base,dim_baser,phonbs,lev,phon,ih_lev,ip_lev,hol_lev,par_lev,i_ph_int,nlev,ind_red,a_mat)
implicit none

include 'types_eqm_hole.inc'

type(level_typ),dimension(:), allocatable :: lev
type(phon_typ), dimension (:), allocatable :: phon
type(phonbase_typ), dimension (:), allocatable :: phonbs
type(rho_typ), dimension(:), allocatable :: ronp,ropp,ronh,roph

integer :: i,ired,j,dim_base,dim_baser,isi,j_tot,iout,ii,iii,ih_lev,ip_lev,nlev
integer :: i_sp,j_sp,i_lam,j_lam,i_lam_old
integer :: ihp_min,ihp_max,ihn_min,ihn_max,ipp_min,ipp_max,ipn_min,ipn_max,isi_max,ndla_max
integer :: jmin, jmax,nronp,nropp,nronh,nroph
integer :: lwork,info
integer :: i1,i2,ji1,ji2
integer :: ifaz

integer, dimension(:,:), allocatable :: i_ph_int
integer, dimension(:), allocatable :: par_lev, hol_lev
integer, dimension (:), allocatable :: ind_red

double precision, dimension (:,:),allocatable :: a_mat
double precision, dimension(:,:,:,:),allocatable :: ronh1
double precision :: rac
double precision, dimension(:), allocatable :: work,e

double precision, dimension(:,:,:,:,:),allocatable ::fp,fpn

integer, dimension (:,:,:), allocatable :: ipozbr
integer, dimension (:,:), allocatable :: ndbr


character*30 namefp,namefpn

allocate(a_mat(dim_baser,dim_base))
a_mat=0.0d0
i_lam_old=-10
isi_max=15

call read_input_par(ihp_min,ihp_max,ihn_min,ihn_max,ipp_min,ipp_max,ipn_min,ipn_max,ndla_max)

namefp='fmat.dat'

!     loads F(p) or F(n) interaction
jmin=0
jmax=isi_max

call readfin(namefp,jmin,jmax,1,nlev,1,nlev,fp)
  do ired=1,dim_baser
!  do i=1,dim_base
     i=ind_red(ired)

     i_sp=phonbs(i)%isp
     i_lam=phonbs(i)%ila

     if (i_lam.ne.i_lam_old) then
        call readro('1f_rp.dat',i_lam,ronp,nronp)
!        call readro('1f_rpp.dat',i_lam,ropp,nropp)
        call readro('1f_rh.dat',i_lam,ronh,nronh)
!        call readro('1f_rph.dat',i_lam,roph,nroph)
        i_lam_old=i_lam
        call redrsum(ndla_max,ronp,nronp,ronh,nronh,ipozbr,ndbr)
     endif

   do j=1,dim_base
     j_sp=phonbs(j)%isp
     j_lam=phonbs(j)%ila


      if (i_sp.eq.j_sp.and.i_lam.eq.j_lam) then 
        a_mat(ired,j)=a_mat(ired,j)-lev(i_sp)%en+phon(i_lam)%enf
      
      endif 

       do iii=1,ndbr(1,j_lam)  !nronp          ! neutron particle
         ii=ipozbr(1,j_lam,iii)
          i1=ronp(ii)%i1
          i2=ronp(ii)%i2
          isi=ronp(ii)%j
          ji1=lev(i1)%j
          ji2=lev(i2)%j

          ifaz=(lev(i_sp)%j+lev(j_sp)%j)/2-isi+(ji1-ji2)/2+phon(j_lam)%j+phon(i_lam)%j
          ifaz=(-1)**ifaz
          rac=racah(2*phon(i_lam)%j,2*isi,j_tot,lev(j_sp)%j,2*phon(j_lam)%j,lev(i_sp)%j)
          a_mat(ired,j)=a_mat(ired,j)+rac*ifaz*(2*isi+1)**0.5d0*ronp(ii)%ro*fp(isi,j_sp,i_sp,i2,i1)
        enddo

       do iii=1,ndbr(2,j_lam) ! nronh          ! neutron hole
         ii=ipozbr(2,j_lam,iii)
          i1=ronh(ii)%i1
          i2=ronh(ii)%i2
          ji1=lev(i1)%j
          ji2=lev(i2)%j
          isi=ronh(ii)%j

          ifaz=(lev(i_sp)%j+lev(j_sp)%j)/2-isi+(ji1-ji2)/2+phon(j_lam)%j+phon(i_lam)%j
          ifaz=(-1)**ifaz
          rac=racah(2*phon(i_lam)%j,2*isi,j_tot,lev(j_sp)%j,2*phon(j_lam)%j,lev(i_sp)%j)
          a_mat(ired,j)=a_mat(ired,j)+rac*ifaz*(2*isi+1)**0.5d0*0.5d0*ronh(ii)%ro*fp(isi,j_sp,i_sp,i2,i1)
        enddo

   enddo
  enddo


   iout=0
     if (iout.eq.1) then
      write(99,*)
      write(99,*)'******** matrix  A *************'
      write(99,*)
      do i=1,dim_base
        write(99,'(1000f10.5)')(a_mat(i,j),j=1,dim_base)
      enddo
     endif


end subroutine amatrix
! **********************************************************************
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

!******************************************************************************
subroutine reduce_rows(nx,mat,ndim,ndim_red)
implicit none
double precision, dimension (:,:),allocatable :: mat,matc
integer, dimension(:), allocatable :: nx
integer i,j,ii,jj,ndim,ndim_red

allocate(matc(ndim,ndim))
matc=mat
deallocate(mat)
allocate(mat(ndim_red,ndim))


ii=0
do i=1,ndim
 if (nx(i).ne.0) then
   ii=ii+1
!   jj=0
   do j=1,ndim
!   if (nx(j).ne.0) then
!       jj=jj+1
       mat(ii,j)=matc(i,j)
!     endif
    enddo

  endif
enddo
deallocate(matc)

end subroutine reduce_rows

!******************************************************************************


 subroutine  read_input_par(ihp_min,ihp_max,ihn_min,ihn_max,ipp_min,ipp_max,ipn_min,ipn_max,ndla_max)

      implicit double precision (a-h,o-z)

      character*30 name1f

      open(1,file='input_tda_coup.dat',status='old',form='formatted')

      read(1,'(30i6)')ia,iz
      read(1,'(30i6)')ihn_min,ihn_max
      read(1,'(30i6)')ihp_min,ihp_max
      read(1,'(30i6)')ipn_min,ipn_max
      read(1,'(30i6)')ipp_min,ipp_max
!      read(1,26)alfa,beta
!      read(1,*)
!      read(1,15)iparmn,iparmx
!      read(1,15)jminn,jmaxn
      close(1)

!      isi_max=jmaxn

      name1f='1phonon/1f_states.dat'
      open (3,file=name1f,status='old',form='unformatted')

      do while (.not.eof(3))
       read(3)i,ipar,ijj,en
      enddo


      ndla_max=i

      close(3)

      return
 end subroutine read_input_par
!***************************************************************************
subroutine normal(idim_ind,idim_base,idim_baser,vr,xamp,ind_red,nx,ijj)

      implicit double precision (a-h,o-z)

!c      include 'chole.inc'

      double precision, dimension(:,:), allocatable :: vr,xamp,xampr
      integer, dimension (:), allocatable :: ind_red,nx

      allocate(xampr(idim_ind,idim_ind))
      
      do j=1,idim_ind
      ii=0
      do i=1,idim_baser
        if (nx(i).eq.1) then
          ii=ii+1
          xampr(ii,j)=xamp(ind_red(i),j)
        endif
      enddo
     enddo

      

      xfact=-1.d0*(dfloat(ijj+1))**0.5d0

      do i=1,idim_ind
        xpom=0.d0
        do j=1,idim_ind
          xpom=xpom+vr(j,i)*xampr(j,i)
         enddo
        do j=1,idim_ind
          vr(j,i)=vr(j,i)/dsqrt(xpom)
        enddo
        do j=1,idim_base
          xamp(j,i)=xamp(j,i)/dsqrt(xpom)
!          xampr(j,i)=xampr(j,i)/dsqrt(xpom)
        enddo
        write(771,*)' norm i =',i,xpom
      enddo

     xamp=xamp*xfact
     do j=1,idim_ind
      ii=0
      do i=1,idim_baser
        if (nx(i).eq.1) then
          ii=ii+1
          xampr(ii,j)=xamp(ind_red(i),j)
        endif
      enddo
     enddo



!     xampr=xampr*xfact
!    test of normalization
   
   do k=1,idim_ind
    do i=1,idim_ind
     xpom=0.0d0
     do j=1,idim_ind
       xpom=xpom+vr(j,k)*xampr(j,i)
     enddo 
     if (dabs(xpom).gt.1.d-10) write(771,*) 'ovrl ',k,i,xpom,xpom/xfact
    enddo
   enddo
    
      return
 end subroutine normal
!**************************************************************

! calculatesm racah w

double precision function racah(a,b,c,d,e,f) 

 implicit none 
 integer, intent(in) :: a,b,c,d,e,f
 integer :: ifaz
 double precision :: sixj_sym

 sixj_sym=sixj(a,b,e,d,c,f)
 ifaz=(a+b+c+d)/2
 ifaz=(-1)**ifaz
 racah=sixj_sym*dfloat(ifaz)

end function racah 



end module dmatr


