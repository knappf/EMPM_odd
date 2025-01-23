module dmatr

use anglib 
use read_densfmat

contains  

subroutine dmatrix(j_tot,dim_base,phonbs,levn,levp,phon,d_mat)


implicit none 

include 'types_eqm_hole.inc'

type(level_typ),dimension(:), allocatable :: levn,levp
type(phon_typ), dimension (:), allocatable :: phon
type(phonbase_typ), dimension (:), allocatable :: phonbs
integer :: i,j,dim_base,nronh,ii,isi,j_tot,iout
integer :: i_sp,j_sp,i_lam,j_lam,i_lam_old
integer :: ihp_min,ihp_max,ihn_min,ihn_max,ipp_min,ipp_max,ipn_min,ipn_max,isi_max,ndla_max
integer :: lwork,info
double precision, dimension (:,:),allocatable :: d_mat,d_matc
double precision, dimension(:,:,:,:),allocatable :: ronh1
double precision :: rac
double precision, dimension(:), allocatable :: work,e


allocate(d_mat(dim_base,dim_base))
d_mat=0.0d0
i_lam_old=-10
isi_max=15

call read_input_par(ihp_min,ihp_max,ihn_min,ihn_max,ipp_min,ipp_max,ipn_min,ipn_max,ndla_max)



  do i=1,dim_base 
     i_sp=phonbs(i)%isp
     i_lam=phonbs(i)%ila 
    
     if (i_lam.ne.i_lam_old) then 
          call readro11('1f_rnh.dat',i_lam,ronh1,ndla_max,ihn_min,ihn_max,ihn_min,ihn_max,isi_max)
          i_lam_old=i_lam
     endif 

   do j=1,dim_base 
     j_sp=phonbs(j)%isp
     j_lam=phonbs(j)%ila
     
     if (i_sp.eq.j_sp.and.i_lam.eq.j_lam) d_mat(i,j)=1.0d0


     do isi=0,isi_max
       d_mat(i,j)=d_mat(i,j)+ronh1(j_lam,isi,i_sp,j_sp)*(2*isi+1)**0.5d0*racah(2*isi,levn(j_sp)%j,2*phon(i_lam)%j,j_tot,levn(i_sp)%j,2*phon(j_lam)%j)
     enddo


   enddo
  enddo

    iout=1
     if (iout.eq.1) then
      write(99,*)
      write(99,*)'******** matrix  D *************'
      write(99,*)
      do i=1,dim_base
        write(99,'(1000f15.10)')(d_mat(i,j),j=1,dim_base)
      enddo
     endif

  allocate(d_matc(dim_base,dim_base))
  d_matc=d_mat

  lwork=26*dim_base
  allocate(work(26*dim_base))
  allocate(e(dim_base))

 call DSYEV('V','L',dim_base,d_matc,dim_base,e,WORK,LWORK,INFO )
 write(98,'(1000f15.10)')(e(i),i=1,dim_base)

 deallocate(work,e)
 deallocate(d_matc)


end subroutine dmatrix
!*************************************************************************
subroutine amatrix(j_tot,dim_base,phonbs,levn,levp,phon,a_mat)

implicit none

include 'types_eqm_hole.inc'

type(level_typ),dimension(:), allocatable :: levn,levp
type(phon_typ), dimension (:), allocatable :: phon
type(phonbase_typ), dimension (:), allocatable :: phonbs
type(rho_typ), dimension(:), allocatable :: ronp,ropp,ronh,roph

integer :: i,j,dim_base,isi,j_tot,iout,ii,iii
integer :: i_sp,j_sp,i_lam,j_lam,i_lam_old
integer :: ihp_min,ihp_max,ihn_min,ihn_max,ipp_min,ipp_max,ipn_min,ipn_max,isi_max,ndla_max
integer :: jmin, jmax,nronp,nropp,nronh,nroph
integer :: lwork,info
integer :: i1,i2,ji1,ji2
integer :: ifaz

double precision, dimension (:,:),allocatable :: a_mat
double precision, dimension(:,:,:,:),allocatable :: ronh1
double precision :: rac
double precision, dimension(:), allocatable :: work,e

double precision, dimension(:,:,:,:,:),allocatable ::fp,fpn

integer, dimension (:,:,:), allocatable :: ipozbr
integer, dimension (:,:), allocatable :: ndbr


character*30 namefp,namefpn

allocate(a_mat(dim_base,dim_base))
a_mat=0.0d0
i_lam_old=-10
isi_max=15

call read_input_par(ihp_min,ihp_max,ihn_min,ihn_max,ipp_min,ipp_max,ipn_min,ipn_max,ndla_max)

namefp='fmat_n.dat'
namefpn='fmat_pn.dat'

!     loads F(p) or F(n) interaction
jmin=0
jmax=isi_max

call readfin(namefp,jmin,jmax,ihn_min,ipn_max,ihn_min,ipn_max,fp)
!     loads F(pn)
call readfin(namefpn,jmin,jmax,ihp_min,ipp_max,ihn_min,ipn_max,fpn)

  do i=1,dim_base
     i_sp=phonbs(i)%isp
     i_lam=phonbs(i)%ila

     if (i_lam.ne.i_lam_old) then
        call readro('1f_rnp.dat',i_lam,ronp,nronp)
        call readro('1f_rpp.dat',i_lam,ropp,nropp)
        call readro('1f_rnh.dat',i_lam,ronh,nronh)
        call readro('1f_rph.dat',i_lam,roph,nroph)
        i_lam_old=i_lam
        call redrsum(ndla_max,ronp,nronp,ropp,nropp,ronh,nronh,roph,nroph,ipozbr,ndbr)
     endif

   do j=1,dim_base
     j_sp=phonbs(j)%isp
     j_lam=phonbs(j)%ila

      if (i_sp.eq.j_sp.and.i_lam.eq.j_lam) then 
        a_mat(i,j)=a_mat(i,j)-levn(i_sp)%en+phon(i_lam)%enf
      
      endif 

       do iii=1,ndbr(1,j_lam)  !nronp          ! neutron particle
         ii=ipozbr(1,j_lam,iii)
          i1=ronp(ii)%i1
          i2=ronp(ii)%i2
          isi=ronp(ii)%j
          ji1=levn(i1)%j
          ji2=levn(i2)%j

          ifaz=(levn(i_sp)%j+levn(j_sp)%j)/2-isi
          ifaz=(-1)**ifaz
          rac=racah(2*phon(i_lam)%j,2*isi,j_tot,levn(j_sp)%j,2*phon(j_lam)%j,levn(i_sp)%j)
          a_mat(i,j)=a_mat(i,j)+rac*ifaz*(2*isi+1)**0.5d0*ronp(ii)%ro*fp(isi,i_sp,j_sp,i1,i2)
        enddo

       do iii=1,ndbr(3,j_lam) ! nronh          ! neutron hole
         ii=ipozbr(3,j_lam,iii)
          i1=ronh(ii)%i1
          i2=ronh(ii)%i2
          ji1=levn(i1)%j
          ji2=levn(i2)%j
          isi=ronh(ii)%j

          ifaz=(levn(i_sp)%j+levn(j_sp)%j)/2-isi
          ifaz=(-1)**ifaz
          rac=racah(2*phon(i_lam)%j,2*isi,j_tot,levn(j_sp)%j,2*phon(j_lam)%j,levn(i_sp)%j)
          a_mat(i,j)=a_mat(i,j)+rac*ifaz*(2*isi+1)**0.5d0*0.5d0*ronh(ii)%ro*fp(isi,j_sp,i_sp,i1,i2)
        enddo



       do iii=1,ndbr(2,j_lam)        !nropp          ! proton particle
         ii=ipozbr(2,j_lam,iii)
          i1=ropp(ii)%i1
          i2=ropp(ii)%i2
          ji1=levp(i1)%j
          ji2=levp(i2)%j
          isi=ropp(ii)%j

          ifaz=(levn(i_sp)%j+levn(j_sp)%j)/2-isi
          ifaz=(-1)**ifaz
          rac=racah(2*phon(i_lam)%j,2*isi,j_tot,levn(j_sp)%j,2*phon(j_lam)%j,levn(i_sp)%j)
          a_mat(i,j)=a_mat(i,j)+rac*ifaz*(2*isi+1)**0.5d0*ropp(ii)%ro*fpn(isi,i_sp,j_sp,i1,i2)
        enddo


       do iii=1,ndbr(4,j_lam)  !nroph          ! proton hole
         ii=ipozbr(4,j_lam,iii)
          i1=roph(ii)%i1
          i2=roph(ii)%i2
          ji1=levp(i1)%j
          ji2=levp(i2)%j
          isi=roph(ii)%j

          ifaz=(levn(i_sp)%j+levn(j_sp)%j)/2-isi
          ifaz=(-1)**ifaz
          rac=racah(2*phon(i_lam)%j,2*isi,j_tot,levn(j_sp)%j,2*phon(j_lam)%j,levn(i_sp)%j)
          a_mat(i,j)=a_mat(i,j)+rac*ifaz*(2*isi+1)**0.5d0*roph(ii)%ro*fpn(isi,i_sp,j_sp,i1,i2)
        enddo

   enddo
  enddo

   iout=1
     if (iout.eq.1) then
      write(99,*)
      write(99,*)'******** matrix  A *************'
      write(99,*)
      do i=1,dim_base
        write(99,'(1000f15.10)')(a_mat(i,j),j=1,dim_base)
      enddo
     endif

  lwork=26*dim_base
  allocate(work(26*dim_base))
  allocate(e(dim_base))

! call DSYEV('V','L',dim_base,a_mat,dim_base,e,WORK,LWORK,INFO )
! write(98,'(1000f15.10)')(e(i),i=1,dim_base)

 deallocate(work,e)

end subroutine amatrix


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


! calculatesm racah w

real function racah(a,b,c,d,e,f) 

 implicit none 
 integer, intent(in) :: a,b,c,d,e,f
 integer :: ifaz
 real :: sixj_sym

 sixj_sym=sixj(a,b,e,d,c,f)
 ifaz=(a+b+c+d)/2
 ifaz=(-1)**ifaz
 racah=sixj_sym*ifaz

end function racah 



end module dmatr


