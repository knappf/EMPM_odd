  module admatr

   contains

   subroutine admat(nf,ipar,jcal,phonbs,ndim,no,phon1,phon2,mxt,myid,numprocs)

   use anglib
   use input_sp

   implicit double precision (a-h,o-z)


   include 'types_eqm.inc'

   integer, dimension (:), allocatable :: phonus, phonmus
   type(phon_typ), dimension (:), allocatable :: phon1,phon2
   type(phonbase_typ), dimension (:), allocatable :: phonbs
   type(rho_typ), dimension(:), allocatable :: rop,roh
   type(level_typ),dimension(:), allocatable :: lev

   double precision, dimension(:,:,:,:),allocatable :: rop1,roh1

   double precision, dimension(:,:,:,:), allocatable :: fph

   double precision, dimension(:,:,:,:,:), allocatable :: csixj,csixjd

   integer, dimension(:), allocatable :: nx,mxt,mxtr,ig_resh

   integer, dimension (:,:,:), allocatable :: ipozbr
   integer, dimension (:,:), allocatable :: ndbr
   double precision, dimension (:,:), allocatable :: rtdaovr
   double precision, dimension (:), allocatable :: a_m,d_m

   character*30 namernp,namerpp,namernh,namerph,namefp,namefh,namernp1,namerpp1,namernh1,namerph1
   character*15 row_number
   integer myid,numprocs


   allocate(rtdaovr(0:5000,0:5000))

   rtdaovr=0.0d0
   ilamcm=0
   iii=0

!    open(881,file='tda_r_overl.dat',status='unknown',form='formatted')
!     do while (.not.eof(881))
!      read(881,*)iii,eee,rr
!    enddo
!    close(881)


!      ilamcm=iii
!      write(*,*)' CM phonon is number ',ilamcm

!     open(881,file='tda_r_overl.dat',status='unknown',form='formatted')
!      do while (.not.eof(881))
!       read(881,*)iii,eee,rr
!       rtdaovr(iii,ilamcm)=rr
!       rtdaovr(ilamcm,iii)=rr
!      enddo
!     close(881)

!      phon1(ilamcm)%enf=rtdaovr(ilamcm,ilamcm)
!      phon2(ilamcm)%enf=rtdaovr(ilamcm,ilamcm)

    open(2,file='mxtr.dat',status='unknown',form='unformatted')
     read(2)idphontr
      allocate(mxtr(idphontr))
     read(2)(mxtr(i),i=1,idphontr)
    close(2)
 
    no=idphontr
    call loadsp(lev,nlev,isi_max)


    open(6,file='a_mat.dat',status='unknown',form='unformatted')
    open(7,file='d_mat.dat',status='unknown',form='unformatted')


      namefp='V_phon_par'
      namefh='V_phon_hol'

      call read_input_par(ihp_min,ihp_max,ihn_min,ihn_max,ipp_min,ipp_max,ipn_min,ipn_max,isi_max,ndla_max)

!      write(*,*)ihp_min,ihp_max,ihn_min,ihn_max,ipp_min,ipp_max,ipn_min,ipn_max,isi_max,ndla_max
      ndla=ndla_max
      nsi=isi_max
      ndla2=ndla_max

      iaaold=-10
      iaold=-10

      jmx=nsi !20

      allocate(csixj(0:jmx,0:jmx,0:jmx,0:jmx,0:jmx))
      csixj=0.0d0

      do i=0,jmx
       do j=0,jmx
        do k=0,jmx
         do l=0,jmx
          do m=0,jmx
           csixj(i,j,k,l,m)=sixj(2*jcal,2*i,2*j,2*k,2*l,2*m)
          enddo
         enddo
        enddo
       enddo
      enddo

     allocate(csixjd(0:jmx,0:jmx,0:jmx,0:jmx,0:jmx))
      csixjd=0.0d0

      do i=0,jmx
       do j=0,jmx
        do k=0,jmx
         do l=0,jmx
          do m=0,jmx
           csixjd(i,j,k,l,m)=sixj(2*i,2*j,2*k,2*jcal,2*l,2*m)
          enddo
         enddo
        enddo
       enddo
      enddo


      allocate(a_m(ndim),d_m(ndim))

      do i=1,idphontr
       iii=mxtr(i)
       iaa=phonbs(iii)%ilap  ! 1phonon index
       ia=phonbs(iii)%ila    ! n-1 phonon index
      enddo

      allocate(ig_resh(idphontr))
      do i=1,idphontr
       ig_resh(i)=i
      enddo

      call rperm(idphontr,ig_resh)

      if (mod(idphontr,numprocs).eq.0) then
            n_seg=idphontr/numprocs
      else
            n_seg=idphontr/numprocs+1
      endif

      if (myid.eq.0) write(*,*) ' size of segment = ',n_seg


      do irs=myid*n_seg+1,min((myid+1)*n_seg,idphontr)
      i=ig_resh(irs)

      write(row_number,'(i15.15)')i
      open(6,file='./scratch/a_mat_'//row_number,status='unknown',form='unformatted')
      open(7,file='./scratch/d_mat_'//row_number,status='unknown',form='unformatted')

!      do i=1,idphontr
       a_m=0.d0
       d_m=0.d0

       iii=mxtr(i)
       iaa=phonbs(iii)%ilap  ! 1phonon index
       ia=phonbs(iii)%ila    ! n-1 phonon index
       jila=phon1(iaa)%j
       jial=phon2(ia)%j

       if (ia.ne.iaold) then
            call readro('1f_rp.dat',ia,rop,nrop)
            call readro('1f_rh.dat',ia,roh,nroh)
            iaold=ia
            call redrsum(ndla,rop,nrop,roh,nroh,ipozbr,ndbr)

        endif

      if (iaa.ne.iaaold) then
      if (.not.allocated(fph)) allocate(fph(ndla,0:nsi,1:nlev,1:nlev))

        fph=0.0d0
        call readfn('V_phon_par',iaa,fph,ndla,1,nlev,nsi)
        call readfn('V_phon_hol',iaa,fph,ndla,1,nlev,nsi)
!        iaaold=iaa
!       endif
!  
!       if (iaa.ne.iaaold) then

      if (.not.allocated(rop1)) allocate(rop1(ndla,0:nsi,1:nlev,1:nlev))
    
        rop1=0.d0
        call readro11('1f_rp.dat',iaa,rop1,ndla,1,nlev,nsi)
        call readro11('1f_rh.dat',iaa,rop1,ndla,1,nlev,nsi)

        iaaold=iaa
       endif

!       call OMP_SET_NUM_THREADS(6)
!$omp parallel default(shared) private(jjj,dd,aa, &
!$omp xsixj,ii,iii,ib,ibb,jlap,jilap,jialp,ibt,i1,i2, &
!$omp isi,faz,itid,ji1,ji2)

!$omp do schedule (dynamic)


       do j=1,ndim

        jjj=j

        aa=0.d0
        dd=0.d0

        ib=phonbs(jjj)%ila
        ibb=phonbs(jjj)%ilap

        jilap=phon1(ibb)%j
        jialp=phon2(ib)%j
  
        tzfact=1.0d0
!       if ((phon1(ibb)%tz+phon2(ib)%tz).eq.(phon1(iaa)%tz+phon2(ia)%tz)) tzfact=1.d0

        if (ia.eq.ib.and.iaa.eq.ibb) aa=aa+phon1(iaa)%enf+phon2(ia)%enf
!       if (ipar.eq.-1.and.jcal.eq.1) then

!          if (ia.eq.ib.and.iaa.eq.ilamcm) then
!              aa=aa+rtdaovr(ibb,iaaa)
!          endif

          if (ia.eq.ib.and.iaa.ne.ibb) then
!              write(*,*)'***'
             aa=aa+rtdaovr(iaa,ibb)
          endif
          if (iaa.eq.ibb.and.ia.ne.ib) then
!              write(*,*)'***'
             aa=aa+rtdaovr(ia,ib)
          endif

        if (nf.eq.2) then

          if (ia.eq.ib.and.iaa.eq.ibb) dd=dd+1.d0
          jlap=phon1(ibb)%j
          jla=phon1(iaa)%j
!          xsixj=sixj(2*jla,0,2*jla,2*jlap,2*jcal,2*jlap)
          xsixj=csixjd(0,jla,jla,jlap,jlap)
!c          faz=dfloat((-1)**((ji1-ji2)/2))
          if ((ia.eq.ibb).and.(ib.eq.iaa)) then
            dd=dd+xsixj*(dfloat(2*jla+1)*(2*jlap+1))**0.5d0
          endif
        endif

         faz1=-1.d0*dfloat((-1)**(jcal+phon1(ibb)%j+phon2(ia)%j))



        do iii=1,ndbr(1,ib)  !nronp          ! neutron particle
         ii=ipozbr(1,ib,iii)
!         ibt=ronp(ii)%ilap
!         if (ibt.eq.ib) then
          i1=rop(ii)%i1
          i2=rop(ii)%i2
          isi=rop(ii)%j
          ji1=lev(i1)%j
          ji2=lev(i2)%j

!          xsixj=sixj(2*jcal,2*jilap,2*jialp,2*isi,2*jial,2*jila)
          xsixj=csixj(jilap,jialp,isi,jial,jila)
!          write(9911,*)xsixj
          faz=dfloat((-1)**(jcal+jial))
          aa=aa+xsixj*faz*rop(ii)%ro*fph(ibb,isi,i1,i2)

          xsixj=csixjd(isi,jilap,jila,jial,jialp)
          faz=dfloat((-1)**((ji1-ji2)/2))
          dd=dd+faz1*faz*xsixj*rop(ii)%ro*rop1(ibb,isi,i2,i1)

!          write(993,*)xsixj,faz,ronp(ii)%ro,fnp(ibb,isi,i1,i2)
!         endif
        enddo


        do iii=1,ndbr(2,ib) ! nronh          ! neutron hole
         ii=ipozbr(2,ib,iii)
!         ibt=ronh(ii)%ilap
!         if (ibt.eq.ib) then
          i1=roh(ii)%i1
          i2=roh(ii)%i2

          ji1=lev(i1)%j
          ji2=lev(i2)%j

          isi=roh(ii)%j
!          xsixj=sixj(2*jcal,2*jilap,2*jialp,2*isi,2*jial,2*jila)
          xsixj=csixj(jilap,jialp,isi,jial,jila)

          faz=dfloat((-1)**(jcal+jial))
          aa=aa+xsixj*faz*roh(ii)%ro*fph(ibb,isi,i1,i2)
!          write(993,*)xsixj,faz,ronh(ii)%ro,fnh(ibb,isi,i1,i2)

          xsixj=csixjd(isi,jilap,jila,jial,jialp)
          faz=dfloat((-1)**((ji2-ji1)/2))
          dd=dd+faz1*faz*xsixj*roh(ii)%ro*rop1(ibb,isi,i2,i1)

!         endif
        enddo

!        write(993,*)'**************'


!        if (dabs(dd).gt.1.d-10) write(999,*)i,j,dd
!        if (dabs(aa).gt.1.d-10) write(6)i,j,aa
!        if (dabs(dd).gt.1.d-10) write(7)i,j,dd
       
        a_m(j)=aa
        d_m(j)=dd*tzfact

       enddo  ! over j

!$omp end do

!$omp end parallel

       
       write(6)(a_m(jjj),jjj=1,ndim)
       write(7)(d_m(jjj),jjj=1,ndim)

      enddo  ! over i

      deallocate(a_m,d_m)

      close(6)
      close(7)



      return

      end subroutine admat

!*****************************************************************************
   subroutine readro(fname,ig,ron,ndgg)

      implicit double precision (a-h,o-z)

      include 'formats_eqm.inc'
      include 'types_eqm.inc'

      type(rho_typ), dimension(:), allocatable :: ron

      character(len=9)fname
      character(len=4)nlam
      logical je_tam

      ifile=33
 
      write(nlam,'(i4.4)')ig

      inquire(file='scratch/'//fname//'_'//nlam,exist=je_tam)

      if (je_tam.eq..FALSE.) then
        ndgg=0
        return
      endif


      open(ifile,file='scratch/'//fname//'_'//nlam,status='unknown',form='unformatted')

      ndro=5000000
      ndgg=0

      if (.not.allocated(ron)) allocate (ron(ndro))

       read(ifile)igg,ndgg

       if (igg.ne.ig) then
               write(*,*)' Loaded Ig does not match !!! '
               stop
       endif

       if (ndgg.gt.ndro) then
                write(*,*)'WARNING: Increase dimension in readro'
                stop
       endif

       read(ifile)(ron(ii)%ilap,ii=1,ndgg)
       read(ifile)(ron(ii)%j,ii=1,ndgg)
       read(ifile)(ron(ii)%i1,ii=1,ndgg)
       read(ifile)(ron(ii)%i2,ii=1,ndgg)
       read(ifile)(ron(ii)%ro,ii=1,ndgg)


       close(ifile)
      return
      end subroutine readro

!***********************************************************************
      subroutine redrsum(ifonmx,rop,nrop,roh,nroh,ipozbr,ndbr)

      implicit double precision (a-h,o-z)

      include 'formats_eqm.inc'
      include 'types_eqm.inc'

      type(rho_typ), dimension(:), allocatable :: rop,roh
      integer, dimension (:,:,:), allocatable :: ipozbr
      integer, dimension (:,:), allocatable :: ndbr

      character(len=30)fname

      ndipo=50000

      if (.not.allocated(ipozbr)) allocate(ipozbr(2,ifonmx,ndipo))
      ipozbr=0

      if (.not.allocated(ndbr)) allocate(ndbr(2,ifonmx))
      ndbr=0



      do i=1,nrop
       ibt=rop(i)%ilap
       ndbr(1,ibt)=ndbr(1,ibt)+1

        if (ndbr(1,ibt).gt.ndipo) then
               write(*,*)' Increase dimension in redrsum'
               stop
            endif
       ipozbr(1,ibt,ndbr(1,ibt))=i

      enddo


      do i=1,nroh
       ibt=roh(i)%ilap
       ndbr(2,ibt)=ndbr(2,ibt)+1

        if (ndbr(2,ibt).gt.ndipo) then
               write(*,*)' Increase dimension in redrsum'
               stop
            endif
       ipozbr(2,ibt,ndbr(2,ibt))=i

      enddo


      return
      end subroutine redrsum
!***********************************************************************


      subroutine readro11(fname,ig,ron,ndla,n1mn,n1mx,nsi)

      implicit double precision (a-h,o-z)

      include 'types_eqm.inc'

      double precision, dimension(:,:,:,:), allocatable :: ron
      type(rho_typ), dimension(:), allocatable :: ronn
      logical je_tam_subor



      character(len=9)fname
      character(len=4)nlam

      ifile=33
      write(nlam,'(i4.4)')ig

      inquire(file='scratch/'//fname//'_'//nlam,exist=je_tam_subor)

      if (je_tam_subor.eq..FALSE.) then
        ndgg=0
        return
      endif


      open(ifile,file='scratch/'//fname//'_'//nlam,status='unknown',form='unformatted')


      ndro=5000000
      ndgg=0

      if (.not.allocated(ronn)) allocate (ronn(ndro))


       read(ifile)igg,ndgg

       if (igg.ne.ig) then
               write(*,*)' Loaded Ig does not match !!! '
               stop
       endif

       if (ndgg.gt.ndro) then
                write(*,*)'WARNING: Increase dimension in readro',ndgg,ndro
                stop
       endif

       read(ifile)(ronn(ii)%ilap,ii=1,ndgg)
       read(ifile)(ronn(ii)%j,ii=1,ndgg)
       read(ifile)(ronn(ii)%i1,ii=1,ndgg)
       read(ifile)(ronn(ii)%i2,ii=1,ndgg)
       read(ifile)(ronn(ii)%ro,ii=1,ndgg)

      close(ifile)

!      if (.not.allocated(ron)) then 
!       allocate(ron(ndla,0:nsi,n1mn:n1mx,n1mn:n1mx))
!       ron=0.d0
!      endif 

       do ii=1,ndgg

        ibt=ronn(ii)%ilap
        isit=ronn(ii)%j
        i1t=ronn(ii)%i1
        i2t=ronn(ii)%i2
        rot=ronn(ii)%ro

        if (ibt.gt.ndla.or.isit.gt.nsi.or.i1t.lt.n1mn.or.i1t.gt.n1mx.or.i2t.lt.n1mn.or.i2t.gt.n1mx) then
         write(*,*)' Small dimensions of ron in read11',ibt,isit,i1t,i2t
         stop
       endif

        ron(ibt,isit,i1t,i2t)=rot

       enddo



      return
      end subroutine readro11


!******************************************************************************
     subroutine readfn(fname,ig,ron,ndla,n1mn,n1mx,nsi)

      implicit double precision (a-h,o-z)

      include 'types_eqm.inc'

      double precision, dimension(:,:,:,:), allocatable :: ron
      type(rho_typ), dimension(:), allocatable :: ronn

      character(len=10)fname
      character(len=4)nlam

      ifile=33

      write(nlam,'(i4.4)')ig
 
      open(ifile,file='scratch/'//fname//'_'//nlam,status='unknown',form='unformatted')

      ndro=50000000
      ndgg=0

      if (.not.allocated(ronn)) allocate (ronn(ndro))

       read(ifile)igg,ndgg

       if (igg.ne.ig) then
               write(*,*)' Loaded Ig does not match !!! '
               stop
       endif

       if (ndgg.gt.ndro) then
                write(*,*)'WARNING: Increase dimension in readfn',ndgg,ndro
                stop
       endif

       read(ifile)(ronn(ii)%ilap,ii=1,ndgg)
       read(ifile)(ronn(ii)%j,ii=1,ndgg)
       read(ifile)(ronn(ii)%i1,ii=1,ndgg)
       read(ifile)(ronn(ii)%i2,ii=1,ndgg)
       read(ifile)(ronn(ii)%ro,ii=1,ndgg)
   
       close(33)


!      if (.not.allocated(ron)) then 
!           allocate(ron(ndla,0:nsi,n1mn:n1mx,n1mn:n1mx))
!           ron=0.d0
!      endif

       do ii=1,ndgg

        ibt=ronn(ii)%ilap
        isit=ronn(ii)%j
        i1t=ronn(ii)%i1
        i2t=ronn(ii)%i2
        rot=ronn(ii)%ro

        if (ibt.gt.ndla.or.isit.gt.nsi.or.i1t.lt.n1mn.or.i1t.gt.n1mx.or.i2t.lt.n1mn.or.i2t.gt.n1mx) then
         write(*,*)' Small dimensions of ron in readfn',ibt,isit,i1t,i2t
         stop
       endif

        ron(ibt,isit,i1t,i2t)=rot

       enddo

       return

       end subroutine readfn

!***********************************************************************

      subroutine readf(ifile,ig,ron,ndla,n1mn,n1mx,n2mn,n2mx,nsi)

      implicit double precision (a-h,o-z)

      include 'types_eqm.inc'

      double precision, dimension(:,:,:,:), allocatable :: ron

      character(len=30)fname


      if (.not.allocated(ron)) allocate(ron(ndla,0:nsi,n1mn:n1mx,n2mn:n2mx))


  111 continue
       ibt=1
       noig=0

       ron=0.d0

  112  read(ifile)igg

       if(igg.eq.0) goto 112


       do while (ibt.ne.0)
       read(ifile)ibt,isit,i1t,i2t,rot

       if (igg.eq.ig) then

       if (ibt.ne.0) then
       if (ibt.gt.ndla.or.isit.gt.nsi.or.i1t.lt.n1mn.or.i1t.gt.n1mx.or.i2t.lt.n2mn.or.i2t.gt.n2mx) then
         write(*,*)' Small dimensions of ron in read1',ibt,isit,i1t,i2t
         stop
       endif

        ron(ibt,isit,i1t,i2t)=rot
       endif

       endif

       enddo

        if (ig.ne.igg) then
!c            write(*,*)'WARNING: Loaded ig does not match '
            if (igg.gt.ig) rewind(ifile)
            goto 111
       endif

      return
      end subroutine readf
!***********************************************************************


      subroutine loadsp(lev,nlev,jmax)

      use input_sp

      implicit double precision (a-h,o-z)

      include 'types_eqm.inc'
      include 'formats_eqm.inc'

      type(level_typ),dimension(:), allocatable :: lev


!c      write(*,*)'Loading of input '


      open(1,file='input_tda_coup.dat',status='old',form='formatted')

      read(1,15)ia,iz
      read(1,15)ihnmn,ihnmx
      read(1,15)ihpmn,ihpmx
      read(1,15)ipnmn,ipnmx
      read(1,15)ippmn,ippmx
!      read(1,26)alfa,beta
!      read(1,*)
!      read(1,15)iparmn,iparmx
!      read(1,15)jminn,jmaxn
      close(1)

      allocate(lev(1000))


      call inp_sp(lev,nlev,jmax)
      

      
      end subroutine loadsp
!
      subroutine  read_input_par(ihp_min,ihp_max,ihn_min,ihn_max,ipp_min,ipp_max,ipn_min,ipn_max,isi_max,ndla_max)


      include 'formats_eqm.inc'

      character*30 name1f


      open(1,file='input_tda_coup.dat',status='old',form='formatted')

      read(1,15)ia,iz
      read(1,15)ihn_min,ihn_max
      read(1,15)ihp_min,ihp_max
      read(1,15)ipn_min,ipn_max
      read(1,15)ipp_min,ipp_max
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


      end subroutine read_input_par

!***********************************************************************
!-----------------------------------------------------------------
subroutine rperm(N, p)

 integer(kind=4), intent(in) :: N
 integer(kind=4), dimension(:), intent(out) :: p

 integer(kind=4) :: i
 integer(kind=4) :: k, j, ipj, itemp, m
 real(kind=4), dimension(100) :: u

p = (/ (i, i=1,N) /)

! Generate up to 100 U(0,1) numbers at a time.
do i=1,N,100
m = min(N-i+1, 100)
call random_number(u)
do j=1,m
ipj = i+j-1
k = int(u(j)*(N-ipj+1)) + ipj
itemp = p(ipj)
p(ipj) = p(k)
p(k) = itemp
end do
end do
return

end subroutine rperm
!--------------------------------------------------------------


      end module admatr

