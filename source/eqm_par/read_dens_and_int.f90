module read_densfmat

contains 


!***********************************************************************

     subroutine readfin(fname,jmin,jmax,imin,imax,kmin,kmax,fpp)

      implicit double precision (a-h,o-z)

!      include 'formats_phon_int.inc'

      double precision, dimension(:,:,:,:,:), allocatable ::fpp

      character(len=30)fname
      integer(kind=1) :: j_f
      integer(kind=2) :: i_i,i_j,i_k,i_l


      allocate(fpp(jmin:jmax,imin:imax,imin:imax,kmin:kmax,kmin:kmax))
      fpp=0.d0

!      open(2,file=fname,status='old',form='formatted')
      open(2,file=fname,status='old',form='unformatted')

      do while (.not.eof(2))
!       read(2,10)itt,ipt,ijt,i,j,k,l,vint
!       read(2)itt,ipt,ijt,i,j,k,l,vint
      read(2)j_f,i_i,i_j,i_k,i_l,vint

      ijt=j_f
      i=i_i
      j=i_j
      k=i_k
      l=i_l

!      read(2)itt,ipt,ijt,i,j,k,l,vint

!c        if (ipt.ne.ipar) goto 11
      if (ijt.gt.jmax.or.ijt.lt.jmin) goto 11

      if (i.gt.imax.or.i.lt.imin.or.j.gt.imax.or.j.lt.imin.or.k.gt.kmax.or.k.lt.kmin.or.l.gt.kmax.or.l.lt.kmin) goto 11


      fpp(ijt,i,j,k,l)=vint

 11   enddo
      close(2)


      return
      end subroutine readfin

!*****************************************************************************
   subroutine readro(fname,ig,ron,ndgg)

      implicit double precision (a-h,o-z)

!      include 'formats_eqm.inc'
      include 'types_eqm_hole.inc'

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
      subroutine redrsum(ifonmx,ronp,nronp,ronh,nronh,ipozbr,ndbr)

      implicit double precision (a-h,o-z)

!      include 'formats_eqm.inc'
      include 'types_eqm_hole.inc'

      type(rho_typ), dimension(:), allocatable :: ronp,ropp,ronh,roph
      integer, dimension (:,:,:), allocatable :: ipozbr
      integer, dimension (:,:), allocatable :: ndbr

      character(len=30)fname

      ndipo=50000

      if (.not.allocated(ipozbr)) allocate(ipozbr(2,ifonmx,ndipo))
      ipozbr=0

      if (.not.allocated(ndbr)) allocate(ndbr(2,ifonmx))
      ndbr=0



      do i=1,nronp
       ibt=ronp(i)%ilap
       ndbr(1,ibt)=ndbr(1,ibt)+1

        if (ndbr(1,ibt).gt.ndipo) then
               write(*,*)' Increase dimension in redrsum'
               stop
            endif
       ipozbr(1,ibt,ndbr(1,ibt))=i

      enddo


      do i=1,nronh
       ibt=ronh(i)%ilap
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
      subroutine readro11(fname,ig,ron,ndla,n1mn,n1mx,n2mn,n2mx,nsi)

      implicit double precision (a-h,o-z)

      include 'types_eqm_hole.inc'

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
      if (.not.allocated(ron)) allocate(ron(ndla,0:nsi,n1mn:n1mx,n2mn:n2mx))

      ron=0.d0

      do ii=1,ndgg

        ibt=ronn(ii)%ilap
        isit=ronn(ii)%j
        i1t=ronn(ii)%i1
        i2t=ronn(ii)%i2
        rot=ronn(ii)%ro

        if (ibt.gt.ndla.or.isit.gt.nsi.or.i1t.lt.n1mn.or.i1t.gt.n1mx.or.i2t.lt.n2mn.or.i2t.gt.n2mx) then
         write(*,*)' Small dimensions of ron in read11',ibt,isit,i1t,i2t
         stop
       endif

        ron(ibt,isit,i1t,i2t)=rot

       enddo

      return
      end subroutine readro11


!******************************************************************************

end module read_densfmat

