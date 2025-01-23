!*     Phoninterac contains routines for calculation of redefined phonon
!*     interaction in p-h J-scheme
!*     last update 10.2.2015

      module phoninterac

      contains
      
      subroutine vinth(ifmx,isimax,lev,nlev,xthrun_min,xthrun_max,ih_lev,ip_lev,par_lev,hol_lev)

      use anglib   ! angular momentum staff
      use read_densfmat

      implicit double precision(a-h,o-z)

      include 'types_phon_int.inc'
      include 'input_phon_int.inc'
      include 'formats_phon_int.inc'


      double precision, dimension(:,:,:,:,:),allocatable ::fp

     double precision, dimension(:,:,:),allocatable :: vint
     integer, dimension(:), allocatable :: par_lev, hol_lev


      integer, dimension(:), allocatable :: jphon,irop,iroh,ndcam,ius

      type(amp_typ), dimension(:,:), allocatable :: camp,camn
      type(rho_typ), dimension(:), allocatable :: rop,roh
      type(level_typ),dimension(*) :: lev
      type(ro_typ),dimension(:), allocatable :: rh

      character*10 namer,namerp,namerh,namev
      character*30 namefp,namefpn,namecph
      character*4 nlam
      
      ndrho=10000000
      allocate(rh(ndrho))

      allocate(ius(ifmx))
      ius=0

!     ius=1
     
      jmin=0
      jmax=isimax

      ndlam=10000      
      ndamp=10000
      ifmxx=10000
      allocate (jphon(ifmxx))
      jphon=0


      ius=0
      open (3,file='1phonon/1f_states.dat',status='old',form='unformatted')

      do while (.not.eof(3))
       read(3)i,ipar,ijj,en
       jphon(i)=ijj
       if (en.le.xthrun_max.and.en.ge.xthrun_min) ius(i)=1
      enddo

      close(3)

      ifmx=i

      icount=0
      do i=1,ifmx
       if (ius(i).eq.1) icount=icount+1
      enddo

      write(*,*)' Number of selected phonons  ',icount

      

     namecph='1phonon/1f_cph.dat'
     namefp='fmat.dat'
     namerp='1f_rp.dat'
     namerh='1f_rh.dat'
     namev='V_phon_hol'

      write(*,*)'Phonons for which interaction is calculated?',1,ifmx
      read(*,*)ig_min,ig_max
      write(*,*)ig_min,ig_max


!c     loads F(p) or F(n) interaction
      call readfin(namefp,jmin,jmax,1,nlev,1,nlev,fp)
!c     loads Cph protons
!      call readcam(namecph,ndlam,ndamp,camp,ndcam)
!c     loads Cph neutrons
!      call readcam(namecn,ndlam,ndamp,camn,ndcamn)
      

 write(*,*)'Calculation of phonon-hole interaction'

    allocate(vint(0:isimax,1:nlev,1:nlev))

      vint=0.0d0


!      do ig=1,ifmx
      do ig=ig_min,ig_max

      iirg=0

!      write(5)ig

      if (ius(ig).ne.0) then


       write(nlam,'(i4.4)')ig

       open(5,file='scratch/'//namev//'_'//nlam,status='unknown',form='unformatted')



       write(997,*)ig,ifmx


       jig=jphon(ig)



!       call readro(namernp,ig,ronp,nronp)
       call readro(namerp,ig,rop,nrop)
       call readro(namerh,ig,roh,nroh)
!       call readro(namerph,ig,roph,nroph)


    
      do ib=1,ifmx
       jib=jphon(ib)

!        isi_min=iabs(jib-jig)
!        isi_max=jib+jig
        isi_min=0
        isi_max=isimax

      
!      do isi=0,isimax 
!      ifaz=(-1)**(jib+isi)
       
!      allocate(ironp(nronp))
!      ironp=0 
      allocate(irop(nrop))
      irop=0 
      allocate(iroh(nroh))
      iroh=0 
!      allocate(iroph(nroph))
!      iroph=0 

      call rosub2(ib,rop,nrop,irop,nrops)
!      call rosub2(ib,ronp,nronp,ironp,nronps)     
!      call rosub2(ib,roph,nroph,iroph,nrophs)
      call rosub2(ib,roh,nroh,iroh,nrohs)     
    
     
!      do ih1=ihmin,ihmax
!       jh1=lev(ih1)%j 
!      do ih2=ihmin,ihmax
!       jh2=lev(ih2)%j
             
        vint=0.d0

        do ii=1,nrops
         i1=rop(irop(ii))%i1
         i2=rop(irop(ii))%i2
         isi=rop(irop(ii))%j
         ifaz=(-1)**(jib+isi)
        do ihh1=1,ih_lev
          ih1=hol_lev(ihh1)

          do ihh2=1,ih_lev
            ih2=hol_lev(ihh2)
    
           vint(isi,ih1,ih2)=vint(isi,ih1,ih2)+dfloat(ifaz)*rop(irop(ii))%ro*fp(isi,i1,i2,ih1,ih2)
          enddo
        enddo
        enddo
 

        do ii=1,nrohs
         i1=roh(iroh(ii))%i1
         i2=roh(iroh(ii))%i2
         isi=roh(iroh(ii))%j
         ifaz=(-1)**(jib+isi)
        do ihh1=1,ih_lev
          ih1=hol_lev(ihh1)

          do ihh2=1,ih_lev
            ih2=hol_lev(ihh2)
     
         vint(isi,ih1,ih2)=vint(isi,ih1,ih2)+0.5d0*dfloat(ifaz)*roh(iroh(ii))%ro*fp(isi,i1,i2,ih1,ih2)
          enddo
        enddo
        enddo


!c        if (dabs(vint).gt.xrotrunc.and.ipart.eq.1) 
!c     *write(997,102)isi,ig,ib,ih1,ih2,vint

!c        if (dabs(vint).gt.xrotrunc.and.ipart.eq.-1) 
!c     *write(998,102)isi,ig,ib,ih1,ih2,vint


!        if (dabs(vint).gt.xrotrunc) 
!     *write(5)ib,isi,ih1,ih2,vint

!       do isi=0,isimax
       do isi=isi_min,isi_max

       do ihh1=1,ih_lev
        ih1=hol_lev(ihh1)
       do ihh2=1,ih_lev
        ih2=hol_lev(ihh2)

       if (dabs(vint(isi,ih1,ih2)).gt.xrotrunc) then
               iirg=iirg+1
               if (iirg.gt.ndrho) then
                  write(*,*)' Increase dimension of ndrho in rop!!'
                  stop
               endif
               rh(iirg)%ib=ib
               rh(iirg)%isi=isi
               rh(iirg)%i1=ih1
               rh(iirg)%i2=ih2
               rh(iirg)%rho=vint(isi,ih1,ih2)
!c        write(33)ib,isi,ip1,ip2,ronp
       endif
       enddo
       enddo
       enddo




!        enddo ! loop ih2
!      enddo ! loop ih1
      deallocate (irop,iroh)
!      enddo ! loop isi
      enddo ! loop ib

!      endif

      write(5)ig,iirg
      if (iirg.gt.0) then
      write(5)(rh(iii)%ib,iii=1,iirg)
      write(5)(rh(iii)%isi,iii=1,iirg)
      write(5)(rh(iii)%i1,iii=1,iirg)
      write(5)(rh(iii)%i2,iii=1,iirg)
      write(5)(rh(iii)%rho,iii=1,iirg)
      endif

!     do iii=1,iirg
!       write(332,'(5i5,f15.10)')ig,rh(iii)%ib,rh(iii)%isi,rh(iii)%i1,rh(iii)%i2,rh(iii)%rho
!      enddo
!      write(332,*)ig,iirg


      close(5)

      endif

!      if (iirg.eq.0) then
!      write(5)iirg
!      write(5)iirg
!      write(5)iirg
 !     write(5)iirg
 !     write(5)dfloat(iirg)
 !     endif



!      write(5)0,0,0,0,0.d0      
      enddo ! loop ig      

      deallocate(rh,jphon,fp,rop,roh)

!      write(5)10000000
!      write(5)0,0,0,0,0.d0 

!      close(33)
!      close(34)
!      close(43)
!      close(44)
 !     close(5)
      return
      end subroutine vinth
!*************************************************************************
      subroutine vintp(ifmx,isimax,lev,nlev,xthrun_min,xthrun_max,ih_lev,ip_lev,par_lev,hol_lev)

      use anglib   ! angular momentum staff
      use read_densfmat

      implicit double precision(a-h,o-z)

      include 'types_phon_int.inc'
      include 'input_phon_int.inc'
      include 'formats_phon_int.inc'

      double precision, dimension(:,:,:,:,:),allocatable ::fp

      double precision, dimension(:,:,:),allocatable :: vint
      integer, dimension(:), allocatable :: par_lev, hol_lev

      integer, dimension(:), allocatable :: jphon,irop,iroh,ndcam,ius

!      type(amp_typ), dimension(:,:), allocatable :: camp,camn
      type(rho_typ), dimension(:), allocatable :: rop,roh
      type(level_typ),dimension(*) :: lev
      type(ro_typ),dimension(:), allocatable :: rh


      character*10 namer,namerp,namerh
      character*30 namefp,namefpn,namecph
      character*10 namev
      character*4 nlam

      ndrho=100000000
      allocate(rh(ndrho))


      allocate(ius(ifmx))
      ius=0
      
!     ius=1

     
      jmin=0
      jmax=isimax

      ndlam=10000      
      ndamp=10000
      ifmxx=10000
      allocate (jphon(ifmxx))
      jphon=0

      ius=0

      open (3,file='1phonon/1f_states.dat',status='old',form='unformatted')

      do while (.not.eof(3))
       read(3)i,ipar,ijj,en
!c       write(991,*)i,ipar,ijj,en
       jphon(i)=ijj
!       if (en.le.xthrun) ius(i)=1
       if (en.le.xthrun_max.and.en.ge.xthrun_min) ius(i)=1

      enddo

      close(3)

      ifmx=i

      icount=0
      do i=1,ifmx
       if (ius(i).eq.1) icount=icount+1
      enddo

      write(*,*)' Number of selected phonons ',icount

      namecph='1phonon/1f_cph.dat'
      namefp='fmat.dat'
      namerp='1f_rp.dat'
      namerh='1f_rh.dat'
      namev='V_phon_par'


      write(*,*)'Phonons for which interaction is calculated?',1,ifmx
      read(*,*)ig_min,ig_max
      write(*,*)ig_min,ig_max


!     loads F(p) or F(n) interaction
      call readfin(namefp,jmin,jmax,1,nlev,1,nlev,fp)
!     loads F(pn) 
!      call readfin(namefpn,jmin,jmax,ihminp,ipmaxp,ihmin,ipmax,fpn)
!     loads Cph protons
!      call readcam(namecp,ndlam,ndamp,camp,ndcamp)
!     loads Cph neutrons
!      call readcam(namecn,ndlam,ndamp,camn,ndcamn)
      

      write(*,*)'Calculation of phonon-particle interaction'


      allocate(vint(0:isimax,1:nlev,1:nlev))
      vint=0.0d0


!      do ig=1,ifmx
       do ig=ig_min,ig_max 
   
       iirg=0 

       jig=jphon(ig)

!      write(5)ig
      if (ius(ig).ne.0) then

!       open(5,file=namev,status='unknown'
!     *,form='unformatted')


       write(nlam,'(i4.4)')ig
       
       open(5,file='scratch/'//namev//'_'//nlam,status='unknown',form='unformatted')

       write(996,*)ig,ifmx

      call readro(namerp,ig,rop,nrop)
      call readro(namerh,ig,roh,nroh)

    
      do ib=1,ifmx
       jib=jphon(ib)
       isi_min=0 !iabs(jib-jig)
       isi_max=isimax !jib+jig

      
!      do isi=0,isimax 
      ifaz=(-1)**(jib+isi)
       
      allocate(irop(nrop))
      irop=0 
      allocate(iroh(nroh))
      iroh=0 

      call rosub2(ib,rop,nrop,irop,nrops)
      call rosub2(ib,roh,nroh,iroh,nrohs)     
    
     
!      do ip1=ipmin,ipmax
!       jp1=lev(ip1)%j 
!      do ip2=ipmin,ipmax
!       jp2=lev(ip2)%j
             
        vint=0.d0

        do ii=1,nrops
         i1=rop(irop(ii))%i1
         i2=rop(irop(ii))%i2
         isi=rop(irop(ii))%j
         ifaz=(-1)**(jib+isi)

        do ipp1=1,ip_lev
         ip1=par_lev(ipp1)
         do ipp2=1,ip_lev
         ip2=par_lev(ipp2)
         vint(isi,ip1,ip2)=vint(isi,ip1,ip2)+0.5d0*dfloat(ifaz)*rop(irop(ii))%ro*fp(isi,i1,i2,ip1,ip2)
         enddo
        enddo

        enddo


        do ii=1,nrohs
         i1=roh(iroh(ii))%i1
         i2=roh(iroh(ii))%i2
         isi=roh(iroh(ii))%j
         ifaz=(-1)**(jib+isi)

        do ipp1=1,ip_lev
         ip1=par_lev(ipp1)
         do ipp2=1,ip_lev
          ip2=par_lev(ipp2)
          vint(isi,ip1,ip2)=vint(isi,ip1,ip2)+dfloat(ifaz)*roh(iroh(ii))%ro*fp(isi,i1,i2,ip1,ip2)
         enddo
        enddo
        enddo



       do isi=isi_min,isi_max

        do ipp1=1,ip_lev
         ip1=par_lev(ipp1)
         do ipp2=1,ip_lev
          ip2=par_lev(ipp2)

       if (dabs(vint(isi,ip1,ip2)).gt.xrotrunc) then
               iirg=iirg+1
               if (iirg.gt.ndrho) then
                  write(*,*)' Increase dimension of ndrho in rop!!'
                  stop
               endif
               rh(iirg)%ib=ib
               rh(iirg)%isi=isi
               rh(iirg)%i1=ip1
               rh(iirg)%i2=ip2
               rh(iirg)%rho=vint(isi,ip1,ip2)
!c        write(33)ib,isi,ip1,ip2,ronp
       endif

         enddo
        enddo

       enddo




!        enddo ! loop ih2
!      enddo ! loop ih1
      deallocate (irop,iroh)
!      enddo ! loop isi
      enddo ! loop ib
  
!      endif

      write(5)ig,iirg
      if (iirg.gt.0) then
      write(5)(rh(iii)%ib,iii=1,iirg)
      write(5)(rh(iii)%isi,iii=1,iirg)
      write(5)(rh(iii)%i1,iii=1,iirg)
      write(5)(rh(iii)%i2,iii=1,iirg)
      write(5)(rh(iii)%rho,iii=1,iirg)
      endif

!      do iii=1,iirg
!       write(331,'(5i5,f15.10)')ig,rh(iii)%ib,rh(iii)%isi,rh(iii)%i1,rh(iii)%i2,rh(iii)%rho
!      enddo
!      write(331,*)ig,iirg


      close(5)

      endif

!      if (iirg.eq.0) then
!      write(5)iirg
!      write(5)iirg
!      write(5)iirg
!      write(5)iirg
!      write(5)dfloat(iirg)
!      endif



!      write(5)0,0,0,0,0.d0 

      enddo ! loop ig      

      deallocate(rh,jphon,fp,rop,roh)

!      write(5)10000000
!      write(5)0,0,0,0,0.d0 

!      close(33)
!      close(34)
!      close(43)
!      close(44)
!      close(5)
      return
      end subroutine vintp


!************************************************************************
      subroutine readcam(fname,ndimi,ndimj,cam,ndcc)

      implicit double precision (a-h,o-z)

      include 'formats_phon_int.inc'
      include 'types_phon_int.inc'

      type (amp_typ), dimension(:,:), allocatable :: cam
      integer, dimension (:), allocatable :: ndcc

      character(len=30)fname

      allocate(cam(ndimi,ndimj))
      allocate(ndcc(ndimi))
  
      open(2,file=fname,status='old',form='unformatted')

      ilam=0

      do while (.not.eof(2))
       ilam=ilam+1
       if (ilam.gt.ndimi) then 
               write(*,*)'Readcam: allocate bigger array in ndimi'
               stop
           endif
       read(2)ipar,ijj,ndc
         if (ndc.gt.ndimj) then 
            write(*,*)'Readcam: allocate bigger array in ndimj'
               stop
           endif
  
       read(2)(cam(ilam,i)%par,cam(ilam,i)%hol,cam(ilam,i)%am,i=1,ndc)
       ndcc(ilam)=ndc
      enddo

      close(2)

      return
      end subroutine readcam

!
!*****************************************************************************
      subroutine readro(fname,ig,ron,ndgg)

      implicit double precision (a-h,o-z)

!      include 'formats_eqm.inc'
      include 'types_phon_int.inc'

      type(rho_typ), dimension(:), allocatable :: ron

      character(len=10)fname
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




      subroutine rosub(j,ilam,rop,nrop,irop,nrops)

      implicit double precision (a-h,o-z)

      include 'types_phon_int.inc'

      type(rho_typ), dimension(:), allocatable :: rop
      integer, dimension(:), allocatable :: irop

      ii=0
      do i=1,nrop
      if (rop(i)%j.eq.j.and.rop(i)%ilap.eq.ilam) then
              ii=ii+1
              irop(ii)=i
             endif
 
      enddo

      nrops=ii

      end subroutine rosub
!***********************************************************************
      subroutine rosub2(ilam,rop,nrop,irop,nrops)

      implicit double precision (a-h,o-z)

      include 'types_phon_int.inc'

      type(rho_typ), dimension(:), allocatable :: rop
      integer, dimension(:), allocatable :: irop

      ii=0
      do i=1,nrop
      if (rop(i)%ilap.eq.ilam) then
              ii=ii+1
              irop(ii)=i
             endif

      enddo

      nrops=ii

      end subroutine rosub2

     
      end module phoninterac 
      
      
