!     
!     Program phon_dens1 computes 1phonon densities in proton-neutron 
!     J-coupled formalism.

!     last update 12.4.2013

      program phon_dens1
 
      use input_sp


      implicit double precision (a-h,o-z)

      include 'input_phon_dens1.inc'
      include 'formats_phon_dens1.inc'
      include 'types_phon_dens1.inc'

      type(level_typ),dimension(:), allocatable :: levn,levp
      
*     loading of input data 

      xrotrunc=1.d-8
      
      write(*,*)'Loading of input '

      open(1,file='input_tda_coup.dat',status='old',form='formatted')
            
      read(1,15)ia,iz
      read(1,15)ihnmn,ihnmx
      read(1,15)ihpmn,ihpmx
      read(1,15)ipnmn,ipnmx
      read(1,15)ippmn,ippmx
      read(1,26)alfa,beta
      read(1,*)
      read(1,15)iparmn,iparmx
      read(1,15)jminn,jmaxn
      close(1)

      allocate(levn(ipnmx+1000),levp(ippmx+1000))

      call inp_sp(levn,levp)


      open (3,file='1phonon/1f_states.dat',
     *status='old',form='unformatted')

      do while (.not.eof(3))
       read(3)i,ipar,ijj,en
c       jphon(i)=ijj
      enddo

      close(3)



      ifmx=i ! treba upravit aby sa to nacitalo 
      jmax=jmaxn

      write(*,*)' Energy threshodl for 1 phonon states?'
      read(*,*)xthres
      write(*,*)xthres
      call rop(-1,ifmx,jmax,levp,xthres)
      call rop(1,ifmx,jmax,levn,xthres) 
      call roh(-1,ifmx,jmax,levp,xthres)
      call roh(1,ifmx,jmax,levn,xthres) 

      
      end 
!     
!     END of the main program 
! 

      
!     subrouutines and functions       
!************************************************************************
      subroutine rop(ipart,ifmx,isimax,lev,xthres)

      use anglib   ! angular momentum staff

      implicit double precision(a-h,o-z)

      include 'types_phon_dens1.inc'
      include 'input_phon_dens1.inc'

      double precision, dimension(:,:,:), allocatable :: 
     *came

      integer, dimension(:), allocatable :: jphon,ius

      type(amp_typ), dimension(:), allocatable :: camp
      type(level_typ),dimension(*) :: lev
      type(ro_typ),dimension(:), allocatable :: rh


      character*10 namer
      character*9  namec

      ndamp=10000
      ndrho=10000000     

      allocate (camp(ndamp))
      allocate (jphon(ifmx))
      jphon=0

      allocate(rh(ndrho))

      if (ipart.eq.1) then 
              namer='1f_rnp.dat'
              namec='1f_cn.dat'
              ipmin=ipnmn
              ipmax=ipnmx
              ihmin=ihnmn
              ihmax=ihnmx
      endif

      if (ipart.eq.-1) then 
              namer='1f_rpp.dat'
              namec='1f_cp.dat'
              ipmin=ippmn
              ipmax=ippmx
              ihmin=ihpmn
              ihmax=ihpmx
      endif

      allocate (came(ifmx,ipmin:ipmax,ihmin:ihmax))
      came=0.d0

      allocate(ius(ifmx))
      ius=0



      open (3,file='0hom_phon.dat',
     *status='old',form='formatted')

      do while (.not.eof(3))
       read(3,*)i,ee
       ius(i)=1
      enddo

      close(3)


      open (3,file='1hom_phon.dat',
     *status='old',form='formatted')

      do while (.not.eof(3))
       read(3,*)i,ee
       ius(i)=1
      enddo

      close(3)

      open (3,file='2hom_phon.dat',
     *status='old',form='formatted')

      do while (.not.eof(3))
       read(3,*)i,ee
       if (ee.le.xthres) ius(i)=1
      enddo

      close(3)


      open (3,file='3hom_phon.dat',
     *status='old',form='formatted')

      do while (.not.eof(3))
       read(3,*)i,ee
       if (ee.le.xthres) ius(i)=1

!       ius(i)=1
      enddo

      close(3)

      ius=1

      icount=0
      do i=1,ifmx
       if (ius(i).eq.1) icount=icount+1
      enddo

      write(*,*)' Number of 0hom+1hom+2hom+3hom phonons ',icount


 
      if (ipart.eq.1)
     *write(*,*)'Calculation of 1-phonon particle neutron densities'

      if (ipart.eq.-1)
     *write(*,*)'Calculation of 1-phonon particle proton densities'


      open(33,file='1phonon/'//namer,status='unknown'
     *,form='unformatted')

      open(11,file='1phonon/'//namec,status='old',
     *form='unformatted')

      ilam=0
      do while (.not.eof(11))
       ilam=ilam+1
       read(11)ipar,ijj,ndc
       jphon(ilam)=ijj
        if (ndc.gt.ndamp) then
                           write(*,*)'small dimension of array camp'
                           stop
                   endif 
  
       read(11)(camp(i)%par,camp(i)%hol,camp(i)%am,i=1,ndc)

        do i=1,ndc
         ipar=camp(i)%par
         ihol=camp(i)%hol
         cc=camp(i)%am
         came(ilam,ipar,ihol)=cc
        enddo

      enddo

      

      do ig=1,ifmx

       iirg=0

c       write(33)ig
c       write(997,*)ig
       jig=jphon(ig)

      if (ius(ig).ne.0) then 

      do ib=1,ifmx
       jib=jphon(ib)

      do isi=0,isimax
 
      do ip1=ipmin,ipmax
       jp1=lev(ip1)%j 
      do ip2=ipmin,ipmax
       jp2=lev(ip2)%j
   
          ronp=0.d0
      
          do ih=ihmin,ihmax
            jh=lev(ih)%j
            ifaz=(-1)**(jig+isi+(jh+jp1)/2)
            xsixj=sixj(2*jib,2*isi,2*jig,jp2,jh,jp1) 
            fact=(dfloat((2*jib+1)*(2*jig+1)*(2*isi+1)))**0.5d0
            ronp=ronp+dfloat(ifaz)*fact*
     *came(ig,ip2,ih)*came(ib,ip1,ih)*xsixj
          enddo
             
       if (dabs(ronp).gt.xrotrunc) then
               iirg=iirg+1
               if (iirg.gt.ndrho) then
                  write(*,*)' Increase dimension of ndrho in rop!!'
                  stop
               endif
               rh(iirg)%ib=ib
               rh(iirg)%isi=isi
               rh(iirg)%i1=ip1
               rh(iirg)%i2=ip2
               rh(iirg)%rho=ronp
c        write(33)ib,isi,ip1,ip2,ronp
       endif

c       if (dabs(ronp).gt.xrotrunc.and.ipart.eq.-1) 
c     *write(997,*)ib,isi,ip1,ip2,ronp

c       if (dabs(ronp).gt.xrotrunc.and.ipart.eq.1) 
c     *write(977,*)ib,isi,ip1,ip2,ronp


  55    continue      
        enddo ! loop ip2
      enddo ! loop ip1
      enddo ! loop isi
      enddo ! loop ib

      endif

      write(33)ig,iirg
      if (iirg.gt.0) then
      write(33)(rh(iii)%ib,iii=1,iirg)
      write(33)(rh(iii)%isi,iii=1,iirg)
      write(33)(rh(iii)%i1,iii=1,iirg)
      write(33)(rh(iii)%i2,iii=1,iirg)
      write(33)(rh(iii)%rho,iii=1,iirg)
      endif

      if (iirg.eq.0) then       
      write(33)iirg
      write(33)iirg
      write(33)iirg
      write(33)iirg
      write(33)dfloat(iirg)
      endif

      write(997,*)ig,iirg

!c      write(977,*)0,0,0,0,0.d0
!      endif
      enddo ! loop ig      
      

      deallocate(came,camp,jphon,rh)
!c      write(33)10000000
!c      write(33)0,0,0,0,0.d0
      close(33)
      return
      end 
!************************************************************************

      subroutine roh(ipart,ifmx,isimax,lev,xthres)

      use anglib   ! angular momentum staff

      implicit double precision(a-h,o-z)

      include 'types_phon_dens1.inc'
      include 'input_phon_dens1.inc'

      double precision, dimension(:,:,:), allocatable :: 
     *came

      integer, dimension(:), allocatable :: jphon,ius

      type(amp_typ), dimension(:), allocatable :: camp
      type(level_typ),dimension(*) :: lev
      type(ro_typ),dimension(:), allocatable :: rh


      character*10 namer
      character*11 namerf
      character*9  namec

      ndamp=10000    
      ndrho=10000000  
      allocate (camp(ndamp))
      allocate (jphon(ifmx))
      jphon=0

      allocate(rh(ndrho))

      if (ipart.eq.1) then 
              namer='1f_rnp.dat'
              namec='1f_cn.dat'
              ipmin=ipnmn
              ipmax=ipnmx
              ihmin=ihnmn
              ihmax=ihnmx
      endif



      if (ipart.eq.1) then 
              namer='1f_rnh.dat'
c              namerf='1f_rnhf.dat'
              namec='1f_cn.dat'
              ipmin=ipnmn
              ipmax=ipnmx
              ihmin=ihnmn
              ihmax=ihnmx
      endif

      if (ipart.eq.-1) then 
              namer='1f_rph.dat'
c              namerf='1f_rphf.dat'
              namec='1f_cp.dat'
              ipmin=ippmn
              ipmax=ippmx
              ihmin=ihpmn
              ihmax=ihpmx
      endif

      allocate (came(ifmx,ipmin:ipmax,ihmin:ihmax))
      came=0.d0

      allocate(ius(ifmx))
      ius=0


      open (3,file='0hom_phon.dat',
     *status='old',form='formatted')

      do while (.not.eof(3))
       read(3,*)i,ee
       ius(i)=1
      enddo

      close(3)


      open (3,file='1hom_phon.dat',
     *status='old',form='formatted')

      do while (.not.eof(3))
       read(3,*)i,ee
       ius(i)=1
      enddo

      close(3)

      open (3,file='2hom_phon.dat',
     *status='old',form='formatted')

      do while (.not.eof(3))
       read(3,*)i,ee
       if (ee.le.xthres) ius(i)=1
!       ius(i)=1
      enddo

      close(3)

      open (3,file='3hom_phon.dat',
     *status='old',form='formatted')

      do while (.not.eof(3))
       read(3,*)i,ee
       if (ee.le.xthres) ius(i)=1

!       ius(i)=1
      enddo

      close(3)


      ius=1

      icount=0
      do i=1,ifmx
       if (ius(i).eq.1) icount=icount+1
      enddo





      write(*,*)' Number of 0hom+1hom+2hom+3hom phonons ',icount


 
      if (ipart.eq.1)
     *write(*,*)'Calculation of 1-phonon hole neutron densities'

      if (ipart.eq.-1)
     *write(*,*)'Calculation of 1-phonon hole proton densities'


      open(33,file='1phonon/'//namer,status='unknown'
     *,form='unformatted')

c      open(433,file='1phonon/'//namerf,status='unknown'
c     *,form='unformatted')


      open(11,file='1phonon/'//namec,status='old',
     *form='unformatted')

      ilam=0
      do while (.not.eof(11))
       ilam=ilam+1
       read(11)ipar,ijj,ndc
       jphon(ilam)=ijj
        if (ndc.gt.ndamp) then
                           write(*,*)'small dimension of array camp'
                           stop
                   endif 
  
       read(11)(camp(i)%par,camp(i)%hol,camp(i)%am,i=1,ndc)

        do i=1,ndc
         ipar=camp(i)%par
         ihol=camp(i)%hol
         cc=camp(i)%am
         came(ilam,ipar,ihol)=cc
        enddo

      enddo

      
      do ig=1,ifmx

       iirg=0

c       write(33)ig
c       write(433)ig
c       write(998,*)ig
       jig=jphon(ig)

       if (ius(ig).ne.0) then


      do ib=1,ifmx
       jib=jphon(ib)
      do isi=0,isimax
      do ih1=ihmin,ihmax
       jh1=lev(ih1)%j 
      do ih2=ihmin,ihmax
       jh2=lev(ih2)%j
   
          ronh=0.d0
      
          do ip=ipmin,ipmax
            jp=lev(ip)%j
            ifaz=(-1)**(jib+(jp+jh2)/2)
            xsixj=sixj(2*jib,2*isi,2*jig,jh1,jp,jh2) 
            fact=(dfloat((2*jib+1)*(2*jig+1)*(2*isi+1)))**0.5d0
            ronh=ronh-dfloat(ifaz)*fact*
     *came(ig,ip,ih1)*came(ib,ip,ih2)*xsixj
          enddo
          
c       if (dabs(ronh).gt.xrotrunc) write(33)ib,isi,ih1,ih2,ronh

       if (dabs(ronh).gt.xrotrunc) then
               iirg=iirg+1
               if (iirg.gt.ndrho) then
                  write(*,*)' Increase dimension of ndrho in roh!!'
                  stop
               endif
               rh(iirg)%ib=ib
               rh(iirg)%isi=isi
               rh(iirg)%i1=ih1
               rh(iirg)%i2=ih2
               rh(iirg)%rho=ronh
c        write(33)ib,isi,ip1,ip2,ronp
       endif


c       if (ib.eq.ig.and.ih1.eq.ih2.and.isi.eq.0) 
c     *ronh=ronh+dfloat((jh1+1)*(2*jig+1))**0.5d0
c       if (dabs(ronh).gt.xrotrunc) write(433)ib,isi,ih1,ih2,ronh
c       if (dabs(ronh).gt.xrotrunc) 
c     *write(998,*)ib,isi,ih1,ih2,ronh

  55    continue      
        enddo ! loop ih2
      enddo ! loop ih1
      enddo ! loop isi
      enddo ! loop ib

      endif

      write(33)ig,iirg
      if (iirg.gt.0) then
      write(33)(rh(iii)%ib,iii=1,iirg)
      write(33)(rh(iii)%isi,iii=1,iirg)
      write(33)(rh(iii)%i1,iii=1,iirg)
      write(33)(rh(iii)%i2,iii=1,iirg)
      write(33)(rh(iii)%rho,iii=1,iirg)
      endif

      if (iirg.eq.0) then
      write(33)iirg
      write(33)iirg
      write(33)iirg
      write(33)iirg
      write(33)dfloat(iirg)
      endif

!      endif
      write(998,*)ig,iirg

!c      write(33)0,0,0,0,0.d0
!c      write(433)0,0,0,0,0.d0
!c      write(998,*)0,0,0,0,0.d0
      enddo ! loop ig      
      
      deallocate(came,camp,jphon,rh)
!c      write(33)10000000
!c      write(33)0,0,0,0,0.d0
      close(33)
!c      write(433)10000000
!c      write(433)0,0,0,0,0.d0
!c      close(433)


      return
      end 
************************************************************************
      
      
      
