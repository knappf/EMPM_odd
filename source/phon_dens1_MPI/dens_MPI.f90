module dens_pn
      

contains 
!     subrouutines and functions       

!************************************************************************
     subroutine rop(ifmx,nlev,isimax,lev,xthres,ih_lev,ip_lev,par_lev,hol_lev,myid,numprocs)

      use anglib   ! angular momentum staff

      implicit double precision(a-h,o-z)

      include 'types_phon_dens1.inc'
      include 'input_phon_dens1.inc'

      double precision, dimension(:,:,:), allocatable :: came

      integer, dimension(:), allocatable :: jphon,ig_resh
      integer, dimension(:), allocatable :: par_lev, hol_lev

      type(amp_typ), dimension(:), allocatable :: camp
      type(level_typ),dimension(*) :: lev
      type(ro_typ),dimension(:), allocatable :: rh


      character*9 namer
      character*10 namec
      character*4 nlam

      xrotrunc=0.0d-12
      ndamp=10000
      ndrho=10000000     

      allocate (camp(ndamp))
      allocate (jphon(ifmx))
      jphon=0

      allocate(rh(ndrho))

      namer='1f_rp.dat'
      namec='1f_cph.dat'

      allocate (came(ifmx,nlev,nlev))
      came=0.d0

      if (myid == 0 ) then 
        write(*,*)'Calculation of 1-phonon particle  densities'
        write(*,*)' # of MPI processes = ',numprocs
      endif 


      open(11,file='1phonon/'//namec,status='old',form='unformatted')

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

      close(11)

      allocate(ig_resh(ifmx))
      
      do i=1,ifmx
        ig_resh(i)=i
      enddo
      icount=ifmx

      call random_permutation(ig_resh,icount)
      
      if (mod(icount,numprocs).eq.0) then
              n_seg=icount/numprocs
      else
              n_seg=icount/numprocs+1
      endif

      if (myid.eq.0) write(*,*) ' size of segment = ',n_seg  

      do igg=myid*n_seg+1,min((myid+1)*n_seg,icount)

      ig=ig_resh(igg)

      write(nlam,'(i4.4)')ig

      open(33,file='scratch/'//namer//'_'//nlam,status='unknown',form='unformatted')

       iirg=0

!c       write(33)ig
!c       write(997,*)ig
       jig=jphon(ig)


      do ib=1,ifmx
       jib=jphon(ib)

       isi_min=iabs(jib-jig)
       isi_max=jib+jig


      do ipp1=1,ip_lev
       ip1=par_lev(ipp1)
       jp1=lev(ip1)%j
      do ipp2=1,ip_lev
       ip2=par_lev(ipp2)
       jp2=lev(ip2)%j

!         if (((jp1+jp2)/2).lt.isi_max) isi_max=(jp1+jp2)/2
!         if ((iabs(jp1-jp2)/2).gt.isi_min) isi_min=iabs(jp1-jp2)/2

         do isi=isi_min,isi_max

         ronp=0.d0

          do ihh=1,ih_lev
            ih=hol_lev(ihh)
            jh=lev(ih)%j
            ifaz=(-1)**(jig+isi+(jh+jp1)/2)
            xsixj=sixj(2*jib,2*isi,2*jig,jp2,jh,jp1)
            fact=(dfloat((2*jib+1)*(2*jig+1)*(2*isi+1)))**0.5d0
            ronp=ronp+dfloat(ifaz)*fact*came(ig,ip2,ih)*came(ib,ip1,ih)*xsixj
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
       endif

        enddo ! loop isi

      enddo ! loop ip2
      enddo ! loop ip1


      enddo ! loop ib


      write(33)ig,iirg
      if (iirg.gt.0) then
!       write(33)(rh(iii),iii=1,iirg)
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

      close(33)
      enddo ! loop ig      
      

      deallocate(came,camp,jphon,rh)
!c      write(33)10000000
!c      write(33)0,0,0,0,0.d0
      close(33)
      return
      end subroutine rop
!************************************************************************

      subroutine roh(ifmx,nlev,isimax,lev,xthres,ih_lev,ip_lev,par_lev,hol_lev,myid,numprocs)

      use anglib   ! angular momentum staff

      implicit double precision(a-h,o-z)

      include 'types_phon_dens1.inc'
      include 'input_phon_dens1.inc'

      double precision, dimension(:,:,:), allocatable :: came

      integer, dimension(:), allocatable :: jphon,ig_resh
      integer, dimension(:), allocatable :: par_lev, hol_lev
      type(amp_typ), dimension(:), allocatable :: camp
      type(level_typ),dimension(*) :: lev
      type(ro_typ),dimension(:), allocatable :: rh


      character*9 namer
      character*11 namerf
      character*10 namec
      character*4 nlam

      xrotrunc=0.0d-14
      ndamp=10000    
      ndrho=10000000  
      allocate (camp(ndamp))
      allocate (jphon(ifmx))
      jphon=0

      allocate(rh(ndrho))

      namer='1f_rh.dat'
      namec='1f_cph.dat'

      allocate (came(ifmx,nlev,nlev))
      came=0.d0

     write(*,*)'Calculation of 1-phonon hole densities'

      open(11,file='1phonon/'//namec,status='old',form='unformatted')

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
      
      close(11)

      allocate(ig_resh(ifmx))
      
      do i=1,ifmx
        ig_resh(i)=i
      enddo
      icount=ifmx

      call random_permutation(ig_resh,icount)
      
      if (mod(icount,numprocs).eq.0) then
              n_seg=icount/numprocs
      else
              n_seg=icount/numprocs+1
      endif

      if (myid.eq.0) write(*,*) ' size of segment = ',n_seg  

      do igg=myid*n_seg+1,min((myid+1)*n_seg,icount)

      ig=ig_resh(igg)
      
!      do ig=1,ifmx

      write(nlam,'(i4.4)')ig

      open(33,file='scratch/'//namer//'_'//nlam,status='unknown',form='unformatted')


       iirg=0

       jig=jphon(ig)


      do ib=1,ifmx
       jib=jphon(ib)
       isi_min=iabs(jib-jig)
       isi_max=jib+jig


!      do isi=0,isimax
      do ihh1=1,ih_lev
       ih1=hol_lev(ihh1)
       jh1=lev(ih1)%j 
      do ihh2=1,ih_lev
       ih2=hol_lev(ihh2)
       jh2=lev(ih2)%j
 
       do isi=isi_min,isi_max

          ronh=0.d0
      
          do ipp=1,ip_lev
            ip=par_lev(ipp)
            jp=lev(ip)%j
            ifaz=(-1)**(jib+(jp+jh2)/2)
            xsixj=sixj(2*jib,2*isi,2*jig,jh1,jp,jh2) 
            fact=(dfloat((2*jib+1)*(2*jig+1)*(2*isi+1)))**0.5d0
            ronh=ronh-dfloat(ifaz)*fact*came(ig,ip,ih1)*came(ib,ip,ih2)*xsixj
          enddo
          
!       if (dabs(ronh).gt.xrotrunc) write(33)ib,isi,ih1,ih2,ronh

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

           endif
        enddo ! loop isi
      enddo ! loop ih2
      enddo ! loop ih1
      enddo ! loop ib

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

      close(33)
      enddo ! loop ig      
      
      deallocate(came,camp,jphon,rh)
!c      write(33)10000000
!c      write(33)0,0,0,0,0.d0
!      close(33)
!c      write(433)10000000
!c      write(433)0,0,0,0,0.d0
!c      close(433)


      return
      end subroutine roh
!************************************************************************

subroutine random_permutation(arr, n)
      USE IFPORT
      implicit none
      integer, intent(inout) :: arr(:), n
      integer :: i, j, temp

      ! Initialize the random number generator
      call random_seed()

      ! Generate a random permutation of the array
      do i = n, 2, -1
         j = 1 + int(i*rand())
         temp = arr(i)
         arr(i) = arr(j)
         arr(j) = temp
      end do

end subroutine random_permutation

end module dens_pn      
      
      
