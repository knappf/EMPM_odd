!     last update 17.4 .2015
module rdens
 contains 

!************************************************************************
      subroutine r_value(lev,ip_lev,par_lev,ih_lev,hol_lev,jamax,ihmn,ipmx,isimax,iphous,iphous2,phonbs,nphon,ns1,ns2,myid,numprocs)

      use anglib   ! angular momentum staff

      implicit double precision(a-h,o-z)

      include 'types_phon_dens.inc'
      include 'formats_phon_dens.inc'

      type (amp2_typ), dimension(:), allocatable :: c_ampl!,cc
      integer(kind=8) :: ndimroc
      type (roc_typ), dimension(:), allocatable :: roc
      type (ro2_typ), dimension(:), allocatable :: r2p,r2h
!      type (ro2_typ), dimension(:), allocatable :: r2ph

      integer, dimension (:), allocatable :: ndx,ndc

      type (phon_typ), dimension (:,:), allocatable :: phonbs
      integer, dimension (:), allocatable :: nphon,iphous,iphous2,ipozl

      type(level_typ),dimension(:), allocatable :: lev
      integer, dimension(:), allocatable :: par_lev, hol_lev
      integer, dimension(:,:), allocatable :: i_ph_int


      double precision, dimension (:,:,:), allocatable :: x_ampl
      double precision, dimension (:,:,:,:,:,:), allocatable :: csixj
      double precision, dimension (:,:), allocatable :: rog

!     double precision, dimension(:,:,:,:,:), allocatable :: rnp,rpp,rnh,rph

      real, dimension(:,:,:,:,:), allocatable :: rph

!      double precision, dimension (:,:,:), allocatable :: rr

      character(len=30) fnamex,fnamec
      character(len=30) fnamer1,fnamer2,fnamer,fnameo
      character(len=9) fnamern,name_new
      character(len=10) myid_name
      integer(kind=8) iroc



      integer myid,numprocs
      integer, dimension(:), allocatable :: ig_resh

      jjmx=jamax/2
      allocate (csixj(0:isimax,0:jjmx,0:jjmx,0:jjmx,0:jamax,0:jamax))
      csixj=0.d0

!      if (myid.eq.0) write(*,*)'CG coef. for jmx =',jjmx

      do i1=0,isimax
      ifaz1=(-1)**(i1)
       do i2=0,jjmx
        do i3=0,jjmx
         do i4=0,jjmx
         ifaz4=(-1)**(i4)
          do i5=0,jamax
           do i6=0,jamax             
             csixj(i1,i2,i3,i4,i5,i6)=ifaz1*ifaz4*sixj(2*i1,2*i2,2*i3,2*i4,2*i5,2*i6)
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo

      write(myid_name,'(i10.10)')myid
      open(63,file='2_phon_dens_calc_myid_'//myid_name,status='unknown',form='formatted')

      ndimroc=600000
      allocate(roc(ndimroc))
      roc%rho=0.0
      roc%ib=0
      roc%is=0
      roc%i1=0
      roc%i2=0
    
      iroc=1

      nff=2

      if (nff.eq.2) then 
       fnamex='2phonon/2f_x.dat'
       fnamec='2phonon/2f_c.dat'
      endif

      allocate(rog(0:isimax,nphon(1)))
      rog=0.d0

      allocate(ipozl(nphon(2)))
      ipozl=0

      call read_ampl(fnamex,x_ampl,ns2,nphon(1),nphon(1),iphous2,ipozl)
      
      write(*,*)'Proces ',myid, ' loading 1-phonon densities'

      call  readro1_all(rph,nphon(1),ihmn,ipmx,isimax,myid)

      ndrh=20000000
      ndrp=50000000
      allocate(r2p(ndrp),r2h(ndrh))

! MPI paralelization 
      allocate(ig_resh(ns1))
      do i=1,ns1
       ig_resh(i)=i
      enddo

      call rperm(ns1, ig_resh)

!      do ig=1,ifmx
!       do ig=ig_min,ig_max

      if (mod(ns1,numprocs).eq.0) then
              n_seg=ns1/numprocs
      else
              n_seg=ns1/numprocs+1
      endif

      if (myid.eq.0) write(*,*) ' size of segment = ',n_seg


do ia_cal=myid*n_seg+1,min((myid+1)*n_seg,ns1)
   ndr2p=0
   ndr2h=0
      
   ialpha=iphous(ig_resh(ia_cal))
   jalpha=phonbs(nff,ialpha)%jj
   itzalpha=phonbs(nff,ialpha)%tz
   facta=(dfloat(2*jalpha+1))**0.5d0

   write(*,*) ' Process #',myid,'  calculating ia = ',ialpha, 'Phonon energy ',phonbs(nff,ialpha)%enf
   call read_amp_ia(ialpha,fnamec,c_ampl,ndc,ipa,ja,nnz_c,nnz_x)
 
   
   do ibeta=1,nphon(nff)
      jbeta=phonbs(nff,ibeta)%jj
      itzbeta=phonbs(nff,ibeta)%tz
      iroc=0

      do ilp=1,nphon(1) ! lambda'

         rog=0.d0
         jlp=phonbs(1,ilp)%jj

         do i=1,nnz_c  ! sum over all nonzero C
            il1=c_ampl(i)%is
            il2=c_ampl(i)%ig
            jl1=phonbs(1,il1)%jj
            jl2=phonbs(1,il2)%jj
            fact=facta*(-1)**(jl1+jbeta)
      
            rog(0:isimax,il1)=rog(0:isimax,il1)+fact*c_ampl(i)%am*x_ampl(ibeta,ilp,il2)*csixj(0:isimax,jlp,jl1,jl2,jalpha,jbeta)

            fact=facta*(-1)**(jlp+jalpha)
            
            rog(0:isimax,il2)=rog(0:isimax,il2)+fact*c_ampl(i)%am*x_ampl(ibeta,il1,ilp)*csixj(0:isimax,jlp,jl2,jl1,jalpha,jbeta)
          enddo ! over i

           do ila=1,nphon(1)
            do isi=0,isimax

             if (dabs(rog(isi,ila)).gt.1.d-12) then
              iroc=iroc+1
              roc(iroc)%rho=rog(isi,ila)
              roc(iroc)%ib=ibeta
              roc(iroc)%is=isi
              roc(iroc)%i1=ila
              roc(iroc)%i2=ilp

              if (iroc.gt.ndimroc) then
                 write(*,*)' Increase dimension of array roc '
                 stop
              endif

!              if ((ialpha == 350).and.(ibeta == 418)) write(991,'(5i5,f15.10)') ialpha,ibeta,isi,ila,ilp,roc(iroc)%rho
!              if ((ialpha == 350).and.(ibeta == 418)) write(991,'(5i5,f15.10)') ialpha,ibeta,isi,ila,ilp,roc(iroc)%rho
             endif

            enddo
          enddo
      enddo  !over il

      if (iroc.ne.0) then
       call rdens2(ialpha,ibeta,jalpha,jbeta,ihmn,ipmx,isimax,roc,iroc,nphon,rph,lev,par_lev,hol_lev,ip_lev,ih_lev,r2p,r2h,ndr2p,ndr2h)     
      endif

      enddo ! cycle over ibbb

if (ndr2p.ne.0) then
 call  write_dens(ialpha,ndr2p,r2p,'2f_rp.dat')
endif 

if (ndr2h.ne.0) then
 call  write_dens(ialpha,ndr2h,r2h,'2f_rh.dat')
endif 

write(63,*)ialpha

enddo ! cycle over ia_cal

!       enddo

!      enddo

return
end subroutine r_value

!************************************************************************
subroutine write_dens(iaaa,iirg,rh,fname)

implicit double precision (a-h,o-z)

include 'types_phon_dens.inc'
include 'formats_phon_dens.inc'

type (ro2_typ), dimension(:), allocatable :: rh
character(len=9) fname
character(len=6) nlam

if (iirg.gt.0) then
  write(nlam,'(i6.6)')iaaa                  
  open(73,file='scratch/'//fname//'_'//nlam,status='unknown',form='unformatted')
  write(73)iaaa,iirg
  write(73)(rh(iii)%ib,iii=1,iirg)
  write(73)(rh(iii)%is,iii=1,iirg)
  write(73)(rh(iii)%i1,iii=1,iirg)
  write(73)(rh(iii)%i2,iii=1,iirg)
  write(73)(rh(iii)%rho,iii=1,iirg)
  close(73)
endif

return 
end subroutine write_dens

subroutine rdens2(ia,ib,ja,jb,imin,imax,isimax,rom,irom,nphon,rph,lev,par_lev,hol_lev,ip_lev,ih_lev,r2p,r2h,idimp,idimh)

use anglib   ! angular momentum staff

implicit double precision(a-h,o-z)

include 'types_phon_dens.inc'
!include 'input_phon_dens.inc'
include 'formats_phon_dens.inc'


type (level_typ), dimension(:),allocatable :: lev
type (ro2_typ), dimension(:) :: r2p,r2h
integer, dimension(:), allocatable :: par_lev, hol_lev

type (roc_typ), dimension(:), allocatable :: roc,rom

integer, dimension (:), allocatable :: nphon
real(kind=4), dimension (:,:,:), allocatable :: rr
real, dimension(:,:,:,:,:), allocatable :: rph
character(len=10) fnamern,name_new
character(len=2) denst
integer*8 irom,ii
character*6 nlam


if (irom.ne.0) then

  nff=2
  ndla=nphon(nff-1)
  if (.not.allocated(rr)) allocate(rr(imin:imax,imin:imax,0:isimax))
  rr=0.d0

  do ii=1,irom
     isir=rom(ii)%is
     iba=rom(ii)%ib
     ilam=rom(ii)%i1
     ilamp=rom(ii)%i2
     rr(imin:imax,imin:imax,isir)=rr(imin:imax,imin:imax,isir)+rph(imin:imax,imin:imax,isir,ilam,ilamp)*rom(ii)%rho
     xpom=rr(1,1,0)
  enddo ! over ii

endif ! irom
!       deallocate(rr)

do isi=0,isimax
!   if(abs(jb-ja).le.isi.and.isi.le.(jb+ja)) then
     do i_1=1,ip_lev
       ind1=par_lev(i_1)
       do i_2=1,ip_lev
        ind2=par_lev(i_2)
!         if(abs(lev(ind1)%j-lev(ind2)%j).le.2*isi.and.2*isi.le.(lev(ind1)%j+lev(ind2)%j)) then
 
          if (abs(rr(ind1,ind2,isi))>1.d-10) then 
          idimp=idimp+1

         if (idimp.gt.size(r2p)) then
           write(*,*)'Increase dimension of rp in rdens2'
           stop
         endif

         r2p(idimp)%ib=iba
         r2p(idimp)%is=int(isi,1)
         r2p(idimp)%i1=int(ind1,2)
         r2p(idimp)%i2=int(ind2,2)
         r2p(idimp)%rho=rr(ind1,ind2,isi)
!         endif
         endif
       enddo ! i_1
     enddo ! i_2
!    end if
enddo

do isi=0,isimax
   if(abs(jb-ja).le.isi.and.isi.le.(jb+ja)) then
     do i_1=1,ih_lev
       ind1=hol_lev(i_1)
       do i_2=1,ih_lev
        ind2=hol_lev(i_2)
         if(abs(lev(ind1)%j-lev(ind2)%j).le.2*isi.and.2*isi.le.(lev(ind1)%j+lev(ind2)%j)) then
 
         if (abs(rr(ind1,ind2,isi))>1.d-10) then 
          idimh=idimh+1

         if (idimh.gt.size(r2h)) then
           write(*,*)'Increase dimension of rp in write_dens'
           stop
         endif

         r2h(idimh)%ib=iba
         r2h(idimh)%is=int(isi,1)
         r2h(idimh)%i1=int(ind1,2)
         r2h(idimh)%i2=int(ind2,2)
         r2h(idimh)%rho=real(rr(ind1,ind2,isi))
         endif
         endif
       enddo ! i_1
     enddo ! i_2
    end if
enddo


return
end subroutine rdens2


      subroutine loadphon(phonbs,nphon)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      character*30 namef

      type (phon_typ), dimension (:,:), allocatable :: phonbs
      integer, dimension (:), allocatable :: nphon


      allocate (nphon(0:4))
      nphon=0
      allocate (phonbs(0:4,1:1000000))

      do nfon=1,2

      if (nfon.eq.1) namef='1phonon/1f_states.dat'
      if (nfon.eq.2) namef='2phonon/2f_states.dat'
      
      open(1,file=namef,status='old',form='unformatted')

      do while (.not.eof(1))
       read(1)i,ipartt,ijj,en,itt
         phonbs(nfon,i)%par=ipartt
         phonbs(nfon,i)%jj=ijj
         phonbs(nfon,i)%enf=en
         phonbs(nfon,i)%tz=itt
      enddo
      nphon(nfon)=i

      close(1)

      enddo

      end subroutine loadphon
!*********************************************************************
      subroutine selphon(nf,phonbs,nphon,iphous,ns1,myid)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      character*30 namef

      type (phon_typ), dimension (:,:), allocatable :: phonbs
      integer, dimension (:), allocatable :: nphon,iphous,list_2ph


      allocate(iphous(nphon(nf)))
      iphous=0


!      xthr_min=0.0d0
!      xthr_max=10000000000000000.0d0

      allocate(list_2ph(nphon(nf)))
      list_2ph=0

      open(881,file='2_phon_dens_list.dat',status='old',form='formatted')
      do while(.not.eof(881))
         read(881,*)ii
         list_2ph(ii)=1
      enddo
      close(881)


      ii=0
      do i=1,nphon(nf)

!      xene=phonbs(nf,i)%enf
!       jjf=phonbs(nf,i)%jj
!        do j=1,nlist
!       if (xene.le.xthr_max.and.xene.gt.xthr_min) then
        if (list_2ph(i).eq.1) then
          ii=ii+1
          iphous(ii)=i
        endif
!          enddo
      enddo

!      write(*,*)'Energy threshold for 2 phon. dens.',xthr_min,xthr_max
      if(myid.eq.0) write(*,*)' Number of selected phonons a)',ii

      ns1=ii

      end subroutine selphon

!**********************************************************************
      subroutine selphon2(nf,phonbs,nphon,iphous,ns2,myid)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      character*30 namef

      type (phon_typ), dimension (:,:), allocatable :: phonbs
      integer, dimension (:), allocatable :: nphon,iphous


      allocate(iphous(nphon(nf)))
      iphous=0



      ii=0
      do i=1,nphon(nf)
       xene=phonbs(nf,i)%enf
       jjf=phonbs(nf,i)%jj

       if (xene.gt.-10.0d0) then
        ii=ii+1
        iphous(i)=1
       endif
      enddo

      if (myid.eq.0) write(*,*)' Number of selected phonons b) ',ii
      ns2=ii

      end subroutine selphon2

!**********************************************************************
    

      subroutine readx(fname,xcc,ndcc,ipcal,jcal,ilamps,no,idphon)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      type (amp2_typ), dimension(:,:), allocatable :: xcc
      integer, dimension (:), allocatable :: ndcc

      character(len=30)fname


      open(2,file=fname,status='old',form='unformatted')

      ilamp=0

      ndimj=100000


      do while (.not.eof(2))

      read(2)ipar,ijj,no,idphon

      if (allocated(xcc).eq..TRUE.) deallocate(xcc,ndcc)

      allocate(xcc(ilamp+1:ilamp+no+1,ndimj))
      allocate(ndcc(ilamp+1:ilamp+1+no))

      xcc%is=0
      xcc%ig=0
      xcc%am=0.d0
      ndcc=0


      do ilam=1,no

       
         if (idphon.gt.ndimj) then 
      write(*,*)'Readcam: allocate bigger array in ndimj',idphon,ndimj
               stop
           endif
  
       read(2)(xcc(ilam+ilamp,i)%ig,xcc(ilam+ilamp,i)%is,xcc(ilam+ilamp,i)%am,i=1,idphon)
       ndcc(ilam+ilamp)=idphon
      enddo

      ilamps=ilamp
      ilamp=ilamp+no


      if (ipar.eq.ipcal.and.ijj.eq.jcal) goto 11


      deallocate(xcc,ndcc)
      no=0
      idphon=0

      enddo

  11  continue     

      close(2)

      return
      end subroutine readx
!************************************************************
      subroutine readx_ia(ia,fname,xcc,ndcc,ipcal,jcal,ilamps,no,idphon)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      type (amp2_typ), dimension(:,:), allocatable :: xcc
      integer, dimension (:), allocatable :: ndcc

      character(len=30)fname

      open(2,file=fname,status='old',form='unformatted')

      ilamp=0
      ndimj=100000
      ii=0

      if (allocated(xcc)) deallocate(xcc)

      do while (.not.eof(2))

        read(2)ipar,ijj,no,idphon
        ii=ii+no

        if (ia.gt.ii) then 

          do ilam=1,no
            read(2)
          enddo

        else


          allocate(xcc(1,no))

          do ilam=1,ia-ii+no-1
            read(2)
          enddo
          
          read(2)(xcc(1,i)%ig,xcc(1,i)%is,xcc(1,i)%am,i=1,idphon)
          close(2) 
          return 

        endif 
      enddo


      return
      end subroutine readx_ia


!******************************************************************************
      subroutine readro1_all(rph,ndla,ihmn,ipmx,isimax,myid)

      implicit double precision (a-h,o-z)

!      include 'types_eqm.inc'
      include 'types_phon_dens.inc'

!      double precision, dimension(:,:,:,:,:), allocatable :: rpp,rnp,rph,rnh
      real, dimension(:,:,:,:,:), allocatable :: rph
      type(rho_typ), dimension(:), allocatable :: ronn
      logical je_tam_subor


      character(len=9)fname
      character(len=4)nlam

      if (.not.allocated(rph)) allocate(rph(ihmn:ipmx,ihmn:ipmx,0:isimax,ndla,ndla))
      rph=0.0
      
      ndro=5000000
      if (.not.allocated(ronn)) allocate (ronn(ndro))
 
      ifile=33


      do ityp=1,2
       if (ityp.eq.1) fname='1f_rp.dat'
       if (ityp.eq.2) fname='1f_rh.dat'

      ifile=33+myid

      do ig=1,ndla

      write(nlam,'(i4.4)')ig

      inquire(file='scratch/'//fname//'_'//nlam,exist=je_tam_subor)

      if (je_tam_subor.eq..FALSE.) then
        ndgg=0
        write(*,*)'File scratch/'//fname//'_'//nlam,'not found! '
        return
      endif


      open(ifile,file='scratch/'//fname//'_'//nlam,status='unknown',form='unformatted')

      ndgg=0

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
!      if (.not.allocated(ron)) allocate(ron(ndla,ndla,0:nsi,n1mn:n1mx,n2mn:n2mx))

!      ron=0.d0

       do ii=1,ndgg

        ibt=ronn(ii)%ilap
        isit=ronn(ii)%j
        i1t=ronn(ii)%i1
        i2t=ronn(ii)%i2
        rot=ronn(ii)%ro

!        if (ibt.gt.ndla.or.isit.gt.nsi.or.i1t.lt.n1mn.or.i1t.gt.n1mx.or.i2t.lt.n2mn.or.i2t.gt.n2mx) then
!         write(*,*)' Small dimensions of ron in read11',ibt,isit,i1t,i2t
!         stop
!       endif
!       if (isit.le.isimax) then
        rph(i1t,i2t,isit,ig,ibt)=real(rot)
!       endif
       enddo
       enddo
       enddo

      write(*,*)' Proces ',myid,' 1 phon. dens. loaded'
      return
      end subroutine readro1_all


!******************************************************************************

      subroutine readc2_part(fname,cc,n1,n2,n3,iphous2,ipozl,myid)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      type (amp2_typ), dimension(:,:), allocatable :: xcc
!      double precision, dimension (:,:,:), allocatable :: cc
      double precision, dimension (:,:,:), allocatable :: cc
      integer, dimension (:), allocatable :: ndcc,iphous2,ipozl
!      real :: xxr
      character(len=30)fname

!      ndimi=2000
!      ndimj=600000

      write(*,*)'Process ',myid,' reading X values'


      allocate (cc(n1,n2,n3))
      cc=0.0
      illl=1

      write(*,*)'Process ',myid,' cc allocated'

      open(2,file=fname,status='old',form='unformatted')

      ilamp=0

      do while (.not.eof(2))

       read(2)ipar,ijj,no,idphon

!      if (ipar.eq.-1.and.ijj.eq.1) then
!     read(2)ipar,ijj,no,idphon
      if (allocated(xcc).eq..TRUE.) deallocate(xcc,ndcc)

      allocate(xcc(1,idphon))
      allocate(ndcc(1))

      do ilam=1,no

!         if (idphon.gt.ndimj) then
!      write(*,*)'Readcam: allocate bigger array in ndimj',idphon,ndimj
!               stop
!           endif


       if (iphous2(ilam+ilamp).eq.1) then

       read(2)(xcc(1,i)%ig,xcc(1,i)%is,xcc(1,i)%am,i=1,idphon)
       ndcc(1)=idphon


       do i=1,idphon
        ii=xcc(1,i)%is
        jj=xcc(1,i)%ig
        xxr=xcc(1,i)%am
        cc(illl,ii,jj)=xxr
        ipozl(ilam+ilamp)=illl
       enddo

       illl=illl+1
       else

       read(2)

       endif
      enddo

 !     endif

      ilamps=ilamp
      ilamp=ilamp+no

      enddo
  11  continue
      close(2)

      deallocate(xcc,ndcc)
      return
      end subroutine readc2_part
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

      subroutine read_ampl(fname,cc,n1,n2,n3,iphous2,ipozl)

            implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

            type (amp2_typ), dimension(:), allocatable :: xcc
            double precision, dimension (:,:,:), allocatable :: cc
            integer, dimension (:), allocatable :: ndcc,iphous2,ipozl
            real :: xxr
            character(len=30)fname

            allocate (cc(n1,n2,n3))
            cc=0.0d0
            illl=1

            open(2,file=fname,status='old',form='unformatted')

            ilamp=0

            do while (.not.eof(2))

             read(2)ipar,ijj,no,idphon
             if (.not.allocated(xcc)) allocate(xcc(idphon))
             read(2)(xcc(i)%ig,xcc(i)%is,i=1,idphon)

            do ilam=1,no

             if (iphous2(ilam+ilamp).eq.1) then

             read(2)(xcc(i)%am,i=1,idphon)

             do i=1,idphon
              cc(illl,xcc(i)%is,xcc(i)%ig)=xcc(i)%am
              ipozl(ilam+ilamp)=illl
             enddo
             illl=illl+1
             else

             read(2)

             endif
            enddo

       !     endif

            ilamps=ilamp
            ilamp=ilamp+no

            deallocate(xcc)

            enddo

            close(2)

            return
            end subroutine read_ampl
!************************************************************
     subroutine read_amp_ia(ia,fname,xcc,ndcc,ipar,ijj,no,idphon)

      implicit double precision (a-h,o-z)
      include 'types_phon_dens.inc'
      
            type (amp2_typ), dimension(:), allocatable :: xcc
            integer, dimension (:), allocatable :: ndcc

            character(len=30)fname

            open(2,file=fname,status='old',form='unformatted')

            ilamp=0
            ndimj=100000
            ii=0



            do while (.not.eof(2))

              read(2)ipar,ijj,no,idphon
      !        if (allocated(xcc)) deallocate(xcc)
      !        allocate(xcc(idphon))
      !        read(2)(xcc(i)%ig,xcc(i)%is,i=1,idphon)

              ii=ii+no

              if (ia.gt.ii) then
                do ilam=1,no+1
                  read(2)
                enddo

              else
              if (allocated(xcc)) deallocate(xcc)
                allocate(xcc(idphon))
                read(2)(xcc(i)%ig,xcc(i)%is,i=1,idphon)  ! ig -correspodns to gamma, is to lambda
              do ilam=1,ia-ii+no-1
                 read(2)
              enddo

                read(2)(xcc(i)%am,i=1,idphon)
                close(2)
               return
              endif

            enddo


            return
            end subroutine read_amp_ia



end module rdens
      
