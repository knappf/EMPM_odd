module full_mat_hole

contains 

subroutine bstot(nfon,ipar,jcal,tzcal,nbs,ndmpho,phonbs,par_lev,hol_lev)

implicit none

include 'types_fullham.inc'

character*30 namef

type (phon_typ), dimension (:,:), allocatable :: phonbs
integer, dimension(:), allocatable :: nphon,ndmpho
integer, dimension(:,:), allocatable :: nbs
integer :: nf_max,nfon,ipar,iparl,jcal,tzcal,i,ii,nlev,jmax,ip_lev,ih_lev,ih_levp,ia,iz
type(level_typ),dimension(:), allocatable :: lev
integer, dimension(:,:), allocatable :: i_ph_int
integer, dimension(:), allocatable :: par_lev, hol_lev


nf_max=50000

if (.not.allocated(phonbs)) allocate(phonbs(0:4,nf_max))

      if(.not.allocated(nphon)) then
         allocate(nphon(0:4))
         nphon=0
      endif

      if (.not.allocated(ndmpho)) then
         allocate(ndmpho(0:4))
         ndmpho=0
      endif

      if (.not.allocated(nbs)) then
         allocate(nbs(0:4,nf_max))
         nbs=0
      endif

 if (nfon.eq.0) then
!              phonbs(0,1)%par=1
!              phonbs(0,1)%jj=0
!              phonbs(0,1)%enf=0.d0
!              nphon(0)=1

allocate(i_ph_int(4,-1:1))
write(*,*)'Loading of input '

open(1,file='input_tda_coup.dat',status='old',form='formatted')

read(1,'(30i6)')ia,iz
read(1,'(30i6)')i_ph_int(1,1),i_ph_int(2,1)
read(1,'(30i6)')i_ph_int(1,-1),i_ph_int(2,-1)
read(1,'(30i6)')i_ph_int(3,1),i_ph_int(4,1)
read(1,'(30i6)')i_ph_int(3,-1),i_ph_int(4,-1)

close(1)



 nlev=300
 allocate(lev(nlev))

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


! if (tzcal.eq.-1) then    
   do i=1,ih_lev
     iparl=(-1)**lev(hol_lev(i))%l
!     if (lev(hol_lev(i)%j.eq.jcal.and.iparl.eq.iparc) then  
     phonbs(0,i)%par=iparl
     phonbs(0,i)%jj=lev(hol_lev(i))%j
     phonbs(0,i)%tz=1
     if (mod(hol_lev(i),2).eq.0) phonbs(0,i)%tz=-1     
     phonbs(0,i)%enf=lev(hol_lev(i))%en
!     endif
   enddo

  nphon(nfon)=ih_lev

 endif 


 if (nfon.gt.0) then

        if (nfon.eq.1) namef='1phon_hole/1f_states.dat'
!        if (nfon.eq.2) namef='2phonon/2f_states.dat'
!        if (nfon.eq.3) namef='3phonon/3f_states.dat'

    open(1,file=namef,status='old',form='unformatted')

      do while (.not.eof(1))
       read(1)i,phonbs(nfon,i)%par,phonbs(nfon,i)%jj,phonbs(nfon,i)%tz,phonbs(nfon,i)%enf
      enddo
      nphon(nfon)=i

    close(1)

 endif

      ii=0

 do i=1,nphon(nfon)

      if (ipar.eq.phonbs(nfon,i)%par.and.jcal.eq.phonbs(nfon,i)%jj.and.tzcal.eq.phonbs(nfon,i)%tz) then
              ii=ii+1
              nbs(nfon,ii)=i
!         if (nfon.eq.0) nbs(nfon,ii)=hol_lev(i)
      endif
 enddo

   ndmpho(nfon)=ii


return
end subroutine bstot


subroutine inp_sp(lev,nlev,jmax)

      implicit none

      include 'types_fullham.inc'

      integer :: i,ii,jmax,nlev
      integer, dimension(:), allocatable:: il_n,il_p
      type(level_typ),dimension(:), allocatable :: lev

      allocate(il_n(0:1000))
      il_n=0

      allocate(il_p(0:1000))
      il_p=0

      jmax=0

      open(1,file='singpart_coup.dat',status='old',form='formatted')

      read(1,*)
!      do i=1,ipnmx
!       read(1,16)nt,lt,jt,erspnt,ersppt
       i=0
       do while (.not.eof(1))
        i=i+1
        read(1,'(4i5,7x,f15.5)')lev(i)%tz,lev(i)%n,lev(i)%l,lev(i)%j,lev(i)%en

        if (lev(i)%j.gt.jmax) jmax=lev(i)%j



      enddo
      write(*,*)' Number of levels loaded ',i
      write(*,*)' Maximal value 2J',jmax
      nlev=i

      close(1)

   end subroutine inp_sp
!-----------------------------------------------------------------------------------------

subroutine inter_hol_phon1(ndmpho,nbs,xamp,ndx,vlam,jcal,par_lev,hol_lev)
implicit none 

include 'types_fullham.inc'

type (amp_typ), dimension(:,:), allocatable :: xamp
integer, dimension (:), allocatable :: ndx
integer, dimension(:), allocatable :: ndmpho
integer, dimension(:,:), allocatable :: nbs
integer, dimension(:), allocatable :: jphon
integer, dimension(:),allocatable :: ndc

type(amp_ph_typ), dimension(:,:), allocatable :: camp
integer :: i,ih,ila,jsp,i_jsp,ic,ipp,ihh,jc,inu,iu,jcal

double precision, dimension(:,:,:,:,:),allocatable ::fp
double precision, dimension(:,:), allocatable :: vlam
integer, dimension(:), allocatable :: par_lev, hol_lev

character*30 namefp


namefp='fmat.dat'
! upravit na dynamicku zmenu dimenzie poli
call readfin(namefp,0,30,1,100,1,100,fp)
call read_cph(camp,ndc,jphon)

allocate (vlam(ndmpho(1),ndmpho(0)))
vlam=0.0d0

do iu=1,ndmpho(1)
   inu=nbs(1,iu)
do jsp=1,ndmpho(0)
   i_jsp=hol_lev(nbs(0,jsp))
   

 
do i=1,ndx(inu)
  ih=xamp(inu,i)%isp
  ila=xamp(inu,i)%ila
  
! do jsp=1,ndmpho(0)
!    i_jsp=nbs(0,jsp)
    
!    vlam=0.d0
    do ic=1,ndc(ila)
      ipp=camp(ila,ic)%par
      ihh=camp(ila,ic)%hol
      jc=jphon(ila)
      vlam(iu,jsp)=vlam(iu,jsp)+0.5d0*camp(ila,ic)%am*fp(jc,ipp,ihh,i_jsp,ih)*xamp(inu,i)%am*dfloat(2*jc+1)**0.5d0/(dfloat(jcal+1))
    enddo
    
 enddo
  
enddo

enddo


return 
end subroutine inter_hol_phon1
!--------------------------------------------------------------------------------------
subroutine fulhamt(ndmpho,nbs,diag_h,vlam)


      implicit double precision (a-h,o-z)

      double precision, dimension(:,:), allocatable :: amtr
      double precision, dimension(:), allocatable :: wr,work

      integer, dimension(:), allocatable :: ndmpho
      integer,  dimension (:,:), allocatable :: nbs
      character*30 :: nameen,namewf,namecon
      double precision, dimension(:,:), allocatable :: vlam
      double precision, dimension(:),allocatable :: diag_h

!c     construction of total H matrix

      ndimtot=ndmpho(0)+ndmpho(1)

      allocate(amtr(ndimtot,ndimtot))
      amtr=0.d0

      do i=1,ndimtot
        amtr(i,i)=diag_h(i)
      enddo

! 0-1 phonon
        
      do i=1,ndmpho(1)
        do j=1,ndmpho(0)
    
        amtr(i+ndmpho(0),j)=vlam(i,j)
       enddo
      enddo
      

      iout=1

      if (iout.eq.1) then

  302 format(10000f13.4)
      do i=1,ndimtot
        write(203,302)(amtr(i,j),j=1,ndimtot)
      enddo

      endif

      allocate(wr(ndimtot))

      wr=0.d0
      lwork=26*ndimtot
      allocate(work(lwork))
      work=0.d0

      call DSYEV('V','L',ndimtot,amtr,ndimtot,wr,WORK,LWORK,INFO)

      write(*,*)' info =', info

      open(23,file='eigenenergies.dat',status='unknown',form='formatted')
      do i=1,ndimtot
        write(23,*)wr(i)
      enddo

      close(23)


      iphomax=1
      if (ndmpho(1).eq.0) iphomax=0

      write(*,*)' Phonmax ',iphomax


      if (iphomax.eq.0) then
             nameen='ef_0f.dat'
             namewf='wf_0f.dat'
             namecon='eigst_cont_0f.dat'
      endif


      if (iphomax.eq.1) then
             nameen='ef_1f.dat'
             namewf='wf_1f.dat'
             namecon='eigst_cont_1f.dat'
      endif


      open(23,file=namewf,status='unknown',form='unformatted')
      write(23)iphomax
      write(23)(ndmpho(ipho),ipho=0,iphomax)
      write(23)ndimtot
      write(23)((nbs(ipho,i),i=1,ndmpho(ipho)),ipho=0,iphomax)
      write(23)((amtr(i,j),i=1,ndimtot),j=1,ndimtot)


      write(201,302)(wr(i),i=1,ndimtot)
      write(201,*)
      do i=1,ndimtot
         write(201,302)(amtr(i,j),j=1,10)
      enddo


      close(23)



      open(23,file=nameen,status='unknown',form='unformatted')
      write(23)ndimtot
      write(23)(wr(i),i=1,ndimtot)
      close(23)


 274  format(i5,6f12.5)
      open(2,file=namecon,status='unknown',form='formatted')
      do i=1,ndimtot
        con0=0.d0
        con1=0.d0

        do k=0,iphomax

        aaa=0.d0

         if (k.eq.0) then
              jj=0
         else

           jj=0
           do kkk=0,k-1
            jj=jj+ndmpho(kkk)
           enddo
         endif

         do j=1+jj,ndmpho(k)+jj
           eee=amtr(j,i)
           aaa=aaa+amtr(j,i)**2.0d0
         enddo

         if (k.eq.0) con0=aaa
         if (k.eq.1) con1=aaa

         enddo

        write(2,274)i,wr(i),con0,con1
      enddo

      close(2)

      deallocate(amtr,wr,work)


      end subroutine fulhamt


!---------------------------------------------------------------------------------------
subroutine read_cph(camp,ndc,jphon)

implicit none 

include 'types_fullham.inc'

type(amp_ph_typ), dimension(:,:), allocatable :: camp
double precision, dimension(:,:,:), allocatable :: came
integer, dimension(:), allocatable :: jphon
integer :: iparc,ijj,ilam,ipar,ihol,i,ndamp
integer, dimension(:),allocatable :: ndc
double precision :: cc
character*10 namec


namec='1f_cph.dat'
open(11,file='1phonon/'//namec,status='old',form='unformatted')

ndamp=2000
allocate(camp(ndamp,ndamp))
allocate(ndc(ndamp),jphon(ndamp))

ilam=0
do while (.not.eof(11))
       ilam=ilam+1
       read(11)iparc,ijj,ndc(ilam)
       jphon(ilam)=ijj
        if (ndc(ilam).gt.ndamp) then
                           write(*,*)'small dimension of array camp'
                           stop
                   endif

       read(11)(camp(ilam,i)%par,camp(ilam,i)%hol,camp(ilam,i)%am,i=1,ndc(ilam))

!        do i=1,ndc
!         ipar=camp(i,ilam)%par
!         ihol=camp(i,ilam)%hol
!         cc=camp(i,ilam)%am
!         came(ilam,ipar,ihol)=cc
!        enddo

enddo

close(11)


end subroutine read_cph

!-------------------------------------------------------------------------------------------

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

!------------------------------------------------------------------------------------------

 subroutine readx(fname,xcc,ndcc,ipcal,jcal,itzcal)

      implicit double precision (a-h,o-z)

      include 'types_fullham.inc'

      type (amp_typ), dimension(:,:), allocatable :: xcc
!      type (ampr_typ), dimension(:,:), allocatable :: xccr

      integer, dimension (:), allocatable :: ndcc

      character(len=30)fname

!c      allocate(xcc(ndimi,ndimj))
!c      allocate(ndcc(ndimi))

      open(2,file=fname,status='old',form='unformatted')

      ilamp=0

!c      ndimj=533000



      do while (.not.eof(2))

      read(2)ipar,ijj,itz,no,idphon

      ndimj=idphon

!c      write(*,*)ipar,ijj,no,idphon

      if (allocated(xcc).eq..TRUE.) deallocate(xcc)
!     if (allocated(xccr).eq..TRUE.) deallocate(xccr)
      if (allocated(ndcc).eq..TRUE.) deallocate(ndcc)

!c      write(*,*)'ccc'
!      allocate(xcc(1,ndimj))
      allocate(xcc(ilamp:ilamp+no,ndimj))
      allocate(ndcc(ilamp:ilamp+no))
!c      write(*,*)'ddd'

!c      stop


      do ilam=1,no

!c      write(*,*)ilam,ilamp


         if (idphon.gt.ndimj) then
      write(*,*)'Readcam: allocate bigger array in ndimj',idphon,ndimj
               stop
           endif
       read(2)(xcc(ilam+ilamp,i)%isp,xcc(ilam+ilamp,i)%ila,xcc(ilam+ilamp,i)%am,i=1,idphon)
!       xccr(ilam+ilamp,:)%isp=xcc(1,:)%isp
!       xccr(ilam+ilamp,:)%ila=xcc(1,:)%ila
!       xccr(ilam+ilamp,:)%am=real(xcc(1,:)%am)

       ndcc(ilam+ilamp)=idphon
!c       read(2)(xcc(ilam+ilamp,i)%ig,xcc(ilam+ilamp,i)%is
!c     *,xcc(ilam+ilamp,i)%am,i=1,idphon)
!c       ndcc(ilam+ilamp)=idphon
      enddo

      ilamp=ilamp+no

      if (ipar.eq.ipcal.and.ijj.eq.jcal.and.itzcal.eq.itz) then 
       close(2)
       return 
      endif 

      deallocate(xcc,ndcc)
!      deallocate(xccr)

      enddo

      return


      end subroutine readx


end module full_mat_hole
