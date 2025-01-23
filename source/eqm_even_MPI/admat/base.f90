!c     last modification 14.4.2010  

      module phonon_base 

      contains

      subroutine phonbase(nf,ipare,j,itzcal,phonus,phonmus,idim1,idim2,idimbs,idphon,idphontr,phonbs,phon1,phon2,mxtr,myid)

      implicit double precision (a-h,o-z)

      include 'types_eqm.inc'

      integer, dimension (:), allocatable :: phonus, phonmus
      type(phon_typ), dimension (:), allocatable :: phon1,phon2
      type(phonbase_typ), dimension (:), allocatable :: phonbs
      integer, dimension(:), allocatable :: mxtr,nxtr

      character*30 name1f,name2f

!      open(23,file='phon_phon_base.dat',status='unknown',form='formatted')

      if (nf.eq.2) name1f='1phonon/1f_states.dat'
      if (nf.eq.2) name2f='1phonon/1f_states.dat'
      
      if (allocated(phon1).eq..TRUE.) then 
       deallocate(phon1,phon2,phonus,phonmus)
       deallocate(phonbs)
      endif

      allocate(phon1(0:idim1),phon2(0:idim2),phonus(0:idim1),phonmus(0:idim2))
      allocate(phonbs(idimbs))

!   ALL PHONONS     
      do i=1,idim1
      phon1(i)%us=1
      enddo

      do i=1,idim2
      phon2(i)%us=1
      enddo


      do i=1,idim1
      phonus(i)=i
      enddo

      do i=1,idim2
      phonmus(i)=i
      enddo

      open (3,file=name1f,status='old',form='unformatted')

      do while (.not.eof(3))
       read(3)i,ipar,ijj,en,itt
       phon1(i)%par=ipar
       phon1(i)%j=ijj
       phon1(i)%enf=en
       phon1(i)%tz=itt
!       phon1(i)%us=1
      enddo

      nnphontr=i

      close(3)

      open (3,file=name2f,status='old',form='unformatted')

      do while (.not.eof(3))
       read(3)i,ipar,ijj,en,itt
       phon2(i)%par=ipar
       phon2(i)%j=ijj
       phon2(i)%enf=en
       phon2(i)%tz=itt
!       phon2(i)%us=1
      enddo

      nnphonmtr=i

      close(3)

      i=0

      do ilamp=1,nnphontr

       do ilam=1,nnphonmtr

        ilamu=phonmus(ilam)    ! used phonons
        iparl=phon2(ilamu)%par
        jl=phon2(ilamu)%j
        itz=phon2(ilamu)%tz
        
!        do ilamp=1,nnphontr

         
        
          ilamup=phonus(ilamp)
          iparlp=phon1(ilamup)%par
          jlp=phon1(ilamup)%j
          itzp=phon1(ilamup)%tz

          ipart=iparl*iparlp

          jmin=abs(jl-jlp)
          jmax=jl+jlp
  
          if ((itz+itzp).eq.itzcal) then
!           if (itz.eq.0.and.itzp.eq.0) then
          if (ipart.eq.ipare) then 
          if (j.le.jmax.and.j.ge.jmin) then 
!          if (j.gt.jmax.or.j.lt.jmin) goto 6
!          if (phon2(ilamu)%us.ne.1) goto 6
!          if (phon1(ilamup)%us.ne.1) goto 6          
 
          i=i+1

          if (i.gt.idimbs) then 
                  write(*,*)'Small dimension of array phonbs'
                  stop
                endif
          
          phonbs(i)%ila=ilamu
          phonbs(i)%ilap=ilamup  

          endif
          endif 
          endif
        end do
        
      end do

      idphon=i

!      write(*,*)' Total dimension ',idphon

     
      allocate(mxtr(idphon),nxtr(idphon))
      mxtr=0
      nxtr=0

!      xetrun=25.0

!      write(*,*)' Energetic truncation? '
!      read(*,*)xetrun
!      write(*,*)xetrun
      xetrun=10000000.0d0

      jj=0
      do ii=1,idphon
         ilamu=phonbs(ii)%ila
         ilamup=phonbs(ii)%ilap
         eila=phon2(ilamu)%enf
         eilap=phon1(ilamup)%enf
         ius1=phon1(ilamup)%us
         ius2=phon2(ilamu)%us
         it1=phon2(ilamu)%tz
         it1p=phon1(ilamup)%tz
!         if (ius1.ne.0.and.ius2.ne.0) then
!         write(991,*)ilamu,eila,ilamup,eilap,ius1,ius2
!         endif

!         if ((it1+it1p).eq.itzcal) then 

!       if ((eila+eilap).le.xetrun.and.ilamu.le.ilamup.and.ius1.eq.1.and.ius2.ne.0.and.(ius1+ius2).le.3) then 
        if (ilamu.le.ilamup) then 
               jj=jj+1
               mxtr(jj)=ii
               nxtr(ii)=jj
!         write(23,'(3i5,f10.5,2i5,f10.5)')jj,ilamu,phon2(ilamu)%tz,eila,ilamup,phon1(ilamup)%tz,eilap
         endif 
!      endif

!  661 continue       

      enddo

      idphontr=jj

!      write(*,*) ' Truncated dimension ',idphontr

      deallocate(nxtr)
 
!c      deallocate(phon1,phon2,phonus,phonmus,phonbs)

!      close(23)
      return

      end subroutine phonbase

      end module phonon_base
