!*     Program FMAT generates F-matrix elements in Tz, J-coupled of two-body interaction
!*     input from HFB code       
!*     last modification 22.1.2025

      program fmatr

      use anglib   ! angular momentum staff

      implicit double precision (a-h,o-z)

 !     include 'formats_fmat.inc'
      include 'types_fmat.inc'

      integer (kind=1) :: ip,ij,it
      integer (kind=2) :: ia,ib,ic,id
      real (kind=8) :: vv,fmat

      character(len=30)file_int,file_sp,file_fp,file_fn,file_fpn,file_fnp

      type(level_typ),dimension(:),allocatable :: lev,levn,levp

      double precision, dimension(:,:,:,:,:), allocatable :: vp,vn,vpn

      neerno=0 ! error info 

      open(3,status='unknown',file='Input_files',form='formatted')

      write(*,*)' Name of interaction file'
      read(3,*)
      read(3,'(A )')file_int
      write(*,*)file_int
      write(*,*)
      write(*,*)' Name of output Fmat  interaction file'
      read(3,*)
      read(3,'(A )')file_fp
      write(*,*)file_fp

      close(3)

      nlevdim=10000
      allocate(lev(nlevdim),levn(nlevdim),levp(nlevdim))
      lev%n=0
      levn%n=0
      levp%n=0
      jmax=0

      open(1,file='HF_p.out',status='unknown',form='formatted')
      read(1,*)
        do while (.not.eof(1))
!         read(1,*) lhfp(ii)%index,lhfp(ii)%l,lhfp(ii)%j2,
!     &    lhfp(ii)%ei,lhfp(ii)%qei,lhfp(ii)%vi

         read(1,*)ii,ll,jj,ee,qei,vi

!        lev(i)%n=nn
        i=2*ii-1
        lev(i)%l=ll
        lev(i)%j=jj
        lev(i)%en=ee
        lev(i)%it=-1

        levp(ii)%l=ll
        levp(ii)%j=jj
        levp(ii)%en=ee
         if (jj.gt.jmax) jmax=jj

        enddo
       close(1)

       nlevmaxp=ii
       write(*,*)' Number of proton  levels',nlevmaxp


       open(1,file='HF_n.out',status='unknown',form='formatted')
       read(1,*)
        do while (.not.eof(1))
!         read(1,*) lhfp(ii)%index,lhfp(ii)%l,lhfp(ii)%j2,
!     &    lhfp(ii)%ei,lhfp(ii)%qei,lhfp(ii)%vi

         read(1,*)ii,ll,jj,ee,qei,vi

!        lev(i)%n=nn
        i=2*ii
        lev(i)%l=ll
        lev(i)%j=jj
        lev(i)%en=ee
        lev(i)%it=1

        levn(ii)%l=ll
        levn(ii)%j=jj
        levn(ii)%en=ee

        if (jj.gt.jmax) jmax=jj

        end do
       close(1)

      nlevmaxn=ii

      write(*,*)' Number of neutron levels',nlevmaxn
      write(*,*)'Maximal J',jmax
      

      nlev=nlevmaxn
      if (nlevmaxp.gt.nlev) nlev=nlevmaxp

      write(*,*)'nlev_max?'
      read(*,*)nlev_max
      write(*,*)nlev_max
      nlev=nlev_max

      open(1,file='singpart_coup.dat',status='unknown',form='formatted')
      write(1,*)'  tz   n    l    2j '
      do i=1,nlev
       write(1,'(4i5,7x,f15.5)')-1,levp(i)%n,levp(i)%l,levp(i)%j,levp(i)%en
       write(1,'(4i5,7x,f15.5)')1,levn(i)%n,levn(i)%l,levn(i)%j,levn(i)%en       
      enddo


      jmax=0
      do i=1,nlev
      if (levp(i)%j.gt.jmax) jmax=levp(i)%j
      enddo

      do i=1,nlev
       if (levn(i)%j.gt.jmax) jmax=levn(i)%j
      enddo
   
      write(*,*)'F matrix J_max',jmax
 
      nlev=2*nlev

      open(2,status='unknown',file=file_int,form='unformatted')
      open(3,status='unknown',file=file_fp,form='unformatted')

      write(*,*)
      jmin=0

      allocate(vp(0:jmax,1:nlev,1:nlev,1:nlev,1:nlev))
      vp=0.d0
 
      itzt=-1   ! Proton-proton interaction
      jtotmax=0
      do while (.not.eof(2))
        read(2)it,ij,ia,ib,ic,id,vv
        if (it.eq.itzt) then 
           if (ij/2.le.jmax) then
 
           if (ij.gt.jtotmax) jtotmax=ij
           if (jtotmax/2.gt.jmax) then 
                   write(*,*)' Jmax larger then in input'
                   stop
           endif

           ifazab=(-1)**((lev(ia)%j+lev(ib)%j-ij)/2)
           ifazcd=(-1)**((lev(ic)%j+lev(id)%j-ij)/2)


           if (ia.le.nlev.and.ib.le.nlev.and.ic.le.nlev.and.id.le.nlev) then 

           factab=1.d0
           factcd=1.d0
           if (ia.eq.ib) factab=dsqrt(2.d0)
           if (ic.eq.id) factcd=dsqrt(2.d0)
           xnorm=factab*factcd
           vv=vv*xnorm

           
           vp(ij/2,ia,ib,ic,id)=vv
           vp(ij/2,ic,id,ia,ib)=vv
           vp(ij/2,ib,ia,ic,id)=dfloat(-1*ifazab)*vv
           vp(ij/2,ic,id,ib,ia)=dfloat(-1*ifazab)*vv
           vp(ij/2,ia,ib,id,ic)=dfloat(-1*ifazcd)*vv
           vp(ij/2,id,ic,ia,ib)=dfloat(-1*ifazcd)*vv
           vp(ij/2,ib,ia,id,ic)=dfloat(ifazcd*ifazab)*vv
           vp(ij/2,id,ic,ib,ia)=dfloat(ifazcd*ifazab)*vv

           endif
          endif
         endif
      enddo

      jtotmax=jtotmax/2 

      write(*,*)'J_max of loaded PP interaction m.e ',jtotmax
      write(*,*)

      rewind(2)

      itzt=1   ! Neutron-neutron interaction
      jtotmax=0
      do while (.not.eof(2))
        read(2)it,ij,ia,ib,ic,id,vv

        if (it.eq.itzt) then
          if (ij/2.le.jmax) then

           if (ij.gt.jtotmax) jtotmax=ij
           if (jtotmax/2.gt.jmax) then
                   write(*,*)' Jmax larger then in input'
                   stop
           endif

           ifazab=(-1)**((lev(ia)%j+lev(ib)%j-ij)/2)
           ifazcd=(-1)**((lev(ic)%j+lev(id)%j-ij)/2)


           if (ia.le.nlev.and.ib.le.nlev.and.ic.le.nlev.and.id.le.nlev) then

           factab=1.d0
           factcd=1.d0
           if (ia.eq.ib) factab=dsqrt(2.d0)
           if (ic.eq.id) factcd=dsqrt(2.d0)
           xnorm=factab*factcd
           vv=vv*xnorm


           vp(ij/2,ia,ib,ic,id)=vv
           vp(ij/2,ic,id,ia,ib)=vv
           vp(ij/2,ib,ia,ic,id)=dfloat(-1*ifazab)*vv
           vp(ij/2,ic,id,ib,ia)=dfloat(-1*ifazab)*vv
           vp(ij/2,ia,ib,id,ic)=dfloat(-1*ifazcd)*vv
           vp(ij/2,id,ic,ia,ib)=dfloat(-1*ifazcd)*vv
           vp(ij/2,ib,ia,id,ic)=dfloat(ifazcd*ifazab)*vv
           vp(ij/2,id,ic,ib,ia)=dfloat(ifazcd*ifazab)*vv

            endif
          endif
         endif
      enddo

      jtotmax=jtotmax/2

      write(*,*)'J_max of loaded NN interaction m.e ',jtotmax
      write(*,*)

      rewind(2)

      itzt=0   ! Proton-neutron interaction
      jtotmax=0
      do while (.not.eof(2))


        read(2)it,ij,ia,ib,ic,id,vv


        if (it.eq.itzt) then
         if (ij/2.le.jmax) then
           if (ij.gt.jtotmax) jtotmax=ij
           if (jtotmax/2.gt.jmax) then
                   write(*,*)' Jmax larger then in input'
                   stop
           endif

           ifazab=(-1)**((lev(ia)%j+lev(ib)%j-ij)/2)
           ifazcd=(-1)**((lev(ic)%j+lev(id)%j-ij)/2)


           if (ia.le.nlev.and.ib.le.nlev.and.ic.le.nlev.and.id.le.nlev) then


           factab=1.d0
           factcd=1.d0
!           if (ia.eq.ib) factab=dsqrt(2.d0)
!           if (ic.eq.id) factcd=dsqrt(2.d0)
           xnorm=factab*factcd
!           vv=-1.d0*vv


           vp(ij/2,ia,ib,ic,id)=vv
           vp(ij/2,ic,id,ia,ib)=vv
           vp(ij/2,ib,ia,ic,id)=dfloat(-1*ifazab)*vv
           vp(ij/2,ic,id,ib,ia)=dfloat(-1*ifazab)*vv
           vp(ij/2,ia,ib,id,ic)=dfloat(-1*ifazcd)*vv
           vp(ij/2,id,ic,ia,ib)=dfloat(-1*ifazcd)*vv
           vp(ij/2,ib,ia,id,ic)=dfloat(ifazcd*ifazab)*vv
           vp(ij/2,id,ic,ib,ia)=dfloat(ifazcd*ifazab)*vv




           endif
          endif
         endif
      enddo

      jtotmax=jtotmax/2

      write(*,*)'J_max of loaded interaction PN  m.e ',jtotmax
      write(*,*)




      write(*,*) 'F matrix calculation'



      do isi=jmin,jmax

       do ii=1,nlev
         write(1234,*)isi,jmax,ii,nlevmax
          jii=lev(ii)%j
        do jj=1,nlev
            jjj=lev(jj)%j
         do kk=1,nlev
             jkk=lev(kk)%j
          do ll=1,nlev
               jll=lev(ll)%j

               ipari=(-1)**(lev(ii)%l+lev(kk)%l)

             fmat=0.d0  
             do j=jmin,jmax

             ifaz=(-1)**((jjj+jkk)/2-isi-j) 
            sixjc=sixj(jii,jjj,2*isi,jll,jkk,2*j) 
            fmat=fmat+dfloat((2*j+1)*ifaz)*sixjc*vp(j,ii,kk,jj,ll)

             enddo
          
!        if (dabs(fmat).gt.1.d-10) write(3)-1,ipari,isi,ii,jj,kk,ll,fmat 
         if (dabs(fmat).gt.1.d-10) write(3)int(isi,1),int(ii,2),int(jj,2),int(kk,2),int(ll,2),fmat
          enddo
         enddo
        enddo
       enddo



      enddo 


      deallocate(vp)


      end
