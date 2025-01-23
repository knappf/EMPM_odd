!     last modification 21.6.2010      
      
      program eqm 

      use phonon_base
      use choleski
      use hami
      use types_eqm

      implicit double precision (a-h,o-z)

!      include 'types_eqm.inc'

      type(phonbase_typ), dimension (:), allocatable :: phonbs
      type(phon_typ), dimension (:), allocatable :: phon1,phon2
      integer, dimension (:), allocatable :: phonus, phonmus

!     choleski arrays
      double precision, dimension(:,:), allocatable :: dd,cq,d1,cdu,xr,vr,h_corr
      double precision, dimension(:), allocatable :: wr
      
      integer, dimension(:), allocatable :: nx,mxt,irow,mxtr
      
      character*30 namex,names,namec

!      CALL OMP_SET_NUM_THREADS(16)

      nf=2
      idim1=10000
      idim2=10000
      idimbs=2000000
            
      nlam=0
      namex='2phonon/2f_x.dat'
      namec='2phonon/2f_c.dat'
      names='2phonon/2f_states.dat'
      
      open(12,file=namex,status='unknown',form='unformatted')
      open(22,file=namec,status='unknown',form='unformatted')
      open(13,file=names,status='unknown',form='unformatted')
      open(99,file='2phon.log',status='unknown',form='formatted')

      jmax=12
      jmin=0

      write(*,*)'ipmin,ipmax?'
      read(*,'(i5)')ipmin,ipmax
      write(*,*)'jmin,jmax?'
      read(*,'(i5)')jmin,jmax
      write(*,*)'tzmin,tzmax?'
      read(*,'(i5)')itmin,itmax

      write(*,*)'2-phonon calculation parity interval: ','<',ipmin,',',ipmax,'     >'
      write(*,*)'                         J  interval: ','<',jmin,',',jmax,'     >'
      write(*,*)'                         Tz interval: ','<',itmin,',',itmax,'     >'


      do ip=ipmin,ipmax,2    
      do jcal=jmin,jmax
      do itzcal=itmin,itmax

      write(*,*)
      write(*,*)'----------------------------------------------'
      write(*,*)'Tz = ',itzcal,' Par = ',ip,'   J = ',jcal

      call phonbase(nf,ip,jcal,itzcal,phonus,phonmus,idim1,idim2,idimbs,idphon,idphontr,phonbs,phon1,phon2,mxtr)

      write(*,*)' Dimension = ',idphon
      write(*,*)' Truncated dimension = ',idphontr

      idphontot=idphon

      if (idphon.gt.0) then 

      open(23,file='AD_J_Pi_Tz.dat',status='unknown',form='formatted')
      write(23,*)ip,jcal,itzcal
      close(23)
      CALL execute_command_line('./run_admat2.sh' )

!  Choleski anal. of spurious subspace
       n_spur=0 ! temporary
       call read_sub_dmat_str(dd,n_spur,idphon)
!       call read_sub_dmat(dd,n_spur)
       write(*,*)
       write(*,*)' Cholesky analysis of spurious subspace'
!       call cholesk_dec(n_spur,dd,nx,mxt)

       call cholesk(1,n_spur,no,noo,dd,nx,mxt)

       open(3,file='chol_spur_subspace.dat',form='formatted',status='unknown')
       write(3,*)n_spur
       do k=1,n_spur
        write(3,*)nx(k)
       enddo
       close(3)

       deallocate(dd,nx,mxt)

       call cholesk(0,idphontr,no,idphontr,dd,nx,mxt)


      call ham_geev(idphon,idphontr,no,nor,ns,irow,wr,xr,vr,ip,jcal,itzcal,mxtr,phonbs,nx,h_corr)
   
      do i=1,no
       write(13)nlam+i,ip,jcal,wr(i),itzcal
      enddo
 
      write(*,*)' X ',idphontr,no
      write(12)ip,jcal,no,idphontot
      write(12)(phonbs(i)%ila,phonbs(i)%ilap,i=1,idphontot)

      mxt=0
      ii=0
      do i=1,idphontr
        if (nx(i).eq.1) then
            ii=ii+1
            mxt(ii)=mxtr(i)
        endif
      enddo


      write(22)ip,jcal,no,no
!      write(22)(phonbs(mxtr(mxt(i)))%ila,phonbs(mxtr(mxt(i)))%ilap,i=1,no)
      write(22)(phonbs(mxt(i))%ila,phonbs(mxt(i))%ilap,i=1,no)


      write(122,*)(mxt(i),i=1,no)
      write(122,*)(mxtr(i),i=1,idphontr)
      write(122,*)(mxtr(mxt(i)),i=1,idphontr)


      write(*,*)' idphontot =',idphontot      
      do j=1,no
       write(12)(xr(i,j),i=1,idphontot)
       write(22)(vr(i,j),i=1,no)
!       write(12)(phonbs(i)%ila,phonbs(i)%ilap,xr(i,j),i=1,idphontot)
!       write(22)(phonbs(mxtr(mxt(i)))%ila,phonbs(mxtr(mxt(i)))%ilap,vr(i,j),i=1,no)
      enddo

      nlam=nlam+no

      write(*,*)'*************'


      deallocate(dd,nx,mxt,wr,xr,vr,h_corr)

      endif

      deallocate(mxtr)

!c      deallocate(phon1,phon2,phonus,phonmus,phonbs)

      write(*,*)'----------------------------------------------'

      CALL execute_command_line('rm ./scratch/a_mat*')
      CALL execute_command_line('rm ./scratch/d_mat*')

      enddo
      enddo  ! cycle over J
      enddo  ! cycle over parity

      write(*,*) ' Number of 2 phonon states ', nlam
      
!c      close(66)
!c      close(33)
      close(99)
      close(12)
      close(22)
    

      end
