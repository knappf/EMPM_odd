!     last modification 21.6.2010      
      
      program eqm 

      use phonon_base
      use admatr
!      use mpi_f08

      implicit double precision (a-h,o-z)

      include 'types_eqm.inc'

      include 'mpif.h'

      type(phonbase_typ), dimension (:), allocatable :: phonbs
      type(phon_typ), dimension (:), allocatable :: phon1,phon2
      integer, dimension (:), allocatable :: phonus, phonmus

!     choleski arrays
      double precision, dimension(:,:), allocatable :: dd,cq,d1,cdu,xr,vr
      double precision, dimension(:), allocatable :: wr
      
      integer, dimension(:), allocatable :: nx,mxt,irow,mxtr, i_resh
      
      character*30 namex,names,namec

      integer :: myid,ierr,numprocs

      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )


      nf=2
      idim1=10000
      idim2=50000
      idimbs=2000000

      open(23,file='AD_J_Pi_Tz.dat',status='old',form='formatted')
      read(23,*)ipcal,jcal,itzcal
      close(23)

      if (myid.eq.0) then
        write(*,*)
        write(*,*)'----------------------------------------------'
        write(*,*)' Parity = ',ipar,'   J = ',jcal, ' Tz = ', itzcal
      endif
      
      call phonbase(nf,ipcal,jcal,itzcal,phonus,phonmus,idim1,idim2,idimbs,idphon,idphontr,phonbs,phon1,phon2,mxtr,myid)

      if (myid.eq.0) then
          write(*,*)' Dimension = ',idphon
          write(*,*)' Truncated dimension = ',idphontr
      endif 

      call admat(nf,ipcal,jcal,phonbs,idphon,no,phon1,phon2,mxt,myid,numprocs)

      call MPI_FINALIZE(ierr)
      
      end
