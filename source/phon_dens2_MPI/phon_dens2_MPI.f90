!*     
!*     Program phon_dens computes phonon densities in proton-neutron 
!*     J-coupled formalism.

!*     last update 12.7.2010

      program phon_dens

      use rdens
      use input_sp
!c      use rdens_al

      implicit double precision (a-h,o-z)

!      include 'input_phon_dens.inc'

      include 'mpif.h'
      include 'formats_phon_dens.inc'
      include 'types_phon_dens.inc'

      type(level_typ),dimension(:), allocatable :: lev

      type (phon_typ), dimension (:,:), allocatable :: phonbs
      integer, dimension (:), allocatable :: nphon,iphous,iphous2
      integer, dimension(:), allocatable :: par_lev, hol_lev
      integer, dimension(:,:), allocatable :: i_ph_int

      character(len=16) fname

      integer :: myid,ierr,numprocs
      
!*     loading of input data 

call MPI_INIT( ierr )
call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )


! for debugging
!numprocs=1
!myid=0

xrotrunc=1.d-10
      
write(*,*)'Loading of input '


allocate(i_ph_int(4,-1:1))

  open(1,file='input_tda_coup.dat',status='old',form='formatted')
    read(1,'(30i6)')ia,iz
    read(1,'(30i6)')i_ph_int(1,1),i_ph_int(2,1)
    read(1,'(30i6)')i_ph_int(1,-1),i_ph_int(2,-1)
    read(1,'(30i6)')i_ph_int(3,1),i_ph_int(4,1)
    read(1,'(30i6)')i_ph_int(3,-1),i_ph_int(4,-1)
  close(1)


 idim_sp=1000 ! dimension of single particle space

 allocate(lev(idim_sp))

 call inp_sp(lev,nlev,jmaxn)

     allocate(par_lev(nlev),hol_lev(nlev))
     par_lev=0
     hol_lev=0

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


      write(*,*)'2* jmax s.p  =',jmaxn
      
      jmax=0


      open (3,file='1phonon/1f_states.dat',status='old',form='unformatted')

      do while (.not.eof(3))
       read(3)i,ipar,ijj,en
        if (ijj > jmax) jmax=ijj
      enddo

      close(3)

      write(*,*) ' jmax 1phon = ', jmax
 
      ifmx=i ! treba upravit aby sa to nacitalo 
!      jmax=jmaxn

      jamax=2*jmax
      isimax=jmaxn


      call loadphon(phonbs,nphon)

      if (myid.eq.0) then 
       write(*,*)' Number of 1phonon states ',nphon(1)
       write(*,*)' Number of 2phonon states ',nphon(2)
      endif 

      nff=2

      call selphon(nff,phonbs,nphon,iphous,ns1,myid)
      call selphon2(nff,phonbs,nphon,iphous2,ns2,myid)

      if (ns1 > 0) then 
      call r_value(lev,ip_lev,par_lev,ih_lev,hol_lev,jamax,1,nlev,isimax,iphous,iphous2,phonbs,nphon,ns1,ns2,myid,numprocs)
 
      endif 

call MPI_BARRIER(  MPI_COMM_WORLD, ierr)

      
call MPI_FINALIZE(ierr)


      end 
!*     
!*     END of the main program 
!* 

      




      
      
