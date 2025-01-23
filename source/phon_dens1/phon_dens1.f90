!     
!     Program phon_dens1 computes 1phonon densities in proton-neutron 
!     J-coupled formalism.

!     last update 12.4.2013

      program phon_dens1
 
      use input_sp 
      use dens_pn


      implicit double precision (a-h,o-z)

      include 'input_phon_dens1.inc'
      include 'formats_phon_dens1.inc'
      include 'types_phon_dens1.inc'

      type(level_typ),dimension(:), allocatable :: lev
      integer, dimension(:,:), allocatable :: i_ph_int
      integer, dimension(:), allocatable :: par_lev, hol_lev

      
!*     loading of input data 

      xrotrunc=1.d-8
      allocate(i_ph_int(4,-1:1))

      
      write(*,*)'Loading of input '

      open(1,file='input_tda_coup.dat',status='old',form='formatted')

      read(1,15)ia,iz
      read(1,15)i_ph_int(1,1),i_ph_int(2,1)
      read(1,15)i_ph_int(1,-1),i_ph_int(2,-1)
      read(1,15)i_ph_int(3,1),i_ph_int(4,1)
      read(1,15)i_ph_int(3,-1),i_ph_int(4,-1)

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

   

      write(*,*)'jmax =',jmax


      open (3,file='1phonon/1f_states.dat',status='old',form='unformatted')

      do while (.not.eof(3))
       read(3)i,ipar,ijj,en
!       jphon(i)=ijj
      enddo

      close(3)



      ifmx=i ! treba upravit aby sa to nacitalo 

!      write(*,*)' Energy threshodl for 1 phonon states?'
!      read(*,*)xthres
!      write(*,*)xthres  
      xthres=100000000.0
      call rop(ifmx,nlev,jmax,lev,xthres,ih_lev,ip_lev,par_lev,hol_lev)
      call roh(ifmx,nlev,jmax,lev,xthres,ih_lev,ip_lev,par_lev,hol_lev)
      
      end 
!     
!     END of the main program 
! 

