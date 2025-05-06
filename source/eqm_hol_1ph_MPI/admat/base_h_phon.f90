module base_h_phon

use types_eqm

contains 

  subroutine odd_phonbase(nf,ipar_bs,ijj_bs,it_bs,dim_phon,lev,phon,isp_min,isp_max,par_lev,dim_base,dim_baser,phonbs,ind_red,n_spur) 
   implicit none 

!   include 'types_eqm_hole.inc'

   integer ::  ipar_bs,ijj_bs,ipar_sp,nf,dim_phon,dim_bs,dim_sp_n,dim_sp_p,jmin,jmax,it_bs,tzi,i_sp,n_spur

   double precision :: en,etrunc,en_CM
   
   type(phonbase_typ), dimension (:), allocatable :: phonbs,phonbs_reor
   type(phon_typ), dimension (:), allocatable :: phon   
   type(level_typ),dimension(:), allocatable :: lev  
   integer :: i,j,ii,isp_min,isp_max,dim_base,dim_baser,ipar,ijj,itt
   integer, dimension (:), allocatable :: ind_red
   integer, dimension(:), allocatable :: par_lev, hol_lev
   integer, dimension(:), allocatable :: mxtr,nxtr

   character*30 name1f
   
    

   if (nf.eq.1) name1f='1phonon/1f_states.dat'

   dim_bs=500000
   allocate(phon(dim_phon))
   allocate(phonbs(dim_bs),phonbs_reor(dim_bs))
   allocate(ind_red(dim_bs))
    

  open(23,file='ethr_1f.dat',status='old',form='formatted')
!   write(*,*)' Energy threshold for 2 phonon states?'
   read(23,*)etrunc
   read(23,*)en_CM
!   write(*,*)etrunc
  close(23) 

   open (3,file=name1f,status='old',form='unformatted')

    do while (.not.eof(3))
      read(3)i,ipar,ijj,en,itt
       phon(i)%par=ipar
       phon(i)%j=ijj
       phon(i)%enf=en
       phon(i)%tz=itt
       phon(i)%spur=0
       if (ipar==-1.and.ijj==1.and.itt==0.and.en <= en_CM) phon(i)%spur=1  ! tag the spurious 1-phonon state
    enddo

   close(3)

   dim_phon=i
!   write(*,*)' Number of phonons = ',dim_phon
!   write(*,*)' Single-particle levels',isp_min,isp_max
 
  phonbs%spur=0 
   
   ii=0
   do i_sp=isp_min,isp_max
    i=par_lev(i_sp)
    ipar_sp=(-1)**lev(i)%l
    tzi=-1 ! neutron hole 
    if (mod(i,2).ne.0) tzi=1 ! proton hole 
    do j=1,dim_phon
      jmin=iabs(2*phon(j)%j-lev(i)%j)
      jmax=2*phon(j)%j+lev(i)%j
      if (phon(j)%par*ipar_sp.eq.ipar_bs) then
      if (mod(ijj_bs,2).ne.0) then
       if (ijj_bs.le.jmax.and.ijj_bs.ge.jmin) then
         if ((2*phon(j)%tz+tzi).eq.it_bs) then 
           ii=ii+1
           phonbs(ii)%ila=j
           phonbs(ii)%isp=i
           phonbs(ii)%tz=2*phon(j)%tz+tzi
           if (phon(j)%spur==1) phonbs(ii)%spur=1
         endif
        endif 
       endif 
      endif
 
    enddo
   enddo
  
  dim_base=ii

 ii=0
 do i=1,dim_base
   if (phonbs(i)%spur.eq.1) then
    ii=ii+1
    phonbs_reor(ii)=phonbs(i)
   endif
 enddo

 n_spur=ii  ! number of selected spurious states

! write(*,*)' Dimension of spurious subspace = ',n_spur

 ii=0
 do i=1,dim_base
   if (phonbs(i)%spur.eq.0) then
    ii=ii+1
    phonbs_reor(ii+n_spur)=phonbs(i)
   endif
 enddo

 phonbs=phonbs_reor

   allocate(mxtr(dim_base),nxtr(dim_base))
   mxtr=0
   nxtr=0



!  write(*,*)' Energy threshold for 1 phonon states?'
!  write(*,*)etrunc
   
   ii=0
   do i=1,dim_base
     if (phonbs(i)%tz.eq.it_bs) then 
      if (lev(phonbs(i)%isp)%tz.eq.-1*it_bs) then 
        if ((phonbs(i)%spur.eq.1).or.(phonbs(i)%spur.eq.0.and.phon(phonbs(i)%ila)%enf.lt.etrunc)) then
          ii=ii+1
          ind_red(ii)=i
          mxtr(ii)=i
          nxtr(i)=ii
       endif
      endif
     endif
   enddo

   dim_baser=ii

!   open(1,file=' phonon_base.dat',status='unknown',form='formatted')
!   write(1,*)'i   2*Tz    i_sp      en_sp      lam      en_phon'
!   write(1,*)'i   spur   2*Tz   i_sp   en_sp    lam   en_phon'
!    do i=1,dim_base
!    write(1,'(4i5,f10.5,i5,f10.5)')i,phonbs(i)%spur,phonbs(i)%tz,phonbs(i)%isp,lev(phonbs(i)%isp)%en,phonbs(i)%ila,phon(phonbs(i)%ila)%enf
!      write(1,'(3i10,f10.5,i5,f10.5)')i,phonbs(i)%tz,phonbs(i)%isp,lev(phonbs(i)%isp)%en,phonbs(i)%ila,phon(phonbs(i)%ila)%enf
!    enddo

!   write(1,*)'spurious basis  states'
!   do i=1,dim_base
!    if (phonbs(i)%spur ==1) write(1,'(4i5,f10.5,i5,f10.5)')i,phonbs(i)%spur,phonbs(i)%tz,phonbs(i)%isp,lev(phonbs(i)%isp)%en,phonbs(i)%ila,phon(phonbs(i)%ila)%enf
!   enddo


!  list of 2 ph densitities to be calculated    
!    open(891,file='2_phon_dens_list.dat',status='unknown',form='formatted',position='append')
  
!   write(1,*)'reduced basis  '

!   do ii=1,dim_baser
!    i=ind_red(ii)
!    write(1,'(3i10,f10.5,i5,f10.5)')i,phonbs(i)%tz,phonbs(i)%isp,lev(phonbs(i)%isp)%en,phonbs(i)%ila,phon(phonbs(i)%ila)%enf
!    write(1,'(4i5,f10.5,i5,f10.5)')i,phonbs(i)%spur,phonbs(i)%tz,phonbs(i)%isp,lev(phonbs(i)%isp)%en,phonbs(i)%ila,phon(phonbs(i)%ila)%enf
!    write(891,*)phonbs(i)%ila
!   enddo
!   close(891)

!   close(1)  

!  open(2,file='mxtr.dat',status='unknown',form='unformatted')
!  write(2)dim_baser
!  write(2)(mxtr(i),i=1,dim_baser)
!  close(2)

!  open(2,file='nxtr.dat',status='unknown',form='unformatted')
!  write(2)dim_base
!  write(2)(nxtr(i),i=1,dim_base)
!  close(2)
    
   
 
end subroutine odd_phonbase

end  module base_h_phon
