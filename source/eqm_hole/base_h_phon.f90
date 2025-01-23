module base_h_phon

contains 

  subroutine odd_phonbase(nf,ipar_bs,ijj_bs,it_bs,dim_phon,lev,phon,isp_min,isp_max,hol_lev,dim_base,dim_baser,phonbs,ind_red) 
   implicit none 

   include 'types_eqm_hole.inc'

   integer ::  ipar_bs,ijj_bs,ipar_sp,nf,dim_phon,dim_bs,dim_sp_n,dim_sp_p,jmin,jmax,it_bs,tzi,i_sp

   double precision :: en 
   
   type(phonbase_typ), dimension (:), allocatable :: phonbs
   type(phon_typ), dimension (:), allocatable :: phon   
   type(level_typ),dimension(:), allocatable :: lev  
   integer :: i,j,ii,isp_min,isp_max,dim_base,dim_baser,ipar,ijj,itt
   integer, dimension (:), allocatable :: ind_red
   integer, dimension(:), allocatable :: par_lev, hol_lev

   character*30 name1f
   
    

   if (nf.eq.1) name1f='1phonon/1f_states.dat'

   dim_bs=30000
   allocate(phon(dim_phon))
   allocate(phonbs(dim_bs))
   allocate(ind_red(dim_bs))
    


   open (3,file=name1f,status='old',form='unformatted')

    do while (.not.eof(3))
      read(3)i,ipar,ijj,en,itt
       phon(i)%par=ipar
       phon(i)%j=ijj
       phon(i)%enf=en
       phon(i)%tz=itt
    enddo

   close(3)

   dim_phon=i
   write(*,*)' Number of phonons = ',dim_phon
   write(*,*)' Single-particle levels',isp_min,isp_max
 
   
   ii=0
   do i_sp=isp_min,isp_max
    i=hol_lev(i_sp)
    ipar_sp=(-1)**lev(i)%l
    tzi=1
    if (mod(i,2).ne.0) tzi=-1
    do j=1,dim_phon
      jmin=iabs(2*phon(j)%j-lev(i)%j)
      jmax=2*phon(j)%j+lev(i)%j
      if (phon(j)%par*ipar_sp.eq.ipar_bs) then
      if (mod(ijj_bs,2).ne.0) then
       if (ijj_bs.le.jmax.and.ijj_bs.ge.jmin) then
         if ((2*phon(j)%tz-tzi).eq.it_bs) then 
!           if (dabs(phon(j)%enf).gt.0.9d0) then 
           ii=ii+1
           phonbs(ii)%ila=j
           phonbs(ii)%isp=i
           phonbs(ii)%tz=2*phon(j)%tz-tzi
!        endif
         endif
        endif 
       endif 
      endif
 
    enddo
   enddo
  
   dim_base=ii

   ii=0
   do i=1,dim_base
     if (phonbs(i)%tz.eq.it_bs) then 
      if (lev(phonbs(i)%isp)%tz.eq.-1*it_bs) then 
        if (phon(phonbs(i)%ila)%enf.gt.0.9d0) then

     ii=ii+1
     ind_red(ii)=i
       endif
      endif
     endif
   enddo

   dim_baser=ii

   open(1,file=' phonon_base.dat',status='unknown',form='formatted')
   write(1,*)'i   2*Tz    i_sp      en_sp      lam      en_phon'
    do i=1,dim_base
     write(1,'(3i5,f10.5,i5,f10.5)')i,phonbs(i)%tz,phonbs(i)%isp,lev(phonbs(i)%isp)%en,phonbs(i)%ila,phon(phonbs(i)%ila)%enf
    enddo


    write(1,*)'reduced basis  '

    do ii=1,dim_baser
    i=ind_red(ii)
     write(1,'(3i5,f10.5,i5,f10.5)')i,phonbs(i)%tz,phonbs(i)%isp,lev(phonbs(i)%isp)%en,phonbs(i)%ila,phon(phonbs(i)%ila)%enf
    enddo


    close(1)  
 
   
 
end subroutine odd_phonbase

end  module base_h_phon
