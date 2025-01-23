!*     Program tda_pn_coup.f computes TDA phonons in  J-coupled proton-neutron 
!*     formalism.
!*
!*     last update 4.2.2015 

      program tda_pn_coup
      
     

 
      use cmconst
      use input_sp

      implicit double precision (a-h,o-z)

 
      include 'formats_tda_cp.inc' !Definiton of used formats 
      include 'types_tda_cp.inc'      
      include 'input_tda_cp.inc'


      integer, dimension(:), allocatable :: ipozi,isosp
      integer, dimension(:,:), allocatable :: i_ph_int

      double precision, dimension (:), allocatable :: work,e 
      double precision, dimension(:,:), allocatable :: amtr,tbase,hami,amtrold 

      double precision, dimension(:,:,:,:,:), allocatable :: f_int

      type(level_typ),dimension(:), allocatable :: lev
      type(ph_typ),dimension(:), allocatable :: ph_base

      character*1, dimension(0:40) :: orbit 
      character(len=30)file_fp,file_fn,file_fpn,file_fpcm,file_fncm,file_fpncm
      character*5 tnf
      character*1 tcm

      integer(kind=1) :: j_f
      integer(kind=2) :: i_i,i_j,i_k,i_l

        
      xrotrunc=0.000000000001d0


     
      orbit(0)='s'
      orbit(1)='p'
      orbit(2)='d'
      orbit(3)='f'
      orbit(4)='g'
      orbit(5)='h'
      orbit(6)='i'
      orbit(7)='j'
      orbit(8)='k'
      orbit(9)='l'
      orbit(10)='m'
      orbit(11)='n'
      orbit(13)='o'

      open(881,file='tda_r_overl.dat',status='unknown',form='formatted')
      open(77,file='tda_pn_cp.log',status='unknown',form='formatted') 

      open(94,file='phon_struct.dat',status='unknown',form='formatted')
      open(96,file='sp_levord.dat',status='unknown',form='formatted')
    
      open(743,file='1phonon/1ph_cont_cp.dat',status='unknown',form='formatted')
                
!*     loading of input data 
      allocate(i_ph_int(4,-1:1))     
 
      write(*,*)'Loading of input '

      open(1,file='input_tda_coup.dat',status='old',form='formatted')      
      
      read(1,15)ia,iz
      read(1,15)i_ph_int(1,1),i_ph_int(2,1)
      read(1,15)i_ph_int(1,-1),i_ph_int(2,-1)
      read(1,15)i_ph_int(3,1),i_ph_int(4,1)
      read(1,15)i_ph_int(3,-1),i_ph_int(4,-1)
      read(1,'(a1)')tcm
      close(1)

      write(*,*) 'CM ort? ',tcm

      nlev=300
      allocate(lev(nlev))
      allocate(ph_base(nlev*nlev))
      allocate(isosp(nlev*nlev))
      

      call inp_sp(lev,nlev,jmax)
      


      open(33,file='1phonon/1f_states.dat',status='unknown',form='unformatted')


      open(2,file='1phonon/1f_cph.dat',status='unknown',form='unformatted')


      call read_fmat(f_int,nlev,jmax)

      nlam=0


      do ipar=-1,1,2   ! cycle over parity
      do ijj=0,2*jmax        ! cycle over J


      call phbase(ipar,ijj,i_ph_int,lev,ph_base,idim_ph)


       write(*,*)'Parity =   ',ipar
       write(*,*)'     J =   ',ijj
       write(*,*)'dimension  = ',idim_ph
                

       if (idim_ph.ne.0) then  


       allocate(amtr(idim_ph,idim_ph))
       amtr=0.d0

        do i=1,idim_ph
          ip=ph_base(i)%par
          ih=ph_base(i)%hol
          do j=1,idim_ph
            ipp=ph_base(j)%par
            ihp=ph_base(j)%hol
 
             amat=0.d0
      
              if (ih.eq.ihp.and.ip.eq.ipp ) amat=amat+lev(ip)%en-lev(ih)%en
      
              ifz=(-1)**((lev(ipp)%j+lev(ihp)%j)/2+ijj)
              amat=amat-f_int(ip,ih,ihp,ipp,ijj)*dfloat(ifz)
      
              amtr(i,j)=amat
            enddo
          enddo

       
          do i=1,idim_ph
            do j=1,idim_ph
              if (dabs(amtr(i,j)-amtr(j,i)).gt.1.d-10) write(*,*)'nsym   ',i,j,amtr(i,j),amtr(j,i)
            enddo
          enddo


          ndmx=idim_ph

!c   control output of A matrix 
          write(999,*)'--------------------------------------'
          do iii=1,ndmx
           write(999,'(100f15.5)')(amtr(iii,jjj),jjj=1,ndmx)
          enddo 
          write(999,*)'--------------------------------------'
     

          lwork=26*ndmx 
          
          allocate(work(26*ndmx))
          allocate(e(ndmx))

           
          write(*,*)' ndmx = ',ndmx 
!  CM orthogonalisation

if (ipar.eq.-1.and.ijj.eq.1.and.tcm.eq.'y') then 

          allocate(hami(ndmx,ndmx))
          hami=0.0d0
          
          allocate(amtrold(ndmx,ndmx))
          amtrold=0.d0
          amtrold=amtr

          call cmcreate(idim_ph,ph_base,tbase)

          call ortog(ndmx,tbase)
          
!  control output od tbase          
          write(777,*)'tbase'          
          write(777,*)
          do i=1,ndmx
          write(777,'(1000f10.5)')(tbase(i,j),j=1,ndmx)
          enddo


          call dgemm('N','N',ndmx,ndmx,ndmx,1.0d0,amtr,ndmx,tbase,ndmx,0.0d0,hami,ndmx)
!          call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
!           call gemm(amtr,tbase,hami,'N','N',1.0d0,0.0d0)

          amtr=0.0d0

          call dgemm('T','N',ndmx,ndmx,ndmx,1.0d0,tbase,ndmx,hami,ndmx,0.0d0,amtr,ndmx)


 321      format(300f15.10)
          write(777,*)' A matica'
          do i=1,ndmx
            write(777,321)(amtr(i,j),j=1,ndmx)
          enddo

!   R decouple
          do i=2,ndmx
           amtr(i,1)=0.0d0
           amtr(1,i)=0.0d0
          enddo


          do j=1,ndmx
            do i=1,j-1
              amtr(i,j)=0.d0
            enddo
          enddo

          amtr(1,1)=10000000.0d0

          write(777,*)' A matica ort. '
          do i=1,ndmx
            write(777,321)(amtr(i,j),j=1,ndmx)
          enddo


endif 

         
 
          call DSYEV('V','l',ndmx,amtr,ndmx,e,WORK,LWORK,INFO )

          write(777,'(100f15.5)')          
          write(777,*)'vl. vektory A'
          do i=1,ndmx
            write(777,'(100f15.5)')(amtr(i,j),j=1,ndmx)
          enddo

          do j=1,ndmx
            do i=1,ndmx
             if (dabs(amtr(i,j)).gt.1.d-10) isosp(j)=ph_base(i)%tz
            enddo

           do i=1,ndmx
             if (dabs(amtr(i,j)).gt.1.d-10.and.ph_base(i)%tz.ne.isosp(j)) write(*,*)' eigenvector Tz mixing !!!!!!!'
            enddo
          enddo
          
          
          write(77,*)'  Parity = ',ipar,'    J = ',ijj
          write(77,792)(e(ii),ii=1,ndmx)
          write(77,*)


      do i=1,ndmx
        write(tnf,101)nlam+i
        write(33)nlam+i,ipar,ijj,e(i),isosp(i)
        write(94,*)
        write(94,*)'________________________________________'
        write(94,*)
        write(94,*)'i     Parity    J        energy       Tz'
        write(94,'(3i5,f10.5,i4)')nlam+i,ipar,ijj,e(i),isosp(i)
        write(94,*)
!c        open(2,file='1phonon/1f_cp.'//tnf//''
!c     *,status='unknown',form='unformatted')


        write(2)ipar,ijj,ndmx


        write(2)(ph_base(j)%par,ph_base(j)%hol,amtr(j,i),j=1,ndmx)
!c        close(2)

        write(94,*)'-  p(h)^-1 -'

        do j=1,ndmx
         if (dabs(amtr(j,i)).gt.1.d-7) write(94,'(2i2,5x,f10.5)')ph_base(j)%par,ph_base(j)%hol,amtr(j,i)
        enddo
      enddo
                             

          deallocate(amtr,work,e)

         
          nlam=nlam+ndmx
          
       endif
          
          
          
        enddo
      enddo

      close(33)
      close(2)
      close(4)

!      close(11)
!      close(12)
!      close(13)
!      close(14)
      close(99)

      write(*,*)'dimension of 1 phonon space  ',  nlam




      
      end 
!*     
!*     END of the main program 
!* 
      


      
      

      
