module read_admat

    use types_eqm
    contains

    subroutine read_sub_dmat_str(dmatr,ndim,ndimt)

        double precision, dimension(:,:), allocatable :: dmatr
        double precision, dimension(:), allocatable :: dmm
        character (len=15) :: row_num  

        if (.not.allocated(dmatr))  allocate(dmatr(ndim,ndim))
        dmatr=0.d0
  
        allocate(dmm(ndimt))
        dmm=0.0d0

        do iii=1,ndim
          write(row_num,'(i15.15)')iii
          open(66,file='./scratch/d_mat_'//row_num,status='old',form='unformatted')
           read(66)(dmm(jjj),jjj=1,ndimt)
           dmatr(iii,1:ndim)=dmm(1:ndim)
          close(66) 
        enddo
  
        deallocate(dmm)

        return

    end subroutine read_sub_dmat_str

    subroutine read_dmat(dmatr)

        implicit double precision (a-h,o-z)
  
  !c      include 'formats_eqm.inc'
  
        double precision, dimension(:,:), allocatable :: dmatr
        double precision, dimension(:), allocatable :: dmm
  
        integer, dimension(:), allocatable :: nxtr
        character (len=15) :: row_num 
  !c      double precision, dimension(:), allocatable :: work,e
  
    
  
        open(2,file='nxtr.dat',status='unknown',form='unformatted')
        read(2)ndim_tot
        allocate(nxtr(ndim_tot))
        nxtr=0
        read(2)(nxtr(i),i=1,ndim_tot)
        close(2)
  
        write(991,*)'ndim_tot =', ndim_tot
        write(991,*)(nxtr(i),i=1,ndim_tot)

        ndim_tr=0
        do i=1,ndim_tot
            if (nxtr(i).ne.0) ndim_tr=ndim_tr+1
        enddo

!        allocate(dmatr(ndim_tr,ndim_tot))
        dmatr=0.d0
  
  !c      iold=-100
  
   
!        read(66)ndimrt,ndimt

        allocate(dmm(ndim_tot))
        dmm=0.0d0
  
  !      if (ndimrt.ne.ndimr.or.ndimt.ne.ndim) then
  !        write(*,*)' Dimensions does not match in d_m file'
  !        stop
  !      endif


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,row_num,ifile,dmm)
!$OMP DO 
        do i=1,ndim_tr
         write(row_num,'(i15.15)')i
         ifile=800+i
         open(ifile,file='./scratch/d_mat_'//row_num,status='old',form='unformatted') 
         read(ifile)(dmm(j),j=1,ndim_tot)
         do j=1,ndim_tot
           if (nxtr(j).ne.0) then
            dmatr(i,nxtr(j))=dmm(j)
           endif
         enddo
         close(ifile)
        enddo
!$OMP END DO
!$OMP END PARALLEL  
 
        deallocate(dmm)
  
  !      do while (.not.eof(66))
  !      read(66)i,j,dd
  !      if (nxtr(j).ne.0) then 
  !      dmatr(i,nxtr(j))=dd
  !      endif
  
  !      enddo
  
  !      close(66)
  
        deallocate(nxtr)
        
        end subroutine read_dmat

        subroutine read_admatr(path,matr,ndim_tr,ndim_tot)

            double precision, dimension(:,:), allocatable :: matr
            double precision, dimension(:), allocatable :: dmm
            character (len=15) :: row_num
            character (len=16) :: path   
    
            if (.not.allocated(matr))  allocate(matr(ndim_tr,ndim_tot))
            matr=0.d0
      
            allocate(dmm(ndim_tot))
            dmm=0.0d0
    
            do i=1,ndim_tr
              write(932,*)i,ndim_tr
              write(row_num,'(i15.15)')i
!              open(66,file='./scratch/d_mat_'//row_num,status='old',form='unformatted')
              open(66,file=path//row_num,status='old',form='unformatted')
               read(66)(dmm(j),j=1,ndim_tot)
               matr(i,1:ndim_tot)=dmm(1:ndim_tot)
              close(66) 
            enddo
      
            deallocate(dmm)
    
            return

        end subroutine read_admatr

        subroutine read_admatr_OMP(path,matr,ndim_tr,ndim_tot)

            double precision, dimension(:,:), allocatable :: matr
            double precision, dimension(:), allocatable :: dmm
            character (len=15) :: row_num
            character (len=16) :: path   
    
            if (.not.allocated(matr))  allocate(matr(ndim_tr,ndim_tot))
            matr=0.d0
      
!            allocate(dmm(ndim_tot))
!            dmm=0.0d0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,row_num,ifile)
!$OMP DO 
            do i=1,ndim_tr
!              write(932,*)i,ndim_tr
              write(row_num,'(i15.15)')i
              ifile=800+i
              open(ifile,file=path//row_num,status='old',form='unformatted')
!               read(ifile)(dmm(j),j=1,ndim_tot)
               read(ifile)(matr(i,j),j=1,ndim_tot)
!               matr(i,1:ndim_tot)=dmm(1:ndim_tot)
              close(ifile) 
            enddo
!$OMP END DO
!$OMP END PARALLEL      
      
!            deallocate(dmm)
    
            return

        end subroutine read_admatr_OMP



end module read_admat
