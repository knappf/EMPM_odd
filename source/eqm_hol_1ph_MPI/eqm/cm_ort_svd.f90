module cm_ort_svd

use types_eqm

contains 

subroutine dmat_spur_set(dim_base,dim_baser,dim_ind,phonbs,d_mat,ind_red,nx,s,vt,ns)

  implicit double precision (a-h,o-z)
  
!  include 'types_eqm_hole.inc'
!  include 'types_eqm.inc'
  
  type(phonbase_typ), dimension (:), allocatable :: phonbs
  integer :: dim_base,dim_baser,dim_spur,dim_spur_all,dim_ind,i,j,k,ii,jj,ns,n_spur
  integer, dimension (:), allocatable :: ind_spur,ind_spur_all 
  double precision, dimension (:,:),allocatable :: d_mat,d_mat_sp,d_mat_sp_full
  double precision, dimension (:), allocatable :: s,work
  double precision, dimension (:,:), allocatable :: u,vt,a
  integer, dimension (:), allocatable :: ind_red,nx,nx_spur
  integer :: lwork,info
  
  double precision :: xpom
  
  allocate(ind_spur(dim_base))
  ind_spur=0
  allocate(ind_spur_all(dim_base))
  ind_spur_all=0

  open(3,file='chol_spur_subspace.dat',form='formatted',status='old')
  read(3,*)n_spur
  allocate(nx_spur(n_spur))
  nx_spur=0
  do k=1,n_spur 
   read(3,*)nx_spur(k)
  enddo
  close(3)
  
  ii=0
  do i=1,dim_base
    if (phonbs(i)%spur.eq.1.and.nx_spur(i).eq.1) then 
      ii=ii+1
      ind_spur(ii)=i
    endif 
  enddo
  
  dim_spur=ii
 
  ii=0
  do i=1,dim_base
    if (phonbs(i)%spur.eq.1) then 
      ii=ii+1
      ind_spur_all(ii)=i
    endif 
  enddo

  dim_spur_all=ii

  write(*,*)' Dims of spur. subspace  (lin. ind. total) '
  write(*,*)dim_spur,dim_spur_all
  
  if (dim_spur.eq.0) then 
    ns=0
    allocate(vt(dim_ind,dim_ind))
    vt=0.d0 
    do i=1,dim_ind
      vt(i,i)=1.0d0
    enddo
  return 
  endif 
  
  
  allocate(d_mat_sp(dim_spur,dim_ind))
  d_mat_sp=0.d0
  
  do i=1,dim_spur
   jj=0
   do j=1,dim_baser
    if (nx(j).eq.1) then 
     jj=jj+1
     d_mat_sp(i,jj)=d_mat(j,ind_spur(i))
    endif
    enddo
  enddo
  
  write(*,*)' Independent rows  ', jj

  allocate(d_mat_sp_full(dim_spur_all,dim_ind))
  d_mat_sp_full=0.d0
  
  do i=1,dim_spur_all
   jj=0
   do j=1,dim_baser
    if (nx(j).eq.1) then 
     jj=jj+1
     d_mat_sp_full(i,jj)=d_mat(j,ind_spur_all(i))
    endif
    enddo
  enddo
  
  !    SVD
  allocate(s(dim_spur))
  allocate(u(dim_spur,dim_spur),vt(dim_ind,dim_ind))
  lwork=10*dim_baser
  allocate(work(lwork))
  
  call dgesvd('A','A', dim_spur, dim_ind, d_mat_sp, dim_spur, s, u, dim_spur, vt, dim_ind, work, lwork, info)
  
  
  write(*,*) 'SVD  info =',info
  
  write(998,*)
  write(998,*)' SVD '
  write(998,*)(s(i),i=1,dim_spur)
  
  ii=0
  do i=1,dim_spur
  if (dabs(s(i)).gt.1.d-4) then
   ii=ii+1
  endif
  enddo
  ns=ii
  write(*,*)' Number of nonzero dg. values ',ns
  if (dim_spur.ne.ns) write(*,*)'***** Number of nonzero SVD egv different from the dimension of sp. subspace !******'
  
  
  allocate(a(dim_spur_all,dim_ind))
  a=0.0d0
  
  call dgemm('N','T',dim_spur_all,dim_ind,dim_ind,1.d0,d_mat_sp_full,dim_spur_all,vt,dim_ind,0.d0,a,dim_spur_all)
  
  write(*,*)'Test SVD'
  do i=1,dim_spur_all
  ! write(999,'(1000f10.5)')(a(i,j),j=1,dim_ind)
   do j=ns+1,dim_ind
     if (dabs(a(i,j)).gt.1.d-7) write(*,*)'Nonzero elem <spur| bas> !',i,j,a(i,j)
   enddo
  enddo
  

!  call read_sub_dmat(d_full,dim_base)

!  call dgemm('N','T',dim_spur,dim_base,dim_base,1.d0,d_full,dim_spur,vt,dim_ind,0.d0,a,dim_spur)
  
  !write(999,*)
  !write(999,*)' SVD '
  !write(999,*)(s(i),i=1,dim_spur)
  !do i=1,dim_ind
  ! write(999,'(1000f10.5)')(vt(i,j),j=1,dim_ind)
  !enddo
  
  
  ! nulify 
  
  !do i=1,ns
  !  vt(i,:)=0.d0
  !enddo
  
  !do i=1,ns
  ! vt(i,i)=1.0
  !enddo
  
  !ii=0
  !jj=0
  !do i=1,dim_baser
  !  if (nx(i).eq.1) then 
  !    ii=ind_red(i)
  !  if (phonbs(ii)%spur.eq.1) then 
  !    jj=jj+1
  !    vt(iii,i)=1.0d0
  !  endif
  ! endif 
  !enddo
  
  !write(*,*)' Number of spurious states check ',jj
  
  
  
  return
  
  end subroutine dmat_spur_set
!***********************************************************************
  

end module cm_ort_svd
