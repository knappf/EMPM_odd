module Tcm_2body
  use technical
  use math
  use geom

  contains

  double precision function T_nn(i,j,k,l,Jp)
    integer :: i,j,k,l,Jp
    integer :: j2i,j2j,j2k,j2l
    
    T_nn = 0.d0
    j2i=levn(i)%j2
    j2j=levn(j)%j2
    j2k=levn(k)%j2
    j2l=levn(l)%j2
    sixj=sixj_real(0.5d0*j2i,0.5d0*j2j,dfloat(Jp),0.5d0*j2l,0.5d0*j2k,1.0d0)
    T_nn= dfloat((-1)**((j2j+j2k)/2+Jp))*sixj*grad_red_me(i,k)*grad_red_me(j,l)
  end function T_nn

  double precision function T_nn_asym(i,j,k,l,Jp)
    integer :: i,j,k,l,Jp
    integer :: j2k,j2l
    double precision :: fact

    T_nn_asym = 0.d0

    j2k=levn(k)%j2
    j2l=levn(l)%j2

    fact=1.d0/dsqrt((1.d0+delta1(i,j))*(1.d0+delta1(k,l)))
    T_nn_asym = fact*(T_nn(i,j,k,l,Jp)-dfloat((-1)**((j2k+j2l)/2-Jp))*T_nn(i,j,l,k,Jp))

  end function T_nn_asym

  double precision function grad_red_me(i,j)

    integer :: i,j,j2i,j2j,li,lj
    double precision :: sixj,phase,ji,jj

    grad_red_me = 0.d0

    j2i=levn(i)%j2
    j2j=levn(j)%j2
    li=levn(i)%l
    lj=levn(j)%l
    ni=levn(i)%nn
    nj=levn(j)%nn
    sixj=sixj_real(0.5d0*j2i,0.5d0*j2j,1.d0,dfloat(lj),dfloat(li),0.5d0)
    fact=dfloat((-1)**(li+(j2j+1)/2))
    fact=fact*dsqrt((j2i+1.d0)*(j2j+1.d0))
    grad_red_me = fact*sixj*(dsqrt((lj+1.d0)*(nj+lj+3.d0/2.d0))*delta1(li,lj+1)*delta1(ni,nj)+dsqrt((lj+1.d0)*nj)*delta1(li,lj+1)*delta1(ni,nj-1)+dsqrt(lj*(nj+lj+1.d0/2.d0))*delta1(li,lj-1)*delta1(ni,nj)+dsqrt(lj*(nj+1.d0))*delta1(li,lj-1)*delta1(ni,nj+1))

  end function grad_red_me


end module Tcm_2body
