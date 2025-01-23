module v3body_no2b

use technical

contains

real function V3Bno(i,j)
implicit none
integer :: i,j,ii,jj,it
integer(kind=8) :: iv3b,arrl
integer(kind=8) :: low,high,middle,rang
integer(kind=8) :: ipoz,ipozjt,pair_ij
real :: diff,var

V3Bno=0.0

arrl=size(V3B_ar,DIM=1,KIND=8)
low=1
high=arrl
rang=high-low
middle=(low+high)/2

!This is only for testing
!diff=0.0

if (lev3(i)%jab.ne.lev3(j)%jab) return

ipoz=0
ipozjt=0
!do it=1,lev3(j)%TT-2,2
!do jj=0,lev3(j)%jab-1
! ipoz=ipoz+(idim3jt(jj,it)*(idim3jt(jj,it)+1)/2)
! ipozjt=ipozjt+idim3jt(jj,it)
!enddo
!enddo


if (lev3(j)%TT==1) then  
 do jj=0,lev3(j)%jab-1
   ipoz=ipoz+(idim3jt(jj,1)*(idim3jt(jj,1)+1)/2)
   ipozjt=ipozjt+idim3jt(jj,1)
 enddo
endif 

if (lev3(j)%TT==3) then  

!         do jj=0,jmax
!          ipoz=ipoz+(idim3jt(jj,1)*(idim3jt(jj,1)+1)/2)
!          ipozjt=ipozjt+idim3jt(jj,1)
!         enddo
  ipoz=ipozt1
  ipozjt=ipozjt1
  do jj=0,lev3(j)%jab-1
    ipoz=ipoz+(idim3jt(jj,3)*(idim3jt(jj,3)+1)/2)
    ipozjt=ipozjt+idim3jt(jj,3)
  enddo
endif 

ii=lev3(i)%ord-ipozjt
jj=lev3(j)%ord-ipozjt

if ((ii.ne.0).and.(jj.ne.0)) then
  if (i >= j ) then
!     diff=V3B_ar2((int(ii,8)*(int(ii,8)-1))/2+int(jj,8)+ipoz)
     pair_ij=(int(ii,8)*(int(ii,8)-1))/2+int(jj,8)+ipoz

     do while((V3B_pair(middle).ne.pair_ij).and.(rang.gt.0))
        if(pair_ij>V3B_pair(middle)) then
           low=middle+1
        else
           high=middle-1
        end if
        
        rang=high-low
        middle=(high+low)/2
     end do
     if(V3B_pair(middle)==pair_ij) then
       V3Bno=V3B_ar(middle)
     end if

  endif
endif 

!if (abs(diff-V3Bno)>1.0e-6) then
!  write(*,*)'Diff at',i,j,V3Bno,diff
!endif 


if (i<j) then 
  write(*,*)'i<j in V3B !'
  stop
endif 

return
end function 

real function V3BNO2_me(i,j,k,l,m,n,Jij,Tij,Tlm,Ttot)

integer :: i,j,k,l,m,n,Jij,Tij,Tlm,Ttot,ii,jj,ipoz,ipozj
integer :: ibra,iket,phase_bra,phase_ket

V3BNO2_me=0.0

if (i >= j ) then
  ibra=lpoint(i,j,k,Jij,Tij,Ttot)
  phase_bra=1
  else 
  ibra=lpoint(j,i,k,Jij,Tij,Ttot)
  phase_bra=(-1)**(Jij+Tij-(lev1pn(i)%j2+lev1pn(j)%j2)/2)
endif

if (l>=m) then
  iket=lpoint(l,m,n,Jij,Tlm,Ttot)
  phase_ket=1
  else
  iket=lpoint(m,l,n,Jij,Tlm,Ttot)
  phase_ket=(-1)**(Jij+Tlm-(lev1pn(l)%j2+lev1pn(m)%j2)/2)
endif



if(ibra.ne.0.and.iket.ne.0) then
 V3BNO2_me=V3Bno(max(ibra,iket),min(ibra,iket))*float(phase_bra*phase_ket)
endif 

return 

end function V3BNO2_me

end module v3body_no2b
