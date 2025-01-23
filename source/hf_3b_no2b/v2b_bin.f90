module NN_bin

! several functions adopted and modified from NuHamil code by T.Myiagi

use technical
use Tcm_2body

contains


subroutine read_2b_bin(e1max,e2max,hbarom,ACM)

integer :: io, runit = 22

real(8), allocatable :: v(:)

integer :: e1max, e2max,lmax, nlev
integer(kind=8) :: nelms
integer :: a, b, c, d, dmax
integer :: la, ja, ea
integer :: lb, jb, eb
integer :: lc, jc, ec
integer :: ld, jd, ed
integer :: J, Jmin, Jmax
integer(8) :: nelm, icnt
integer :: ap, bp, cp, dp
integer :: an, bn, cn, dn
integer :: phase
real(8) :: me_00, me_pp, me_10, me_nn, fact
real(8) :: hbarom, ACM
    

nelms=count_2bme(e1max, e2max)

write(*,*) "Calculated number of 2b MEs: ", nelms


allocate(v(nelms))
open(runit, file='2belem.bin', action='read',iostat=io, form='unformatted', access='stream')
read(runit) v
close(runit)

write(*,*) "Loaded number of 2b MEs: ", nelms

nlev=(e1max+1)*(e1max+2)/2  ! number of single-particle levels

    icnt = 0
    do a = 1, nlev
      la = lev1pn(a)%l
      ja = lev1pn(a)%j2
      ea = lev1pn(a)%N
      if(ea > e1max) cycle
       ap=2*a-1 
       an=2*a
!      ap = sps_me2j%iso2pn(sps,a,-1)
!      an = sps_me2j%iso2pn(sps,a, 1)
      do b = 1, a
        lb = lev1pn(b)%l
        jb = lev1pn(b)%j2
        eb = lev1pn(b)%N
        bp=2*b-1
        bn=2*b
!        bp = sps_me2j%iso2pn(sps,b,-1)
!        bn = sps_me2j%iso2pn(sps,b, 1)
        if(ea + eb > e2max) cycle
        do c = 1, a
          lc = lev1pn(c)%l
          jc = lev1pn(c)%j2
          ec = lev1pn(c)%N
          cp=2*c-1
          cn=2*c
!          cp = sps_me2j%iso2pn(sps,c,-1)
!          cn = sps_me2j%iso2pn(sps,c, 1)
          dmax = c
          if(a == c) dmax = b
          do d = 1, dmax
            ld = lev1pn(d)%l
            jd = lev1pn(d)%j2
            ed = lev1pn(d)%N
            dp=2*d-1
            dn=2*d
!            dp = sps_me2j%iso2pn(sps,d,-1)
!            dn = sps_me2j%iso2pn(sps,d, 1)
            if(ec + ed > e2max) cycle
            if( (-1)**(la+lb) /= (-1)**(lc+ld) ) cycle
            Jmin = max(abs(ja-jb), abs(jc-jd))/2
            Jmax = min(   (ja+jb),    (jc+jd))/2
            if(Jmin > Jmax) cycle
            do J = Jmin, Jmax
              me_00 = v(icnt+1)
              me_nn = v(icnt+2)
              me_10 = v(icnt+3)
              me_pp = v(icnt+4)
              icnt = icnt + 4

!              if(a == b .and. mod(J,2) == 1 .and. abs(me_nn) > 1.d-8) call me2j_read_warning(a,b,c,d,J,me_nn)
!              if(c == d .and. mod(J,2) == 1 .and. abs(me_nn) > 1.d-8) call me2j_read_warning(a,b,c,d,J,me_nn)
!              if(a == b .and. mod(J,2) == 1 .and. abs(me_pp) > 1.d-8) call me2j_read_warning(a,b,c,d,J,me_nn)
!              if(c == d .and. mod(J,2) == 1 .and. abs(me_pp) > 1.d-8) call me2j_read_warning(a,b,c,d,J,me_nn)

!              if(ap == 0) cycle
!              if(bp == 0) cycle
!              if(cp == 0) cycle
!              if(dp == 0) cycle

!              if(sps%orb(ap)%e + sps%orb(bp)%e > ms%e2max) cycle
!              if(sps%orb(cp)%e + sps%orb(dp)%e > ms%e2max) cycle



              fact = 1.d0
              phase = 1
              if(a == b) fact = fact / dsqrt(2.d0)
              if(c == d) fact = fact / dsqrt(2.d0)

              if(abs(1.d0-fact) < 1.d-4 .or. mod(J,2) == 0) then
                 call store_2bme(a,b,c,d,J,1,me_nn,hbarom,ACM)
                 call store_2bme(a,b,c,d,J,-1,me_pp,hbarom,ACM)
              endif 

!                write(933,'(9i5,3f10.5)')ap,bp,cp,dp,an,bn,cn,dn,J,0.5d0*(me_10+me_00),0.5d0*(me_10-me_00),dabs(me_10-me_00)-dabs(me_10+me_00)
!                write(933,'(5i5,2f10.5)')a,b,c,d,J,0.5d0*(me_10+me_00),0.5d0*(me_10-me_00)
 
               if(abs(1.d0-fact) < 1.d-4 .or. mod(J,2) == 0 ) then              
!                 phase=1
                 call store_2bme(a,b,c,d,J,0,0.5d0*me_10,hbarom,ACM)
                 if (c/=d) then
                  phase=(-1)**((lev1pn(c)%j2+lev1pn(d)%j2)/2+J)  
                  call store_2bme(a,b,d,c,J,0,phase*(-0.5d0*me_10),hbarom,ACM)
                 endif 
                 if (a/=b .and. c/=d) then
                  phase = (-1)**((lev1pn(a)%j2+lev1pn(b)%j2+lev1pn(c)%j2+lev1pn(d)%j2)/2)  
                  call store_2bme(b,a,d,c,J,0,phase*(0.5d0*me_10),hbarom,ACM)
                  endif 
                 if (a/=b .and. (a/=c .or. b/=d)) then 
                  phase = (-1)**((lev1pn(a)%j2+lev1pn(b)%j2)/2+J)
                  call store_2bme(b,a,c,d,J,0,phase*(-0.5d0*me_10),hbarom,ACM)
                 endif
                endif

              if(abs(1.d0-fact) < 1.d-4 .or. mod(J,2) == 1 ) then              
                  call store_2bme(a,b,c,d,J,0,0.5d0*me_00,hbarom,ACM)
                 if (c/=d) then
                  phase=(-1)**((lev1pn(c)%j2+lev1pn(d)%j2)/2+J)  
                  call store_2bme(a,b,d,c,J,0,phase*(0.5d0*me_00),hbarom,ACM)
                 endif 
                 if (a/=b .and. c/=d) then
                  phase = (-1)**((lev1pn(a)%j2+lev1pn(b)%j2+lev1pn(c)%j2+lev1pn(d)%j2)/2)  
                  call store_2bme(b,a,d,c,J,0,phase*(0.5d0*me_00),hbarom,ACM)
                  endif 
                 if (a/=b .and. (a/=c .or. b/=d)) then 
                  phase = (-1)**((lev1pn(a)%j2+lev1pn(b)%j2)/2+J)
                  call store_2bme(b,a,c,d,J,0,phase*(0.5d0*me_00),hbarom,ACM)
                 endif
              endif





                 
!                call two%SetTwBME(ap,bp,cp,dp,J,me_pp*fact)
!                call two%SetTwBME(an,bn,cn,dn,J,me_nn*fact)
!                call two%AddToTwBME(ap,bn,cp,dn,J,0.5d0*me_10)
!                if(c/=d) call two%AddToTwBME(ap,bn,cn,dp,J,0.5d0*me_10)
!                if(a/=b .and. c/=d) &
!                    &    call two%AddToTwBME(an,bp,cn,dp,J,0.5d0*me_10)
!                if(a/=b .and. (a/=c .or. b/=d)) &
!                    &    call two%AddToTwBME(an,bp,cp,dn,J,0.5d0*me_10)
!              end if

!              if(abs(1.d0-fact) < 1.d-4 .or. mod(J,2) == 1) then
!                call two%AddToTwBME(ap,bn,cp,dn,J,0.5d0*me_00)
!                if(c/=d) call two%AddToTwBME(ap,bn,cn,dp,J,-0.5d0*me_00)
!                if(a/=b .and. c/=d) &
!                    &    call two%AddToTwBME(an,bp,cn,dp,J, 0.5d0*me_00)
!                if(a/=b .and. (a/=c .or. b/=d)) &
!                    &    call two%AddToTwBME(an,bp,cp,dn,J,-0.5d0*me_00)
!              end if
            end do
          end do
        end do
      end do
    end do
    deallocate(v)
!    call sps_me2j%fin()
    return

end subroutine read_2b_bin


function count_2bme(e_cut,e2max) result(r)
!    type(OrbitsIsospin), intent(in) :: sps
    integer, intent(in) :: e_cut,e2max
    integer(8) :: r
    integer :: a, b, c, d, dmax
    integer :: la, ja, ea
    integer :: lb, jb, eb
    integer :: lc, jc, ec
    integer :: ld, jd, ed
    integer :: J, Jmin, Jmax

    r = 0
    nlev=(e_cut+1)*(e_cut+2)/2  ! number of single-particle levels

    do a = 1, nlev
      la = lev1pn(a)%l
      ja = lev1pn(a)%j2
      ea = lev1pn(a)%N
      if(ea>e_cut) cycle
      do b = 1, a
        lb = lev1pn(b)%l
        jb = lev1pn(b)%j2
        eb = lev1pn(b)%N
        if(ea + eb > e2max) cycle
        do c = 1, a
          lc = lev1pn(c)%l
          jc = lev1pn(c)%j2
          ec = lev1pn(c)%N
          dmax = c
          if(a == c) dmax = b
          do d = 1, dmax
            ld = lev1pn(d)%l
            jd = lev1pn(d)%j2
            ed = lev1pn(d)%N
            if(ec + ed > e2max) cycle
            if( (-1)**(la+lb) /= (-1)**(lc+ld) ) cycle
            Jmin = max(abs(ja-jb), abs(jc-jd))/2
            Jmax = min(   (ja+jb),    (jc+jd))/2
            if(Jmin > Jmax) cycle
            do J = Jmin, Jmax
              r = r + 4
            end do
          end do
        end do
      end do
    end do

  end function count_2bme

  subroutine store_2bme(ia,ib,ic,id,J,Tz,me,hbarom,ACM)

    integer :: a, b, c, d, J, Tz, phase
    integer :: ia, ib, ic, id
    real(kind=4) :: me4
    real(8) :: me
    real(8) :: fact, hbarom, ACM,fact_ab,fact_cd

    a=lp1(ia)
    b=lp1(ib)
    c=lp1(ic)
    d=lp1(id)
    me4=real(me,4)
!    2-body part of CM Hamiltonian
!    be careful with the factor of dsqrt(2.d0) 
!    fact_ab=1.d0
!    if(a == b) fact_ab=dsqrt(2.d0)
!    fact_cd=1.d0
!    if(c == d) fact_cd=dsqrt(2.d0)
    
!    if (Tz==0)  me=me+ hbarom/ACM*T_nn(a,b,c,d,J)
!    if (Tz==1 .or. Tz == -1) me=me+hbarom/ACM*T_nn_asym(a,b,c,d,J)*fact_ab*fact_cd
  
    if (Tz==1) then
       Vnn(a,b,c,d,J)=me4
       Vnn(c,d,a,b,J)=me4
!       phase=(-1)**((lev1pn(a)%j2+lev1pn(b)%j2)/2+J)
       phase=(-1)**((levn(a)%j2+levn(b)%j2)/2+J) 
       Vnn(b,a,c,d,J)=-1*phase*me4
       Vnn(c,d,b,a,J)=-1*phase*me4
!       phase=(-1)**((lev1pn(c)%j2+lev1pn(d)%j2)/2+J)
       phase=(-1)**((levn(c)%j2+levn(d)%j2)/2+J)
       Vnn(a,b,d,c,J)=-1*phase*me4
       Vnn(d,c,a,b,J)=-1*phase*me4
!       phase=(-1)**((lev1pn(a)%j2+lev1pn(b)%j2+lev1pn(c)%j2+lev1pn(d)%j2)/2)
       phase=(-1)**((levn(a)%j2+levn(b)%j2+levn(c)%j2+levn(d)%j2)/2)
       Vnn(b,a,d,c,J)=phase*me4
       Vnn(d,c,b,a,J)=phase*me4
    endif

    if (Tz==-1) then
       Vpp(a,b,c,d,J)=me4
       Vpp(c,d,a,b,J)=me4
!       phase=(-1)**((lev1pn(a)%j2+lev1pn(b)%j2)/2+J)
       phase=(-1)**((levp(a)%j2+levp(b)%j2)/2+J)
       Vpp(b,a,c,d,J)=-1*phase*me4
       Vpp(c,d,b,a,J)=-1*phase*me4
!       phase=(-1)**((lev1pn(c)%j2+lev1pn(d)%j2)/2+J)
       phase=(-1)**((levp(c)%j2+levp(d)%j2)/2+J)
       Vpp(a,b,d,c,J)=-1*phase*me4
       Vpp(d,c,a,b,J)=-1*phase*me4
!       phase=(-1)**((lev1pn(a)%j2+lev1pn(b)%j2+lev1pn(c)%j2+lev1pn(d)%j2)/2)
       phase=(-1)**((levp(a)%j2+levp(b)%j2+levp(c)%j2+levp(d)%j2)/2)
       Vpp(b,a,d,c,J)=phase*me4
       Vpp(d,c,b,a,J)=phase*me4
    endif

    if (Tz==0) then
       Vpn(a,b,c,d,J)=Vpn(a,b,c,d,J)+me4
       Vpn(c,d,a,b,J)=Vpn(a,b,c,d,J)
    endif

  end subroutine store_2bme


end module NN_bin