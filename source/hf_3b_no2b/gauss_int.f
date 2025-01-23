      double precision function gauss_int(lam,ni,li,nj,lj,b1)

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       double precision :: dlfc1(0:170), dlfc2(0:170)
       double precision :: dlgm1(0:170), gi(0:170)

       dlfc1(0) = 0.d0
       do i = 1,170
        dlfc1(i) = dlog(dble(i)) + dlfc1(i-1)  
       enddo

       dlfc2(0) = 0.d0
       dlfc2(1) = 0.d0
       do i = 2,170
        dlfc2(i) = dlog(dble(i)) + dlfc2(i-2)
       enddo

       dlgm1(0) = 0.d0
       dlgm1(1) = dlog(dsqrt(pi))
       dlgm1(2) = 0.d0
       do i=3,170
        dlgm1(i) = dlog(dble(i)/2.d0-1.d0) + dlgm1(i-2)
       enddo

       gi(0)=dsqrt(pi)/2.d0
       gi(1)=0.5d0
       do i=2,170
        if(mod(i,2).eq.0) gi(i)=dsqrt(pi)*dexp(dlfc2(i-1))
     &                                            /(2.d0**dble(1+i/2))
        if(mod(i,2).eq.1) gi(i)=dexp(dlfc2(i-1))
     &                                        /(2.d0**dble(1+(i-1)/2))
       enddo

       gauss_int=0.d0

       dd1=4.d0/(dsqrt(pi)*b1**lam)
       dd2=dsqrt(2.d0**dble(ni+nj+li+lj))
     &              *dsqrt(dexp(dlfc1(ni)+dlfc1(nj)
     &                         -dlfc2(2*ni+2*li+1)-dlfc2(2*nj+2*lj+1)))
       ss1=0.d0
       do mi=0,ni
        do mj=0,nj
         ss1=ss1+dble((-1)**(mi+mj))*gi(lam+2+li+lj+2*mi+2*mj)
     &     *dexp(dlgm1(2*ni+2*li+3)+dlgm1(2*nj+2*lj+3)-dlfc1(mi)
     &     -dlfc1(mj)-dlfc1(ni-mi)-dlfc1(nj-mj)-dlgm1(2*li+2*mi+3)
     &                                          -dlgm1(2*lj+2*mj+3))
        enddo
       enddo

       gauss_int=dd1*dd2*ss1

       return
      end
