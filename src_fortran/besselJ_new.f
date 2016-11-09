      recursive subroutine besselJNew(g, alp, valueX, valueKin, n, absN, 
     *funJ)
      real g, alp, valueX, valueKin
      integer n, absN
      real sinX, val, fact
      real funJ(absN+2)
      double precision fun(absN+2)
      double precision arg, alpha
C
      alpha = 0.0
      sinX = valueX/sqrt(1.0+valueX*valueX)
      val = sqrt(g*g-1)*valueKin*sinX*sin(alp)
      val = abs(val)
C      PRINT 2,val
C    2 FORMAT(3E16.7)
C
      arg = val
      call rjbesl( arg, alpha, absN+2, fun, ncalc )
C
C
      do i=1, absN+2
         funJ(i) = fun(i)
      end do
C
      if( n < 0 ) then
         do i=1, absN+2
           funJ(i) = funJ(i)*(-1.0)**(i-1)
         end do
      end if
C
      if( valueKin*sinX*sin(alp) > 0 ) then
         do i=1, absN+2
           funJ(i) = funJ(i)*(-1.0)**(i-1)
         end do
      end if
C
      end

C      if( n .eq. -4 .and. i .eq. absN ) then
C      print 6,  val, fact*(-val/2.0)**(i-1), F1(absN),  funJ(absN)
C    6 format(3E15.7)
C      end if
