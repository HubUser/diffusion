      subroutine besselJ(g, alp, valueX, valueKin, n, absN, funJ)
      real g, alp, valueX, valueKin
      integer n, absN
      real sinX, val, fact
      real funJ(absN+2), F1(absN+2), F2(absN+2), WORK(absN+2)
C
      sinX = valueX/sqrt(1.0+valueX*valueX)
      val = sqrt(g*g-1)*valueKin*sinX*sin(alp)
      val = abs(val)
C      PRINT 2,val
C    2 FORMAT(3E16.7)
C
      if( val < 1 ) then
         NMAX = max(absN+2.,2*val+4)
         call SF33R(val,absN+2,NMAX,F1,F2,WORK,IERR)
C
         fact=1
         do i=1, absN+2
              if( i .NE. 1 ) then
                fact=fact*(i-1)
              end if
             funJ(i)=F1(i)*1.0/fact*(val/2.0)**(i-1)
         enddo
C
      else
         funJ(1) = SF24R(val)
         funJ(2) = SF25R(val)
         do i=2, absN+2
           funJ(i+1)=funJ(i)*2*(i-1)/val-funJ(i-1)
         enddo
      end if
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
