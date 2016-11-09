      subroutine funNOld(s, numX, arrayX, arrayG, valueKin, valueWin, va
     *lueN)
      real s, valueKin, valueWin
      integer numX
      real arrayX(numX), arrayG(numX)
      real valueN
      real dX, zn, c, n2, F, X, Y2
C
      valueN = 0
      X=s*s/valueWin**2
      Y2=1.0/valueWin/valueWin
      n2=valueKin**2/valueWin**2
      do i = 1, numX
!      A1 = 4*valueKin**3*co**2*s**2+(0-4*valueKin**3-2*co**2*valueKin*s*
!     *s-2*valueKin*s**2-4*valueKin*s**4-4*valueKin**3*s**2)*valueWin**2
!      A2 = 2*(0-valueKin**4-co**2*valueKin**2*s**2-valueKin**2*s**2-2*va
!     *lueKin**2*s**4-s**6-valueKin**4*s**2)*valueWin
C
      F=n2*(n2-X)/X
!      A1 = 4*valueKin**3*co**2-(4*s**2*valueKin+4*valueKin**3*(1+d/s/s)+
!     *a*2*valueKin*(1+co**2))*valueWin**2
!      A2 = 2*( s**4 + 2*s**2*valueKin**2 +valueKin**4*(1+d/s/s)+a*valueK
!     *in**2*(1+co**2))*valueWin
      c = valueKin/valueWin*(0.5*F+n2)/n2*0.5/(3.14159)**2
      c = abs(c)
         if( i < numX ) then
            dX = arrayX(i+1)-arrayX(i)
         else
            dX = arrayX(numX)-arrayX(numX-1)
         end if
         zn = (1+arrayX(i)*arrayX(i))*sqrt(1+arrayX(i)*arrayX(i))
         valueN = valueN + c*valueKin**2*arrayG(i)*arrayX(i)*dX/zn
      enddo
C
      end