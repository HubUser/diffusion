      subroutine funNTest(s, numX, arrayX, arrayG, valueKin, valueWin, v
     *alueN, iX)
      real s, valueKin, valueWin
      integer numX, iX
      real arrayX(numX), arrayG(numX)
      real valueN
      real dX, zn, c, co2, si2, n2, F, dZ, Z, X, Y2, Sq
C
      valueN = 0
      dX = abs(arrayX(numX)-arrayX(1))
      X=s*s/valueWin**2
      Y2=1.0/valueWin/valueWin
      n2=valueKin**2/valueWin**2
C
      co2 = 1.0/(1+arrayX(iX)*arrayX(iX))
      si2 = arrayX(iX)*arrayX(iX)/(1+arrayX(iX)*arrayX(iX))
!      A1 = 4*valueKin**3*co**2*s**2+(0-4*valueKin**3-2*co**2*valueKin*s*
!     *s-2*valueKin*s**2-4*valueKin*s**4-4*valueKin**3*s**2)*valueWin**2
!      A2 = 2*(0-valueKin**4-co**2*valueKin**2*s**2-valueKin**2*s**2-2*va
!     *lueKin**2*s**4-s**6-valueKin**4*s**2)*valueWin
C

      Sq = sqrt( (0.5*Y2*si2)**2+(1-X)**2*Y2*co2 )
      Z=1-X-0.5*Y2*si2+Sq
      dZ=2*X+Y2*si2-(0.5*(Y2*si2)**2+(3*X-1)*(X-1)*Y2*co2)/Sq
      F=-2*X*(2*X-1)/Z+X*(1-X)*dZ/Z**2
!      A1 = 4*valueKin**3*co**2-(4*s**2*valueKin+4*valueKin**3*(1+d/s/s)+
!     *a*2*valueKin*(1+co**2))*valueWin**2
!      A2 = 2*( s**4 + 2*s**2*valueKin**2 +valueKin**4*(1+d/s/s)+a*valueK
!     *in**2*(1+co**2))*valueWin
      c = valueKin/valueWin*(0.5*F+n2)/n2*0.5/(3.14159)**2
      c = abs(c)
      zn = (1+arrayX(iX)*arrayX(iX))*sqrt(1+arrayX(iX)*arrayX(iX))
      valueN = c*valueKin**2*arrayX(iX)*dX/zn
C
      end