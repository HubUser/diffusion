      recursive subroutine funQNew(g, s, alp, valueX, valueKin, 
     *valueWin, valueQ)
      real g, s, alp, valueX
      real valueKin, valueWin
      real valueQ, co2, si2, n2, F, dZ, Z, X, Y2, Sq, vII
C
      co2 = 1.0/(1+valueX*valueX)
      si2 = valueX*valueX/(1+valueX*valueX)
C
      X=s*s/valueWin**2
      Y2=1.0/valueWin/valueWin
      n2=valueKin**2/valueWin**2
C
      Sq = sqrt( (0.5*Y2*si2)**2+(1-X)**2*Y2*co2)
      Z=1-X-0.5*Y2*si2+Sq
      dZ=2*X+Y2*si2-(0.5*(Y2*si2)**2+(3*X-1)*(X-1)*Y2*co2)/Sq
      F=-2*X*(2*X-1)/Z+X*(1-X)*dZ/Z**2
C
      vII=valueWin/valueKin*n2/(0.5*F+n2)/sqrt(co2)
C
      valueQ = sqrt(1-1.0/g/g)*cos(alp)-vII
C
!      valueQ = sqrt(1-1.0/g/g)*cos(alp)-2*valueKin/valueWin
C
      valueQ = abs(valueQ)
      if( valueQ <0.0001 ) then
        valueQ = 0
      else
        valueQ = 1.0/valueQ
      end if
      end
