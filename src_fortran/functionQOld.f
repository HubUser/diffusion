      subroutine funQOld(g, s, alp, valueX, valueKin, valueWin, valueQ)
       real g, s, alp, valueX
      real valueKin, valueWin
      real valueQ, co2, n2, F, X, Y2, vII
C
      co2 = 1.0/(1+valueX*valueX)
C
      X=s*s/valueWin**2
      Y2=1.0/valueWin/valueWin
      n2=valueKin**2/valueWin**2
C
      F=n2*(n2-X)/X
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