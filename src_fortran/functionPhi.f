      recursive subroutine funPhi(g, s, alp, valueX, valueKin, valueWin, 
     *n, absN, phi)
      real g, s, alp
      real valueX, valueKin, valueWin
      integer n, absN
      real phi
      real jN(absN+2)!, jN1(absN+2)
      real sinA, cosA, sinX, cosX
      real Rf, Lf, Pf, Sf, Df, znp, znm, w, mu
      real phi1
C
!      call besselJ(g, alp, valueX, valueKin, n, absN, jN)
       call besselJNew(g, alp, valueX, valueKin, n, absN, jN)
C      call besselJnew(g, alp, valueX, valueKin, n, absN, jN1)
C      print 5,  jN(1), jN(2), jN(5), valueKin*valueX*sin(alp)
C      print 5,  jN1(1), jN1(2), jN1(5), valueKin*valueX*sin(alp)
C    5 format(4E15.7)
C      pause
      
      sinA = sin(alp)
      cosA = cos(alp)
      sinX = valueX/sqrt(1+valueX*valueX)
      cosX = 1.0/sqrt(1+valueX*valueX)
C
      w = 1.0/valueWin
      znp = 1.0/(1.0+w)
      znm = 1.0/(1.0-w)
      mu = valueKin*valueKin*w*w
      Pf = -s*s*w*w
      Rf = 1+Pf*znm
      Lf = 1+Pf*znp
      Sf = 0.5*(Rf+Lf)
      Df = 0.5*(Rf-Lf)
      Pf = 1 + Pf
C
      phi1 = Df*Df/(mu-Sf)/(mu-Sf)
      phi2 = (Pf-mu*sinX*sinX)*(Pf-mu*sinX*sinX)/mu/mu
      phi3 = Pf*Pf*cosX*cosX/mu/mu
C
      if( n .eq. 0) then
        phi4 = -cosA*sinX*cosX*jN(1)/sinA
        phi5 = (1+Df/(mu-Sf))*jN(2) - (1-Df/(mu-Sf))*jN(2)
      else if( n > 0 ) then
        phi4 = -cosA*sinX*cosX*jN(absN+1)/sinA
        phi5 = (1+Df/(mu-Sf))*jN(absN+2) + (1-Df/(mu-Sf))*jN(absN)
      else
        phi4 = -cosA*sinX*cosX*jN(absN+1)/sinA
        phi5 = (1+Df/(mu-Sf))*jN(absN) + (1-Df/(mu-Sf))*jN(absN+2)
      end if
C
      phi6 = phi4 + 0.5*(mu*sinX*sinX-Pf)/mu*phi5
      phi = phi6*phi6/(phi1*phi2+phi3)
C
C      if( n .eq. -4 ) then
C      print 6,  phi, jN(absN),  jN(absN+2)
C    6 format(3E15.7)
C      end if
C
      end
      
