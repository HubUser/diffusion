      recursive subroutine funNTet2(s, numT, arrayT, arrayG, valueKin, 
     *valueWin, valueN, Nmax)
      real s, valueKin, valueWin
      integer numT
      real arrayT(numT), arrayG(numT)
      real valueN, Nmax
      real dT, zn, c, co2, si2, n2, F, dZ, Z, X, Y2, Sq
      real cosTetR, Tetr
C
      valueN = 0
      X=s*s/valueWin**2
      Y2=1.0/valueWin/valueWin
!      n2=valueKin**2/valueWin**2
C
      cosTetR = valueWin*sqrt(1.0+(1-valueWin**2)/s/s)
      Tetr = acos(cosTetR)
C
      do i = 1, numT
      co2 = cos(arrayT(i))**2
      si2 = sin(arrayT(i))**2
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
      n2 = 1-X*(1-X)/Z
      c = (0.5*F+n2)*sqrt(n2)*0.5/(3.14159)**2/Y2
C      n2=valueKin**2/valueWin**2
C      c = valueKin**2*valueKin/valueWin*(0.5*F+n2)/n2*0.5/(3.14159)**2
      c = abs(c)
         if( i < numT ) then
            dT = arrayT(i+1)-arrayT(i)
         else
            dT = arrayT(numT)-arrayT(numT-1)
         end if
         if(arrayT(i)<Tetr-0.01*3.14159/180.0 .and. n2>0
!     *   .and. tan(arrayT(i))**2<0.99*(1.0/cosTetR**2-1)
     *   .and. sqrt(n2)<Nmax) then
         zn = 1.0 !sin(arrayT(i))
         valueN = valueN + c*arrayG(i)*dT*zn
         end if
      enddo
C
      end
