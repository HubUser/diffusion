      subroutine funNTet2LyonsTest(s, numT, arrayT, arrayG, valueKin, va
     *lueWin, valueN)
      real s, valueKin, valueWin
      integer numT
      real arrayT(numT), arrayG(numT)
      real valueN(2)
      real dT, zn, c, co2, si2, n2, F, dZ, Z, X, Y2, Sq, mu2
      real cosTetR, Tetr
C
      valueN(1) = 0.0
      valueN(2) = 0.0
C
      X=s*s/valueWin**2
      Y2=1.0/valueWin/valueWin
!      n2=valueKin**2/valueWin**2
C
!      cosTetR = valueWin
      cosTetR = valueWin*sqrt(1.0+(1-valueWin**2)/s/s)
      Tetr = acos(cosTetR)
C
      do i = 1, numT
      mu2 = 1.0/(cos(arrayT(i))-valueWin)*s*s/valueWin
      co2 = cos(arrayT(i))**2
      si2 = sin(arrayT(i))**2
C
      Sq = sqrt( (0.5*Y2*si2)**2+(1-X)**2*Y2*co2 )
      Z=1-X-0.5*Y2*si2+Sq
      dZ=2*X+Y2*si2-(0.5*(Y2*si2)**2+(3*X-1)*(X-1)*Y2*co2)/Sq
      F=-2*X*(2*X-1)/Z+X*(1-X)*dZ/Z**2
      n2 = 1-X*(1-X)/Z
C
!      if( arrayT(i)<Tetr-0.01*3.14159/180.0 .and. n2<0 ) then
!      print*, "!!!!!"
!      pause
!      end if
      
!      c = cos(arrayT(i))*mu2**(2.5)*valueWin**3/2.0/s/s
      c = (0.5*F+n2)*sqrt(n2)/Y2
!      c = valueKin/valueWin*(0.5*F+n2)/n2
         if( i < numT ) then
            dT = arrayT(i+1)-arrayT(i)
         else
            dT = arrayT(numT)-arrayT(numT-1)
         end if
         if(arrayT(i)<Tetr-0.01*3.14159/180.0 .and. sqrt(n2)<300 ) then
         zn = sin(arrayT(i))
         valueN(1) = valueN(1) + abs(c)*dT*zn*0.5/(3.14159)**2*arrayG(i)
         end if
C
         if(arrayT(i)<Tetr-0.01*3.14159/180.0 .and. valueN(2)<arrayT(i)
     * .and. sqrt(n2)<300) then
           valueN(2)=arrayT(i)
         end if
      enddo
C
      end