      subroutine localDMSW(gam, s, alp, localDaa, localDee, valueXm, val
     *ueXw, b, valueWm, valuedW)
      real gam, s, alp
      integer numX
      parameter(numX=1000)
      real arrayX(numX), arrayG(numX)
      real arrayKin(2), arrayWin(2)
      real phi, valueQ, valueN, valueB2, aa, valueXr
C
      real vXmin, vXmax, deltaX, valueXm, valueXw, dX
C
      real vWmin, vWmax, valueWm, valuedW
C
      real denDaa,  denDee, daa, b
      real localDaa, localDee
      real eta, nHe, nO
      parameter(nHe=0.2, nO=0.2)
      parameter(eta=(nHe/4.0+nO/16.0+(1-nHe-nO))/1836.0)
C
      vWmax=4.4/1000.0!valueWm+1.5*valuedW
      vWmin=2.6/1000.0!valueWm-1.5*valuedW
!      vWmax=valueWm+2.0*valuedW
!      vWmin=valueWm-2.0*valuedW
C
!      vXmax = getX(gam, s, alp, vWmax, eta)
!      vXmin = getX(gam, s, alp, vWmin, eta)
!      deltaX =  vXmax-vXmin
      vXmax = tan(89.2*3.1415926/180.0)!vXmax + 5*deltaX
      vXmin = tan(88.67*3.1415926/180.0)!vXmin - 3*deltaX
C
      dX = (vXmax-vXmin)/(numX*1.0)
      do i = 1, numX
         arrayX(i) = vXmin + dX*i
      enddo
      call getGMSW( numX, arrayX, valueXm, valueXw, arrayG)
C
C
      localDaa=0
      localDee=0
      do iX = 1, numX
         denDaa = 0
         denDea = 0
         denDee = 0
C
C
      call kwMSW(gam, alp, arrayKin, arrayWin, s, arrayX(iX), eta)
C
      do i = 1, 2
!         valueXr = 1.0/arrayWin(iX,in,i)**2-1.0
!        valueXr = 1.0/(arrayWin(i)**2)-1.0
!         if( arrayWin(iX,in,i) >0 .and. arrayX(iX)**2<0.8*valueXr) then
!        arrayX(iX)=tan(89.193*3.14159/180.0)
!        arrayWin(i) = sqrt(1.0/1836.0/sqrt(1.0-1.0/gam/gam)/cos(alp)*sqr
!     *t(1+arrayX(iX)**2)-s**2)
        valueB2=getB2MSW(arrayWin(i)*b, valueWm, valuedW, vWmax, vWmin)
!      print 5, sqrt(1.0-1.0/gam/gam)*cos(alp)
!    5 format(13E15.7)
!      pause
!      print 5, arrayWin(i), arrayX(iX), valueB2
!    5 format(13E15.7)
!
         if( arrayWin(i) >0 .and. valueB2>0
!     *.and. arrayX(iX)**2<1./(1./1836.0-arrayWin(i)**2)
     *) then

!
      call funNMSW(s, numX, arrayX, arrayG, arrayKin(i), arrayWin(i),
     * valueN, eta)
      call funPhiMSW(gam, s, alp, arrayX(iX), arrayKin(i), arrayWin(i),
     * phi, eta)
      call funQMSW(gam, s, alp, arrayX(iX), arrayKin(i), arrayWin(i),val
     *ueQ, eta)
!
      aa = sin(alp)**2
C
        daa = arrayWin(i)**2/(1+arrayX(iX)*arrayX(iX))/valueN*valueQ*phi
     **valueB2*arrayG(iX)
C
      denDaa = denDaa + daa/cos(alp)**2*aa**2
      denDee = denDee + daa*sin(alp)**2

C
!      print 9,  daa,  valueQ
!    9 format(13E15.7)
C      if( harm(in) .eq. -4 ) then
C      end if
         end if
C
      enddo ! end i=1,3
C

         localDaa = localDaa + denDaa*arrayX(iX)*dX/(4*3.141
     *59)/(gam*gam-1)
         localDee = localDee + denDee*arrayX(iX)*dX/(4*3.141
     *59)/(gam*gam-1)*(1+1./gam)*(1+1./gam)!/gam/(gam-1)/gam*(gam+1)
         enddo ! end X
!      pause
C
C
      end
C
C
C
C
C     distribution of X
      subroutine getGMSW( numX, arrayX, valueXm, valueXw, arrayG)
      integer numX
      real arrayX(numX), arrayG(numX)
      real valueXm, valueXw
      real c
C
      do i = 1, numX
        c = (arrayX(i)-valueXm)/valueXw
        arrayG(i)=exp(-c*c)
      enddo
      end
C

C     distribution of omega
      function getB2MSW( valueWin, valueWm, valuedW, vWmax, vWmin)
      external ERF
      real valueWin
      real valueWm, valuedW, vWmax, vWmin
      real getB2MSW, A2, erf1, erf2
C
      getB2MSW = 0
      if( valueWin < vWmax .AND. valueWin > vWmin) then
          erf1 = ERF((valueWm-vWmin)/valuedW)
          erf2 = ERF((vWmax-valueWm)/valuedW)
          A2=2.0/sqrt(3.14159)/valuedW/(erf1+erf2)
          getB2MSW = A2*exp(-(valueWin-valueWm)**2/valuedW/valuedW)
      end if
!      if( valueWin < 0.2 )then
!        getB2 = 0
!      end if
      return
      end
C     distribution of X
C
C
C
      subroutine kwMSW(gam, alp, arrayKin, arrayWin, s, valueX, eta)
      real arrayKin(2), arrayWin(2)
      real s, alp, gam, valueX, a1, eta, cost, u, d, fp, fm
      cost = 1.0/sqrt(1.0+valueX**2)
      u = sqrt(1.0-1.0/gam/gam)*cos(alp)*cost
      a1 = cost**2+eta
      do i = 1, 2
      arrayKin(i)=0
      arrayWin(i)=0
      end do
C
      d=a1**2-4*(s*u*cost)**2
      fp = 0.5*(a1+sqrt(d))/u/u
      fm = 0.5*(a1-sqrt(d))/u/u
      if(d>0) then
        if(fp>s**2) then
          arrayKin(1)=sqrt(fp-s**2)
          arrayWin(1)=arrayKin(1)*u
        end if
        if(fm>s**2) then
          arrayKin(2)=sqrt(fm-s**2)
          arrayWin(2)=arrayKin(2)*u
        end if
      end if
C
      end
C
C
C
      subroutine funQMSW(g, s, alp, valueX, valueKin, valueWin, valueQ,
     *eta)
      real g, s, alp, valueX, eta
      real valueKin, valueWin
      real valueQ, Teh, co2, n2, vII
C
      co2 = 1.0/(1+valueX*valueX)
      Teh = 1.0/(1+s**2/valueKin/valueKin)
      n2 = 2*Teh**3*co2+Teh**2*eta
C
C
      vII = s**2/(valueWin*valueKin**3)*n2/sqrt(co2)
C
      valueQ = sqrt(1-1.0/g/g)*cos(alp)-vII
!      print 5, valueWin, valueKin, n2, Teh
!    5 format(13E15.7)
!!        valueQ = sqrt(1-1.0/g/g)*cos(alp)*valueWin**2*1836.0
C
!      valueQ = sqrt(1-1.0/g/g)*cos(alp)-2*valueKin/valueWin
C
      valueQ = abs(valueQ)
      if( valueQ <0.00001 ) then
        valueQ = 0
      else
        valueQ = 1.0/valueQ
      end if
      end
C
C
C
      subroutine funNMSW(s, numX, arrayX, arrayG, valueKin, valueWin, va
     *lueN, eta)
      real s, valueKin, valueWin, eta
      integer numX
      real arrayX(numX), arrayG(numX)
      real valueN
      real dX, zn, n2, c, co2, Teh
C
      valueN = 0
      Teh = 1.0/(1+s**2/valueKin/valueKin)
      do i = 1, numX
      co2 = 1.0/(1+arrayX(i)*arrayX(i))
      n2 = 2*Teh**3*co2+Teh**2*eta
C
      c = valueWin*valueKin**3/s/s/n2*0.5/(3.14159)**2
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
C
C
C
      function getX(g, s, alp, valueW, eta)
      real g, s, alp, valueW, eta
      real u, getX, d
C
      u = sqrt(1.0-1.0/g/g)*cos(alp)
      d = sqrt(eta-valueW**2)/s/u
      getX = 0
      if( d< 1.) then
      getX = sqrt(1.0/d/d-1.0)
      end if
      
      return
      end
C
C
C
      subroutine funPhiMSW(g, s, alp, valueX, valueKin, valueWin, phi, e
     *ta)
      real g, s, alp
      real valueX, valueKin, valueWin, eta
      real phi
      real jN(2)!, jN1(absN+2)
      real sinA, cosA, sinX, cosX
      real Rf, Lf, PfE, PfI, Sf, Df, znpE, znmE, znpI, znmI, w, mu
      real phi1, phi2, phi3, phi4, phi5, phi6
C
!      call besselJ(g, alp, valueX, valueKin, n, absN, jN)
       call besselJNew(g, alp, valueX, valueKin, 0, 0, jN)
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
      znpE = 1.0/(1.0+w)
      znmE = 1.0/(1.0-w)
      znpI = 1.0/(1.0+w*eta)
      znmI = 1.0/(1.0-w*eta)
      mu = valueKin*valueKin*w*w
      PfE = -s*s*w*w
      PfI = -s*s*w*w*eta
      Rf = 1+PfE*znpE+PfI*znpI
      Lf = 1+PfE*znmE+PfI*znmI
      Sf = 0.5*(Rf+Lf)
      Df = 0.5*(Rf-Lf)
      Pf = 1 + PfE+PfI
C
      phi1 = Df*Df/(mu-Sf)/(mu-Sf)
      phi2 = (Pf-mu*sinX*sinX)*(Pf-mu*sinX*sinX)/mu/mu
      phi3 = Pf*Pf*cosX*cosX/mu/mu
C
      phi4 = -cosA*sinX*cosX*jN(1)/sinA
      phi5 = (1+Df/(mu-Sf))*jN(2) - (1-Df/(mu-Sf))*jN(2)
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