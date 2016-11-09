      recursive subroutine localDTet2Dst(gam, s, alp, numN, localDaa, lo
     *calDea, localDee, valueTm, valueTw, Bw1, Bw2, testVmax, b, 
     *valueWm, testNumber)
C> input parameter variable calculated from particle energy in keV
      real gam
      real s, alp
      integer numT, numN, testNumber
!      parameter(numT=200+800)
!      parameter(numT=1000+2000)
      parameter(numT=1000+1000)
      integer harm(2*numN-1), harmTest(1)
      real arrayT(numT), arrayX(numT), arrayX1(1), arrayG(numT)
      real arrayKin(1,1,8), arrayWin(1,1,8)
      real phi, valueQ, valueN, valueB2, aa, valueXr, nn2, cosTetR, Tetr
      real funAA(numT), funEE(numT), funEA(numT)
C
      real vTmin, vTmax, valueTm(2), valueTw(2), Bw1, Bw2, dT1, dT2, dT,
     *testVmax, testSign, coefTan
      parameter(vTmin=0, vTmax=80*3.14159/180.0)
C
      real vWmin, vWmax, valueWm, valuedW
      parameter(valuedW=0.15)
C
      real denDaa, denDea, denDee, daa, b
      real localDaa(2*numN+5), localDea(2*numN+5), localDee(2*numN+5)

C
      vWmax=valueWm+1.5*valuedW
      vWmin=valueWm-1.5*valuedW
      if(vWmin<0.05) then
        vWmin=0.05
      end if
C
C
      harm(1)=0
      do i = 2, numN
         harm(i) = (i-1)
         harm(numN-1+i) = -(i-1)
      enddo
C

!      dT2 = 10.0/2000.0*(3.14159/180.0)
!      dT1 = (vTmax-vTmin)/(numT*1.0-1.0-2000)
!      dT2 = 10.0/800.0*(3.14159/180.0)
!      dT1 = (vTmax-vTmin)/(numT*1.0-1.0-800)
      dT2 = 10.0/1000.0*(3.14159/180.0)
      dT1 = (vTmax-vTmin)/(numT*1.0-1.0-1000)
      arrayT(1) = 0.0
      arrayX(1) = 0
      do i = 2, numT
         if( arrayT(i-1)+dT2<80*3.1415/180) then
            arrayT(i) = arrayT(i-1) + dT1
         else
           arrayT(i) = arrayT(i-1) + dT2
         end if
          arrayX(i) = tan(arrayT(i))
      enddo
      call getGTet22( numT, arrayT, valueTm, valueTw, Bw1, Bw2, arrayG,
     *testVmax)
!      call  getGtest( numT, arrayX, 0.0, 0.577, arrayG)
C
C      print 55, Bw1, Bw2
C      pause
C      open (unit = 77, file = 'testG.dat')
C      do iX=1,numT
C      write(77,55) arrayT(iX), arrayG(iX)
C      end do
C   55  format(5E15.7)
C      close(77)
C      pause
C
C

C
      localDaa(1)=alp*180.0/3.14159
      localDea(1)=alp*180.0/3.14159
      localDee(1)=alp*180.0/3.14159
C
      do in = 1, 2*numN-1
         localDaa(in+1)=0
         localDea(in+1)=0
         localDee(in+1)=0
      end do
!      call kwNew(gam, s, alp, numN, harm, arrayX, numX, arrayKin,
!     * arrayWin)
C

      do in = 1, 2*numN-1
C
      harmTest(1) = harm(in)
      do iT = 1, numT
      funAA(iT) = 0.0
      funEA(iT) = 0.0
      funEE(iT) = 0.0
C
C
      if( arrayT(iT)<testVmax .and. arrayG(iT)>0 ) then
C
      do i = 1, 8
      arrayKin(1, 1, i )=0
      arrayWin(1, 1, i )=0
      end do
C
      arrayX1(1) = arrayX(iT)
      call kwFull(gam, s, alp, 1, harmTest, arrayX1, 1, arrayKin,
     * arrayWin)
C
      denDaa = 0
      denDea = 0
      denDee = 0
      do i = 1, 8
C
C
      valueXr = (1.0/arrayWin(1,1,i)**2-1.0)/(1.0+1./s/s)
      cosTetR = arrayWin(1,1,i)*sqrt(1.0+(1-arrayWin(1,1,i)**2)/s/s)
      Tetr = acos(cosTetR)
C
      nn2 =  abs( arrayKin(1,1,i)/arrayWin(1,1,i) )
      valueB2=getB2Tet22(arrayWin(1,1,i)*b, valueWm, valuedW, vWmax, v
     *Wmin)
!
!
      coefTan = 1.0!-s*s/650.0;
!
      if( arrayWin(1,1,i) >0 .and.
     * arrayX(iT)**2<coefTan*valueXr .and.
     * arrayT(iT)<Tetr-dT2
     *.and. valueB2>0 .and. nn2 <300.0  ) then!.and. arrayKin(1,in,i) >0) then
!
      testSign = -1.0
      if( arrayWin(1,1,i)+harm(in)/gam<0)then
        testSign = 1.0
      end if
!
      call funNTet2(s, numT, arrayT, arrayG, arrayKin(1,1,i), arrayWin
     *(1,1,i), valueN, 300.0, coefTan)
!      call funNNew(s, numT, arrayX, arrayG, arrayKin(iT,in,i), arrayWin(
!     *iT,in,i), valueN)
C
      call funPhi(gam, s, alp, arrayX(iT),testSign*abs(arrayKin(1,1,i))
     *,arrayWin(1,1,i), harm(in), abs(harm(in)), phi)
      call funQNew(gam, s, alp, arrayX(iT), arrayKin(1,1,i), arrayWin(
     *1, 1,i), valueQ)
C
      aa = (harm(in)/gam/arrayWin(1,1,i)+sin(alp)**2)
C
      daa = arrayWin(1,1,i)**2/valueN*valueQ*phi*valueB2*arrayG(iT)
C
      denDaa = denDaa + daa/cos(alp)**2*aa**2
      denDea = denDea - daa*tan(alp)*aa
      denDee = denDee + daa*sin(alp)**2
C
      end if    ! end if theta<thetaR!
      enddo ! end roots!
C
C
      funAA(iT) =  denDaa/(4*3.14159)/(gam*gam-1)/cos(arrayT(iT))
      funEA(iT) = denDea/(4*3.14159)*(gam+1)/gam**2/(gam-1)/cos(ar
     *rayT(iT))
      funEE(iT) =  denDee/(4*3.14159)/gam/(gam-1)/cos(arrayT(iT))
C
      end if  ! end if theta<thetaMax!
      enddo ! end theta!
C
C
      do iT= 1, numT-1
      localDaa(in+1) = localDaa(in+1) + 0.5*(funAA(iT)+funAA(iT+1))*(ar
     *rayT(iT+1)-arrayT(iT))
      localDea(in+1) = localDea(in+1) + 0.5*(funEA(iT)+funEA(iT+1))*(arr
     *ayT(iT+1)-arrayT(iT))
      localDee(in+1) = localDee(in+1) + 0.5*(funEE(iT)+funEE(iT+1))*(arr
     *ayT(iT+1)-arrayT(iT))
      end do
      enddo ! end N
C
      localDaa(2*numN+1) = 0
      localDea(2*numN+1) = 0
      localDee(2*numN+1) = 0
      do in = 1, 2*numN-1
        localDaa(2*numN+1) = localDaa(2*numN+1) + localDaa(in+1)
        localDea(2*numN+1) = localDea(2*numN+1) + localDea(in+1)
        localDee(2*numN+1) = localDee(2*numN+1) + localDee(in+1)
      end do
C
      localDaa(2*numN+2) = 0
      localDea(2*numN+2) = 0
      localDee(2*numN+2) = 0
      do in = 1, 6
        localDaa(2*numN+2) = localDaa(2*numN+2) + localDaa(in+1) + local
     *Daa(numN+in+1)
        localDea(2*numN+2) = localDea(2*numN+2) + localDea(in+1) + local
     *Dea(numN+in+1)
        localDee(2*numN+2) = localDee(2*numN+2) + localDee(in+1) + local
     *Dee(numN+in+1)
      end do
C
      if(numN>6) then
      localDaa(2*numN+3) = 0
      localDea(2*numN+3) = 0
      localDee(2*numN+3) = 0
      do in = 1, 11
        localDaa(2*numN+3) = localDaa(2*numN+3) + localDaa(in+1) + local
     *Daa(numN+in+1)
        localDea(2*numN+3) = localDea(2*numN+3) + localDea(in+1) + local
     *Dea(numN+in+1)
        localDee(2*numN+3) = localDee(2*numN+3) + localDee(in+1) + local
     *Dee(numN+in+1)
      end do
C
      localDaa(2*numN+4) = 0
      localDea(2*numN+4) = 0
      localDee(2*numN+4) = 0
      do in = 1, 21
        localDaa(2*numN+4) = localDaa(2*numN+4) + localDaa(in+1) + local
     *Daa(numN+in+1)
        localDea(2*numN+4) = localDea(2*numN+4) + localDea(in+1) + local
     *Dea(numN+in+1)
        localDee(2*numN+4) = localDee(2*numN+4) + localDee(in+1) + local
     *Dee(numN+in+1)
      end do
C
      localDaa(2*numN+5) = 0
      localDea(2*numN+5) = 0
      localDee(2*numN+5) = 0
      do in = 1, 41
        localDaa(2*numN+5) = localDaa(2*numN+5) + localDaa(in+1) + local
     *Daa(numN+in+1)
        localDea(2*numN+5) = localDea(2*numN+5) + localDea(in+1) + local
     *Dea(numN+in+1)
        localDee(2*numN+5) = localDee(2*numN+5) + localDee(in+1) + local
     *Dee(numN+in+1)
      end do
      end if
C
      end
C
C
C
C
C     distribution of X
      recursive subroutine getGTet22(numT,arrayT, valueTm, valueTw, Bw1, 
     *Bw2, arrayG, testVmax)
      integer numT
      real arrayT(numT), arrayG(numT)
      real valueTm(2), valueTw(2), Bw1, Bw2
      real c1, c2, testVmax
C
      do i = 1, numT
        arrayG(i) = 0
        if(arrayT(i)<testVmax) then
         c1 = (arrayT(i)-valueTm(1))/valueTw(1)
          c2 = (arrayT(i)-valueTm(2))/valueTw(2)
          arrayG(i)=Bw2*Bw2*exp(-c2*c2)
          if(arrayT(i)<0.785) then
           arrayG(i)=arrayG(i) + Bw1*Bw1*exp(-c1*c1)
          end if
        end if
      enddo
      end
C

C     distribution of omega
      recursive function getB2Tet22( valueWin, valueWm, valuedW, vWmax, 
     *vWmin)
      external ERF
      real valueWin
      real valueWm, valuedW, vWmax, vWmin
      real getB2Tet22, A2, erf1, erf2
C
      getB2Tet22 = 0
      if( valueWin < vWmax .AND. valueWin > vWmin) then
          erf1 = ERF((valueWm-vWmin)/valuedW)
          erf2 = ERF((vWmax-valueWm)/valuedW)
          A2=2.0/sqrt(3.14159)/valuedW/(erf1+erf2)
          getB2Tet22 = A2*exp(-(valueWin-valueWm)**2/valuedW/valuedW)
      end if
!      if( valueWin < 0.2 )then
!        getB2 = 0
!      end if
      return
      end
