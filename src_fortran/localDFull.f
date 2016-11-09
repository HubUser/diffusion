      recursive subroutine localDFull(gam, s, alp, numN, localDaa, 
     *localDea, localDee, valueXm, valueXw, vXmax, b, valueWm, valuedW, 
     *vWmin, vWmax)
      real gam, s, alp
      integer numX, numN
      parameter(numX=100)
      integer harm(2*numN-1)
      real arrayX(numX), arrayG(numX), arrayGNew(numX)
      real arrayKin(numX,2*numN-1,8), arrayWin(numX,2*numN-1,8)
      real phi, valueQ, valueN, valueB2, aa, valueXr, nn2
C
      real vXmin, vXmax, valueXm, valueXw, dX, tanG, factor, testSign
      parameter(vXmin=0, factor=1.0)
C
      real vWmin, vWmax, valueWm, valuedW
C
      real denDaa, denDea, denDee, daa, b
      real localDaa(2*numN+5), localDea(2*numN+5), localDee(2*numN+5)
C
!      vWmax=valueWm+1.5*valuedW
!      vWmin=valueWm-1.5*valuedW
C
C
      harm(1)=0
      do i = 2, numN
         harm(i) = (i-1)
         harm(numN-1+i) = -(i-1)
      enddo
C
      dX = (vXmax-vXmin)/(numX*1.0)
      do i = 1, numX
         arrayX(i) = vXmin + dX*i
      enddo
      call getGFull( numX, arrayX, valueXm, valueXw, arrayG)
C
C
C

C
      do in = 1, 2*numN-1
      do iX = 1, numX
      do i = 1,8
         arrayWin(iX,in,i) = 0
         arrayKin(iX,in,i) = 0
      end do
      end do
      end do

      call kwFull(gam, s, alp, numN, harm, arrayX, numX, arrayKin,
     * arrayWin)

!      call kwNew(gam, s, alp, numN, harm, arrayX, numX, arrayKin,
!     * arrayWin)
C
         localDaa(1)=alp*180.0/3.14159
         localDea(1)=alp*180.0/3.14159
         localDee(1)=alp*180.0/3.14159
         do in = 1, 2*numN-1
         localDaa(in+1)=0
         localDea(in+1)=0
         localDee(in+1)=0
         do iX = 1, numX
         denDaa = 0
         denDea = 0
         denDee = 0
         do i = 1, 8
C
C
!         valueXr = 1.0/arrayWin(iX,in,i)**2-1.0
        valueXr = (1.0/arrayWin(iX,in,i)**2-1.0)/(1.0+1./s/s)
         nn2 =  abs( arrayKin(iX,in,i)/arrayWin(iX,in,i) )
!         if( arrayWin(iX,in,i) >0 .and. arrayX(iX)**2<0.8*valueXr) then
      valueB2=getB2Full(arrayWin(iX,in,i)*b, valueWm, valuedW, vWmax, vW
     *min)
!
!      call funNTest(s, numX, arrayX, arrayG, arrayKin(iX,in,i), arrayWin
!     *(iX,in,i), valueN, iX)
!
!
         if( arrayWin(iX,in,i) >0 .and. arrayX(iX)<0.995*sqrt(valueXr)
     *.and. valueB2>0! ) then
     *.and. nn2 < 300) then
!
      tanG = 10.0!sqrt(1.0+0.25/arrayWin(iX,in,i)/arrayWin(iX,in,i))
      do ixx = 1, numX
        arrayGNew(ixx)=arrayG(ixx)
        if( arrayX(ixx) > 0.75*tanG ) then
          arrayGNew(ixx)=arrayG(ixx)/factor;
        end if
      enddo
!
!
      call funNNew(s, numX, arrayX, arrayGNew, arrayKin(iX,in,i), arrayW
     *in(iX,in,i), valueN)
!      call funNOld(s, numX, arrayX, arrayG, arrayKin(iX,in,i), arrayWin
!     *(iX,in,i), valueN)


      testSign = 1.0
      if( arrayWin(iX,in,i)+harm(in)/gam<0)then
        testSign = -1.0
      end if

      call funPhi(gam, s, alp, arrayX(iX),
     *testSign*abs(arrayKin(iX,in,i)),
     * arrayWin(iX,in,i), harm(in), abs(harm(in)), phi)
      call funQNew(gam, s, alp, arrayX(iX), arrayKin(iX,in,i), arrayWin(
     *iX, in,i),valueQ)
!      call funQOld(gam, s, alp, arrayX(iX), arrayKin(iX,in,i), arrayWin(
!     *iX, in,i),valueQ)
      aa = (harm(in)/gam/arrayWin(iX,in,i)+sin(alp)**2)
C
        daa = arrayWin(iX,in,i)**2/(1+arrayX(iX)*arrayX(iX))/valueN*valu
     *eQ*phi*valueB2*arrayGNew(iX)
C
      denDaa = denDaa + daa/cos(alp)**2*aa**2
      denDea = denDea - daa*tan(alp)*aa
      denDee = denDee + daa*sin(alp)**2

C

C      if( harm(in) .eq. -4 ) then
C      print 5,  arrayX(iX), phi
C    5 format(13E15.7)
C      end if
         end if
C
         enddo ! end i=1,3
C

         localDaa(in+1) = localDaa(in+1) + denDaa*arrayX(iX)*dX/(4*3.141
     *59)/(gam*gam-1)
         localDea(in+1) = localDea(in+1) + denDea*arrayX(iX)*dX/(4*3.141
     *59)*(gam+1)/gam**2/(gam-1)
         localDee(in+1) = localDee(in+1) + denDee*arrayX(iX)*dX/(4*3.141
     *59)/gam/(gam-1)
         enddo ! end X
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
C
      end
C
C
C
C
C     distribution of X
      recursive subroutine getGFull( numX, arrayX, valueXm, valueXw, 
     *arrayG)
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
      recursive function getB2Full( valueWin, valueWm, valuedW, vWmax, 
     *vWmin)
      external ERF
      real valueWin
      real valueWm, valuedW, vWmax, vWmin
      real getB2Full, A2, erf1, erf2
C
      getB2Full = 0
      if( valueWin < vWmax .AND. valueWin > vWmin) then
          erf1 = ERF((valueWm-vWmin)/valuedW)
          erf2 = ERF((vWmax-valueWm)/valuedW)
          A2=2.0/sqrt(3.14159)/valuedW/(erf1+erf2)
          getB2Full = A2*exp(-(valueWin-valueWm)**2/valuedW/valuedW)
      end if
!      if( valueWin < 0.2 )then
!        getB2 = 0
!      end if
      return
      end
