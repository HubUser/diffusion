      recursive subroutine kwFull(g, s, alp, numN, harm, arrayX, numX, 
     *arrayKin, arrayWin)
      real g, s, alp
      real arrayX(numX)
      integer numX, numN, harm(2*numN-1)
      real arrayKin(numX,2*numN-1, 8), arrayWin(numX,2*numN-1, 8)
      real ci, u, ar, nn, uc
      dimension coef(9), ROOT(8), DA(9), DZ(9)
      dimension coef3(4), ROOT3(3), DA3(4), DZ3(4)
      complex ROOT, ROOT3
C
      u = sqrt(1-1.0/g/g)*cos(alp)
C
      do n = 1, 2*numN-1
      do i = 1, numX
         ci = 1.0/sqrt(1+arrayX(i)*arrayX(i))
         uc = u*ci
         nn = harm(n)*1.0/g
C

      if( harm(n) .eq. 0) then
         coef3(4) = uc**4+uc**8-2*uc**6
         coef3(3) = -uc**6+2*uc**4-uc**2-uc**2*s**2+4*uc**4*s**2-3*uc**6
     **s**2
         coef3(2) = -s**2*ci**2*uc**2+uc**4*s**2+3*uc**4*s**4-2*uc**2*s*
     **4+s**2*ci**2-uc**2*s**2
         coef3(1) =-uc**2*s**6
C
      call ZP10R(3,coef3,ROOT3,DA3,DZ3,IERR)

C
      do j=1, 8
        arrayKin(i,n,j) = 0
        arrayWin(i,n,j) = 0
        if( j < 4 ) then
        ar = abs(sqrt(ROOT3(j)))
         if( ar .eq. abs(real(sqrt(ROOT3(j)))) .and. ar .ne. 0 ) then
          if( real(sqrt(ROOT3(j)))*u*ci > 0 ) then
             arrayWin(i,n,j) = real(sqrt(ROOT3(j)))*uc
             arrayKin(i,n,j) = real(sqrt(ROOT3(j)))
          end if
         end if
C
        end if
        if( j>3 .and. j<7 ) then
        ar = abs(-sqrt(ROOT3(j-3)))
         if( ar .eq. abs(real(-sqrt(ROOT3(j-3)))) .and. ar .ne. 0 ) then
          if( real(-sqrt(ROOT3(j-3)))*u*ci > 0 ) then
             arrayWin(i,n,j) = real(-sqrt(ROOT3(j-3)))*uc
             arrayKin(i,n,j) = real(-sqrt(ROOT3(j-3)))
          end if
         end if
C
        end if
C
C
      enddo
C
      else
         coef(9) =  uc**4+uc**8-2*uc**6
         coef(8) = -8*nn*uc**7+12*nn*uc**5-4*nn*uc**3
         coef(7) = 2*uc**4-uc**2+6*nn**2*uc**2+4*uc**4*s**2+28*nn**2*uc*
     **6-uc**2*s**2-3*uc**6*s**2-30*nn**2*uc**4-uc**6
         coef(6) = -4*nn**3*uc-8*nn*uc**3+18*nn*uc**5*s**2+6*nn*uc**5-56
     **nn**3*uc**5+2*nn*uc*s**2+40*nn**3*uc**3+2*nn*uc-16*nn*uc**3*s**2
         coef(5) =12*nn**2*uc**2+nn**4-2*uc**2*s**4-45*nn**2*uc**4*s**2+
     *24*nn**2*uc**2*s**2-30*nn**4*uc**2+3*uc**4*s**4-uc**2*s**2-15*nn**
     *2*uc**4+uc**4*s**2+70*nn**4*uc**4-nn**2-nn**2*s**2+s**2*ci**2-s**2
     **ci**2*uc**2
         coef(4) = 4*nn*uc*s**4-16*nn**3*uc*s**2+2*s**2*ci**2*nn*uc-8*nn
     ***3*uc+2*nn*uc*s**2+60*nn**3*uc**3*s**2+12*nn**5*uc-12*nn*uc**3*s*
     **4+20*nn**3*uc**3-56*nn**5*uc**3-4*nn*uc**3*s**2
         coef(3) = 18*nn**2*uc**2*s**4-uc**2*s**6-s**2*ci**2*nn**2-45*nn
     ***4*uc**2*s**2-2*nn**2*s**4+28*nn**6*uc**2+2*nn**4-15*nn**4*uc**2-
     *2*nn**6+6*nn**2*uc**2*s**2-nn**2*s**2+4*nn**4*s**2
         coef(2) = 2*nn*uc*s**6-12*nn**3*uc*s**4-8*nn**7*uc-4*nn**3*uc*s
     ***2+18*nn**5*uc*s**2+6*nn**5*uc
         coef(1) =-nn**6+3*nn**4*s**4-nn**2*s**6+nn**4*s**2+nn**8-3*nn**
     *6*s**2
C
!      if( harm(n) .eq. 8 ) then
!      print 5, coef(9)/coef(1), coef(8)/coef(1), coef(7)/coef(1), coef(6
!     *)/coef(1), coef(5)/coef(1), coef(4)/coef(1), coef(2)/coef(1)
!    5  format(7E10.3)
!      end if
!      if( n > numN ) then
!      pause
!      end if

      call ZP10R(8,coef,ROOT,DA,DZ,IERR)
C
      do j=1, 8
        arrayKin(i,n,j) = 0
        arrayWin(i,n,j) = 0
        ar = abs(ROOT(j))
C        if( ar .eq. abs(real(ROOT(j))) .and. ar .ne. 0 ) then
        if( ar .eq. abs(real(ROOT(j))) .and. ar .ne. 0 ) then
          if( real(ROOT(j))*u*ci-harm(n)*1.0/g > 0 ) then
             arrayWin(i,n,j) = real(ROOT(j))*u*ci-harm(n)*1.0/g
             arrayKin(i,n,j) = real(ROOT(j))
          end if
        end if
      enddo
C
       end if
      enddo
      enddo
      end
