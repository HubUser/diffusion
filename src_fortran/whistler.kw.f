      subroutine kw(g, s, alp, numN, harm, arrayX, numX, arrayKin, array
     *Win)
      real g, s, alp
      real arrayX(numX)
      integer numX, numN, harm(2*numN-1)
      real arrayKin(numX,2*numN-1, 6), arrayWin(numX,2*numN-1, 6)
      real ci, si2, u, ar, a, dd
      dimension coef(7), ROOT(6), DA(7), DZ(7)
      complex ROOT
C
      u = sqrt(1-1.0/g/g)*cos(alp)
C
      do n = 1, 2*numN-1
      do i = 1, numX
         ci = 1.0/sqrt(1+arrayX(i)*arrayX(i))
         si2 = arrayX(i)*arrayX(i)/(1+arrayX(i)*arrayX(i))
         a = 1.0
C
         coef(7) = u*u*ci*ci*(1+si2/s/s)
         coef(6) = -2*harm(n)*1.0/g*u*ci*(1+si2/s/s)
         coef(5) = harm(n)*harm(n)*1.0/g/g*(1+si2/s/s)+u*u*ci*ci*(a+s*s)
     **(2+si2/s)-ci*ci
         coef(4) = -2*harm(n)*1.0/g*u*ci*(a+s*s)*(2+si2/s)
         coef(3) = harm(n)*harm(n)*1.0/g/g*(a+s*s)*(2+si2/s)+u*u
     **ci*ci*(a+s*s)**2
         coef(2) = -2*harm(n)*1.0/g*u*ci*(a+s*s)**2
         coef(1) = harm(n)*harm(n)*1.0/g/g*(a+s*s)**2
         if( harm(n) .eq. 0) then
           ROOT(1) = 0
           ROOT(2) = 0
           if( coef(5)**2> 4*coef(7)*coef(3) ) then
             dd = sqrt(coef(5)**2-4*coef(7)*coef(3))
             if(dd-coef(5)> 0 ) then
             ROOT(1) = sqrt( (dd-coef(5))/2.0/coef(7) )
             end if
             if(-dd-coef(5)> 0 ) then
             ROOT(2) = sqrt( (-dd-coef(5))/2.0/coef(7) )
             end if
           end if
C
           do j=3, 6
             ROOT(j) = 0
           end do
         else
C      print 6, coef
C    6 FORMAT(7E15.7)
C      pause
          call ZP10R(6,coef,ROOT,DA,DZ,IERR)
         end if
C

C      print 5, ROOT
C    5 FORMAT(6E15.7)
C      pause
C      print 6, n
C    6 FORMAT(I5)
C      print 5, ROOT
C    5 FORMAT(4E15.7)
C
C
      do j=1, 6
        arrayKin(i,n,j) = 0
        arrayWin(i,n,j) = 0
        ar = abs(ROOT(j))
C        if( ar .eq. abs(real(ROOT(j))) .and. ar .ne. 0 ) then
        if( ar .eq. abs(real(ROOT(j))) .and. ar .ne. 0 ) then
          if( real(ROOT(j))*u*ci-harm(n)*1.0/g > 0 ) then
             arrayWin(i,n,j) = real(ROOT(j))*u*ci-harm(n)*1.0/g
             arrayKin(i,n,j) = real(ROOT(j))
          end if
C             arrayKin(i,n,j) = real(ROOT(j))
C             arrayWin(i,n,j) = arrayKin(i,n,j)*fgam*ci*ca-harm(n)*1.0/g
C             print 2, arrayKin(i,n,j), arrayWin(i,n,j)
C    2        FORMAT(4E15.7)
        end if
      enddo
C
      enddo
      enddo
C
C
      end