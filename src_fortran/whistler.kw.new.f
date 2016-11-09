      subroutine kwNew(g, s, alp, numN, harm, arrayX, numX, arrayKin, ar
     *rayWin)
      real g, s, alp
      real arrayX(numX)
      integer numX, numN, harm(2*numN-1)
      real arrayKin(numX,2*numN-1, 8), arrayWin(numX,2*numN-1, 8)
      real ci, u, ar, dd, nn
      dimension coef(4), ROOT(3), DA(4), DZ(4)
      complex ROOT
C
      
C
      do n = 1, 2*numN-1
      do i = 1, numX
         ci = 1.0/sqrt(1+arrayX(i)*arrayX(i))
         u = sqrt(1-1.0/g/g)*cos(alp)*ci
         nn = harm(n)/g;
C
         coef(4) = u
         coef(3) = -(ci+nn)
         coef(2) = u*s*s
         coef(1) = -nn*s*s
C
!      u*s^2-(nn+co)*k+u*k^2
         if( nn .eq. 0) then
           ROOT(1) = 0
           ROOT(2) = 0
           ROOT(3) = 0
           if( coef(3)**2> 4*coef(4)*coef(2) ) then
             dd = sqrt(coef(3)**2-4*coef(4)*coef(2))
             if(dd-coef(3)> 0 ) then
             ROOT(1) = (dd-coef(3))/2.0/coef(4)
             end if
             if(-dd-coef(3)> 0 ) then
             ROOT(2) = (-dd-coef(3))/2.0/coef(4)
             end if
           end if
         else
C      print 6, coef
C    6 FORMAT(7E15.7)
C      pause
          call ZP10R(3,coef,ROOT,DA,DZ,IERR)
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
      do j=1, 8
        arrayKin(i,n,j) = 0
        arrayWin(i,n,j) = 0
        if( j < 4) then
          ar = abs(ROOT(j))
          if( ar .eq. abs(real(ROOT(j))) .and. ar .ne. 0 ) then
            if( real(ROOT(j))*u-nn > 0 ) then
               arrayWin(i,n,j) = real(ROOT(j))*u-nn
               arrayKin(i,n,j) = real(ROOT(j))
            end if
          end if
        end if
      enddo  ! j
C
      enddo
      enddo
C
C
      end