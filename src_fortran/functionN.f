      subroutine funN(s, numX, arrayX, arrayG, valueKin, valueN)
      real s, valueKin
      integer numX
      real arrayX(numX), arrayG(numX)
      real valueN
      real dX, zn, c
C
      c = valueKin*(valueKin**2+s*s+1)**2/(2*3.14159)**2/(1+s*s)
      c = abs(c)
      valueN = 0
      do i = 1, numX
         if( i < numX ) then
            dX = arrayX(i+1)-arrayX(i)
         else
            dX = arrayX(numX)-arrayX(numX-1)
         end if
         zn = 1+arrayX(i)*arrayX(i)
         valueN = valueN + c*arrayG(i)*arrayX(i)*dX/zn
      enddo
C
      end