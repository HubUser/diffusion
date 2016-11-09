      subroutine funQ(g, s, alp, valueX, valueKin, valueWin, valueQ)
      real g, s, alp, valueX
      real valueKin, valueWin
      real valueQ, f, df, si2, co2
C
C      si2 = valueX**2/(1+valueX**2)
C      f = (valueKin*valueKin+1+s*s)**2+valueKin*valueKin/s*(valueKin*val
C     *ueKin/s+1+s*s)*si2
C      df = 2*(valueKin*valueKin+1+s*s)+1.0/s*(2*valueKin*valueKin/s+1+s*
C     *s)*si2
C
C      valueQ = sqrt(1-1.0/g/g)*cos(alp)-valueKin/sqrt(f)*(2-valueKin*val
C     *ueKin*df/f)
!      valueQ = sqrt(1-1.0/g/g)*cos(alp)-2*valueKin*s*s/(valueKin*valueKi
!     *n+s*s)/(valueKin*valueKin+s*s)
      co2 = 1.0/(1+valueX*valueX)
!      valueQ = sqrt(1-1.0/g/g)*cos(alp)-valueKin/(valueKin*valueKi
!     *n+s*s)/(valueKin*valueKin+s*s)*(valueKin**2*(1-2*co2)+s*s+co2)
C
      valueQ = sqrt(1-1.0/g/g)*cos(alp)-2*(s*s+1)/valueKin**3/(1+(s*s+1)
     */valueKin**2)**2
C
!      valueQ = sqrt(1-1.0/g/g)*cos(alp)-2*valueKin/valueWin
C
      valueQ = abs(valueQ)
      if( valueQ <0.01 ) then
        valueQ = 0
      else
        valueQ = 1.0/valueQ
      end if
      end