      include 'UTZP10.FOR'
      include 'ZP10R.FOR'
      include 'ZP14C.FOR'
      include 'ZP14R.FOR'
      include 'SF34R.FOR'
C     bessel function
      include 'SF33R.FOR'
      include 'UTSF10.FOR'
      include 'SF24R.FOR'
      include 'SF25R.FOR'
      include 'specfun.FOR'
C     subroutine
      ! include 'whistler.kw.f'
      ! include 'whistler.kw.old.f'
      ! include 'whistler.kw.new.f'
      include 'whistler.kw.full.f'
      include 'besselJ.f'
      include 'besselJ_new.f'
      include 'functionPhi.f'
      include 'functionQ.f'
      ! include 'functionQOld.f'
      include 'functionQNew.f'
      ! include 'functionN.f'
      include 'functionNOld.f'
      include 'functionNNew.f'
      include 'functionNTest.f'
      include 'functionNTet2.fast.f'
      include 'functionNTet2.Lyons.Test.f'
      include 'functionNTet2.Lyons.f'
      include 'localDFull.f'
      include 'localDTet2Dst.longTet.New.Fast.f'
      include 'localDMSW.f'

      recursive subroutine radbelInputTet2DstDoubleQ(gam, s, nameS, 
     *nameG, lM, alphaL, kp)
C> input parameter string containing energy
      character(len=6)::nameG
C> input parameter string containing wpe/wce
      character(len=4)::nameS
C> input parameter variable calculated from particle energy in keV
      real gam
C> input parameter ratio wpe/wce
      real s 
C> input parameter K_p index
      real kp
C> input parameter 
      real alphaL, lM

      real alp
      character(len=5)::nA
      character(len=4)::nameB
      integer numN
      parameter(numN=61)
C
      integer nalp, alphaMax
      parameter(nalp=100, alphaMax=90.0)
      real alpEQ(nalp)
C
      real coefDaaLoc(2*numN+5), coefDeeB(2*numN+5),
     *coefDaaB(2*numN+5), coefDaaSB(2*numN+5), coefDeeSB(2*numN+5)
      real coefDeaLoc(2*numN+5), coefDeeLoc(2*numN+5)
C
      real dl, Tl, lam, ldeg, b, valueWm
      parameter(dl=0.01)
C
      real valueTm(2), valueTw(2)
      real lM0, valueTm1, valueTm01, a1Tm1, a2Tm1, a3Tm1, a4Tm1, a5T
     *m1, valueTm2, valueTm02, a1Tm2, a2Tm2, a3Tm2, a4Tm2, a5Tm2
      parameter(lM0=3.14159)
C
      parameter(valueTm01=11.5372, a1Tm1=14.3157, a2Tm1=-8.18009,
     * a3Tm1=1.38532, a4Tm1=0.0, a5Tm1=0.0)  !L 4-5
      parameter(valueTm02=66.0, a1Tm2=1.0, a2Tm2=0.0,
     *a3Tm2=0.0, a4Tm2=0.0, a5Tm2=0.0)  !L 4-5
C
C
      real valueTw1, valueTw01, a1Tw1, a2Tw1, a3Tw1, a4Tw1, a5Tw1, valu
     *eTw2, valueTw02, a1Tw2, a2Tw2, a3Tw2, a4Tw2, a5Tw2
C
      parameter(valueTw01=5.68019, a1Tw1=4.62479, a2Tw1=3.04756, a3T
     *w1=-5.05697, a4Tw1=1.83305, a5Tw1=-2.02501e-001)  !L 4-5
      parameter(valueTw02=5.68019, a1Tw2=4.62479, a2Tw2=3.04756,
     *a3Tw2=-5.05697, a4Tw2=1.83305, a5Tw2=-2.02501e-001)  !L 4-5
C
C
      real factorQ, Bw, Bw15, coefA(4,4), coefB(4,4)
C
      real norm, L, ne0
!      parameter (L=4.5, norm = 330.0/(L/4.5)**3)
!      parameter (L=5.5, norm = 330.0/(L/4.5)**3)
      parameter (L=5.0, norm = 184.0/(L/4.5)**3)     !L-shell
C
      integer now(3), nownew(3), dTime
C
      write (Unit=nA, FMT="(F5.2)") alphaL
      write (Unit=nameB, FMT="(F4.0)") kp
C
      coefB(1,1) = 0.002280709100887
      coefB(1,2) = 0.000007740005458
      coefB(1,3) = 0.000000043666208
      coefB(1,4) = 0.000000011799739
C
      coefB(2,1) = -0.000006133966053
      coefB(2,2) = 0.000004421397080
      coefB(2,3) = 0.000000357794391
      coefB(2,4) = -0.000000001681957
C
      coefB(3,1) = 0.000003074870165
      coefB(3,2) = 0.000000691596767
      coefB(3,3) = -0.000000052153258
      coefB(3,4) = 0.000000000679970
C
      coefB(4,1) = -0.000000065216319
      coefB(4,2) = -0.000000019691972
      coefB(4,3) = 0.000000001242614
      coefB(4,4) = -0.000000000019907
C
C
C
C
      coefA(1,1) = -2.962275028228760
      coefA(1,2) = 0.005432456266135
      coefA(1,3) = 0.000125454942463
      coefA(1,4) = -0.000001982993581
C
      coefA(2,1) =  -0.022489553317428
      coefA(2,2) = 0.025674106553197
      coefA(2,3) = -0.000471087114420
      coefA(2,4) = 0.000001404062914
C
      coefA(3,1) = 0.004958000034094
      coefA(3,2) = -0.001571639673784
      coefA(3,3) = 0.000026058032745
      coefA(3,4) = -0.000000054605785
C
      coefA(4,1) = -0.000074383700849
      coefA(4,2) = 0.000023944545319
      coefA(4,3) = -0.000000390142418
      coefA(4,4) = 0.000000000685272

C
      open (unit = 7, file = nameG//'_ee_mag('//nameS//')_'//nA//'_Kp'//
     *nameB//'.dat')
      write(7,*)  'alp n0 n1 n2 n3 n4 n5 n6 n7 n8 n9 n10 n11 n12 n13 n14
     * n15 n16 n17 n18 n19 n20 n21 n22 n23 n24 n25 n26 n27 n28 n29 n30 n
     *31 n32 n33 n34 n35 n36 n37 n38 n39 n40 n41 n42 n43 n44 n45 n46 n47
     * n48 n49 n50 n51 n52 n53 n54 n55 n56 n57 n58 n59 n60 n-1 n-2 n-3 n
     *-4 n-5 n-6 n-7 n-8 n-9 n-10 n-11 n-12 n-13 n-14 n-15 n-16 n-17 n-1
     *8 n-19 n-20 n-21 n-22 n-23 n-24 n-25 n-26 n-27 n-28 n-29 n-30 n-31
     * n-32 n-33 n-34 n-35 n-36 n-37 n-38 n-39 n-40 n-41 n-42 n-43 n-44
     *n-45 n-46 n-47 n-48 n-49 n-50 n-51 n-52 n-53 n-54 n-55 n-56 n-57 n
     *-58 n-59 n-60 sum sum5 sum10 sum20 sum40'
C
      open (unit = 8, file = nameG//'_aa_mag('//nameS//')_'//nA//'_Kp'//
     *nameB//'.dat')
      write(8,*)  'alp n0 n1 n2 n3 n4 n5 n6 n7 n8 n9 n10 n11 n12 n13 n14
     * n15 n16 n17 n18 n19 n20 n21 n22 n23 n24 n25 n26 n27 n28 n29 n30 n
     *31 n32 n33 n34 n35 n36 n37 n38 n39 n40 n41 n42 n43 n44 n45 n46 n47
     * n48 n49 n50 n51 n52 n53 n54 n55 n56 n57 n58 n59 n60 n-1 n-2 n-3 n
     *-4 n-5 n-6 n-7 n-8 n-9 n-10 n-11 n-12 n-13 n-14 n-15 n-16 n-17 n-1
     *8 n-19 n-20 n-21 n-22 n-23 n-24 n-25 n-26 n-27 n-28 n-29 n-30 n-31
     * n-32 n-33 n-34 n-35 n-36 n-37 n-38 n-39 n-40 n-41 n-42 n-43 n-44
     *n-45 n-46 n-47 n-48 n-49 n-50 n-51 n-52 n-53 n-54 n-55 n-56 n-57 n
     *-58 n-59 n-60 sum sum5 sum10 sum20 sum40'
C
      open (unit = 9, file = nameG//'_aa_mgs('//nameS//')_'//nA//'_Kp'//
     *nameB//'.dat')
      write(9,*)  'alp n0 n1 n2 n3 n4 n5 n6 n7 n8 n9 n10 n11 n12 n13 n14
     * n15 n16 n17 n18 n19 n20 n21 n22 n23 n24 n25 n26 n27 n28 n29 n30 n
     *31 n32 n33 n34 n35 n36 n37 n38 n39 n40 n41 n42 n43 n44 n45 n46 n47
     * n48 n49 n50 n51 n52 n53 n54 n55 n56 n57 n58 n59 n60 n-1 n-2 n-3 n
     *-4 n-5 n-6 n-7 n-8 n-9 n-10 n-11 n-12 n-13 n-14 n-15 n-16 n-17 n-1
     *8 n-19 n-20 n-21 n-22 n-23 n-24 n-25 n-26 n-27 n-28 n-29 n-30 n-31
     * n-32 n-33 n-34 n-35 n-36 n-37 n-38 n-39 n-40 n-41 n-42 n-43 n-44
     *n-45 n-46 n-47 n-48 n-49 n-50 n-51 n-52 n-53 n-54 n-55 n-56 n-57 n
     *-58 n-59 n-60 sum sum5 sum10 sum20 sum40'
C
      open (unit =10, file = nameG//'_ee_mgs('//nameS//')_'//nA//'_Kp'//
     *nameB//'.dat')
      write(10,*) 'alp n0 n1 n2 n3 n4 n5 n6 n7 n8 n9 n10 n11 n12 n13 n14
     * n15 n16 n17 n18 n19 n20 n21 n22 n23 n24 n25 n26 n27 n28 n29 n30 n
     *31 n32 n33 n34 n35 n36 n37 n38 n39 n40 n41 n42 n43 n44 n45 n46 n47
     * n48 n49 n50 n51 n52 n53 n54 n55 n56 n57 n58 n59 n60 n-1 n-2 n-3 n
     *-4 n-5 n-6 n-7 n-8 n-9 n-10 n-11 n-12 n-13 n-14 n-15 n-16 n-17 n-1
     *8 n-19 n-20 n-21 n-22 n-23 n-24 n-25 n-26 n-27 n-28 n-29 n-30 n-31
     * n-32 n-33 n-34 n-35 n-36 n-37 n-38 n-39 n-40 n-41 n-42 n-43 n-44
     *n-45 n-46 n-47 n-48 n-49 n-50 n-51 n-52 n-53 n-54 n-55 n-56 n-57 n
     *-58 n-59 n-60 sum sum5 sum10 sum20 sum40'
C
      do ia = 1, nalp
        call itime(now)
C        alpEQ(ia) = (5+84*(ia-1.0)/(nalp-1.0) )*3.14159/180.0
C nalp point in range [5, 90] but in radians
        alpEQ(ia) = (5.0+(alphaMax-5.0)*(ia-1.0)/(nalp-1.0) )*3.14159/18
     *0.0
        Tl = 1.3 - 0.56*sin(alpEQ(ia))
C
        coefDeeB(1)=alpEQ(ia)*180.0/3.14159
        coefDaaB(1)=alpEQ(ia)*180.0/3.14159
        coefDaaSB(1)=alpEQ(ia)*180.0/3.14159
        coefDeeSB(1)=alpEQ(ia)*180.0/3.14159
        do i=2,2*numN+5
          coefDeeSB(i) = 0
          coefDeeB(i) = 0
          coefDaaB(i) = 0
          coefDaaSB(i) = 0
        enddo
C
      lam = 0.0
      alp = alpha( lam, alpEQ(ia) )
      do while (alp > 0)
C         print 5,  alp, lam
      b = sqrt(1+3*sin(lam)**2)/cos(lam)**6
      ne0 = cos(lam)**(-2*alphaL)

      ldeg = lam*180/3.14159/10.0
C
        valueTm1=(valueTm01+a1Tm1*ldeg+a2Tm1*ldeg**2+a3Tm1*ldeg**3+a4Tm1
     **ldeg**4+a5Tm1*ldeg**5)*3.14159/180.0
        valueTm2=(valueTm02+a1Tm2*ldeg+a2Tm2*ldeg**2+a3Tm2*ldeg**3+a4Tm2
     **ldeg**4+a5Tm2*ldeg**5)*3.14159/180.0
        valueTw1=sqrt(2.0)*(valueTw01+a1Tw1*ldeg+a2Tw1*ldeg**2+a3Tw1*lde
     *g**3+a4Tw1*ldeg**4+a5Tw1*ldeg**5)*3.14159/180.0
        valueTw2=sqrt(2.0)*(valueTw02+a1Tw2*ldeg+a2Tw2*ldeg**2+a3Tw2*lde
     *g**3+a4Tw2*ldeg**4+a5Tw2*ldeg**5)*3.14159/180.0
C
      Bw = 0
      do j1 = 1, 4
        do i1 = 1, 4
         Bw = Bw + coefB(j1, i1)*(kp)**(i1-1)*(ldeg*10)**(j1-1)
        end do
      end do
      Bw = Bw*1000.0
C
      if( kp>50 .and. ldeg>3.0 .and. Bw<3.0) then
      Bw = 3.0
      end if
C
      if( Bw<1.0 ) then
      Bw = 1.0
      end if
!     Bw = 100.0
!      print 5,  Bw
!      pause
C
      factorQ = 0.0
      if( kp<65 ) then
      do j1 = 1, 4
        do i1 = 1, 4
         factorQ = factorQ + coefA(j1, i1)*(kp)**(i1-1)
     **(ldeg*10)**(j1-1)
        end do
      end do
C      factorQ = 0.2+1.0*log10(Bw)
      factorQ = 10.0**(factorQ)
      if( factorQ<0.005 ) then
        factorQ = 0.0
      end if
      end if

      valueWm = 0.35!  wave frequency

C      print 5, factorQ, Bw
C      pause
      if( lam < lM ) then
C
!      call localD(gam, s/b, alp, numN, coefDaaLoc, coefDeaLoc, coefDeeLo
!     *c, valueXm, valueXw, 20.0, b)
C
C       print 5, alp, lam, ldeg
      valueTm(1) = valueTm1
      valueTm(2) = valueTm2
      valueTw(1) = valueTw1
      valueTw(2) = valueTw2
C
C      call localDTet2Dst(gam, s*sqrt(ne0)/b, alp, numN, coefDaaLoc, coef
C     *DeaLoc, coefDeeLoc, valueTm, valueTw, 1.0, factorQ, 1.521, b,
C     * 0.35)
     
      call localDTet2Dst(gam, s*sqrt(ne0)/b, alp, numN, coefDaaLoc, coef
     *DeaLoc, coefDeeLoc, valueTm, valueTw, 1.0, factorQ, 1.571, b,
     * valueWm, numT)
C
         do i=2,2*numN+5
            coefDaaB(i) = coefDaaB(i) + coefDaaLoc(i)*cos(alp)/cos(alpEQ
     *(ia))**2*cos(lam)**7*dl/Tl/norm*Bw*Bw/100.0/100.0
            coefDeeB(i) = coefDeeB(i) + coefDeeLoc(i)/cos(alp)*sqrt(1+3*
     *sin(lam)**3)*cos(lam)*dl/Tl/norm*Bw*Bw/100.0/100.0
         enddo
C
      end if
C
      if( lam < lM0 ) then
C
C       print 5, alp, lam, ldeg
      valueTm(1) = 0.0
      valueTm(2) = 0.0
      valueTw(1) = 0.523
      valueTw(2) = 1.0
C
!      call localDTet2Dst(gam, s*sqrt(ne0)/b, alp, numN, coefDaaLoc, coef
!     *DeaLoc, coefDeeLoc, valueTm, valueTw, 1.0, 0.0, 0.785, b, 0.35)
C
      call localDTet2Dst(gam, s*sqrt(ne0)/b, alp, numN, coefDaaLoc, coef
     *DeaLoc, coefDeeLoc, valueTm, valueTw, 1.0, 0.0, 0.785, b,
     * valueWm, numT)
C
         do i=2,2*numN+5
          coefDaaSB(i) = coefDaaSB(i) + coefDaaLoc(i)*cos(alp)/cos(alpEQ
     *(ia))**2*cos(lam)**7*dl/Tl/norm*Bw*Bw/100.0/100.0
          coefDeeSB(i) = coefDeeSB(i) + coefDeeLoc(i)/cos(alp)*sqrt(1+3*
     *sin(lam)**3)*cos(lam)*dl/Tl/norm*Bw*Bw/100.0/100.0
         enddo
      end if
C
         lam = lam + dl
         alp = alpha( lam, alpEQ(ia) )
      end do
C
C      write (7,5) alpEQ(ia)*180/3.14159, coefDaa
       write(7,6)  (coefDeeB(n), n=1,2*numN+5 )
       write(8,6)  (coefDaaB(n), n=1,2*numN+5 )
       write(9,6)  (coefDaaSB(n), n=1,2*numN+5 )
       write(10,6)  (coefDeeSB(n), n=1,2*numN+5 )
    6  format(308E15.7)
C
      call itime(nownew)
      dTime = (nownew(1)-now(1))*3600 + (nownew(2)-now(2))*60 + (nownew(
     *3)-now(3))
C
!       print *,'!'
       print 5,  coefDaaB(1), coefDaaB(2*numN+1), coefDaaSB(2*numN+1),
     *coefDeeB(2*numN+1), dTime*(nalp-ia)/60.0
    5  format(5E15.7)
      enddo  ! end alpha
      close(7)
      close(8)
      close(9)
      close(10)
C
C
      end
C
C
C     alpha
      recursive function alpha( lam, alpEQ )
      real lam, alpEQ
      real dipole, fun
      real alpha
C
      dipole = cos(lam)**6/sqrt(1+3*sin(lam)**2)
      dipole = 1.0/dipole
      fun = sin(alpEQ)*sqrt(dipole)
      alpha = -10
      if ( fun<1 ) then
        alpha = asin( fun )
      end if
      return
      end
C
C     lambdaMax
      subroutine lambdaMax( alpEQ, lambdaM )
      real alpEQ
      real lambdaM
      dimension coef(7), ROOT(6), DA(7), DZ(7)
      complex ROOT
C
      coef(7) = 1.0
      coef(6) = 0 !5
      coef(5) = 0 !4
      coef(4) = 0 !3
      coef(3) = 0.0001 !2
      coef(2) = 3*sin(alpEQ)**4 !1
      coef(1) = -4*sin(alpEQ)**4
C
      call ZP10R(6,coef,ROOT,DA,DZ,IERR)
      do j = 1, 6
        if( abs(ROOT(j)) .eq. real(ROOT(j)) .and. real(ROOT(j))>0 .and.
     *abs(ROOT(j))<=1 ) then
           lambdaM = acos(sqrt(real(ROOT(j))))
        end if
      end do
      end
