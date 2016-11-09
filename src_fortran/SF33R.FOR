      SUBROUTINE SF33R(X,N,NMAX,FUN1,FUN2,WORK,IERR)
      DIMENSION FUN1(NMAX),FUN2(NMAX),WORK(NMAX)
      INTEGER N,NMAX,N1,N2,NA,NB,NC1,NC2,I,J,K,INT,
     1MAX0,MAX1
      REAL X,FUN1,FUN2,WORK,P,P0,P1,Q0,Q1,Q2,R,R1,R2,R3,R4,
     1R5,U,U1,U2,U3,U4,U5,V,V1,V2,V3,V4,V5,V6,V7,SYS039,
     1SYS015,SYS005,ALOG,SQRT,COS,SIN
      DATA  SYS039/0.577215665/,SYS015/0.63661977/,SYS005/
     11.57079633/
      IF(N.LT.2) GO TO 26
      NC1=MAX0(N,2*INT(X)+4)
      IF(NMAX.LT.NC1) GO TO 28
      IERR=0
      IF(X) 1,1,2
    1 IERR=67
      GO TO 27
    2 NC1=MAX1(X+1.,3.)
      NC2=2*NC1+1
      IF(N-NC1) 3,3,4
    3 N2=N
      N1=N
      GO TO 7
    4 N1=NC1
      IF(N-NC2) 5,5,6
    5 N2=NC2+1
      GO TO 7
    6 N2=N
    7 IF(X.GE.10.) GO TO 14
      IF(X.GT.1.) GO TO 12
      R1=0.5*X
      R2=-R1**2
      U1=1.
      U2=1.
      U3=1.
      U4=1.
      P0=1.
      P1=1.
      Q0=-(ALOG(R1)+SYS039)
      V1=Q0
      Q1=V1+0.5
      V2=Q1
      DO 8 I=1,8
      U1=U1*R2/FLOAT((N-1+I)*I)
      P0=P0+U1
      U2=U2*R2/FLOAT((N-2+I)*I)
      P1=P1+U2
      R3=1.E0/FLOAT(I)
      R4=1.E0/FLOAT(I+1)
      U3=U3*R2*R3**2
      V1=V1+R3
      Q0=Q0+U3*V1
      U4=U4*R2*R3*R4
      V2=V2+0.5*(R3+R4)
      Q1=Q1+U4*V2
    8 CONTINUE
      FUN1(N)=P0
      FUN1(N-1)=P1
      R3=-SYS015
      FUN2(1)=R3*Q0
      FUN2(2)=R3*(1.E0-2.E0*R2*Q1)
      IF(N.EQ.2) GO TO 11
      FUN2(3)=FUN2(2)+R2*FUN2(1)*2
      DO 9 I=4,N
      FUN2(I)=FUN2(I-1)+FUN2(I-2)*R2/FLOAT((I-2)*(I-3))
    9 CONTINUE
      NA=N-2
      DO 10 J=1,NA
      I=NA-J+1
      FUN1(I)=FUN1(I+1)+FUN1(I+2)*R2/FLOAT(I*(I+1))
   10 CONTINUE
   11 GO TO 25
   12 R1=0.5*X
      R2=-R1**2
      U1=1.
      U2=1.
      P0=1.
      P1=1.
      Q0=-(ALOG(R1)+SYS039)
      V1=Q0
      Q1=V1+0.5
      V2=Q1
      DO 13 I=1,20
      R3=1.E0/FLOAT(I)
      R4=1.E0/FLOAT(I+1)
      U1=U1*R2*R3**2
      P0=P0+U1
      V1=V1+R3
      Q0=Q0+U1*V1
      U2=U2*R2*R3*R4
      P1=P1+U2
      V2=V2+0.5*(R3+R4)
      Q1=Q1+U2*V2
   13 CONTINUE
      FUN1(1)=P0
      FUN1(2)=R1*P1
      R2=-SYS015
      FUN2(1)=R2*Q0
      FUN2(2)=(0.5/R1+R1*Q1)*R2
      GO TO 16
   14 R1=0.125/X
      R2=-R1**2
      P0=1.
      Q0=1.
      P1=1.
      Q1=1.
      U1=1.
      U2=1.
      U3=1.
      U4=1.
      DO 15 I=1,9
      V1=FLOAT(2*I)
      V2=2.E0*V1
      V3=((V2-3.)*(V2-1.))**2
      V4=V2**2-1.
      V5=V4*(V2-3.)
      V6=V5*(V2+3.)
      V5=V5*(V2-5.)
      V4=V4**2
      V7=R2/((V1-1.)*V1)
      V2=R2/((V1+1.)*V1)
      U1=U1*V3*V7
      P0=P0+U1
      U2=U2*V4*V2
      Q0=Q0+U2
      U3=U3*V5*V7
      P1=P1+U3
      U4=U4*V6*V2
      Q1=Q1+U4
   15 CONTINUE
      Q0=R1*Q0
      Q1=3.E0*R1*Q1
      U1=SYS005
      U2=SQRT(U1*X)
      U3=X+0.5*U1
      U4=SIN(U3)/U2
      U5=COS(U3)/U2
      FUN1(1)=P0*U4-Q0*U5
      FUN1(2)=Q1*U4-P1*U5
      FUN2(1)=-(P0*U5+Q0*U4)
      FUN2(2)=-(P1*U4+Q1*U5)
   16 R1=2.E0/X
      DO 17 I=3,N1
      R2=R1*FLOAT(I-2)
      FUN1(I)=R2*FUN1(I-1)-FUN1(I-2)
      FUN2(I)=R2*FUN2(I-1)-FUN2(I-2)
   17 CONTINUE
      IF(N1.EQ.N) GO TO 25
      R1=0.5*X
      R2=R1**2
      Q2=1.
      U=1.
      V=1.
      R3=R1
      NA=N2-2
      DO 18 I=1,NA
      V=V+1.E0/FLOAT(I+1)
      U=U*R2/FLOAT(I*(NA+1-I))
      Q2=Q2+U
      R3=R3*R1/FLOAT(I)
   18 CONTINUE
      U=1.
      P=1.
      R2=-R2
      IF(X-10.) 19,20,20
   19 K=20
      GO TO 21
   20 K=NC2-1
   21 Q1=.5*V-ALOG(R1)-SYS039
      DO 22 I=1,K
      R4=1.E0/FLOAT(I)
      R5=1.E0/FLOAT(N2-1+I)
      U=U*R2*R4*R5
      P=P+U
      V=V+0.5*(R4+R5)
      Q1=Q1+U*V
   22 CONTINUE
      Q2=0.5*Q2/R3
      R3=R3/FLOAT(N2-1)
      P=P*R3
      Q1=Q1*R3
      R2=-SYS015
      FUN1(N2)=P
      FUN2(N2)=R2*(Q1+Q2)
      R=2./X
      WORK(N1)=0.E0
      NA=N1+1
      NB=N2-1
      DO 23 I=NA,NB
      WORK(I)=1.E0/(R*FLOAT(I-1)-WORK(I-1))
      FUN1(I)=WORK(I)*FUN1(I-1)
      FUN2(I)=WORK(I)*FUN2(I-1)
   23 CONTINUE
      DO 24 J=NA,NB
      I=NA+NB-J
      FUN1(I)=WORK(I)*FUN1(I+1)+FUN1(I)
      FUN2(I)=WORK(I)*FUN2(I+1)+FUN2(I)
   24 CONTINUE
   25 RETURN
   26 IERR=65
   27 CALL UTSF10(IERR,33)
      GO TO 25
   28 IERR=66
      GO TO 27
      END
