      SUBROUTINE UTSF10(IERR,N)
      IF(N.EQ.28.AND.IERR.EQ.1)GO TO 28
      IF(N.EQ.29.AND.IERR.EQ.1)GO TO 29
      I=IERR-64
      IF(N.EQ.10)GOTO 10
      IF(N.EQ.11)GO TO 11
      IF(N.EQ.13) GOTO 13
      IF(N.EQ.14) GO TO 14
      IF(N.EQ.15) GO TO 15
      IF (N. EQ. 16) GO TO 16
      IF (N .EQ. 17) GO TO 17
      IF (N .EQ. 18) GO TO 18
      IF (N.EQ.19) GO TO 19
      IF(N.EQ.20)GOTO 200
      IF(N.EQ.21)GOTO 21
      IF(N.EQ.22)GO TO 22
      IF(N.EQ.23)GO TO 23
      IF(N.EQ.28)GOTO 282
      IF(N.EQ.29)GO TO 292
      IF(N.EQ.31) GO TO 31
      IF (N .EQ. 32) GOTO 320
      IF(N.EQ.33)GO TO 330
      IF(N.EQ.34) GO TO 340
      IF(N.EQ.35) GO TO 350
      RETURN
   10 PRINT 101
      GO TO (105,106,107),I
  105 PRINT102
      RETURN
  106 PRINT 103
      RETURN
  107 PRINT 104
      RETURN
   11 PRINT 111
      GO TO (105,106,115),I
  115 PRINT 112
      RETURN
   13 IF (IERR .NE. 1) GOTO 135
      PRINT 131
  131 FORMAT(' OTEKA HB MY,YHK SF13R:HEATAHA OKA')
      PRINT 132
  132 FORMAT(' N 1 - APYMEHT AAH OTPATEHM HAEHEM '/
     *'       L=2,BCEHE POOEH C COOBAHEM ACOTHOO',
     *' HAEH APYMEHTA')
      GOTO 9005
  135 I=IERR-65
      PRINT 136
  136 FORMAT(' OTEKA HB MY,YHK SF13R:ATAHA OKA')
      GOTO(137,138,139,1391),I
  137 PRINT 1392
 1392 FORMAT(' N 2 - HAEHE L MEHE 1  OE 3.'/
     *'       HAEHE SF13R OOEHO PABHM 3.4E38')
      GOTO 9005
  138 PRINT 1393
 1393 FORMAT(' N 3 - HAEHE APYMEHTA PABHO 0.HAEHE SF13R OOE',
     *'HO PABHM'/7X,'3.4E38,EC',
     *' L=2, -3.4E38,EC L=1  3')
      GOTO 9005
  139 PRINT 1394
 1394 FORMAT(' N 4 - P L=1 HAEHE APYMEHTA PEBOCXOT MAKCMA',
     *'HO OYCTMOE 43.66.'/
     *'       HAEHE SF13R OOEHO PABHM 3.4E38')
      GOTO 9005
 1391 PRINT 1395
 1395 FORMAT(' N 5 - HAEHE APYMEHTA MEHE -4.118E01 P L=1  ',
     *'2.'/
     *'       HAEHE SF13R OAAETC PABHM 0.')
      GO TO 9005
   14 PRINT141
      GO TO(142,143,144,145,146),I
  142 PRINT 147
      RETURN
  143 PRINT148
      RETURN
  144 PRINT 149
      RETURN
  145 PRINT 1490
      RETURN
  146 PRINT 1491
      RETURN
   15 PRINT 151
      GOTO (142,143,152),I
  152 PRINT 153
      RETURN
   16 IF (IERR .NE. 1) GO TO 161
C  162 PRINT 163
      PRINT 164
      GOTO 9005
  161 PRINT 165
      PRINT 166
      GOTO 9005
   17 IF (IERR .NE. 1) GOTO 171
C  172 PRINT 173
      PRINT 174
      GOTO 9005
  171 PRINT 175
      GOTO 9005
   18 IF (IERR .NE. 1) GO TO 181
C  182 PRINT 183
      PRINT 184
      GOTO 9005
  181 PRINT 185
      PRINT 186
 9005 RETURN
   19 PRINT 191
      GOTO(195,196,197),I
  195 PRINT 192
      RETURN
  196 PRINT 193
      RETURN
  197 PRINT 194
      RETURN
  200 PRINT 201
      GOTO(202,203),I
  202 PRINT192
      RETURN
  203 PRINT 193
      RETURN
   21 PRINT 210
      GOTO(211,212,211),I
  211 PRINT 213
      RETURN
  212 PRINT 214
      RETURN
   22 PRINT 220
      GOTO(221,224),I
  224 PRINT 222
      RETURN
  221 PRINT 223
      RETURN
   23 PRINT 230
      GO TO (231,335,232),I
  231 PRINT 233
      RETURN
  232 PRINT 234
      RETURN
   28 PRINT 281
      RETURN
   29 PRINT 291
      RETURN
  282 PRINT 283
      RETURN
  292 PRINT 293
      RETURN
   31 PRINT 311
      RETURN
  320 PRINT 321
      GO TO (322,323),I
  322 PRINT 324
      GOTO 326
  323 PRINT 325
  326 RETURN
  330 PRINT 331
      GO TO(332,333,334,335),I
  332 PRINT 338
      RETURN
  333 PRINT 336
      RETURN
  334 PRINT 337
      RETURN
  335 PRINT 339
      RETURN
  101 FORMAT (' OTEKA HB MY YHK SF10R: ATAHA OKA')
  111 FORMAT (' OTEKA HB MY YHK SF11R: ATAHA OKA')
  141 FORMAT (' OTEKA HB MY YHK SF14R: ATAHA OKA')
  151 FORMAT (' OTEKA HB MY YHK SF15R: ATAHA OKA')
  102 FORMAT(' N 1 - YCOBHOE CO L HE PABHO 1  HE PABHO 2.',
     *'HAEHE YHK'/7X,'OOEHO PABHM 3.4E38')
  103 FORMAT(' N 2 - HAEHE APYMEHTA MEHE  PABHO 0.',
     *'HAEHE YHK'/7X,'OOEHO PABHM 3.4E38')
  104 FORMAT(' N 3 - HAEHE APYMEHTA OE 43,39.',
     *' HAEHE YHK'/7X,'OOEHO PABHM 0.')
  112 FORMAT(' N 3 - HAEHE APYMEHTA OE 44.36.',
     *'HAEHE YHK'/7X,'OOEHO PABHM 0.')
  147 FORMAT(' N 1 - HAEHE L MEHE 1  OE 3.',
     *'HAEHE YHK'/7X,'OOEHO PABHM 3.4E38')
  148 FORMAT(' N 2 - HAEHE APYMEHTA MEHE 0 P L=1.',
     *'HAEHE YHK'/7X,'OOEHO PABHM 3.4E38')
  149 FORMAT(' N 3 - HAEHE APYMEHTA OE  PABHO 1',
     *' P L=1  L=2.'/7X,
     *'HAEHE YHK OOEHO PABHM 3.4E38')
 1490 FORMAT(' N 4 - HAEHE APYMEHTA OE 1 P L=3.',
     *' HAEHE YHK'/7X,'OOEHO PABHM 3.4E38')
 1491 FORMAT(' N 5 - HAEHE APYMEHTA MEHE  PABHO 0 P L=3.'/
     *7X,'HAEHE YHK OOEHO PABHM 3.4E38')
  153 FORMAT(' N 3 - HAEHE APYMEHTA OE 1 P L=1,2,3.',
     *'HAEHE'/7X,'YHK OOEHO PABHM 3.4E38')
C  163 FORMAT(' OTEKA HB MY,OPOPAMMA SF16R:',
C     *'HEATAHA OKA')
  164 FORMAT (' N 1 - HAEHE APYMEHTA MEHE 0.'/
     *'       HAEH AKER  AKEI OOEH PABHM 3.4E38')
  165 FORMAT(' OTEKA HB MY,OPOPAMMA SF16R:',
     *' ATAHA OKA')
  166 FORMAT(' N 2 - ACOTHOE HAEHE APYMEHTA OE 57.58.',
     *' HAEH BER  BEI'/7X,'OOEH PABHM 0.',
     *'B CYAE KOA APYMEHT HEOTPATEEH,'/
     *'       AKER  AKEI TAKE OOEH PABHM 0,',
     *'HAE AKER  AKEI'/7X,'OAATC PABHM 3.4E38')
C  173 FORMAT(' OTEKA HB MY,OPOPAMMA SF17R:',
C     *'HEATAHA OKA')
  174 FORMAT(' N 1 - HAEHE APYMEHTA MEHE 0.'/
     *'       HAEH AKERD  AKEID OOEH PABHM 3.4E38')
  175 FORMAT(' N 2 - ACOTHOE HAEHE APYMEHTA OE 57.58.'/
     *'       HAEH BERD  BEID OOEH PABHM 0,',
     *'B CYAE,'/7X,'KOA APYMEHT HEOTPATEEH,',
     *' AKERD  AKEID TAKE'/7X,'OOEH PABHM 0,',
     *'HAE AKERD  AKEID OAATC,'/7X,'PABHM 3.4E38')
C  183 FORMAT(' OTEKA HB MY, OPOPAMMA SF18R:',
C     *'HEATAHA OKA')
  184 FORMAT(' N 1 - HAEHE APYMEHTA MEHE 0.'/
     *'       HAEH AKER1  AKEI1 OOEH PABHM 3.4E38')
  185 FORMAT(' OTEKA HB MY,OPOPAMMA SF18R:',
     *'ATAHA OKA')
  186 FORMAT(' N 2 - ACOTHOE HAEHE APYMEHTA OE 57.58.'/
     *'       HAEH BER1  BEI1 OOEH PABHM 0.',
     *'B CYAE, KOA'/7X,'APYMEHT HEOTPATEEH,',
     *' AKER1  AKEI1 TAKE OOEH PABHM 0,'/7X,
     *'HAE AKER1  AKEI1 OAATC PABHM 3.4E38')
  191 FORMAT(' OTEKA HB MY, YHK SF19R: ATAHA OKA')
  201 FORMAT(' OTEKA HB MY, YHK SF20R: ATAHA OKA')
  210 FORMAT(' OTEKA HB MY,YHK SF21R: ATAHA OKA')
  220 FORMAT(' OTEKA HB MY, YHK SF22R:ATAHA OKA')
  192 FORMAT(' N 1 - HEBEPHO AAH AKTECK APAMETP RK2: 0<RK2<1'/
     *7X,'HAEHE YHK OOEHO PABHM 3.4E38')
  193 FORMAT(' N 2 - HEBEPHO AAH AKTECK APAMETP T: 0<T<1.5707',
     *'9632679'/
     1'       HAEHE YHK OOEHO PABHM 3.4E38')
  194 FORMAT(' N 3 - P RK2=1  T=1.57079632679 HAMEHATE ',
     *'OHTEPAHO YHK'/7X,'PABEH HY',
     2'HAEHE YHK OOEHO PABHM 3.4E38')
  213 FORMAT(' N 1 - AMMA-YHK HE OPEEEHA P X PABHOM HY ',
     *' OMY'/7X,'EOMY OTPATEHOMY CY.',
     1' HAEHE SF21R OOEHO PABHM 3.4E38')
  214 FORMAT(' N 2 - HAEHE APYMEHTA OE 21.43. HAEHE',
     1' YHK'/7X,'OOEHO PABHM 3.4E38')
  223 FORMAT(' N 1 - TA-YHK HE OPEEEHA P X, Y,  X+Y ',
     1'PABHX HY '/7X,
     *' OMY EOMY OTPATEHOMY CY.',
     2' HAEHE YHK'/7X,'OOEHO PABHM 3.4E38')
  222 FORMAT(' N 2 - HAEHE APYMEHTA MEHE -20.4',
     *' HAEHE YHK'/7X,'OOEHO PABHM 3.4E38')
  230 FORMAT(' OTEKA HB MY,O',
     *'POPAMMA SF23C: ATAHA OKA')
  233 FORMAT(' N 1 - BEECTBEHHA ACT ',
     *'APYMEHTA HE MOET T OE 21.4' )
  234 FORMAT(' N 2 - AMMA-YHK ',
     *'HE OPEEEHA P X=0  OMY ',
     *'EOMY'/7X,'OTPATEHOMY CY. HAEHE ',
     1' YHK OOEHO PABHM 3.4E38')
  281 FORMAT (' OTEKA HB MY, YHK SF28R: HEATAHA OKA'/
     *' N 1 - HAEHE APYMEHTA MEHE 0;',
     *'  YHK BCETC '/7X,'ACOTHOO HAEH APYMEHTA')
  283 FORMAT (' OTEKA HB MY, YHK SF28R: ATAHA OKA'/
     *' N 2 - HAEHE APYMEHTA OE 46.499'/7X,
     *'HAEHE YHK OOEHO PABHM 3.4E38')
  291 FORMAT (' OTEKA HB MY,YHK SF29R: HEATAHA OKA'/
     *' N 1 - HAEHE APYMEHTA MEHE 0;',
     *'YHK BCETC '/7X,'ACOTHOO HAEH APYMEHTA')
  293 FORMAT (' OTEKA HB MY,YHK SF29R: ATAHA OKA'/
     *' N 2 - HAEHE APYMEHTA OE 46.509',
     *' HAEHE YHK OOEHO'/7X,'PABHM 3.4E38')
  311 FORMAT(' OTEKA HB MY, YHK SF31R: ATAHA OKA'/
     *' N 1 - HTEPAH KOCHYC HE OPEEEH  HEOTPATEHX',
     */7X,'APYMEHTOB HAEHE YHK OOEHO PABHM 3.4E38')
  321 FORMAT (' OTEKA HB MY,OPOPAMMA SF32R:ATAHA OKA'
     */)
  324 FORMAT (' N 1 - HAEHE APYMEHTA OTPATEHO.')
  325 FORMAT (' N 2 - AAHH OPOK CEPECKO YHK',
     *' ECCE'/7X,'HE PHAET HTEPBAY (0,29).')
  331 FORMAT(' OTEKA HB MY, ',
     *'OPOPAMMA SF33R: ATAHA OKA')
  336 FORMAT(' N 2 - NMAX (PAMEPHOCT',
     *' MACCBOB FUN1  FUN2) > MAX',
     *'(N,2*[X]+2)'/)
  337 FORMAT(' N 3 - HAEHE APYME',
     *'HTA X OHO T OE 0')
  338 FORMAT(' N 1 - N (CO HAEH ',
     *'OPKA,  KOTOPX HYHO'/7X,'OC',
     *'TAT YHK) > 2')
  339 FORMAT(' N 4 - HAEHE ',
     *'BCEMO YHK PEBOCX',
     *'OT HAOEE CO,'/7X,'PE',
     *'CTABMOE B MAHE (CM. COH',
     *'TEH APEC)'/)
  340 PRINT 341
  341 FORMAT('   , SFG4R:',
     *' ')
      GO TO 342
  342 GO TO (343,344,345),I
  343 PRINT 346
  346 FORMAT(' N 1 -   X  0,  ',
     *' '/7X,'    ',
     *' ')
      RETURN
  344 PRINT 347
  347 FORMAT(' N 2 -       0, ',
     *''/7X,'   ',
     *' '/7X,'  ')
      RETURN
  345 PRINT 348
  348 FORMAT(' N 3 -    ',
     *'  '/7X,'  A  EPS, ',
     *'  '/7X,'  ',
     *'   ')
      RETURN
  350 PRINT 351
  351 FORMAT('   , SFG5R:',
     *' ')
      GO TO 352
  352 GO TO (353,354,355),I
  353 PRINT 356
  356 FORMAT(' N 1 -   X  0,  ',
     *' '/7X,'    ',
     *' ')
      RETURN
  354 PRINT 357
  357 FORMAT(' N 2 -   A    0, ',
     *''/7X,'   ',
     *' '/7X,'  ')
      RETURN
  355 PRINT 358
  358 FORMAT(' N 3 -    ',
     *'  '/7X,'  A  EPS, ',
     *'  '/7X,'  ',
     *'   ')
      RETURN
      END
