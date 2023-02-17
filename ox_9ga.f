C

	SUBROUTINE OX(DTIMEB, IXNTIME)

c
C   New Ox routine that uses the AER chemistry (diurnally avgs) 
C      with the Stolarski odd oxygen counting scheme                 - EF  2/24/03
C      also uses long lived constituents that have 
C      just been transported
C
C  NOTE: THe C array for Ox here has been transported, so update 
C        the C array with chemistry below to become the CN array
C
C
C   FOR J and K reactions which are products of 2 or more diurnally varying things, we need to 
C      sum up over the 24-hour cycle first, and then average - use JDC(PH$,18,L$,Z$) and CNDC(S$,18+3,L$,Z$)
C
C      PURE Ox DAY/NIGHT EXTENSION FOR 95 KM MODEL (7/14/88)
C
C
C   this is Ox_9be.f = Ox prod/loss computed diagnostically


        include "com2d.h"


        INTEGER IXNTIME

        COMMON/OXCC/OXTR(L$,Z$)
        REAL OXVAL(L$,Z$), OXCOL(L$), OXTEN(13,L$)

        
        REAL POXJO2(L$,Z$), POXOTH(L$,Z$), 
     >       BALOX(L$,Z$), BALHOX(L$,Z$), BALNOX(L$,Z$), BALCLX(L$,Z$),
     >       BALBRX(L$,Z$), BALCHX(L$,Z$), BALCFC(L$,Z$), 
     >       BALRAI(L$,Z$), BALTOT(L$,Z$)


        DOUBLE PRECISION DTIMEB(L$,IXNTIME)
        REAL*8 jsum10, jsum11, jsum13, jsum17, jsum36, ksum(175)
        REAL*8 ksumh2, ksumh3, ksumh5, ksumh6, ksumh7


        SAVE


C
C  **********************************************************************************************
c
c  start new Stolarski counting scheme here:
C

       DO 1000 ik=1,Z$
       DO 1000 ij=1,L$

C
C
C  Ox PRODUCTION:
C  --------------
C
C    first compute diurnally varying terms
C
C
      jsum10 = 0.D0   
      jsum11 = 0.D0   
      jsum13 = 0.D0   
      jsum17 = 0.D0   
      jsum36 = 0.D0   

      do 212 iip=1,175
 212	 ksum(iip) = 0.D0

      ksumh2 = 0.D0
      ksumh3 = 0.D0
      ksumh5 = 0.D0
      ksumh6 = 0.D0
      ksumh7 = 0.D0

      DO 200 it=1,ixntime-1
        jsum10 = jsum10 + 0.5D0*(JDC(10,it,ij,ik)*CNDC(23,it,ij,ik) +
     >                           JDC(10,it+1,ij,ik)*CNDC(23,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

        jsum11 = jsum11 + 0.5D0*(JDC(11,it,ij,ik)*CNDC(23,it,ij,ik) +
     >                           JDC(11,it+1,ij,ik)*CNDC(23,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

        jsum13 = jsum13 + 0.5D0*(JDC(13,it,ij,ik)*CNDC(24,it,ij,ik) +
     >                           JDC(13,it+1,ij,ik)*CNDC(24,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

        jsum17 = jsum17 + 0.5D0*(JDC(17,it,ij,ik)*CNDC(7,it,ij,ik) +
     >                           JDC(17,it+1,ij,ik)*CNDC(7,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

        jsum36 = jsum36 + 0.5D0*(JDC(36,it,ij,ik)*CNDC(61,it,ij,ik) +
     >                           JDC(36,it+1,ij,ik)*CNDC(61,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0


       ksum(2) = ksum(2) + 0.5D0*(CNDC(1,it,ij,ik)*CNDC(4,it,ij,ik) 
     >                     +    CNDC(1,it+1,ij,ik)*CNDC(4,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(5) = ksum(5) + 0.5D0*(CNDC(14,it,ij,ik)*CNDC(4,it,ij,ik) 
     >                   +       CNDC(14,it+1,ij,ik)*CNDC(4,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(6) = ksum(6) + 0.5D0*(CNDC(28,it,ij,ik)*CNDC(14,it,ij,ik) 
     >                   +    CNDC(28,it+1,ij,ik)*CNDC(14,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(15) = ksum(15) + 0.5D0*(CNDC(5,it,ij,ik)*CNDC(22,it,ij,ik) 
     >                     +    CNDC(5,it+1,ij,ik)*CNDC(22,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(21) = ksum(21) + 0.5D0*(CNDC(1,it,ij,ik)*CNDC(23,it,ij,ik) 
     >                     +    CNDC(1,it+1,ij,ik)*CNDC(23,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(22) = ksum(22) + 0.5D0*(CNDC(14,it,ij,ik)*CNDC(22,it,ij,ik) 
     >                     +    CNDC(14,it+1,ij,ik)*CNDC(22,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(25) = ksum(25) + 0.5D0*(CNDC(1,it,ij,ik)*CNDC(28,it,ij,ik) 
     >                     +    CNDC(1,it+1,ij,ik)*CNDC(28,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(27) = ksum(27) + 0.5D0*(CNDC(13,it,ij,ik)*CNDC(29,it,ij,ik) 
     >                     +    CNDC(13,it+1,ij,ik)*CNDC(29,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(29) = ksum(29) + 0.5D0*(CNDC(13,it,ij,ik)*CNDC(16,it,ij,ik) 
     >                     +    CNDC(13,it+1,ij,ik)*CNDC(16,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(40) = ksum(40) + 0.5D0*(CNDC(13,it,ij,ik)*CNDC(14,it,ij,ik)  
     >                     +    CNDC(13,it+1,ij,ik)*CNDC(14,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(41) = ksum(41) + 0.5D0*(CNDC(13,it,ij,ik)*CNDC(1,it,ij,ik) 
     >                     +     CNDC(13,it+1,ij,ik)*CNDC(1,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(42) = ksum(42) + 0.5D0*(CNDC(14,it,ij,ik)*CNDC(1,it,ij,ik)  
     >                     +     CNDC(14,it+1,ij,ik)*CNDC(1,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(43) = ksum(43) + 0.5D0*(CNDC(6,it,ij,ik)*CNDC(1,it,ij,ik) 
     >                     +     CNDC(6,it+1,ij,ik)*CNDC(1,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(51) = ksum(51) + 0.5D0*(CNDC(13,it,ij,ik)*CNDC(23,it,ij,ik) 
     >                     +    CNDC(13,it+1,ij,ik)*CNDC(23,it+1,ij,ik))
     >                                       *DTIMEB(ij,it)/86400.D0

       ksum(53) = ksum(53) + 0.5D0*(CNDC(14,it,ij,ik)*CNDC(27,it,ij,ik) 
     >                     +    CNDC(14,it+1,ij,ik)*CNDC(27,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(55) = ksum(55) + 0.5D0*(CNDC(13,it,ij,ik)*CNDC(38,it,ij,ik) 
     >                     +    CNDC(13,it+1,ij,ik)*CNDC(38,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(58) = ksum(58) + 0.5D0*(CNDC(13,it,ij,ik)*CNDC(24,it,ij,ik) 
     >                     +    CNDC(13,it+1,ij,ik)*CNDC(24,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(59) = ksum(59) + 0.5D0*(CNDC(13,it,ij,ik)*CNDC(13,it,ij,ik) 
     >                     +    CNDC(13,it+1,ij,ik)*CNDC(13,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(61) = ksum(61) + 0.5D0*(CNDC(13,it,ij,ik)*CNDC(28,it,ij,ik) 
     >                     +    CNDC(13,it+1,ij,ik)*CNDC(28,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(62) = ksum(62) + 0.5D0*(CNDC(13,it,ij,ik)*CNDC(25,it,ij,ik) 
     >                     +    CNDC(13,it+1,ij,ik)*CNDC(25,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(63) = ksum(63) + 0.5D0*(CNDC(27,it,ij,ik)*CNDC(23,it,ij,ik) 
     >                     +    CNDC(27,it+1,ij,ik)*CNDC(23,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(64) = ksum(64) + 0.5D0*(CNDC(14,it,ij,ik)*CNDC(14,it,ij,ik) 
     >                     +    CNDC(14,it+1,ij,ik)*CNDC(14,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(71) = ksum(71) + 0.5D0*(CNDC(12,it,ij,ik)*CNDC(14,it,ij,ik) 
     >                     +    CNDC(12,it+1,ij,ik)*CNDC(14,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(81) = ksum(81) + 0.5D0*(CNDC(9,it,ij,ik)*CNDC(6,it,ij,ik) 
     >                     +     CNDC(9,it+1,ij,ik)*CNDC(6,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(90) = ksum(90) + 0.5D0*(CNDC(1,it,ij,ik)*CNDC(1,it,ij,ik) 
     >                     +     CNDC(1,it+1,ij,ik)*CNDC(1,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(92) = ksum(92) + 0.5D0*(CNDC(14,it,ij,ik)*CNDC(45,it,ij,ik) 
     >                     +    CNDC(14,it+1,ij,ik)*CNDC(45,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(93) = ksum(93) + 0.5D0*(CNDC(44,it,ij,ik)*CNDC(28,it,ij,ik) 
     >                     +    CNDC(44,it+1,ij,ik)*CNDC(28,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(94) = ksum(94) + 0.5D0*(CNDC(44,it,ij,ik)*CNDC(44,it,ij,ik) 
     >                     +    CNDC(44,it+1,ij,ik)*CNDC(44,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(95) = ksum(95) + 0.5D0*(CNDC(13,it,ij,ik)*CNDC(46,it,ij,ik) 
     >                     +    CNDC(13,it+1,ij,ik)*CNDC(46,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(107) = ksum(107) + 0.5D0*(CNDC(44,it,ij,ik)*CNDC(1,it,ij,ik) 
     >                       +   CNDC(44,it+1,ij,ik)*CNDC(1,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

      ksum(108) = ksum(108) + 0.5D0*(CNDC(44,it,ij,ik)*CNDC(14,it,ij,ik) 
     >                       +  CNDC(44,it+1,ij,ik)*CNDC(14,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

      ksum(109) = ksum(109) + 0.5D0*(CNDC(45,it,ij,ik)*CNDC(23,it,ij,ik) 
     >                       +  CNDC(45,it+1,ij,ik)*CNDC(23,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

      ksum(140) = ksum(140) + 0.5D0*(CNDC(28,it,ij,ik)*CNDC(41,it,ij,ik) 
     >                       +  CNDC(28,it+1,ij,ik)*CNDC(41,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(141) = ksum(141) + 0.5D0*(CNDC(1,it,ij,ik)*CNDC(38,it,ij,ik) 
     >                       +   CNDC(1,it+1,ij,ik)*CNDC(38,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

      ksum(150) = ksum(150) + 0.5D0*(CNDC(13,it,ij,ik)*CNDC(42,it,ij,ik) 
     >                       +  CNDC(13,it+1,ij,ik)*CNDC(42,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

      ksum(151) = ksum(151) + 0.5D0*(CNDC(2,it,ij,ik)*CNDC(4,it,ij,ik) 
     >                      +     CNDC(2,it+1,ij,ik)*CNDC(4,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(153) = ksum(153) + 0.5D0*(CNDC(1,it,ij,ik)*CNDC(7,it,ij,ik) 
     >                       +    CNDC(1,it+1,ij,ik)*CNDC(7,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(158) = ksum(158) + 0.5D0*(CNDC(28,it,ij,ik)*CNDC(7,it,ij,ik) 
     >                       +   CNDC(28,it+1,ij,ik)*CNDC(7,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

      ksum(160) = ksum(160) + 0.5D0*(CNDC(28,it,ij,ik)*CNDC(28,it,ij,ik) 
     >                       +  CNDC(28,it+1,ij,ik)*CNDC(28,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(162) = ksum(162) + 0.5D0*(CNDC(7,it,ij,ik)*CNDC(7,it,ij,ik) 
     >                        +   CNDC(7,it+1,ij,ik)*CNDC(7,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

       ksum(163) = ksum(163) + 0.5D0*(CNDC(62,it,ij,ik)*CNDC(1,it,ij,ik) 
     >                       +   CNDC(62,it+1,ij,ik)*CNDC(1,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

      ksum(164) = ksum(164) + 0.5D0*(CNDC(62,it,ij,ik)*CNDC(13,it,ij,ik) 
     >                       +  CNDC(62,it+1,ij,ik)*CNDC(13,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

      ksum(175) = ksum(175) + 0.5D0*(CNDC(27,it,ij,ik)*CNDC(61,it,ij,ik) 
     >                      +   CNDC(27,it+1,ij,ik)*CNDC(61,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0


        ksumh2 = ksumh2 + 0.5D0*(KH(2,it,ij,ik)*CNDC(30,it,ij,ik) + 
     >                         KH(2,it+1,ij,ik)*CNDC(30,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksumh3 = ksumh3 + 0.5D0*(KH(3,it,ij,ik)*CNDC(8,it,ij,ik) + 
     >                         KH(3,it+1,ij,ik)*CNDC(8,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksumh6 = ksumh6 + 0.5D0*(KH(6,it,ij,ik)*CNDC(47,it,ij,ik) + 
     >                         KH(6,it+1,ij,ik)*CNDC(47,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksumh5 = ksumh5 + 0.5D0*(KH(5,it,ij,ik)*CNDC(25,it,ij,ik) +
     >                         KH(5,it+1,ij,ik)*CNDC(25,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksumh7 = ksumh7 + 0.5D0*(KH(7,it,ij,ik)*CNDC(68,it,ij,ik) +
     >                        KH(7,it+1,ij,ik)*CNDC(68,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0
 200  CONTINUE

                               ! load in duplicates
       ksum(104) = ksum(94)
       ksum(105) = ksum(93)
       ksum(148) = ksum(58)

C
C    POXJO2 and POXOTH are in num den/sec
C
C
	POXJO2(IJ,IK) = (J(1,IJ,IK) + J(46,IJ,IK))*C(3,IJ,IK)*2.


	POXOTH(IJ,IK) = (J(4,IJ,IK) + J(25,IJ,IK))*c(15,IJ,IK)
     >                + (J(12,IJ,IK) + J(41,IJ,IK))*c(20,IJ,IK)
     >         + J(14,IJ,IK)*c(11,IJ,IK) + jsum13*2.D0                  ! J(13,IJ,IK)*c(24,IJ,IK)*2.
     >         + J(29,IJ,IK)*c(49,IJ,IK)*2.
     >         + K(3,IJ,IK)*C(12,IJ,IK)*C(3,IJ,IK)*2. 
     >         + K(14,IJ,IK)*C(18,IJ,IK)*C(13,IJ,IK)
     >         + K(15,IJ,IK)*ksum(15)*2.D0                              ! C(22,IJ,IK)*C(5,IJ,IK)
     >         + K(20,IJ,IK)*C(9,IJ,IK)*C(3,IJ,IK)*2.
     >         + K(26,IJ,IK)*C(18,IJ,IK)*C(27,IJ,IK)*2.
     >         + K(39,IJ,IK)*C(15,IJ,IK)*C(2,IJ,IK)
     >         + K(45,IJ,IK)*C(11,IJ,IK)*C(2,IJ,IK)
     >         + K(49,IJ,IK)*C(18,IJ,IK)*C(2,IJ,IK)*2.
     >         + K(124,IJ,IK)*C(49,IJ,IK)*C(2,IJ,IK)*2.
     >         + K(140,IJ,IK)*ksum(140)*M(IJ,IK)*2.D0                  ! C(28,IJ,IK)*C(41,IJ,IK)*M(IJ,IK)*2.
                                                                       ! per AER code, M NOT included in K140
     >         + ksumh2                                                ! HET, ClONO2+H2O
     >         + ksumh3                                                ! HET, N2O5+H2O
     >         + ksumh6                                                ! HET, BrONO2+H2O

cccccccccc     >         + K(79,IJ,IK)*C(8,IJ,IK)*C(15,IJ,IK)  -  gas phase N2O5+H2O = 0.0


C                                                              PROD is in num den./sec
cox	PROD = POXJO2(IJ,IK) + POXOTH(IJ,IK)


C  Load in Ox prod into arrays for write out in MAIN;  CHMPOX(2,L$,Z$), 1=PJO2;  2=POTH;  IN MIX RAT/SEC

        chmpox(1,ij,ik) = POXJO2(IJ,IK)/M(IJ,IK)
        chmpox(2,ij,ik) = POXOTH(IJ,IK)/M(IJ,IK)

        CTPROD(1,IJ,IK) = (POXJO2(IJ,IK) + POXOTH(IJ,IK))/M(IJ,IK)


CC *********************************************************************************************
C
C
C
C  Ox LOSS:
C  --------
C             BALOX, BALHOX, BALNOX, BALCLX, BALBRX, BALCHX, BALCFC, BALRAI are in num den/second
C           
C

        BALOX(IJ,IK) = K(2,IJ,IK)*ksum(2)*2.D0                          ! C(1,IJ,IK)*C(4,IJ,IK)*2.
     >               + K(90,IJ,IK)*ksum(90)*M(IJ,IK)*2.D0               ! C(1,IJ,IK)*C(1,IJ,IK)*2.   
     >               + K(151,IJ,IK)*ksum(151)*2.D0                      ! C(2,IJ,IK)*C(4,IJ,IK)*2.



        BALHOX(IJ,IK) = K(5,IJ,IK)*ksum(5)*2.D0                         ! C(14,IJ,IK)*C(4,IJ,IK)*2.
     >                + K(29,IJ,IK)*ksum(29)                            ! C(13,IJ,IK)*C(16,IJ,IK)
     >                + K(30,IJ,IK)*C(17,IJ,IK)*C(13,IJ,IK)
     >                + K(40,IJ,IK)*ksum(40)*3.D0                       ! C(13,IJ,IK)*C(14,IJ,IK)*3.
     >                + K(41,IJ,IK)*ksum(41)*2.D0                       ! C(13,IJ,IK)*C(1,IJ,IK)*2.
     >                + K(42,IJ,IK)*ksum(42)*2.D0                       ! C(14,IJ,IK)*C(1,IJ,IK)*2.
     >                + K(59,IJ,IK)*ksum(59)                            ! C(13,IJ,IK)*C(13,IJ,IK)
     >                + K(64,IJ,IK)*ksum(64)*2.D0               ! K64 here is k64 + k65*M,  C(14,IJ,IK)*C(14,IJ,IK)*2.
     >                + K(71,IJ,IK)*ksum(71)*2.D0                       ! C(12,IJ,IK)*C(14,IJ,IK)*2.
     >                + K(73,IJ,IK)*ksum(71)                            ! C(12,IJ,IK)*C(14,IJ,IK) -> H2O+O


	BALNOX(IJ,IK) = jsum17*2.D0                                    ! J(17,IJ,IK)*C(7,IJ,IK)*2.
     >                + K(37,IJ,IK)*C(10,IJ,IK)*C(13,IJ,IK)
     >                + K(43,IJ,IK)*ksum(43)*2.D0                      ! C(6,IJ,IK)*C(1,IJ,IK)*2.
     >                + K(55,IJ,IK)*ksum(55)*3.D0                      ! C(13,IJ,IK)*C(38,IJ,IK)*3.
     >                + K(81,IJ,IK)*ksum(81)                           ! C(9,IJ,IK)*C(6,IJ,IK)
     >                + K(141,IJ,IK)*ksum(141)*2.D0                    ! C(38,IJ,IK)*C(1,IJ,IK)*2.  
     >                + K(150,IJ,IK)*ksum(150)                         ! C(42,IJ,IK)*C(13,IJ,IK)
     >                + K(153,IJ,IK)*ksum(153)*2.D0                    ! C(7,IJ,IK)*C(1,IJ,IK)*2.
     >                + K(162,IJ,IK)*ksum(162)*2.D0                    ! C(7,IJ,IK)*C(7,IJ,IK)*2.


	BALCLX(IJ,IK) = jsum36*2.D0                                     ! J(36,IJ,IK)*C(61,IJ,IK)*2.
     >                + K(6,IJ,IK)*ksum(6)*2.D0                         ! C(28,IJ,IK)*C(14,IJ,IK)*2.
     >                + K(25,IJ,IK)*ksum(25)*2.D0                       ! C(28,IJ,IK)*C(1,IJ,IK)*2.
     >                + K(27,IJ,IK)*ksum(27)                            ! C(29,IJ,IK)*C(13,IJ,IK)
     >                + K(53,IJ,IK)*ksum(53)*2.D0                       ! C(27,IJ,IK)*C(14,IJ,IK)*2.
     >                + K(61,IJ,IK)*ksum(61)*2.D0                       ! C(28,IJ,IK)*C(13,IJ,IK)*2.
     >                + K(62,IJ,IK)*ksum(62)                            ! C(25,IJ,IK)*C(13,IJ,IK)
     >                + K(158,IJ,IK)*ksum(158)*2.D0                     ! C(7,IJ,IK)*C(28,IJ,IK)*2.
     >                + K(160,IJ,IK)*ksum(160)*2.D0                     ! C(28,IJ,IK)*C(28,IJ,IK)*2.
     >                + K(161,IJ,IK)*ksum(160)*2.D0                     ! C(28,IJ,IK)*C(28,IJ,IK)*2.
     >                + K(163,IJ,IK)*ksum(163)*2.D0                     ! C(1,IJ,IK)*C(62,IJ,IK)*2.
     >                + K(164,IJ,IK)*ksum(164)*2.D0                     ! C(13,IJ,IK)*C(62,IJ,IK)*2.
     >                + K(175,IJ,IK)*ksum(175)*2.D0                     ! C(27,IJ,IK)*C(61,IJ,IK)*2.
     >                + ksumh5                                          ! HET, HOCl+HCl

ccccc     >                + K(139,IJ,IK)*C(28,IJ,IK)*C(14,IJ,IK)*2.     Ks = 0. here
ccccc     >                + K(143,IJ,IK)*C(30,IJ,IK)*C(13,IJ,IK)*2.     Ks = 0. here
ccccc     >                + K(144,IJ,IK)*C(28,IJ,IK)*C(10,IJ,IK)*2.     Ks = 0. here
ccccc     >                + K(173,IJ,IK)*C(25,IJ,IK)*C(46,IJ,IK)        Ks = 0. here



        BALBRX(IJ,IK) = K(92,IJ,IK)*ksum(92)*2.D0                       ! C(45,IJ,IK)*C(14,IJ,IK)*2.
     >                + K(93,IJ,IK)*ksum(93)*2.D0                       ! C(44,IJ,IK)*C(28,IJ,IK)*2.
     >                + K(94,IJ,IK)*ksum(94)*2.D0                       ! C(44,IJ,IK)*C(44,IJ,IK)*2.
     >                + K(95,IJ,IK)*ksum(95)                            ! C(13,IJ,IK)*C(46,IJ,IK)
     >                + K(104,ij,ik)*ksum(104)*2.D0                     ! c(44,IJ,IK)*c(44,IJ,IK)*2.
     >                + K(105,IJ,IK)*ksum(105)*2.D0                     ! C(44,IJ,IK)*C(28,IJ,IK)*2.
     >                + K(107,IJ,IK)*ksum(107)*2.D0                     ! C(44,IJ,IK)*C(1,IJ,IK)*2.
     >                + K(108,IJ,IK)*ksum(108)*2.D0                     ! C(44,IJ,IK)*C(14,IJ,IK)*2.
     >                + ksumh7                                          ! HET, HOBr+HCl

ccccc     >                + K(137,IJ,IK)*C(44,IJ,IK)*C(14,IJ,IK)*2.     Ks = 0. here
ccccc     >                + K(145,IJ,IK)*C(44,IJ,IK)*C(4,IJ,IK)*2.
ccccc     >                + K(170,IJ,IK)*C(44,IJ,IK)*C(13,IJ,IK)*2.
ccccc     >                + K(172,IJ,IK)*C(47,IJ,IK)*C(1,IJ,IK)*2.


C                                                                       ! BALCHx is the Ox loss from CHx

       BALCHX(IJ,IK) = jsum10 + jsum11                               ! J(10,IJ,IK) + J(11,IJ,IK))*c(23,IJ,IK)
     >               + K(21,IJ,IK)*ksum(21)                                 ! C(23,IJ,IK)*C(1,IJ,IK)
     >               + K(22,IJ,IK)*ksum(22)    ! *2.D0                      ! C(22,IJ,IK)*C(14,IJ,IK)  ! *2.
     >               + K(51,IJ,IK)*ksum(51)*2.D0                            ! C(23,IJ,IK)*C(13,IJ,IK)*2.
     >               + K(58,IJ,IK)*ksum(58)                                 ! C(13,IJ,IK)*C(24,IJ,IK)
     >               + K(63,IJ,IK)*ksum(63)                                 ! C(23,IJ,IK)*C(27,IJ,IK)
     >               + K(109,IJ,IK)*ksum(109)                               ! C(45,IJ,IK)*C(23,IJ,IK)
     >               + K(148,IJ,IK)*ksum(148)                               ! C(24,IJ,IK)*C(13,IJ,IK)




C                                                                        ! BALCFC is Ox loss from CFC breakup
       BALCFC(IJ,IK) =  K(16,IJ,IK)*C(13,IJ,IK)*C(37,IJ,IK)
     >                + K(54,IJ,IK)*C(36,IJ,IK)*C(2,IJ,IK)
     >                + K(66,IJ,IK)*c(2,IJ,IK)*c(57,IJ,IK)   
     >                + K(75,IJ,IK)*C(13,IJ,IK)*C(40,IJ,IK)
     >                + K(80,IJ,IK)*C(35,IJ,IK)*C(2,IJ,IK)
     >                + K(83,IJ,IK)*C(34,IJ,IK)*C(2,IJ,IK)
     >                + K(97,IJ,IK)*C(49,IJ,IK)*C(13,IJ,IK)
     >                + K(98,IJ,IK)*C(13,IJ,IK)*C(52,IJ,IK)
     >                + K(99,IJ,IK)*c(53,IJ,IK)*C(2,IJ,IK)
     >                + K(100,IJ,IK)*c(54,IJ,IK)*C(2,IJ,IK)
     >                + K(101,IJ,IK)*c(55,IJ,IK)*C(2,IJ,IK)
     >                + K(112,IJ,IK)*C(51,IJ,IK)*C(2,IJ,IK)
     >                + K(113,IJ,IK)*C(50,IJ,IK)*C(2,IJ,IK)
     >                + K(114,IJ,IK)*c(69,ij,ik)*c(2,ij,ik)
     >                + K(115,IJ,IK)*c(69,ij,ik)*c(13,ij,ik)
     >                + K(117,IJ,IK)*c(70,ij,ik)*c(2,ij,ik)
     >                + K(118,IJ,IK)*c(70,ij,ik)*c(13,ij,ik)
     >                + K(120,IJ,IK)*c(71,ij,ik)*c(2,ij,ik)
     >                + K(121,IJ,IK)*c(71,ij,ik)*c(13,ij,ik)
     >                + K(123,IJ,IK)*c(72,ij,ik)*c(2,ij,ik)
     >               + K(36,IJ,IK)*C(19,IJ,IK)*C(13,IJ,IK)
     >               + K(78,IJ,IK)*C(11,IJ,IK)*C(2,IJ,IK)
     >               + K(82,IJ,IK)*N2(IJ,IK)*C(2,IJ,IK)
     >               + K(111,IJ,IK)*C(1,IJ,IK)*C(19,IJ,IK)*M(IJ,IK) 



C                                          BALRAI is Ox loss from rainout of HNO3, H2O2, CH2O, CH3OOH

        BALRAI(IJ,IK) = (c(10,IJ,IK)*3. + c(16,IJ,IK)*2. 
     >                +  c(23,IJ,IK)    + c(24,IJ,IK)*2.)*c(59,ij,ik)




        BALTOT(IJ,IK) = BALOX(IJ,IK) + BALHOX(IJ,IK) + BALNOX(IJ,IK)
     >		      + BALCLX(IJ,IK) + BALBRX(IJ,IK) + BALCHX(IJ,IK)
     >                + BALCFC(IJ,IK) + BALRAI(IJ,IK)



C  BALOX, BALHOX, BALNOX, BALCLX, BALBRX, BALCHX, BALCFC, BALRAI, BALTOT are in num den/second
C
C                                                                        ! LOSS is in 1/second
cox        LOSS = BALTOT(IJ,IK)/C(79,IJ,IK)


c  update Ox at all levels - cn(79), here is just computed as a stand alone diagnostic

C  for surface dry deposition - OXDEP(L$,360) from HARVARD model are the daily values in m/sec,   
C  use DELTAZ(L$,Z$X) in KM (in COMMON) which is at the model grid points,   so DEPLOSS is in 1/sec

cox       IF (IK .eq. 1) THEN
cox
cox         deploss = OXDEP(ij,iday360)/(DELTAZ(ij,1)*1000.)
cox         cn(79,IJ,IK) = (c(79,IJ,IK) + prod*dt)/(1. + (loss+deploss)*dt)
C
cox       ELSE
cox
cox         cn(79,IJ,IK) = (c(79,IJ,IK) + prod*dt)/(1. + loss*dt)
cox
cox       ENDIF

cox    
c
C  Load in Ox loss into arrays for write out in MAIN;  
C   CHMLOX(8,L$,Z$): 1=LOx;  2=LHOx;  3=LNOx; 4=LClx; 5=LBrx; 6=LCHx; 7=LCFC;  8=Lrain,  ALL IN MIX RAT/SEC

       chmlox(1,ij,ik) = BALOX(IJ,IK)/M(IJ,IK)
       chmlox(2,ij,ik) = BALHOX(IJ,IK)/M(IJ,IK)
       chmlox(3,ij,ik) = BALNOX(IJ,IK)/M(IJ,IK)
       chmlox(4,ij,ik) = BALCLX(IJ,IK)/M(IJ,IK)
       chmlox(5,ij,ik) = BALBRX(IJ,IK)/M(IJ,IK)
       chmlox(6,ij,ik) = BALCHX(IJ,IK)/M(IJ,IK)
       chmlox(7,ij,ik) = BALCFC(IJ,IK)/M(IJ,IK)
       chmlox(8,ij,ik) = BALRAI(IJ,IK)/M(IJ,IK)
C                                                          ! CTLOSS(1,IJ,IK) is in MIX RATIO/SEC
       CTLOSS(1,IJ,IK) = BALTOT(IJ,IK)/M(IJ,IK)


1000   CONTINUE



CC *********************************************************************************************
C
C
C   NOW do loop to compute the column Ox tendencies due to PROD, LOSS, and transport:  
C     OXTEN(13,L$), OXVAL(L$,Z$), OXCOL(L$),    DELTAZ(L$,Z$X) is real*8
C
C   initialize:
C
        do 600 ij=1,L$
        do 600 iq=1,13
 600	   OXTEN(iq,ij) = 0.0



C  first load in updated Ox value for use in column integrations, OXVAL(L$,Z$)

       DO 1500 ik=1,Z$
       DO 1500 ij=1,L$
 1500	  oxval(ij,ik) = cn(39,ij,ik)


C  first integrate column OX - FIRST DEFINE TOTAL COLUMN ABOVE TOP MODEL LEVEL

        DO 602 IJ=1,L$
          OXCOL(ij) = OXVAL(IJ,Z$)*5.e5
c                                                                 integrate top-down
          DO 604 IK = Z$,1,-1
 604        OXCOL(ij) = OXCOL(ij) + OXVAL(IJ,IK)*DELTAZ(IJ,IK)*1.D5

            oxten(13,ij) = OXCOL(ij)
 602	CONTINUE


       iq = 1
       DO 701 IJ = 1,L$
       DO 701 IK = Z$,1,-1
         dq = DELTAZ(IJ,IK)*1.D5
       oxten(iq,ij)= oxten(iq,ij) + POXJO2(IJ,IK)*dq*OXVAL(IJ,IK)*dq
 701   CONTINUE


       iq = 2
       DO 702 IJ = 1,L$
       DO 702 IK = Z$,1,-1
         dq = DELTAZ(IJ,IK)*1.D5
       oxten(iq,ij)= oxten(iq,ij) + POXOTH(IJ,IK)*dq*OXVAL(IJ,IK)*dq
 702   CONTINUE


       iq = 3
       DO 703 IJ = 1,L$
       DO 703 IK = Z$,1,-1
         dq = DELTAZ(IJ,IK)*1.D5
       oxten(iq,ij) = oxten(iq,ij) + BALOX(IJ,IK)*dq*OXVAL(IJ,IK)*dq
 703   CONTINUE


       iq = 4
       DO 704 IJ = 1,L$
       DO 704 IK = Z$,1,-1
         dq = DELTAZ(IJ,IK)*1.D5
       oxten(iq,ij)= oxten(iq,ij) + BALHOX(IJ,IK)*dq*OXVAL(IJ,IK)*dq
 704   CONTINUE


       iq = 5
       DO 705 IJ = 1,L$
       DO 705 IK = Z$,1,-1
         dq = DELTAZ(IJ,IK)*1.D5
       oxten(iq,ij)= oxten(iq,ij) + BALNOX(IJ,IK)*dq*OXVAL(IJ,IK)*dq
 705   CONTINUE


       iq = 6
       DO 706 IJ = 1,L$
       DO 706 IK = Z$,1,-1
         dq = DELTAZ(IJ,IK)*1.D5
       oxten(iq,ij)= oxten(iq,ij) + BALCLX(IJ,IK)*dq*OXVAL(IJ,IK)*dq
 706   CONTINUE


       iq = 7
       DO 707 IJ = 1,L$
       DO 707 IK = Z$,1,-1
         dq = DELTAZ(IJ,IK)*1.D5
       oxten(iq,ij)= oxten(iq,ij) + BALBRX(IJ,IK)*dq*OXVAL(IJ,IK)*dq
 707   CONTINUE


       iq = 8
       DO 708 IJ = 1,L$
       DO 708 IK = Z$,1,-1
         dq = DELTAZ(IJ,IK)*1.D5
       oxten(iq,ij)= oxten(iq,ij) + BALCHX(IJ,IK)*dq*OXVAL(IJ,IK)*dq
 708   CONTINUE


       iq = 9
       DO 709 IJ = 1,L$
       DO 709 IK = Z$,1,-1
         dq = DELTAZ(IJ,IK)*1.D5
       oxten(iq,ij)= oxten(iq,ij) + BALCFC(IJ,IK)*dq*OXVAL(IJ,IK)*dq
 709   CONTINUE


       iq = 10
       DO 710 IJ = 1,L$
       DO 710 IK = Z$,1,-1
         dq = DELTAZ(IJ,IK)*1.D5
       oxten(iq,ij)= oxten(iq,ij) + BALRAI(IJ,IK)*dq*OXVAL(IJ,IK)*dq
 710   CONTINUE


C  for surface dry deposition - OXDEP(L$,360) from HARVARD model are the daily values in m/sec,   
C  use DELTAZ(L$,Z$X) in KM (in COMMON) which is at the model grid points to convert to 1/sec,   
C  then multiply by OX num den, so DEPLOSS is in num den/second

       iq = 11
       DO 711 IJ = 1,L$
         dq = DELTAZ(IJ,1)*1.D5
         deploss = OXDEP(ij,iday360)/(DELTAZ(ij,1)*1000.)*OXVAL(IJ,1)
         oxten(iq,ij) = oxten(iq,ij) + deploss*dq*OXVAL(IJ,1)*dq
 711   CONTINUE


C                    transport term, OXTR(L$,Z$) is in num den/sec
       iq = 12
       DO 712 IJ = 1,L$
       DO 712 IK = Z$,1,-1
         dq = DELTAZ(IJ,IK)*1.D5
      oxten(iq,ij) = oxten(iq,ij) + OXTR(IJ,IK)*dq*OXVAL(IJ,IK)*dq
 712   CONTINUE

                              ! normalize all PROD/LOSS and transport to total column density
        do 650 ij=1,L$
        do 650 iq=1,12
 650	   OXTEN(iq,ij) = OXTEN(iq,ij)/OXCOL(ij)



C    and write out every day:  OXTEN(13,L$)


       WRITE (213) OXTEN


       RETURN
       END

