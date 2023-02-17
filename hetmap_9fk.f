C
C 
	SUBROUTINE HETMAP(RAER,nkr,IJ,IK,IANGG,gfaer, gfnat, gfice)
C

        include "com2d.h"


        INTEGER NKR
        REAL*8 RAER(NKR)

        REAL gfaer, gfnat, gfice

        SAVE


C   LOAD in REACTION RATES from AER array R/RAER(207) into GSFC array , khaer(207,18,L$,Z$)
C        GFAEROUT(3,L$,Z$) for aerosols, PSCs 
C      - this is done at the end of each time step, and at the new NOONTIME time step (NTIME=18)
C                                                     also need to account for the H2O factor in KH(2,3,6)
C
C  NOTE: khaer(207,18,L$,Z$) appears to only be used to output diurnally varying reaction rates.....

ccoutput        do 100 ikh=1,nkr
ccoutput 100	   khaer(ikh,iangg,ij,ik) = RAER(ikh)


        GFAEROUT(1,ij,ik) = gfaer
        GFAEROUT(2,ij,ik) = gfnat
        GFAEROUT(3,ij,ik) = gfice


c  also load into array KH(RH$=14,18,L$,Z$) (REAL*8) in COMMON, for SOLVER/OX

      KH(1,iangg,ij,ik) = RAER(91) + RAER(161) + RAER(166)
      KH(2,iangg,ij,ik)= (RAER(93) + RAER(160) + RAER(165))  ! /cn(15,ij,ik)         ! need to divide by H2O
      KH(3,iangg,ij,ik)= (RAER(92) + RAER(158) + RAER(163))  ! /cn(15,ij,ik)         ! need to divide by H2O
      KH(4,iangg,ij,ik) =            RAER(159) + RAER(164)
      KH(5,iangg,ij,ik) = RAER(90) + RAER(162) + RAER(167)
      KH(6,iangg,ij,ik)= (RAER(122)            + RAER(169))  ! /cn(15,ij,ik)         ! need to divide by H2O
      KH(7,iangg,ij,ik) = RAER(124)            + RAER(168)

      KH(8,iangg,ij,ik) =                        RAER(285)
      KH(9,iangg,ij,ik) =                        RAER(275)
      KH(10,iangg,ij,ik)=            RAER(277) + RAER(276)
      KH(11,iangg,ij,ik)=                        RAER(278)
      KH(12,iangg,ij,ik)=            RAER(279) + RAER(280)
      KH(13,iangg,ij,ik)= RAER(282)            + RAER(281)
      KH(14,iangg,ij,ik)= RAER(284)            + RAER(283)

      RETURN
      END
