C-*-Fortran-*-  -  this is used for emacs/Fortran mode capability, or just do M-x fortran-mode
C
C     CONSTANTS AND SWITCHES COMMON BLOCKS


      COMMON /CONST1_1/ DZ,DY,DT,DTH,DY1,DY2,DZ1,DZ2,DELT,C00,GZ
     1 ,DYX,DZX,XNORM,DAYL,YP(N$),YY(N$),ZP(M$),ZZ(M$),RHO(M$)
     2 ,RSTAR(M$),h25,YPP(NP$),ZPP(MP$),XEP(MP$),XN2ST(MP$)
     3 ,EZ2H(M$),S(N$),C(N$),DTHDZ(M$),THG(M$),TG(M$),TTOTH(M$)
     4 ,CF(N$),BETA(N$),XMM(N$,M$),PRNDTL,CST(NP$),XN2(MP$),zbreak
     5 ,W1,W2,ZRES,SST(NP$),CFST(NP$),THERMYN,THERMCR,RHOST(MP$)
     6 ,TGX(MP$), THGX(MP$), TTOTHX(MP$), rfric(NP$)

      COMMON /CONST2_1/ PI,R,H,CP,TOMEG,A,XKAPPA,IPRINT,NDBUG

      COMMON /SWITC0/ ISW(NSWTCHS), LH_SW(12), KZZ_SW(15)

      COMMON /PRANDT/ KYYPR(15), IYYPR(15), KZZPR(15), IZZPR(15)
     2               ,DMOLPR(15)
      REAL            KYYPR,     IYYPR,     KZZPR,     IZZPR
