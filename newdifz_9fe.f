c***************************************************************
c 
c this subroutine will calculate the diffusion terms using a 
c mass conserving differencing scheme
c***************************************************************
C
c   ******  FOR  BASE9  (new grid resolution)   **********
C
C    everything in REAL*8  -  LOOP THROUGH EACH LATITUDE  
c

       SUBROUTINE NEWDIFZ


       include "com2d.h"


       REAL*8 mk(Z$X1), qk(Z$X1), hsc, DPH(Z$X1)
       REAL*8 term1, term11, term12
       REAL*8 tempin1(Z$X), tempk(Z$X1), kzzk(Z$X1)


       SAVE


C  DPH(Z$X1) is the DELTAP between grid points, defined at box edges, where top and bottom = 0 (not used)
C     PRESS8(Z$X) is REAL*8 in COMMON
C
	DO 108 IK=1,Z$X1-2
 108		DPH(IK+1) = DLOG(press8(ik+1)/press8(ik))

        DPH(1) = 0.0
        DPH(Z$X1) = 0.0
C

       DO 1000 IJ=1,L$

C
C  first need to interpolate Temperature to the BOX EDGES, then compute M and Q at the BOX EDGES
C    ZALTE8(Z$X1), PRESSE8(Z$X1) are REAL*8,  INPUT: TEMP(L$,Z$X), ZALT8(Z$X) is REAL*8
C
          do 801 ik=1,Z$X
 801         tempin1(ik) = TEMP(ij,ik)

                                                                           ! tempk(Z$X1)
            CALL LINTERP8(ZALT8, Z$X, TEMPIN1, ZALTE8, Z$X1, 0, TEMPK)


          DO 802 IK=1,Z$X1
              kzzk(ik) = ekzz(ij,ik)
              mk(ik) = PRESSE8(IK)*1.d3/(WTA*RD*TEMPK(IK)*1.66d-24)
              HSC = 29.3*TEMPK(IK)*1.D-3
 802          QK(IK) = -1.d0/(HSC*1.d5)


C  ST(T$,L$,Z$X), DP(Z$X), Q(L$,Z$X), difnz(T$,L$,Z$X) all are REAL*8,  
C       KZZK(Z$X1) is REAL*8 defined at the box edges


       DO 600 IT=1,ITRANS

        term1 = q(ij,1)/dp(1)*mk(2)*kzzk(2)*qk(2)*
     >               (st(IT,IJ,2)-st(IT,IJ,1))/dph(2)
        difnz(it,ij,1) = -1.d0*(term1)

       do 800 IK=2,Z$X-1
              term11 = q(ij,ik)/dp(ik)*mk(IK+1)*kzzk(IK+1)*qk(IK+1)
     >                       *(st(IT,IJ,ik+1) - st(IT,IJ,ik))/DPH(IK+1)

              TERM12 = q(ij,ik)/dp(ik)*mk(IK)*kzzk(IK)*qk(IK)
     >                       *(st(IT,IJ,ik) - st(IT,IJ,ik-1))/dph(ik)

              difnz(it,ij,IK) = -1.d0*(term11-term12)
800     CONTINUE

              term1=-1.*q(ij,Z$X)/dp(Z$X)*mk(Z$X)*kzzk(Z$X)*qk(Z$X)
     >                       *(st(IT,IJ,Z$X) - st(IT,IJ,Z$X-1))/dph(Z$X)

              difnz(it,ij,Z$X) = -1.d0*(term1)

 600  CONTINUE
 1000 CONTINUE


	RETURN
	END

