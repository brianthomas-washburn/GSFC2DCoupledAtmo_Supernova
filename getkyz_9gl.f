C
C
          SUBROUTINE GETKYZ
c
C
c  GETKYZ_9ap.f - routine TO compute Kyz from input Kyy, Kzz, and TEMPs,
C     Kyz computed at the gridbox corners for use in AER adapted diffusive flux routiune (NEWDIFYZ)
C
C  This ensures that the diffusion matrix is always positive (diffusive) 
C        such that, Kyy*Kzz > Kyz*Kyz       (EF, 8/13/03)
C                                                                                            

        include "com2d.h"


 	SAVE


        REAL dptdz(L$+1,Z$X1), dptdy(L$+1,Z$X1), tslope(L$+1,Z$X1)
        REAL kzzmin(L$+1,Z$X1), pt(L$,Z$X), KYZ1(L$+1,Z$X1)
        REAL EKYYC(L$+1,Z$X1), EKYYC1(L$+1,Z$X1), EKZZC(L$+1,Z$X1)



C  first compute potential temperature, TEMP(L$,Z$X) is REAL*8, PRESS(Z$X)

      do 33 ik=1,Z$X
      do 33 ij=1,L$
 33      pt(ij,ik) = TEMP(IJ,IK)*(1013./press(ik))**.286
C
C
C
C  *******************************     KYZ   CALCULATION   *********************************************
C
C
cc *** pot. temp. gradients calculated below (K/meters),     Kyz = -Kyy*(dptdy/dptdz) 
c
C    should compute dz using actual scale height in meters (not just 7000 m), ie, thickness of layer
C    assume the temperature at that level is the mean temperature in the layer
C
C    Kyz is evaluated at grid box corners, so take average theta grads (in both y and z)
C       Kyz = 0. at sides (poles) and bottom and top levels (so DONT compute GRADS at sides/bottom/top)
C
C     dptdz(L$+1,Z$X1), dptdy(L$+1,Z$X1), TEMP(L$,Z$X), PT(L$,Z$X), DZZ1, DZZ2 are in METERS
C
      do 567 ik=2,Z$X
      do 567 ij=2,L$
        dzz1 = 29.3*(TEMP(ij-1,ik-1) + TEMP(ij-1,ik))/2. 
     >             *ALOG(press(ik-1)/press(ik))

        dzz2 = 29.3*(TEMP(ij,ik-1) + TEMP(ij,ik))/2. 
     >             *ALOG(press(ik-1)/press(ik))

        dptdz(ij,ik) = ((pt(ij-1,ik) - pt(ij-1,ik-1))/dzz1 + 
     >                  (pt(ij,ik) - pt(ij,ik-1))/dzz2)/2.
 567  CONTINUE
c
C   delyc is in METERS (defined in STREAMF) for the constituent grid
C     need to interpolate dptdy to proper grid box corner with irregular grid, ZALT(Z$X), ZALTE(Z$X1)
C     don't need to do thickness dz, since we're multiplying back dz, so it should be only a very small error
C
       do 211 ik=2,Z$X
       do 211 ij=2,L$
          dpdy2  = (pt(ij,ik) - pt(ij-1,ik))/delyc
          dpdy1 =  (pt(ij,ik-1) - pt(ij-1,ik-1))/delyc

          dptdy(ij,ik) = (dpdy2 - dpdy1)/(zalt(ik) - zalt(ik-1))*
     >                   (zalte(ik) - zalt(ik-1)) + dpdy1
 211   CONTINUE
c
C
c
c   load into output arrays for DYNOUT;   DPTDYALL(L$,Z$X,NMON$), DPTDZALL(L$,Z$X,NMON$) are in COMMON
C
c        DO 522 IK=1,Z$X
c	DO 522 IJ=1,L$
c          DPTDYALL(ij,ik,im324) = dptdy(ij,ik) 
c 522      DPTDZALL(ij,ik,im324) = dptdz(ij,ik) 
c
C
C  *****  Kyz  CALCULATIONS  -- they are dependent on Kyy, Initialize first   *******************  
C       Kyz = -Kyy*(dptdy/dptdz) ,   KYZALL(L$,Z$X,NMON$)
C
cef	DO 322 IK=1,Z$X
cef	DO 322 IJ=1,L$
cef 322       KYZALL(ij,ik,im324) = 1.e5
c
c
C
c  since we're dividing by dpdz, set limit so that Kyz doesn't blow up here, set to .001 K/m 
C
C  The Kyz calculations are from the  equation, Kyz = Kzy = -Kyy*[(dO/dy)/(dO/dz)] from 
C  Newman et al. [1988], where O=potential temp. There may be an uncertainty of 30% in slope to
C  reduce the Kyz and make sure diffusion matrix is positive (diffusive) since the 
C  relation Kyy*Kzz/Kyz^2 > 1 must hold to avoid NUMERICAL PROBLEMS;  the uncertainty in the mixing slope 
C  is probably around +-30% or so  (P. Newman, 1/97), since mixing is APPROXIMATELY along isentropes.
C  Therefore, the Kzz (computed below) in lower strat won't have to be so large 
C  and H2O will be better simulated
C
C   IN PRACTICE, I found that reducing the slope by at least 25% to as much as 50% gave better results 
C    in CH4, H2O, O3 at 15-25 km. - EF, 8/2005
C
C  dptdz(L$+1,Z$X1), dptdy(L$+1,Z$X1), tslope(L$+1,Z$X1), EKYZ(L$+1,Z$X1)
C  Kyz is evaluated at BOX CORNERS, so need to interpolate Kyy to same 
C  load into array EKYYC(L$+1,Z$X1),  EKYY(L$+1,Z$X)
C  also do the same for Kzz to check that the diffusion matrix is positive, EKZZC(L$+1,Z$X1), EKZZ(L$,Z$X1)
C
C
       DO 4242 IK=2,Z$X
       DO 4242 IJ=2,L$
         if (dptdz(ij,ik) .le. .001) dptdz(ij,ik) = .001
c
C                            reduce slope to 70% of orginal gives better results
c
         tslope(ij,ik) = dptdy(ij,ik)/dptdz(ij,ik)*.7  ! *.75

C  need to properly interpolate Kyy for irregular vertical grid

         EKYYC(ij,ik) = (EKYY(ij,ik) - EKYY(ij,ik-1))/
     >  (zalt(ik) - zalt(ik-1))*(zalte(ik) - zalt(ik-1)) + EKYY(ij,ik-1)

         EKYYC1(ij,ik) = EKYYC(ij,ik)

         EKZZC(ij,ik) = (EKZZ(ij-1,ik) + EKZZ(ij,ik))/2.
 4242   CONTINUE

C
C  Compute Kyz:  if Kzz > Kyy*slope^2, then compute Kyz as usual,
C  if NOT, then adjust 1) slope (up to 50%), 2) Kyy,  BUT DON'T CHANGE KZZ
C  in troposphere, adjust slope up to 0% (ie, no Kyz) to maintain large Kyy
C  for coupled model tropospheric temps which can give near vertical theta (July 2009)
C
C    adjust slope to 20% of original - WCONV131
C
       DO 477 IK=2,Z$X
       DO 477 IJ=2,L$
         kyyslp1 = EKYYC(ij,ik)*(tslope(ij,ik))**2

         slff = .2
         if (zalt(ik) .le. 9.) slff = .1
         if (zalt(ik) .le. 8.) slff = 0.

         tslope1 = slff*tslope(ij,ik)
C                                            first adjust slope up to 20% of original - WCON131
         xslope = tslope(ij,ik)
         if (EKZZC(ij,ik) .LT. kyyslp1) then
             xslope = SQRT(EKZZC(ij,ik)/EKYYC(ij,ik))*.95
             xslope = SIGN(xslope, tslope(ij,ik))                     ! ensure SIGN of original tslope(ij,ik)
             if (ABS(xslope) .LE. ABS(tslope1)) xslope = tslope1
         endif
c                                                     if still a problem, reduce Kyy appropriately
         kyyslp2 = EKYYC(ij,ik)*(xslope**2)
         if (EKZZC(ij,ik) .LT. kyyslp2) 
     >          EKYYC(ij,ik) = .95*EKZZC(ij,ik)/(xslope**2)
                                      

         EKYZ(ij,ik) = -xslope*EKYYC(ij,ik)
                                            ! load final xslope into COMMON array THSLP(L$+1,Z$X1) for output
         THSLP(ij,ik) = xslope
 477   CONTINUE


C  set Kyz = 0. at sides (poles) and bottom and top levels , EKYZ(L$+1,Z$X1)

       do 457 ik=1,Z$X1
          EKYZ(1,ik) = 0.0
 457      EKYZ(L$+1,ik) = 0.0

       do 467 ij=1,L$+1
          EKYZ(ij,1) = 0.0
 467      EKYZ(ij,Z$X1) = 0.0


C  store original Kyz values for use in adjustment below, KYZ1(L$+1,Z$X1)
C
       do 476 ik=1,Z$X1
       do 476 ij=1,L$+1
 476     KYZ1(ij,ik) = EKYZ(ij,ik)

C
C  If Kyy has been modified, interpolate back to get EKYY(L$+1,Z$X) from the REDUCED array EKYYC(L$+1,Z$X1)
C    and original EKYYC1(L$+1,Z$X1)  -  levels 1 and Z$X can just stay the same
c
C   this is OK for irregular vertical grid, since we're interpolating from grid box top/bottom to center
C                                            which is always equidistant by definition
C
       DO 377 IK=2,Z$X-1
       DO 377 IJ=2,L$
          kfact1 = EKYYC(ij,ik)/EKYYC1(ij,ik)
          kfact2 = EKYYC(ij,ik+1)/EKYYC1(ij,ik+1)
            if (kfact1 .lt. 1.  .or.  kfact2 .lt. 1.) then
               kfact = AMIN1(kfact1, kfact2)
               EKYY(ij,ik) = kfact*(EKYYC1(ij,ik) + EKYYC1(ij,ik+1))/2.
            endif
 377   CONTINUE


C
C  Final check: at box centers, check that Kyy*Kzz > Kyz^2 (positive diffusive matrix at grid box centers)
C     if not, then reduce Kyz (ie, the slope) at the surrounding corners
C
C   EKYY(L$+1,Z$X), EKZZ(L$,Z$X1), EKYZ(L$+1,Z$X1),  EKYYCE, EKZZCE, EKYZCE  - values at box centers
C     this is OK for irregular vertical grid, since we're interpolating from grid box edges/corners to center
C

       DO 577 IK=1,Z$X
       DO 577 IJ=1,L$
          EKYYCE = (EKYY(ij,ik) + EKYY(ij+1,ik))/2.
          EKZZCE = (EKZZ(ij,ik) + EKZZ(ij,ik+1))/2.
          EKYZCE = (EKYZ(ij,ik) + EKYZ(ij+1,ik) + EKYZ(ij,ik+1) 
     >                   + EKYZ(ij+1,ik+1))/4.

          difb = EKYYCE*EKZZCE - EKYZCE**2
C                                                      reduce Kyy*Kzz magnitude by 75% to really nail it
          if (difb .le. 0.) then 
              kyzmag = .75*SQRT(EKYYCE*EKZZCE)
c                                                              ! now reduce Kyz at the 4 surrounding corners
              EKYZ(ij,ik) = AMIN1( ABS(EKYZ(ij,ik)), 
     >                  ABS(kyzmag/ekyzce*KYZ1(ij,ik)) )
              EKYZ(ij,ik) = SIGN(EKYZ(ij,ik), KYZ1(ij,ik))            ! reset to sign of original Kyz

              EKYZ(ij+1,ik) = AMIN1( ABS(EKYZ(ij+1,ik)), 
     >                  ABS(kyzmag/ekyzce*KYZ1(ij+1,ik)) )
              EKYZ(ij+1,ik) = SIGN(EKYZ(ij+1,ik), KYZ1(ij+1,ik))            ! reset to sign of original Kyz

              EKYZ(ij,ik+1) = AMIN1( ABS(EKYZ(ij,ik+1)), 
     >                  ABS(kyzmag/ekyzce*KYZ1(ij,ik+1)) )
              EKYZ(ij,ik+1) = SIGN(EKYZ(ij,ik+1), KYZ1(ij,ik+1))            ! reset to sign of original Kyz

              EKYZ(ij+1,ik+1) = AMIN1( ABS(EKYZ(ij+1,ik+1)), 
     >                  ABS(kyzmag/ekyzce*KYZ1(ij+1,ik+1)) )
              EKYZ(ij+1,ik+1) = SIGN(EKYZ(ij+1,ik+1), KYZ1(ij+1,ik+1))        ! reset to sign of original Kyz
          END IF

 577   CONTINUE


C  reset Kyz = 0. at sides (poles) and bottom and top levels , EKYZ(L$+1,Z$X1), LATEG(L$+1), ZALTE(Z$X1)

       do 757 ik=1,Z$X1
          EKYZ(1,ik) = 0.0
 757      EKYZ(L$+1,ik) = 0.0

       do 767 ij=1,L$+1
          EKYZ(ij,1) = 0.0
 767      EKYZ(ij,Z$X1) = 0.0

c  also need to reduce Kyz in troposphere below ~8 km at high lats (factor of 2 works) 
C       - to avoid Cly->0 and NaNs

       do 777 ik=1,Z$X1
       do 777 ij=1,L$+1
          if (ABS(LATEG(ij)) .ge. 70.  .and.  ZALTE(ik) .le. 8.5) 
     >                 EKYZ(ij,ik) = EKYZ(ij,ik)/5.
 777   CONTINUE

c
C  Now need to check that the relation, Kyy*Kzz/Kyz^2 > 1 holds true. Since Kyz = ~-Kyy*slope, then
C    Kzz must be GE Kyy*slope^2, but to ensure that the matrix is positive (diffusive), we'll add an 
C    uncertainty factor of 5% to Kzz. This seemed to be necessary for H2O which has 
C    sharp vertical grads at 15-20 km in the tropics where the Kzz is very small.  
C    Don't do in region above level 46 where slope is questionable.  KZZMIN, KZZTOT are in cm2/sec
C
C    Since max slope is +-1.e-3 in upper troposphere/lower stratosphere, and Kzz > Kyy*slope^2
C    but if Kyy is ~1.e9 or 1.e10 cm2/sec as it is at mid-latitudes, then Kzz ~.1 - 1 m2/sec which
C    gives anomolous transport of H2O across tropopause, because of this, we have set Kyz=0 
C    (and get better total ozone). 
C
C    EKZZ(L$,Z$X1) is the total Kzz in cm2/sec,  kzzmin(L$+1,Z$X1), tslope(L$+1,Z$X1)
c
c
ckzz        DO 425 IK=2,Z$X
ckzz        DO 425 IJ=2,L$
ckzz           kzzmin(ij,ik) = EKYYC(ij,ik)*1.05*(tslope(ij,ik))**2
ckzz 425    CONTINUE
ckzz
ckzz
C                             don't need to reset the ground or top, since Kyz=0. there
ckzz       DO 430 IK=2,Z$X
ckzz
ckzz         DO 434 IJ=2,L$-1
ckzz             kzzmina = (kzzmin(ij,ik) + kzzmin(ij+1,ik))/2.
ckzz             if (EKZZ(ij,ik) .lt. kzzmina) EKZZ(ij,ik) = kzzmina 
ckzz 434     CONTINUE
c                                                                          do poles separately here
ckzz         if (EKZZ(1,ik) .lt. kzzmin(2,ik)) EKZZ(1,ik) = kzzmin(2,ik) 
ckzz         if (EKZZ(L$,ik) .lt. kzzmin(L$,ik)) EKZZ(L$,ik) = kzzmin(L$,ik)
ckzz
ckzz 430   CONTINUE
ckzz

	RETURN
	END

