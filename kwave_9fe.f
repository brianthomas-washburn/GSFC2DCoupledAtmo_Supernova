          SUBROUTINE KWAVE(im10)
c
c   FILE KWAVE_9AA_9200.F  FILE  TO  COMPUTE  momentum deposition from thermally
c      damped Kelvin wave for generation of SAO residual circulation/tracer transport
c      (ie, double peak structure) in 2D model following Gray and Pyle (1987)  (EF 6/13/96)
c
C   ALL CALCULTIONS DONE IN MKS UNITS  ==>  AND DOUBLE PRECISION, for Z$S LEVELS
C
C    now includes QBO forcing, with latitudinal variations,
C       adapted from kwave_base8wa_9200_lat.pro (7/00)


       include "com2d.h"


       REAL*8 aa3, kk3, cc3, vv3, betaf(L$S), yf(L$S), yyl, ydis
       REAL*8 xlatr(L$S), fwave(L$S,Z$S)


       SAVE


C    ik0 = bottom boundary for calcs, set to around 16 km (~100 mb), level 9 for original grid
C                                                        (need to add 1 since we're now starting at 0 km)
       ik0 = INT(16./116.*Z$S)+1


c                                      convert degrees latitude to radians, initialize yf array
       do 101 ij=1,L$S
          yf(ij) = 0.0
          betaf(ij) = 2.*om1*dcosd(latst(ij))/rad
 101      xlatr(ij) = latst(ij)*qp1

c                                  !beta is 0 at poles (cos = 0) - we now will not have points right at pole
CEF          betaf(1) = 0.d0
CEF          betaf(37) = 0.d0


C  *****************************************************************************************
C
C     FIRST DO FAST KELVIN WAVE 
C 
c  Define vertical momentum flux at tropopause, zonal wave# (k=1), and phase speed, from base8zg model
c
c       aa3 = 1.5*7.0d-3
c       kk3 = 1./rad
c       cc3 = 50.d0
c       vv3 = kk3*cc3

c  Define vertical momentum flux at tropopause, zonal wave# (k=1), and phase speed, 
C      these values are from GP89 - gives better SAO signature than base8zg model

       aa3 = 3.5d-3
       kk3 = 1./rad
       cc3 = 60.d0
       vv3 = kk3*cc3

cc  - define latitudinal weighting factor, ydis is distance from equator in m
CC    we assume there is only forcing at 30S-30N, poleward of this fwave=0.0, yf(L$S)

       do 440 ij=1,L$S
          IF (LATST(ij) .ge. -30. .and. LATST(ij) .le. 30.) then
             ydis = qp1*rad*DABS(latst(ij))
             yyl = DSQRT(2.*vv3/(kk3*betaf(ij)))
             yf(ij) = DEXP(-(ydis**2)/(yyl**2))
             if (yf(ij) .lt. 1.d-20) yf(ij) = 0.0d0
          ENDIF
 440   CONTINUE


       CALL MOMDEP(aa3, kk3, cc3, yf, ik0, fwave)


cc  load in fwave(L$S,Z$S) array , also if just using equator/symmetric about equator, do HERE, FK(L$S,Z$S) 
cc      load in above 16 km 
cc     and set limit of < 43.2 m/s/day (0.0005 m/sec^2) to keep it from getting really large above 100 km
C

       do 500 ik=ik0,Z$S
       do 500 ij=1,L$S
          if (LATST(ij) .ge. -30. .and. LATST(ij) .le. 30.) then
              fk(ij,ik) = fwave(ij,ik)
ccccc                   fk(ij,ik) = fwave(19,ik)*yf(ij)
              if (fk(ij,ik) .ge. 0.0005) fk(ij,ik) = 0.0005
          endif
 500   CONTINUE


cc GP87 way of defining latitude dependence
Cc
cc       do 450 ik=8,59
cc       do 450 ij=1,37
cc          if (ij .ne. 19) then
cc             yf = DEXP(-2.*om1*rad*(xlatr(ij)**2)/(cc3 - ubar1(19,ik)))
cc             if (yf .lt. 1.d-20) yf = 0.0d0
cc             fk(ij,ik) = fk(19,ik)*yf
cc          endif
cc 450   continue



C  *****************************************************************************************
C
C     NOW SLOW KELVIN WAVE #1,  from Dunkerton, 1997  (amp*5, following GP89)
C
c  Define vertical momentum flux at tropopause, zonal wave# (k=2), and phase speed - from D97

       aa3 = 12.5d-3
       kk3 = 2./rad
       cc3 = 30.d0
       vv3 = kk3*cc3

cc  - define latitudinal weighting factor, ydis is distance from equator in m
CC    we assume there is only forcing at 30S-30N, poleward of this fwave=0.0, yf(L$S)

       do 540 ij=1,L$S
          IF (LATST(ij) .ge. -30. .and. LATST(ij) .le. 30.) then
             ydis = qp1*rad*DABS(latst(ij))
             yyl = DSQRT(2.*vv3/(kk3*betaf(ij)))
             yf(ij) = DEXP(-(ydis**2)/(yyl**2))
             if (yf(ij) .lt. 1.d-20) yf(ij) = 0.0d0
          ENDIF
 540   CONTINUE


       CALL MOMDEP(aa3, kk3, cc3, yf, ik0, fwave)


cc  load in fwave(L$S,Z$S) array , also if just using equator/symmetric about equator, do HERE, FKQ(L$S,Z$S)
cc     and set limit of < 43.2 m/s/day (0.0005 m/sec^2) to keep it from getting really large above 100 km
CC     where it might hit a critical line

       do 600 ik=ik0,Z$S
       do 600 ij=1,L$S
          if (LATST(ij) .ge. -30. .and. LATST(ij) .le. 30.) then
              fkq(ij,ik) = fwave(ij,ik)
ccccc                  fkq(ij,ik) = fwave(19,ik)*yf(ij)
              if (fkq(ij,ik) .ge. 0.0005) fkq(ij,ik) = 0.0005
           endif
 600   CONTINUE



C
C  *****************************************************************************************
C
C     NOW SLOW KELVIN WAVE #2,  from Dunkerton, 1997  (use 2X amp of Kelvin wave #1 above)
C
c  Define vertical momentum flux at tropopause, zonal wave# (k=2), and phase speed - from D97

       aa3 = 25.d-3
       kk3 = 2./rad
       cc3 = 20.d0
       vv3 = kk3*cc3

cc  - define latitudinal weighting factor, ydis is distance from equator in m
CC    we assume there is only forcing at 30S-30N, poleward of this fwave=0.0, yf(L$S)

       do 740 ij=1,L$S
          IF (LATST(ij) .ge. -30. .and. LATST(ij) .le. 30.) then
             ydis = qp1*rad*DABS(latst(ij))
             yyl = DSQRT(2.*vv3/(kk3*betaf(ij)))
             yf(ij) = DEXP(-(ydis**2)/(yyl**2))
             if (yf(ij) .lt. 1.d-20) yf(ij) = 0.0d0
          ENDIF
 740   CONTINUE


       CALL MOMDEP(aa3, kk3, cc3, yf, ik0, fwave)


cc  load in fwave(L$S,Z$S) array , also if just using equator/symmetric about equator, do HERE, FKQ2(L$S,Z$S)
cc     and set limit of < 43.2 m/s/day (0.0005 m/sec^2) to keep it from getting really large above 100 km
CC     where it might hit a critical line

       do 617 ik=ik0,Z$S
       do 617 ij=1,L$S
          if (LATST(ij) .ge. -30. .and. LATST(ij) .le. 30.) then
              fkq2(ij,ik) = fwave(ij,ik)
ccccc                  fkq2(ij,ik) = fwave(19,ik)*yf(ij)
              if (fkq2(ij,ik) .ge. 0.0005) fkq2(ij,ik) = 0.0005
           endif
 617   CONTINUE



C  *****************************************************************************************
C
C     ROSSBY-GRAVITY WAVE, from Dunkerton, 1997  (amp*6, following GP89)
C
c  Define vertical momentum flux at tropopause, zonal wave# (k=4), and phase speed - from D97
c     need to make the momentum flux at tropopause negative for westward moving waves

       aa3 = -15.d-3
       kk3 = 4./rad
       cc3 = -40.d0
       vv3 = kk3*cc3

cc  - define latitudinal weighting factor, ydis is distance from equator in m
CC    we assume there is only forcing at 30S-30N, poleward of this fwave=0.0, yf(L$S)

       do 820 ij=1,L$S
          IF (LATST(ij) .ge. -30. .and. LATST(ij) .le. 30.) then
            ydis = qp1*rad*DABS(latst(ij))
            yyl = DSQRT(2.*vv3/(betaf(ij)*(betaf(ij)/vv3 - kk3)))
            yf(ij) = DEXP(-(ydis**2)/(yyl**2))
            if (yf(ij) .lt. 1.d-20) yf(ij) = 0.0d0
          ENDIF
 820   CONTINUE


       CALL MOMDEP(aa3, kk3, cc3, yf, ik0, fwave)


cc   load in fwave(L$S,Z$S) array , also if just using equator/symmetric about equator, do HERE, FRQ(L$S,Z$S)
cc     and set limit of > -130 m/s/day (-0.0015 m/sec^2) to keep it from getting really large above 65 km
CC     where it might hit a critical line
CC
       do 700 ik=ik0,Z$S
       do 700 ij=1,L$S
          if (LATST(ij) .ge. -30. .and. LATST(ij) .le. 30.) then
              frq(ij,ik) = fwave(ij,ik)
ccccc                  frq(ij,ik) = fwave(19,ik)*yf(ij)
              if (frq(ij,ik) .le. -0.0015) frq(ij,ik) = -0.0015
          endif
 700   CONTINUE

       RETURN
       END


c
c   ROutine to compute momentum deposition from tropical waves
C      based on Gray and Pyle, 1989
C
C   ALL CALCULTIONS DONE IN MKS UNITS  ==>  AND DOUBLE PRECISION, for Z$S LEVELS
C

       SUBROUTINE MOMDEP(aa3, kk3, cc3, yf, ik0, fwave)


       include "com2d.h"


       REAL*8 ab1, ab2, aa3, kk3, cc3, sum1, xa, nnn(L$S,Z$S)
       REAL*8 yf(L$S), icr0(L$S), succ(Z$S), fwave(L$S,Z$S)
       REAL*8 rr0(L$S,Z$S), pp0(L$S,Z$S)

       REAL*4 thd59(59), thd(Z$S)


       SAVE


c  Slow damping rates as defined in Dunkerton, 1979
c          ab1 = 7.5257d-4
c          ab2 = .69897d0
c
c  Fast damping rates as defined in Dunkerton, 1979
c       ab1 = 1.7474d-3
c       ab2 = .30103d0

c  Modified Fast damping rates (Not quite as fast as the fast rates above)
       ab1 = 1.5d-3
       ab2 = .43103d0


c define thermal damping profile  (in 1/sec), set to zero in troposphere               ;; zz59 in km
c
c   NOTE: THERE are some small differences between the FK, FKQ, and FRQ values
C      computed with this FORTRAN CODE vs. the  kwave_base8wa_9200_lat.pro IDL routine
C      at levels 16-59 (in FORTRAN). This is likely due to the thd array 
C      definition here for levels 16-59 which is a complicated function of EXPs and such
C      so there are probably some difference here between IDL and FORTRAN. However for 
C      levels 1-15, thd is just a constant, and the IDL and FORTRAN Values are identical 
C      to within 1.e-9 (just machine accuracy for I/O).

       do 120 ik=1,5
 120      thd59(ik) = 0.0

       do 121 ik=6,15
 121      thd59(ik) = 0.4e-6 + 0.8e-6*((zz59(ik)-17.)/13.)

       do 122 ik=16,59
 122     thd59(ik) = 1./(86400.*EXP(2.3*(ab1*(zz59(ik)-50.)**2 + ab2)))

c                        re-define below 30 km to be constant value, as defined in GP89
       do 123 ik=1,15
 123      thd59(ik) = 0.35e-6

C  
C  interpolate (logarithmically) to proper grid, ZSTR4(Z$S) is REAL*4  as are THD(Z$S), thd59(59)
C
       CALL LINTERP(ZZ59, 59, THD59, ZSTR4, Z$S, 1, THD)  

c
c  Seems like better results in HF FOR ALL SEASONS are achieved with a modified FAST damping rate,
c     and 1.5 times the Gray and Pyle tropopause mom fluxes (theirs was 7.e-3), as listed below
c     this seems a bit better overall compared with the Base8qk and Base8ql runs
C     - THIS WAS FOR THE SAO ONLY 
c
c  if COMPUTing R (rr0)  and  P (pp0) AT EQUATOR only, and damping out in latitude, need to 
C              rescale below


       do 150 ik=1,Z$S
       do 150 ij=1,L$S
          rr0(ij,ik) = 0.0
          pp0(ij,ik) = 0.0
          icr0(ij) = Z$S
          fwave(ij,ik) = 0.0
          nnn(ij,ik) = DSQRT(xn2(ij,ik))
 150   continue

C                                                 UBAR1(L$S,Z$S)
       DO 350 ij=1,L$S 
          do 351 ik=1,Z$S
 351         succ(ik) = 0.0

          do 352 ik=1,Z$S
             if (ubar1(ij,ik)-cc3 .lt. 0.) succ(ik) = -1.
             if (ubar1(ij,ik)-cc3 .ge. 0.) succ(ik) = 1.
 352      CONTINUE

         do 353 ik=ik0,Z$S-1
           if (DABS(ubar1(ij,ik)-cc3) .lt. 1.0  .or.  
     >                          succ(ik) .ne. succ(ik+1)) then
             icr0(ij) = ik  
             GO TO 350
           endif
 353     continue
 350   CONTINUE


       do 201 ij=1,L$S
       do 201 ik=1,icr0(ij)
 201      rr0(ij,ik) = thd(ik)*nnn(ij,ik)/(kk3*((ubar1(ij,ik)-cc3)**2))


       do 300 ij=1,L$S
       do 300 ik=ik0,icr0(ij)
           sum1 = 0.0d0
              do 305 ii = ik0-1,ik-1
                 xa = (rr0(ij,ii) + rr0(ij,ii+1))/2.
                 sum1 = sum1 + xa*delz
 305          continue
          pp0(ij,ik) = sum1
 300  continue 

C          ;;  FWAVE(L$S,Z$S) is du/dt (in m/sec2), eq (6) in GP87, also define the 14 km index for zstr
C
       ik14 = INT(14./116.*Z$S)+1

       do 400 ij=1,L$S
          if (LATST(ij) .ge. -30. .and. LATST(ij) .le. 30.) then
             do 405 ik=ik0,icr0(ij)
               fwave(ij,ik) = aa3*DEXP((zstr(ik)-zstr(ik14))*1000./hh) 
     >                        *rr0(ij,ik)*DEXP(-pp0(ij,ik))*yf(ij)
 405         CONTINUE
          ENDIF
 400   CONTINUE


	RETURN
	END


