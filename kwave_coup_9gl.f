C
          SUBROUTINE KWAVE_COUP(NP$, MP$, yp, ypp, zp, zpp, ttin, uuin, 
     >                          FK, FKQ, FKQ2, FRQ, FEGWC)

c
c   FILE KWAVE_COUP.F  FILE  TO  COMPUTE  momentum deposition from thermally
c      damped Kelvin wave for generation of SAO residual circulation/tracer transport
c      (ie, double peak structure) in 2D model following Gray and Pyle (1987)  (EF 6/13/96)
c
C    now includes QBO forcing, with latitudinal variations,
C       adapted from kwave_base8wa_9200_lat.pro (7/00)
C
C    also includes Dunkerton's PARAMs for equatorial gravity waves
C
C     
C     THIS  is adopted for the coupled model (converted from cm/sec) to MKS
C
C     ALL CALCULTIONS DONE IN MKS UNITS  ==>  AND DOUBLE PRECISION
C

CCCCCC       include "com2d.h"

!        USE degree_trig

        INTEGER NP$, MP$, N$, M$
        REAL*4 yp(NP$-1), ypp(NP$), zp(MP$-1), zpp(MP$)
        REAL*4 ttin(NP$-1,MP$-1), tt(NP$,MP$), uuin(NP$,MP$)

        REAL*8 delz, rad, om1, rr1, hh, kap, xpi, qp1, tt8(NP$,MP$)
        REAL*8 aa3, kk3, cc3, vv3, betaf(NP$), yf(NP$), yyl, ydis
        REAL*8 dtdz(NP$,MP$), xn2(NP$,MP$), ubar1(NP$,MP$)
        REAL*8 xlatr(NP$), fwave(NP$,MP$)
 
        REAL*4 fk(NP$,MP$), fkq(NP$,MP$), fkq2(NP$,MP$), frq(NP$,MP$)
        REAL*4 fegwc(NP$,MP$)


        SAVE

C                                         get delz in Meters
        delz = (zpp(2)- zpp(1))*1000.d0

        rad = 6371000.d0
        om1 = 7.292d-5
        rr1 = 287.d0
        hh = 7000.d0
        kap = 2./7.d0 
c               qp1=pi/180. in REAL*8, dthst=latitude resolution in deg for STREAMF; dely, delz are in meters

        xpi = 4.*DATAN(1.d0)
        qp1 = xpi/180.d0


C    ik0 = index of bottom boundary for calcs, set to around 16 km (~100 mb)
C
        ik0 = 8

         do 877 ik=1,MP$
             if (zpp(ik) .lt. 16.) ik0 = ik
 877     CONTINUE



C  interpolate INPUTE MODEL temperatures to expanded grid:  tt(NP$,MP$), get DT/DZ, XN2  dtdz(NP$,MP$)

          N$ = NP$-1
          M$ = MP$-1

          CALL BINTERP(yp, N$, zp, M$, TTIN, 
     >                 ypp, NP$, zpp, MP$, 0, 0, TT)

C                                                           convert to REAL*8 for DERV4, XN2
          do 110 ik=1,MP$
          do 110 ij=1,NP$
 110         tt8(ij,ik) = tt(ij,ik)


          CALL DERV4(1, TT8, dtdz, delz, NP$, MP$, 2)


          do 130 ik=1,MP$
          do 130 ij=1,NP$
            xn2(ij,ik) = rr1/hh*(dtdz(ij,ik) + kap*tt8(ij,ik)/hh)
              if (xn2(ij,ik) .lt. 1.d-4)  xn2(ij,ik) = 1.d-4
 130  continue


c                                      convert degrees latitude to radians, initialize yf array
       do 101 ij=1,NP$
          yf(ij) = 0.0
          betaf(ij) = 2.*om1*COSD(ypp(ij))/rad
 101      xlatr(ij) = ypp(ij)*qp1

c                                  !beta is 0 at poles (cos = 0) - we now will not have points right at pole
CEF          betaf(1) = 0.d0
CEF          betaf(37) = 0.d0


C  initialize final momentum arrays,  CONVERT  UBAR -> m/sec,  ubar1(NP$,MP$)

          do 717 ik=1,MP$
          do 717 ij=1,NP$
             ubar1(ij,ik) = uuin(ij,ik)/100.

             fk(ij,ik) = 0.
             fkq(ij,ik) = 0.
             fkq2(ij,ik) = 0.
             frq(ij,ik) = 0.
             fegwc(ij,ik) = 0.
 717  CONTINUE


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

       aa3 = 1.d-3         ! 3.5d-3
       kk3 = 1./rad
       cc3 = 60.d0
       vv3 = kk3*cc3

cc  - define latitudinal weighting factor, ydis is distance from equator in m
CC    we assume there is only forcing at 30S-30N, poleward of this fwave=0.0, yf(NP$)

       do 440 ij=1,NP$
          IF (YPP(ij) .ge. -30. .and. YPP(ij) .le. 30.) then
             ydis = qp1*rad*ABS(ypp(ij))
             yyl = DSQRT(2.*vv3/(kk3*betaf(ij)))
             yf(ij) = DEXP(-(ydis**2)/(yyl**2))
             if (yf(ij) .lt. 1.d-20) yf(ij) = 0.0d0
          ENDIF
 440   CONTINUE


       CALL MOMDEPC(NP$, MP$, ypp, zpp, xn2, ubar1, delz,
     >              aa3, kk3, cc3, betaf, yf, ik0, fwave)


cc  load in fwave(NP$,MP$) array , also if just using equator/symmetric about equator, do HERE, FK(NP$,MP$) 
cc      load in above 16 km 
cc     and set limit of < 43.2 m/s/day (0.0005 m/sec^2) to keep it from getting really large above 100 km
C

       do 500 ik=ik0,MP$
       do 500 ij=1,NP$
          if (YPP(ij) .ge. -30. .and. YPP(ij) .le. 30.) then
              fk(ij,ik) = fwave(ij,ik)
ccccc                   fk(ij,ik) = fwave(19,ik)*yf(ij)
cckw              if (fk(ij,ik) .ge. 0.0005) fk(ij,ik) = 0.0005
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

       aa3 = 1.d-3    ! 2.5d-3
       kk3 = 2./rad
       cc3 = 30.d0
       vv3 = kk3*cc3

cc  - define latitudinal weighting factor, ydis is distance from equator in m
CC    we assume there is only forcing at 30S-30N, poleward of this fwave=0.0, yf(NP$)

       do 540 ij=1,NP$
          IF (YPP(ij) .ge. -30. .and. YPP(ij) .le. 30.) then
             ydis = qp1*rad*ABS(ypp(ij))
             yyl = DSQRT(2.*vv3/(kk3*betaf(ij)))
             yf(ij) = DEXP(-(ydis**2)/(yyl**2))
             if (yf(ij) .lt. 1.d-20) yf(ij) = 0.0d0
          ENDIF
 540   CONTINUE


       CALL MOMDEPC(NP$, MP$, ypp, zpp, xn2, ubar1, delz,
     >              aa3, kk3, cc3, betaf, yf, ik0, fwave)



cc  load in fwave(NP$,MP$) array , also if just using equator/symmetric about equator, do HERE, FKQ(NP$,MP$)
cc     and set limit of < 43.2 m/s/day (0.0005 m/sec^2) to keep it from getting really large above 100 km
CC     where it might hit a critical line

       do 600 ik=ik0,MP$
       do 600 ij=1,NP$
          if (YPP(ij) .ge. -30. .and. YPP(ij) .le. 30.) then
              fkq(ij,ik) = fwave(ij,ik)
ccccc                  fkq(ij,ik) = fwave(19,ik)*yf(ij)
cckw              if (fkq(ij,ik) .ge. 0.0005) fkq(ij,ik) = 0.0005
           endif
 600   CONTINUE



C
C  *****************************************************************************************
C
C     NOW SLOW KELVIN WAVE #2,  from Dunkerton, 1997  (use 2X amp of Kelvin wave #1 above)
C
c  Define vertical momentum flux at tropopause, zonal wave# (k=2), and phase speed - from D97

       aa3 = 1.d-3    ! 5.d-3
       kk3 = 2./rad
       cc3 = 20.d0
       vv3 = kk3*cc3

cc  - define latitudinal weighting factor, ydis is distance from equator in m
CC    we assume there is only forcing at 30S-30N, poleward of this fwave=0.0, yf(NP$)

       do 740 ij=1,NP$
          IF (YPP(ij) .ge. -30. .and. YPP(ij) .le. 30.) then
             ydis = qp1*rad*ABS(ypp(ij))
             yyl = DSQRT(2.*vv3/(kk3*betaf(ij)))
             yf(ij) = DEXP(-(ydis**2)/(yyl**2))
             if (yf(ij) .lt. 1.d-20) yf(ij) = 0.0d0
          ENDIF
 740   CONTINUE


       CALL MOMDEPC(NP$, MP$, ypp, zpp, xn2, ubar1, delz,
     >              aa3, kk3, cc3, betaf, yf, ik0, fwave)


cc  load in fwave(NP$,MP$) array , also if just using equator/symmetric about equator, do HERE, FKQ2(NP$,MP$)
cc     and set limit of < 43.2 m/s/day (0.0005 m/sec^2) to keep it from getting really large above 100 km
CC     where it might hit a critical line

       do 617 ik=ik0,MP$
       do 617 ij=1,NP$
          if (YPP(ij) .ge. -30. .and. YPP(ij) .le. 30.) then
              fkq2(ij,ik) = fwave(ij,ik)
ccccc                  fkq2(ij,ik) = fwave(19,ik)*yf(ij)
cckw              if (fkq2(ij,ik) .ge. 0.0005) fkq2(ij,ik) = 0.0005
           endif
 617   CONTINUE



C  *****************************************************************************************
C
C     ROSSBY-GRAVITY WAVE, from Dunkerton, 1997  (amp*6, following GP89)
C
c  Define vertical momentum flux at tropopause, zonal wave# (k=4), and phase speed - from D97
c     need to make the momentum flux at tropopause negative for westward moving waves

       aa3 = -2.d-3
       kk3 = 4./rad
       cc3 = -40.d0
       vv3 = kk3*cc3

cc  - define latitudinal weighting factor, ydis is distance from equator in m
CC    we assume there is only forcing at 30S-30N, poleward of this fwave=0.0, yf(NP$)

       do 820 ij=1,NP$
          IF (YPP(ij) .ge. -30. .and. YPP(ij) .le. 30.) then
            ydis = qp1*rad*ABS(ypp(ij))
            yyl = DSQRT(2.*vv3/(betaf(ij)*(betaf(ij)/vv3 - kk3)))
            yf(ij) = DEXP(-(ydis**2)/(yyl**2))
            if (yf(ij) .lt. 1.d-20) yf(ij) = 0.0d0
          ENDIF
 820   CONTINUE


       CALL MOMDEPC(NP$, MP$, ypp, zpp, xn2, ubar1, delz,
     >              aa3, kk3, cc3, betaf, yf, ik0, fwave)


cc   load in fwave(NP$,MP$) array , also if just using equator/symmetric about equator, do HERE, FRQ(NP$,MP$)
cc     and set limit of > -130 m/s/day (-0.0015 m/sec^2) to keep it from getting really large above 65 km
CC     where it might hit a critical line
CC
       do 700 ik=ik0,MP$
       do 700 ij=1,NP$
          if (YPP(ij) .ge. -30. .and. YPP(ij) .le. 30.) then
              frq(ij,ik) = fwave(ij,ik)
ccccc                  frq(ij,ik) = fwave(19,ik)*yf(ij)
cckw              if (frq(ij,ik) .le. -0.0015) frq(ij,ik) = -0.0015
          endif
 700   CONTINUE

C
C
C  now do small scale graivty waves following Dunkerton - adapted from offline IDL code
C
C 
cckw       CALL MOMDEPC_EQGW(NP$, MP$, ypp, zpp, xn2, ubar1, delz, 
cckw     >                   qp1, rad, ik0, fegwc)
cckw          

C                                   convert momentum dep back to cm/sec^2 before returning
          do 757 ik=1,MP$
          do 757 ij=1,NP$
             fk(ij,ik)   = fk(ij,ik)*100.
             fkq(ij,ik)  = fkq(ij,ik)*100.
             fkq2(ij,ik) = fkq2(ij,ik)*100.
             frq(ij,ik)  = frq(ij,ik)*100.
             fegwc(ij,ik) = fegwc(ij,ik)*100.
 757  CONTINUE


       RETURN
       END


c
c   ROutine to compute momentum deposition from tropical waves
C      based on Gray and Pyle, 1989
C
C   ALL CALCULTIONS DONE IN MKS UNITS  ==>  AND DOUBLE PRECISION, for MP$ LEVELS
C

        SUBROUTINE MOMDEPC(NP$, MP$, ypp, zpp, xn2, ubar1, delz,
     >                     aa3, kk3, cc3, betaf, yf, ik0, fwave)


cccccccccc       include "com2d.h"


        COMMON/CC59/zz59(59), thd59(59)


        INTEGER NP$, MP$
        REAL*8 delz, hh, ab1, ab2, aa3, kk3, cc3, sum1, xa, nnn(NP$,MP$)
        REAL*8 yf(NP$), icr0(NP$), succ(MP$), fwave(NP$,MP$), betaf(NP$)
        REAL*8 rr0(NP$,MP$),pp0(NP$,MP$), xn2(NP$,MP$),ubar1(NP$,MP$),bb

        REAL*4 thd(MP$), ypp(NP$), zpp(MP$)

 
        SAVE


        hh = 7000.d0

C                           define zz59 in km
        do 109 ik=1,59 
 109	  zz59(ik) = 7.*.2844*(ik-1.)


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
C  interpolate (logarithmically) to proper grid, ZPP(MP$) is REAL*4  as are THD(MP$), thd59(59)
C
       CALL LINTERP(ZZ59, 59, THD59, ZPP, MP$, 1, THD)  

c
c  Seems like better results in HF FOR ALL SEASONS are achieved with a modified FAST damping rate,
c     and 1.5 times the Gray and Pyle tropopause mom fluxes (theirs was 7.e-3), as listed below
c     this seems a bit better overall compared with the Base8qk and Base8ql runs
C     - THIS WAS FOR THE SAO ONLY 
c
c  if COMPUTing R (rr0)  and  P (pp0) AT EQUATOR only, and damping out in latitude, need to 
C              rescale below


       do 150 ik=1,MP$
       do 150 ij=1,NP$
          rr0(ij,ik) = 0.0
          pp0(ij,ik) = 0.0
          icr0(ij) = MP$
          fwave(ij,ik) = 0.0
          nnn(ij,ik) = DSQRT(xn2(ij,ik))
 150   continue

C                                                 UBAR1(NP$,MP$)
       DO 350 ij=1,NP$ 
          do 351 ik=1,MP$
 351         succ(ik) = 0.0

          do 352 ik=1,MP$
             if (ubar1(ij,ik)-cc3 .lt. 0.) succ(ik) = -1.
             if (ubar1(ij,ik)-cc3 .ge. 0.) succ(ik) = 1.
 352      CONTINUE

         do 353 ik=ik0,MP$-1
           if (DABS(ubar1(ij,ik)-cc3) .lt. 1.0  .or.  
     >                          succ(ik) .ne. succ(ik+1)) then
             icr0(ij) = ik  
             GO TO 350
           endif
 353     continue
 350   CONTINUE

C  
C                             rr0 is different for Kelvin and R-G waves
      if (aa3 .ge. 0.) then

         do 201 ij=1,NP$
         do 201 ik=1,icr0(ij)
 201       rr0(ij,ik) = thd(ik)*nnn(ij,ik)/(kk3*((ubar1(ij,ik)-cc3)**2))

      else

         do 202 ij=1,NP$
         do 202 ik=1,icr0(ij)
           bb = betaf(ij)/(kk3*kk3*(ubar1(ij,ik)-cc3)) - 1.
        rr0(ij,ik) = thd(ik)*nnn(ij,ik)/(kk3*((ubar1(ij,ik)-cc3)**2))*bb
 202     CONTINUE

      endif



       do 300 ij=1,NP$
       do 300 ik=ik0,icr0(ij)
           sum1 = 0.0d0
              do 305 ii = ik0-1,ik-1
                 xa = (rr0(ij,ii) + rr0(ij,ii+1))/2.
                 sum1 = sum1 + xa*delz
 305          continue
          pp0(ij,ik) = sum1
 300  continue 

C          ;;  FWAVE(NP$,MP$) is du/dt (in m/sec2), eq (6) in GP87, also define the 14 km index for zpp(MP$)

         ik14 = 7

         do 777 ik=1,MP$
             if (zpp(ik) .lt. 14.) ik14 = ik
 777     CONTINUE


       do 400 ij=1,NP$
          if (YPP(ij) .ge. -30. .and. YPP(ij) .le. 30.) then
             do 405 ik=ik0,icr0(ij)
               fwave(ij,ik) = aa3*DEXP((zpp(ik)-zpp(ik14))*1000./hh) 
     >                        *rr0(ij,ik)*DEXP(-pp0(ij,ik))*yf(ij)
 405         CONTINUE
          ENDIF
 400   CONTINUE


	RETURN
	END


c
c   ROutine to compute momentum deposition from tropical gravity waves, based on Dunkerton, 1997
C      adapted from offline IDL code,   - computes gravity wave spectrum of Dunkerton, 1997
C
C
C   ALL CALCULTIONS DONE IN MKS UNITS  ==>  
C
c
C    THIS IS IN SINGLE PRECISION: - results are nearly identical to DOUBLE PRECISION 
C                                   (just machine accuracy differences, ~e-6 to e-8)


        SUBROUTINE MOMDEPC_EQGW(NP$, MP$, ypp, zpp, xn2, ubar1, delz, 
     >                          qp1, rad, ik0, fegw)


        COMMON/CC59/zz59(59), thd59(59)


        INTEGER NP$, MP$
        REAL*8 delz, qp1, rad, nnn(MP$), xn2(NP$,MP$), ubar1(NP$,MP$)

        REAL*4 thd(MP$), ypp(NP$), zpp(MP$), zdunk(60), thdunk(60)
        REAL*4 ccc(121), aaa(121), ygg(NP$), zgg(MP$), thd1(39)
        REAL*4 ucc(MP$), ucc2(MP$), fgw1(MP$)
        REAL*4 zexp(MP$), rr0(MP$), pp0(MP$), fwave(MP$,12)
        REAL*4 fegw1(NP$,MP$), fegw(NP$,MP$), kk3


        DATA thd1/.05,  .05,  .05,  .05,  .05,  .049, .045, .039, .036, 
     >            .037, .039, .042, .047, .054, .062, .072, .081, .091, 
     >            .10, .109,  .115, .120, .122, .124, .126, .1275, .126, 
     >            .124, .122, .120, .115, .109, .10,  .091, .081, .072, 
     >            .062, .05, .035/

 
        SAVE


C                           find index for 55 km, or 63 km, as in Dunk, 97
        ik55 = 27
         do 677 ik=1,MP$
             if (zpp(ik) .lt. 63.) ik55 = ik
ccccc             if (zpp(ik) .lt. 55.) ik55 = ik
 677     CONTINUE


         do 777 ik=1,MP$
 777        zexp(ik) = EXP((zpp(ik)-zpp(ik0))/7.)


C   initialize final arrays:

          do 577 ik = 1,MP$
          do 577 ij = 1,NP$
             fegw1(ij,ik) = 0.0
 577         fegw(ij,ik) =  0.0



C  note: B0 = 25.e-3 is DUNK, 1997

        b0 = 10.e-3   ! 25.e-3    ! /2.5, reset Dunk, 1997
                                                            ! sign func: result = abs(1st)*sgn(2nd)
        do 101 i=1,121
           ccc(i) = (i-1) - 60.
           c1 = SIGN(1., ccc(i) )
           aaa(i) = 35.*c1/20.*ABS(ccc(i))/10.*EXP(-ABS(ccc(i))/10.)
 101    CONTINUE


C    define latitude, altitude masks:

        do 110 ij=1,NP$
          yydis = qp1*6371.*ABS(ypp(ij))
          ygg(ij) = EXP(-(yydis/2000.)**2)
          if (ABS(ypp(ij)) .gt. 40.) ygg(ij) = 0.0
 110    CONTINUE


        do 120 ik=1,MP$
           zgg(ik) = 1.
           if (zpp(ik) .ge. 35.) zgg(ik) = EXP((35.-zpp(ik))/3.49999)
           if (zpp(ik) .ge. 63.) zgg(ik) = 0.0
C                                                            reducing above 25 km doesn't seem to help much
ccc            if (zpp(ik) .ge. 25.) zgg(ik) = EXP((25.-zpp(ik))/7.)
ccc            if (zpp(ik) .ge. 50.) zgg(ik) = EXP((50.-zpp(ik))/3.49999)
cccc           if (zpp(ik) .ge. 55.) zgg(ik) = 0.0
 120    CONTINUE



C from Dunk, 1997, Fig. 1, middle curve, in 1/days, convert to 1/sec
C    just do mirror image for 52-76 km, then add on from above for >76 km
C    NOTE: the Dunk, 1997 THD profile gives a bit weaker wave driving, but the same patterns.....
C

        do 130 ik=1,39
           zdunk(ik) = (ik-1)*2.
 130       thdunk(ik) = thd1(ik)/86400.
  
        do 131 ik=40,60
           zdunk(ik) = (ik-1)*2.-1.
 131       thdunk(ik) = thd59(ik-1)
        
C                                                                      thd(MP$)
        CALL LINTERP(ZDUNK, 60, THDUNK, ZPP, MP$, 0, THD)  


C                                      loop through latitudes, only do 40S-40N:
       DO 100 ij=1,NP$

          if (ABS(ypp(ij)) .le. 40.) then

                                                   ! define NNN, initialize fgw1(MP$)
          do 150 ik=1,MP$
             nnn(ik) = DSQRT(xn2(ij,ik))*thd(ik)
 150         fgw1(ik) = 0.0


C  for each phase speed, define cc3, and aa3 (Bottom amplitude), for c=0, set aa3 to the value for c=1 m/s
C    Note that the factor of 35 is needed to get Dunk's (1997) non-dim profile in Figure 3
C     then apply his bottom flux, B0 = 25e-3/2.5
C
C
       do 200 ICC = -60,60

C                                              !  initialize fwave(MP$,12)
           do 225 ikk=1,12
           do 225 ik=1,MP$
 225          fwave(ik,ikk) = 0.0

                                                                    ! sign func: result = abs(1st)*sgn(2nd)
           cc3 = REAL(icc)
           c1 = SIGN(1., cc3)

           aa3 = b0*35.* c1/20.*ABS(cc3)/10.*EXP(-ABS(cc3)/10.)   
cccc           if (icc eq 0) aa3 = b0*35./20.*1./10.*EXP(-1./10.)    ! don't use, as in Dunk97

C                                                               ;  define U-C, (if 0, set to small positive)
           do 210 ik=1,MP$
              ucc(ik) = UBAR1(ij,ik) - cc3
              if (ucc(ik) .eq. 0.) ucc(ik) = .01
              ucc2(ik) = ucc(ik)**2
 210       CONTINUE


C  find critical level at each latitude,
C  fast method: load into uu(16:116), find the minimum of: 
C    1) lowest level where ABS(u<1)
C    2) lowest level where sgn changes w/ that of level above
C    3) top level (in case there's no critical level)
C
      
               i1 = MP$
               icr = MP$

               do 221 ik = ik0,MP$
                  if (ABS(ucc(ik)) .lt. 1.) then
                    icr = ik
                    go to 220
                  endif
 221           CONTINUE

               do 222 ik = MP$-1,ik0,-1
                  s1 = SIGN(1., ucc(ik))
                  s2 = SIGN(1., ucc(ik+1))
                  if (s1 .ne. s2) i1 = ik
 222           CONTINUE

                icr = MIN(i1, MP$)
 220         CONTINUE


C  now wave number loop: kk3 = 1, 2, 4, 8, 16,...., 1024, 2048

C  for CPU speed, just compute rr0, pp0 at ALL POINTS, then confine to below critical levels in FGW1 array
C                              rr0(MP$), pp0(NP$,MP$), delz is in meters, just go up to 55 km


           do 310 IKKN=1,12

              do 315 ik=1,MP$
 315             pp0(ik) = 0.0


              kk3 = (2.**(ikkn-1.))/rad

              do 320 ik=1,ik55
 320             rr0(ik) = nnn(ik)/(kk3*ucc2(ik))

              do 330 ik=ik0, ik55
 330             pp0(ik) = pp0(ik-1) + (rr0(ik-1) + rr0(ik))/2.*delz


C  fwave is du/dt (in m/sec2), eq (6) in GP87;  distribute flux (aa3) equally among the 12 zonal wavenumbers
C                                                      ! fwave(MP$,12)
              do 340 ik=ik0, ik55 
 340             fwave(ik,ikkn) = aa3/12.*zexp(ik)*rr0(ik)*EXP(-pp0(ik))

 310       CONTINUE


C  now sum up over all wavenumbers, and limit by critical level, which is a function of phase speed, cc3
C                                                           fgw1(MP$) is ALSO summed up over ALL PHASE SPEEDS
            do 400 ik = ik0, icr
                 do 410 ikkn=1,12
 410                fgw1(ik) = fgw1(ik) + fwave(ik,ikkn)
 400        CONTINUE

C                              end loop over PHASE SPEEDS
 200     CONTINUE


C  load into latitude array, apply masks, fegw1, fegw(NP$,MP$)
C
            do 450 ik = ik0,ik55
 450           fegw1(ij,ik) = fgw1(ik)

          ENDIF

 100    CONTINUE
C                   end latitude loop


          do 700 ik = 1,MP$
          do 700 ij = 1,NP$
 700         fegw(ij,ik) = fegw1(ij,ik)*ygg(ij)*zgg(ik)      ! *2.5/4., reset to Dunk, 1997


	RETURN
	END
