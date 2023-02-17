          SUBROUTINE STREAMF(IM10, IIYR, IM324) 
C
C   FILE STREAMF_9AL.F   FORTRAN ROUTINE TO  COMPUTE  2D MODEL RESID CIRC FROM TEMPS, HEAT RATES,
C     BY SOLVING STREAM FUNCTION AS IN THE GARCIA-SOLOMON MODEL, AND THE 2D COUPLED MODEL
C     SOLVES BY  SUCCESSIVE OVER-RELAXATION, SOLVES FOR 72 5-DAY MEANS
C
C   THIS ALSO INCLUDES EXTRA X (CHI) TERM NOT CONTAINED IN GS83 FORMULATION, AND ALSO SMALL MODIFICATIONS
C   TO THE CZZ, AND CZ COEFFICIENTS. HOWEVER, THESE CHANGES MADE ONLY VERY SMALL DIFFERENCES IN THE 
C   CHI, W, AND V FIELDS. NOTE ALSO THAT THE SOLUTION CONVERGED A LITTLE FASTER WITH THESE CHANGES
C
C  NOTE: USING CONSTANT N2 HAS VERY LITTLE EFFECT EXCEPT IN TROPICAL TROPOSPHERE WHERE
C   CIRC CELLS ARE 2-3X BIGGER WITH N2 FROM TEMPS (WITH OLD HEATING RATES)
C
C   *****  NOTE  THE  COPIOUS  DOCUMENTATION  ****************
C
C   ALL CALCULTIONS DONE IN MKS UNITS  ==>   V*, W*, KYY, KZZ CONVERTED TO CM(CM2)/SEC
C    ***   USES DELY AND DELZ, AND RHOST, IDENTICAL TO THAT USED IN TRANSPORT ROUTINES ***
C    ***  INCLUDES ROUTINE KWAVE FOR SAO ******
C    ****  ADD IN  EDDY HEAT FLUX TERM FOR PWAVES IN THERMODYNAMIC EQUATION  *********
C
C    *** EVERYTHING IS IN DOUBLE PRECISION (REAL*8)   ***********
C
C    *****  NOW  FOR  VARIABLE  RESOLUTION   ****************


       include "com2d.h"


       INTEGER ITERC    !!, nt1
C                                                    IF IT DOESN'T CONVERGE AFTER 2000 ITERATIONS, IT STOPS
       PARAMETER(ITERC=2000)


       REAL*8 kyy(L$+1,Z$X)

       REAL*8 srr, sf1, sf2, sf3, sf4, sf6, sf7, sf8, sf9, sf0, sf5
       REAL*8 adiffc, dtest, kk1, kk2
       REAL*8 ff1(L$S,Z$S), ff2(L$S,Z$S), kru(L$S,Z$S), qvad
       REAL*8 kry(L$S,Z$S), kry1(L$S,Z$S), df1dy(L$S,Z$S), kyysh
       REAL*8 heat1(L$S,Z$S), heat1s(L$S,Z$S), lh1(L$S,Z$S)
       REAL*8 temp1(L$S,Z$S), dtt1(L$S,Z$S), dt1y(L$S,Z$S), ehf(L$S,Z$S)
       REAL*8 ddh(L$S,Z$S), ddhs(L$S,Z$S), ehfz(L$S,Z$S)
       REAL*8 fxall(L$S,Z$S), dfxdz(L$S,Z$S), qv(L$S,Z$S)
       REAL*8 qtot(L$S,Z$S), qtots(L$S,Z$S), qtots2(L$S,Z$S)
       REAL*8 czz(L$S,Z$S), cz(L$S,Z$S), czy(L$S,Z$S), cy(L$S,Z$S)
       REAL*8 cyy(L$S,Z$S), c0(L$S,Z$S), chif(L$S,Z$S)
       REAL*8 cff(L$S,Z$S), cffs(L$S,Z$S), chi0(L$S,Z$S,2)
       REAL*8 dcdy(L$S,Z$S), dcdz(L$S,Z$S), diffc(L$S,Z$S)
       REAL*8 dudy(L$S,Z$S), dtdy(L$S,Z$S), dtdz(L$S,Z$S), lhex(L$S,Z$S)
       REAL*8 dt1z(L$S,Z$S), dt1yy(L$S,Z$S), dhdy(L$S,Z$S)                      !! , wchi(L$,Z$S)
       REAL*8 bbc2(L$S,10), bbcy(L$S,10), wbot(L$S), wbot2(L$S),costot37

       REAL*4 zzink(L$S,Z$S), zzoutk(L$,Z$X1), fegw(L$S,Z$S), yyg(L$S)
       REAL*4 zz91(L$S,Z$S), zzout(L$,Z$X), zzout1(L$+1,Z$X)

       REAL*4 CHI_REAN(L$S,21,72,33), LAT7(L$S), ZZZ7(21)


C                                        ! WHACK72 are in K/sec (WACCM HACK latent heating)
       COMMON/CCHACK/WHACK72(91,117,72)

C
c
C srr=successive over-relaxation parameter: 1= stan (no over-relxtion) 1<srr<2; BEST SRR SEEMS TO BE 1.7-1.75
C
       DATA srr/1.7d0/

       SAVE


C
C  Read in tropospheric X* (BBC) from NCEP data for stream function;  READ IN 1ST TIME THROUGH ONLY
C    already on model streamfunction latitude grid, every 0.5 km for 0-10 km    
c    CHI_REAN(L$S=91,21,72,33), LAT7(L$S=91), ZZZ7(21) are defined above (NOT in COMMON)  - in m^2/sec 
c                                       iy=1,32 is for 1979-2010;   33 = CLIMATOLOGY
       IF (im324 .eq. 1) then 

          OPEN (75, FILE = '/misc/kah04/fleming/bbc05/chi_7910_j.xdr',
     >                      convert="big_endian", form = 'binary', 
     >                      access = 'sequential', status = 'old')

            READ (75) CHI_REAN
            READ (75) LAT7
            READ (75) ZZZ7
          CLOSE (75)

       ENDIF



C       INITIALIZE PHYSICAL CONSTANTS for Streamfunction routines, all in COMMON

       rad = 6371000.d0
       om1 = 7.292d-5
       hh = 7000.d0
       kap = 2./7.d0 
c               qp1=pi/180. in REAL*8, dthst=latitude resolution in deg for STREAMF; dely, delz are in meters

       xpi = 4.*DATAN(1.d0)
       qp1 = xpi/180.d0

       dely = rad*qp1*dthst
       dely2 = dely*dely
C                              delz is in COMMON, now defined as a PARAMETER (in meters)
       delz2 = delz*delz

c                            delyc is the DELTA-Y (in meters) for the constituent grid in COMMON - key on dth
       DELYC = rad*qp1*dth

c                                          ! yyg(L$S) is for the equatorial gravity wave momentum source
       do 103 ij=1,L$S
          yydis = qp1*6371.*DABS(latst(ij))
          yyg(ij) = EXP(-(yydis/2000.)**2)
          if (DABS(latst(ij)) .gt. 40.) yyg(ij) = 0.0
ccccc          yyg(ij) = 1.0
c                                                            ! set yyg=1 everywhere if not using mask
          cose(ij) = dcosd(latst(ij))
          tana(ij) = dtand(latst(ij))/rad  
 103      ff0(ij) = 2.*om1*dsind(latst(ij))

c  redefine COSE and TANA at poles, although they're NEVER USED below, 

          cose(1) = 0.d0
          cose(L$S) = 0.d0

          tana(1) = dtand(-89.9d0)/rad   
          tana(L$S) = dtand(89.9d0)/rad   

C
C  now load in/convert dynamics arrays to proper units:  epin in units of m/sec2 => Dont need to convert
C                                                        heatin in units of deg K/day => convert to K/sec
C                           ehfin, ehfzin (add together here) in units of deg K/sec => DON'T CONVERT
C                     LHGCM(L$S,Z$S,72) daily values from GCM in units of deg K/day => convert to K/sec
C                                                  Kyyin in units of m2/sec = Don't need to convert
C                                             double check current limit of 5.e10 cm2/s (this can be reduced)
C                                              set minimum Kyy to 1.e8 cm2/s (1.e6 seems to small)
C                                                but note that only 1.e10 is the Kyy limit at the top levels
c                                                  NOTE: Kyy and Kzz are used here only to compute ddh term
c
         do 105 ik=1,Z$S
         do 105 ij=1,L$S
            temp1(ij,ik) = temp5in(ij,ik,im10) 
            ubar1(ij,ik) = ubarin(ij,ik,im10)   !use zonal mean:  UBARin(L$S,Z$S,72)

            qv(ij,ik) = epin(ij,ik,im10)
c                                                           Make Delf>0.0 == 0.0
              if (qv(ij,ik) .gt. 0.) qv(ij,ik) = 0.0
c               if (qv(ij,ik) .gt. 0.) qv(ij,ik) = qv(ij,ik)*cose(ij)

C                            ! lhgcm(ij,ik,im10)/86400. - for MERRA latent heating included in HEATIN
C                            ! but add in WACCM HACK latent heating, WHACK72(91,117,72) is in K/sec
            lh1(ij,ik) = WHACK72(ij,ik,im10)

            ddh(ij,ik) = 0.0
            ehf(ij,ik) = ehfin(ij,ik,im10) + ehfzin(ij,ik,im10)
            heat1(ij,ik) = heatin(ij,ik,im10)/86400.            !! HEATIN(L$S,Z$S,72) is in K/day

c                                                     FEGWIN(L$S,Z$S,72) in COMMON (m/sec^2) - apply yyg(L$S)

            fegw(ij,ik) = 0.    ! fegwin(ij,ik,im10)*yyg(ij)   - SET FEGW = 0. for 9GTB run
 105     continue


C
C  interpolate Kyy to chemistry grid box sides, LATST4(L$S), ZSTR4(Z$S) are in COMMON,  ZZ91(L$S,Z$S)
C    KYYIN(L$S,Z$S,72) is on STREAMF grid, 
C    KYY(L$+1,Z$X) is defined at the chemistry grid box sides and vertical levels, eg, 58 levs
C                                        LATEG4(L$+1), ZALT(Z$X),  ZZOUT1(L$+1,Z$X)

            do 117 ik=1,Z$S 
            do 117 ij=1,L$S
 117	       zz91(ij,ik) = KYYIN(ij,ik,im10)
                                                                             ! ZALT(Z$X)
            CALL BINTERP(LATST4, L$S, ZSTR4, Z$S, zz91, 
     >                   LATEG4, L$+1, ZALT, Z$X, 0, 0, ZZOUT1)

C                                                      load into KYY array, limits are now set in TROPKZZ
            DO 875 ik=1,Z$X
            DO 875 ij=1,L$+1
    	      kyy(ij,ik) = zzout1(ij,ik)
 875        CONTINUE



C   
C   interpolate TEMPERATURE to the CHEMISTRY GRID - TEMPALL(L$,Z$X,74)
C   TEMP5in(L$S,Z$S,72) is on STREAMF grid -   zz91(L$S,Z$S), LAT4(L$), ZALT(Z$X), ZZOUT(L$,Z$X)

            do 517 ik=1,Z$S 
            do 517 ij=1,L$S
 517	       zz91(ij,ik) = TEMP5in(ij,ik,im10) 

            CALL BINTERP(LATST4, L$S, ZSTR4, Z$S, zz91, 
     >                   LAT4, L$, ZALT, Z$X, 0, 0, ZZOUT)

            DO 675 ik=1,Z$X
            DO 675 ij=1,L$
 675	       TEMPALL(ij,ik,im10+1) = zzout(ij,ik)


c
C   linearly interpolate DELF between 20S-20N in tropical troposphere (levels 1-4) 
C     where there is no HRDI data, blend over level 4 to get a smooth transition. 
C     Don't worry about consistency with QBARY, since this is just the TROPOSPHERE!
C     (Kyy in the troposphere now adjusted in TROPKZZ routine)
C
C
CEF - NEED TO RE-DO THIS WITH VARIABLE RESOLUTION??????
CEF         do 112 ik=1,4
CEF         do 112 ij=16,22
CEF            qvad = (qv(23,ik)-qv(15,ik))/8.*(ij-15) + qv(15,ik)
CEF              if (ik .le. 3) qv(ij,ik) = qvad
CEF              if (ik .eq. 4) qv(ij,ik) = (qv(ij,ik) + qvad)/2.
CEF 112     continue
CEF
c

c
c add in extra latent heating mainly in sub-tropics, so Hadley cell is not as strong, and more spread out
c                                           NO, DON'T DO THIS WITH PPM
c         do 390 ik=1,Z$S
c         do 390 ij=1,L$S
c 390          lhex(ij,ik) = 0.0
c
c         do 391 ij=7,31
c 391          lhex(ij,3) = 2.0/86400.*dcosd(3.*(dabs(latst(ij))-30.))
c
c         do 392 ij=9,29
c 392          lhex(ij,4) = 2.0/86400.*dcosd(3.6*(dabs(latst(ij))-25.))
c
c         do 393 ij=11,27
c 393          lhex(ij,5) = 1.35/86400.*dcosd(4.5*(dabs(latst(ij))-20.))
c
c         do 394 ij=13,25
c 394          lhex(ij,6) = .65/86400.*dcosd(6.*(dabs(latst(ij))-15.))
cc                                                                        add to latent heat
c         do 398 ik=1,
c         do 398 ij=1,L$S
c 398        lh1(ij,ik) = lh1(ij,ik) + lhex(ij,ik)
c

C background Rayleigh fric w/ latitude dependence only, 1/RF = ~100 days in tropics, 
C   50 days at mid-high lats, decrease to 1/300 days in troposphere everywhere since we already 
C   have momentum source from qv from the NMC data - convert to 1/sec  
CEF
CEF - Rayleigh friction not used now - set to zero below
CEF
CEF         do 106 ik=1,Z$S
CEF         do 106 ij=1,L$S
CEF 106        kry1(ij,ik) = 150.-100*dexp(-dcosd(latst(ij))**4)
CEF
CEF        do 415 ik=1,7
CEF        do 415 ij=1,L$S 
CEF 415       kry1(ij,ik) = 300.
CEF
CEF        do 416 ij=1,L$S
CEF           kry1(ij,8) = .75*kry1(ij,7) + .25*kry1(ij,11)
CEF           kry1(ij,9) = .5*kry1(ij,7) + .5*kry1(ij,11)
CEF 416       kry1(ij,10) = .25*kry1(ij,7) + .75*kry1(ij,11)
CEF
CEF        do 419 ik=1,Z$S
CEF        do 419 ij=1,L$S
CEF 419       kry(ij,ik) = 1./kry1(ij,ik)/86400.
CEF

c
C   ; Now compute derivatives of  temp1, ubar1, dudy, dudz, dtdy, dtdz = (L$S,Z$S)

       call derv4(1, ubar1, dudy, dely, L$S, Z$S, 1)
       call derv4(1, ubar1, dudz, delz, L$S, Z$S, 2)
       call derv4(1, temp1, dtdy, dely, L$S, Z$S, 1)
       call derv4(1, temp1, dtdz, delz, L$S, Z$S, 2)


C  Brunt-Visala freq as defined in Lindzen (1981), is for a variable scale height, but we
C    just need it from dtheta/dz, in which H is a constant scale height for log p altitude
C    so we should use Randel's definition as below, w/ H=7 km
C
C    ALSO, the 3.2d-4 limit is NO longer necessary since the Cyy coefficient is now just N^2
C          which is always positive (dtheta/dz >0) (N^2 > 2 day-1), but this N^2 change has an 
C          impact on GWD and diffusion which are generally larger, 
C
C     now w/ new MLS data in mesosphere, need to ensure N^2 (xn2) is always positive so that Cyy is 
C          always positive, otherwise streamfunction blows up during summer - EF, June 04
C
       do 130 ik=1,Z$S 
       do 130 ij=1,L$S
          xn2(ij,ik) = rr1/hh*(dtdz(ij,ik) + kap*temp1(ij,ik)/hh)
              if (xn2(ij,ik) .lt. 1.d-4)  xn2(ij,ik) = 1.d-4
 130  continue


C                            call gravity wave routine to compute Kzzt (m2/sec)  and  Fxtt (m/sec2)
C                                                 after dudz, xn2 have been determined

            CALL GWAVE(im10, iiyr, im324) 


C                                 initialize Kelvin/Rossby-Grav. wave arrays here
            do 348 ik=1,Z$S
            do 348 ij=1,L$S 
               fk(ij,ik) = 0.0
               fkq(ij,ik) = 0.0
               fkq2(ij,ik) = 0.0
 348           frq(ij,ik) = 0.0

c                               call Kelvin wave routine to get FK(L$S.Z$S) (m/sec2) momentum forcing for SAO
c                             and FKQ(L$S.Z$S), FKQ2(L$S.Z$S), FRQ(L$S,Z$S) (m/sec2) momentum forcing for QBO
C                                 note: im10 is NOT used in KWAVE (it's just a dummy argument here)
            CALL KWAVE(im10) 

c                              to test the sensitivity of FKQ/FRQ, set = 0.0
cc            do 338 ik=1,Z$S
cc            do 338 ij=1,L$S 
cc               fkq(ij,ik) = 0.0
cc 338           frq(ij,ik) = 0.0




C 
C
C  NOW COMPUTE FORCING TERMS:
C
C  First the Dh term (eq 7) in GS83, eddy diffusion of heat by GWs,  this was redone following
C   Schoeberl et al., 1983  and Huang and Smith, 1991, DDH = 1/rho*d/dz(rho*Kzz*d@/dz) in terms of theta
C   which has been converted to Temp -- this ends up being the same as the DH term in GS83 
C   they had the correct sign after all. UNITS are in (DEG K/sec), as before, NO TROPOSPHERIC DDH -- Sept 97
C   also now includes the PW eddy heating, read in off-line (the EHF+EHFZ array).
C
        do 151 ik=1,Z$S
        do 151 ij=1,L$S 
 151      dtt1(ij,ik)= xkzzt(ij,ik)*(.286/hh*temp1(ij,ik) + dtdz(ij,ik))

       CALL DERV4(1, dtt1, dt1z, delz, L$S, Z$S, 2)

         ik16 = INT(16./116.*Z$S)+1
c                                                only compute DDH above tropopause, NO TROPOSPHERIC DDH 
         do 155 ik=ik16,Z$S
         do 155 ij=1,L$S
 155        ddh(ij,ik) = (.286-1.)/hh*dtt1(ij,ik) + dt1z(ij,ik)

ccc
ccc   OLD Version of DDH which includes Kyy, Kzz, and rho, pressure terms
ccc
ccc        do 151 ik=1,Z$S
ccc        do 151 ij=1,L$S 
ccc            dtt1(ij,ik)= rhost(ik)*xkzzt(ij,ik)
ccc     c                    *((pres59(1)/pres59(ik))**.286)
ccc     c                    *(.286/hh*temp1(ij,ik) + dtdz(ij,ik))
ccc 151    continue
ccc
ccc       call derv4(1, dtt1, dt1z, delz, L$S1, Z$S, 2)
ccc
ccc        do 152 ik=1,Z$S
ccc        do 152 ij=1,L$S
ccc 152       dt1y(ij,ik) = cose(ij)*kyy(ij,ik)*dtdy(ij,ik)
ccc
ccc       call derv4(1, dt1y, dt1yy, dely, L$S1, Z$S, 1)
ccc
ccc         do 155 ik=10,Z$S
ccc         do 155 ij=1,L$S
ccc       ddh(ij,ik) = ((pres59(ik)/pres59(1))**.286)/rhost(ik)*dt1z(ij,ik)
ccc         ddh(ij,ik) = -dt1yy(ij,ik)/cose(ij)                   ! old version includes DDH from Kyy, Kzz
ccc c      + ((pres59(ik)/pres59(1))**.286)/rhost(ik)*dt1z(ij,ik)
ccc 155     continue
ccc
cccc      CALL SMOOTH5(ddh, ddhs, L$S, Z$S)
CCCC
CCCC      CALL SMOOTH5(heat1, heat1s, L$S, Z$S)          ! New heating rates already smoothed 1X offline
c                                            
C
C  TOTAL HEATING  ;; qtot, heat1, lh1, ehf, ddh, gwh1, gwh2  are in K/sec
C     ehf is plan. wave heating (EHF+EHFZ);  ddh, gwh1, gwh2 are gravity wave HEATING;  gwh1, gwh2 (L$S,Z$S)
C
         do 175 ik=1,Z$S
         do 175 ij=1,L$S
           qtot(ij,ik) = heat1(ij,ik) + lh1(ij,ik) + ehf(ij,ik) + 
     >        ddh(ij,ik) + gwh1(ij,ik) + gwh2(ij,ik) 
 175     CONTINUE
C
c                                                       smooth total heating rates 2X
cc          CALL SMOOTH5(qtot, qtots, L$S, Z$S)
cc          CALL SMOOTH5(qtots, qtots2, L$S, Z$S)
c
C
C  TOTAL MOMENTUM FORCING, all terms in  m/sec2:  and include EQ gravity wave momentum source (offline) 
C                                                 - FEGW(L$S,Z$S)
      do 180 ik=1,Z$S 
      do 180 ij=1,L$S
C       kru(ij,ik) = kry(ij,ik)*ubar1(ij,ik)
        kru(ij,ik) = 0.0                                                     ! set Rayleigh Friction = 0 

        fxall(ij,ik) = -kru(ij,ik) + qv(ij,ik) + fk(ij,ik)
     >                  + fkq(ij,ik) + fkq2(ij,ik) + frq(ij,ik)
     >                  + fxtt(ij,ik) + dmom(ij,ik) + fegw(ij,ik)
 180  CONTINUE

C  Compute derivatives of Forcing
C
      call derv4(1, qtot, dhdy, dely, L$S, Z$S, 1)
      call derv4(1, fxall, dfxdz, delz, L$S, Z$S, 2)

C
C  ff1 = ff0 + 2*Ubar*tan(lat)/a   ; ff2 = ff0 + ubar*tan(lat)/a - du/dy    ;;  ff0 = f 
C       df1dy is in analytic form after expanding in spherical coords,  
C
       do 190 ik=1,Z$S
       do 190 ij=1,L$S
          ff1(ij,ik) = ff0(ij) + 2.*ubar1(ij,ik)*tana(ij) 
          ff2(ij,ik) = ff0(ij) + ubar1(ij,ik)*tana(ij) - dudy(ij,ik)
          df1dy(ij,ik) = 2.*om1*cose(ij)/rad + 2.*tana(ij)*dudy(ij,ik)
     >                 + 2.*ubar1(ij,ik)/(rad*rad*cose(ij)*cose(ij))
 190   continue

C    set f's, dfdy at poles so they don't blow up, doesn't matter since they're not used in X* calculation

       do 191 ik=1,Z$S
          ff1(1,ik) = ff0(1)
          ff1(L$S,ik) = ff0(L$S)
          ff2(1,ik) = ff0(1)
          ff2(L$S,ik) = ff0(L$S)
          df1dy(1,ik) = 0.0
 191      df1dy(L$S,ik) = 0.0

C
C  NOW SET UP COEFFICIENTS FOR STREAMFUNCTION EQUATION, note: forcing term includes cos(lat) factor
C      these have been corrected - June 1998 
C  ALSO, CZZ can only be VERY slightly negative or else X* blows up w/ fine resolution - set limit (Feb 02)
C                                              ! actually for pure elliptic eqn, Czz>0
      do 200 ik=1,Z$S
      do 200 ij=1,L$S
        czz(ij,ik) = ff1(ij,ik)*ff2(ij,ik)
          if (czz(ij,ik) .le. -7.e-11) czz(ij,ik) = -7.e-11 

        cz(ij,ik) = -ff1(ij,ik)*ff2(ij,ik)/hh + df1dy(ij,ik)*dudz(ij,ik)
     >              + 2.*ff1(ij,ik)*dudz(ij,ik)*tana(ij)

        czy(ij,ik) = 2.*ff1(ij,ik)*dudz(ij,ik)

        cy(ij,ik) = (xn2(ij,ik) - 2.*dudz(ij,ik)*dudz(ij,ik))*tana(ij) 
     >               - ff1(ij,ik)/hh*dudz(ij,ik)*(1. + kap)
                                                                      !! XN2 and Cyy are now always > 0. 6/04
        cyy(ij,ik) = xn2(ij,ik)

        c0(ij,ik) = -df1dy(ij,ik)*dudz(ij,ik)/hh 
     >              - 2.*ff1(ij,ik)*dudz(ij,ik)*tana(ij)/hh

        cff(ij,ik) = (ff1(ij,ik)*dfxdz(ij,ik) + rr1/hh*dhdy(ij,ik))*
     c                  cose(ij) 
 200   continue
c                                  reset cy at poles (prop. to tan(lat)), but its not used below anyway
        do 203 ik=1,Z$S
            cy(1,ik) = cy(2,ik)
 203        cy(L$S,ik) = cy(L$S-1,ik)

c                                   Don't smooth forcing term, elliptic equation will take care of this
ccc          CALL SMOOTH5(cff, cffs, L$S, Z$S)

        do 798 ik=1,Z$S
        do 798 ij=1,L$S
 798       cffs(ij,ik) = cff(ij,ik)

c damp out forcing for 110-116 km, 0 at top (level Z$S), doesn't help to smooth 
c  can just interpolate using the indicies (Instead of ZSTR) since STREAMF uses an equadistant grid
c
        i110 = INT(110./116.*Z$S)+1

        do 801 ij=1,L$S        
           cffs(ij,Z$S) = 0.0
        do 801 ik=i110+1,Z$S-1
 801       cffs(ij,ik) = cffs(ij,i110)/(i110-Z$S)*(ik-Z$S)


cjune
c         write (77) czz, cz, czy, cy, cyy, c0, cff, cffs
c         write (77) ff1, ff2, df1dy      

c
C              
C  BOUNDARY CONDITIONS, RESET FOR ALL TIME STEPS:     Set X* = 0 at poles, TOP and bottom for rigid walls
C       use G92 eq(9) for LEVEL 2, interpolate between +-15 deg, and adjust for mass balance
C   USE ONLY BBC = 0 AT LEVEL 1, LEVEL 2 IS specified, NOTE: Make BBC2 a 2D array for DERIV4 call
c     fxall(L$S,Z$S), rhost(Z$S), ff0(L$S)
C

      DO 220 ij=2,L$S-1
        sum=0.0
           do 221 ik=2,Z$S
221     sum = sum + rhost(ik)*fxall(ij,ik)*delz

      if (ff0(ij) .ne. 0.) bbc2(ij,1)= -cose(ij)/(ff0(ij)*rhost(2))*sum   
220    CONTINUE                                                         
C                                                     interpolate between 15S-15N
       i15s = INT((-15.-(-90.))/180.*L$S)+1
       i15n = INT((15.-(-90.))/180.*L$S)+1

       do 222 ij=i15s+1,i15n-1   
        bbc2(ij,1) = (bbc2(i15n,1) - bbc2(i15s,1))/(i15n-i15s)*(ij-i15s) 
     >                             + bbc2(i15s,1)
 222   CONTINUE

        bbc2(1,1) = 0.d0
        bbc2(L$S,1) = 0.d0
                                       ! Reduce X* to get a weaker circ., and load in ik=2,10 for DERV4 call
                                       ! only reduce by .8 for 9200r run
        do 780 ik=1,10                 
        do 780 ij=1,L$S
 780       bbc2(ij,ik) = bbc2(ij,1)*.8
C
C  Now adjust for mass balance, first compute w*,  bbc2, bbcy (L$S,10), wbot(L$S), wbot2(L$S) all in meters
C                                                                              dely= in meters

       CALL DERV4(1, bbc2, bbcy, dely, L$S, 10, 1)

        do 781 ij=2,L$S-1
 781       wbot(ij) = bbcy(ij,1)/cose(ij)
c                                            
        wbot(1) = wbot(2)
        wbot(L$S) = wbot(L$S-1)

        costot37 = 0.d0
        sum = 0.0
        do 782 ij=1,L$S
           sum = sum + wbot(ij)*cose(ij)
 782       costot37 = costot37 + cose(ij)

        do 784 ij=1,L$S
 784       wbot2(ij) = wbot(ij) - sum/costot37

C finally integrate to get X*, everything is in meters, seconds, X*=0.0 at poles, of course

        bbc2(1,1) = 0.d0
        bbc2(L$S,1) = 0.d0

        do 785 ij=2,L$S-1
            bbc2(ij,1) = bbc2(ij-1,1) + 
     >         (wbot2(ij-1)*cose(ij-1) + wbot2(ij)*cose(ij))/2.*dely
 785    CONTINUE

       do 223 ij=1,L$S
 223      bbcadj(ij,im10) = bbc2(ij,1)
                                                                  ! REAL*8 BBCADJ(L$S,72) in COMMON

C   Initialize  chi  field, start with zero's everywhere except bottom boundary (ie, level 2) 
C                       chi0(L$S,Z$S,2)

       do 225 itt=1,2
       do 225 ik=1,Z$S 
       do 225 ij=1,L$S
 225      chi0(ij,ik,itt) = 0.d0


C    ITERATION  PROCEDURE, iterate until change from between iterations is <1% at all grid points
C                          Set boundary conditionS for first step, bottom, top = 0.0, specified for lev 2
C
C  - use X* from NCEP data - CHI_REAN(L$S=91,21,72,33) defined ABOVE,  in m^2/sec 
C     iiyr = 1,32 for 1979-2010,  for  CHI_REAN, 1=1979, ...., 32=2010;  33=climatology

          DO 1000 nt1 = 1,iterc

          do 230 ij=1,L$S
             chi0(ij,Z$S,1) = 0.d0
             chi0(ij,1,1) = 0.d0
                                                             ! chi0(L$S,Z$S,2), set here for 1-3 km
             do 231 ik=2,4
 231            chi0(ij,ik,1) = CHI_REAN(ij,ik+ik-1,im10,iiyr)
 230      CONTINUE


cccccc 230         chi0(ij,2,1) = bbcadj(ij,im10)


c                            DON'T solve for X* at poles, top, or bottom, or levels 2-4, these are pre-set
          do 400 ik=5,Z$S-1
          do 400 ij=2,L$S-1

      sf1=chi0(ij-1,ik-1,2)*(czy(ij,ik)*dely*delz)
      sf2=chi0(ij,ik-1,2)*(4.*czz(ij,ik)*dely2 -2.*cz(ij,ik)*dely2*delz)
      sf3=chi0(ij+1,ik-1,2)*(-czy(ij,ik)*dely*delz)
      sf4=chi0(ij-1,ik,2)*(4.*cyy(ij,ik)*delz2 -2.*cy(ij,ik)*dely*delz2)
      sf6=chi0(ij+1,ik,1)*(2.*dely*delz2*cy(ij,ik) +4.*delz2*cyy(ij,ik))
      sf7=chi0(ij-1,ik+1,1)*(-dely*delz*czy(ij,ik))
      sf8=chi0(ij,ik+1,1)*(4.*dely2*czz(ij,ik) +2.*dely2*delz*cz(ij,ik))
      sf9=chi0(ij+1,ik+1,1)*(dely*delz*czy(ij,ik))
      sf0 = 4.*dely2*delz2*cffs(ij,ik)
      sf5 = 8.*dely2*czz(ij,ik) + 8.*delz2*cyy(ij,ik)
     c      - 4.*dely2*delz2*c0(ij,ik)

        chi0(ij,ik,2) = chi0(ij,ik,1) + srr/sf5*(sf1 + sf2 + sf3 + sf4
     c          + sf6 + sf7 + sf8 + sf9 - sf0 - sf5*chi0(ij,ik,1))

 400  continue

C
C  SET Boundaries = 0 for all time steps

       do 260 ik=1,Z$S
           chi0(1,ik,2) = 0.0  
 260       chi0(L$S,ik,2) = 0.0

       do 270 ij=1,L$S
cc         chi0(ij,Z$S,2) = chi0(ij,Z$S-1,2)

           chi0(ij,Z$S,2) = 0.d0
           chi0(ij,1,2) = 0.d0
                                                       ! chi0(L$S,Z$S,2), CHI_REAN(L$S=91,21,72,33)
             do 271 ik=2,4
 271            chi0(ij,ik,2) = CHI_REAN(ij,ik+ik-1,im10,iiyr)
 270   CONTINUE



ccccccccccc 270       chi0(ij,2,2) = bbcadj(ij,im10)

 
C  Now check difference between iterations, after 100 iterations, first initialize each time
C    if maximum difference is <0.1%, end iteration loop, and goto w, v calculation

          if (nt1 .gt. 100) then

       do 280 ik=2,Z$S-1
       do 280 ij=2,L$S-1
              diffc(ij,ik) = 0.0d0

        if (chi0(ij,ik,1) .ne. 0.0) 
     >    diffc(ij,ik)= 100.*(chi0(ij,ik,2)-chi0(ij,ik,1))/chi0(ij,ik,1)
 280  continue

      dtest = 0.0

         do 285 ik=2,Z$S-1 
         do 285 ij=2,L$S-1
            adiffc = DABS(diffc(ij,ik))
 285        dtest = DMAX1(adiffc, dtest)

            if (dtest .lt. 0.1) go to 99
ccccc            if (dtest .lt. 1.0) go to 99
          endif

c                                              if continuing, update chi0 for next iteration, chi0(L$S,Z$S,2)
          do 295 ik=1,Z$S
          do 295 ij=1,L$S
 295         chi0(ij,ik,1) = chi0(ij,ik,2)

cjune
c         write (77) nt1, chi0

 1000  CONTINUE

                                    ! chif(L$S,Z$S) is the FINAL STREAMFUNCTION ARRAY
 99       do 305 ik=1,Z$S 
          do 305 ij=1,L$S
 305         chif(ij,ik) = chi0(ij,ik,2)


c  Compute  v* from stream function (A-grid), just for diagnostics walla8, valla(L$S,Z$S) are IN COMMON

       CALL DERV4(1, chif, dcdz, delz, L$S, Z$S, 2)

          do 307 ik=1,Z$S 
          do 307 ij=2,L$S-1
            valla(ij,ik) = -(dcdz(ij,ik) - chif(ij,ik)/hh)/cose(ij)*100.           ! also convert to cm/sec
 307      continue
                                                                            ! v* is zero at the poles by def.
          do 309 ik=1,Z$S
             valla(1,ik) = 0.0
 309         valla(L$S,ik) = 0.0


c  Compute w* from stream function, A-grid = walla8(L$S,Z$S) is REAL*8   
c
       CALL DERV4(1, chif, dcdy, dely, L$S, Z$S, 1)

          do 310 ik=2,Z$S-1
          do 310 ij=2,L$S-1
 310         walla8(ij,ik) = dcdy(ij,ik)/cose(ij)*100.              ! convert to cm/sec here

          do 311 ij=1,L$S
             walla8(ij,1) = 0.0d0                                    ! top and bottom w* = 0, 
 311         walla8(ij,Z$S) = 0.0d0
                                                     ! at poles, set w*=1 point inside poles, for diagnostics
          do 312 ik=1,Z$S 
             walla8(1,ik) = walla8(2,ik)   
 312         walla8(L$S,ik) = walla8(L$S-1,ik)

C
ccccc       do 315 ik=1,Z$S 
ccccc       do 315 ij=1,L$
ccccc 315    wchi(ij,ik) = (chif(ij+1,ik) - chif(ij,ik))/dely/cosc(ij)*100.d0 
cccc
C 
c  Now interpolate w* from streamfunction grid to C-grid, for constituent transport,  WALLA8(L$S,Z$S) 
C     and at top and bottom edges of boxes corresponding to ZSTR(Z$S), LATST(L$S)
C     ie, interpolate to the constituent grid box edges (L$,Z$X1) before adjusting,
C      ! chif(L$S,Z$S),  and convert to CM/SEC
C     in REAL*8   LAT(L$),  ZALTE8(Z$X1), W1(L$,Z$X1) in COMMON - this interpolation has been checked - OK   

            CALL BINTERP8(LATST, L$S, ZSTR, Z$S, WALLA8, 
     >                    LAT, L$, ZALTE8, Z$X1, 0, 0, W1)

c
c  subroutine vwmass does mass adjustment at each p-level, and computes v* by simple mass continuity 
c    NOTE:  THIS IS DONE HERE JUST FOR THE OUTPUT FILES (WALL, VALL). VWMASS is then called again 
C           every day in SETDAILY, after w* is interpolated to daily values
C
C  NOW, (in case of irregular grid), you need to also give it:
C  1) the w* field;  2) the altitude grid you're using;  3) the number of levels (# latitudes is ALWAYS L$)
C  and get back: the adjusted/computed w* and v* values: 
C     REAL*8 WBAR(2,Z$X1), W1(L$,Z$X1), W(L$,Z$X1), V(L$+1,Z$X) in COMMON
C

          CALL VWMASS 


          print *, 'im324,  streamf iterations = ', im324, nt1

cc      if (mod(im10-2,9) .eq. 0.0) then
         do 199 ik=1,Z$X1
ccccccccccc       write(6,299) ik, WBAR(1,ik), WBAR(2,ik), w1(13,ik), w(13,ik)
 299        format(1x, i5, 2x, 1P4D14.4)
 199      continue
cc              print *, ' '
cc      endif


c  now load arrays for use in model (SETDAILY), need to do WALL and VALL separately
c     VALL here is just used for output, it's updated every day in model, so DON'T NEED to store for each
C        time period:     WALL(L$,Z$X1,74), VALL(L$+1,Z$X) are in COMMON as REAL*4
C
           do 600 ik=1,Z$X1
           do 600 ij=1,L$
 600          wall(ij,ik,im10+1) = w(ij,ik)

           do 610 ik=1,Z$X
           do 610 ij=1,L$+1
 610          vall(ij,ik) = v(ij,ik)

C
c  interpolate Kzz to constituent box edges for NEWDIFZ, convert to cm2/sec, XKZZT(L$S,Z$S) is REAL*8, 
C  REAL*4 ZZINK(L$S,Z$S), LATST4(L$S), ZSTR4(Z$S), zzoutk(L$,Z$X1), LAT4(L$), ZALTE(Z$X1),
C    - probably best to interpolate logarithmically
C
           do 701 ik=1,Z$S
           do 701 ij=1,L$S
 701          zzink(ij,ik) = xkzzt(ij,ik)
       
            CALL BINTERP(LATST4, L$S, ZSTR4, Z$S, ZZINK, 
     >                   LAT4, L$, ZALTE, Z$X1, 0, 0, ZZOUTK)
c                                                               KZZGW(L$,Z$X1) in COMMON (REAL*4)
           do 501 ik=1,Z$X1
           do 501 ij=1,L$
 501          kzzgw(ij,ik) = zzoutk(ij,ik)*1.e4


c   convert Kyy to cm2/sec,  KYY(L$+1,Z$X) ,  KYYALL(L$+1,Z$X,74)

           do 502 ik=1,Z$X
           do 502 ij=1,L$+1
 502          kyyall(ij,ik,im10+1) = kyy(ij,ik)*1.d4



c  now load other dynamics arrays NOT USED IN MODEL, into full annual arrays, 
C       TEMPALLS(L$S,Z$S) is for STREAMF/TROPKZZ

           do 507 ik=1,Z$S
           do 507 ij=1,L$S
               ubarall(ij,ik) = ubar1(ij,ik)   
               heatall(ij,ik) = heat1(ij,ik)   
               chiall(ij,ik) = chif(ij,ik)   
               fxtall(ij,ik) = fxtt(ij,ik)         
               fkall(ij,ik) = fk(ij,ik)         
               fkqall(ij,ik) = fkq(ij,ik)         
               fkq2all(ij,ik) = fkq2(ij,ik)         
               frqall(ij,ik) = frq(ij,ik)
               fegwall(ij,ik) = fegw(ij,ik)

               cffall(ij,ik) = cffs(ij,ik)   
               qvall(ij,ik) = qv(ij,ik)         
               kruall(ij,ik) = kru(ij,ik)         
               ehfall(ij,ik) = ehf(ij,ik)   
               qtotall(ij,ik) = qtot(ij,ik)   
               lhall(ij,ik) = lh1(ij,ik)   
               ddhall(ij,ik) = ddh(ij,ik)   
               tempalls(ij,ik) = temp1(ij,ik)   

cc      ! load in REAL*4 array for output below, walla8(L$S,Z$S) is REAL*8, walla(L$S,Z$S) is REAL*4
               walla(ij,ik) = walla8(ij,ik)   
 507    continue
            

C
C  ****  Write out dynamics  arrays  to fort.37, fort.38, fort.57  for  new  CIRCULATION  if requested   ****
C        
C
      IF (LDYNOUT) THEN
                                                 ! just write out dimensions, etc. 1st time through
        if (im324 .eq. 1) then 
ccdyn           WRITE (37) L$, L$S, Z$S, Z$X, itdyn, NMON$
ccdyn           WRITE (37) LAT
ccdyn           WRITE (37) ZALT
ccdyn           WRITE (37) PRESS
ccdyn           WRITE (37) ZSTR
ccdyn           WRITE (37) PRST
ccdyn           WRITE (37) LATST
ccdyn           WRITE (37) ZALTE
ccdyn           WRITE (37) PRESSE
ccdyn           WRITE (37) LATEG

           WRITE (57) L$, L$S, Z$S, Z$X, itdyn, NYR$
           WRITE (57) LAT
           WRITE (57) ZALT
           WRITE (57) PRESS
           WRITE (57) ZSTR
           WRITE (57) PRST
           WRITE (57) LATST
           WRITE (57) ZALTE
           WRITE (57) PRESSE
           WRITE (57) LATEG
         endif


        if (im324 .eq. 1) then 
           WRITE (38) L$, L$S, Z$S, Z$X, itdyn, NYR$
           WRITE (38) LAT
           WRITE (38) ZALT
           WRITE (38) PRESS
           WRITE (38) ZSTR
           WRITE (38) PRST
           WRITE (38) LATST
        endif

           WRITE (38) tempalls
           WRITE (38) walla
           WRITE (38) valla
           WRITE (38) ubarall
           WRITE (38) heatall
           WRITE (38) chiall
           WRITE (38) fxtall
           WRITE (38) fkall
           WRITE (38) fkqall
           WRITE (38) fkq2all
           WRITE (38) frqall
           WRITE (38) fegwall

           WRITE (38) cffall
           WRITE (38) qvall
           WRITE (38) kruall
           WRITE (38) ehfall
           WRITE (38) qtotall
           WRITE (38) lhall
           WRITE (38) ddhall
           WRITE (38) dmom
           WRITE (38) gwh1
           WRITE (38) gwh2
      ENDIF


        RETURN
	END


cdbc  The following subroutine is copied from the 2d model written by
cdbc  Mark Schoeberl. (copied 10/22/91) It has been modified so 
cdbc  that boundaries are not set to zero.   -- TAKEN FROM OLD KYY ROUTINE

      SUBROUTINE DERV4(TYP, XIN, OUT, DEL, N, M, SW)

C      DO A FOURTH OR SECOND ORDER FIRST DERIVATIVE

C
C       IF TYP = 1 THEN SECOND ORDER
C       IF TYP = 2 THEN FOURTH ORDER
C

C      IF SW=1 THEN INNER INDEX
C      IF SW=2 THEN OUTER INDEX


       REAL*8 xin(N,M), out(N,M), del
       REAL*8 ay1, ay2, ay3

       INTEGER N,M
       INTEGER TYP,SW


      GO TO (1,2),TYP
      STOP ' WRONG TYP IN DERV4'

2      AY1=1./(12.*DEL)
       AY2=2./(3.*DEL)
       AY3=1./(2.*DEL)
      GOTO 8
1       AY1=0.
      AY2=1./(2.*DEL)
      AY3=AY2
8      GO TO (10,20),SW
      STOP 'ILLEGAL CALL TO DERV4'      

10      N1=N-1
      N2=N-2

C  - F. Vitt correction put in 5/19/94

      DO 51 IK=1,M
      DO 50 IJ=3,N2

      OUT(IJ,IK)=(-AY1*XIN(IJ+2,IK)+AY2*XIN(IJ+1,IK)-AY2*XIN(IJ-1,IK)
     1  +AY1*XIN(IJ-2,IK))

50      CONTINUE

      OUT(2,IK)=(XIN(3,IK)-XIN(1,IK))*AY3
      OUT(n1,IK)=(XIN(n,IK)-XIN(n2,IK))*AY3
      OUT(1,IK)=(xin(2,ik)-xin(1,ik))/del
      OUT(N,IK)=(xin(n,ik)-xin(n1,ik))/del
51      CONTINUE

      RETURN

20    M1=M-1
      M2=M-2
      DO 61 IJ=1,N
         if (ij .eq. N+1) print *,' IN DERV4: ', ij, n
      DO 60 IK=3,M2
      OUT(IJ,IK)=(-AY1*XIN(IJ,IK+2)+AY2*XIN(IJ,IK+1)-AY2*XIN(IJ,IK-1)
     1+AY1*XIN(IJ,IK-2))
60    CONTINUE
      OUT(IJ,2)=(XIN(IJ,3)-XIN(IJ,1))*AY3
      OUT(IJ,m1)=(XIN(IJ,m)-XIN(IJ,m2))*AY3
      OUT(IJ,1)=(xin(ij,2)-xin(ij,1))/del
      OUT(IJ,M)=(xin(ij,m)-xin(ij,m1))/del
61    CONTINUE

      RETURN
      END

