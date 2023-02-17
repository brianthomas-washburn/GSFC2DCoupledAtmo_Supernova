        SUBROUTINE GWAVE(im10, iiyr, im324) 
C
C    program GWAVE_9BA.F   (EF -- Aug-Sept. 95), re-done to include longitudinal variations June 97
C                              also includes eddy heat flux effects as in Huang and Smith, 1991 (Aug 97)
C
C   Gravity waves 'r us......Program to generate gravity wave model using the parameterizations
C     of Lindzen(1981), Holton(1982), and Holton and Zhu (1984), assuming an isotropic 
C     distribution of gwaves in upper troposphere, see also Garcia & Solomon(1985)
C     THIS EXTENDS to Z$S LEVS (116 km), since 90 km level can be affected by waves breaking above
C    
C  USES  WAVE SPEEDS OF 0, +-10, +-20, +-30, +-40, +-50, W/ efficiency factors, and
C     ZONAL WAVELENGTHS defined as functions of phase speed, also ASSUME k=K. l=0, so
C     cos(angle) = k/K = 1, K=sqrt(l^2 + k^2), as described in HZ84, for simplicity, and
C     as used in other 2D modelling studies, eg, GS85, and Brasseur et al, 1990.
C
C; NOTE: differences w/ GS85, obs.: GS85 has summer transistion hgt of Ubar higher than in winter
C;   which is opposite to obs. We have a sharper summer peak than in winter, consistent w/ obs,
C;   and GS85. Our winter Fx, Kzz, and v are more spread out in altitude than in summer (GOOD). 
C;   You would expect the Fx and v to correspond to the Ubar transistion hgt -- this is so in
C;   summer in obs, but not necessarily in winter since things are spread out more in altitude.
C;  So our computed Fx and v seem to be OK, except that the winter peak is a bit larger in magnitude
C;  than summer, GS85 shows the opposite, and Lindzen shows a bit larger Fx in summer, but
C;  a bit larger Kzz in winter than summer; Radar obs. from Manson et al., 1991 don't show any
C;  significant seasonal differences in v magnitudes, but the altitudes of maxima certainly have
C;  a seasonal dependence. Note finally that tweeking the k (using 300-2000 km) and/or w0 in our
C;   calculations did NOT qualitatively change the seasonal variations of the altitude of largest
C;   Kzz/Fx, so we'll just keep them at 200 km and standard w0.
C;  ALSO, HB90 show no such seasonal variation in the altitudes of max Kzz, Fx, but they do show
C;    the largest Kzz and Fx in the winter hemisphere, but with a much greater difference than we do.
C
C ****  NOTE:  Kzz is limited to  </= 100 m2/sec, in accordance w/ 1/8 day (3 hr) time step used for vertical
C       diffusion in MAIN, Fx also reduced accordingly by modifying w0 and k (??), base efficiency
C       factor, beff set to .025, as before, anything bigger gives a Fx, and w* and v* which seem too 
C       big, even .03 gives Fx of ~250 m/s/day - too big. Also, Kzz w/ .03 seems OK, but maybe too big
C       at winter mid-high lats, as H2O at 80 km is a bit too wet there compared to HALOE v18 ************
C
C     ALL CALCULTIONS DONE IN MKS UNITS
C
C  **** NOTE: THIS IS done in SINGLE PRECISION, since it's too messy to convert to REAL*8,
C       and it doesn't matter anyway since it's just a parameterization, so just 
C       put the final (COMMON BLOCK) arrays in REAL*8   
C     
C   Longitudinal variations NOT INCLUDED in this routine- UBARIN(L$S,Z$S,72), (in COMMON), 
C       longitudinal BC's read in TEMPIN, array GWBC(13,L$S,72) defined in COMMON
C
C   NEW - REVAMPED computations for Zbr, sorting, Diffusion, and Diff-spor., should now
C         be robust, (same algorithm as HZ84), and NO ITERATING!      (EF, April 98)


       include "com2d.h"


C REAL*8 XN2(L$S,Z$S), UBAR1(L$S,Z$S), DUDZ(L$S,Z$S), XKZZT(L$S,Z$S), FXTT(L$S,Z$S)  ALL IN COMMON
C   zstr(Z$S), prst(Z$S), rhost(Z$S), latst(L$S), cose(L$S), tana(L$S), ff0(L$S)  ALL IN COMMON



C  ITN = Number of iterations to compute Zbr and DIffusion;; 

       INTEGER inw, itn, ik0

       PARAMETER(inw=11, itn=1)

       INTEGER icr(inw), iuns(inw+1), iccbr(inw+1)
       INTEGER iwbr(inw+1), iwbrf(inw+1), ibrlev(inw+1)
       INTEGER ibr(inw+1,inw,itn+1), isort(inw+1,inw,itn+1)

       REAL*8 UUU(L$S,Z$S), UUZ(L$S,Z$S), DUDZG(13,L$S,Z$S)
       REAL*8 UG(13,L$S,Z$S)

cccc       REAL FXTLON(13,L$S,Z$S,72), KZZLON(13,L$S,Z$S,72)
       REAL FXTLON(13,L$S,Z$S), KZZLON(13,L$S,Z$S) 

       REAL nc1(58), nci(Z$S), nc(Z$S)
       REAL ang1(1), skk0(inw), eff1(inw), kk1(inw)
       REAL cc0(inw), kk0(inw), wb0(inw), w0(inw), effg(inw), eff0(inw)
       REAL ucc(inw,Z$S), usqg(inw), succ(inw,Z$S), zxy(inw,Z$S,itn+1)

C   difft = total diffusion from all waves
C   diffw = total diffusion by each wave sorted by brkg level, ie, has diffusion in indicies in the 
C           order of breaking, and DOES NOT zero-out above the Zbr of the next wave, whereas 
C           DIFFA has total diffusion for each wave via the ICC indicies, and does zero-out
C   diffa = total diffusion by wave (ICC), used only at end where it's summed to get Kzz
C   fxsp is exponential decay below zb due to sporadic brkng  ;; fx43 is from eq 43 in HZ84 at zb, NOT USED


       REAL difft(Z$S,itn), diff(inw,Z$S,itn), diffw(inw+1,Z$S,itn)   
       REAL nn(Z$S), dmm(Z$S), diffa(inw,Z$S,itn) 

       REAL fxg(inw,Z$S), fxsp(inw,Z$S), fx43(inw), lamiw(inw,Z$S)
       REAL lamrw(inw,Z$S), lamr(Z$S), lami(Z$S), fxt(Z$S), kzz(Z$S)

       REAL kzzl(13,L$S,Z$S), fxtl(13,L$S,Z$S)
       REAL kzzdi(L$S,Z$S), fxtdi(L$S,Z$S)
       REAL gwh1a(inw,Z$S), gwh1b(13,L$S,Z$S), gwh2b(13,L$S,Z$S)
       REAL*8 kz1(L$S,Z$S), kz2(L$S,Z$S), fx1(L$S,Z$S), fx2(L$S,Z$S)
       REAL*8 kz0(L$S,Z$S), fx0(L$S,Z$S)
       REAL*8 dm1(12,Z$S), dmz(12,Z$S), dm2(L$S,Z$S), dmz2(L$S,Z$S)

       DATA ang1/0.0/

       DATA nc1/.12, .16, .2, .24, .28, .31, .34, .4, .5, .63, .76, .87,
     c  .99, 1.1, 1.3, 1.5, 1.7, 2.0, 2.2, 2.4, 2.6, 2.8, 3., 3.3, 3.6, 
     c  3.8, 3.9, 4., 4.1, 4.3, 4.5, 4.7, 4.9, 5., 5.1, 5.3, 5.5, 5.7, 
     c  5.75, 5.72, 5.7, 5.65, 5.47, 5., 4.3, 3.3, 2.3, 1.6, 1.6, 2.15,
     c  2.9, 3.75, 5.05, 7.25, 10., 12.5, 14.3, 15.7/

C
C  DEFINE  WAVE  PARAMETERS  HERE  (HZ84), and DEFINE k as a function of c, for tropical gwaves, see G92
C    ; phase speeds, and zonal wavelengths, Lx, defined to be 200 km for +-40 m/s, and 800 km for +,-50m/s
C
C  Lx of 200-500 km, seems best, sometimes 200 km is better than 500 km, and vice-versa,
C      probably just see how residual circ. and constituent distr. are.
C      KEEP AT 200 km for now, w/ exponential dependence of w0 on c, and sin(lat), as in GS85
C  Also, define base efficiency factor, beff

       DATA cc0/-50, -40, -30, -20, -10, 0., 10, 20, 30, 40, 50/   !, ik0/4/
       DATA kk1/800., 9*200., 800./
       DATA eff1/11*1./, beff/.025/, pran/1./              
C                                                                      7
C               -50, -40,  -30,  -20,  -10,  0,  10,  20,  30,  40, 50 2
       DATA wb0/.02,  1.,   2. ,  4,    8.,  9,  10., 10,  7,   4,  .02/

CCCold DATA wb0/.02, 2.03, 4.43, 7.71, 10.8, 12.03, 10.8, 7.71, 4.43, 2.03, .02/


         SAVE


C    ik0 = bottom boundary for calcs, set to around 8 km (~325 mb), level 4 for original grid
C                                                        (need to add 1 since we're now starting at 0 km)
       ik0 = INT(8./116.*Z$S)+1

C
C  define indicies for diurnal tide computations at upper levels
C

       i60 = INT(60./116.*Z$S)+1
       i80 = INT(80./116.*Z$S)+1
       i86 = INT(86./116.*Z$S)+1
       i100 = INT(100./116.*Z$S)+1

       i25 = INT(25./116.*Z$S)+1
       i35 = INT(35./116.*Z$S)+1
       i45 = INT(45./116.*Z$S)+1


       do 201 ii=1,inw
          kk0(ii) = 2.*xpi/(kk1(ii)*1000.)
            if (kk0(ii)*cos(qp1*ang1(1)) .ge. 0.) skk0(ii) = 1.
            if (kk0(ii)*cos(qp1*ang1(1)) .lt. 0.) skk0(ii) = -1.
 201      eff0(ii) = eff1(ii)*beff


C  newtonian cooling from HZ84, convert to proper units, nc1(58) on zz58(58) grid, 
C     first interpolate (logarithmically) to proper grid,   NCI(Z$S),  ZSTR4(Z$S) is REAL*4
C
       CALL LINTERP(ZZ58, 58, NC1, ZSTR4, Z$S, 1, NCI)  

       do 103 ik=1,Z$S
          nc(ik) = nci(ik)*1.e-6
 103      dmm(ik) = 200.*exp((zstr(ik) - 110.)/7.)
c                                                  molecular diffusion (dmm) used in zbr and lami calculation
C                                                  initialize kzzt, fxtt at begining of each 10-day period

       do 312 ik=1,Z$S
       do 312 ij=1,L$S
       do 312 ix=1,13
          kzzl(ix,ij,ik) = 0.0
 312      fxtl(ix,ij,ik) = 0.0

       do 313 ik=1,Z$S
       do 313 ij=1,L$S
          xkzzt(ij,ik) = 0.0
 313      fxtt(ij,ik) = 0.0


C define du/dz for each lon/lat; UBARIN(L$S,Z$S,72); REAL*8 UUU, UUZ(L$S,Z$S), DUDZG(13,L$S,Z$S)
C    also load in UBAR for this particular im10 period into UG(13,L$S,Z$S), IX=13 = zonal mean
C
C    for 9BA, just use ZONAL MEAN UBAR for ALL LONGITUDES in UG
C
       do 315 ix=1,13
           do 317 ik=1,Z$S
           do 317 ij=1,L$S 
              uuu(ij,ik)  =  ubarin(ij,ik,im10)        !! ubar10(ix,ij,ik,im10,iiyr)
 317          ug(ix,ij,ik) = ubarin(ij,ik,im10)        !! ubar10(ix,ij,ik,im10,iiyr)

              CALL DERV4(1, uuu, uuz, delz, L$S, Z$S, 2)

           do 318 ik=1,Z$S
           do 318 ij=1,L$S 
 318          dudzg(ix,ij,ik) = uuz(ij,ik)
 315   CONTINUE


c  SWITCH for longitudinally varying Gwaves, igwlon=1  (igwlon=0=zonal avg)
C
C  NOTE: for interannual variability, we currently have yearly variations for UBAR, 
C    but the longitudinal variations are CLIMATOLOGICAL (see /misc/kah02/newmod/ukmok99/kyy99k_year.pro)
C    so just compute the ZONAL MEAN GWAVE FORCING here for ALL QBO INTERANNUAL VARIABILITY RUNS
C 
C    but for doing interannual variability outside of tropics, it's best to put in the yearly 
C    longitudinal zonal wind variations and compute the longitudinal GWAVE effects (as per Siskind Proposal)
C
C   if  IGWLON=2 it does longitudinal variations of u and the zonal mean u loaded into kzzlon, fxtlon
C      but does NOT zonally avg the values (ix=13 in these is the values based on the zonal mean U)
c     also, IGWLON determines loop bounds if zonal mean or longitudinal variations

ccc       igwlon = 1
       igwlon = 0

          ix1 = 1
          ix2 = 12

       if (igwlon .eq. 0) then 
          ix1 = 13
          ix2 = 13
       endif

       if (igwlon .eq. 2) ix2 = 13


      DO 5000 IJT=1,L$S 

C           print *, 'ijt, latst = ', ijt, latst(ijt)

           do 202 ik=1,Z$S
 202          nn(ik) = sqrt(abs(xn2(ijt,ik)))


        DO 7000 IXT=ix1,ix2 


        do 104 ii=1,inw
           effg(ii) = eff0(ii)
 104       w0(ii) = wb0(ii)*1.e-3*3.0*gwbc(ixt,ijt,im10)

cccc      gwbc(ixt,ijt,im10) = exp(-cos(latst(ijt)*qp1*2)**6)  !original BCs, simple latitude-dependence only
cccc 104       w0(ii) = wb0(ii)*1.e-3*3.0*exp(-cos(latst(ijt)*qp1*2)**6)

cccccc   *******  For longitude variations, use the 104 w0(... statement here with  ! gwbc(13,L$S,72)
CCCCCC            read in from D. Wu's MLS data w/ proper scaling, this applies the cc dependence in wb0
CCCCCC            ALSO change special CC = +-50 m/s tropical waves below??????


C    re-define efficiencies to be lower for NH 30-90N winter (Nov-March), as described in Shine (1989)
C    *****   we don't need this out with a 3D wind field  so that effg = eff0 (.10) above
CC           for all phase speeds for all latitudes, longitudes, and time periods  
CC           Note: the effg(inw) array is only needed for the NH winter  **********
CC 
c       Try taking thsi out even with zonal mean U, since Ubar will be weaker in NH winter anyway,
C         and this should be reflected in the Kzz, Fxt
C
CCCCCC      if (igwlon .eq. 0) then 
CCCCCC         if (im10 .ge. 31 .or. im10 .le. 9) then
CCCCCC           if (ijt .ge. 25) then
CCCCC                do 105 ii=1,inw
CCCCC 105            effg(ii) = eff0(ii) - eff0(ii)*.5*
CCCCC     c                exp(-abs((latst(ijt)-60.)/10.))*cos(qp1*10.*im10)
CCCCC            endif
CCCCC          endif
CCCCC      endif

C
C  redefine w0 for c=+-50 m/s from G92 (eq 23) for 30S-30N, reduced to 0.2e-3 since kzz,Fx were
C   too large otherwise; Note: even at uw0=.2e-3, these waves still have a significant impact, 
C   but any larger of a uw0>.2e-3 would give too large of a Kzz, Fx
C   We won't mess with observational/longitudinal BCs for these tropical waves, just use, UG(13,L$S,Z$S)
C
      if (LATST(ijt) .ge. -30. .and. LATST(ijt) .le. 30.) then
        uw0 = 0.2e-3*exp(-(LATST(ijt)*qp1)**2/.02)
        w0(1) = sqrt(uw0*2.*kk0(1)*abs(ug(ixt,ijt,ik0)-cc0(1))/nn(ik0))
       w0(11)= sqrt(uw0*2.*kk0(11)*abs(ug(ixt,ijt,ik0)-cc0(11))/nn(ik0))
      endif 
              
       do 106 ii=1,inw
 106      icr(ii) = Z$S+1
                                                        
       do 108 ik=1,itn+1
       do 108 ij=1,inw  
       do 108 ii=1,inw+1
          isort(ii,ij,ik) = 0
 108      ibr(ii,ij,ik) = 1

       do 119 itt=1,itn+1                                   
       do 119 ik=1,Z$S
       do 119 ii=1,inw
 119      zxy(ii,ik,itt) = 0.

c  difft is the total diffusion, backgrd = .1 m2/sec everywhere, tropospheric Kzz in TROPKZZ routine
C     also initialize diff, difft, diffa, diffw arrays
C
       do 109 itt=1,itn
       do 110 ik=1,Z$S
 110      difft(ik,itt) = .1 
CEF       do 111 ik=4,6
CEF 111      difft(ik,itt) = 1.*exp(5.-ik) 
 109   continue

       do 801 itt=1,itn
       do 801 ik=1,Z$S
       do 801 ii=1,inw
           diff(ii,ik,itt) = 0.0
 801       diffa(ii,ik,itt) = 0.0

       do 802 itt=1,itn
       do 802 ik=1,Z$S
       do 802 ii=1,inw+1
 802       diffw(ii,ik,itt) = 0.0

C
C ********  Now loop through all waves to compute usqg, and icr  ***********************

       DO 2000 ICC=1,11
                                                                                 ! UG(13,L$S,Z$S) 
          do 112 ik=1,Z$S
             ucc(icc,ik) = ug(ixt,ijt,ik)*cos(qp1*ang1(1)) - cc0(icc)
               if (ucc(icc,ik) .lt. 0.0) succ(icc,ik) = -1.
               if (ucc(icc,ik) .ge. 0.0) succ(icc,ik) = 1.
 112      continue
             usqg(icc) = ((w0(icc)*nn(ik0)/kk0(icc))**.666667)/
     c                (abs(ucc(icc,ik0))**.333333)

C  now loop through in height to find critical level, and breaking height of wave
C  Note that for large ucc, the breaking level can get very high (> 90 km), ICR always is from ik0+2->Z$S+1
C
        do 115 ik = ik0+2, Z$S

           if (abs(ucc(icc,ik)) .lt. 1.0 .or. 
     c                          succ(icc,ik) .ne. succ(icc,ik-1)) then
                icr(icc) = ik
CC                               print, 'icr, zz = ', icr(icc), zstr(icr(icc))
                GO TO 2000
           endif
 115     continue

 2000  CONTINUE

C
C   Now loop through all waves, find breaking heights by eq (34), and compute
C   diffusion by eq (33);  Note that sporadic breaking from components with higher Zb
C   is not included in the sum1 integral for the Zb computation (see discussion of 
C   eq 34 in HZ84), HOWEVER, MOLECULAR DIFFUSION IS NOW INCLUDED IN THE ZBR CALCULATION.
C   Diffusion not reduced by effg here either, only at very end in Kzz calc.
C   When searching for the breaking levels, need to start at ik0+1 (for zxy), 
C       but when integrating, start at z0 (ik0), use trapazoidal rule for integrals
c
c   PROBLEM WITH ITERATING - GET OSCILLATORY BEHAVIOR, I.E., USING PREVIOUSLY COMPUTED DIFFT
C      TO THEN RE-COMPUTE ZBR - ONLY WORKS IF DIFFUSION IS MUCH SMALLER (E.G., USE INTERMITANCY),
C      BUT THIS IS NO GOOD. - JUST GO THROUGH ONE TIME ONLY - ITN=1, include Dsp for each wave
C      ONLY after Zbr and D for the wave have been determined - don't make Zbr and D and Dsp dependent
C        on one another - it'll oscillate also.
C
C    DON"T NEED TO ITERATE!!! ONE TIME IS FINE, THE ANSWER DOESN'T CHANGE WITH MORE ITERATIONS!!!!!

       DO 3000 ITER = 1,ITN

          do 138 iq = 1, inw+1
            isx = 0
            iuns(iq) = Z$S+1
            iccbr(iq) = 0
 138        iwbr(iq) = 1


          DO 150 IWV = 1,inw

c  now loop through ALL WAVES, and compute the Zbr, with difft = total diffusion from lower waves
c        initialize ibr and isort, also check for any waves that previously did not break, where isort=0
C        NOTE:   ICR always is from ik0+2-> Z$S+1

       do 232 icc=1,inw
          ibr(iwv,icc,iter) = 1
 232      isort(iwv,icc,iter) = 0

         
        DO 160 ICC=1,inw
c                                                              set up so it correctly does 1st wave
          if (iwv .eq. 1) isortp = 1
          if (iwv .ge. 2) isortp = isort(iwv-1,icc,iter)

          IF (isortp .ge. IWV  .or.  isortp .eq. 0) THEN

             do 155 ik = ik0+1,icr(icc)-1
                     sum=0.0                                        
                        do 153 ii = ik0,ik-1
                          xa2 = (nn(ii)*nc(ii)/(ucc(icc,ii)**2) + 
     c                          nn(ii+1)*nc(ii+1)/(ucc(icc,ii+1)**2))/2. 
                          sum = sum + xa2*delz
 153                    continue

                sum1 = 0.0
        do 154 ii = ik0,ik-1
         xa3 = ((nn(ii)**3)*(difft(ii,iter)+dmm(ii))/(ucc(icc,ii)**4) +
     c   (nn(ii+1)**3)*(difft(ii+1,iter)+dmm(ii))/(ucc(icc,ii+1)**4))/2.
                sum1 = sum1 + xa3*delz
 154   continue

         zxy(icc,ik,iter) = zstr(ik0)*1000. + hh/kk0(icc)*sum +
     c      2.*hh/kk0(icc)*sum1 + 3.*hh*alog(abs(ucc(icc,ik)/usqg(icc)))
     c      - zstr(ik)*1000. 
 155     continue

c                                                          breaking level should be at or above ik0+2
         do 157 ik = ik0+2,icr(icc)-1
            if (zxy(icc,ik,iter) .le. 0.0 .and. 
     c                               zxy(icc,ik-1,iter) .gt. 0.0) then
                                                                              ! updated zbr
               ibr(iwv,icc,iter) = ik
C                                                  ;; print, 'ibr, zz = ', ibr(iwv,icc,1), zz(ibr(iwv,icc,1))
c                                                   ; If Zb < 30 km, continue searching for a higher Zb
               if (zstr(ik) .gt. 30.) go to 160
            endif
 157     continue

         ENDIF

 160     CONTINUE

c
C   Now sort the breaking levels for the REMAINING waves (all waves if on first go-around of iwv)
c       Note: just loop through as usual, for waves that have already broken, ibr=1, so it won't be sorted
C                                             ! ikst actually can (should) be ik0+2, but it shouldn't matter
          if (iwv .eq. 1) ikst = ik0+1  
          if (ikst .ge. Z$S) ikst = Z$S 

          iii = iwv

       do 165 ik = ikst,Z$S
 
           do 334 ixyz=1,12
 334          ibrlev(ixyz) = 0

          is = 0
          do 166 icc = 1,11
             if (ibr(iwv,icc,iter) .eq. ik) then 
                isort(iwv,icc,iter) = iii
                if (iii .eq. iwv) ikst = ik + 1
                                              !for the next iwv go-around, start sorting at the next level up
                is = is + 1
                ibrlev(is) = icc
             endif
 166        continue
c                                      ! IS counts how many waves have Zbr at the current ik, store in ibrlev
            if (is .ge. 12) then
               print *, 'IN GWAVE, IS ge 12, STOP!!!'
               STOP
            endif
c                                           ! this part is for 2 or more waves with same Zbr - see note below
C
        if (is .ge. 2  .and.  iii .eq. iwv) then
           CALL SORTOUT(ucc, Z$S, is, ibrlev, ik, index)
c                                                        ! index = index (icc) of wave with the largest u-c
           ijcount = 1
           do 266 ipp=1,is
             icc = ibrlev(ipp)
              if (icc .ne. index) then
                  isort(iwv,icc,iter) = iii + ijcount
                  ijcount = ijcount + 1
              endif
 266       continue
        endif          

          iii = iii + is
 165   continue

c Note: if two or more waves have the same Zbr and isort, keep the one with the greatest abs(ubar-c), since
C       this will give the greater diffusion, and presumably would cause the wave with the
C       smaller ubar-c to become unsaturated, then the smaller ubar-c becomes isort+1
C        -- only worry about this for the current IWV, since upper IWV's will get modified anyway
C       Note, the greatest U-C will depend on which waves have Zbr at the level in question

c
c   Now compute diffusion D for current IWV

       DO 260 ICC=1,11
           IF (isort(iwv,icc,iter) .eq. IWV) THEN

              iwbr(iwv) = ibr(iwv,icc,iter) 
              iccbr(iwv) = icc
                                 ! at end of IWV loop, ISX is the total number of waves that break, as before
              isx = iwv

C set difft=0 for wv1 at & above wv2 brk lev, since wv1 is no longer saturated above Zb2, but set to 0 + Mdif
C     also do for diffa, ie just zero these for each IWV
C     I dont think zeroing diffw does anything here, since the diffusion for iwv hasn't been computed yet,
C     but it does matter for diffa, since we loop through ALL waves (and of course, difft is affected also)

                   do 158 ik = iwbr(iwv),Z$S
                       diffw(iwv,ik,iter) = 0.0  
 158                   difft(ik,iter) = 0.0
   
                   do 159 ik = iwbr(iwv),Z$S
                   do 159 icc1 = 1,11
 159                  diffa(icc1,ik,iter) = 0.0 

  
c Note: by definition, iwbr(iwv) LE icr(icc)-1  ALWAYS!!    (see 157 loop above)         ! DUDZG(13,L$S,Z$S)
                                                              
         do 161 ik = iwbr(iwv), icr(icc)-1
            diff(icc,ik,iter) = kk0(icc)/2.*(ucc(icc,ik)**2)/
     c          (nn(ik)**3)*((ucc(icc,ik)**2)/hh - 3.*dudzg(ixt,ijt,ik)*
     c          cos(qp1*ang1(1))*ucc(icc,ik) - nc(ik)*nn(ik)/kk0(icc))

C  If diff goes negative, means that Lami<0, so assume wave becomes unsaturated below Zcr and can grow
         
             if (diff(icc,ik,iter) .lt. 0.0) then
                  diff(icc,ik,iter) = 0.0
                  iuns(iwv) = ik
                  go to 370
             endif
c                                                        for iwv=1, just adding to the .1 backgrnd for difft
C
         diffw(iwv,ik,iter) = diffw(iwv,ik,iter) + diff(icc,ik,iter)
             difft(ik,iter) = difft(ik,iter)     + diff(icc,ik,iter)
         diffa(icc,ik,iter) = diffa(icc,ik,iter) + diff(icc,ik,iter)
 161        CONTINUE

C   diffw = total diffusion by each wave sorted by brkg level, ie, has diffusion in indicies in the 
C           order of breaking, and DOES NOT zero-out above the Zbr of the next wave, whereas 
C           DIFFA has total diffusion for each wave via the ICC indicies, and does zero-out

C   add on here for sporadic breaking below zb(j) level, down to zb(j-1) level, also smooths discontinuities 
c    -- this is the only place where diffw array is used, and in sporadic computing of Fx (near 187 loop)

 370    if (iwv .eq. 1  .and.  ik0+1 .le. iwbr(1)-1) then
             do 142 ik = ik0+1, iwbr(1)-1      
               diff(icc,ik,iter) = diff(icc,iwbr(1),iter)*
     c                              exp((zstr(ik)-zstr(iwbr(1)))/7.)

              difft(ik,iter) = difft(ik,iter) + diff(icc,ik,iter)
              diffw(1,ik,iter) = diffw(1,ik,iter) + diff(icc,ik,iter)

             diffa(icc,ik,iter) = diffa(icc,ik,iter) + diff(icc,ik,iter)
 142         continue
        endif
C                 note, the limits here should be rigorous now, such that iwbr(iwv-1) LE iwbr(iwv)-1 ALWAYS
C   ! diffw(12,Z$S,5)  ,    also check that diff(zb-curr) of current wave is > diff(zb-curr) of previous wave


      if (iwv .ge. 2) then
        diffwd = diffw(iwv,iwbr(iwv),iter) - diffw(iwv-1,iwbr(iwv),iter)
            if (diffwd .le. 0.) diffwd = diffw(iwv,iwbr(iwv),iter) 

           do 162 ik = iwbr(iwv-1), iwbr(iwv)-1
              diffsp = diffwd*exp((zstr(ik)-zstr(iwbr(iwv)))/7.)
  
              difft(ik,iter)     = difft(ik,iter) + diffsp
              diffa(icc,ik,iter) = diffa(icc,ik,iter) + diffsp
              diffw(iwv,ik,iter) = diffw(iwv,ik,iter) + diffsp
 162        continue
      endif


        ENDIF 

 260   CONTINUE
c                                                         260 ends the diffusion calculation

 150  CONTINUE
c                                                         150 ends the IWV loop

 3000  CONTINUE
c                                                3000 ends the Iteration loop - not really used, since ITN=1
c 
C  Note, on final loop through of IWV (150 loop):  
C       iwbr(1-12) = zbr of the waves in order of breaking, 1->ISX   (ISX = total # of waves that break)
C              iccbr(1-12) = icc indicies of the waves in order of breaking, 1->ISX
C
C
C  *****************   NOW COMPUTE MOMENTUM FLUX CONVERGENCE  *********************************
C
C   keep lami always positive below Zb for wave damping by diffusion and Newtonian cooling,
C   otherwise the exponential will go to infinity if lami is negative, as expected
C   lamr must be of opposite sign of little k for upward energy propagation
C
C     Need to initialize Fx arrays here first

       do 805 ik=1,Z$S
       do 805 ii=1,inw
           gwh1a(ii,ik) = 0.0   
           fxg(ii,ik) = 0.0
           fxsp(ii,ik) = 0.0
           fx43(ii) = 0.0
           lamiw(ii,ik) = 0.0
 805       lamrw(ii,ik) = 0.0


       DO 4000 ICC = 1,INW
C
C                    !iwv here = sorting number of wave for icc index, iwbr(12) = Zbr in order of breaking
           iwv = 0
           do 335 imm = 1,11
 335          if (iccbr(imm) .eq. icc) iwv = imm

         IF (IWV .GE. 1) THEN 
              
              do 171 ik=1,Z$S
                 lamr(ik) = 0.0
 171             lami(ik) = 0.0

              do 172 ik = ik0,Z$S
 172             lamr(ik) = abs(nn(ik)/ucc(icc,ik))*(-skk0(icc))
c                                                                               ;; (eq 16) for lami
              do 173 ik = ik0,iwbr(1)-1
 173             lami(ik) = nn(ik)*nc(ik)/(2.*kk0(icc)*(ucc(icc,ik)**2))
c                                                                               ;; (eq 29) for lami
              do 174 ik = iwbr(1),icr(icc)-1
 174             lami(ik) = abs(lamr(ik)/(kk0(icc)*ucc(icc,ik))*
     c              ((difft(ik,itn)+dmm(ik))*(lamr(ik)**2) + nc(ik)/2.))

C   compute lami in saturated regions -eq (32), just write over what was computed above by eq (29)

            if (iwv .lt. isx) then
                ib2 = iwbr(iwv+1)                              
                if (icr(icc) .lt. iwbr(iwv+1)) ib2 = icr(icc)           
            endif                                    
                 if (iwv .eq. isx) ib2 = icr(icc)
c                                                           if wave goes unsat, compute lami by (29), w/ D=0
            if (iuns(iwv) .lt. ib2) ib2 = iuns(iwv)

C
C  NOTE:  PREVIOUS PROBLEMS (both in the diffusion calc. and here) due to loop bounds with ib1 and ib2
C         going backwards (e.g., do 176 ik=50,40) which the compiler couldn't handle, 
C         so it went haywire. But, the loops here now with ib1, ib2, ib3 should all now be robust, 
C         since ib1 and ib2 are the lowest of the ICR(ICC), IWBR(IWV+1), or IUNS(IWV), and these 
C         should all be above IWBR(IWV). Problems before were with iuns, which should be fixed now
C         but any problems in the future could point to IUNS. And BEWARE of computing loop bounds
C         interactively within the routine like this - BE VERY CAREFUL!!!  --   DUDZG(13,L$S,Z$S)
C         

         do 176 ik = iwbr(iwv),ib2-1
             lami(ik) = 1./(2.*hh) 
     c           - 3./2.*dudzg(ixt,ijt,ik)*cos(qp1*ang1(1))/ucc(icc,ik)
 176     continue

C                                           Fx and diffusion start at ik0+1, they = 0 at ik0 (z0)
       do 180 ik = ik0+1, iwbr(iwv)-1
          sum2 = 0.0
             do 181 ii = ik0,ik-1
 181            sum2 = sum2 + (lami(ii) + lami(ii+1))/2.*delz

             if (sum2 .gt. 25.) sum2 = 500.  

C           ;; if integral gets high (strong damping), this means that Fx=>0, so set sum2 to a 
C            ;;   large number so the exp[..] = 0.0  (below) ;; compute Fx below breaking level, eq(40) 
C 
         fxg(icc,ik) = -kk0(icc)*cos(qp1*ang1(1))*
     c        ((w0(icc)*nn(ik0)/kk0(icc))**2)*lami(ik)*abs(ucc(icc,ik))/ 
     c         (nn(ik)*abs(ucc(icc,ik0))*ucc(icc,ik))*
     c          exp((zstr(ik)-zstr(ik0))*1000./hh - 2.*sum2)

 180   continue
C                                                                 compute Fx in saturated regions by eq (44)
           if (iwv .lt. isx) then
               ib1 = iwbr(iwv+1)
                 if (icr(icc) .lt. iwbr(iwv+1)) ib1 = icr(icc)
           endif
C                                                               ;; if wave goes unsat, DO NOT Compute Fx here
           if (iwv .eq. isx) ib1 = icr(icc)
           if (iuns(iwv) .lt. ib1) ib1 = iuns(iwv) 
                                                                             ! DUDZG(13,L$S,Z$S)
          do 182 ik = iwbr(iwv), ib1-1
            fxg(icc,ik) = -cos(qp1*ang1(1))*kk0(icc)/(2.*nn(ik))
     c       *(ucc(icc,ik)**3)*
     c       (1./hh - 3.*dudzg(ixt,ijt,ik)*cos(qp1*ang1(1))/ucc(icc,ik))
 182       continue

C  also compute heat flux convergence in breaking regions only from eqn. 9 in Huang and Smith, 1991
C    term 1 is here, term 2 is dependent only on the total diffusion, below see also Schoeberl et al., 1983.
C                                          ! gwh1a(inw,Z$S), this is in terms of HEATING (i.e., negative)
          do 1129 ik = iwbr(iwv), ib1-1
 1129        gwh1a(icc,ik) = -717./(2.*rr1*1004.)*nc(ik)*ucc(icc,ik)**2

C                                                     use eq 40 in HZ84 to compute Fx for zb(j+1) < z < zc
           if (iwv .lt. isx) then
             ib3 = iwbr(iwv+1)
             if (iuns(iwv) .lt. iwbr(iwv+1)) ib3 = iuns(iwv) 

                if (icr(icc) .gt. ib3) then
c                                                            this loop is OK because of condition icr > ib3
                  do 185 ik = ib3, icr(icc)-1
                     sum2 = 0.0
                       do 186 ii = ik0,ik-1
 186                     sum2 = sum2 + (lami(ii) + lami(ii+1))/2.*delz

                         if (sum2 .gt. 25.) sum2 = 500. 
c                                                        ; if sum2 large, need to set it so exp below is 0.0
             fxg(icc,ik) = -kk0(icc)*cos(qp1*ang1(1))*
     c               ((w0(icc)*nn(ik0)/kk0(icc))**2)*lami(ik)*
     c          abs(ucc(icc,ik))/(nn(ik)*abs(ucc(icc,ik0))*ucc(icc,ik))
     c               *exp((zstr(ik)-zstr(ik0))*1000./hh - 2.*sum2)

 185                continue

                endif
           endif

C  compute Fx due to sporadic breaking below zb, fx43 at zb in which the exponential is 0.0 
C     Note: using lami at Zb from eq (32) gives fx43 at Zb identical to Fx from eq (44), and lami
C           from (29) gave strange results, so just specify sporadic breaking below Zb for the given
C           wave based on the ratio of the diffusion as in eq (46) and specified below. So by adding
C           fxsp+fx to get total drag, you are just adding the sporadic breaking for each wave to
C           the total ONE time only, which seems more realistic ;;  NOTE: FX43 NOT USED NOW

         ijbr = iwbr(iwv)                            
C                                     Fx decay proportional to diffusive decay,  seems to work OK
         r4 = 0.0 
           if (diffw(iwv,ijbr,itn) .gt. 0.0) then
               if (iwv .eq. 1) dbw = 0.0
               if (iwv .ge. 2) dbw = diffw(iwv-1,ijbr,itn)

                   r4 = fxg(icc,ijbr)*dbw/diffw(iwv,ijbr,itn) 
                   if (dbw .ge. diffw(iwv,ijbr,itn)) r4 = 0.
           endif
c                           this loop is OK, since iwbr(iwv) GE ik0+2  ALWAYS, so iwbr(iwv)-1 GE ik0+1

               do 187 ik = ik0+1, iwbr(iwv)-1
                  fxsp(icc,ik) = (fxg(icc,ijbr) - r4)
     c                            *exp((zstr(ik)-zstr(ijbr))*1000./hh)
 187           continue


                do 191 ik1 = 1,Z$S
                   lamrw(icc,ik1) = lamr(ik1)
 191               lamiw(icc,ik1) = lami(ik1)


         ENDIF

 4000  CONTINUE

C  Sum up Fx for all wave components to obtain the total drag  ;;   fxg, fxsp = fltarr(inw,Z$S)
C     Sum up diffusion in the same way,  ADD IN  EFF FACTOR HERE for BOTH
C              factor in Prandtl number here, so that Kzz = Dh (diffusion for heat and constituents)

        do 195 ik = 1,Z$S
          fxt(ik) = 0.0
           sum4 = 0.0
             do 196 ICC=1,INW
 196            sum4 = sum4 + (fxg(icc,ik) + fxsp(icc,ik))*effg(icc)
 195    fxt(ik) = sum4
C                     

        do 205 ik = 1,Z$S
          kzz(ik) = 0.0
           sum5 = 0.0
              do 206 ICC=1,INW
 206             sum5 = sum5 + diffa(icc,ik,itn)*effg(icc)/pran
 205    kzz(ik) = sum5   
                                                                       
C                                        kzz, kzzl in m2/sec  ;;;    fxt, fxtl in m/sec2            
           do 210 ik=1,Z$S
              kzzl(ixt,ijt,ik) = 0.0
              fxtl(ixt,ijt,ik) = 0.0

              kzzl(ixt,ijt,ik) = kzz(ik)            
 210          fxtl(ixt,ijt,ik) = fxt(ik)            


C  also sum wave heat flux convergence for term 1 which is wave dependent, gwh1a(inw,Z$S), gwh1b(13,L$S,Z$S)
C     include efficiency factor since term is caused by wave breaking
C
          do 1115 ik = 1,Z$S
             gwh1b(ixt,ijt,ik) = 0.0
                sum6 = 0.0
                do 1116 ICC=1,INW
 1116              sum6 = sum6 + gwh1a(icc,ik)*effg(icc)
 1115     gwh1b(ixt,ijt,ik) = sum6   


 7000   CONTINUE
 5000   CONTINUE


C  Now compute zonal mean of kzzl and fxtl,  kzzl, fxtl(13,L$S,Z$S), if igwlon=1  ;  

      if (igwlon .eq. 1) then 

        do 900 ik = 1,Z$S
        do 900 ij = 1,L$S
           sumz = 0.0
           do 902 ix=1,12  
 902           sumz = sumz + kzzl(ix,ij,ik)
           kzzl(13,ij,ik) = sumz/12.

           sumf = 0.0
           do 903 ix=1,12  
 903           sumf = sumf + fxtl(ix,ij,ik)
           fxtl(13,ij,ik) = sumf/12.

 900    CONTINUE
      ENDIF

C
C  NOW SPECIFY Fx and DIFFusion from breaking diurnal tide in upper mesosphere
C  adapted from Lieberman and Hays, 1994, combination of estimates from HRDI data
C   for (1,1) mode only, and estimates from Forbes (1982), including (1,1) and (1,-2) modes,
C   with contributions to both vertical momentum flux convergence (for above 90 km),
C      and horiz mom flux conv (important < 90 km), just compute for 45S-45N

        do 215 ik=1,Z$S 
        do 215 ij=1,L$S
           fxtdi(ij,ik) = 0.0
 215       kzzdi(ij,ik) = 0.0

C semi-annual variation in diurnal tide; for im=1,12 do dta(im-1)=-15.+10.*cos(im*60.*qp1) for monthly values
C   NOTICE dta (and fxtdi) are PUT INTO m/sec2  ;;;   and  dtak is in m2/sec for diffusion
C    note that the largest kzzdi is ~24 m2/sec for model at 90 km, eq. for March, so just keep as is
C    this has been adjusted for 72 time periods per year (5-day avgs)
C
          dta = (-15. + 10.*cos(qp1*(2+im10)*2.*360/72.))/86400.

      do 217 ik=i60,Z$S 
      do 217 ij=1,L$S
        if (LATST(ij) .ge. -45. .and. LATST(ij) .le. 45.) fxtdi(ij,ik) = 
     >  dta*cos(latst(ij)*6.*qp1)*exp(-abs((zstr(ik) - zstr(i100))/5.))
     > +dta*.6*cos(latst(ij)*6.*qp1)*exp(-abs((zstr(ik)-zstr(i80))/3.5))
     > -dta*.3*cos(latst(ij)*6.*qp1)*exp(-abs(zstr(ik) - zstr(i86)))
 217     continue


       dtak = 60. + 40.*cos(qp1*(20+im10)*2.*360./72.)

         do 218 ik=i60,i100
         do 218 ij=1,L$S
       if (LATST(ij) .ge. -45. .and. LATST(ij) .le. 45.)  kzzdi(ij,ik) = 
     >     dtak*exp((zstr(ik)-zstr(i100))/7.)*exp(-((latst(ij)/20.)**2))
 218     continue

         do 219 ik=i100+1,Z$S
         do 219 ij=1,L$S
             if (LATST(ij) .ge. -45. .and. LATST(ij) .le. 45.) 
     >           kzzdi(ij,ik) = kzzdi(ij,i100)
 219     CONTINUE 

C
C   kzzl is in m2/sec   ;;;  fxtl is in m/sec2,  FXTLON(13,L$S,Z$S), KZZLON(13,L$S,Z$S) defined above
C
         DO 250 ik=1,Z$S
         DO 250 ij=1,L$S
         DO 250 ix=1,13
             kzzlon(ix,ij,ik) = kzzl(ix,ij,ik) + kzzdi(ij,ik)
 250         fxtlon(ix,ij,ik) = fxtl(ix,ij,ik) + fxtdi(ij,ik)


C  Now smooth values 2X, using bi-linear smoother, re-load into KZZLON, FXTLON arrays for return
C     REAL*8 kz0(L$S,Z$S), fx0(L$S,Z$S), do for each longitude

       DO 280 ix=1,13

          do 283 ik=1,Z$S
          do 283 ij=1,L$S
            kz0(ij,ik) = kzzlon(ix,ij,ik)
 283        fx0(ij,ik) = fxtlon(ix,ij,ik)

              CALL SMOOTH5(kz0, kz1, L$S, Z$S)
              CALL SMOOTH5(kz1, kz2, L$S, Z$S)

              CALL SMOOTH5(fx0, fx1, L$S, Z$S)
              CALL SMOOTH5(fx1, fx2, L$S, Z$S)

            do 600 ik=1,Z$S
            do 600 ij=1,L$S
                kzzlon(ix,ij,ik) = kz2(ij,ik)
 600            fxtlon(ix,ij,ik) = fx2(ij,ik)
 280   CONTINUE


C  Now in LOWER STRATOSPHERE, if level 16 (12 mb) Kzz is greater than .1 m2/sec (1000 cm2/sec), 
C    then re-interpolate Kzz exponentially from the value at level 20 (4 mb - ~39 km) (ie, interpolate 
C    over levels 17-19) to .1 m2/sec at level 16. Set levels 1-15 (16 mb) to 0.0. These will 
C    be specified in TROPKZZ.  Also scale the Fxt array by the diffusion changes.
C  
C    new for Sept. 2005, use calculated values down to 35 km (use i35, defined above), 
C        then interpolate to 0 below 25 km;   zstr(Z$S), prst(Z$S);  FXTLON(13,L$S,Z$S), KZZLON(13,L$S,Z$S) 
C        best to interpolate in log, get a more natural taper off, 
C        assume MINUTE values at 25 km:  Kzz=.01 m2/sec, and Fxt = .01 m/sec/day = 1.157e-7 m/sec2
C
C    UPDATED, May 2006:
C        use calculated values down to 45 km (use i45, defined above), then interpolate to 0 below 35 km
C        assume MINUTE values at 35 km:  Kzz=.01 m2/sec, and Fxt = .01 m/sec/day = 1.157e-7 m/sec2


      fxt35 = .01/86400.

      do 620 ij=1,L$S
      do 620 ix=1,13

        do 625 ik=i35+1,i45-1
          fxt0 = EXP( ALOG( ABS(fxtlon(ix,ij,i45))/fxt35 ) /
     >     (zstr(i45) - zstr(i35))*(zstr(ik) - zstr(i35)) + ALOG(fxt35))
                                                                            ! reset to sign of original fxt
          fxtlon(ix,ij,ik) = SIGN(fxt0, fxtlon(ix,ij,i45))

          kzzlon(ix,ij,ik) = EXP( ALOG( ABS(kzzlon(ix,ij,i45))/.01) /
     >      (zstr(i45) - zstr(i35))*(zstr(ik) - zstr(i35)) + alog(.01) )
 625         continue

             do 627 ik=1,i35
                 fxtlon(ix,ij,ik) = 0.0
                 kzzlon(ix,ij,ik) = 0.0
 627         continue

C                                                don't do this part, just use Kzz, Fxtt as is
C       if (kzzlon(ix,ij,16,im10) .gt. .1) then
C          do 622 ik=17,19
C              kzzint = exp((alog(kzzlon(ix,ij,20,im10)) - 
C     c                        (-2.3))/4.*(ik-16.)-2.3)
C         
C              fxtlon(ix,ij,ik,im10) = fxtlon(ix,ij,ik,im10)*kzzint
C     c                                   /kzzlon(ix,ij,ik,im10)
C 622          kzzlon(ix,ij,ik,im10) = kzzint
C
C              fxtlon(ix,ij,16,im10) = fxtlon(ix,ij,16,im10)*.1
C     c                                  /kzzlon(ix,ij,16,im10)
C              kzzlon(ix,ij,16,im10) = .1
C       endif
C                            only set to zero at and below ik0 (bottom boundary, ~325 mb), originally level 4
ccc             do 624 ik=1,ik0
ccc                fxtlon(ix,ij,ik) = 0.0
ccc 624            kzzlon(ix,ij,ik) = 0.0

 620   CONTINUE


C  need to set max Kzz (ZONAL MEAN ONLY) to 100 m2/sec for use in model, 
C      LOAD ZONAL MEANS INTO XKZZT, FXTT arrays 
C                                       
Cc         go to 550
c

      do 410 ik=1,Z$S 
      do 410 ij=1,L$S

      do 412 ix=13,13
 412     if (kzzlon(ix,ij,ik) .gt. 100.) kzzlon(ix,ij,ik) = 100.d0

        xkzzt(ij,ik) = kzzlon(13,ij,ik) 
        fxtt(ij,ik) = fxtlon(13,ij,ik) 
cc                                             - and set limits of +- 216 m/s/day on FXTT (+-.0025 m/sec2)

        if (fxtt(ij,ik) .le. -.0025) fxtt(ij,ik) = -.0025
        if (fxtt(ij,ik) .ge. .0025)  fxtt(ij,ik) = .0025
 410  CONTINUE


C  Now compute eddy heat flux convergence, based on smoothed Kzz's, and don't need to smooth term 1
C    Note: if just zonal mean, (igwlon=0), then ix=1,12 are 0.0, and ix=13 is the zonal mean
C         gwh1b(13,L$S,Z$S), gwh2b(13,L$S,Z$S)

      do 880 ik=1,Z$S 
      do 880 ij=1,L$S
      do 880 ix=1,13
 880     gwh2b(ix,ij,ik) = -kzzlon(ix,ij,ik)*xn2(ij,ik)*
     >                          717./(2.*rr1*1004.)
C
C    if longitudinal variations turned on, compute zonal mean of eddy heat flux convergence, term1, term2  
C     
      if (igwlon .eq. 1) then 
        do 890 ik=1,Z$S 
        do 890 ij=1,L$S
           sumh1 = 0.0
           do 904 ix=1,12  
 904          sumh1 = sumh1 + gwh1b(ix,ij,ik)
           gwh1b(13,ij,ik) = sumh1/12.

           sumh2 = 0.0
           do 905 ix=1,12  
 905          sumh2 = sumh2 + gwh2b(ix,ij,ik)
           gwh2b(13,ij,ik) = sumh2/12.
 890    CONTINUE
      endif

C  load heat flux conver terms into arrays for output, use in STREAMF;  gwh1, gwh2 (L$S,Z$S) in DYNEXTRA

      do 950 ik=1,Z$S 
      do 950 ij=1,L$S
         gwh1(ij,ik) = gwh1b(13,ij,ik)
 950     gwh2(ij,ik) = gwh2b(13,ij,ik)
ccccccc 950     gwh2(ij,ik,im324) = gwh2b(13,ij,ik)

C
C Finally, do momentum forcing from turbulent diffusion of mean flow, 1/rho*d/dz(rho*Dm*du/dz),
C    if longitudinal variations are on (igwlon=1), then compute at each longitude first, then zonally avg
C    based on smoothed Kzz's, NOTE: need to use momentum diffusion, Dm = Kzz*Prandtl number  (Kzz or Dh)
                                   ! rhost(Z$S), DUDZG(13,L$S,Z$S), dmom(L$S,Z$S), kzzlon(13,L$S,Z$S)
                                   ! REAL*8 dm1(12,Z$S), dmz(12,Z$S)
      if (igwlon .eq. 1) then                                          

        do 770 ij = 1,L$S

         do 771 ik = 1,Z$S
         do 771 ix = 1,12
 771   dm1(ix,ik) = rhost(ik)*kzzlon(ix,ij,ik)*pran*dudzg(ix,ij,ik)

         CALL DERV4(1, dm1, dmz, delz, 12, Z$S, 2)

        do 772 ik = 1,Z$S
           sumom = 0.0
                do 773 ix = 1,12
 773            sumom = sumom + dmz(ix,ik)/rhost(ik)
 772     dmom(ij,ik) = sumom/12.
cccccccccccccccccc 772     dmom(ij,ik,im324) = sumom/12.

 770    CONTINUE
      ENDIF
C                              Now for zonal mean U, ! rhost(Z$S), DUDZG(13,L$S,Z$S), kzzlon(13,L$S,Z$S)
C                                                REAL*8 dm2(L$S,Z$S), dmz2(L$S,Z$S);  REAL dmom(L$S,Z$S)
      if (igwlon .eq. 0) then    

          do 780 ik = 1,Z$S
          do 780 ij = 1,L$S
 780   dm2(ij,ik) = rhost(ik)*kzzlon(13,ij,ik)*pran*dudzg(13,ij,ik)

         CALL DERV4(1, dm2, dmz2, delz, L$S, Z$S, 2)

           do 782 ik = 1,Z$S
           do 782 ij = 1,L$S
 782           dmom(ij,ik) = dmz2(ij,ik)/rhost(ik)
cccccccc 782           dmom(ij,ik,im324) = dmz2(ij,ik)/rhost(ik)

      ENDIF



CC 550        print *, ' END  OF  GRAVITY WAVE ROUTINE  '

         RETURN
         END


C
C   **************************   SMOOTHING  SUBROUTINE  ****************************************
C

        SUBROUTINE SMOOTH5(u1, u2, m, n)
 
c
c  subroutine to smooth 2D array,  u1(m,n),  and  return  smooooothed array u2(m,n)
C     smooth in the identical manner as smooth5.pro, adapted from smooth5.pro  (EF, 7/27/95)
C
C  *****  weights are:  .05, .075, .5
C     since, eg  x(-1,-1) is sqrt(2) further from x(0,0) than is x(-1,0)
C

        integer m, n
        REAL*8 u1(m,n), u2(m,n)


C  FIRST   SMOOTH all INTERIOR points

       do 200 ik=2,n-1
       do 200 ij=2,m-1

       u2(ij,ik) = .05*(u1(ij-1,ik-1) + u1(ij-1,ik+1) + u1(ij+1,ik-1) +
     c     u1(ij+1,ik+1)) + .075*(u1(ij-1,ik) + u1(ij,ik+1) +
     c     u1(ij,ik-1) + u1(ij+1,ik)) + .5*u1(ij,ik)
 200       continue


C  *****Compute smoooooooothed boundaries  -- new  (EF, 6/26/95)
C     use more weight to surrounding points for more smoothing, here and for corner points
  
      do 300 ik=2,n-1
          u2(1,ik) = .1*(u1(2,ik-1) + u1(2,ik+1)) + .35*u1(1,ik) +
     c              .15*(u1(1,ik-1) + u1(1,ik+1) + u1(2,ik))
 300  continue

      do 305 ik=2,n-1
          u2(m,ik) = .1*(u1(m-1,ik-1) + u1(m-1,ik+1)) + .35*u1(m,ik) +
     c              .15*(u1(m,ik-1) + u1(m,ik+1) + u1(m-1,ik))
 305  continue


      do 320 ij=2,m-1
          u2(ij,1) = .1*(u1(ij-1,2) + u1(ij+1,2)) + .35*u1(ij,1) +
     c              .15*(u1(ij-1,1) + u1(ij+1,1) + u1(ij,2))
 320  continue

      do 325 ij=2,m-1
          u2(ij,n) = .1*(u1(ij-1,n-1) + u1(ij+1,n-1)) + .35*u1(ij,n) +
     c              .15*(u1(ij-1,n) + u1(ij+1,n) + u1(ij,n-1))
 325  continue


c  *****Compute smoooooooooothed corner points

      u2(1,1) = .16*u1(2,2) + .24*(u1(1,2) + u1(2,1)) + .36*u1(1,1)
      u2(m,1) = .16*u1(m-1,2) + .24*(u1(m,2) + u1(m-1,1)) + .36*u1(m,1)

      u2(1,n) = .16*u1(2,n-1) + .24*(u1(2,n) + u1(1,n-1)) + .36*u1(1,n)
      u2(m,n)= .16*u1(m-1,n-1) +.24*(u1(m-1,n) + u1(m,n-1)) +.36*u1(m,n)
  

        RETURN
      END


C
C   **************************   SMOOTHING  SUBROUTINE - for COUPLED MODEL    *******************************
C                                                         less smoothing than SMOOTH5 above

        SUBROUTINE SMOOTH5C(u1, u2, m, n)
 
c
c  subroutine to smooth 2D array,  u1(m,n),  and  return  smooooothed array u2(m,n)
C
C  *****  weights are:  .03, .045, .7
C     since, eg  x(-1,-1) is sqrt(2) further from x(0,0) than is x(-1,0)
C

        integer m, n
        REAL*8 u1(m,n), u2(m,n)


C  FIRST   SMOOTH all INTERIOR points

       do 200 ik=2,n-1
       do 200 ij=2,m-1

       u2(ij,ik) = .03*(u1(ij-1,ik-1) + u1(ij-1,ik+1) + u1(ij+1,ik-1) +
     c     u1(ij+1,ik+1)) + .045*(u1(ij-1,ik) + u1(ij,ik+1) +
     c     u1(ij,ik-1) + u1(ij+1,ik)) + .7*u1(ij,ik)
 200       continue


C  *****Compute smoooooooothed boundaries  -- new  (EF, 6/26/95)
C     use more weight to surrounding points for more smoothing, here and for corner points
  
      do 300 ik=2,n-1
          u2(1,ik) = .06*(u1(2,ik-1) + u1(2,ik+1)) + .61*u1(1,ik) +
     c               .09*(u1(1,ik-1) + u1(1,ik+1) + u1(2,ik))
 300  continue

      do 305 ik=2,n-1
          u2(m,ik) = .06*(u1(m-1,ik-1) + u1(m-1,ik+1)) + .61*u1(m,ik) +
     c               .09*(u1(m,ik-1) + u1(m,ik+1) + u1(m-1,ik))
 305  continue


      do 320 ij=2,m-1
          u2(ij,1) = .06*(u1(ij-1,2) + u1(ij+1,2)) + .61*u1(ij,1) +
     c               .09*(u1(ij-1,1) + u1(ij+1,1) + u1(ij,2))
 320  continue

      do 325 ij=2,m-1
          u2(ij,n) = .06*(u1(ij-1,n-1) + u1(ij+1,n-1)) + .61*u1(ij,n) +
     c               .09*(u1(ij-1,n) + u1(ij+1,n) + u1(ij,n-1))
 325  continue


c  *****Compute smoooooooooothed corner points

      u2(1,1) = .1*u1(2,2)    + .15*(u1(1,2) + u1(2,1)) +     .6*u1(1,1)
      u2(m,1) = .1*u1(m-1,2)  + .15*(u1(m,2) + u1(m-1,1)) +   .6*u1(m,1)

      u2(1,n) = .1*u1(2,n-1)  + .15*(u1(2,n) + u1(1,n-1)) +   .6*u1(1,n)
      u2(m,n) = .1*u1(m-1,n-1) +.15*(u1(m-1,n) + u1(m,n-1)) + .6*u1(m,n)
  

        RETURN
      END




C   **************************   SORTING  SUBROUTINE  ****************************************
C
C
       SUBROUTINE SORTOUT(ucc, Z$S, is, ibrlev, ik, index)
C
C        routine to figure out index of maximum value,  ucc(11,Z$S), ibrlev(12), 
C       is = number of waves at the given Zbr, ik=current level
C
C    OUTPUT:>  index = index of icc needed (i.e., -50,-40,-30,....40,50, phase speed waves)
C

       INTEGER ibrlev(12), Z$S
       REAL ucc(11,Z$S)

       utest = abs(ucc(ibrlev(1),ik))
       index = ibrlev(1)

       ij = 2

 1000  if (utest .lt. abs(ucc(ibrlev(ij),ik))) then 
             utest = abs(ucc(ibrlev(ij),ik))
             index = ibrlev(ij)
       endif

       ij = ij + 1
       if (ij .le. is) go to 1000


       RETURN
       END

