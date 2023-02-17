c
c
      SUBROUTINE TPCORE_DYN(np$, mp$, zres, dt, xc0, wrbr, xcnew)
c
c
c *******************************************************************************
c    Routine to do advection on 2D model grid, EVERYTHING IN DOUBLE PRECISION
c      -- E. Fleming (Aug. 96) 
C            - updated Nov. 2001 to handle an irregular grid in the vertical
cc ******************************************************************************
C
C Purpose: advects constituents from given input v* and w* on a C-grid.
c          The current code assumes that the individual Courant numbers
c          are never greater than 1 which is well within the limits of the 
c          2D model circulation for a 12 hour time step.  Inclusion of
c          the semi-lagrangian code from TPCORE.F is needed if Courant
c          numbers are greater than 1. This is relatively easy to do.
c
C Schemes: Piecewise parabolic method (PPM), with cross-derivative 
c          terms to account for multidirectional advection 
c          without the operator splitting that the Prather scheme imposes
C          (see Lin and Rood 1994, 1996; also see Colella and Woodward, 1984, 
c          and Carpenter et al., 1990).
C
C This is adopted from the TPCORE.F routine (version 3.5) for 3D transport
c written by Shian-Jiann Lin and Richard B. Rood for the GSFC 3D CTM, 
c the GEOS-GCM and the GSFC Data Assimilation System (GEOS-DAS). 1D and
c 2D (x-y) advection codes utilizing the PPM written by S.J. Lin were also used.
c
c This code can accomodate the following options on the PPM:
C
C  IORD= 2: uses Improved full monotonicity constraint (Lin, 2000) (otherwise same as IORD=3)
C                this is significantly different from IORD=3, is more diffusive 
C                (doesn't maintain gradients as well).
C        3: monotonic PPM* (slightly improved PPM of Collela & Woodward 1984) -
C           uses 4th order reference slopes
C        4: semi-monotonic PPM (same as 3, but overshoots are allowed)
C        5: positive-definite PPM (constraint on the subgrid distribution is
C           only strong enough to prevent generation of negative values;
C           both overshoots & undershootes are possible).
C        6: un-constrained PPM (nearly diffusion free; slightly faster but
C           positivity not quaranteed. Use this option only when the fields
C           and winds are very smooth).
C
C Very little difference in the 2D model simulations of CH4 was found between
c these 4 options, since standard 2D model fields are rather smooth. Note,
c however, that the un-constrained option did show a small area of negative
c mixing ratios. Because of this, and since the numerical diffusion seems to be
c quite small (actually smaller than Prather, probably because of the operator 
c splitting), Lin recommends just using the fully montonic
c constraint (IORD=2 or 3), since this guarantees exact tracer correlations. 
C Best to use IORD=3 - it maintains gradients better than IORD=2 
C (unless I'm not doing IORD=2 correctly).
C
C Note that the implicit numerical diffusion decreases as IORD increases.
C The last two options (IORD=5, 6) should only be used when there is
C significant explicit diffusion (such as a turbulence parameterization). You
C might get dispersive results otherwise.
C
C PPM is 4th order accurate when grid spacing is uniform (y & z) as
C     with the 2D model; 3rd order accurate for non-uniform grid.
C Again, no space or time splitting is performed.
C ********************************************************************
C
c
c   XC0(NP$,MP$) for theta, ang mom., REAL*8 ,  THINGS ACT ON MIXING RATIO
C
C
C   WBAR(2,MP$+1), W1(MP$+1), W(NP$,MP$+1), V(NP$+1,MP$), COSC(NP$), cose(NP$+1)
C
C
C  *************************************************************************************************
C    
C  THIS HAS BEEN MODIFIED TO ADVECT COUPLED MODEL THETA AND ANGULAR MOMENTUM, treated as mixing ratio
C  (otherwise it won't work if COS=1 and MS=1)
C  IT'S ALSO BEEN AUTOMATED TO HANDLE ANY INPUT GRID BASED ON NP$, MP$, and ZRES (EF, March 2008)
C
C  *************************************************************************************************
C
C
ccccc       include "com2d.h"
ccccc       INCLUDE 'PARAM.INC'
ccccc       INCLUDE 'COMMONC.INC'
ccccc       INCLUDE 'COMMONT.INC'

!       USE degree_trig

       integer NP$, MP$, npe, mpe
       real*4 dt, zres, wrbr(NP$,MP$)

       real*8 w1(MP$), ww(MP$+1), XC0(NP$,MP$), XCNEW(NP$,MP$)
       real*8 dtdy1, dtdze(MP$), dzz(MP$), ype(NP$+1), zpe(MP$+1)
       real*8 w8(NP$,MP$+1), wbar(2,MP$+1), w(NP$,MP$+1), v(NP$+1,MP$)
       real*8 cosc(NP$), cose(NP$+1), tcos, press8(MP$), presse8(MP$+1)
       real*8 ypp(NP$), zpp(MP$), dy8, dt8, dth, dzres

       real*8 XMR0(NP$,MP$), mse(MP$+1), ms(MP$), DP(MP$), deltazp(MP$)

       real*8 dqv(NP$,MP$), dqw(NP$,MP$), dqt(NP$,MP$), ffb, fft
       real*8 px(NP$), pcr(NP$), ptt(NP$), rx(MP$), rcr(MP$), rtt(MP$)
       real*8 px1(0:NP$+1), rx1(0:MP$+1), ady(NP$,MP$), adz(NP$,MP$)  
       real*8 vc(NP$+1), ymass(NP$+1), vav(NP$)
       real*8 wc(MP$+1), zmass(MP$+1), wav(MP$)
       real*8 fymr(NP$+1), ffy(NP$+1), fzmr(MP$+1), ffz(MP$+1)
       real*8 sum1, sum2, sumd, zzzz

       LOGICAL LMASS  

ccc       DATA IORD/2/
       DATA IORD/3/

       SAVE


C  define latitudes and Z-alt (km) at box centers and edges as appropriate for COUPLED MODEL dynamics grid
C  also define cosines at centers and edges:   follow as done in SETCON in TWODLIB

        dt8 = DBLE(dt)
        dzres = DBLE(zres)*7.d0

        dth = 180.d0/NP$
        dy8 = 6371000.d0*100.d0*3.1415927d0/NP$


c  DTDY1 is in sec/cm, DY8 (CM), dt8 is in seconds

        DTDY1 = dt8/dy8


        tcos = 0.d0
        do 505 ij=1,NP$ 
           ypp(ij) = -90.d0 + (ij - .5d0)*dth
           cosc(ij) = DCOSD(ypp(ij))
           tcos = tcos + cosc(ij)
 505    CONTINUE


        do 507 ij=1,NP$+1 
           ype(ij) = -90.d0 + (ij - 1.d0)*dth
 507       cose(ij) = DCOSD(ype(ij))

           cose(1) = 0.d0
           cose(NP$+1) = 0.d0



        do 510 ik=1,MP$
 510       zpp(ik) = (ik - .5d0)*dzres

        do 511 ik=1,MP$+1
 511       zpe(ik) = (ik - 1.d0)*dzres


C
c  define MSE array at box edges using PRESSE8(MP$+1), number density in 1/cm3 at constant T = 239.27K, 
C    also define MS array and deltazp (MP$) at box centers using PRESS8(MP$), DP(MP$)
C
        do 103 ik=1,MP$+1
           presse8(ik) = 1000.d0*DEXP(-zpe(ik)/7.d0)
 103       mse(ik) = presse8(ik)*1.d3/(28.97*2.87d6*239.27*1.66d-24)

        do 104 ik=1,MP$
           press8(ik) = 1000.d0*DEXP(-zpp(ik)/7.d0)
           ms(ik) = press8(ik)*1.d3/(28.97*2.87d6*239.27*1.66d-24)
           dp(ik) = -DBLE(zres)      ! DLOG(presse8(ik+1)/presse8(ik))
 104       deltazp(ik) = -29.3*239.27*1.D-3*DP(ik)



C  define dz and dtdz between BOX EDGES in DTDZE(MP$), DZZ(MP$), ZPE(MP$+1), convert to CM

        do 110 ik=1,MP$
           DZZ(ik) = (zpe(ik+1) - zpe(ik))*1000.d0*100.d0
 110       DTDZE(ik) = dt8/dzz(ik)



C  interpolate w* field to top and bottom edges of grid boxes, reload to REAL*8, ww(MP$+1), w8(NP$,MP$+1)
C                        zpp(MP$),  set to 0.0 at bottom and top edge - W* is in CM/SEC

           NPE = NP$+1
           MPE = MP$+1


           do 420 ij=1,NP$

              do 221 ik=1,MP$
 221               w1(ik) = WRBR(ij,ik)

           CALL LINTERP8(zpp, MP$, W1, zpe, MPE, 0, WW)

              do 222 ik=2,MP$
 222             w8(ij,ik) = ww(ik)

                 w8(ij,1) = 0.d0
                 w8(ij,MP$+1) = 0.d0
 420       CONTINUE



C     Currant condition w/ current coupled model grid (~5 degrees x ~2 km) is: max W* = 28.2717 cm/sec
C                                               cosc(NP$),  wbar(2,MP$+1),     max V* = 7532.56 cm/sec
         DO 300 ik=1,MP$+1 

               WBAR(1,ik) = 0.0
            do 304 ij=1,NP$
 304           WBAR(1,IK) = WBAR(1,ik) + w8(IJ,IK)*cosc(ij)

            do 305 ij=1,NP$
 305           w(ij,ik) = w8(ij,ik) - WBAR(1,ik)/tcos
                                                                  ! check adjusted w* for mass balance
               WBAR(2,ik) = 0.0
            do 307 ij=1,NP$
 307           WBAR(2,ik) = WBAR(2,ik) + w(IJ,ik)*cosc(IJ)
 300     CONTINUE

c                                       ensure that w* at bottom and top levels = 0.0 for adjusted array
         DO 404 ij=1,NP$
            w(ij,1) = 0.0
            w(ij,MP$+1) = 0.0
 404     continue

C                                            set up V array,  v* = 0 at poles by definition, v(NP$+1,MP$)
         do 400 ik=1,MP$                                   
               v(1,ik) = 0.d0  
 400           v(NP$+1,ik) = 0.d0

C  just use the exp(-z/H) factor for the density variation, since the rho(0) factors cancel 
C     - EVERYTHING IS CONVERTED TO  CM  HERE!!!! dy8 is in CM,  cosc(NP$), cose(NP$+1), zpe(MP$+1), zpp(MP$)

         do 450 ik=1,MP$ 
         do 450 ij=2,NP$
           v(ij,ik) = (v(ij-1,ik)*cose(ij-1) - 
     >                (w(ij-1,ik+1)*DEXP(-zpe(ik+1)/7.d0)  
     >               - w(ij-1,ik)*DEXP(-zpe(ik)/7.d0))
     > /(DEXP(-zpp(ik)/7.d0)*(zpe(ik+1)-zpe(ik))*1000.d0*100.d0)
     >                *cosc(ij-1)*dy8)/cose(ij)    
 450  CONTINUE



c Check  Mass conservation (total number of molecules) before (and after) advection every 60 days, CO2 here 
c  note - checking mass conservation uses M*DELTAZ quantity, in which TEMPERATURE cancels, 
C    so just use constant temp for both for simplicity, also, total area of earth not included (constant)
c    XC0(NP$,MP$)

cc        LMASS = .TRUE.
        LMASS = .FALSE.

       IF (LMASS) THEN
ccccccc        if (mod(iday360-15,60) .eq. 0.0) then
         sum1 = 0.d0
           do 990 ik=1,MP$
           do 990 ij=1,NP$
 990        sum1 = sum1 + xc0(ij,ik)*ms(ik)*deltazp(ik)*cosc(ij)/tcos
ccccccc        endif
       ENDIF

c                                        ****************  loop through dynamical species *********
cc      DO 1000 INS=1,2     ! ITRANS

c                            first load mixing ratios into temporary array
        do 150 ik=1,MP$
        do 150 ij=1,NP$
 150       XMR0(ij,ik) = XC0(ij,ik)


C                         ***** Y  SETUP  --  COMPUTE X-TERMS  ****************
      DO 160 ik=1,MP$ 
             do 162 ij=1,NP$
 162            px1(ij) = xmr0(ij,ik)                
c                                     extend array one grid point each side, Lin says just use constant value
                  px1(0) = px1(1)
                  px1(NP$+1) = px1(NP$)            

             do 164 ij=1,NP$+1 
 164            vc(ij) = v(ij,ik)*dtdy1                        ! vc(NP$+1), Courant number

        do 166 ij=1,NP$                                        ! vav(NP$), box-averged Courant number
 166       vav(ij) = (vc(ij) + vc(ij+1))*.5

            CALL CROSST1(NP$, VAV, PX1, PCR, PTT, IORD)   !px1(0:NP$+1), PCR, PTT(NP$)=total (PX1+PCR) x-terms

             do 168 ij=1,NP$          
 168            ady(ij,ik) = ptt(ij)                     
 160  CONTINUE


C                            ***** Z  SETUP  -- COMPUTE X-TERMS  ****************
C
C  with irregular grid, use the UPWIND BOX to define the Courant number across the box edge
C     i.e., the mass that is moving across the box edge, and W, WC = 0 at bottom and top 
C
C  w(NP$,MP$+1), DTDZE(MP$), WC(MP$+1) 
C
       DO 170 ij=1,NP$               
              do 172 ik=1,MP$
 172             rx1(ik) = xmr0(ij,ik)
c                                     extend array one grid point each side, Lin says just use constant value
                 rx1(0) = rx1(1)
                 rx1(MP$+1) = rx1(MP$)            
 
         do 174 ik=2,MP$
             if (w(ij,ik) .ge. 0.d0) wc(ik) = w(ij,ik)*DTDZE(ik-1) 
             if (w(ij,ik) .lt. 0.d0) wc(ik) = w(ij,ik)*DTDZE(ik) 
 174     continue
                                            ! wc(MP$+1), Courant number across box edges
         wc(1) = 0.d0
         wc(MP$+1) = 0.d0

        do 175 ik=1,MP$                                           ! wav(MP$), box-averged Courant number
 175       wav(ik) = (wc(ik) + wc(ik+1))*.5

           CALL CROSST1(MP$, WAV, RX1, RCR, RTT, IORD)    !rx1(0:MP$+1), RCR,RTT(MP$) = total(RX1+RCR) x-terms

              do 177 ik=1,MP$         
 177            adz(ij,ik) = rtt(ik)   
 170  CONTINUE


C                 ******** Do  Y ADVECTION  [conservative (liberal) transport]  ******

      DO 210 ik=1,MP$                                                   !!  v(NP$+1,MP$), adz, dqv(NP$,MP$)
            do 212 ij=1,NP$+1
              vc(ij) = v(ij,ik)*DTDY1    
 212          ymass(ij) = vc(ij)*ms(ik)*cose(ij)          ! horizontal mass flux; ymass, vc, cose(NP$+1)
                                                          !    use ms(MP$) at grid points (box centers)
            do 215 ij=1,NP$
 215           px(ij) = adz(ij,ik)                       ! px(NP$), use adz array here because of cross terms

               CALL FYPPM1(NP$, VC, PX, fymr, IORD)

         do 600 ij=1,NP$+1                                      ! constit. mass flux across box edge
 600        ffy(ij) = fymr(ij)*ymass(ij)                       ! ymass=0 and so fx=0 at poles (90s, 90N)

         do 610 ij=1,NP$                                        ! constit mass flux convergence at box center
 610       dqv(ij,ik) = -(ffy(ij+1) - ffy(ij))/cosc(ij)         
 210  CONTINUE
C
c  Note, get a compile warning, "out of bounds array reference" on the COSC array, 610 statement above
c  this is gremlin-like, I checked this THOROUGHLY and COULDN'T figure out why the warning occurs, the
C array is properly dimensioned; get the same error with NP$ or 18, the problem is only with -O3 optimization
C
C
C
C                         ******** Z (vertical)  ADVECTION   ****************
C
C  with irregular grid, use the UPWIND BOX to define the Courant number across the box edge
C     i.e., the mass that is moving across the box edge, and W, WC at bottom and top = 0
C
C w(NP$,MP$+1), ady, dqw(NP$,MP$), DTDZE(MP$),mse, WC(MP$+1), zmass(MP$+1)=vertical mass flux across box edge
C
      DO 220 ij=1,NP$                                        

         do 223 ik=2,MP$
             if (w(ij,ik) .ge. 0.d0) wc(ik) = w(ij,ik)*DTDZE(ik-1) 
             if (w(ij,ik) .lt. 0.d0) wc(ik) = w(ij,ik)*DTDZE(ik) 
 223     continue
       
         wc(1) = 0.d0
         wc(MP$+1) = 0.d0
                                  ! zmass should use the M at the box EDGE (not center) here since it gets
                                  ! multiplied by the mixing ratio at the edge (M and q are always together) 

         do 224 ik=1,MP$+1
 224         zmass(ik) = w(ij,ik)*mse(ik)
                                             ! at boundaries, set zmass=0 since w=0
             zmass(1) = 0.d0
             zmass(MP$+1) = 0.d0


            do 225 ik=1,MP$
 225           rx(ik) = ady(ij,ik)      ! rx(MP$), use ady array here for mixing ratio because of cross terms


               CALL FZPPM1(MP$, DZZ, WC, RX, fzmr, IORD)


         do 700 ik=1,MP$+1                                 ! constit. mass flux across box edge - top/bottom
 700        ffz(ik) = fzmr(ik)*zmass(ik)*dt8               ! zmass=0 and so fx=0 at top and bottom boundary

         do 710 ik=1,MP$                                   ! constit mass flux convergence at box center
 710        dqw(ij,ik) = -(ffz(ik+1) - ffz(ik))/dzz(ik)    ! need to scale it by the size of the box, 
                                                           ! per eqn. 1.13 in C&W, 1984

C  to compute constit mass flux convergence at box center, need to 
C     account for the size of the upwind box relative to the box in question, 
C     DZZ=fltarr(MP$)  - YES!!! THIS IS THE CORRECT WAY TO DO IT, dqw(NP$,MP$)
C 
C  need to do top and bottom boxes separately
C
C  NEW NOTE:  THIS TURNS OUT TO BE IDENTICAL TO  -(ffz(ik+1) - ffz(ik))*dt8/dzz(ik) - the 710 loop ABOVE
C
CC       if (ffz(2) .ge. 0.d0) fft = ffz(2)*dzz(1)
CC       if (ffz(2) .lt. 0.d0) fft = ffz(2)*dzz(2)
CC
CC       dqw(ij,1) = -(fft - 0.d0)/dzz(1)
CC
CC       do 710 ik=2,MP$-1
CC          if (ffz(ik+1) .ge. 0.d0) fft = ffz(ik+1)*dzz(ik)
CC          if (ffz(ik+1) .lt. 0.d0) fft = ffz(ik+1)*dzz(ik+1)
CC
CC          if (ffz(ik) .ge. 0.d0) ffb = ffz(ik)*dzz(ik-1)
CC          if (ffz(ik) .lt. 0.d0) ffb = ffz(ik)*dzz(ik)
CC
CC          dqw(ij,ik) = -(fft - ffb)/dzz(ik)
CC 710     CONTINUE
CC
CC       if (ffz(MP$) .ge. 0.d0) ffb = ffz(MP$)*dzz(MP$-1)
CC       if (ffz(MP$) .lt. 0.d0) ffb = ffz(MP$)*dzz(MP$)
CC
CC       dqw(ij,MP$) = -(0.d0 - ffb)/dzz(MP$)

 220  CONTINUE


c  update new mixing ratio at box center, take old mixing ratio and add on change
c   in number density in both horiz and vertical, converted back to mixing ratio 
c  (don't divide by cosine here, it's already been done)  --  dqv, dqw, dqt (NP$,MP$)  
c
                                                            !! need to scale by the ms at box center
        do 330 ik=1,MP$
        do 330 ij=1,NP$
 330       dqt(ij,ik) = (dqv(ij,ik) + dqw(ij,ik))/ms(ik)


        do 440 ik=1,MP$
        do 440 ij=1,NP$
           XCNEW(ij,ik) = XMR0(ij,ik) + dqt(ij,ik)                      ! UPDATE w/ mixing ratio change
c                                                            sometimes get negatives in solid HNO3, H2O
c                                                            because of sharp density changes, readjust here
cc           zzzz = XCNEW(ij,ik)*m(ij,ik)  
cc           if (zzzz .le. 1.D-12) XCNEW(ij,ik) = 1.D-12/m(ij,ik)  
 440    continue


C   Now load in advective tendency terms here for output in MAIN;  
C               ADVY, ADVZ(T$,NP$,MP$) are in MIX RAT/SEC (DT8 is in seconds)

cc        do 444 ik=1,MP$
cc        do 444 ij=1,NP$
cc           advy(ij,ik) = dqv(ij,ik)/ms(ik)/DT8
cc 444       advz(ij,ik) = dqw(ij,ik)/ms(ik)/DT8


ccc 1000   CONTINUE                                                              ! end species loop


c Check  Mass conservation (total number of molecules) (before) and after advection every 60 days, SPECIES 2

       IF (LMASS) THEN
ccccccc        if (mod(iday360-15,60) .eq. 0.0) then
          sum2 = 0.d0
           do 992 ik=1,MP$
           do 992 ij=1,NP$
 992       sum2 = sum2 + xcnew(ij,ik)*ms(ik)*deltazp(ik)*cosc(ij)/tcos
             sumd = (sum2-sum1)/sum1
             write(6,995) iday360, iyr, sum1, sum2, sumd
 995         format(1x, 2I3, 3x, 1P3D24.15)
ccccccc        endif
       ENDIF
    

ccccccc       write (1991) NP$, MP$, dt8, dy8, dtdy1
ccccccc       write (1991) ypp, zpp, ype, zpe, dp, press8, presse8, dtdze
ccccccc       write (1991) dzz, deltazp, cosc, cose, tcos, wrbr, w8, wbar, w, v
ccccccc       write (1991) ms, mse, xc0, xcnew, sum1, sum2, sumd


	RETURN
	END

CC *************************  SUBROUTINES **************************************************************


          subroutine FYPPM1(L$, VC, PX, fluxt, IORD)

          integer L$
          real*8 R3, R23, R24
          parameter (R3 = 1.d0/3., R23 = 2.d0/3., R24 = 1./24.d0)

          real*8 VC(L$+1), px(L$), DC(L$)
          real*8 fluxt(L$+1), AL(L$), AR(L$), A6(L$), tmp, pmax, pmin 


          ILMT = IORD - 3


       do 100 ij=2,L$-1
          if (ij .ge. 3  .and. ij .le. L$-2) then
            tmp = (8.d0*(px(ij+1) - px(ij-1)) + px(ij-2) - px(ij+2))*R24   !4th order difference for interior
          else
            tmp = .25*(px(ij+1) - px(ij-1))                                !2nd order difference for 1 pt on
          endif                                                                    ! e.g., 75S,75N

          Pmax = DMAX1(px(ij-1), px(ij), px(ij+1)) - px(ij)
          Pmin = px(ij) - DMIN1(px(ij-1), px(ij), px(ij+1))
 100      DC(ij) = DSIGN(DMIN1(DABS(tmp), Pmax, Pmin), tmp)           ! sign func: result = abs(1st)*sgn(2nd)
 
          DC(1) = 0.d0
          DC(L$) = 0.d0

cc  Compute values at box edges, for poles, extrapolate 1/2 grid point
cc     to get first guess for AL(1) and AR(L$) - these will be modified in LMTPPM
cc     just use constant gradient (ratio) for extrapolation, need sqrt since its 1/2 grid point
c
        do 400 ij=2,L$
 400       AL(ij) = 0.5d0*(px(ij-1)+px(ij)) + (DC(ij-1) - DC(ij))*R3
              AL(1) = DSQRT(px(1)/px(2))*px(1)

        do 402 ij=1,L$-1
 402       AR(ij) = AL(ij+1)
              AR(L$) = DSQRT(px(L$)/px(L$-1))*px(L$)  

        do 404 ij=1,L$
 404       A6(ij) = 3.d0*(px(ij)+px(ij)  - (AL(ij)+AR(ij)))


       if (ILMT .le. 2) call LMTPPM1(L$, PX, DC, A6, AR, AL, ILMT)          ! adjusts AL, AR, A6 arrays


       do 450 ij = 2,L$
          IF (VC(ij) .GT. 0.) then
             fluxt(ij) = AR(ij-1) + 0.5d0*VC(ij)*(AL(ij-1) - AR(ij-1) + 
     c                     A6(ij-1)*(1.-R23*VC(ij)) )
          else
             fluxt(ij) = AL(ij) - 0.5d0*VC(ij)*(AR(ij) - AL(ij) + 
     c                     A6(ij)*(1.+R23*VC(ij)))
          ENDIF
 450   continue
c                                                  at poles, VC = 0 and YMASS = 0, so it doesn't matter
c                                                  what fluxt(poles) is, just set to zero to define
       fluxt(1) = 0.d0 
       fluxt(L$+1) = 0.d0                         
                                                  
      RETURN                                      
      END                                         



          SUBROUTINE FZPPM1(Z$X, DZZ, WC, PX, fluxt, IORD)

          integer Z$X
          real*8 R3, R23, R24
          parameter(R3= 1.d0/3., R23 = 2.d0/3., R24 = 1./24.d0)

          real*8 WC(Z$X+1), px(Z$X), DC(Z$X), tmp, pmax, pmin
          real*8 fluxt(Z$X+1), AL(Z$X), AR(Z$X), A6(Z$X), dzz(Z$X)

                                             !dzz here is the dz between BOX EDGES 

          ILMT = IORD - 3


       do 100 ij=2,Z$X-1
Cef                                         these centered difference formulas for regular grid not used here
cef       if (ij .ge. 3  .and. ij .le. Z$58-2) then
cef         tmp = (8.d0*(px(ij+1) - px(ij-1)) + px(ij-2) - px(ij+2))*R24  !4th order difference for levs 3-56
cef          else
cef            tmp = .25*(px(ij+1) - px(ij-1))                            !2nd order difference for levs 2,57
cef       endif
Cef
C instead, do centered differencing for irregular grid from eqn 1.7 in C&W, 1984,  (tmp=centered difference)
C   dzz=fltarr(Z$X), px=fltarr(Z$X), 
C   NOTE, we have altered the DC array computation to match eqn 1.8 in C&W, whereas the previous code carried
C         a factor of 1/2 through all computations that used the DC array. Here, we will be consistent 
C         w/ the equations in C&W throughout the code (ie, for all uses of the DC array).
C         (although the factor of 2 made very LITTLE difference in the simple test profile).
C
C  also, PMAX and PMIN are the left and right differences, and DC=minimum of the left, right, and center 
C        differences, and DC=0 if px(ij) is an extremum - see Fig. 2 in Carpenter et al., and eqn. 1.8 in C&W
C

         tmp = dzz(ij)/(dzz(ij-1) + dzz(ij) + dzz(ij+1))
     >         *( (2.*dzz(ij-1) + dzz(ij))/(dzz(ij+1) + dzz(ij))
     >         *(px(ij+1) - px(ij))
     >       + (dzz(ij) + 2.*dzz(ij+1))/(dzz(ij-1) + dzz(ij))
     >         *(px(ij) - px(ij-1)) )

          Pmax = DMAX1(px(ij-1), px(ij), px(ij+1)) - px(ij)
          Pmin = px(ij) - DMIN1(px(ij-1), px(ij), px(ij+1))
 100      DC(ij) = DSIGN(DMIN1(DABS(tmp), 2.d0*Pmax, 2.d0*Pmin), tmp)  !sign func: result = abs(1st)*sgn(2nd)
 
          DC(1) = 0.d0
          DC(Z$X) = 0.d0


cc  Compute values at box edges, for top/bottom boundary, extrapolate 1/2 grid point
cc     to get first guess for AL(0) and AR(Z$X-1) - these will be modified in LMTPPM
cc     just use constant gradient (ratio) for extrapolation, need sqrt since its 1/2 grid point
c
c for irregular grid, use eqn 1.6 in C&W, 1984 for AL,   dzz=fltarr(Z$X), px=fltarr(Z$X) 
c
      do 400 ij=2,Z$X-2
           AL(ij+1) = px(ij) + dzz(ij)/(dzz(ij+1) + dzz(ij))
     >                        *(px(ij+1) - px(ij)) 
     >       +  1./(dzz(ij-1) + dzz(ij) + dzz(ij+1) + dzz(ij+2))
     >        *( 2.*dzz(ij+1)*dzz(ij)/(dzz(ij) + dzz(ij+1))
     >        *((dzz(ij-1) + dzz(ij))/(2.*dzz(ij) + dzz(ij+1)) 
     > - (dzz(ij+2)+dzz(ij+1))/(2.*dzz(ij+1)+dzz(ij)))*(px(ij+1)-px(ij)) 
     > - dzz(ij)*(dzz(ij-1) + dzz(ij))/(2.*dzz(ij) + dzz(ij+1))*DC(ij+1)
     > + dzz(ij+1)*(dzz(ij+1)+dzz(ij+2))/(dzz(ij)+2.*dzz(ij+1))*DC(ij) )
 400   CONTINUE
C
C  now do levels 2 and Z$X using the eqn on p. 589 in Carpenter et al. which is for a regular
C  grid, since these levels will always have a regular grid, 
C   i.e, the 2nd level above the ground, and the top level 
C   (which is in the 92-116 km sponge layer which will always have a regular grid)
C
        AL(2) = (px(1) + px(2))/2.d0 - (DC(2) - DC(1))/6.d0
        AL(Z$X)= (px(Z$X-1) + px(Z$X))/2.d0 - (DC(Z$X) - DC(Z$X-1))/6.d0

        AL(1) = DSQRT(px(1)/px(2))*px(1)


        do 402 ij=1,Z$X-1
 402       AR(ij) = AL(ij+1)
           AR(Z$X) = DSQRT(px(Z$X)/px(Z$X-1))*px(Z$X)  


        do 404 ij=1,Z$X
 404       A6(ij) = 3.d0*(px(ij)+px(ij)  - (AL(ij)+AR(ij)))


       if (ILMT .le. 2) CALL LMTPPM1(Z$X, PX, DC, A6, AR, AL, ILMT)          ! adjusts AL, AR, A6 arrays


C these fluxt calculations have all been checked to be correct, 
C   NOTE that the fluxt value at a box edge is computed from the WC at that edge
C
       do 450 ij = 2,Z$X
          IF (WC(ij) .GT. 0.) then
             fluxt(ij) = AR(ij-1) + 0.5d0*WC(ij)*(AL(ij-1) - AR(ij-1) + 
     c                     A6(ij-1)*(1.-R23*WC(ij)) )
          else
             fluxt(ij) = AL(ij) - 0.5d0*WC(ij)*(AR(ij) - AL(ij) + 
     c                     A6(ij)*(1.+R23*WC(ij)))
          ENDIF
 450   continue
c                                             at top/bottom, WC = 0 and ZMASS = 0, so it doesn't matter
c                                             what fluxt there is, just set to zero to define
       fluxt(1) = 0.d0 
       fluxt(Z$X+1) = 0.d0  

      RETURN
      END   


      subroutine LMTPPM1(IM, PX, DCP, A6, AR, AL, ILMT) 
C
C A6 =  CURVATURE OF THE TEST PARABOLA
C AR =  RIGHT EDGE VALUE OF THE TEST PARABOLA
C AL =  LEFT  EDGE VALUE OF THE TEST PARABOLA
C DCP =  0.5 * MISMATCH
C P  =  CELL-AVERAGED VALUE
C IM =  VECTOR LENGTH
C
C OPTIONS:
C
C ILMT =-1: Improved full monotonicity constraint (Lin, 2000) (IORD=2)
C ILMT = 0: FULL MONOTONICITY
C ILMT = 1: SEMI-MONOTONIC CONSTRAINT (NO UNDERSHOOTS)
C ILMT = 2: POSITIVE-DEFINITE CONSTRAINT
C                                  
       real*8 R12                  
       parameter ( R12 = 1./12. )  
       real*8 A6(IM), AR(IM), AL(IM), px(IM), DCP(IM)
       real*8 da1, da2, A6DA, fmin, dll, drr
C

      if (ILMT .eq. 0) then                             ! Full constraint
c                                          ! these constraints check out exactly as in eqn. 1.10 in C&W, 1984
      do 100 i=1,IM
      if (DCP(i) .eq. 0.) then
            AR(i) = px(i)
            AL(i) = px(i)
            A6(i) = 0.d0
      else
          da1  = AR(i) - AL(i)
          da2  = da1**2
          A6DA = A6(i)*da1
        if (A6DA .lt. -da2) then
            A6(i) = 3.*(AL(i)-px(i))
            AR(i) = AL(i) - A6(i)
        elseif (A6DA .gt. da2) then
            A6(i) = 3.*(AR(i)-px(i))
            AL(i) = AR(i) - A6(i)
        endif
      endif
100   continue


      elseif (ILMT .eq. -1) then        !!  Improved Full monotonicity constraint (Lin, 2000)
                                        !!  FORTRAN sign func: result = abs(1st)*sgn(2nd)
      do 323 i=1,im 
          da1 = px(i) + px(i)
          dll  = DSIGN( DMIN1( DABS(da1), DABS(AL(i)-px(i))), da1)
          drr  = DSIGN( DMIN1( DABS(da1), DABS(AR(i)-px(i))), da1)
 
          AR(i) = px(i) + drr
          AL(i) = px(i) - dll
          A6(i) = 3.*(dll - drr)
 323   CONTINUE


      elseif (ILMT .eq. 1) then                               ! Semi-monotonic constraint

      do 150 i=1,IM
      if (DABS(AR(i)-AL(i)) .GE. -A6(i)) go to 150
      if (px(i) .lt. AR(i) .and. px(i) .lt. AL(i)) then
            AR(i) = px(i)
            AL(i) = px(i)
            A6(i) = 0.d0
      elseif (AR(i) .gt. AL(i)) then
            A6(i) = 3.*(AL(i)-px(i))
            AR(i) = AL(i) - A6(i)
      else
            A6(i) = 3.*(AR(i)-px(i))
            AL(i) = AR(i) - A6(i)
      endif
150   continue


      elseif (ILMT .eq. 2) then                               ! Positive Definite constraint

      do 250 i=1,IM
      if (DABS(AR(i)-AL(i)) .GE. -A6(i)) go to 250
      fmin = px(i) + 0.25*(AR(i)-AL(i))**2/A6(i) + A6(i)*R12
      if (fmin .ge. 0.) go to 250
      if (px(i) .lt. AR(i) .and. px(i) .lt. AL(i)) then
            AR(i) = px(i)
            AL(i) = px(i)
            A6(i) = 0.
      elseif (AR(i) .gt. AL(i)) then
            A6(i) = 3.*(AL(i)-px(i))
            AR(i) = AL(i) - A6(i)
      else
            A6(i) = 3.*(AR(i)-px(i))
            AL(i) = AR(i) - A6(i)
      endif
250   continue

      ENDIF

      RETURN
      END


       subroutine CROSST1(IN, UA0, PX, PCR, PTT, IORD)                   ! px(0:IN+1), UA0(IN), pcr,ptt(in)

c  this computes the cross-derivative (SPLITTING) term which arise from the sequential
c  splitting process, and is critical to the stability of the scheme. It can be interpreted 
c  as the convergence of a secondary oblique flux. It arises since the physical transport 
c  can neither be perfectly modeled by the indepedent applications of the orthogonal w*
c  and v* fields, nor can it be perfectly modeled by a sequential process of w* advection
c  and then v* advection as was done by  the Prather scheme. This term vanishes for a
c  uniform mixing ratio field.  See Lin and Rood (1996) for a complete discussion of this.
c  Also, inclusion of the splitting terms are necessary for stability, ie, MAX(Cy, Cz) <1,
c  instead of Cy+Cz <1 (Cy, Cz are the Courant nums in y,z), the larger the Cy and Cz, the
c  larger the cross terms, UA0 = box averaged Courant number using WC of the UPWIND boxes
c
CEF  We DONT NEED TO adjust ANYTHING HERE for the irregular grid, SINCE THE GRID SPACING IS ALREADY 
CEF   TAKEN INTO ACCOUNT IN THE DTDZ COURANT NUMBER, 
CEF     (if there is a problem, then need to somehow normalize the terms px(iiu-1) - px(iiu) 
CEF           and px(iiu) - px(iiu+1) by the grid spacing.....)


        real*8 ua0(in), px(0:IN+1), ptt(in), pcr(in), rut


             do 120 ij=1,IN
                 iu = IDINT(ua0(ij))
                 rut = ua0(ij) - iu
                 iiu = ij - iu

             if (ua0(ij) .ge. 0.0) then
                pcr(ij) = px(iiu) + rut*(px(iiu-1) - px(iiu))
             else
                pcr(ij) = px(iiu) + rut*(px(iiu) - px(iiu+1))
             endif
 120       continue

       do 300 ij=1,IN
 300      ptt(ij) = 0.5*(pcr(ij) + px(ij))

       RETURN 
      END

