C
        SUBROUTINE KZZADJ(IIDAY360, TIMEKZ, xlatz, zaltz, KZZ0, 
     >                    KZZTH, KZZDAY)

!        USE degree_trig

C
C   THIS ROUTINE adjusts current day's tropospheric Kzz for fixed and coupled models
C        
C
C   common for the vertical gradients of the DU/km ozone climatology, for Kzz adjustment
C
        COMMON/CDUKZ/dukz(14,45,89)


        INTEGER IIDAY360

        REAL*4 KZZ0(45,89), KZZDAY(45,89), KZZTH(45,89), xlatz(45)
        REAL*4 ttime(1), ttout(1), g4ss(14), timekz(14), zaltz(89)


        ttime(1) = IIDAY360*1.


C
C  use Kzz from fixed model for current day : KZZ0(45,89)
C   modify troposphere, also used for coupled model dynamics
C   set max of 30.e4 in lower trop, 100.e4 elsewhere  -  this seems to work well
C   set minimum of .001*1.e4 for DYNAMICS - this also done in XCOUPLED for CHEMISTRY
C   first find indicies closest to 5 km and 16 km of current dynamics grid, KZZ0(45,89)
C
C  WCONV268: KZZTH(45,89) is used for coupled model dynamics, xlatz(45), zaltz(89)
C     KZZDAY(45,89) is used for the chemistry (July 2011)
C

        i5z = 6
        i12z = 12
        i15z = 16
        i16z = 16

        do 877 ik=1,89
           if (zaltz(ik) .lt. 5.2)   i5z = ik
           if (zaltz(ik) .lt. 12.2) i12z = ik
           if (zaltz(ik) .lt. 15.2) i15z = ik
           if (zaltz(ik) .lt. 16.2) i16z = ik
 877    CONTINUE



      DO 8517 IK=1,89
      DO 8517 IJ=1,45

         xxlt = ABS(xlatz(ij))
         zzk  = zaltz(ik)


         xkzmx = 100.e4
         if (zzk .le. 4.5) xkzmx = 10.e4

CWCONV200 - enhance trop Kzz at 50S-50N, largest at 20S-20N, theta Kzz done separately (TRCRLIB)

         if (xxlt .le. 60.  .and.  xxlt .gt. 50.) then
            if (zzk .gt. 4.5  .and.  zzk .le. 15.2) then 
              xkzz5 = 10.e4
              if (KZZ0(ij,i5z) .lt. 10.e4) xkzz5 = KZZ0(ij,i5z)

         xkzmx = EXP( ALOG(KZZ0(ij,i15z)/xkzz5)/(zaltz(i15z)-zaltz(i5z))
     >             *(zzk-zaltz(i5z)) + ALOG(xkzz5) )
            endif
         endif

         if (xxlt .le. 50.) then
            xxkz = 10.e4 + 20.e4*COSD((xxlt-20.)*3.)
            if (xxlt .le. 20.) xxkz = 30.e4

            if (zzk .le. 4.5) xkzmx = xxkz

            if (zzk .gt. 4.5  .and.  zzk .le. 15.2) then 
              xkzz5 = xxkz
              if (KZZ0(ij,i5z) .lt. xxkz) xkzz5 = KZZ0(ij,i5z)

         xkzmx = EXP( ALOG(KZZ0(ij,i15z)/xkzz5)/(zaltz(i15z)-zaltz(i5z))
     >             *(zzk-zaltz(i5z)) + ALOG(xkzz5) )
            endif
         endif

cckzz
cckzz
cckzz    for testing, set larger Kzz below 13 km.....
cckzz
cckzz         if (zzk .le. 13.) xkzmx = 10.e4
cckzz
cckzz   for testing, define Kzz based on latent heating - didn't seem to work as well
cckzz      this needs more work if it is to be used......
cckzz
cckzz         if (zzk .le. 11.) then
cckzz              xkzzlh = ( WLHDAY(ij+ij,ik) + HACKDAY(ij+ij,ik) )*7.*1.e4
cckzz              if (xkzzlh .le. .15*1.e4) xkzzlh = .15*1.e4
cckzz              KZZ0(ij,ik) = xkzzlh
cckzz         endif
cckzz

         if (KZZ0(ij,ik) .ge. xkzmx) KZZ0(ij,ik) = xkzmx

         xkzmn = .001*1.e4
         if (KZZ0(ij,ik) .le. xkzmn) KZZ0(ij,ik) = xkzmn


C  set tropospheric values here:  Minima for 0-7 km;  High lat Maxima above 8 km

         if (zzk .le. 5.1) then
            if (KZZ0(ij,ik) .le. 3.e4) KZZ0(ij,ik) = 3.e4
         endif

         if (zzk .gt. 5.1  .and.  zzk .le. 6.1) then
            if (KZZ0(ij,ik) .le. 2.e4) KZZ0(ij,ik) = 2.e4
         endif

         if (zzk .gt. 6.1  .and.  zzk .le. 7.1) then
            if (KZZ0(ij,ik) .le. 1.e4) KZZ0(ij,ik) = 1.e4
         endif


         if (zzk .gt. 7.1  .and.  zzk .le. 8.1) then
            if (xxlt .ge. 60.  .and.  xxlt .lt. 65.) then
              if (KZZ0(ij,ik) .ge. 1.e4) KZZ0(ij,ik) = 1.e4
            endif

            if (xxlt .ge. 65.) then
              if (KZZ0(ij,ik) .ge. .5e4) KZZ0(ij,ik) = .5e4
            endif
         endif

         if (zzk .gt. 8.1  .and.  zzk .le. 9.1) then
            if (xxlt .ge. 60.  .and.  xxlt .lt. 65.) then
              if (KZZ0(ij,ik) .ge. 1.e4) KZZ0(ij,ik) = 1.e4
            endif

            if (xxlt .ge. 65.) then
              if (KZZ0(ij,ik) .ge. .4e4) KZZ0(ij,ik) = .4e4
            endif
         endif

         if (zzk .gt. 9.1  .and.  zzk .le. 10.1) then
            if (xxlt .ge. 60.  .and.  xxlt .lt. 65.) then
              if (KZZ0(ij,ik) .ge. .7e4) KZZ0(ij,ik) = .7e4
            endif

            if (xxlt .ge. 65.) then
              if (KZZ0(ij,ik) .ge. .4e4) KZZ0(ij,ik) = .4e4
            endif
         endif

         if (zzk .gt. 10.1  .and.  zzk .le. 14.1) then
            if (xxlt .ge. 60.) then
              if (KZZ0(ij,ik) .ge. .4e4) KZZ0(ij,ik) = .4e4
            endif
         endif


C  reduce SUBTROP/MIDLAT Kzz in UT/LS-Midstrat for chemistry to get better AGE/Cly
C    this was the Extra turbulent Kzz added into the fixed model
C    This will SLIGHTLY effect Kzz for UBAR - Kzz for theta gets overwritten in TRCRLIB
C    xlatz(45), zaltz(89); kzz0(45,89)

CCcc           xkzzeq = KZZ0(23,ik)
cccc             if (KZZ0(ij,ik) .ge. xkzzeq) KZZ0(ij,ik) = xkzzeq

          if (xxlt .ge. 20.) then
             if (zzk .ge. 15.  .and.  zzk .le. 25.) then

c   first define minima (before KZZ0 is reset), based on polar values:

                if (xlatz(ij) .lt. 0.) ijj = 1
                if (xlatz(ij) .ge. 0.) ijj = 45
                xkzn = KZZ0(ijj,ik)*(SIND(xlatz(ij))**2)


                xkzz0 = KZZ0(ij,ik)

                if (zzk .le. 16.  .or.  zzk .ge. 24.) then
                    KZZ0(ij,ik) = .7*xkzz0
                else
                    KZZ0(ij,ik) = .3*xkzz0 
                endif

                if (KZZ0(ij,ik) .le. xkzn) KZZ0(ij,ik) = xkzn
                if (KZZ0(ij,ik) .le. .01e4) KZZ0(ij,ik) = .01e4
             endif
          endif


C   if Kzzeq LE .01 (ie, above ~19-20 km), set extratropical min

cccc             if (xkzzeq .le. .01e4) then
cccc      xkzn = 2.*xkzzeq + .01e4 - (xkzzeq + .01e4)*(COSD(xlatz(ij))**2)
cccc
cccc                  if (xkzn .le. .001e4) xkzn = .001e4
cccc                  if (KZZ0(ij,ik) .le. xkzn) KZZ0(ij,ik) = xkzn
cccc              endif


C
C  WCONV255 - adjust Kzz in UT/LS based on du/km ozone clim vertical gradients
C        dukzd is for current day
C
C   define maximum Kzz based on dukzd, for 7-16 km, don't do 60-90N,S for WCONV256
C
C   first interpolate to current day (already on Kzz grid) -> dukzd
C      dukz(14,45,89) -> dukzd ;     g4ss(14), timekz(14)
C

      if (zzk .ge. 7.  .and.  zzk .le. 16.) then

            do 4575 iii=1,14
 4575          g4ss(iii) = dukz(iii,ij,ik)

               CALL LINTERP(timekz, 14, g4ss, ttime, 1, 0, ttout)

               dukzd = ttout(1)


        if (xxlt .le. 60) then
          if (dukzd .ge. 0.  .and.  dukzd .le. 2.) then
            xkzs = 1.e4 
            if (xxlt .le. 35) xkzs = 2.e4 
            if (xxlt .le. 30) xkzs = 3.e4 
            if (xxlt .le. 25) xkzs = 5.e4 
            if (xxlt .le. 20) xkzs = 7.e4 

            xkzmax = EXP( ALOG(.3e4/xkzs)/(2.-0.)*dukzd + ALOG(xkzs) )


C  WCONV257 - don't do during SH Summer/Fall (Jan 1 - May 1), 30S-60S, ramp up/down

            if (xlatz(ij) .le. -30.  .and.  IIDAY360 .le. 120) then
                xkzmax0 = xkzmax
                xzf = 0.

                if (IIDAY360 .le. 30.) xzf = COSD(3.*IIDAY360)**2
                if (IIDAY360 .ge. 90.) xzf = SIND(3.*IIDAY360 + 90.)**2

                xkzmax = xzf*xkzmax0 + (1.-xzf)*KZZ0(ij,ik)
            endif


            if (KZZ0(ij,ik) .gt. xkzmax) KZZ0(ij,ik) = xkzmax
          endif
        endif

      endif

C  WCONV268 - store Kzz for dynamics calculations

         KZZTH(IJ,IK) = KZZ0(IJ,IK)



C   WCONV271 - further adjust tropical tropospheric Kzz for chemstry ONLY
C       set min of 30.e4 below 9 km, ramp down, interpolate over 12-16 km
C
         if (zzk .le. 9.1  .and.  xxlt .le. 20.) then
             if (KZZ0(ij,ik) .le. 30.e4) KZZ0(ij,ik) = 30.e4
         endif

         if (zzk .gt. 9.1  .and.  zzk .le. 10.1) then
       if (xxlt .le. 20. .and. KZZ0(ij,ik) .le. 20.e4) KZZ0(ij,ik)=20.e4
         endif

         if (zzk .gt. 10.1  .and.  zzk .le. 11.1) then
       if (xxlt .le. 20. .and. KZZ0(ij,ik) .le. 10.e4) KZZ0(ij,ik)=10.e4
         endif

         if (zzk .gt. 11.1  .and.  zzk .le. 12.1) then
       if (xxlt .le. 20. .and. KZZ0(ij,ik) .le. 5.e4) KZZ0(ij,ik) = 5.e4
         endif

         if (zzk .gt. 12.1  .and.  zzk .le. 15.1) then
            if (xlatz(ij) .gt. -30.  .and.  xlatz(ij) .lt. 33.) then
         KZZ0(ij,ik) = EXP( ALOG(KZZ0(ij,i16z)/KZZ0(ij,i12z))/
     > (zaltz(i16z)-zaltz(i12z))*(zzk-zaltz(i12z))+alog(KZZ0(ij,i12z)) )
            endif
         endif


         KZZDAY(IJ,IK) = KZZ0(IJ,IK)
 8517  CONTINUE


        RETURN
	END

        function cosd(dgr_argument)
        real(4) cosd
        real(4) dgr_argument
        real(16), parameter ::
     &     quadpi = 3.141592653589793238462643383279502884197Q0
        real(16), parameter :: dgr_to_rad = (quadpi/180Q0)
        cosd = cos(dgr_to_rad * dgr_argument)
        end function

        function sind(dgr_argument)
        real(4) sind
        real(4) dgr_argument
        real(16), parameter ::
     &     quadpi = 3.141592653589793238462643383279502884197Q0
        real(16), parameter :: dgr_to_rad = (quadpi/180Q0)
        sind = sin(dgr_to_rad * dgr_argument)
        end function

