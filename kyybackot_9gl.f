C
        SUBROUTINE KYYBACK(latyp, zp)

C
c  file KYYBACK_9GL.f - FORTRAN routine TO  PRESCRIBE  BACKGROUND Kyys for coupled model
C                       similar to code used in FIXED MODEL
C
C                                   include PARAM.INC here for N$, M$
        INCLUDE 'PARAM.INC'


        !INTEGER N$, M$, NP$, MP$
        REAL LATYP(N$), ZP(M$), xkot


        COMMON/CKYY1/ BGKYY0(N$,M$), BGKYY(N$,M$), EPOTX(NP$,MP$), 
     >                XKYYFX(N$,M$)


C
C  BGKYY(N$,M$) is in cm2/sec ;  first set MINIMUM KYY to 2.e8 EVERYWHERE
C
        do 7756 ik=1,M$
        do 7756 ij=1,N$
 7756      bgkyy(ij,ik) = 2.e8


C
C  Specify minimum Kyy below 50 km: this is based on fixed model clim Kyy
C

      do 730 ik=1,M$

         if (zp(ik) .le. 50.) then

      do 732 ij=1,N$
C
c                                              first set tropics/midlats - go up to 50N,S:
       IF (ABS(latyp(ij)) .le. 50.) then

        if (zp(ik) .ge. 40.7) bgkyy(ij,ik) = 20.e8

        if (zp(ik) .ge. 31.7 .and.  zp(ik) .lt. 40.7) bgkyy(ij,ik)=20.e8

        if (zp(ik) .ge. 30.7 .and.  zp(ik) .lt. 31.7) then
             bgkyy(ij,ik)= 20.e8
             if (ABS(latyp(ij)) .le. 15.) bgkyy(ij,ik) = 20.e8
             if (ABS(latyp(ij)) .le. 10.) bgkyy(ij,ik) = 15.e8
             if (ABS(latyp(ij)) .le.  6.) bgkyy(ij,ik) = 10.e8
        endif

        if (zp(ik) .ge. 29.7 .and.  zp(ik) .lt. 30.7) then
             bgkyy(ij,ik)= 15.e8
             if (ABS(latyp(ij)) .le. 15.) bgkyy(ij,ik) = 15.e8
             if (ABS(latyp(ij)) .le. 10.) bgkyy(ij,ik) = 10.e8
             if (ABS(latyp(ij)) .le.  6.) bgkyy(ij,ik) =  7.e8
        endif

        if (zp(ik) .ge. 27.7 .and.  zp(ik) .lt. 29.7) then
             bgkyy(ij,ik)= 10.e8
             if (ABS(latyp(ij)) .le. 15.) bgkyy(ij,ik) = 10.e8
             if (ABS(latyp(ij)) .le. 10.) bgkyy(ij,ik) =  7.e8
        endif

        if (zp(ik) .ge. 26.7 .and.  zp(ik) .lt. 27.7) then
             bgkyy(ij,ik)=  7.e8
             if (ABS(latyp(ij)) .le. 15.) bgkyy(ij,ik) =  7.e8
             if (ABS(latyp(ij)) .le. 10.) bgkyy(ij,ik) =  5.e8
             if (ABS(latyp(ij)) .le.  6.) bgkyy(ij,ik) =  4.e8
        endif

        if (zp(ik) .ge. 20.7 .and.  zp(ik) .lt. 26.7) then
             bgkyy(ij,ik)=  5.e8
             if (ABS(latyp(ij)) .le. 15.) bgkyy(ij,ik) = 4.e8
             if (ABS(latyp(ij)) .le. 10.) bgkyy(ij,ik) = 3.e8
             if (ABS(latyp(ij)) .le.  6.) bgkyy(ij,ik) = 3.e8
        endif

        if (zp(ik) .ge. 19.7 .and.  zp(ik) .lt. 20.7) then
             bgkyy(ij,ik)=  7.e8
             if (ABS(latyp(ij)) .le. 15.) bgkyy(ij,ik) =  5.e8
             if (ABS(latyp(ij)) .le. 10.) bgkyy(ij,ik) =  3.e8
             if (ABS(latyp(ij)) .le.  6.) bgkyy(ij,ik) =  3.e8
        endif

        if (zp(ik) .ge. 18.7 .and.  zp(ik) .lt. 19.7) then
             bgkyy(ij,ik)=  7.e8
             if (ABS(latyp(ij)) .le. 15.) bgkyy(ij,ik) =  5.e8
             if (ABS(latyp(ij)) .le. 10.) bgkyy(ij,ik) =  4.e8
             if (ABS(latyp(ij)) .le.  6.) bgkyy(ij,ik) =  4.e8
        endif

        if (zp(ik) .ge. 17.7 .and.  zp(ik) .lt. 18.7) then
             bgkyy(ij,ik)=  7.e8
             if (ABS(latyp(ij)) .le. 15.) bgkyy(ij,ik) =  7.e8
             if (ABS(latyp(ij)) .le. 10.) bgkyy(ij,ik) =  5.e8
             if (ABS(latyp(ij)) .le.  6.) bgkyy(ij,ik) =  5.e8
        endif

        if (zp(ik) .ge. 16.7 .and. zp(ik) .lt. 17.7) then
             bgkyy(ij,ik)= 10.e8
             if (ABS(latyp(ij)) .le. 15.) bgkyy(ij,ik) = 10.e8
        endif

        if (zp(ik) .ge. 15.7 .and. zp(ik) .lt. 16.7) then
             bgkyy(ij,ik)= 20.e8
             if (ABS(latyp(ij)) .le. 15.) bgkyy(ij,ik) = 20.e8
        endif

        if (zp(ik) .ge. 13.7 .and. zp(ik) .lt. 15.7) then
             bgkyy(ij,ik)= 20.e8
             if (ABS(latyp(ij)) .le. 30.) bgkyy(ij,ik) = 30.e8
        endif

        if (zp(ik) .ge. 9.7 .and. zp(ik) .lt. 13.7) then
             bgkyy(ij,ik)= 20.e8
             if (ABS(latyp(ij)) .le. 30.) bgkyy(ij,ik) = 30.e8
             if (ABS(latyp(ij)) .le. 25.) bgkyy(ij,ik) = 40.e8
        endif

        if (zp(ik) .ge. 7.7 .and. zp(ik) .lt. 9.7) then
             bgkyy(ij,ik)= 30.e8
             if (ABS(latyp(ij)) .le. 30.) bgkyy(ij,ik) = 40.e8
             if (ABS(latyp(ij)) .le. 25.) bgkyy(ij,ik) = 50.e8
        endif

        if (zp(ik) .ge. 6.7 .and.  zp(ik) .lt. 7.7) bgkyy(ij,ik)= 75.e8

        if (zp(ik) .ge. 5.7  .and.  zp(ik) .lt. 6.7) bgkyy(ij,ik)=100.e8

        if (zp(ik) .ge. 4.7  .and.  zp(ik) .lt. 5.7) bgkyy(ij,ik)=150.e8

        if (zp(ik) .lt. 4.7) bgkyy(ij,ik) = 200.e8

      ENDIF
                              ! end tropics/subtropics loop

ccmid 
ccmid   -  MIDLATITUDE LOOP - just use as transition zone, 51-55N,S
 
      IF (ABS(latyp(ij)) .gt. 50.  .and. ABS(latyp(ij)) .le. 55.) then

        if (zp(ik) .ge. 40.7) bgkyy(ij,ik)= 15.e8
 
        if (zp(ik) .ge. 29.7 .and. zp(ik) .lt. 40.7) bgkyy(ij,ik)= 15.e8

        if (zp(ik) .ge. 17.7 .and. zp(ik) .lt. 29.7) bgkyy(ij,ik)= 10.e8

        if (zp(ik) .ge. 13.7 .and. zp(ik) .lt. 17.7) bgkyy(ij,ik)= 15.e8

        if (zp(ik) .ge. 11.7 .and. zp(ik) .lt. 13.7) bgkyy(ij,ik)= 20.e8

        if (zp(ik) .ge. 9.7 .and.  zp(ik) .lt. 11.7) bgkyy(ij,ik)= 20.e8

        if (zp(ik) .ge. 7.7  .and. zp(ik) .lt.  9.7) bgkyy(ij,ik)= 30.e8

        if (zp(ik) .ge. 6.7 .and.  zp(ik) .lt. 7.7) bgkyy(ij,ik)= 75.e8

        if (zp(ik) .ge. 5.7  .and.  zp(ik) .lt. 6.7) bgkyy(ij,ik)=100.e8

        if (zp(ik) .ge. 4.7  .and.  zp(ik) .lt. 5.7) bgkyy(ij,ik)=150.e8

        if (zp(ik) .lt. 4.7) bgkyy(ij,ik) = 200.e8
 
       ENDIF
                                ! end mid-latitude loop


C  for polar region, ramp down Kyy faster w/ height than at midlatitudes to avoid blowing up the vortex. 
C
       IF (ABS(latyp(ij)) .gt. 55.) then

        if (zp(ik) .ge. 40.7) bgkyy(ij,ik) = 15.e8

        if (zp(ik) .ge. 35.7 .and. zp(ik) .lt. 40.7) bgkyy(ij,ik)= 10.e8

        if (zp(ik) .ge. 15.7 .and. zp(ik) .lt. 35.7) bgkyy(ij,ik)= 7.5e8

        if (zp(ik) .ge. 13.7 .and. zp(ik) .lt. 15.7) bgkyy(ij,ik)= 10.e8

        if (zp(ik) .ge. 11.7 .and. zp(ik) .lt. 13.7) bgkyy(ij,ik)= 20.e8

        if (zp(ik) .ge. 9.7 .and.  zp(ik) .lt. 11.7) bgkyy(ij,ik)= 20.e8

        if (zp(ik) .ge. 7.7  .and.  zp(ik) .lt. 9.7) bgkyy(ij,ik)= 30.e8

        if (zp(ik) .ge. 6.7 .and.  zp(ik) .lt. 7.7)  bgkyy(ij,ik)= 75.e8

        if (zp(ik) .ge. 5.7 .and. zp(ik) .lt. 6.7) bgkyy(ij,ik)= 100.e8

        if (zp(ik) .ge. 4.7  .and.  zp(ik) .lt. 5.7) bgkyy(ij,ik)=150.e8

        if (zp(ik) .lt. 4.7) bgkyy(ij,ik) = 200.e8

      ENDIF
                             ! end high latitude loop

 732       CONTINUE

        endif      

c                             ! ramp up MIN Kyy=20.e8  AT ALL LATITUDES above 50  km


       if (zp(ik) .gt. 50.) then
         do 744 ij=1,N$
 744        bgkyy(ij,ik) = 20.e8
       endif


 730   CONTINUE


C                                  ! store off these original values
        do 5677 ik=1,M$
        do 5677 ij=1,N$
 5677      bgkyy0(ij,ik) = bgkyy(ij,ik)


	RETURN
	END



C
        SUBROUTINE SKYYOT(latyp, zp, IDOY360, QQY)
c
C
c     this loads in Kyyot, and adjusts EP Flx divergence accordingly for coupled model
C                                  NOTE: EP FLx div correction DOES NOT WORK!!!!! get screwy easterlies

        INCLUDE 'PARAM.INC'


        !INTEGER N$, M$, NP$, MP$, IDOY360
        INTEGER IDOY360
        REAL LATYP(N$), ZP(M$), QQY(N$,M$), EPOT(N$,M$), xkot
        REAL XKYYIN(36,46), xkday(N$,M$)
C                                            common for KYYOT on coupled model dynamics grid (from TEMPIN)
C                                                         XKYYOTC, BGKYY are in cm2/sec,  EPOTX is in cm/sec2
        COMMON/CKYYOT/ XKYYOTC(36,46,360), ypot(36), zpot(46)

        COMMON/CKYY1/ BGKYY0(N$,M$), BGKYY(N$,M$), EPOTX(NP$,MP$),
     >                XKYYFX(N$,M$)


C                                ! load in orginal BGKYY, and initialize EP flux adjustment
        do 700 ik=1,M$
        do 700 ij=1,N$
           bgkyy(ij,ik) = bgkyy0(ij,ik)
 700       epot(ij,ik) = 0.0



C   interpolate KYYOTs to current coupled model dynamics grid for current day:
C                       XKYYOTC(36,46,360) => xkday(N$,M$),  latyp(N$), zp(M$)
C
        do 705 ik=1,46
        do 705 ij=1,36
 705       xkyyin(ij,ik) = xkyyotc(ij,ik,idoy360)

        CALL BINTERP(ypot, 36, zpot, 46, XKYYIN,
     >               latyp, N$, zp, M$, 0, 0, xkday)



c                                                 use KYYOTC for 50S-50N, 7-19 km
      do 1730 ik=1,M$

      if (zp(ik) .ge. 7.  .and.  zp(ik) .le. 19.) then

        do 1732 ij=1,N$


       IF (ABS(latyp(ij)) .gt. 20.  .and.  ABS(latyp(ij)) .le. 50.) then
C                                                                ! blend over several levels
           xkot = xkday(ij,ik)

           if (zp(ik) .ge. 15.7)xkot = .9*bgkyy(ij,ik) + .1*xkday(ij,ik)
   
           if (zp(ik) .ge. 12.7 .and. zp(ik) .lt. 15.7)
     >            xkot = .8*bgkyy(ij,ik) + .2*xkday(ij,ik)

           if (zp(ik) .ge. 10.7 .and. zp(ik) .lt. 12.7) 
     >            xkot = .7*bgkyy(ij,ik) + .3*xkday(ij,ik)

           if (zp(ik) .ge. 9.7 .and. zp(ik) .lt. 10.7) 
     >            xkot = .6*bgkyy(ij,ik) + .4*xkday(ij,ik)

ccc           if (zp(ik) .lt. 14.7) xkot = xkday(ij,ik)

                if (bgkyy(ij,ik) .lt. xkot) then

                      bgkyy(ij,ik) = xkot
                                                                    ! EPOT is in cm/sec2
ccelf                      epot(ij,ik) = -bgkyy(ij,ik)*qqy(ij,ik)
ccelf                      if (epot(ij,ik) .ge. 0.) epot(ij,ik) = 0.0
                 endif
         ENDIF


         IF (ABS(latyp(ij)) .le. 20.) then
C                                                                ! blend over several levels
           xkot = xkday(ij,ik)

           if (zp(ik) .ge. 17.7)xkot = .9*bgkyy(ij,ik) + .1*xkday(ij,ik)
   
           if (zp(ik) .ge. 16.7 .and. zp(ik) .lt. 17.7) 
     >            xkot = .7*bgkyy(ij,ik) + .3*xkday(ij,ik)

           if (zp(ik) .ge. 15.7 .and. zp(ik) .lt. 16.7) 
     >            xkot = .5*bgkyy(ij,ik) + .5*xkday(ij,ik)

           if (zp(ik) .ge. 14.7 .and. zp(ik) .lt. 15.7) 
     >            xkot = .25*bgkyy(ij,ik) + .75*xkday(ij,ik)

ccc           if (zp(ik) .lt. 14.7) xkot = xkday(ij,ik)

                if (bgkyy(ij,ik) .lt. xkot) then

                      bgkyy(ij,ik) = xkot
                                                                    ! EPOT is in cm/sec2
ccelf                      epot(ij,ik) = -bgkyy(ij,ik)*qqy(ij,ik)
ccelf                      if (epot(ij,ik) .ge. 0.) epot(ij,ik) = 0.0
                 endif

         ENDIF


 1732      CONTINUE

        ENDIF      

 1730   CONTINUE

                                        ! convert EPOT to proper grid
        CALL MRS2RBR( EPOT , EPOTX )


ccccccc        write (755) N$, M$, NP$, MP$, IDOY360
ccccccc        write (755) latyp, zp
ccccccc        write (755) QQY
ccccccc        write (755) BGKYY     
ccccccc        write (755) EPOT
ccccccc        write (755) EPOTX


	RETURN
	END
