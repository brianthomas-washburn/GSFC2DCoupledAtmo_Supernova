C
         SUBROUTINE BINTERP8(xxin, ijn, yyin, ikn, zzin, xxout, ijo, 
     >                       yyout, iko, iyl, izl, ZZOUT)
C
C   routine to do bilinear interpolation (latitude-height) on model input files, (EF, 8/01)
C
C    SAME AS BINTERP, except everything in REAL*8
C
C
C example CALLING sequence:
C
C        CALL BINTERP8(LAT, 37, ZZ, 59, TEMP, LATI, 73, ZZI, 150, 1, 1, TEMPI)
C
C    
C    INPUT:
C            XXIN - 1-D array of X-axis values of dimension IJN (usually latitudes) corresponding to 
C                   the 1st dimension of ZZIN
C
C            YYIN - 1-D array of Y-axis values of dimension IKN (usually Pressure levels) corresponding to 
C                   the 2nd dimension of ZZIN
C
C            ZZIN - 2-D array of values to be interpolated, of dimension  (IJN,IKN) 
C
C           XXOUT - 1-D array of X-axis values (usually latitudes) to be interpolated to, of dimension IJO
C
C          YYOUT - 1-D array of Y-axis values (usually pressure levs) to be interpolated to, of dimension IKO
C
C           IYL - 0/1 for linear/logarithmic interpolation in latitude 
C
C           IZL - 0/1 for linear/logarithmic interpolation in altitude 
C
C   OUTPUT: 
C
C           ZZOUT - 2-D array of interpolated values, of dimension (IJO,IKO)
C
C  
C   NOTE: If you are doing either linear or log interpolation in both directions, the order of 
C         interpolation DOES NOT MATTER, i.e., latitude then height vs. height then latitude
C         HOWEVER, if you are mixing them, e.g., linear in latitude and logarithmic in height, 
C         the order DOES make a small difference - on the order of 1-2% in extreme cases, 
C         probably less in more normal cases.
C
C
C don't think we need
CEF#include "com2d.h"


       INTEGER ijn, ikn, ijo, iko, iyl, izl

       REAL*8 xxin(ijn), yyin(ikn), zzin(ijn,ikn), q1, q2
       REAL*8 xxout(ijo), yyout(iko), zzout(ijo,iko)



       do 100 ik=1,IKO
       do 100 ij=1,IJO

C                                            find latitude indicies of input file
          ij0=1  
          if (xxout(ij) .lt. xxin(1)) ij0=0
 500      if (ij0 .lt. ijn) then
            if (xxout(ij) .gt. xxin(ij0+1)) then
               ij0 = ij0 + 1
               go to 500
            endif 
          endif 
          if (xxout(ij) .gt. xxin(ijn)) ij0=ijn+1

C                                            find altitude indicies of input file
          ik0=1  
          if (yyout(ik) .lt. yyin(1)) ik0=0
 600      if (ik0 .lt. ikn) then
            if (yyout(ik) .gt. yyin(ik0+1)) then
               ik0 = ik0 + 1
               go to 600
            endif 
          endif 
          if (yyout(ik) .gt. yyin(ikn)) ik0=ikn+1

c  interpolate in x-axis direction, then y-direction, fix up edges and corners as special cases

      IF (ik0 .eq. 0) THEN
        if (ij0 .eq. 0) zzout(ij,ik) = zzin(1,1)

        if (ij0 .ge. 1  .and.  ij0 .le. ijn) then
          if (iyl .eq. 0) then
            zzout(ij,ik) = (zzin(ij0+1,1)-zzin(ij0,1))/
     >       (xxin(ij0+1)-xxin(ij0))*(xxout(ij)-xxin(ij0)) + zzin(ij0,1)
          else
            zzout(ij,ik) = DEXP(DLOG(zzin(ij0+1,1)/zzin(ij0,1))/
     > (xxin(ij0+1)-xxin(ij0))*(xxout(ij)-xxin(ij0))+ DLOG(zzin(ij0,1)))
          endif
        endif

        if (ij0 .ge. ijn+1) zzout(ij,ik) = zzin(ijn,1)
      ENDIF


      IF (ik0 .ge. ikn+1) THEN
        if (ij0 .eq. 0) zzout(ij,ik) = zzin(1,ikn)

        if (ij0 .ge. 1  .and.  ij0 .le. ijn) then
          if (iyl .eq. 0) then
           zzout(ij,ik) = (zzin(ij0+1,ikn)-zzin(ij0,ikn))/
     >     (xxin(ij0+1)-xxin(ij0))*(xxout(ij)-xxin(ij0)) + zzin(ij0,ikn)
          else
           zzout(ij,ik) = DEXP(DLOG(zzin(ij0+1,ikn)/zzin(ij0,ikn))/
     >(xxin(ij0+1)-xxin(ij0))*(xxout(ij)-xxin(ij0))+DLOG(zzin(ij0,ikn)))
          endif
        endif

        if (ij0 .ge. ijn+1) zzout(ij,ik) = zzin(ijn,ikn)
      ENDIF


      IF (ik0 .ge. 1  .and.  ik0 .le. ikn) THEN
        if (ij0 .eq. 0) then 
           q1 = zzin(1,ik0)
           q2 = zzin(1,ik0+1)
        endif

        if (ij0 .ge. 1  .and.  ij0 .le. ijn) then
          if (iyl .eq. 0) then
            q1 = (zzin(ij0+1,ik0)-zzin(ij0,ik0))/(xxin(ij0+1)-xxin(ij0))
     >            *(xxout(ij)-xxin(ij0)) + zzin(ij0,ik0)

          q2=(zzin(ij0+1,ik0+1)-zzin(ij0,ik0+1))/(xxin(ij0+1)-xxin(ij0))
     >            *(xxout(ij)-xxin(ij0)) + zzin(ij0,ik0+1)
          ELSE
          q1 = DEXP(DLOG(zzin(ij0+1,ik0)/zzin(ij0,ik0))/(xxin(ij0+1)-
     >           xxin(ij0))*(xxout(ij)-xxin(ij0)) + DLOG(zzin(ij0,ik0)))

          q2= DEXP(DLOG(zzin(ij0+1,ik0+1)/zzin(ij0,ik0+1))/(xxin(ij0+1)-
     >         xxin(ij0))*(xxout(ij)-xxin(ij0)) + DLOG(zzin(ij0,ik0+1)))
          endif
        endif

        if (ij0 .ge. ijn+1) then 
           q1 = zzin(ijn,ik0)
           q2 = zzin(ijn,ik0+1)
        endif

          if (izl .eq. 0) then
             zzout(ij,ik) = (q2-q1)/(yyin(ik0+1)-yyin(ik0))*
     >              (yyout(ik)-yyin(ik0)) + q1
          else
             zzout(ij,ik) = DEXP(DLOG(q2/q1)/(yyin(ik0+1)-yyin(ik0))*
     >            (yyout(ik)-yyin(ik0)) + DLOG(q1))
          endif
      ENDIF

 100	CONTINUE

      RETURN
      END

