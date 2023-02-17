C
         SUBROUTINE LINTERP8(xxin, ijn, zzin, xxout, ijout, iyl, ZZOUT)
C
C   routine to do linear interpolation in one dimenstion on model input files, (EF, 8/01)
C
C    SAME AS LINTERP, except everything in REAL*8
C
C
C example CALLING sequence:
C
C        CALL LINTERP8(LAT, 37, TEMP, LATI, 73, 0, TEMPI)
C
C    
C    INPUT:
C            XXIN - 1-D array of X-axis values of dimension IJN corresponding to the dimension of ZZIN
C
C            ZZIN - 1-D array of values to be interpolated, of dimension  (IJN) 
C
C           XXOUT - 1-D array of X-axis values to be interpolated to, of dimension IJOUT
C
C           IYL - 0/1 for linear/logarithmic interpolation
C
C   OUTPUT: 
C
C           ZZOUT - 1-D array of interpolated values, of dimension IJOUT
C
C don't think we need
CEF#include "com2d.h"


       INTEGER ijn, iyl

       REAL*8 xxin(ijn), xxout(ijout), zzin(ijn), zzout(ijout)


       do 100 ij=1,ijout
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


c  interpolate, fix up edges as special cases

        if (ij0 .eq. 0) zzout(ij) = zzin(1)

        if (ij0 .ge. 1  .and.  ij0 .le. ijn) then
          if (iyl .eq. 0) then
            zzout(ij) = (zzin(ij0+1)-zzin(ij0))/
     >       (xxin(ij0+1)-xxin(ij0))*(xxout(ij)-xxin(ij0)) + zzin(ij0)
          else
            zzout(ij) = DEXP(DLOG(zzin(ij0+1)/zzin(ij0))/
     >   (xxin(ij0+1)-xxin(ij0))*(xxout(ij)-xxin(ij0))+ DLOG(zzin(ij0)))
          endif
        endif

        if (ij0 .ge. ijn+1) zzout(ij) = zzin(ijn)

 100	CONTINUE

      RETURN
      END

