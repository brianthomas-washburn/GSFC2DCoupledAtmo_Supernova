
       SUBROUTINE PDFINT(JJ0, J181)

c
c  This routine interpolates Js LINEARLY to 1 degree latitude grid -
C     LAT(L$) is REAL*8 in COMMON -> lat181(181)
c      JX(PH$,L$,Z$) -->  J181(181), use ijind(181) which is 1->L$-1 EVERYWHERE
C      so this ensures that you're using a constant grad for latitudes outside of the original range
C      ie, to lats > 85N, it uses the 75N-85N gradient, also check for location of terminator using DECD
C
C      NOTE:  dellat is always positive, dth is defined in COMMON, eg, dth = 10 deg, 4 deg,..
C
!       USE degree_trig

       include "com2d.h"


       REAL*8 JJ0(L$), J181(181), JJ1, JJ2, dellat


       SAVE

 
C  find latitude of terminators
   
        laters = -90. + decd
        latern = 90. + decd


C  sometimes there are J-values just beyond the terminator, 
C     ie, CHID slightly greater than 90 degrees (~90.5), so adjust laters, latern
C
C   SH:  find highest southern latitude w/ J values, adjust laters:  Also,
C        ensure that laters minus any lat(L$) is not 0, otherwise divide by 0 below

        do 115 ij=1,L$
            if (lat(ij) .lt. laters .and. JJ0(ij) .gt. 0.d0)
     >          laters = lat(ij) - .5

            xdlat = ABS(laters - lat(ij) )
            if (xdlat .lt. .1) laters = lat(ij) - .5
 115    CONTINUE


C  NH:  find highest northern latitude w/ J values, adjust latern:     Also,
C       ensure that latern minus any lat(L$) is not 0, otherwise divide by 0 below

        do 116 ij=1,L$
            if (lat(ij) .gt. latern .and. JJ0(ij) .gt. 0.d0) 
     >          latern = lat(ij) + .5

            xdlat = ABS(latern - lat(ij) )
            if (xdlat .lt. .1) latern = lat(ij) + .5
 116    CONTINUE

              

        do 250 ij=1,181
C                             for latitudes outside 1-L$ range, just set to 1 and L$ - DON'T USE
C
c              if (LAT181(ij) .le. LAT(1)) then
c                   J181(ij) = JJ0(1)
c                   go to 250
c              endif
c
c              if (LAT181(ij) .ge. LAT(L$)) then
c                   J181(ij) = JJ0(L$)
c                   go to 250
c              endif


           JJ1 = JJ0(ijind(ij))
           JJ2 = JJ0(ijind(ij)+1)

           dellat = lat181(ij) - lat(ijind(ij))
           J181(ij) = (JJ2 - JJ1)/dth*dellat + JJ1


C  Now check for terminator, if so, then RE-INTERPOLATE from terminator, and set latitudes poleward = 0.0
C
C  SH:
           if (lat(ijind(ij)) .lt. laters) then
             dellat = lat181(ij) - laters
             J181(ij) = (JJ2 - 0.d0)/(lat(ijind(ij)+1)-laters)*dellat     ! + 0.d0 is implied here
           endif
          
           if (lat181(ij) .lt. laters) J181(ij) = 0.d0

C  NH:
           if (lat(ijind(ij)+1) .gt. latern) then
             dellat = lat181(ij) - lat(ijind(ij))
             J181(ij) = (0.d0 - JJ1)/(latern-lat(ijind(ij)))*dellat +JJ1
           endif
          
           if (lat181(ij) .gt. latern) J181(ij) = 0.d0
 250   CONTINUE



C   if interpolated  J is NEGATIVE beyond boundary, ie < LAT(1) or > LAT(L$), then the falloff is too rapid, 
C   so just reduce by COS(lat), use COS(89.9) for poles,  ie, cos(89.9)/cos(89.)=.1;     start at 86S, 86N

       do 260 ij=5,1,-1 
           if (J181(ij) .lt. 0.d0 .and. LAT181(ij) .lt. LAT(1)) then
               if (ij .eq. 1)  J181(1) = J181(2)*.1d0

               if (ij .gt. 1)  J181(ij) = J181(ij+1)*
     >                             DCOSD(lat181(ij))/DCOSD(lat181(ij+1))
           endif
                                                                      ! ensure NON-NEGATIVE J181
           if (J181(ij) .lt. 0.d0) J181(ij) = 0.d0
 260   CONTINUE

 
       do 270 ij=177,181
           if (J181(ij) .lt. 0.d0 .and. LAT181(ij) .gt. LAT(L$)) then
               if (ij .eq. 181)  J181(181) = J181(180)*.1d0

               if (ij .lt. 181)  J181(ij) = J181(ij-1)*
     >                             DCOSD(lat181(ij))/DCOSD(lat181(ij-1))
           endif
                                                                      ! ensure NON-NEGATIVE J181
           if (J181(ij) .lt. 0.d0) J181(ij) = 0.d0
 270   CONTINUE


       RETURN
       END
