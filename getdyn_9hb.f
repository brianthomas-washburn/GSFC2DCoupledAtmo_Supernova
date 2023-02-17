
	SUBROUTINE GETDYN


C     THIS routine reads in transport arrays for proper YEAR. 
C         This routine is called ON FIRST DAY OF THE YEAR ONLY,
C
C      INPUT arrays are on the 45x76 grid, need to INTERPOLATE to PROPER grid
C

       include "com2d.h"


       REAL w45(45,89,74), kzz45(45,89,74), kyy45(46,88,74)
       REAL temp45(45,88,74)

       REAL trint1(45,89), trint2(46,88), trint3(45,88)
       REAL lattr(45), lattre(46), zalttr(88), zalttre(89)

       REAL trout1(L$,Z$X1), trout2(L$+1,Z$X), trout3(L$,Z$X)


       CHARACTER*45 CCF
       CHARACTER*4 CCYR



       SAVE

C
C
C  9JG run:  use MERRA Interannual data for 1979-2010, NCEP for 1958-1978, Clim otherwise;
C     now use internal write/read to convert INTEGER year to CHAR (Jan 2012)

         ccyr = 'clim'

         iyx = IYR + 1935
         IF (IYR .ge. 23  .and.  IYR .le. 75) write(ccyr, '(I4.4)' ) iyx

                                          ! set here if climatology is set for ALL YEARS in MAIN
         if (iclim .eq. 1) ccyr = 'clim'



C   read in transport for proper year (on 45x76 grid):
C     REAL w45(45,89,74), kzz45(45,89,74), kyy45(46,88,74), temp45(45,88,74) defined above
C

        ccf = '/misc/kah05/fleming/9gtrans/dyn_9gtp_'//ccyr//'.xdr'

        OPEN (75, FILE = ccf, convert="big_endian",
     >    form = 'binary', access = 'sequential', status = 'old')
              READ (75) w45
              READ (75) kyy45
              READ (75) kzz45
              READ (75) temp45
        CLOSE (75)


C
C   define latitudes/altitudes of input arrays: lattr(45), lattre(46), zalttr(88), zalttre(89)
C   
       do 301 ij=1,45 
 301      lattr(ij) = (ij-1)*4. - 88.

       do 304 ij=1,46
 304	  lattre(ij) = (ij-1)*4. - 90.


          zalttre(1) = 0.0 

          do 310 ik=2,61
 310	     zalttre(ik) = zalttre(ik-1) + 1.

          do 311 ik=62,77
 311	     zalttre(ik) = zalttre(ik-1) + 2.

          do 312 ik=78,89
 312	     zalttre(ik) = zalttre(ik-1) + 2.


          do 315 ik=1,88
 315	     zalttr(ik) = (zalttre(ik) + zalttre(ik+1))/2.


C  now interpolate to current model grid:  
C   w45(45,89,74), kzz45(45,89,74), kyy45(46,88,74), temp45(45,88,74) defined above
C
C   trint1(45,89), trint2(46,88), trint3(45,88),   lattr(45), lattre(46), zalttr(88), zalttre(89)
C
C   WALL(L$,Z$X1,74), KZZTOT(L$,Z$X1,74), KYYALL(L$+1,Z$X,74), TEMPALL(L$,Z$X,74) all in COMMON
C                                                   LAT4(L$), LATEG4(L$+1),  ZALT(Z$X), ZALTE(Z$X1)
C   trout1(L$,Z$X1), trout2(L$+1,Z$X), trout3(L$,Z$X)
C
C
      DO 5000 iim=1,74

           do 551 ik=1,89
           do 551 ij=1,45
 551	      trint1(ij,ik) = w45(ij,ik,iim)

           CALL BINTERP(lattr, 45, zalttre, 89, trint1,
     >                  LAT4, L$, ZALTE, Z$X1, 0, 0, TROUT1)

           do 751 ik=1,Z$X1
           do 751 ij=1,L$
 751	      WALL(ij,ik,iim) = trout1(ij,ik)



           do 553 ik=1,89
           do 553 ij=1,45
 553	      trint1(ij,ik) = kzz45(ij,ik,iim)

           CALL BINTERP(lattr, 45, zalttre, 89, trint1,
     >                  LAT4, L$, ZALTE, Z$X1, 0, 0, TROUT1)

           do 753 ik=1,Z$X1
           do 753 ij=1,L$
 753	      KZZTOT(ij,ik,iim) = trout1(ij,ik)

 

           do 555 ik=1,88
           do 555 ij=1,46
 555	      trint2(ij,ik) = kyy45(ij,ik,iim)

           CALL BINTERP(lattre, 46, zalttr, 88, trint2,
     >                  LATEG4, L$+1, ZALT, Z$X, 0, 0, TROUT2)

           do 755 ik=1,Z$X
           do 755 ij=1,L$+1
 755	      KYYALL(ij,ik,iim) = trout2(ij,ik)



           do 557 ik=1,88
           do 557 ij=1,45
 557	      trint3(ij,ik) = temp45(ij,ik,iim)

           CALL BINTERP(lattr, 45, zalttr, 88, trint3,
     >                  LAT4, L$, ZALT, Z$X, 0, 0, TROUT3)

           do 757 ik=1,Z$X
           do 757 ij=1,L$
 757	      TEMPALL(ij,ik,iim) = trout3(ij,ik)

 5000   CONTINUE


	RETURN
	END
