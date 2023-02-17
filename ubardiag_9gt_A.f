
        SUBROUTINE UBARDIAG 


C  read in yearly files for UBAR, for diagnostic output


        include "com2d.h"


        REAL*4 xjunk(L$S,Z$S,72), zz91(L$S,Z$S), zzzx(L$,Z$X)
        REAL*4 ubarin1(L$S,Z$S,72), epin1(L$S,Z$S,72), qyin1(L$S,Z$S,72)


c  common for UBAR climatology, used in both fixed and COUPLED models, read in TEMPIN

        COMMON/CUBAR/ UBARCL(91,117,74)


C  common for fixed model dynamical fields from MERRA data, clim avg for 1979-2010

        COMMON/CLMERRA/ xkyyclm(91,117,74), epclm(91,117,74), 
     >           tempclm(91,117,74), ubarclm(91,117,74),
     >           ehfclm(91,117,74), ehfzclm(91,117,74), qyclm(91,117,74)



        CHARACTER*47 CCF
        CHARACTER*4 CCYR(32)
 
       DATA CCYR/'1979', '1980', '1981', '1982', '1983', '1984', '1985',
     >           '1986', '1987', '1988', '1989', '1990', '1991', '1992',
     >           '1993', '1994', '1995', '1996', '1997', '1998', '1999',
     >           '2000', '2001', '2002', '2003', '2004', '2005', '2006',
     >           '2007', '2008', '2009', '2010'/


 	SAVE


C  load in UBAR, EP, QBARY MERRA Climatology (read in TEMPIN) -
C     UBARCLM (m/sec), EPCLM (m/sec/day), QYCLM(91,117,74) (1.e11 1/m/sec)
C  then if doing an INterannual run, read in appropriate yearly arrays from MERRA/MLS for 1979-2010
C           generated in IDL Linux, so use "big_endian"

        do 5000 iim=1,72
        do 5000 ik=1,Z$S
        do 5000 ij=1,L$S
           ubarin1(ij,ik,iim) = ubarclm(ij,ik,iim+1)
           epin1(ij,ik,iim) = epclm(ij,ik,iim+1)
           qyin1(ij,ik,iim) = qyclm(ij,ik,iim+1)
 5000   CONTINUE


      IF (iclim .ne. 1) then 
         IF (IYR .ge. 44  .and.  IYR .le. 75) then

           itp = iyr - 43

      ccf = 'kyy_merra_'//ccyr(itp)//'.xdr'

           OPEN (75, FILE = ccf, convert="big_endian",
     >         form = 'unformatted', access = 'sequential', 
     >         status = 'old')
             READ (75) xjunk
             READ (75) epin1
             READ (75) xjunk
             READ (75) ubarin1
             READ (75) xjunk
             READ (75) xjunk
             READ (75) qyin1
           CLOSE (75)

         ENDIF
      ENDIF

c
c  Interpolate to chemistry grid, expand to 74 time periods for SETDAILY
c    LATST4(L$S), ZSTR4(Z$S), LAT4(L$),ZALT(Z$X) all REAL*4 in COMMON,  zz91(L$S,Z$S), zzzx((L$,Z$X)
C
c    ubarin1(L$S,Z$S,72) -> UBARCG(L$,Z$X,74) all in COMMON:
c    epin1(L$S,Z$S,72)   -> EPCG(L$,Z$X,74)
c    qyin1(L$S,Z$S,72)   -> QYCG(L$,Z$X,74)
C
C
       DO 7000 iim=1,72

            do 517 ik=1,Z$S 
            do 517 ij=1,L$S
 517	       zz91(ij,ik) = ubarin1(ij,ik,iim) 

            CALL BINTERP(LATST4, L$S, ZSTR4, Z$S, zz91, 
     >                   LAT4, L$, ZALT, Z$X, 0, 0, ZZZX)

            DO 675 ik=1,Z$X
            DO 675 ij=1,L$
 675	       UBARCG(ij,ik,iim+1) = zzzx(ij,ik)


            do 617 ik=1,Z$S 
            do 617 ij=1,L$S
 617	       zz91(ij,ik) = epin1(ij,ik,iim) 

            CALL BINTERP(LATST4, L$S, ZSTR4, Z$S, zz91, 
     >                   LAT4, L$, ZALT, Z$X, 0, 0, ZZZX)

            DO 775 ik=1,Z$X
            DO 775 ij=1,L$
 775	       EPCG(ij,ik,iim+1) = zzzx(ij,ik)


            do 717 ik=1,Z$S 
            do 717 ij=1,L$S
 717	       zz91(ij,ik) = qyin1(ij,ik,iim) 

            CALL BINTERP(LATST4, L$S, ZSTR4, Z$S, zz91, 
     >                   LAT4, L$, ZALT, Z$X, 0, 0, ZZZX)

            DO 776 ik=1,Z$X
            DO 776 ij=1,L$
 776	       QYCG(ij,ik,iim+1) = zzzx(ij,ik)

 7000    CONTINUE


         DO 677 ik=1,Z$X
         DO 677 ij=1,L$
            UBARCG(ij,ik,1) = UBARCG(ij,ik,73)
            UBARCG(ij,ik,74) = UBARCG(ij,ik,2)

            EPCG(ij,ik,1) = EPCG(ij,ik,73)
            EPCG(ij,ik,74) = EPCG(ij,ik,2)

            QYCG(ij,ik,1) = QYCG(ij,ik,73)
            QYCG(ij,ik,74) = QYCG(ij,ik,2)
 677     CONTINUE


       RETURN
       END
