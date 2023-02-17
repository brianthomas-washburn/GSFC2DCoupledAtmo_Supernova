	SUBROUTINE AEROSOL_DBC

cdbc 8/8/96: This version of aerosols has been written for the chlorine Monte
c            Carlo study. It defines aerosol_all as a 264 x 18 x 30 array
c            and reads file climsad_79-85_v0.dat (or v1, etc). 
c            If the year is less than 1979 (iyr=44) then background conditions 
c            should be used.  For years past 1995, use background.
c       From /home/dbc/newmod/monte/aerosols_clmonte.f
c       Read '/misc/dbc03/dbc/climsad_79-95_v4.dat'
C
C  updated Oct. 2001, using David's new SAD file for 1979-2000 (EF)
C
C
        include "com2d.h"


        COMMON/CSAD/sadbk7(L$,Z$,14), sad7(L$,Z$,492), times7(492)

	common/aer/aer_dbc(264,L$,Z$)


	DIMENSION aer_dbcin(264,18,30), aeroin(18,46), aerout(L$,Z$)

        REAL sadbkr(36,80,14), sadr(36,80,492), lats7(36), zzs7(80)
        REAL sadin(36,80)


        OPEN (unit=80, 
     >          NAME='/misc/kah03/chj26/base9/modsad_79to00_new.dat',
     >		type='OLD', form='unFORMATTED', convert="big_endian")
           READ (80) aer_dbcin
        CLOSE(unit=80)


C 
C   now interpolate AER_DBCIN(264,18,30) array to new model grid, LATIN(18), ZZ46(46) are REAL*4 in COMMON 
C                                        LAT4(L$), ZALT90(Z$) are REAL*4 in COMMON

	do 500 ii=1,264

           do 601 ij=1,18
             do 602 ik=1,30
 602		aeroin(ij,ik) = aer_dbcin(ii,ij,ik)
             do 603 ik=31,46
 603		aeroin(ij,ik) = 0.0
 601       CONTINUE


            CALL BINTERP(LATIN, 18, zz46, 46, aeroin, 
     >                   LAT4, L$, ZALT90, Z$, 0, 0, AEROUT)


	do 200 ik=1,Z$
	do 200 ij=1,L$
    	   aer_dbc(ii,ij,ik) = aerout(ij,ik)
           if (aer_dbc(ii,ij,ik) .lt. 0.0) aer_dbc(ii,ij,ik) = 0.0         ! ensure non-negative aerosol
 200	CONTINUE

 500	CONTINUE

c  reset to 0.0 above 58 km (original level 30), aer_dbc(264,L$,Z$), ZALT90(Z$)

      DO 110 IK=1,Z$
         IF (zalt90(ik) .GT. 58.) THEN
		DO 111 IJ=1,L$
		DO 111 II=1,264
 111		   aer_dbc(ii,ij,ik) = 0.0
         ENDIF
 110		CONTINUE


              do 10 ik=1,Z$
              do 10 ij=1,L$
	      do 10 iq=1,264
                 aer_dbc(iq,ij,ik) = aer_dbc(iq,ij,ik)*1.e-8
 10	      continue

C
C  9JA - read in new SAD data set from CCMVal website, for 1963-2003, plus background clim
C    for Chapter 5 TD simulations  -   written out in IDL on Linux, so use "little_endian"
C
C    sadbkr(36,80,14), sadr(36,80,492), times7(492), lats7(36), zzs7(80)

         OPEN (177, FILE = '/misc/kah01/splife/chap5/sad_refc1.xdr',
     >       convert="little_endian", form = 'binary', 
     >       ACCESS = 'SEQUENTIAL', STATUS = 'OLD')
           READ (177) sadbkr
           READ (177) sadr
           READ (177) times7
           READ (177) lats7
           READ (177) zzs7
         CLOSE (177)


C  interpolate to current model grid -> SADIN(36,80), AEROUT(L$,Z$) => sadbk7(L$,Z$,14)

        do 700 ii=1,14

           do 701 ik=1,80
           do 701 ij=1,36
 701          sadin(ij,ik) = sadbkr(ij,ik,ii)


            CALL BINTERP(lats7, 36, zzs7, 80, SADIN, 
     >                   LAT4, L$, ZALT90, Z$, 0, 0, AEROUT)

C                         ! ensure non-negative aerosol, set=0 above 39.5 km
        do 705 ik=1,Z$
        do 705 ij=1,L$
    	   sadbk7(ij,ik,ii) = aerout(ij,ik)*1.e-8

           if (sadbk7(ij,ik,ii) .lt. 0.) sadbk7(ij,ik,ii) = 0.
           if (zalt90(ik) .GT. 39.5) sadbk7(ij,ik,ii) = 0.
 705	CONTINUE

 700	CONTINUE


C                           sadr(36,80,492) => sad7(L$,Z$,492)
        do 710 ii=1,492

           do 711 ik=1,80
           do 711 ij=1,36
 711          sadin(ij,ik) = sadr(ij,ik,ii)


            CALL BINTERP(lats7, 36, zzs7, 80, SADIN, 
     >                   LAT4, L$, ZALT90, Z$, 0, 0, AEROUT)

C                         ! ensure non-negative aerosol, set=0 above 39.5 km
        do 715 ik=1,Z$
        do 715 ij=1,L$
           sad7(ij,ik,ii) = aerout(ij,ik)*1.e-8

           if (sad7(ij,ik,ii) .lt. 0.) sad7(ij,ik,ii) = 0.
           if (zalt90(ik) .GT. 39.5) sad7(ij,ik,ii) = 0.
 715	CONTINUE

 710	CONTINUE


	RETURN
	END
