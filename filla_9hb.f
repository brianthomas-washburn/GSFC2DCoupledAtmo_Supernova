C
	SUBROUTINE FILLA

C
C  THIS ROUTINE LOADS OUR SPECIES ARRAY INTO THE AER SPECIES ARRAY
C

        include "com2d.h"

        include "com_aerd.h"
        include "com_aerg.h"
        include "com_aers.h"


        REAL*8 frat, BRXF, CLXF, NOXF, noyst, xspecb, xspecc, xspecn
        REAL*8 bryo, xspecbo, brxfo, clyo, xspecco, clxfo, noyo, xspecno
        REAL*8 noxfo, snoyo, noyto

        COMMON/CAER1/ISPMAP(NFSP)

        DOUBLE PRECISION ALLFSP, ALLSSP, ALLPLQ, ALLJR
        COMMON /SPEC/ALLFSP(NLT,NHT,NFSP), ALLSSP(NLT,NHT,NSSP),
     $     ALLPLQ(NLT,NHT,NPLQ), ALLJR(NLT,NHT,NJR)

        COMMON /ALTITUDES/ZSTAR(NHT),ZZALT(NLT,NHT),ALAT(NLT)
        COMMON /AIR/TMP(NLT,NHT),PR(NLT,NHT),AIRVD(NLT,NHT)


C   WCONV276 - extra tropospheric ozone production (transferred from SOLVER)

        COMMON/COZPR/ozpr(360,L$,10)


        SAVE 



        DO 1000 IK=1,Z$
        DO 1000 IJ=1,L$

C
C  FAST SPECIES: ALLFSP(L$,Z$,NFSP): use noontime ratios from the day before - NOONRAT(S$,L$,Z$) in COMMON
C
C  NEW METHOD WITH EVERYTHING BEING TRANSPORTED, just use NOONRAT and UPDATED 24-hr average radical values
C                                                                              which have been transported
C     for Ox, HOx, CHx - nothing more needs to be done 
C
C  Ox
       ALLFSP(ij,ik,1) = NOONRAT(ispmap(1),ij,ik)*cn(ispmap(1),ij,ik)
       ALLFSP(ij,ik,2) = NOONRAT(ispmap(2),ij,ik)*cn(ispmap(2),ij,ik)
       ALLFSP(ij,ik,4) = NOONRAT(ispmap(4),ij,ik)*cn(ispmap(4),ij,ik)

       ALLFSP(ij,ik,3) = NOONRAT(ispmap(3),ij,ik)*cn(ispmap(3),ij,ik)
C                                      !O2(1D) use ratio Noon/diur avg*diur avg

C  HOx
       ALLFSP(ij,ik,5) = NOONRAT(ispmap(5),ij,ik)*cn(ispmap(5),ij,ik)
       ALLFSP(ij,ik,6) = NOONRAT(ispmap(6),ij,ik)*cn(ispmap(6),ij,ik)
       ALLFSP(ij,ik,7) = NOONRAT(ispmap(7),ij,ik)*cn(ispmap(7),ij,ik)
       ALLFSP(ij,ik,8) = NOONRAT(ispmap(8),ij,ik)*cn(ispmap(8),ij,ik)



C  CHx
       ALLFSP(ij,ik,9)  = NOONRAT(ispmap(9),ij,ik)*cn(ispmap(9),ij,ik)
       ALLFSP(ij,ik,10) = NOONRAT(ispmap(10),ij,ik)*cn(ispmap(10),ij,ik)
       ALLFSP(ij,ik,11) = NOONRAT(ispmap(11),ij,ik)*cn(ispmap(11),ij,ik)



c  Ix - not used

        ALLFSP(IJ,IK,12) = 0.0E0
	ALLFSP(IJ,IK,13) = 0.0E0
	ALLFSP(IJ,IK,14) = 0.0E0
	ALLFSP(IJ,IK,15) = 0.0E0
	ALLFSP(IJ,IK,16) = 0.0E0
	ALLFSP(IJ,IK,17) = 0.0E0
	ALLFSP(IJ,IK,18) = 0.0E0
	ALLFSP(IJ,IK,19) = 0.0E0
	ALLFSP(IJ,IK,20) = 0.0E0


C
C  for Bry, Cly, NOy - radical diurnal averages have already bee updated to the current 
C        diurnal average family value (from  previous time step in SOLVER)
C        except x-species (ClONO2, ClNO2, BrONO2, BrNO2, BrCl) have NOT been updated by family updates
C
C        so also rescale NOONTIME values to the diurnal average family, w/o the x-species, as done in SOLVER
C        this also checks for problems w/ using ratios from previous time step, DON'T rescale if RAT LE 0
C
C  Bry -
C

       ALLFSP(ij,ik,41) = NOONRAT(ispmap(41),ij,ik)*cn(ispmap(41),ij,ik)
       ALLFSP(ij,ik,42) = NOONRAT(ispmap(42),ij,ik)*cn(ispmap(42),ij,ik)
       ALLFSP(ij,ik,43) = NOONRAT(ispmap(43),ij,ik)*cn(ispmap(43),ij,ik)
       ALLFSP(ij,ik,44) = NOONRAT(ispmap(44),ij,ik)*cn(ispmap(44),ij,ik)
       ALLFSP(ij,ik,45) = NOONRAT(ispmap(45),ij,ik)*cn(ispmap(45),ij,ik)
       ALLFSP(ij,ik,46) = NOONRAT(ispmap(46),ij,ik)*cn(ispmap(46),ij,ik)
       ALLFSP(ij,ik,47) = NOONRAT(ispmap(47),ij,ik)*cn(ispmap(47),ij,ik)
       ALLFSP(ij,ik,48) = NOONRAT(ispmap(48),ij,ik)*cn(ispmap(48),ij,ik)


       brxf = ALLFSP(ij,ik,41) + ALLFSP(ij,ik,42) + ALLFSP(ij,ik,43)
     >      + ALLFSP(ij,ik,44) + 2.*ALLFSP(ij,ik,47)

       xspecb = ALLFSP(ij,ik,45) + ALLFSP(ij,ik,46) + ALLFSP(ij,ik,48)

       frat = (CN(48,ij,ik) - xspecb)/brxf
       if (frat .le. 0.) then
          bryo    = CN(48,ij,ik)/M(IJ,IK)*1.d12
          xspecbo = xspecb/M(IJ,IK)*1.d12
          brxfo   = brxf/M(IJ,IK)*1.d12

          write(1575,225) iyr, iday360, ij, ik,
     >            frat, bryo, xspecbo, brxfo
 225      format('IN FILLA, Br-frat le 0 : ', 4I4, 2x, 1P4D13.4)

          frat = 1.d0
       endif

C         ! this ensures the total of the NOON Bry species = CN(48)

       ALLFSP(ij,ik,41) = ALLFSP(ij,ik,41)*frat
       ALLFSP(ij,ik,42) = ALLFSP(ij,ik,42)*frat
       ALLFSP(ij,ik,43) = ALLFSP(ij,ik,43)*frat
       ALLFSP(ij,ik,44) = ALLFSP(ij,ik,44)*frat
       ALLFSP(ij,ik,47) = ALLFSP(ij,ik,47)*frat



C  Cly - 


       ALLFSP(ij,ik,21) = NOONRAT(ispmap(21),ij,ik)*cn(ispmap(21),ij,ik)
       ALLFSP(ij,ik,22) = NOONRAT(ispmap(22),ij,ik)*cn(ispmap(22),ij,ik)
       ALLFSP(ij,ik,23) = NOONRAT(ispmap(23),ij,ik)*cn(ispmap(23),ij,ik)
       ALLFSP(ij,ik,24) = NOONRAT(ispmap(24),ij,ik)*cn(ispmap(24),ij,ik)
       ALLFSP(ij,ik,25) = NOONRAT(ispmap(25),ij,ik)*cn(ispmap(25),ij,ik)
       ALLFSP(ij,ik,26) = NOONRAT(ispmap(26),ij,ik)*cn(ispmap(26),ij,ik)
       ALLFSP(ij,ik,27) = NOONRAT(ispmap(27),ij,ik)*cn(ispmap(27),ij,ik)
       ALLFSP(ij,ik,28) = NOONRAT(ispmap(28),ij,ik)*cn(ispmap(28),ij,ik)
       ALLFSP(ij,ik,29) = NOONRAT(ispmap(29),ij,ik)*cn(ispmap(29),ij,ik)
       ALLFSP(ij,ik,30) = NOONRAT(ispmap(30),ij,ik)*cn(ispmap(30),ij,ik)
       ALLFSP(ij,ik,31) = NOONRAT(ispmap(31),ij,ik)*cn(ispmap(31),ij,ik)

       clxf = ALLFSP(ij,ik,21) + ALLFSP(ij,ik,22) + ALLFSP(ij,ik,23)
     >      + ALLFSP(ij,ik,24) + ALLFSP(ij,ik,26)
     >      + ALLFSP(ij,ik,27) + ALLFSP(ij,ik,28)
     >   + 2.*ALLFSP(ij,ik,29) + 2.*ALLFSP(ij,ik,30)

       xspecc = ALLFSP(ij,ik,25) + ALLFSP(ij,ik,31) + ALLFSP(ij,ik,48)

       frat = (CN(33,ij,ik) - xspecc)/clxf
       if (frat .le. 0.) then
          clyo    = CN(33,ij,ik)/M(IJ,IK)*1.d9
          xspecco = xspecc/M(IJ,IK)*1.d9
          clxfo   = clxf/M(IJ,IK)*1.d9

          write(1576,226) iyr, iday360, ij, ik,
     >            frat, clyo, xspecco, clxfo
 226      format('IN FILLA, Cl-frat le 0 : ', 4I4, 2x, 1P4D13.4)

          frat = 1.d0
       endif


C         ! this ensures the total of the NOON Cly species = CN(33)

       ALLFSP(ij,ik,21) = ALLFSP(ij,ik,21)*frat
       ALLFSP(ij,ik,22) = ALLFSP(ij,ik,22)*frat
       ALLFSP(ij,ik,23) = ALLFSP(ij,ik,23)*frat
       ALLFSP(ij,ik,24) = ALLFSP(ij,ik,24)*frat
       ALLFSP(ij,ik,26) = ALLFSP(ij,ik,26)*frat
       ALLFSP(ij,ik,27) = ALLFSP(ij,ik,27)*frat
       ALLFSP(ij,ik,28) = ALLFSP(ij,ik,28)*frat
       ALLFSP(ij,ik,29) = ALLFSP(ij,ik,29)*frat
       ALLFSP(ij,ik,30) = ALLFSP(ij,ik,30)*frat



C  NOy - DONT include SHNO3 -  ! load in NOy w/o SHNO3 here, and load in ALLSSP(29) = AER NOx below
C

       noyst = cn(31,ij,ik) - cn(66,ij,ik)
       if (noyst .le. 0.) then
          noyto = cn(31,ij,ik)/M(IJ,IK)*1.d9
          snoyo = cn(66,ij,ik)/M(IJ,IK)*1.d9

          write(1577,227) iyr, iday360, ij, ik,
     >            noyto, snoyo, noyst
 227      format('IN FILLA, NOy le 0 : ', 4I4, 2x, 1P3D13.4)

          noyst = 1.E-12
       endif


       ALLFSP(ij,ik,32) = NOONRAT(ispmap(32),ij,ik)*cn(ispmap(32),ij,ik)
       ALLFSP(ij,ik,33) = NOONRAT(ispmap(33),ij,ik)*cn(ispmap(33),ij,ik)
       ALLFSP(ij,ik,34) = NOONRAT(ispmap(34),ij,ik)*cn(ispmap(34),ij,ik)
       ALLFSP(ij,ik,35) = NOONRAT(ispmap(35),ij,ik)*cn(ispmap(35),ij,ik)
       ALLFSP(ij,ik,36) = NOONRAT(ispmap(36),ij,ik)*cn(ispmap(36),ij,ik)
       ALLFSP(ij,ik,37) = NOONRAT(ispmap(37),ij,ik)*cn(ispmap(37),ij,ik)
       ALLFSP(ij,ik,38) = NOONRAT(ispmap(38),ij,ik)*cn(ispmap(38),ij,ik)
       ALLFSP(ij,ik,39) = NOONRAT(ispmap(39),ij,ik)*cn(ispmap(39),ij,ik)
       ALLFSP(ij,ik,40) = NOONRAT(ispmap(40),ij,ik)*cn(ispmap(40),ij,ik)

       noxf = ALLFSP(ij,ik,32) + ALLFSP(ij,ik,33) + ALLFSP(ij,ik,34)
     >      + ALLFSP(ij,ik,35) + 2.*ALLFSP(ij,ik,36) + ALLFSP(ij,ik,37)
     >      + ALLFSP(ij,ik,38) + ALLFSP(ij,ik,39) + ALLFSP(ij,ik,40)

       xspecn = ALLFSP(ij,ik,25) + ALLFSP(ij,ik,31)
     >        + ALLFSP(ij,ik,45) + ALLFSP(ij,ik,46)

       frat = (noyst - xspecn)/noxf
       if (frat .le. 0.) then
          noyo    = noyst/M(IJ,IK)*1.d9
          xspecno = xspecn/M(IJ,IK)*1.d9
          noxfo   = noxf/M(IJ,IK)*1.d9

          write(1578,228) iyr, iday360, ij, ik,
     >            frat, noyo, xspecno, noxfo
 228      format('IN FILLA, NOx-frat le 0 : ', 4I4, 2x, 1P4D13.4)

          frat = 1.d0
       endif


C      this ensures the total of the NOON NOy species = CN(31)-CN(66)

       ALLFSP(ij,ik,32) = ALLFSP(ij,ik,32)*frat
       ALLFSP(ij,ik,33) = ALLFSP(ij,ik,33)*frat
       ALLFSP(ij,ik,34) = ALLFSP(ij,ik,34)*frat
       ALLFSP(ij,ik,35) = ALLFSP(ij,ik,35)*frat
       ALLFSP(ij,ik,36) = ALLFSP(ij,ik,36)*frat
       ALLFSP(ij,ik,37) = ALLFSP(ij,ik,37)*frat
       ALLFSP(ij,ik,38) = ALLFSP(ij,ik,38)*frat
       ALLFSP(ij,ik,39) = ALLFSP(ij,ik,39)*frat
       ALLFSP(ij,ik,40) = ALLFSP(ij,ik,40)*frat




c                                 hydrocarbons not used
	ALLFSP(IJ,IK,49) = 0.0E0
	ALLFSP(IJ,IK,50) = 0.0E0
	ALLFSP(IJ,IK,51) = 0.0E0
	ALLFSP(IJ,IK,52) = 0.0E0
	ALLFSP(IJ,IK,53) = 0.0E0
	ALLFSP(IJ,IK,54) = 0.0E0



C  SLOW SPECIES:

	ALLSSP(IJ,IK,1) = cn(11,IJ,IK)
	ALLSSP(IJ,IK,2) = cn(18,IJ,IK)
	ALLSSP(IJ,IK,3) = 0.0E0
	ALLSSP(IJ,IK,4) = cn(17,IJ,IK)
	ALLSSP(IJ,IK,5) = cn(19,IJ,IK)
	ALLSSP(IJ,IK,6) = cn(37,IJ,IK)
	ALLSSP(IJ,IK,7) = cn(36,IJ,IK)
	ALLSSP(IJ,IK,8) = cn(34,IJ,IK)
	ALLSSP(IJ,IK,9) = cn(35,IJ,IK)
	ALLSSP(IJ,IK,10) = cn(40,IJ,IK)
	ALLSSP(IJ,IK,11) = cn(57,IJ,IK)
	ALLSSP(IJ,IK,12) = cn(49,IJ,IK)
	ALLSSP(IJ,IK,13) = cn(76,IJ,IK)
	ALLSSP(IJ,IK,14) = cn(77,IJ,IK)
	ALLSSP(IJ,IK,15) = cn(52,IJ,IK)
	ALLSSP(IJ,IK,16) = cn(53,IJ,IK)
	ALLSSP(IJ,IK,17) = cn(54,IJ,IK)
	ALLSSP(IJ,IK,18) = cn(55,IJ,IK)
	ALLSSP(IJ,IK,19) = cn(69,IJ,IK)
	ALLSSP(IJ,IK,20) = cn(70,IJ,IK)
	ALLSSP(IJ,IK,21) = cn(71,IJ,IK)
	ALLSSP(IJ,IK,22) = cn(51,IJ,IK)
	ALLSSP(IJ,IK,23) = cn(50,IJ,IK)
	ALLSSP(IJ,IK,24) = cn(75,IJ,IK)
	ALLSSP(IJ,IK,25) = cn(72,IJ,IK)

	ALLSSP(IJ,IK,26) = 0.0E0          ! CH3I
	ALLSSP(IJ,IK,27) = 0.0E0          ! CF3I

                                          ! HFCs:
	ALLSSP(IJ,IK,28) = cn(81,IJ,IK)
	ALLSSP(IJ,IK,29) = cn(82,IJ,IK)
	ALLSSP(IJ,IK,30) = cn(83,IJ,IK)
	ALLSSP(IJ,IK,31) = cn(84,IJ,IK)
	ALLSSP(IJ,IK,32) = cn(85,IJ,IK)
	ALLSSP(IJ,IK,33) = cn(86,IJ,IK)
	ALLSSP(IJ,IK,34) = cn(87,IJ,IK)
	ALLSSP(IJ,IK,35) = cn(88,IJ,IK)

	ALLSSP(IJ,IK,36) = cn(39,IJ,IK)
	ALLSSP(IJ,IK,37) = noyst

	ALLSSP(IJ,IK,38) = cn(33,IJ,IK)
	ALLSSP(IJ,IK,39) = cn(48,IJ,IK)
	ALLSSP(IJ,IK,40) = cn(56,IJ,IK)

	ALLSSP(IJ,IK,41) = 0.0E0
	ALLSSP(IJ,IK,42) = 0.0E0

C  initialize NOy (N + NO + NO2 + NO3 + 2*N2O5 + HO2NO2 + HONO + ClONO2 + ClONO + BrONO2) - NOW ALREADY DONE
c
ccc	ALLSSP(IJ,IK,30) =  ALLFSP(IJ,IK,21) + ALLFSP(IJ,IK,22) +
ccc     >  ALLFSP(IJ,IK,23)+ ALLFSP(IJ,IK,24) + 2.*ALLFSP(IJ,IK,25) +
cc     >  ALLFSP(IJ,IK,26) + ALLFSP(IJ,IK,27) 
cc     >  + ALLFSP(IJ,IK,15) + ALLFSP(IJ,IK,20) + ALLFSP(IJ,IK,32)

        ALLSSP(IJ,IK,43) = cn(66,IJ,IK)
	ALLSSP(IJ,IK,44) = cn(15,IJ,IK)
	ALLSSP(IJ,IK,45) = cn(67,IJ,IK)

	ALLSSP(IJ,IK,46) = 0.0E0
	ALLSSP(IJ,IK,47) = 0.0E0
	ALLSSP(IJ,IK,48) = 0.0E0
	ALLSSP(IJ,IK,49) = 0.0E0
	ALLSSP(IJ,IK,50) = 0.0E0
	ALLSSP(IJ,IK,51) = 0.0E0
	ALLSSP(IJ,IK,52) = 0.0E0
	ALLSSP(IJ,IK,53) = 0.0E0
	ALLSSP(IJ,IK,54) = 0.0E0
	ALLSSP(IJ,IK,55) = 0.0E0
	ALLSSP(IJ,IK,56) = 0.0E0

	ALLSSP(IJ,IK,57) = cn(20,IJ,IK)
	ALLSSP(IJ,IK,58) = cn(31,IJ,IK)
	ALLSSP(IJ,IK,59) = cn(33,IJ,IK)
	ALLSSP(IJ,IK,60) = cn(3,IJ,IK)
	ALLSSP(IJ,IK,61) = N2(IJ,IK)
	ALLSSP(IJ,IK,62) = m(IJ,IK)


C   THIS SUBROUTINE FILLS THE COMMON BLOCKS WITH
C   THE NECESSARY VARIABLES FOR USING THE CHEMISTRY
C   SUBROUTINES.

c       LATP=LAT(IJ)
       ZZALT(IJ,IK)=ZKM(IJ,IK)
       TMP(IJ,IK)=TEMP(IJ,IK)
       PR(IJ,IK)=PRESS(IK)
       AIRVD(IJ,IK)=M(IJ,IK)

c	type *,'latp,zp,ts,prp',latp,zp,ts,prp
c	DO 2 JS=1,2
c2       NCOL(JS)=NCOLGD(JS,IJ,IK)



C
C  WCONV276 - add lower trop ozone production, use deficit from climatology w/ .5 day time scale
C   ie, use only if model ozone is too low (due to biomass burning, industrial pollution, etc)
C
C  IYRCT, OZBC(360,L$,10) in COMMON (in ppmv) interpolated to lowest 10 levels of current grid
C   load into ozpr(360,L$,10) (in COMMON above)
C
C   If TD run, do for first 3 years only (moved from SOLVER) - use diurnal average ozone
C                    if SState run, do for all years - LSTSTATE logical in COMMON

        IF (ik .le. 10) THEN
          if (.NOT. LSTSTATE  .and.  IYRCT .le. 3  .or.  LSTSTATE) then
             ozpr(iday360,ij,ik) = 
     >    (OZBC(iday360,ij,ik)*1.E-6*M(ij,ik) - cn(4,ij,ik))/(.5*86400.)

             if (ozpr(iday360,ij,ik) .le. 0.) ozpr(iday360,ij,ik) = 0.
          endif
        ENDIF


 1000  CONTINUE


       RETURN
       END
