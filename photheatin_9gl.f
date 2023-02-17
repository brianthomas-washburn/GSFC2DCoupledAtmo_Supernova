C
        SUBROUTINE PHOTHEATIN

C
C  for coupled model, initialize the final INSOLATION CORRECTED photolysis array JZ
C     and the corrected XSRF array here (for heating rate calc consistent w/ photolysis)
C     this should be done BEFORE photolysis calculations so that values are 0 at night 
C     ihr=8,11 - (EF, Sept 2008)
C
C  common for O3, O2, NO2  cross sections*RFLUX for consistent photolysis/heating rates: 
C                  1=Ozone (UV);  2=O3-vis;   3=O2;  4=NO2;   XSRF is in photons/nm/second


        include "com2d.h"


        REAL*8 JZ, XSRF

        COMMON/JCORRECT/JZ(18,PH$,L$,Z$)
        COMMON/CXSRF/XSRF(18,4,L$,Z$)
 

        SAVE


        do 1000 ik=1,Z$
        do 1000 ij=1,L$
        do 1000 iph=1,PH$
        do 1000 itt2=1,18
 1000        JZ(itt2,iph,ij,ik) = 0.0


        do 2000 ik=1,Z$
        do 2000 ij=1,L$
        do 2000 iph=1,4
        do 2000 itt2=1,18
 2000        xsrf(itt2,iph,ij,ik) = 0.0

	RETURN
	END
