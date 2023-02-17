       SUBROUTINE FILL_COL(IJC,IKC,ICLOCK)

       include "com2d.h"

       SAVE 

C   THIS SUBROUTINE FILLS THE COMMON BLOCK CLOC WITH
C   THE NECESSARY VARIABLES FOR USING THE CHEMISTRY
C   SUBROUTINES.

       LATP=LAT(IJC)
       ZP=ZKM(IJC,IKC)
       TS=TEMP(IJC,IKC)
       PRP=PRESS(IKC)
c	type *,'latp,zp,ts,prp',latp,zp,ts,prp
	DO 2 JS=1,2
2       NCOL(JS)=NCOLGD(JS,IJC,IKC)

C  overwrite with diurnally varying Ozone column here - NCOLD(2,18+3,L$,Z$) in COMMON

        NCOL(2) = NCOLD(1,ICLOCK,IJC,IKC)

       RETURN
       END
