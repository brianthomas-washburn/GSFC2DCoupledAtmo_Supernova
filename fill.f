       SUBROUTINE FILL(IJC,IKC)

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

       RETURN
       END
