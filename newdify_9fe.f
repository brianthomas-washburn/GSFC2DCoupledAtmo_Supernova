c*************************************************************** 
c this subroutine will calculate the diffusion terms using a 
c mass conserving differencing scheme
c***************************************************************
c
c   ******  FOR  BASE9  (new grid resolution)  GENERALIZED  ********
C
C    everything in REAL*8 
C

       SUBROUTINE NEWDIFY


       include "com2d.h"


       REAL*8 mhalfy(L$+1,Z$X), kyyh(L$+1,Z$X)
       REAL*8 term2, term21, term22

       SAVE


C  radd, dphid, phid(L$), phihd(L$+1) are all in COMMON, REAL*8

                                           !xpi is in COMMON and defined in STREAMF
       radd = 6371000.d0*100.d0
       dphid = xpi/180.d0*dth
                                           ! double precision latitude angles in Radians
       DO 105 IJ=1,L$
 105       PHID(IJ) = LAT(IJ)*xpi/180.d0

       do 100 IJ=1,L$+1
 100       phihd(IJ) = (-90.d0 + dth*(IJ-1))*xpi/180.d0


c set up half step diffusion parameters, EKYY(L$+1,Z$X) is NOW defined at grid box sides
c this should happen whenever k or m changes 
C
C     INCLUDE DELTAZ (IN COMMON, and updated every day along with M) here
C     in the MHALFY calculation, then divide by DELTAZ in DIFNY below
C     this elimates long term drift in TOTAL MASS due to temp-dependent M changes
C     (since M~1/temp, and DELTAZ~temp), M(L$,Z$X), DELTAZ(L$,Z$X)
C     and you just DON'T get mass conservation very well without it (precision of ~10e-4 for REAL*8)

      do 330 IK=1,Z$X
      do 330 IJ=1,L$-1
            mhalfy(ij+1,ik) = (m(ij,ik)*deltaz(ij,ik) + 
     >                         m(ij+1,ik)*deltaz(ij+1,ik))/2.
            kyyh(ij+1,ik) = ekyy(ij+1,ik)
cccc     kyyh(ij+1,ik) = 3.d10                               !set to large Kyy everywhere just for testing
330     CONTINUE


!  initialize difny(T$,L$,Z$X) array for every time step

       do 441 IK=1,Z$X
       do 441 IJ=1,L$
       do 441 IT=1,T$
 441      difny(it,ij,ik) = 0.d0

C                                      maybe it's better just to do isentropic mixing say from 7-25 km here??
ccccc       i8 = INT(8./116.*Z$X)+1      
ccccc       i8 = INT(10./116.*Z$X)+1      
ccccc       i109 = INT(109./116.*Z$X)+1  
cccc       i109 = INT(25./116.*Z$X)+1  
                                     !if IKYYTH=0 as defined in CONTROL, then do Kyy on all pressure surfaces
cef       IF (IKYYTH .eq. 0) i8 = Z$X                                   ! (isentropic Kyy NOT done anywhere)


c  for new method w/ Kyy done on isentropes, only do this standard Kyy below 8 km and above 109 km
C    where isentropic mixing is not done
C    NOTE: need to TURN ON (UN-COMMENT OUT) the IF statement below (CEF) if using ISENTROPIC MIXING
C
                                                                ! difny(T$,L$,Z$X)
       DO 600 IT=1,ITRANS

       do 900 IK=1,Z$X
cef        if (ik .le. i8  .or.  ik .ge. i109) then
              term2=1.d0/(radd*DCOS(phid(1))*dphid)*
     c              mhalfy(2,IK)*kyyh(2,IK)/radd*DCOS(phihd(2))
     c              *(st(IT,2,IK)-st(IT,1,IK))/dphid

              difny(it,1,IK) = -1.d0*term2/deltaz(1,ik)

       do 1000 IJ=2,L$-1
              term21=1.d0/(radd*DCOS(phid(IJ))*dphid)*
     c              mhalfy(IJ+1,IK)*kyyh(IJ+1,IK)/radd*DCOS(phihd(IJ+1))
     c              *(st(IT,IJ+1,IK)-st(IT,IJ,IK))/dphid

              term22=-1.d0/(radd*DCOS(phid(IJ))*dphid)*
     c              mhalfy(IJ,IK)*kyyh(IJ,IK)/radd*DCOS(phihd(IJ))
     c              *(st(IT,IJ,IK)-st(IT,IJ-1,IK))/dphid

              term2=term21+term22
              difny(it,IJ,IK) = -1.d0*term2/deltaz(ij,ik)
1000    CONTINUE

              term2=-1.d0/(radd*DCOS(phid(L$))*dphid)*
     c              mhalfy(L$,IK)*kyyh(L$,IK)/radd*DCOS(phihd(L$))
     c              *(st(IT,L$,IK)-st(IT,L$-1,IK))/dphid

              difny(it,L$,IK) = -1.d0*term2/deltaz(L$,ik)
cef         endif
900     CONTINUE

600     CONTINUE

	RETURN
	END

