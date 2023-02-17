	SUBROUTINE COLDEN

C  new version can integrate up to 10 species -- 
C      define species index in array INSP, # OF SPECS = NSPC
C
C   NOW, just compute DU in each grid box, using the species value at the grid point, 
C         and the DELTAZ(L$,Z$X) (in KM) value at the species grid points, which is the 
C         difference between the geopotential height of the top edge of the box minus 
C         the geop. hgt. of the bottom edge of the box - EF 2/21/02
C

        include "com2d.h"

	INTEGER NSPC

!        PARAMETER (NSPC=3)
        PARAMETER (NSPC=7) !see notes below

	DIMENSION NOBO2(L$),NOBSP(NSPC,L$), INSP(NSPC)

!BThomas Nov2017: additional variables for HNO3 rainout (see below)
	REAL RAINCOL
	INTEGER IK1
!
!        DATA INSP/3,4,5/
!BThomas Nov2017 - *** NOTE *** : One can add more species for column density calculations here,
!                                 BUT! be careful to leave 3,4,5 in the first 3 slots
!                                      in the INSP array, 
!                                 otherwise calculations elsewhere will FAIL
!                                 because other code ASSUMES that (e.g.) O3 is in slot 2 
!                                        ( that is, in ncolgd(2,ij,ik) )
!
	DATA INSP/3,4,5,6,10,31,59/
! 3,  4,  5,  6,   10,   31,  59
! O2, O3, NO, NO2, HNO3, NOy, H2O precip
!--

ccccccc        DATA INSP/3,39/
C                          if doing model tracer tests, integrate Ox (since ozone is not computed) 
        SAVE 


C  MAIN COMMON ASSEMBLY FOR SPEC TRANS 2-D MODEL --  INITIALIZE/UPDATE COLUMN DENSITIES 
C INEFFICIENT VERSION;RECALCULATES ENTIRE GRID - NOTE CRUDE TREATMENT OF OVERBURDEN OF O2 AND O3?
C   WORKING DOWN IN EACH LAT. BAND
C  first do O2 column,               PRESS(Z$X) in COMMON (pressures at the box centers)

      DO 100 IJ=1,L$
        NOBO2(IJ) = .21*(PRESS(Z$)*1000./980.)*(6.02E23/28.964)
        NCOLGD(1,IJ,Z$)=NOBO2(IJ)
          DO 200 IK1=2,Z$
          IK=Z$+1-IK1
 200	  NCOLGD(1,IJ,IK)=.21*(PRESS(IK)*1000./980.)*(6.02E23/28.964)
 100	  CONTINUE

C  Now loop through other requested species, FIRST DEFINE TOTAL COLUMN ABOVE TOP MODEL LEVEL
C    also load in DU at each level into CN(80,L$,Z$), now use DELTAZ(L$,Z$X) (in KM) 

        DO 1000 IJ=1,L$
          NOBSP(2,IJ) = C(INSP(2),IJ,Z$)*5.e5
         NOBSP(3,IJ) = C(INSP(3),IJ,Z$)*0.0e0

         DO 3000 IP=2,NSPC
            SPCOL = NOBSP(IP,IJ)
c                                  integrate top-down - SPCOL is the column at the CENTER of the box
C                              ie, DELTAZ/2., store 1/2 box in NCOLGD, then add in the other 1/2 to SPCOL
          DO 400 IK = Z$,1,-1
            colcen = C(INSP(IP),IJ,IK)*DELTAZ(IJ,IK)*1.D5/2.
            SPCOL = SPCOL + colcen

            NCOLGD(IP,IJ,IK) = SPCOL
 400        SPCOL = SPCOL + colcen

CCC            SPCOL = SPCOL + C(INSP(IP),IJ,IK)*DELTAZ(IJ,IK)*1.D5
ccc
cccc                                                      OLD METHOD - integrate by log-denisty 
cccc       DO 400 IK1 = 1,Z$-1
cccc          IK = Z$-IK1
cccc          SPCOL = SPCOL + 1.E5*(ZKM(IJ,IK+1)-ZKM(IJ,IK))*
cccc     c    DEXP((DLOG(C(INSP(IP),IJ,IK+1)) + DLOG(C(INSP(IP),IJ,IK)))*.5)
cccc 400	  NCOLGD(IP,IJ,IK) = SPCOL

 3000	CONTINUE
C                      now compute DU/km in each grid box, and add on TOTAL COLUMN ABOVE TOP MODEL LEVEL
       do 410 ik=1,Z$
 410	  cn(80,ij,ik) = C(4,IJ,IK)*DELTAZ(IJ,IK)*1.D5/DELTAZ(IJ,IK)

          cn(80,ij,Z$) = cn(80,ij,Z$) + C(4,IJ,Z$)*5.d5  

!BThomas Nov2017 --- need "rainout" of HNO3.  Taking this block from my previous
!                    version of the GSFC 2D model.
!                    Computes the vertically integrated "wet deposition" of HNO3
!                    using product of HNO3 concentration and rainout (constit. 59)

	RAINCOL = 0.0
	do ik=1,Z$
	   hno3rain(ij,ik)=0.0
	end do
	DO  IK1 = 1,Z$-1
	   IK = Z$-IK1 !computes going *down* from top to ground
	   HNO3RAIN(IJ,IK) = C(59,IJ,IK)*C(10,IJ,IK)
	   RAINCOL = RAINCOL + 1.D5*(ZKM(IJ,IK+1)-ZKM(IJ,IK))*
     >	    DEXP((DLOG(HNO3RAIN(IJ,IK+1))+DLOG(HNO3RAIN(IJ,IK)))*.5)
	   COLDENHNO3(IJ,IK) = RAINCOL
!debug:	   if(ij.eq.1) then
!	      print *, c(59,ij,ik),c(10,ij,ik),
!     >        hno3rain(ij,ik),hno3rain(ij,ik+1),
!     >        zkm(ij,ik),zkm(ij,ik+1),
!     > 	      raincol, coldenhno3(ij,ik)
!---	   end if
	end do !ik1 loop


 1000	CONTINUE !latitude loop


C  Now DO ozone determined from transported Ox - CN(79), FIRST DEFINE TOTAL COLUMN ABOVE TOP MODEL LEVEL
cc79
cc79        DO 1100 IJ=1,L$
cc79           SPCOL = C(4,IJ,Z$)*5.e5/C(39,ij,Z$)*C(79,ij,Z$) 
c  cc79                                                                    integrate top-down
cc79          DO 1400 IK = Z$,1,-1
cc79      SPCOL=SPCOL+ C(4,ij,ik)*DELTAZ(ij,ik)*1.D5/C(39,ij,ik)*C(79,ij,ik) 
cc79 1400      NCOLGD(3,ij,ik) = SPCOL
cc79
cc79 1100	CONTINUE
cc79


      RETURN
      END
