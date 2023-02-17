C
       SUBROUTINE VWMASS
C
C
C  SUBROUTINE  VWMASS_9AA.F --  ROUTINE TO  ADJUST W VERTICAL VELOCITY FIELD FOR MASS BALANCE 
c     THEN COMPUTE V FIELD BY MASS CONTINUITY for PPM transport
C
C     now for non-hardwired latitude and irregular vertical grids 
C
C  INPUT (REAL*4 IN COMMON):
C
C        1) W1 (L$,Z$X1) = array of input vertical velocities, unadjusted, 
C                          defined at ZALTE (Z$X1), vector of altitudes of the box edges in KM
C
C 
C  OUTPUT (REAL*8 IN COMMON):
C
C        1) WBAR = (2,Z$X1) sum of the global cosine-mass weighted vert vel. at each pres lev, before/after
C        2) W = (L$,Z$X1) ARRAY OF adjusted vertical velocities
C        3) V = (L$+1,Z$X) array of meridional velocities, computed by standard mass continuity for a box
C
C    ALL VALUES HERE ARE IN CM !!!!
C

       include "com2d.h"


       SAVE

c
C  MEAN W AT EACH PRESSURE LEVEL MUST BE ZERO, so calculate the cosine weighted wbar
C     and ensure the global mean of w = 0. by the method of Shine (1989) which gives a better
C     correction than the original method.   NOTE THAT THE streamfunction formulation ensures
C     mass balance is close anyway. But it seems like the adjustment does make a slight difference
C     which can be seen when advecting a uniform mixing ratio. This is probably because the
C     streamfunction is computed for 37 latitudes, and then we reduce it to L$ latitudes, so
c     that the global integral may not still be zero. Accuracy is within the machine precision,
C     ~E-16 factor smaller than w* itself, for the global integral of w* along a pressure surface.
C
C                                    ! COSC(L$) in COMMON - cosine of box centers, TCOS = total (in COMMON)

         DO 200 ik=1,Z$X1 
c                                                         ! set limits on pre-adjusted w* at +- 5 cm/sec, 
c                                                         ! proper DT for Courant condition < 1 set in MAIN
              do 303 ij=1,L$
                 if (w1(ij,ik) .gt. 5.)  w1(ij,ik) = 5.d0
                 if (w1(ij,ik) .lt. -5.) w1(ij,ik) = -5.d0
cc                  w(ij,ik) = w1(ij,ik)
 303          continue


               WBAR(1,ik) = 0.0
            do 204 ij=1,L$
 204           WBAR(1,IK) = WBAR(1,ik) + w1(IJ,IK)*cosc(ij)

            do 205 ij=1,L$
 205           w(ij,ik) = w1(ij,ik) - WBAR(1,ik)/tcos

               WBAR(2,ik) = 0.0
            do 207 ij=1,L$
 207           WBAR(2,ik) = WBAR(2,ik) + w(IJ,ik)*cosc(IJ)       ! check adjusted w* for mass balance
 200     CONTINUE

C                           ensure that w* at bottom and top levels = 0.0 for unadjusted and adjusted arrays
         DO 404 ij=1,L$
            w1(ij,1) = 0.0
            w1(ij,Z$X1) = 0.0

            w(ij,1) = 0.0
            w(ij,Z$X1) = 0.0
 404     continue

c
c
c  Now do the simple integration (by mass continuity) to get v* for each box (C-grid), 
C    Note: You get the same answer if you start the integration from the NPole or SPole
C    THis is now generalized for an irregular vertical grid - EF, 8/01
C
C  V(L$+1,Z$X), w(L$,Z$X1) ;  COSC(L$), COSE(L$+1/L$S) are in COMMON,  
C      ZALTE8(Z$X1), ZALT8(Z$X) are in KM  ;  DELYC (meters) is for the constituent grid (in COMMON)
C
C                                     ! v* = 0 at poles by definition
            do 400 ik=1,Z$X1-1                                   
               v(1,ik) = 0.d0  
 400           v(L$+1,ik) = 0.d0

C  just use the exp(-z/H) factor for the density variation, since the rho(0) factors cancel 
C     - EVERYTHING IS CONVERTED TO  CM  HERE!!!!   

         do 450 ik=1,Z$X1-1 
         do 450 ij=2,L$
           v(ij,ik) = (v(ij-1,ik)*cose(ij-1) - 
     >                (w(ij-1,ik+1)*DEXP(-zalte8(ik+1)/7.d0)  
     >               - w(ij-1,ik)*DEXP(-zalte8(ik)/7.d0))
     > /(DEXP(-zalt8(ik)/7.d0)*(zalte8(ik+1)-zalte8(ik))*1000.d0*100.d0)
     >                *cosc(ij-1)*delyc*100.d0)/cose(ij)    
 450  CONTINUE

	RETURN
	END

