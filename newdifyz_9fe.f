c***************************************************************************************
c this subroutine will calculate the Kyz diffusion terms using a  
c mass conserving differencing scheme  -  adapted from the AER 2D model  (EF, 8/13/03)
c***************************************************************************************
c
c   ******  FOR  BASE9  (new grid resolution)  GENERALIZED  ********
C
C    everything in REAL*8 
C

       SUBROUTINE NEWDIFYZ


       include "com2d.h"



       REAL*8 HFYZ(T$,L$+1,Z$X), HHEY(L$+1,Z$X), MEZ(L$,Z$X1)
       REAL*8 VFYZ(T$,L$,Z$X1), EDDYRT(T$,L$,Z$X), dzz1, dzz2, dzz3


       SAVE



C   ST(T$,L$,Z$X) is MIXING RATIO in COMMON

C  define scale height on box sides, just take average temperature, don't need values on on 1st and 
C    last sides (the poles), since the flux there is zero, but define w/ constant value so its non-zero
C      TEMP(L$,Z$X), HHEY(L$+1,Z$X) in CM
C
         do 301 ik=1,Z$X
         do 302 ij=2,L$
 302        hhey(ij,ik) = 29.3*1.d2*(TEMP(ij-1,ik) + TEMP(ij,ik))/2.d0

            hhey(1,ik) = hhey(2,ik) 
            hhey(L$+1,ik) = hhey(L$,ik) 
 301     CONTINUE


c  compute MEZ(L$,Z$X+1) at box bottoms and tops, again, don't need values at levels 1 and Z$X+1, 
C     but define it as something to avoid zeros, PRESSE8(Z$X1)
C
          do 315 ij=1,L$
          do 317 ik=2,Z$X
             mez(ij,ik) = PRESSE8(ik)*1.d3/
     >              (WTA*RD* (TEMP(ij,ik-1)+TEMP(ij,ik))*.5d0 *1.66d-24)
 317      CONTINUE

           mez(ij,1) = PRESSE8(1)*1.d3/(WTA*RD*TEMP(ij,1)*1.66d-24)
        mez(ij,Z$X1) = PRESSE8(Z$X1)*1.d3/(WTA*RD*TEMP(ij,Z$X)*1.66d-24)
 315   CONTINUE  



c   KYZ:      initialize  DIFNYZ(T$,L$,Z$X)


        do 570 ik=1,Z$X
        do 570 ij=1,L$
        do 570 it=1,T$
 570       DIFNYZ(it,ij,ik) = 0.d0


                                                               
       DO 600 IT=1,ITRANS

c
c  HORIZONTAL FLUX ....HFYZ(T$,L$+1,Z$X), ST(T$,L$,Z$X) is mixing ratio in COMMON, EKYZ(L$+1,Z$X1)
C      HHEY(L$+1,Z$X), DP(Z$X) is negative, set Left and right edges = 0.0
C
C      also account for grid box size weighting via ZALTE8(Z$X1), ZALT8(Z$X) - REAL*8
C            but don't do temperature dependent dz
C
C  NOTE: accounting for the variable vertical grid with this method changes the contour plots 
C    of H2O and CH4 ONLY a NEGLIGIBLE amount across levels w/ changing DZ (you barely notice the difference)
C    (ie, 2 km -> 0.5 km, and vice versa), versus assuming that the grid is always equidistant
C    and it's hard to say which method (w/ or w/o) is better. 
C    But we'll leave it in, since it seems the proper way to do it (MASS IS CONSERVED EITHER WAY)
C                                                  and everything cancels if the grid is equidistant.
C  
        do 150 ik=1,Z$X
           HFYZ(it,1,ik) = 0.d0
 150       HFYZ(it,L$+1,ik) = 0.d0


      do 200 ij=2,L$
        HFYZ(it,IJ,Z$X) = -EKYZ(IJ,Z$X)*(ST(IT,IJ-1,Z$X) + ST(IT,IJ,Z$X)
     >                           - ST(IT,IJ-1,Z$X-1) - ST(IT,IJ,Z$X-1))/
     >                             (zalt8(Z$X) - zalt8(Z$X-1))/
     >        (-DP(Z$X)*4.d0*HHEY(ij,Z$X))*(zalte8(Z$X1) - zalte8(Z$X))

        HFYZ(it,IJ, 1) = -EKYZ(IJ, 2)*(ST(IT,IJ-1, 2) + ST(IT,IJ, 2) 
     >                               - ST(IT,IJ-1, 1) - ST(IT,IJ, 1))/
     >                                     (zalt8(2) - zalt8(1))/
     >                  (-DP(1)*4.d0*HHEY(ij,1))*(zalte8(2) - zalte8(1))

      do 205 IK=2,Z$X-1
        dzz1 = zalt8(ik) - zalt8(ik-1)
        dzz2 = zalt8(ik+1) - zalt8(ik)

        dzz3 = zalte8(ik+1) - zalte8(ik)

        HFYZ(it,IJ,IK) = - (EKYZ(IJ,IK)*(ST(IT,IJ-1,IK) + ST(IT,IJ,IK) 
     >                         - ST(IT,IJ-1,IK-1) - ST(IT,IJ,IK-1))/dzz1 
     >  + EKYZ(IJ,IK+1)*(ST(IT,IJ-1,IK+1) + ST(IT,IJ,IK+1) 
     >                 - ST(IT,IJ-1,IK)   - ST(IT,IJ,IK))/dzz2)/
     >              (-DP(ik)*4.d0*HHEY(ij,ik))*dzz3
 205  CONTINUE
 200  CONTINUE

C
C
C   VERTICAL FLUX ..........VFYZ(T$,L$,Z$X1), ST(T$,L$,Z$X) is mixing ratio in COMMON
C       set bottom and top edges = 0.0,  EKYZ(L$+1,Z$X1), DELYC is for the constituent grid (in Meters)
C
C  
        do 170 ij=1,L$
           VFYZ(it,ij,1) = 0.d0
 170       VFYZ(it,ij,Z$X1) = 0.d0


C  don't need to scale by grid box size here for VFYZ, it's accounted for below in EDDYRT via the deltaz term

        do 410 ik=2,Z$X

           VFYZ(it,L$,IK) = -EKYZ(L$,IK)*(ST(IT,L$,IK) + ST(IT,L$,IK-1) 
     >        - ST(IT,L$-1,IK) - ST(IT,L$-1,IK-1))/(4.d0*delyc*100.d0)

           VFYZ(it, 1,IK) = -EKYZ( 2,IK)*(ST(IT,2,IK) + ST(IT,2,IK-1)
     >        - ST(IT,1,IK) - ST(IT,1,IK-1))/(4.d0*delyc*100.d0)

         do 415 ij=2,L$-1

           VFYZ(it,IJ,IK) = -(EKYZ(IJ,IK)*(ST(IT,IJ,IK) + ST(IT,IJ,IK-1) 
     >                    - ST(IT,IJ-1,IK) - ST(IT,IJ-1,IK-1))
     >       + EKYZ(IJ+1,IK)*(ST(IT,IJ+1,IK) + ST(IT,IJ+1,IK-1) 
     >            - ST(IT,IJ,IK) - ST(IT,IJ,IK-1)))/(4.d0*delyc*100.d0)

 415  CONTINUE
 410  CONTINUE


C
C   CALCULATION OF EDDYRT .............  EDDYRT(T$,L$,Z$X), M(L$,Z$X), DELTAZ(L$,Z$X)=-H*DP(Z$X) in KM
C        HFYZ(T$,L$+1,Z$X), VFYZ(T$,L$,Z$X1), MEZ(L$,Z$X1), 
C
C   LATEG(L$+1), COSEG(L$+1) are the latitudes and cosines of chemistry box edges (sides), REAL*8 in COMMON
C      COSC(L$) - cosines of box centers (chemistry grid points)
C

      do 100 ik=1,Z$X
      do 100 ij=1,L$

        EDDYRT(it,ij,ik) = (HFYZ(it,ij,ik)*coseg(ij) 
     >         - HFYZ(it,ij+1,ik)*coseg(ij+1))/(delyc*100.d0*cosc(ij))
c
     >     + (VFYZ(it,ij,ik)*mez(ij,ik) - VFYZ(it,ij,ik+1)*mez(ij,ik+1))
     >                      /(DELTAZ(ij,ik)*1.d5*M(ij,ik))                       ! convert DELTAZ to CM
                                                           
C
C  define DIFNYZ(T$,L$,Z$X) array here as -EDDYRT*M (convert to number density)
C               and add to the C array in XTRANS as usual
C
C  note: Getting the wrong sign for the Kyz contribution here (ie, DIFNYZ)
C     results in contours oppositely oriented to the isentropes, as opposed to parallel to the isentropes
C

        DIFNYZ(it,ij,ik) = -EDDYRT(it,ij,ik)*M(ij,ik)

cef
cef     DIFNYZ(it,ij,ik) = 0.0d0

 100   CONTINUE               


600     CONTINUE


cef
cef      if (iday360.ge.155) write(27) ekyz, st, HFYZ, VFYZ, EDDYRT, DIFNYZ


       RETURN
       END

