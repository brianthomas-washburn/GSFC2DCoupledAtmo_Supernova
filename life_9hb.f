      SUBROUTINE LIFETIME

      include "com2d.h"

      DIMENSION CLDEN(32,L$,Z$X),FLUXN2O(L$)
      REAL tlife(30,10,4), tlifez(2,30,10,Z$), COLTZ(Z$), aaburd(30)

      SAVE


C  9HB - THIS is for the SPARC LIFETIME ASSESSMENT


       DO 350 IJ=1,L$
          DO 350 IK=1,Z$
             CLDEN(1,IJ,IK)=C(34,IJ,IK)

             CLDEN(2,IJ,IK)=C(35,IJ,IK)

             CLDEN(3,IJ,IK)=C(36,IJ,IK)

             CLDEN(4,IJ,IK)=C(40,IJ,IK)

             CLDEN(5,IJ,IK)=C(52,IJ,IK)

             CLDEN(6,IJ,IK)=C(11,IJ,IK)

             CLDEN(7,IJ,IK)=C(18,IJ,IK)


             CLDEN(8,IJ,IK)=C(51,IJ,IK)

             CLDEN(9,IJ,IK)=C(50,IJ,IK)

             CLDEN(10,IJ,IK)=C(53,IJ,IK)

             CLDEN(11,IJ,IK)=C(55,IJ,IK)

             CLDEN(12,IJ,IK)=C(81,IJ,IK)

             CLDEN(13,IJ,IK)=C(82,IJ,IK)

             CLDEN(14,IJ,IK)=C(83,IJ,IK)


             CLDEN(15,IJ,IK)=C(54,IJ,IK)

             CLDEN(16,IJ,IK)=C(69,IJ,IK)
                  
             CLDEN(17,IJ,IK)=C(70,IJ,IK)

             CLDEN(18,IJ,IK)=C(37,IJ,IK)

             CLDEN(19,IJ,IK)=C(49,IJ,IK)

             CLDEN(20,IJ,IK)=C(75,IJ,IK)

             CLDEN(21,IJ,IK)=C(72,IJ,IK)

             CLDEN(22,IJ,IK)=C(84,IJ,IK)

             CLDEN(23,IJ,IK)=C(85,IJ,IK)

             CLDEN(24,IJ,IK)=C(86,IJ,IK)

             CLDEN(25,IJ,IK)=C(87,IJ,IK)

             CLDEN(26,IJ,IK)=C(88,IJ,IK)


             CLDEN(27,IJ,IK)=C(19,IJ,IK)

             CLDEN(28,IJ,IK)=C(17,IJ,IK)

             CLDEN(29,IJ,IK)=C(71,IJ,IK)

             CLDEN(30,IJ,IK)=C(76,IJ,IK)

             CLDEN(31,IJ,IK)=C(77,IJ,IK)

             CLDEN(32,IJ,IK)=C(20,IJ,IK)
 350   CONTINUE

C CO2
       DO 1350 IJ=1,L$
       DO 1350 IK=Z$+1,Z$X
 1350     CLDEN(32,IJ,IK)=C58(19,IJ,IK)


C
C  LOOP for the 26+4 SPARC compounds (INTEGRATED THROUGHOUT THE YEAR):
C
C
C  xlife(30,6,L$,Z$), xlifej(30,5,L$,Z$) loss rates (in COMMON)
C
C  2nd index (loss):   1=total loss;  2=J (tot);  3=O1D;   4=OH;   5=Cl
C
C   xlifej - J's for 5 wavelength bins: J39(PH$,5,L$,Z$)
C    1=Lyman alpha;   2=1695-1905A;   3=1905-2299A;  4=2299-2857A;   5=2857-7350A
C
C   CLOSS(40,10,L$), CLOSSST(40,10,L$), CLOSSTR(40,10,L$) are in COMMON:
C    also CLOSSM(40,10,L$) for mesosphere only is defined above
C
C   and store vertical profile in clossz(40,10,L$,Z$), columnz(40,L$,Z$) in COMMON
C
C   1st index is constituent
C   2nd index:  1=total loss;  2=J (tot);  3=O1D;   4=OH;   5=Cl
C    6=Lyman alpha;   7=1695-1905A;   8=1905-2299A;  9=2299-2857A;   10=2857-7350A
C
C    initialize arrays on 1st day of each year
C
      if (iday360 .eq. 1) then
        do 100 ij=1,L$
        do 100 ils=1,10
        do 100 iu=1,40
           column(iu,ij) = 0.0E0
           closs(iu,ils,ij) = 0.0E0
           closstr(iu,ils,ij) = 0.0E0
           clossst(iu,ils,ij) = 0.0E0
           clossm(iu,ils,ij) = 0.0E0

           do 105 ik=1,Z$
             columnz(iu,ij,ik) = 0.0E0
             clossz(iu,ils,ij,ik) = 0.0E0
 105       continue

 100    continue
      endif


      DO 145 IU=1,30

       DO 150 IJ=1,L$
       DO 110 IK=1,Z$

C  COLUMN is total column for entire atmosphere

         COLUMN(IU,IJ) = COLUMN(IU,IJ) + 
     >      (DELTAZ(IJ,IK)*CLDEN(IU,IJ,IK))*1.E5

         columnz(iu,ij,ik) = columnz(iu,ij,ik) + 
     >       DELTAZ(IJ,IK)*CLDEN(IU,IJ,IK)*1.E5

         do 170 ils=1,5
            CLOSS(IU,ils,IJ) = CLOSS(IU,ils,IJ) + 
     >        (DELTAZ(IJ,IK)*xlife(IU,ils,IJ,IK)*CLDEN(IU,IJ,IK))*1.E5

            CLOSSZ(iu,ils,ij,ik) = CLOSSZ(iu,ils,ij,ik) +
     >         xlife(iu,ils,ij,ik)*CLDEN(IU,IJ,IK)*DELTAZ(IJ,IK)*1.E5
 170     CONTINUE

         do 175 ils=6,10
           CLOSS(IU,ils,IJ) = CLOSS(IU,ils,IJ) + 
     >       (DELTAZ(IJ,IK)*xlifej(IU,ils-5,IJ,IK)*CLDEN(IU,IJ,IK))*1.E5

           CLOSSZ(iu,ils,ij,ik) = CLOSSZ(iu,ils,ij,ik) +
     >       DELTAZ(IJ,IK)*xlifej(IU,ils-5,IJ,IK)*CLDEN(IU,IJ,IK)*1.E5
 175     CONTINUE


C  also separate troposphere, stratosphere, and mesosphere
C    CLOSSTR(40,10,L$), CLOSSST(40,10,L$) are in COMMON
C    use ITROP360(L$,360) (in COMMON), ie, the tropopause FORTRAN indicies
C     for the current Z$ grid, defined in TROPKZZ
C   ie, the troposphere is at and below the ITROP360 level (defined in TROPKZZ)
C
C   CLOSSM(40,10,L$) is for the mesosphere (above 1 mbar, 49 km), ZALT(Z$X) in COMMON
C
        IF (ik .le. ITROP360(ij,iday360)) then
            do 171 ils=1,5
              CLOSSTR(IU,ils,IJ) = CLOSSTR(IU,ils,IJ) + 
     >          (DELTAZ(IJ,IK)*xlife(IU,ils,IJ,IK)*CLDEN(IU,IJ,IK))*1.E5
 171        CONTINUE

            do 176 ils=6,10
              CLOSSTR(IU,ils,IJ) = CLOSSTR(IU,ils,IJ) + 
     >       (DELTAZ(IJ,IK)*xlifej(IU,ils-5,IJ,IK)*CLDEN(IU,IJ,IK))*1.E5
 176        CONTINUE
        ENDIF


        IF (ik .gt. ITROP360(ij,iday360)  .and.  zalt(ik) .le. 49.) then
            do 172 ils=1,5
              CLOSSST(IU,ils,IJ) = CLOSSST(IU,ils,IJ) + 
     >          (DELTAZ(IJ,IK)*xlife(IU,ils,IJ,IK)*CLDEN(IU,IJ,IK))*1.E5
 172        CONTINUE

            do 177 ils=6,10
              CLOSSST(IU,ils,IJ) = CLOSSST(IU,ils,IJ) + 
     >       (DELTAZ(IJ,IK)*xlifej(IU,ils-5,IJ,IK)*CLDEN(IU,IJ,IK))*1.E5
 177        CONTINUE
        ENDIF


        IF (zalt(ik) .gt. 49.) then
            do 173 ils=1,5
              CLOSSM(IU,ils,IJ) = CLOSSM(IU,ils,IJ) + 
     >          (DELTAZ(IJ,IK)*xlife(IU,ils,IJ,IK)*CLDEN(IU,IJ,IK))*1.E5
 173        CONTINUE

            do 178 ils=6,10
              CLOSSM(IU,ils,IJ) = CLOSSM(IU,ils,IJ) + 
     >       (DELTAZ(IJ,IK)*xlifej(IU,ils-5,IJ,IK)*CLDEN(IU,IJ,IK))*1.E5
 178        CONTINUE
        ENDIF

110    CONTINUE
150    CONTINUE
145    CONTINUE


C  do lifetime calculation on last day of year, and write out
C    NOTE: AREA(L$) is the fraction of the globe at each latitude
C    no need to include 5.1e18 cm2 factor here since it cancels out
C    in the lifetime calc, ie, CLIF = COLT/CLST; COLT = global avg column
C
C    get annual avg burden (# molecules) from COLT - aaburd(30) (w/o 5.1e18 cm2 fac)
C
C    tlife(30,10,4):  1=total;  2=troposphere;  3=stratosphere;  4=mesosphere

      IF (IDAY360 .EQ. 360) THEN

       DO 5300 IU=1,30

          COLT=0.0
          DO 5177 IJ=1,L$
 5177        COLT = COLT + COLUMN(IU,IJ)*AREA(IJ)

          aaburd(iu) = COLT/360.


C  COLTZ(Z$), columnz(40,L$,Z$)

          DO 5400 ik=1,Z$
             COLTZ(ik) = 0.0

             DO 5405 IJ=1,L$
 5405          COLTZ(ik) = COLTZ(ik) + COLUMNZ(IU,IJ,IK)*AREA(IJ)
 5400     CONTINUE


C  loop through the different loss processes:

       do 5050 ils=1,10

           CLST = 0.0
           CLSTTR = 0.0
           CLSTST = 0.0
           CLSTM = 0.0

           DO 5100 IJ=1,L$
             CLST   = CLST   + CLOSS(IU,ils,IJ)*AREA(IJ)
             CLSTTR = CLSTTR + CLOSSTR(IU,ils,IJ)*AREA(IJ)
             CLSTST = CLSTST + CLOSSST(IU,ils,IJ)*AREA(IJ)
             CLSTM  = CLSTM  + CLOSSM(IU,ils,IJ)*AREA(IJ)
5100       CONTINUE

           tlife(iu,ils,1) = COLT/CLST/3.1536e7
           tlife(iu,ils,2) = COLT/CLSTTR/3.1536e7
           tlife(iu,ils,3) = COLT/CLSTST/3.1536e7
           tlife(iu,ils,4) = COLT/CLSTM/3.1536e7


C  do vertical profiles - clossz(40,10,L$,Z$), COLTZ(Z$)
C    tlifez(2,30,10,Z$) is for :
C        1 = total atmospheric burden
C        2 = burden at particular altitude
C
          DO 5200 IK=1,Z$
             CLSTZ = 0.0

             DO 5205 IJ=1,L$
 5205           CLSTZ = CLSTZ + CLOSSZ(IU,ils,ij,ik)*AREA(IJ)

             tlifez(1,iu,ils,ik) = COLT/CLSTZ/3.1536e7
             tlifez(2,iu,ils,ik) = COLTZ(ik)/CLSTZ/3.1536e7
5200     CONTINUE

5050   CONTINUE

5300   CONTINUE

C  output to fort.151 on last day of each year

        write (151) iyr, tlife, tlifez, aaburd

      ENDIF   ! end IDAY360 = 360 loop

      RETURN
      END
