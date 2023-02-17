C
        SUBROUTINE SOLVER58(IJ,IK)

C
C   ROUTINE TO SOLVE only for transported species for levels Z$-Z$X, 
C      DON'T use loss rates from level Z$, just set all production and loss to zero for 
C      conservation purposes
c
C   EXCEPT for cycling between CO and CO2, and H2 and H2O, which are self-consistently calculated.
C      The calculations use J-coefficients and reaction rates from level Z$ - these are only
C      functions of temp, pressure, O2 and O3, so there is no dependence with CO, CO2, H2, H2O.
C      The calculations also need M which is based on temp, and an [O] and [OH] density which we
C      specify. This is OK as long as we use the same values for all calculations. For [O], we
C      linearly interpolate the mixing ratio between that computed at level Z$, and a number 
C      density of 2.e11 (cm-3) at level Z$X (115 km) taken from p. 267 of the Whitten and
C      Poppoff Aeronomy book. For [OH], the computed mixing ratio at 80-90 km indicate a drop 
C      off of ~.75 per level, so we specify [OH] for levels Z$+1-Z$X by using this value (constant) 
C      and the value at level Z$.
C

       include "com2d.h"


       REAL*8 aodens, ohdens


       SAVE

C
C  For upper levels w/ AER FAST chemistry  -  ONLY done for a few species with CHEMISTRY
C
C        the C58 array HAS BEEN TRANSPORTED
C
C         NO DIURNAL VARIATIONS HERE !!!!!!
C

c first specify [O] and [OH] densities
c
ccccccc        aodens = ((2.d11/m(ij,Z$X) - c(1,ij,Z$)/m(ij,Z$))/12.*(ik-Z$) +
ccccccc    c             c(1,ij,Z$)/m(ij,Z$))*M(IJ,IK)
ccccccc
ccccccc        ohdens = c(13,ij,Z$)/m(ij,Z$)*(.75**(ik-Z$))*M(IJ,IK)
ccccccc
C CO 
ccccccc        COPR = (J(41,IJ,Z$) + J(12,IJ,Z$))*C58(12,IJ,IK)
ccccccc       COLS = K(111,IJ,Z$)*aodens + K(36,IJ,Z$)*ohdens
ccccccc
ccccccc       CN58(11,IJ,IK) = (C58(11,IJ,IK) + COPR*DT)/(1. + COLS*DT)
ccccccc
ccccccc
C  CO2 
ccccccc	CO2P = K(36,IJ,Z$)*ohdens*C58(11,IJ,IK) +
ccccccc    >        K(111,IJ,Z$)*aodens*C58(11,IJ,IK)
ccccccc	CO2LOSS = j(12,ij,Z$) + j(41,ij,Z$)
ccccccc
ccccccc   	CN58(12,IJ,IK) = (C58(12,IJ,IK) + CO2P*DT)/(1. + co2loss*DT)
ccccccc
ccccccc
C H2
ccccccc        h2prod = j(25,IJ,Z$)*c58(25,IJ,IK)
ccccccc        h2loss = K(30,IJ,Z$)*ohdens
ccccccc
ccccccc       	cn58(10,IJ,IK) = (c58(10,IJ,IK) + h2prod*dt)/(1. + h2loss*dt)
ccccccc
ccccccc
C H2O
ccccccc        h2op = K(30,ij,Z$)*ohdens*c58(10,ij,ik)
ccccccc        h2ol = j(25,ij,Z$)
ccccccc
ccccccc       cn58(25,ij,ik) = (c58(25,ij,ik) + h2op*dt)/(1. + h2ol*dt)
ccccccc
C
C  Since the C58 array has been transported, need to update the CN58 array for the 
C       species that have NOT been updated with chemistry above
C       the IK loop here is for ik=Z$+1,Z$X
C  

	 do 970 is=1,ITRANS
ccccccc          if (is .ne. 10  .and.  is .ne. 11  .and.  is .ne. 12  .and. 
ccccccc    >        is .ne. 25) cn58(is,ij,ik) = c58(is,ij,ik) 

          cn58(is,ij,ik) = c58(is,ij,ik) 
 970     CONTINUE



C H2OA - HALOE interannual tropical trop. BC,  OLD cn(75)
C
CCCC       cn58(38,ij,ik) = (cn58(38,ij,ik) + h2op*dt)/(1. + h2ol*dt)

C H2OB - HALOE climatological tropical trop. BC, OLD cn(76)
C
CCCC       cn58(39,ij,ik) = (cn58(39,ij,ik) + h2op*dt)/(1. + h2ol*dt)


c
c  Set limit to 1.E-12 number density
C    (m at 115 km ~1.e12, and ~3.e19 at the ground, so the minimum mixing ratio will be ~3.E-31)

       DO 7382 III=1,ITRANS
       IF (CN58(III,IJ,IK) .LT. 1.D-12) CN58(III,IJ,IK) = 1.D-12
7382       CONTINUE


      RETURN
      END
