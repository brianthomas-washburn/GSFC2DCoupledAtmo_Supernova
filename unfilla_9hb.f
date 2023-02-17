C
C  Routine to load updated diurnally avged fast species from the CNDCA array into C and CN arrays
C    also load in HOx and CHx families for transport (NO additional chemistry) in C array ONLY
C
C  alos load noontime ratios into array NOONRAT(S$,L$,Z$) for FAST CHEM initialization on next day
C

	SUBROUTINE UNFILLA(CNDCA)


        include "com2d.h"

        include "com_aerd.h"
        include "com_aerg.h"
        include "com_aers.h"


        COMMON/CAER1/ISPMAP(NFSP)


C   COMMON for SPARC OH field in (#/cm3)
C
        COMMON/CCOH/ohm(14,L$,Z$), ohday(L$,Z$)


C  CNDCA = CN array for 18+3 diurnal time steps for up to 1 daily iterations, CNDCA is in Num dens., NFSP=54
C                                                CNDCA is NOT in COMMON

        REAL*8 CNDCA(S$,ntime+3,L$,Z$), CNDC0(S$), ohrat
        REAL*8 bsum, clsum, nsum, bxrat, clrat, nxrat
        REAL*8 brcl0, brono2a, brno2a, dbrcl, dbrono2, dbrno2
        REAL*8 clof, chfac, clono2a, clno2a, dclono2, dclno2
        REAL*8 noxf, xspecs, xnoy, xcly, xbry

        REAL*8 clofo, dclono2o, dclno2o, dbrono2o, dbrno2o, noxfo,dbrclo
        REAL*8 brclo, clyo, clsumo, noyo, nsumo, xspeco, shno3o


        SAVE 


C
C  AFTER FAST CHEMISTRY:  first save CNDC array from previous day into array CNDC0(S$), then update 
C    array CNDC(S$,ntime+3,L$,Z$) (in COMMON) with current day's values from array CNDCA(S$,ntime+3,L$,Z$)
C
C    then, need to rescale Bry, Cly, NOy families (this is NOT done in NWT2DD/JACOB), also 
C    sometimes get non-conservation of Cly when transporting everything separately
C    this happens in winter spring SH polar region and appears to be related to when HCL -> 0
C    so need to ensure that total family is the same as before fast chemistry (ie, reset to CN(31), CN(33)...
C   
C    Bry doesn't seem to have a problem, but adjust it anyway to be sure - do first
C      rescale X species here, subtracting out the change from the total family
C      ie, do Bry first, rescale BrCl subtract the change in BrCl from the Cly family
C                                                                                      CNDC(S$,ntime+3,L$,Z$)
C
C      if there are problems with rescaling (ie, scaling factor goes .LE. 0) 
C              then set to previous days value for that grid point (done at the end) using the CNDC0 array 
C
C
      DO 250 IK=1,Z$
      DO 250 IJ=1,L$
      DO 250 ITT=1,NTIME+3

         do 210 is=1,NFSP
             if (ISPMAP(is) .ne. 0) then
               CNDC0(ISPMAP(is)) = CNDC(ISPMAP(is),itt,ij,ik)
               CNDC(ISPMAP(is),itt,ij,ik) = CNDCA(ISPMAP(is),itt,ij,ik)
             endif
 210	 continue


C                                initialize IIBAD - if there are problems w/ rescaling, set IIBAD = 1
       iibad = 0


C   Bry -  save initial BrCl, BrONO2, BrNO2

       brcl0   = cndc(60,itt,ij,ik)
       brono2a = cndc(47,itt,ij,ik)
       brno2a  = cndc(79,itt,ij,ik)
    
       bsum = cndc(44,itt,ij,ik)+ cndc(45,itt,ij,ik)+ cndc(46,itt,ij,ik)  
     >      + cndc(47,itt,ij,ik)+ cndc(60,itt,ij,ik)+ cndc(68,itt,ij,ik) 
     >      + cndc(79,itt,ij,ik) + 2.*cndc(43,itt,ij,ik)

       bxrat = cn(48,ij,ik)/bsum

       cndc(43,itt,ij,ik) = cndc(43,itt,ij,ik)*bxrat
       cndc(44,itt,ij,ik) = cndc(44,itt,ij,ik)*bxrat
       cndc(45,itt,ij,ik) = cndc(45,itt,ij,ik)*bxrat
       cndc(46,itt,ij,ik) = cndc(46,itt,ij,ik)*bxrat
       cndc(47,itt,ij,ik) = cndc(47,itt,ij,ik)*bxrat
       cndc(60,itt,ij,ik) = cndc(60,itt,ij,ik)*bxrat
       cndc(68,itt,ij,ik) = cndc(68,itt,ij,ik)*bxrat
       cndc(79,itt,ij,ik) = cndc(79,itt,ij,ik)*bxrat

                                                   ! save change in X-specs, these will be applied to Clx,NOx
       dbrcl   = cndc(60,itt,ij,ik) - brcl0
       dbrono2 = cndc(47,itt,ij,ik) - brono2a
       dbrno2  = cndc(79,itt,ij,ik) - brno2a


C   *******************************************************************************
C
C  Cly - first apply opposite change in BrCl to the ClOx species only
C    so total Cly family is UNCHANGED due to BrCl rescaling
C

       clof = cndc(26,itt,ij,ik)+ cndc(27,itt,ij,ik)+ cndc(28,itt,ij,ik)
     >   + 2.*cndc(61,itt,ij,ik) + cndc(62,itt,ij,ik) 
     >   + 2.*cndc(64,itt,ij,ik) + cndc(65,itt,ij,ik)

       chfac = 1.d0 - dbrcl/clof

       if (chfac .gt. 0.) then
            cndc(26,itt,ij,ik) = cndc(26,itt,ij,ik)*chfac
            cndc(27,itt,ij,ik) = cndc(27,itt,ij,ik)*chfac
            cndc(28,itt,ij,ik) = cndc(28,itt,ij,ik)*chfac
            cndc(61,itt,ij,ik) = cndc(61,itt,ij,ik)*chfac
            cndc(62,itt,ij,ik) = cndc(62,itt,ij,ik)*chfac
            cndc(64,itt,ij,ik) = cndc(64,itt,ij,ik)*chfac
            cndc(65,itt,ij,ik) = cndc(65,itt,ij,ik)*chfac
       else
           dbrclo = dbrcl/m(ij,ik)*1.d9
           clofo  = clof/m(ij,ik)*1.d9

           write(1571,221) iyr, iday360, itt, ij, ik, 
     >            chfac, dbrclo, clofo
 221	   format('BrCl chfac le 0 : ', 5I4, 2x, 1P3D13.4)

           iibad = 1
       endif




C  now rescale Cly species to total family  excluding X-species BrCl,  as done in SOLVER
C      also save initial ClONO2, ClNO2


       clono2a = cndc(30,itt,ij,ik)
       clno2a  = cndc(63,itt,ij,ik)

       clsum = cndc(25,itt,ij,ik)+ cndc(26,itt,ij,ik)+cndc(27,itt,ij,ik)  
     >      + cndc(28,itt,ij,ik)+ cndc(29,itt,ij,ik)+ cndc(30,itt,ij,ik) 
     >   + 2.*cndc(61,itt,ij,ik) + cndc(62,itt,ij,ik) 
     > + cndc(63,itt,ij,ik) + 2.*cndc(64,itt,ij,ik) + cndc(65,itt,ij,ik)

       clrat = (cn(33,ij,ik) - cndc(60,itt,ij,ik))/clsum

       if (clrat .gt. 0.) then
           cndc(25,itt,ij,ik) = cndc(25,itt,ij,ik)*clrat
           cndc(26,itt,ij,ik) = cndc(26,itt,ij,ik)*clrat
           cndc(27,itt,ij,ik) = cndc(27,itt,ij,ik)*clrat
           cndc(28,itt,ij,ik) = cndc(28,itt,ij,ik)*clrat
           cndc(29,itt,ij,ik) = cndc(29,itt,ij,ik)*clrat
           cndc(30,itt,ij,ik) = cndc(30,itt,ij,ik)*clrat
           cndc(61,itt,ij,ik) = cndc(61,itt,ij,ik)*clrat
           cndc(62,itt,ij,ik) = cndc(62,itt,ij,ik)*clrat
           cndc(63,itt,ij,ik) = cndc(63,itt,ij,ik)*clrat
           cndc(64,itt,ij,ik) = cndc(64,itt,ij,ik)*clrat
           cndc(65,itt,ij,ik) = cndc(65,itt,ij,ik)*clrat
       else
           clyo   = cn(33,ij,ik)/m(ij,ik)*1.d9
           brclo  = cn(60,ij,ik)/m(ij,ik)*1.d9
           clsumo = clsum/m(ij,ik)*1.d9

           write(1572,222) iyr, iday360, itt, ij, ik, 
     >            clrat, clyo, brclo, clsumo
 222      format('IN UNFILLA, CLRAT le 0 : ', 5I4, 2x, 1P4D13.4)

           iibad = 1
       endif

                                                   ! save change in X-specs, these will be applied to NOx
       dclono2 = cndc(30,itt,ij,ik) - clono2a
       dclno2  = cndc(63,itt,ij,ik) - clno2a


C   *******************************************************************************
C
C  NOy - first apply opposite change in xspecs ClONO2, ClNO2, BrONO2, BrNO2 to the NOx species only
C    so total NOy family is UNCHANGED due to xspecs rescaling above
C

       noxf = cndc(5,itt,ij,ik) + cndc(6,itt,ij,ik) + cndc(7,itt,ij,ik)
     >   + 2.*cndc(8,itt,ij,ik) + cndc(9,itt,ij,ik) 

       chfac = 1.d0 - (dclono2 + dclno2 + dbrono2 + dbrno2)/noxf

       if (chfac .gt. 0.) then
           cndc(5,itt,ij,ik) = cndc(5,itt,ij,ik)*chfac
           cndc(6,itt,ij,ik) = cndc(6,itt,ij,ik)*chfac
           cndc(7,itt,ij,ik) = cndc(7,itt,ij,ik)*chfac
           cndc(8,itt,ij,ik) = cndc(8,itt,ij,ik)*chfac
           cndc(9,itt,ij,ik) = cndc(9,itt,ij,ik)*chfac
       else
           dclono2o = dclono2/m(ij,ik)*1.d9
           dclno2o  = dclno2/m(ij,ik)*1.d9
           dbrono2o = dbrono2/m(ij,ik)*1.d9
           dbrno2o  = dbrno2/m(ij,ik)*1.d9
           noxfo    = noxf/m(ij,ik)*1.d9

           write(1573,223) iyr, iday360, itt, ij, ik, 
     >            chfac, dclono2o, dclno2o, dbrono2o, dbrno2o, noxfo
 223	   format('ClONO2 chfac le 0 : ', 5I4, 2x, 1P6D13.4)

           iibad = 1
       endif



C  now rescale NOy species to total family  excluding X-species,  as done in SOLVER

       nsum = cndc(5,itt,ij,ik) + cndc(6,itt,ij,ik) + cndc(7,itt,ij,ik)
     >   + 2.*cndc(8,itt,ij,ik) + cndc(9,itt,ij,ik) + cndc(10,itt,ij,ik)
     >    + cndc(38,itt,ij,ik) + cndc(42,itt,ij,ik) + cndc(89,itt,ij,ik)


       xspecs = cndc(30,itt,ij,ik) + cndc(63,itt,ij,ik) 
     >        + cndc(47,itt,ij,ik) + cndc(79,itt,ij,ik) 


       nxrat = (cn(31,ij,ik) - cn(66,ij,ik) - xspecs)/nsum

       if (nxrat .gt. 0.) then 
           cndc(5,itt,ij,ik) = cndc(5,itt,ij,ik)*nxrat
           cndc(6,itt,ij,ik) = cndc(6,itt,ij,ik)*nxrat
           cndc(7,itt,ij,ik) = cndc(7,itt,ij,ik)*nxrat
           cndc(8,itt,ij,ik) = cndc(8,itt,ij,ik)*nxrat
           cndc(9,itt,ij,ik) = cndc(9,itt,ij,ik)*nxrat
           cndc(10,itt,ij,ik) = cndc(10,itt,ij,ik)*nxrat
           cndc(38,itt,ij,ik) = cndc(38,itt,ij,ik)*nxrat
           cndc(42,itt,ij,ik) = cndc(42,itt,ij,ik)*nxrat
           cndc(89,itt,ij,ik) = cndc(89,itt,ij,ik)*nxrat
       else
           noyo   = cn(31,ij,ik)/m(ij,ik)*1.d9
           shno3o = cn(66,ij,ik)/m(ij,ik)*1.d9
           xspeco = xspecs/m(ij,ik)*1.d9
           nsumo  = nsum/m(ij,ik)*1.d9

           write(1574,224) iyr, iday360, itt, ij, ik, 
     >            nxrat, noyo, shno3o, xspeco, nsumo
 224      format('IN UNFILLA, NXRAT le 0 : ', 5I4, 2x, 1P5D13.4)

           iibad = 1
       endif


C  if IIBAD = 1 (ie, if ANY problems occur from rescaling, then reset CNDC array to 
C              previous day for the current grid point - CNDC0(S$),   NFSP=54

       IF (iibad .EQ. 1) then
          DO 225 IS=1,NFSP
             if (ISPMAP(is) .ne. 0)
     >          CNDC(ISPMAP(is),itt,ij,ik) = CNDC0(ISPMAP(is))
 225      CONTINUE
       ENDIF


 250  CONTINUE



C  TO RE-SET tropospheric OH (which has been transported) with TRANSCOM climatology:
C    (IFIXOH = 1, set in CONTROL.dat)
C
C   reset diurnal cycle ratio of 24-hour average OH climatology/model OH
C    OHDAY(L$,Z$) is in #/cm3 for current day;    ITROP360(L$,360) (in COMMON) are the
C    tropopause FORTRAN indicies for the current Z$ grid, defined in TROPKZZ
C    CNDC(S$,ntime+3,L$,Z$) is the Diurnal avg, load model 24-hr avg OH into CN(90)
                              !  do smooth blend around tropopause
      IF (IFIXOH .eq. 1) then

        DO 260 IJ=1,L$
C                                       get tropopause index for current latitude
          iktrop = ITROP360(ij,iday360)
          if (ABS(LAT4(ij)) .ge. 27.  .and.  iktrop .gt. 14) iktrop = 14
          if (ABS(LAT4(ij)) .ge. 30.  .and.  iktrop .gt. 13) iktrop = 13
          if (ABS(LAT4(ij)) .ge. 35.  .and.  iktrop .gt. 12) iktrop = 12
          if (ABS(LAT4(ij)) .ge. 65.) iktrop = 11
cc                 if (ABS(LAT4(ij)) .ge. 75.  .and.  iktrop .gt. 10) iktrop = 10
cc                 if (ABS(LAT4(ij)) .ge. 80.  .and.  iktrop .gt.  9) iktrop = 9

          DO 265 IK=1,Z$
            C(90,ij,ik) = CNDC(13, NTIME+3, ij,ik)
           CN(90,ij,ik) = CNDC(13, NTIME+3, ij,ik)

           ohrat = 1.d0

         if (ik .eq. iktrop) ohrat = 0.5d0*
     >  (OHDAY(ij,ik) + CNDC(13, NTIME+3, ij,ik))/CNDC(13,NTIME+3,ij,ik)

         if (ik .lt. iktrop) ohrat = OHDAY(ij,ik)/CNDC(13,NTIME+3,ij,ik)

           DO 270 ITT=1,NTIME+3
 270	     CNDC(13,itt,ij,ik) = CNDC(13,itt,ij,ik)*ohrat

 265    CONTINUE
 260    CONTINUE

      ENDIF


C
C  LOAD in DIURNAL AVG radicals into CN and C(S$,L$,Z$) are REAL*8,  CNDC and C, CN are in Number density
C    the Diurnal avg is in index NTIME+3   ; 35 fast species are computed in the AER code

C   also load in DAYTIME and NITETIME avgs in array CDN(2,S$,L$,Z$);   DAYAVG=NTIME+1, NITAVG=NTIME+2

        DO 105 IK=1,Z$
        DO 105 IJ=1,L$
        DO 105 IS=1,NFSP
          if (ISPMAP(is) .ne. 0) then
             C(ISPMAP(is),ij,ik) = CNDC(ISPMAP(is), NTIME+3, ij,ik)
            CN(ISPMAP(is),ij,ik) = CNDC(ISPMAP(is), NTIME+3, ij,ik)

            CDN(1,ISPMAP(is),ij,ik) = CNDC(ISPMAP(is), NTIME+1, ij,ik)
            CDN(2,ISPMAP(is),ij,ik) = CNDC(ISPMAP(is), NTIME+2, ij,ik)
          endif
 105	CONTINUE


C
CMESO
C     load in original Noon time ratios for mesosphere, using updated families from SOLVER/TRANSPORT
C       except for Clx and Brx, which we use constant partitioning > 60 km (since the family is constant)
CMESO
CMESO         DO 800 ik=IKMES,Z$
CMESO         DO 800 ij=1,L$
C  Ox
CMESO           cn(1,ij,ik) = NOONRAT(1,ij,ik)*cn(39,ij,ik)
CMESO           cn(2,ij,ik) = NOONRAT(2,ij,ik)*cn(39,ij,ik)
CMESO           cn(4,ij,ik) = NOONRAT(4,ij,ik)*cn(39,ij,ik)
CMESO
C  HOx
CMESO           cn(12,ij,ik) = NOONRAT(12,ij,ik)*cn(74,ij,ik)
CMESO           cn(13,ij,ik) = NOONRAT(13,ij,ik)*cn(74,ij,ik)
CMESO           cn(14,ij,ik) = NOONRAT(14,ij,ik)*cn(74,ij,ik)
CMESO           cn(16,ij,ik) = NOONRAT(16,ij,ik)*cn(74,ij,ik)
CMESO
C  CHx
CMESO           cn(22,ij,ik) = NOONRAT(22,ij,ik)*cn(73,ij,ik)
CMESO           cn(23,ij,ik) = NOONRAT(23,ij,ik)*cn(73,ij,ik)
CMESO           cn(24,ij,ik) = NOONRAT(24,ij,ik)*cn(73,ij,ik)
CMESO
C  NOz
CMESO           cn(5,ij,ik) = NOONRAT(5,ij,ik)*cn(32,ij,ik)
CMESO           cn(6,ij,ik) = NOONRAT(6,ij,ik)*cn(32,ij,ik)
CMESO           cn(7,ij,ik) = NOONRAT(7,ij,ik)*cn(32,ij,ik)
CMESO           cn(8,ij,ik) = NOONRAT(8,ij,ik)*cn(32,ij,ik)
CMESO           cn(9,ij,ik) = NOONRAT(9,ij,ik)*cn(32,ij,ik)
CMESO           cn(38,ij,ik) = NOONRAT(38,ij,ik)*cn(32,ij,ik)
CMESO           cn(42,ij,ik) = NOONRAT(42,ij,ik)*cn(32,ij,ik)
CMESO
C  Clx                                                                    ! use constant partitioning > 60 km
CMESO         cn(25,ij,ik) = cn(25,ij,ikmes-1)/cn(33,ij,ikmes-1)*cn(33,ij,ik)
CMESO         cn(26,ij,ik) = cn(26,ij,ikmes-1)/cn(33,ij,ikmes-1)*cn(33,ij,ik)
CMESO         cn(27,ij,ik) = cn(27,ij,ikmes-1)/cn(33,ij,ikmes-1)*cn(33,ij,ik)
CMESO         cn(28,ij,ik) = cn(28,ij,ikmes-1)/cn(33,ij,ikmes-1)*cn(33,ij,ik)
CMESO         cn(29,ij,ik) = cn(29,ij,ikmes-1)/cn(33,ij,ikmes-1)*cn(33,ij,ik)
CMESO         cn(30,ij,ik) = cn(30,ij,ikmes-1)/cn(33,ij,ikmes-1)*cn(33,ij,ik)
CMESO         cn(61,ij,ik) = cn(61,ij,ikmes-1)/cn(33,ij,ikmes-1)*cn(33,ij,ik)
CMESO         cn(62,ij,ik) = cn(62,ij,ikmes-1)/cn(33,ij,ikmes-1)*cn(33,ij,ik)
CMESO         cn(63,ij,ik) = cn(63,ij,ikmes-1)/cn(33,ij,ikmes-1)*cn(33,ij,ik)
CMESO         cn(64,ij,ik) = cn(64,ij,ikmes-1)/cn(33,ij,ikmes-1)*cn(33,ij,ik)
CMESO
ccc           cn(25,ij,ik) = NOONRAT(25,ij,ik)*cn(33,ij,ik)
ccc           cn(26,ij,ik) = NOONRAT(26,ij,ik)*cn(33,ij,ik)
ccc           cn(27,ij,ik) = NOONRAT(27,ij,ik)*cn(33,ij,ik)
ccc           cn(28,ij,ik) = NOONRAT(28,ij,ik)*cn(33,ij,ik)
ccc           cn(29,ij,ik) = NOONRAT(29,ij,ik)*cn(33,ij,ik)
ccc           cn(30,ij,ik) = NOONRAT(30,ij,ik)*cn(33,ij,ik)
ccc           cn(61,ij,ik) = NOONRAT(61,ij,ik)*cn(33,ij,ik)
ccc           cn(62,ij,ik) = NOONRAT(62,ij,ik)*cn(33,ij,ik)
ccc           cn(63,ij,ik) = NOONRAT(63,ij,ik)*cn(33,ij,ik)
ccc           cn(64,ij,ik) = NOONRAT(64,ij,ik)*cn(33,ij,ik)
CMESO
C  Brx                                                                    ! use constant partitioning > 60 km
CMESO         cn(43,ij,ik) = cn(43,ij,ikmes-1)/cn(48,ij,ikmes-1)*cn(48,ij,ik)
CMESO         cn(44,ij,ik) = cn(44,ij,ikmes-1)/cn(48,ij,ikmes-1)*cn(48,ij,ik)
CMESO         cn(45,ij,ik) = cn(45,ij,ikmes-1)/cn(48,ij,ikmes-1)*cn(48,ij,ik)
CMESO         cn(46,ij,ik) = cn(46,ij,ikmes-1)/cn(48,ij,ikmes-1)*cn(48,ij,ik)
CMESO         cn(47,ij,ik) = cn(47,ij,ikmes-1)/cn(48,ij,ikmes-1)*cn(48,ij,ik)
CMESO         cn(60,ij,ik) = cn(60,ij,ikmes-1)/cn(48,ij,ikmes-1)*cn(48,ij,ik)
CMESO         cn(68,ij,ik) = cn(68,ij,ikmes-1)/cn(48,ij,ikmes-1)*cn(48,ij,ik)
CMESO
ccc           cn(43,ij,ik) = NOONRAT(43,ij,ik)*cn(48,ij,ik)
ccc           cn(44,ij,ik) = NOONRAT(44,ij,ik)*cn(48,ij,ik)
ccc           cn(45,ij,ik) = NOONRAT(45,ij,ik)*cn(48,ij,ik)
ccc           cn(46,ij,ik) = NOONRAT(46,ij,ik)*cn(48,ij,ik)
ccc           cn(47,ij,ik) = NOONRAT(47,ij,ik)*cn(48,ij,ik)
ccc           cn(60,ij,ik) = NOONRAT(60,ij,ik)*cn(48,ij,ik)
ccc           cn(68,ij,ik) = NOONRAT(68,ij,ik)*cn(48,ij,ik)
CMESO
C                                                           ! also update C array here for radicals
CMESO        DO 801 IS=1,35
CMESO 801	   C(ISPMAP(is),ij,ik) = CN(ISPMAP(is),ij,ik) 
CMESO
CMESO 800     CONTINUE
CMESO

CFAM
CFAMC   Sum up DIURNALLY AVERAGED family members here for transport for CHx and HOx - No chemistry 
CFAMC      load both the C and CN arrays here for the families, since these should be the same for
CFAMC      the long lived species until transport/solver
CFAMC  
CFAM
CFAM        DO 100 ik=1,Z$
CFAM	DO 100 ij=1,L$
CFAMC
CFAMC  Ox  -  C(39)  (O + O(1D) + O3)
CFAMC        
CFAM          C(39,ij,ik) = c(1,ij,ik) + c(2,ij,ik) + c(4,ij,ik) 
CFAM          CN(39,ij,ik) = C(39,ij,ik) 
CFAM
CFAM
CFAMC
CFAMC  CHx  -  C(73)  (CH3O2 + CH2O + CH3OOH)
CFAMC        
CFAM          C(73,ij,ik) = c(22,ij,ik) + c(23,ij,ik) + c(24,ij,ik) 
CFAM          CN(73,ij,ik) = C(73,ij,ik) 
CFAM
CFAM
CFAMc
CFAMC  HOx  -  C(74)  (H + OH + HO2 + 2*H2O2)
CFAMC        
CFAM          C(74,ij,ik) = c(12,ij,ik) + c(13,ij,ik) + c(14,ij,ik) 
CFAM     >                 + 2.*c(16,ij,ik)
CFAM          CN(74,ij,ik) = C(74,ij,ik) 
CFAM
CFAM
CFAMC
CFAMC  NOz  -  C(32)
CFAM
CFAM          CN(32,ij,ik) = cn(9,ij,ik) + cn(5,ij,ik) + cn(6,ij,ik)  
CFAM     >                 + cn(7,ij,ik) + 2.*cn(8,ij,ik) + cn(38,ij,ik) 
CFAM     >                 + cn(42,ij,ik) 
CFAM     >                 + cn(30,ij,ik) + cn(63,ij,ik)  + cn(47,ij,ik) 
CFAM
CFAM
CFAMC  Cly  -  CN(33)
CFAM
CFAM          CN(33,ij,ik) = cn(25,ij,ik) + cn(26,ij,ik) + cn(27,ij,ik)  
CFAM     >                 + cn(28,ij,ik) + cn(29,ij,ik) + cn(30,ij,ik) 
CFAM     >                 + 2.*cn(61,ij,ik) + cn(62,ij,ik) + cn(63,ij,ik)  
CFAM     >                 + 2.*cn(64,ij,ik) + cn(60,ij,ik)
CFAM
CFAM
CFAMC  Bry  -  CN(48)
CFAM
CFAM          CN(48,ij,ik) = 2.*cn(43,ij,ik) + cn(44,ij,ik) + cn(45,ij,ik)  
CFAM     >                 + cn(46,ij,ik) + cn(47,ij,ik) + cn(68,ij,ik) 
CFAM     >                 + cn(60,ij,ik)
CFAM
CFAM 100	CONTINUE
CFAM

cc      print *, '   '
cc      print *, '   '
cc      print *, ' NOz, Cly, Bry  at 30 km, EQ.  BEFORE, AFTER  FASTCHEM:'
cc      print *,   c(32,9,15), c(33,9,15), c(48,9,15)
cc      print *,   cn(32,9,15), cn(33,9,15), cn(48,9,15)
cc      print *, '   '
cc      print *, '   '



C  load in ratio of Noontime radical concentrations to the diurnal average of each constituent
C      use the DIURNALLY AVERAGES since this is what is transported
C      So the sum of the noontime ratios for a particular family can be > or < 1.
C
C      (most families DO NOT have a diurnal cycle, but CHx and HOx DO!! - these are summed ABOVE)
C      NOONRAT(S$,L$,Z$), CNDC(S$,ntime+3,L$,Z$) 
C    - except for O2(1D), which is just started from previous day's value
C
C
        DO 700 ik=1,Z$ 
        DO 700 ij=1,L$

C  Ox                                                  
           NOONRAT(1,ij,ik) = CNDC(1, NTIME, ij,ik)/CN(1,ij,ik)                  ! CN(39,ij,ik)
           NOONRAT(2,ij,ik) = CNDC(2, NTIME, ij,ik)/CN(2,ij,ik)                  ! CN(39,ij,ik)
           NOONRAT(4,ij,ik) = CNDC(4, NTIME, ij,ik)/CN(4,ij,ik)                  ! CN(39,ij,ik)

           NOONRAT(41,ij,ik)= CNDC(41, NTIME, ij,ik)/CN(41,ij,ik)           ! O2(1D) use ratio Noon/diur. avg

C  HOx
           NOONRAT(12,ij,ik) = CNDC(12, NTIME, ij,ik)/CN(12,ij,ik)
           NOONRAT(13,ij,ik) = CNDC(13, NTIME, ij,ik)/CN(13,ij,ik)
           NOONRAT(14,ij,ik) = CNDC(14, NTIME, ij,ik)/CN(14,ij,ik)
           NOONRAT(16,ij,ik) = CNDC(16, NTIME, ij,ik)/CN(16,ij,ik)


C  CHx
           NOONRAT(22,ij,ik) = CNDC(22, NTIME, ij,ik)/CN(22,ij,ik) 
           NOONRAT(23,ij,ik) = CNDC(23, NTIME, ij,ik)/CN(23,ij,ik) 
           NOONRAT(24,ij,ik) = CNDC(24, NTIME, ij,ik)/CN(24,ij,ik) 


C  NOz
           NOONRAT(5,ij,ik)  = CNDC(5, NTIME, ij,ik)/CN(5,ij,ik) 
           NOONRAT(6,ij,ik)  = CNDC(6, NTIME, ij,ik)/CN(6,ij,ik) 
           NOONRAT(7,ij,ik)  = CNDC(7, NTIME, ij,ik)/CN(7,ij,ik) 
           NOONRAT(8,ij,ik)  = CNDC(8, NTIME, ij,ik)/CN(8,ij,ik) 
           NOONRAT(9,ij,ik)  = CNDC(9, NTIME, ij,ik)/CN(9,ij,ik) 
           NOONRAT(10,ij,ik) = CNDC(10, NTIME, ij,ik)/CN(10,ij,ik) 
           NOONRAT(38,ij,ik) = CNDC(38, NTIME, ij,ik)/CN(38,ij,ik) 
           NOONRAT(42,ij,ik) = CNDC(42, NTIME, ij,ik)/CN(42,ij,ik) 
           NOONRAT(89,ij,ik) = CNDC(89, NTIME, ij,ik)/CN(89,ij,ik) 


C  Clx
           NOONRAT(25,ij,ik) = CNDC(25, NTIME, ij,ik)/CN(25,ij,ik)
           NOONRAT(26,ij,ik) = CNDC(26, NTIME, ij,ik)/CN(26,ij,ik)
           NOONRAT(27,ij,ik) = CNDC(27, NTIME, ij,ik)/CN(27,ij,ik)
           NOONRAT(28,ij,ik) = CNDC(28, NTIME, ij,ik)/CN(28,ij,ik)
           NOONRAT(29,ij,ik) = CNDC(29, NTIME, ij,ik)/CN(29,ij,ik)
           NOONRAT(30,ij,ik) = CNDC(30, NTIME, ij,ik)/CN(30,ij,ik)
           NOONRAT(61,ij,ik) = CNDC(61, NTIME, ij,ik)/CN(61,ij,ik)
           NOONRAT(62,ij,ik) = CNDC(62, NTIME, ij,ik)/CN(62,ij,ik)
           NOONRAT(63,ij,ik) = CNDC(63, NTIME, ij,ik)/CN(63,ij,ik)
           NOONRAT(64,ij,ik) = CNDC(64, NTIME, ij,ik)/CN(64,ij,ik)
           NOONRAT(65,ij,ik) = CNDC(65, NTIME, ij,ik)/CN(65,ij,ik)


C  Brx
           NOONRAT(43,ij,ik) = CNDC(43, NTIME, ij,ik)/CN(43,ij,ik)
           NOONRAT(44,ij,ik) = CNDC(44, NTIME, ij,ik)/CN(44,ij,ik)
           NOONRAT(45,ij,ik) = CNDC(45, NTIME, ij,ik)/CN(45,ij,ik)
           NOONRAT(46,ij,ik) = CNDC(46, NTIME, ij,ik)/CN(46,ij,ik)
           NOONRAT(47,ij,ik) = CNDC(47, NTIME, ij,ik)/CN(47,ij,ik)
           NOONRAT(60,ij,ik) = CNDC(60, NTIME, ij,ik)/CN(60,ij,ik)
           NOONRAT(68,ij,ik) = CNDC(68, NTIME, ij,ik)/CN(68,ij,ik)
           NOONRAT(79,ij,ik) = CNDC(79, NTIME, ij,ik)/CN(79,ij,ik)


C    and check to ensure that NOONRAT does NOT get EXTREMELY SMALL or else get NANs
C                             NOONRAT(S$,L$,Z$) is REAL*8 in COMMON 

           DO 705 is=1,S$
             IF (NOONRAT(is,ij,ik) .LT. 1.D-30) NOONRAT(is,ij,ik)=1.D-30
 705	   CONTINUE


CCtfam - also load in total Cly - CN(32);  and Bry - CN(58) - 
C        species for transport, use diurnal avgs:
C
C        include BrCl in Bry, OK to include ClNO2 and ClONO2 in Cly


          C(32,ij,ik) = c(25,ij,ik) + c(26,ij,ik) + c(27,ij,ik)  
     >                + c(28,ij,ik) + c(29,ij,ik) + c(30,ij,ik) 
     >                + 2.*c(61,ij,ik) + c(62,ij,ik) + c(63,ij,ik)
     >                + 2.*c(64,ij,ik) + c(65,ij,ik) 


          C(58,ij,ik) = 2.*c(43,ij,ik) + c(44,ij,ik) + c(45,ij,ik)  
     >                   + c(46,ij,ik) + c(47,ij,ik) + c(60,ij,ik)
     >                   + c(68,ij,ik) + c(79,ij,ik)


 700     CONTINUE


      RETURN
      END
