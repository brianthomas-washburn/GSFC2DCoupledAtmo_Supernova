
	SUBROUTINE PDFINSOL(ICLOCK)


C  This routine takes the zonal mean J-coefficients for the current ICLOCK computed
C    on the standard model grid (in SPDR), JX(PH$,L$,Z$), and interpolates exponentially
C    to a 1 degree latitude grid, and then applies the latitudinal probability
C    distribution obtained from the offline parcel model for each model latitude,
C    to get an array of corrected J-coefficients for input into the AER FASTCHEMISTRY in JMAP
C
C   Also apply correction for the heating rates in the same way - XSRFU(4,L$,Z$)
C 
C   base calculations now done in separate routine - PDFINT
C

       include "com2d.h"


       REAL*8 JX, JZ, XSRFU, XSRF, xlsum, JJ0(L$), J181(181), XJINT
       REAL*8 JX39, JZ39, J39

       COMMON/JALLTIME/JX(PH$,L$,Z$)
       COMMON/JCORRECT/JZ(18,PH$,L$,Z$)

       COMMON/CXSRFU/XSRFU(4,L$,Z$)
       COMMON/CXSRF/XSRF(18,4,L$,Z$)

ccccccc       COMMON/CJINT/XJINT(25,50,L$,Z$)

C  wavelength dependent J-arrays (for diagnostics ONLY)

       COMMON/CJX39/JX39(PH$,5,L$,Z$), JZ39(18,PH$,5,L$,Z$),
     >               J39(PH$,5,L$,Z$)


       SAVE

C                                first do J's:  JX(PH$,L$,Z$) -> JZ(18,PH$,L$,Z$)
       DO 100 ik=1,Z$
       DO 100 ipp=1,PH$

            do 101 ij=1,L$
 101           JJ0(ij) = JX(ipp,ij,ik)

C                                      ! interpolation to 181 latitudes now done in PDFINT
            CALL PDFINT(JJ0, J181)

C
C  now with J181(181), apply prob distribution,  LPROB(181,L$,Z$) to get CORRECTED  JZ(18,PH$,L$,Z$)
C                                                also use indicies ILP1(L$,Z$), ILP2(L$,Z$)

            do 500 ij=1,L$
               xlsum = 0.d0

               do 640 ilf=ILP1(ij,ik), ILP2(ij,ik)
     		  xlsum = xlsum + J181(ilf)*LPROB(ilf,ij,ik)
 640	       CONTINUE

               JZ(ICLOCK,ipp,ij,ik) = xlsum
 500	    CONTINUE
c                               ! end  PH$, ik loop
 100   CONTINUE



C                       wavelength dependent J's:  JX39(PH$,5,L$,Z$) -> JZ39(18,PH$,5,L$,Z$)
       DO 200 ik=1,Z$
       DO 200 ilx=1,5
       DO 200 ipp=1,PH$

            do 201 ij=1,L$
 201           JJ0(ij) = JX39(ipp,ilx,ij,ik)

C                                      ! interpolation to 181 latitudes now done in PDFINT
            CALL PDFINT(JJ0, J181)

C
C  now with J181(181), apply prob distribution,  LPROB(181,L$,Z$) to get CORRECTED  JZ(18,PH$,L$,Z$)
C                                                also use indicies ILP1(L$,Z$), ILP2(L$,Z$)

            do 510 ij=1,L$
               xlsum = 0.d0

               do 650 ilf=ILP1(ij,ik), ILP2(ij,ik)
     		  xlsum = xlsum + J181(ilf)*LPROB(ilf,ij,ik)
 650	       CONTINUE

               JZ39(ICLOCK,ipp,ilx,ij,ik) = xlsum
 510	    CONTINUE
c                               ! end  PH$, ik loop
 200   CONTINUE



C                             now for heating rates: XSRFU(4,L$,Z$) -> XSRF(18,4,L$,Z$)
       DO 700 ik=1,Z$
       DO 700 ipp=1,4

            do 701 ij=1,L$
 701           JJ0(ij) = XSRFU(ipp,ij,ik)

C                                      ! interpolation to 181 latitudes now done in PDFINT
            CALL PDFINT(JJ0, J181)

C
C   now with J181(181), apply prob distribution,  LPROB(181,L$,Z$) to get CORRECTED XSRF(18,4,L$,Z$)
C                                                 also use indicies ILP1(L$,Z$), ILP2(L$,Z$)

            do 600 ij=1,L$
               xlsum = 0.d0

               do 645 ilf=ILP1(ij,ik), ILP2(ij,ik)
     		  xlsum = xlsum + J181(ilf)*LPROB(ilf,ij,ik)
 645	       CONTINUE

               XSRF(ICLOCK,ipp,ij,ik) = xlsum
 600	    CONTINUE
c                               ! end  ipp, ik loop
 700   CONTINUE



C        if (iday360 .eq. 180) then
c        if (iday360 .eq. 78) then
Cc          write (117) L$, decd, laters, latern
Cc          write (117) LAT181, LAT, ijind
Cc          write (117) JX, JZ, XJINT
C          write (117) ILP1, ILP2, LPROB
C        endif

       RETURN
       END
