C
        SUBROUTINE PHOTHEAT(DTIMEB)

C
C  for coupled model, this computes heating rates consistent w/ photolysis, 
C    use xsect*rflux arrays computed in SPDR with INSOLATION correction applied in PDFINSOL
C    and updated diurnal variations of ozone, NO2, adjusted in UNFILLA.
C    these will be used in the call to the COUPLED model on the following day:  (EF, Sept 2008)
C
C    NOTE: for diurnal averaging, follow methodology from TIMEAV.f, 
C         SUNRIS and SUNSET are ALWAYS 12 and 7, respectively (at ALL latitudes, days of the year)
C

        include "com2d.h"
        include "com_aerd.h"


C  XSRF is common for O3, O2, NO2  cross sections*RFLUX for consistent photolysis/heating rates: 
C                         1=Ozone UV;  2=O3-vis;  3=O2;  4=NO2  -  XSRF is in photons/second/nm

        COMMON/CXSRF/XSRF(18,4,L$,Z$)

        COMMON/CPHEAT/PHEAT(5,L$,Z$), SORADHEAT(15,L$,Z$)

        DOUBLE PRECISION DTIMEB(L$,NTIME)
        REAL heatd(18,4,L$,Z$)
        REAL*8 x1, o3den, o2den, no2den, rhocp, XSRF


        SAVE

c                                 ! converts from milliWatts-sec to Watts-sec
        ffac = 1.9863e-13/1000.
        cpx = 1004.


C                             loop over grid point for ALL calculations
        do 1000 ik=1,Z$
        do 1000 ij=1,L$


C  ONLY FOR POLAR DAY IS THERE DAYLIGHT AT SUNSET/RIS TIMES, so set = 0 otherwise
c    SUNSET=7, SUNRIS=12 ALWAYS!!  SUNRIS, SUNSET, NTIME are in COM_AERD.H (ICLOCK=8,11 are 0.0)

            IF (TFD(ij) .LT. 1.0) then
               do 555 isss=1,4
                  xsrf(SUNSET,isss,ij,ik) = 0.
 555              xsrf(SUNRIS,isss,ij,ik) = 0.
            END IF

C
C  compute energy in Joule/sec, Integration over wavelength bins is already done in XSRF 
C     use REAL*8 CNDC(S$,18+3,L$,Z$), CN(S$,L$,Z$)  for O3, NO2, O2 number densities - in COMMON
C
C   XSRF(18,4,L$,Z$); 1=Ozone (UV);  2=O3-vis;  3=O2;  4=NO2;  XSRF is in photons/sec/nm
C
C   convert to heating rate (K/sec), ie, 1/(rho*Cp); load into heatd(18,4,L$,Z$)
C       ; 1=O3 UV;  2= O3 Vis;  3=O2;  4=NO2
C
C             ! rho in kg/cm^3, use TEMP(L$,Z$X), PRESS(Z$X) (in COMMON), include Cp

          rhocp = press(ik)*100./(287.*TEMP(ij,ik))/1.e6 * cpx


          do 500 itt = 1,NTIME

             o3den  = cndc(4,itt,ij,ik)
             no2den = cndc(6,itt,ij,ik)
             o2den  = cn(3,ij,ik)

C  Ozone UV:
             heatd(itt,1,ij,ik) = xsrf(itt,1,ij,ik)*o3den*ffac/rhocp

C  O3-Vis:
             heatd(itt,2,ij,ik) = xsrf(itt,2,ij,ik)*o3den*ffac/rhocp

C  O2:
             heatd(itt,3,ij,ik) = xsrf(itt,3,ij,ik)*o2den*ffac/rhocp

C  NO2:
             heatd(itt,4,ij,ik) = xsrf(itt,4,ij,ik)*no2den*ffac/rhocp

 500     CONTINUE

 1000    CONTINUE
                          ! end latitude-altitude loop


C
C  DIURNAL AVERAGING - use methodology in TIMEAV.f- only compute daytime values, nitetime values=0;  DTIMEB(L$,NTIME)
C   SUNRIS and SUNSET are ALWAYS 12 and 7, respectively (at ALL latitudes, days of the year)
C   heatd(18,4,L$,Z$) -> load into PHEAT(5,L$,Z$) (K/sec) in COMMON;  1=O3 UV;  2= O3 Vis;  3=O2;  4=NO2;  5=Total
C

        DO 700 ik=1,Z$
        DO 700 ij=1,L$
        DO 700 ihh=1,4

             htday = 0.

             do 710 it=1,SUNSET-1
                htday = htday + 
     >   (heatd(it,ihh,ij,ik) + heatd(it+1,ihh,ij,ik))*0.5*DTIMEB(ij,it)
 710    CONTINUE

             do 720 it=SUNRIS,NTIME-1
                htday = htday + 
     >   (heatd(it,ihh,ij,ik) + heatd(it+1,ihh,ij,ik))*0.5*DTIMEB(ij,it)
 720    CONTINUE

          PHEAT(ihh,ij,ik) = htday/86400.
 700    CONTINUE


        DO 770 ik=1,Z$
        DO 770 ij=1,L$
            PHEAT(5,ij,ik) = PHEAT(1,ij,ik) + PHEAT(2,ij,ik) 
     >                     + PHEAT(3,ij,ik) + PHEAT(4,ij,ik)
 770    CONTINUE


	RETURN
	END
