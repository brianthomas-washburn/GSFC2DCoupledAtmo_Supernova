C
	SUBROUTINE TEMPLON(PSICR, PSICI, THCE, UBX, IDOY360)

C  PROGRAM TO Generate longitudinal temperature array from pwave GEOP fields
C
C  PSICR, PSICI are in cm^2/sec^2,   dz is in cm in COMMONC.INC
C       THCE is current zonal mean model theta on extended grid
C
!	USE degree_trig

	INCLUDE 'PARAM.INC'
	INCLUDE 'COMMONC.INC'

        COMMON/CTXX/ txx(72,N$,M$)

        COMMON/CEHFC/ ehfc(NP$,MP$), epzz(NP$,MP$)


        REAL PSICR(N$,M$,MMAX$), PSICI(N$,M$,MMAX$), ggg(72,n$,m$)
        REAL thxx(72,N$,M$), vxx(72,N$,M$), THCE(NP$,MP$), THC(N$,M$)

        REAL*8 ggg1(N$,M$), ggzz(N$,M$), ggyy(N$,M$)
        REAL*8 ggg2(74,M$), ggxx(74,M$), ttemp(N$,M$), ttempz(N$,M$)
        REAL*8 THC8(N$,M$), THCY(N$,M$), THCZ(N$,M$), delx, dym, dzm
        REAL*8 vth(N$,M$), tterm(N$,M$), ttermz(N$,M$), ggxz(74,M$)
        REAL*8 wth(N$,M$), wup(N$,M$), wthz(N$,M$), wupz(N$,M$)

        REAL ehfc1(N$,M$), delxo(N$), gxz(72,N$,M$), ggy(72,N$,M$)
        REAL ubx(NP$,MP$), ubb(N$,M$), xxn2(N$,M$)
        REAL wpp(72,N$,M$), uppp(72,N$,M$), ehfzz(N$,M$), epzz1(N$,M$)


        do 100 ik=1,M$
        do 100 ij=1,N$
        do 100 ix=1,72        
C                                  first initialize, ggg is in m^2/sec^2,  ggg is PHI'
           ggg(ix,ij,ik) = 0.0
           xlon = (ix-1)*5.

           do 200 ih=1,MMAX$
              ggg(ix,ij,ik) = ggg(ix,ij,ik) + 
     >                   (-PSICI(ij,ik,ih)*SIND(360./360.*ih*xlon) +  
     >                     PSICR(ij,ik,ih)*COSD(360./360.*ih*xlon))/1.e4
 200       CONTINUE
 100    CONTINUE
       

C  compute T' and theta',   dym, dzm are in meters,  NOTE: DERV4 routine needs values in REAL*8 

        dym = dy/100.d0
        dzm = dz/100.d0
        hhr = 7000./287.
 
        do 400 ix=1,72

          do 405 ik=1,M$
          do 405 ij=1,N$        
 405         ggg1(ij,ik) = ggg(ix,ij,ik)

          CALL DERV4(1, ggg1, ggzz, dzm, N$, M$, 2)

          do 407 ik=1,M$
          do 407 ij=1,N$
             txx(ix,ij,ik) = hhr*ggzz(ij,ik)
             thxx(ix,ij,ik) = txx(ix,ij,ik)*TTOTH(ik)
 407      CONTINUE


          CALL DERV4(1, ggg1, ggyy, dym, N$, M$, 1)

          do 707 ik=1,M$
          do 707 ij=1,N$
 707	     ggy(ix,ij,ik) = ggyy(ij,ik)

 400    CONTINUE


C  compute v';  DERV4 routine needs values in REAL*8 - ggg2, ggxx(74,M$)
C       C(N$) is COS(yp), CF(N$) is f, both in COMMONC.INC,  delx is in meters, vxx is m/sec

        do 500 ij=1,N$

          delx = 5.d0*6371000.*3.14159/180.*C(ij)
          delxo(ij) = delx

          do 502 ik=1,M$
             ggg2(1,ik) = ggg(72,ij,ik)
             ggg2(74,ik) = ggg(1,ij,ik)

          do 505 ix=1,72
 505         ggg2(ix+1,ik) = ggg(ix,ij,ik)
 502	  CONTINUE

          CALL DERV4(1, ggg2, ggxx, delx, 74, M$, 1)

          do 507 ik=1,M$
          do 507 ix=1,72
 507	     vxx(ix,ij,ik) = ggxx(ix+1,ik)/cf(ij)


C   get d/dx(d/dz) Phi' - ggxz(74,M$) is in 1/sec^2
C 
          CALL DERV4(1, ggxx, ggxz, dzm, 74, M$, 2)

          do 607 ik=1,M$
          do 607 ix=1,72
 607	     gxz(ix,ij,ik) = ggxz(ix+1,ik)

 500    CONTINUE



C  Interpolate THCE to interior dynamics grid, compute dth/dy, dthdz, dym, dzm are in meters

          CALL BINTERP(ypp, NP$, zpp, MP$, THCE,
     >                 YP, N$, ZP, M$, 0, 1, THC)

C                            convert to real*8 for DERV4
          do 417 ik=1,M$
          do 417 ij=1,N$
 417	     thc8(ij,ik) = thc(ij,ik)


          CALL DERV4(1, thc8, thcy, dym, N$, M$, 1)
          CALL DERV4(1, thc8, thcz, dzm, N$, M$, 2)


C  get v'th'-bar  ; vth(N$,M$), tterm(N$,M$) is K-m/sec, THCY(N$,M$), THCZ(N$,M$) are both in K/m

          do 517 ik=1,M$
          do 517 ij=1,N$

             sum = 0.
             do 527 ix=1,72
 527		sum = sum + vxx(ix,ij,ik)*thxx(ix,ij,ik)

             vth(ij,ik) = sum/72.
 
             tterm(ij,ik) = vth(ij,ik)*thcy(ij,ik)/thcz(ij,ik)
 517	  CONTINUE

c                                                        ttermz(N$,M$) is in K/sec (theta heating)
          CALL DERV4(1, tterm, ttermz, dzm, N$, M$, 2)


          do 617 ik=1,M$
          do 617 ij=1,N$
 617	     ehfc1(ij,ik) = tterm(ij,ik)/7000. - ttermz(ij,ik)


C  Interpolate to extended dynamics grid - ehfc(NP$,MP$) in K/sec (theta heating)

          CALL BINTERP(yp, N$, zp, M$, EHFC1,
     >                 YPP, NP$, ZPP, MP$, 0, 0, EHFC)




C  get paramertized w' (-UBAR/N^2*Phi'zx ;AHL, p.181) and eddy heating, EP flux contributions:
C        interpolate UBX (cm/sec) to proper grid, compute N^2 -  THC(N$,M$) -> ttemp(N$,M$)
C        ttempz(N$,M$), xxn2(N$,M$)

        CALL BINTERP(ypp, NP$, zpp, MP$, UBX, YP, N$, ZP, M$, 0, 0, UBB)


        do 701 ik=1,M$
        do 701 ij=1,N$
 701	     ttemp(ij,ik) = thc(ij,ik)/TTOTH(ik)


          CALL DERV4(1, ttemp, ttempz, dzm, N$, M$, 2)


        do 702 ik=1,M$
        do 702 ij=1,N$
 702	   xxn2(ij,ik)=287./7000.*(ttempz(ij,ik) 
     >                 + .286*ttemp(ij,ik)/7000.)


C   wpp(72,N$,M$), uppp(72,N$,M$), wup, wth(N$,M$), ggy(72,N$,M$) is in m/sec^2
C      gxz(72,N$,M$) is in 1/sec^2; convert UBB to m/sec, wpp, uppp are in m/sec
C
        do 710 ik=1,M$
        do 710 ij=1,N$

          do 712 ix=1,72
             uppp(ix,ij,ik) = -ggy(ix,ij,ik)/cf(ij)
            wpp(ix,ij,ik) = -(ubb(ij,ik)/100.)/xxn2(ij,ik)*gxz(ix,ij,ik)
 712	  CONTINUE


          wth(ij,ik) = 0.d0
          do 717 ix=1,72
 717	     wth(ij,ik) = wth(ij,ik) + wpp(ix,ij,ik)*thxx(ix,ij,ik)/72.


          wup(ij,ik) = 0.d0
          do 727 ix=1,72
 727	     wup(ij,ik) = wup(ij,ik) + wpp(ix,ij,ik)*uppp(ix,ij,ik)/72.

 710     CONTINUE

C   wth is in K-m/sec;  wthz is in K/sec;  wup is in m^2/sec^2, wupz is in m/sec^2


         CALL DERV4(1, wth, wthz, dzm, N$, M$, 2)
         CALL DERV4(1, wup, wupz, dzm, N$, M$, 2)


C   ehfzz(N$,M$) is K/sec; epzz1(N$,M$) is converted to cm/sec^2
C   NOTE: ehfzz term here is NEGLIGIBLE - ~+-1.e-13 K/sec or 1.e-8 K/day (in the machine noise)
C   and the patterns just look like noise which would indicate that this term is ZERo
C   in the QG Pwave system. The EP flux term (w'u') looks reasonable, so maybe it should be used?


        do 730 ik=1,M$
        do 730 ij=1,N$
            ehfzz(ij,ik) = wth(ij,ik)/7000. - wthz(ij,ik)
 730	     epzz1(ij,ik) = (wup(ij,ik)/7000. - wupz(ij,ik))*100.


C  Interpolate EP flux to extended dynamics grid - epzz(NP$,MP$) in cm/sec^2

          CALL BINTERP(yp, N$, zp, M$, EPZZ1,
     >                 YPP, NP$, ZPP, MP$, 0, 0, EPZZ)


 
ccc        if (idoy360 .eq. 30  .or.  idoy360 .eq. 60) then
ccc             write (444) N$, M$, np$, mp$, dym, dzm, idoy360
ccc             write (444) yp, ypp, zp, zpp, c, cf, delxo
ccc             write (444) txx, thxx, vxx, thce, thc, thcy, thcz
ccc             write (444) vth, tterm, ttermz, ehfc1, ehfc
ccc
ccc             write (444) ubx, ubb, ttemp, ttempz, xxn2
ccc             write (444) uppp, wpp, wth, wthz, wup, wupz, ehfzz, epzz1
ccc        endif


	RETURN
	END
