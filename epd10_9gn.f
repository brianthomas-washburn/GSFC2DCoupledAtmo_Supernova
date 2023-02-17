C
	SUBROUTINE EPD10(PSICR, PSICI, UBW1, THC, IDOY360, DELF10)

C  PROGRAM TO Generate u', v', T', W' from pwave GEOP fields
C          then compute EPFlux div, eddy heating terms
C
C  PSICR, PSICI are in cm^2/sec^2,   dz is in cm in COMMONC.INC
C       THC is current zonal mean model theta 
C
!        USE degree_trig

	INCLUDE 'PARAM.INC'
	INCLUDE 'COMMONC.INC'


        COMMON/CTXX/ txx(72,N$,M$)

        COMMON/CEHFC/ ehfc(NP$,MP$), epzz(NP$,MP$)


        REAL PSICR(N$,M$,MMAX$), PSICI(N$,M$,MMAX$), ggg(72,n$,m$)
        REAL ggx(72,N$,M$), ggy(72,N$,M$), ggz(72,N$,M$), ubw1(N$,M$)

        REAL uiter(10,72,N$,M$), viter(10,72,N$,M$)
        REAL fff1(N$,M$), fff2(N$,M$), dudx(72), dvdx(72)
        REAL uppp(72,N$,M$), vppp(72,N$,M$), tv(N$,M$), uv(N$,M$)

        REAL thxx(72,N$,M$), vxx(72,N$,M$), THC(N$,M$)

        REAL*8 ggg1(N$,M$), ggzz(N$,M$), ggyy(N$,M$)
        REAL*8 ggg2(74,M$), ggxx(74,M$), ttemp(N$,M$), ttempz(N$,M$)
        REAL*8 THC8(N$,M$), THCY(N$,M$), THCZ(N$,M$), delx, dym, dzm
        REAL*8 vth(N$,M$), tterm(N$,M$), ttermz(N$,M$), ggxz(74,M$)
        REAL*8 wth(N$,M$), wup(N$,M$), wthz(N$,M$), wupz(N$,M$)
        REAL*8 ubw8(N$,M$), uby(N$,M$), ubz(N$,M$)
        REAL*8 term1(N$,M$), term2(N$,M$), term1y(N$,M$), term2z(N$,M$)

        REAL ehfc1(N$,M$), delxo(N$), gxz(72,N$,M$)
        REAL xxn2(N$,M$), delf10(N$,M$)
        REAL wpp(72,N$,M$), ehfzz(N$,M$), epzz1(N$,M$)


        do 100 ik=1,M$
        do 100 ij=1,N$
        do 100 ix=1,72        
C                                  first initialize, ggg is in m^2/sec^2,  ggg is PHI'
           ggg(ix,ij,ik) = 0.0
           xlon = (ix-1)*5.

           do 200 ih=1,MMAX$
              ggg(ix,ij,ik) = ggg(ix,ij,ik) + 
     >                   (-PSICI(ij,ik,ih)*SIND(360./360.*ih*xlon)
     >                    +PSICR(ij,ik,ih)*COSD(360./360.*ih*xlon))/1.e4
 200       CONTINUE
 100    CONTINUE
       

C  compute T' and theta',   dym, dzm are in meters,  NOTE: DERV4 routine needs values in REAL*8 

        dym = dy/100.d0
        dzm = dz/100.d0
        hhr = 7000./287.
        cvtt = 3.14159/180.

 
        do 400 ix=1,72

          do 405 ik=1,M$
          do 405 ij=1,N$        
 405         ggg1(ij,ik) = ggg(ix,ij,ik)


          CALL DERV4(1, ggg1, ggzz, dzm, N$, M$, 2)

          do 407 ik=1,M$
          do 407 ij=1,N$
 407	     ggz(ix,ij,ik) = ggzz(ij,ik)


          CALL DERV4(1, ggg1, ggyy, dym, N$, M$, 1)

          do 707 ik=1,M$
          do 707 ij=1,N$
 707	     ggy(ix,ij,ik) = ggyy(ij,ik)

 400    CONTINUE


C   DERV4 routine needs values in REAL*8 - ggg2, ggxx(74,M$) -> ggx(72,N$,M$)
C       C(N$) is COS(yp), CF(N$) is f, both in COMMONC.INC,  delx is in meters, vxx is m/sec

        do 500 ij=1,N$

          delx = 5.d0*6371000.*cvtt*C(ij)
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
 507	     ggx(ix,ij,ik) = ggxx(ix+1,ik)


C   get d/dx(d/dz) Phi' - ggxz(74,M$) is in 1/sec^2
C 
          CALL DERV4(1, ggxx, ggxz, dzm, 74, M$, 2)

          do 607 ik=1,M$
          do 607 ix=1,72
 607	     gxz(ix,ij,ik) = ggxz(ix+1,ik)

 500    CONTINUE



C  get du/dy, dudz - ubw1(N$,M$) is in cm/sec, ubw8(N$,M$) is m/sec, uby, ubz(N$,M$) all real*8

          do 510 ik=1,M$
          do 510 ij=1,N$
 510	     ubw8(ij,ik) = ubw1(ij,ik)/100.

          CALL DERV4(1, ubw8, uby, dym, N$, M$, 1)
          CALL DERV4(1, ubw8, ubz, dzm, N$, M$, 2)


C   fff1(N$,M$), fff2(N$,M$) per Randel, 1987,  CF(N$) is f

          do 511 ik=1,M$
          do 511 ij=1,N$
       	     utana = ubw8(ij,ik)/6371000.*TAN(yp(ij)*cvtt)

             fff1(ij,ik) = cf(ij) - uby(ij,ik) + utana
             fff2(ij,ik) = cf(ij) + 2.*utana
 511	  CONTINUE



C  get u', v' - iterate, starting w/ geostrophic values, uiter(10,72,N$,M$), viter(10,72,N$,M$)
C                             ubw8(N$,M$), CF(N$) is f;   ggx, ggy, ggz(72,N$,M$) are in m/sec^2
         do 517 ik=1,M$
         do 517 ij=1,N$
         do 517 ix=1,72

            dlx = 2.*5.d0*6371000.*cvtt*C(ij)

            uiter(1,ix,ij,ik) = -ggy(ix,ij,ik)/cf(ij)
            viter(1,ix,ij,ik) =  ggx(ix,ij,ik)/cf(ij)

            do 530 ii=2,4

               dudx(1) = (uiter(ii-1,2,ij,ik)-uiter(ii-1,72,ij,ik))/dlx
               dvdx(1) = (viter(ii-1,2,ij,ik)-viter(ii-1,72,ij,ik))/dlx

               dudx(72) = (uiter(ii-1,1,ij,ik)-uiter(ii-1,71,ij,ik))/dlx
               dvdx(72) = (viter(ii-1,1,ij,ik)-viter(ii-1,71,ij,ik))/dlx


            if (ix .ge. 2  .and. ix .le. 71) then
          dudx(ix) = (uiter(ii-1,ix+1,ij,ik)-uiter(ii-1,ix-1,ij,ik))/dlx
          dvdx(ix) = (viter(ii-1,ix+1,ij,ik)-viter(ii-1,ix-1,ij,ik))/dlx
            endif

            viter(ii,ix,ij,ik) = 
     >                (ggx(ix,ij,ik) + ubw8(ij,ik)*dudx(ix))/fff1(ij,ik)
                    
            uiter(ii,ix,ij,ik) = 
     >               (-ggy(ix,ij,ik) - ubw8(ij,ik)*dvdx(ix))/fff2(ij,ik)
 530    CONTINUE
 517	CONTINUE


c  load in winds from final iteration - uppp(72,N$,M$), vppp(72,N$,M$)
C        use geostrophic for < +-5 degrees
C

         do 522 ik=1,M$
         do 522 ij=1,N$
         do 522 ix=1,72
            uppp(ix,ij,ik) = uiter(4,ix,ij,ik)
            vppp(ix,ij,ik) = viter(4,ix,ij,ik)

            if (ABS(yp(ij)) .le. 5.) then
               uppp(ix,ij,ik) = uiter(1,ix,ij,ik)
               vppp(ix,ij,ik) = viter(1,ix,ij,ik)
            endif
 522	  CONTINUE


CC                       ttempz(N$,M$), xxn2(N$,M$)
        do 701 ik=1,M$
        do 701 ij=1,N$
 701	     ttemp(ij,ik) = thc(ij,ik)/TTOTH(ik)


          CALL DERV4(1, ttemp, ttempz, dzm, N$, M$, 2)


        do 702 ik=1,M$
        do 702 ij=1,N$
 702	   xxn2(ij,ik)=287./7000.*(ttempz(ij,ik) + 
     >                 .286*ttemp(ij,ik)/7000.)


C
C   get EP flux divergence:  tv(N$,M$), uv(N$,M$), ggz(72,N$,M$) is in m/sec^2
C                            term1(N$,M$), term2(N$,M$) are REAL*8

          do 617 ik=1,M$
          do 617 ij=1,N$

             tv(ij,ik) = 0.
             uv(ij,ik) = 0.

             do 527 ix=1,72
 527		tv(ij,ik) = tv(ij,ik) + vppp(ix,ij,ik)*ggz(ix,ij,ik)/72.

             do 537 ix=1,72
 537		uv(ij,ik) = uv(ij,ik) + vppp(ix,ij,ik)*uppp(ix,ij,ik)/72.


             term1(ij,ik) = ubz(ij,ik)*tv(ij,ik)/xxn2(ij,ik) - uv(ij,ik)
             term2(ij,ik) = fff1(ij,ik)*tv(ij,ik)/xxn2(ij,ik)
 617      CONTINUE


          CALL DERV4(1, term1, term1y, dym, N$, M$, 1)
          CALL DERV4(1, term2, term2z, dzm, N$, M$, 2)


          do 717 ik=1,M$
          do 717 ij=1,N$
             delf10(ij,ik) = -2.*TAN(yp(ij)*cvtt)/6371000.*term1(ij,ik)
     >             + term1y(ij,ik) + term2z(ij,ik) - term2(ij,ik)/7000.
 717	  CONTINUE





C  get paramertized w' (-UBAR/N^2*Phi'zx ;AHL, p.181) and eddy heating, EP flux contributions:
C        interpolate UBX (cm/sec) to proper grid, compute N^2 -  THC(N$,M$) -> ttemp(N$,M$)
C        ttempz(N$,M$), xxn2(N$,M$)

c        CALL BINTERP(ypp, NP$, zpp, MP$, UBX, YP, N$, ZP, M$, 0, 0, UBB)



C   wpp(72,N$,M$), uppp(72,N$,M$), wup, wth(N$,M$), ggy(72,N$,M$) is in m/sec^2
C      gxz(72,N$,M$) is in 1/sec^2; convert UBB to m/sec, wpp, uppp are in m/sec
C
c        do 710 ik=1,M$
c        do 710 ij=1,N$

c          do 712 ix=1,72
c             uppp(ix,ij,ik) = -ggy(ix,ij,ik)/cf(ij)
c            wpp(ix,ij,ik) = -(ubb(ij,ik)/100.)/xxn2(ij,ik)*gxz(ix,ij,ik)
c 712	  CONTINUE


c          wth(ij,ik) = 0.d0
c          do 717 ix=1,72
c 717	     wth(ij,ik) = wth(ij,ik) + wpp(ix,ij,ik)*thxx(ix,ij,ik)/72.


c          wup(ij,ik) = 0.d0
c          do 727 ix=1,72
c 727	     wup(ij,ik) = wup(ij,ik) + wpp(ix,ij,ik)*uppp(ix,ij,ik)/72.
c
c 710     CONTINUE

C   wth is in K-m/sec;  wthz is in K/sec;  wup is in m^2/sec^2, wupz is in m/sec^2


c         CALL DERV4(1, wth, wthz, dzm, N$, M$, 2)
c         CALL DERV4(1, wup, wupz, dzm, N$, M$, 2)


C   ehfzz(N$,M$) is K/sec; epzz1(N$,M$) is converted to cm/sec^2
C   NOTE: ehfzz term here is NEGLIGIBLE - ~+-1.e-13 K/sec or 1.e-8 K/day (in the machine noise)
C   and the patterns just look like noise which would indicate that this term is ZERo
C   in the QG Pwave system. The EP flux term (w'u') looks reasonable, so maybe it should be used?


c        do 730 ik=1,M$
c        do 730 ij=1,N$
c            ehfzz(ij,ik) = wth(ij,ik)/7000. - wthz(ij,ik)
c 730	     epzz1(ij,ik) = (wup(ij,ik)/7000. - wupz(ij,ik))*100.


C  Interpolate EP flux to extended dynamics grid - epzz(NP$,MP$) in cm/sec^2

c         CALL BINTERP(yp, N$, zp, M$, EPZZ1,
c    >                 YPP, NP$, ZPP, MP$, 0, 0, EPZZ)


 
ccccc        if (idoy360 .eq. 30  .or.  idoy360 .eq. 60) then
ccccc             write (446) N$, M$, np$, mp$, dym, dzm, idoy360
ccccc             write (446) yp, ypp, zp, zpp, c, cf, delxo
cccc
ccccc             write (446) ubw8, uby, fff1, fff2, uiter, viter
ccccc             write (446) uppp, vppp, ttemp, xxn2, tv, uv, ggz
ccccc             write (446) term1, term2, term1y, term2z, delf10
ccccc        endif


	RETURN
	END
