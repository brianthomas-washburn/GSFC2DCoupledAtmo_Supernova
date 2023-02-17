      SUBROUTINE GETRAD(YIN, ZIN)
c      SUBROUTINE GETRAD(flag1)
 
      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONR.INC'
      
      common /pscpass/ psctotal(N$,45),pscdyn(N$,M$)
 
      logical flag1
ccelf      include 'test_value.inc'
ccelf      logical nantest
      
      dimension temphr(91),tempcr(91),d4(m$),d5(m$),tempsc(45)
      INTEGER IDAYS(13) !, NS$, MS$
      real yin(ns$),zin(ms$)
      DATA IDAYS/1,32,60,91,121,152,182,213,244,274,305,335,365/


C     FIRST GET TEMPERATURE FIELD ONTO RADIATION COORDINATES

       MZTOT=MZRAD+MZRAD+1
 8000 FORMAT( 'IN GETRAD , MZRAD,MZTOT ',I4,1X,I4)
c      WRITE(6,8000) MZRAD,MZTOT
 
      do 5 j=1,n$
      DO 6 K=1,M$
6     d3(K)=T(J,K)
      CALL INTER(TEMPhR,ZJRREV,MZTOT,D3,ZP,M$)
      DO 7 K=1,MZTOT
7        TTOT(J,K)=TEMPhR(MZTOT+1-K) !REVERSE ORDER OF TEMPS
5     CONTINUE
 
c      print *,'before radcal'
      CALL RADCAL(DECL, YIN, ZIN)
c      CALL RADCAL(DECL,flag1)
c      print *,'after radcal'
 
C   REVERSE BACK AND INTERPOLATE TO MRS GRID

      daylx=1./dayl
 
      DO 50 J=1,N$
      DO 30 K=1,MZRAD
      TEMPCR(K)=COOLJR(J,MZRAD+1-K)*daylx

c       nantest=test_nan( tempcr(k) )
c       if(nantest) then
c         PRINT*, "#2 TEMPCR NAN in GETRAD:"
c         PRINT*, "       J = ", J
c         PRINT*, "       K = ", K
c         PRINT*, "   DAYLX = ", DAYLX
c         PRINT*, "  COOLJR = ", COOLJR(J, MZRAD+1-K)
c         STOP
c         ENDIF


      TEMPHR(K)=HEATAV(J,MZRAD+1-K)*daylx

c       nantest=test_nan( temphr(k) )
c       if(nantest) then
c         PRINT*, "#2 TEMPHR NAN in GETRAD:"
c         PRINT*, "       J = ", J
c         PRINT*, "       K = ", K
c         PRINT*, "   DAYLX = ", DAYLX
c         PRINT*, "  HEATAV = ", HEATAV(J, MZRAD+1-K)
c         STOP
c         ENDIF


cjer  do psc's also
      
      tempsc(k)=psctotal(j,mzrad+1-k)
c      if(j.eq.3) write(6,*) 'in getrad: j,k,tempsc = ',j,k,tempsc(k)      

30	continue
 
      CALL INTER(D3,ZP,M$,TEMPHR,ZJRmid,MZRAD)
      CALL INTER(D4,ZP,M$,TEMPCR,ZJRmid,MZRAD)
cjer

c interpolate near-ir psc optical depth (written out just to look at)

      call inter(d5,zp,m$,tempsc,zjrmid,mzrad)
    
C	USE NEWTONIAN COOLING ABOVE zp0x KM
 
      ZP0X=85.
      DO 35 K=1,M$
      IF (ZP(K).LE.ZP0X) WGT=1.0
      IF (ZP(K).GT.ZP0X) WGT=EXP(-(ZP(K)-ZP0X)/5.)
      COOL(J,K)=D4(K)*WGT+(1.-WGT)*VNC(K)*(T(J,K)-THG(K)/TTOTH(K))

c      IF (COOL(J,K) .NE. COOL(J,K)) THEN
c       nantest=test_nan( cool(j,k) )
c       if(nantest) then
c         PRINT*, "#1 COOL NAN in GETRAD:"
c         PRINT*, "       J = ", J
c         PRINT*, "       K = ", K
c         PRINT*, "     WGT = ", WGT
c         PRINT*, "     VNC = ", VNC
c         PRINT*, "      D4 = ", D4(K)
c         PRINT*, "       T = ", T(J, K)
c         PRINT*, "     THG = ", THG(K)
c         PRINT*, "   TTOTH = ", TTOTH(K)
c         STOP
c         ENDIF


      HEAT(J,K)=D3(K)

c      IF (HEAT(J,K) .NE. HEAT(J,K)) THEN
c       nantest=test_nan( heat(j,k) )
c       if(nantest) then
c         PRINT*, "#1 HEAT NAN in GETRAD:"
c         PRINT*, "       J = ", J
c         PRINT*, "       K = ", K
c         PRINT*, "      D3 = ", D3(K)
c         STOP
c         ENDIF

      if (d5(k) .lt. 1.e-6) d5(k)=0.0
c      if(j.eq.3) write(6,*) 'in getrad: j,k,zp,d5 = ',j,k,zp(k),d5(k)
c now psc's
      pscdyn(j,k)=d5(k)

  35  continue
 
  50  CONTINUE

c    two types of tropospheric heat sources are implemented
c    below is the cunnold source 980
c    starting at 990 is the dopplick source

c      goto 990    ! SKIP TO DOPPLICK


      IF(ISW(28).EQ.1) goto 991  !SKIP TROPO. KLUGES ALTOGETHER - YES



C
C     ADD TROPOSPHERIC HEAT SOURCE FROM CUNNOLD ET AL.
980      ZP00=1.3*H

      DX=real(IDAYX)/(IDAYS(MONTH+1)-IDAYS(MONTH))
      DO 111 J=1,N$
      DO 111 K=1,M$
  111 DUM1(J,K)=TFROP(J,K,MONTH) + 
     1           (TFROP(J,K,MONTH+1)-TFROP(J,K,MONTH))*DX

c	PRINT 777,MONTH,DX

      CALL GETGLM(DUM1,GHET,C,N$,M$)







c               ------ DUM1 NOW CONTAINS THE INTERPOLATED
C               ------ (FROM MONTHLY MEANS) NMC CAC TEMPERATURE
C               ------ FOR THE CURRENT MODEL DAY

  
      ENHFAC=ISW(39)/100.
      ENHLAT=ISW(40)*1.00
      ENHHGT=ISW(41)*1.00

      DO K=1,M$
      DO J=1,N$

         EQENH=(1.-ABS(YP(J))/ENHLAT)*ENHFAC
         EQENH=AMAX1(EQENH,0.000 )

         EQENH=(1. -   ZP(K) /ENHHGT)*EQENH
         EQENH=AMAX1(EQENH,0.000 )

         DUM1(J,K)=(1.00+EQENH)*DUM1(J,K)
      END DO
      END DO





      HCFAC=ISW(38)/100.

      F1=1.0*HCFAC
      F2=1.0*HCFAC

      WRITE(6,*) '  USING CUNNOLDIAN TROPO IN JR RAD  (FACTOR)=',
     2              ISW(38)


      DO 210 K=1,M$
      IF (ZZ(K).GT.ZP00+1.0*H) GO TO 210
      VNC2=VNC(K)
      WGT=1.0
      IF (ZZ(K).GT.ZP00) WGT=EXP(-(ZZ(K)-ZP00)/H)
      DO 200 J=1,N$
      HEAT(J,K)=WGT*VNC2*DUM1(J,K)*F2+(1.-WGT)*HEAT(J,K)
      COOL(J,K)=(1.-WGT)*COOL(J,K)+F1*WGT*VNC2*(T(J,K)-THG(K)/TTOTH(K))
200   CONTINUE

210   CONTINUE

c       write(*,*) ' after loop 210 in getrad2 ' 
c       do 2222 k=1,m$
c        write(*,*) 'heat , cool       ',heat(3,k),cool(3,k)
c        write(*,*) ' .***.. '
c2222   continue 
c       write(*,*) ' .................... '
c       write(*,*) ' .................... '
c------------------------------------------------------
c------------------------------------------------------
C	END CUNNOLD HEAT SOURCE FOR TROPOSPHERE
c------------------------------------------------------
c------------------------------------------------------

	GOTO 991
c------------------------------------------------------
c------------------------------------------------------




c------------------------------------------------------
c------------------------------------------------------
C	FIXED HEAT SOURCE FOR TROPOSPHERE
c	plus a type of cunnold heating
c------------------------------------------------------
990      ZP00=1.5*H
c------------------------------------------------------

	DAYM=IDAYX
	DX=1./(IDAYS(MONTH+1)-IDAYS(MONTH))
	DX2=real(IDAYX)/(IDAYS(MONTH+1)-IDAYS(MONTH))
        DO 2111 K=1,M$
        DO 2111 J=1,N$
        DUM3(J,K)=TFROP(J,K,MONTH) + 
     1           (TFROP(J,K,MONTH+1)-TFROP(J,K,MONTH))*DX2
2111   DUM2(J,K)=(TFROP(J,K,MONTH+1)-TFROP(J,K,MONTH))*DX

C	DUM3 IS THE NMC TROPOSPHERIC TEMPERATURES FOR THAT DAY
C	DUM2 IS D(DUM3)/DT IN DEG/DAY

      CALL GETGLM(DUM2,GHET,C,N$,M$)
      CALL GETGLM(DUM3,GHET,C,N$,M$)

      MNN=MONTH-1
      MNP=MONTH+1

      IF (MONTH.EQ.12) MNP=1
      IF (MONTH.EQ.1) MNN=12

	IF (IDAYX.LT.15) THEN 
	DAYM=IDAYX+15
	M1=MNN
	M2=MONTH
	ENDIF

	IF (IDAYX.GE.15) THEN
	DAYM=IDAYX-15
	M1=MONTH
	M2=MNP
	ENDIF

	DX=real(DAYM)/30.

      DO 112 K=1,M$
      DO 112 J=1,N$
112   DUM1(J,K)=0.

      DO 113 J=1,N$
      DO 113 K=1,10
113   DUM1(J,K)=TEHE(J,K,M1) + 
     1           (TEHE(J,K,M2)-TEHE(J,K,M1))*DX

C	DUM1 IS THE JACKMAN HEATING RATES
C	DUM3 IS THE NMC TROPOSPHERIC TEMPERATURES FOR THAT DAY
C	DUM2 IS D(DUM3)/DT IN DEG/DAY

c	PRINT 777,MONTH,DX
777	FORMAT(' FROM RADLIB2, MONTH, DAYX',I5,F10.2)

      CALL GETGLM(DUM1,GHET,C,N$,M$)
	dr=1./86400.
      DO 2100 K=1,M$
      IF (ZZ(K).GT.ZP00+1.5*H) GO TO 2100
      WGT=1.0
      IF (ZZ(K).GT.ZP00) WGT=EXP(-(ZZ(K)-ZP00)/H)
C	DUM1 IS THE JACKMAN HEATING RATES
C	DUM3 IS THE NMC TROPOSPHERIC TEMPERATURES FOR THAT DAY
C	DUM2 IS D(DUM3)/DT IN DEG/DAY
      VNC2=VNC(K)
      DO 2000 J=1,N$
      COOL(J,K)=(1.-WGT)*COOL(J,K)+WGT*VNC2*(T(J,K)-THG(K)/TTOTH(K))
      HEAT(J,K)=WGT*dr*DUM1(J,K)+(1.-WGT)*HEAT(J,K)
     1        + WGT*VNC2*dum3(j,k)
c     2        + WGT*dum2(j,k)*dr
2000   CONTINUE



2100   CONTINUE

c------------------------------------------------------
c------------------------------------------------------
c  end fixed heat source 
c------------------------------------------------------

991   CONTINUE


C     CALL OUTA(1,NSTEP,NDAY,NYEAR,YP,ZP,HEAT,N$,M$,0)
C     CALL OUTA(2,NSTEP,NDAY,NYEAR,YP,ZP,COOL,N$,M$,0)
      
c      CALL GETGLM(HEAT,GHET,C,N$,M$)
c      CALL GETGLM(COOL,GCOL,C,N$,M$)
 
c      CALL OUTA(1,NSTEP,NDAY,NYEAR,YP,ZP,HEAT,N$,M$,0)
c      CALL OUTA(2,NSTEP,NDAY,NYEAR,YP,ZP,COOL,N$,M$,0)
 
 
      DO 2 K=1,M$
      DO 2 J=1,N$
      HEAT(J,K)=HEAT(J,K)*TTOTH(K)
      COOL(J,K)=COOL(J,K)*TTOTH(K)

c      IF (HEAT(J,K) .NE. HEAT(J,K)) THEN
c       nantest=test_nan( heat(j,k) )
c       if(nantest) then
c         PRINT*, "HEAT NAN in GETRAD:"
c         PRINT*, "       J = ", J
c         PRINT*, "       K = ", K
c         PRINT*, "   TTOTH = ", TTOTH(K)
c         STOP
c         ENDIF

c      IF (COOL(J,K) .NE. COOL(J,K)) THEN
c       nantest=test_nan( cool(j,k) )
c       if(nantest) then
c         PRINT*, "COOL NAN in GETRAD:"
c         PRINT*, "       J = ", J
c         PRINT*, "       K = ", K
c         PRINT*, "   TTOTH = ", TTOTH(K)
c         STOP
c         ENDIF

2     CONTINUE
 
      RETURN
      END
