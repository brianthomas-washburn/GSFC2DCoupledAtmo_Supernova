
C     NOTE: THIS WHOLE MODULE NEEDS TO BE REVISED TO BE SELF-CONSISTENT AND
C     DO TRACER INITIALIZATIONS IN A MORE RATIONAL MANNER.  THE TR_TABLE.DAT
C     FILE PROBABLY COULD USE REVISION AS WELL.   PEM, 9/17/96.

 
      SUBROUTINE THE_TR

C     THIS ROUTINE SETS UP THETA AS CONSTITUENT DISTRIBUTION
C     IN XC(*,*,1). ALSO INITIALIZES XC(*,*,2) FOR UB ADV.

!      USE degree_trig

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'

C     DUM1X already declared in COMMOND.INC
C      DIMENSION TOT_TH(N$,M$), DUM1X(NP$,MP$), TGP(MP$), DTGP(MP$)
      DIMENSION TOT_TH(N$,M$), TGP(MP$), DTGP(MP$)
    

C********************************************
C SETUP INITIAL POTENTIAL TEMP. DISTRIBUTION
C
      DO 1 K=1,M$
      DO 1 J=1,N$
C          print*, J, K, IT2
          TOT_TH(J,K)=THG(K)+TH(J,K,IT2)
 1    CONTINUE


      CALL Regrid(DUM1X,YPP,ZPP,NP$,MP$,TOT_TH,
     1  YP,ZP,N$,M$,0)

      DO 2021 K=1,MP$
      DO 2021 J=1,NP$
         XC(J,K,1)=DUM1X(J,K)
         X(J,K,1)=XC(J,K,1)*RHO0(J,K)
         X0(J,K,1)=X(J,K,1)
2021  CONTINUE


C********************************************
C SETUP INITIAL ANGULAR MOMENTUM OR ZONAL 
C WIND DISTRIBUTION
C
      OMEG=TOMEG/2.0         ! 2.*PI/DAYL

      DO 2022 K=1,MP$
      DO 2022 J=1,NP$
c         XC(J,K,2)=0.                    ! UB is tracer
         XC(J,K,2)=A*OMEG*CST(J)*CST(J)  !Angular Mom. is tracer
         X(J,K,2)=XC(J,K,2)*RHO0(J,K)
         X0(J,K,2)=X(J,K,2)
2022  CONTINUE
C********************************************
C WRITE OUT INITIAL DISTRIBUTIONS
C
c      CALL OUTA(181,NSTEP,NDAY,NYEAR,YPP,ZPP,XC(1,1,1),NP$,MP$,0)
c      CALL OUTA(182,NSTEP,NDAY,NYEAR,YPP,ZPP,XC(1,1,2),NP$,MP$,0)
      WRITE(*,*)' .... SET UP THETA AS TRACER  '


      RETURN
      END


      SUBROUTINE TRCINIT(NDYN)

C     THIS ROUTINE SETS UP INITIAL TRACER DISTRIBUTIONS

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'
      DIMENSION OZON(M$),OZONP(MP$),XN2OIN(18,30),XLATS(18),
     1PXIN(30),ZXIN(30),XN2OPH(18,30)

      if(ndyn .ge. ncon) go to 24 
      DO 23 N=NDYN+1,NCON

      DO 2 K=1,MP$
      DO 2 J=1,NP$

      X0(J,K,N)=0.0
      XC(J,K,N)=0.
      X(J,K,N)=0.
      XX(J,K,N)=0.0
      XY(J,K,N)=0.0
      XXX(J,K,N)=0.0
      XYY(J,K,N)=0.0
      XXY(J,K,N)=0.0
2     CONTINUE

C     SET UP A TRACER DISTRIBUTION
c     x() is mass distribution
c     xc() is mixing ratio
c     x0() is initial distribution (mass)

      if (n.eq.(NDYN+1)) then
C********** N2O DISTRIBUTION
      READ(85,10) NLAT
10    FORMAT(I5)
      READ(85,11) (XLATS(J),J=1,NLAT)
11    FORMAT(6F12.2)
      READ(85,10) NPRESS
      READ(85,11) (PXIN(J),J=1,NPRESS)

C     READ IN N2O

      DO 333 I=1,NLAT
      READ(85,10) ILAT
      READ(85,12) (XN2OIN(ILAT,J),J=1,NPRESS)
12    FORMAT(6E12.2)
333   CONTINUE

C     READ IN PHOTOALYSIS

      DO 334 I=1,NLAT
      READ(85,10) ILAT
      READ(85,12) (XN2OPH(ILAT,J),J=1,NPRESS)
334   CONTINUE
      CLOSE(85)

C     NOW LOGRITHMICALLY INTERPOLATE DATA
      PRINT 556
556   FORMAT(' N2O READ IN COMPLETED')

      DO 433 K=1,NPRESS
      ZXIN(K)=1.0E-5*H*ALOG(1000./PXIN(K))
      DO 433 J=1,NLAT
      XN2OIN(J,K)=ALOG10(XN2OIN(J,K))
      XN2OPH(J,K)=ALOG10(XN2OPH(J,K))
433   CONTINUE

      CALL Regrid(DUM1X,YPP,ZPP,NP$,MP$,XN2OIN,
     1  XLATS,ZXIN,NLAT,NPRESS,0)
      CALL Regrid(DUM2X,YPP,ZPP,NP$,MP$,XN2OPH,
     1  XLATS,ZXIN,NLAT,NPRESS,0)

      DO 2021 K=1,MP$
      DO 2021 J=1,NP$
      XN2OPHX(J,K)=10.**DUM2X(J,K)
      X0(J,K,N)=10.**DUM1X(J,K)
      XC(J,K,N)=X0(J,K,N)
      X(J,K,N)=X0(J,K,N)*RHO0(J,K)
      X0(J,K,N)=X(J,K,N)

2021   CONTINUE
 
c      CALL OUTA(24,NSTEP,NDAY,NYEAR,YPP,ZPP,XN2OPHX,Np$,Mp$,0)
      endif



C     **** EQUATORIAL STRATO. BLOB ****
      IF (N.EQ.(NDYN+10)) THEN
      DO K=1,MP$
      DO J=1,NP$
         X0(J,K,N)=0.0
         IF( ABS(YPP(J)).LE.30.0) THEN 
         IF((ZPP(K).LE.20.).AND.(ZPP(K).GE.10.))THEN
           X0(J,K,N)=1.0
         ENDIF
         ENDIF
         XC(J,K,N)=X0(J,K,N)
         X(J,K,N)=X0(J,K,N)*RHO0(J,K)
         X0(J,K,N)=X(J,K,N)
      END DO
      END DO 
      ENDIF

C     **** LITTLE VERTICAL COLUMNS ****
      IF (N.EQ.(NDYN+200)) THEN
      DO 201 K=1,MP$
      DO 201 J=1,NP$,4
       X0(J,K,N)=1.0
C      X0(J,K,N)=OZONP(K)
      XC(J,K,N)=X0(J,K,N)
      X(J,K,N)=X0(J,K,N)*RHO0(J,K)
      X0(J,K,N)=X(J,K,N)
201   CONTINUE
      ENDIF

C     ****  ****
      IF (N.EQ.(NDYN+300)) THEN
      DO K=1,MP$
      DO J=1,NP$
       X0(J,K,N)=1.0
       IF(ZPP(K).GT.50.) X0(J,K,N)=2.0
C      X0(J,K,N)=OZONP(K)
      XC(J,K,N)=X0(J,K,N)
      X(J,K,N)=X0(J,K,N)*RHO0(J,K)
      X0(J,K,N)=X(J,K,N)
      END DO
      END DO
      END IF

C     ****  ****
      IF (N.EQ.(NDYN+2)) THEN
      DO K=1,MP$
      DO J=1,NP$
       X0(J,K,N)=0.0
       IF((ZPP(K).GE.50.).AND.(ZPP(K).LE.60.)) X0(J,K,N)=2.0
       XC(J,K,N)=X0(J,K,N)
       X(J,K,N)=X0(J,K,N)*RHO0(J,K)
       X0(J,K,N)=X(J,K,N)
      END DO
      END DO
      END IF

C     ****  ****  ****  ****  ****  ****  **** *** ***
C**     TRACER INITIALLY ALIGNED WITH HEIGHT FIELD
C
      IF (N.EQ.(NDYN+3)) THEN
      DO K=1,MP$
      DO J=1,NP$
       X0(J,K,N)=0.0
c       IF((ZPP(K).GE.60.).AND.(ZPP(K).LE.70.)) X0(J,K,N)=2.0
c       IF((ZPP(K).GE.70.).AND.(ZPP(K).LE.80.)) X0(J,K,N)=4.0
     
cc       X0(J,K,N)=INT(ZPP(K)/10.)*10.
       X0(J,K,N)= ZPP(K)

       XC(J,K,N)=X0(J,K,N)
       X(J,K,N)=X0(J,K,N)*RHO0(J,K)
       X0(J,K,N)=X(J,K,N)
      END DO
      END DO
      END IF
c************************************************

C     ****  ****  ****  ****  ****  ****  **** *** ***
C**     ``CHECKER BOARD'' TRACER
C
      IF (N.EQ.(NDYN+4)) THEN
      DO K=1,MP$
      DO J=1,NP$
       X0(J,K,N)=0.0
c       IF((ZPP(K).GE.60.).AND.(ZPP(K).LE.70.)) X0(J,K,N)=2.0
c       IF((ZPP(K).GE.70.).AND.(ZPP(K).LE.80.)) X0(J,K,N)=4.0
     
       X0(J,K,N)= ZPP(MP$)-(ZPP(K)/10.)*10.+ 
     > 1.*( 90.*COS( YPP(J)*PI/180. ) ) 

       XC(J,K,N)=X0(J,K,N)
       X(J,K,N)=X0(J,K,N)*RHO0(J,K)
       X0(J,K,N)=X(J,K,N)
      END DO
      END DO
      END IF
c************************************************

      IF (N.EQ.(NDYN+30)) THEN
C     **** A LINEAR DISTRIBUTION ****
      DO 202 K=1,MP$
      DO 202 J=1,NP$
      X0(J,K,N)=MP$-K
      XC(J,K,N)=X0(J,K,N)
      X(J,K,N)=X0(J,K,N)*RHO0(J,K)
      X0(J,K,N)=X(J,K,N)
202   CONTINUE
      ENDIF

      IF (N.EQ.(NDYN+400)) THEN
C     **** constant field ****
         K=1
         DO 203 K=1,MP$
         DO 203 J=1,NP$
         X0(J,K,N)=1.
         XC(J,K,N)=X0(J,K,N)
         X(J,K,N)=X0(J,K,N)*RHO0(J,K)
         X0(J,K,N)=X(J,K,N)
203      CONTINUE
      ENDIF


23    CONTINUE
24    continue

      CALL TABINIT




      RETURN
      END

      SUBROUTINE TABINIT

C     THIS ROUTINE SETS UP INITIAL TRACER DISTRIBUTIONS
c     dbc 2/18/94: y dimension of trtab limits number of tracers to le 20


C     Preliminary Revisions 9/17/96 by PEM to use new timing (timer.f).

      include "timer.inc"


      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'
      DIMENSION OZON(M$),OZONP(MP$),XN2OIN(18,30),XLATS(18),
     1PXIN(30),ZXIN(30),XN2OPH(18,30),trtab(10,20)
      CHARACTER*64 HEAD$

      LOGICAL TABFLAG
      DATA    TABFLAG /.FALSE./

 100  FORMAT ( I3,I3,9F8.2 )
 101  FORMAT (A64)


C     Modified 9/16/96 by PEM for new timing (timer.f).  Table should be
C     read the first time SUBROUTINE TABINIT is called, but not thereafter.
C     In previous incarnations, NSTEP counted throughout the model run and
C     was used here to determine whether or not this was the first call to
C     SUBROUTINE TABINIT.  NSTEP now counts through the dynamic timesteps in
C     a single model day and LOGICAL TABFLAG added to keep track of whether
C     the routine has been called before.
   
C      IF (NSTEP.EQ.0) THEN

      IF (.NOT.TABFLAG) THEN

      READ(43,101) HEAD$      

      DO N=1,ISW(7) 
         READ(43,100) NYDO, NTR, YC0, ZC0, RY, RZ, RYY, RZZ, RYZ
     1               ,RVALC, RRADIUS


         TRTAB( 1,N)=NYDO*1.000
         TRTAB( 2,N)=YC0
         TRTAB( 3,N)=ZC0
         TRTAB( 4,N)=RY
         TRTAB( 5,N)=RZ
         TRTAB( 6,N)=RYY
         TRTAB( 7,N)=RZZ
         TRTAB( 8,N)=RYZ
         TRTAB( 9,N)=RVALC
         TRTAB(10,N)=RRADIUS

         TABFLAG = .TRUE.

      END DO
      CLOSE(UNIT=43)
      END IF


      DO N=1,ISW(7) 

         RYX =0.00
         RZX =0.00
         RYYX=0.00
         RZZX=0.00
         RYZX=0.00

         NYDO=INT(TRTAB(1,N))

C        NYDO specifies what day of the model run a time-dependent tracer
C        should be initialized on.  It is the first number in the
C        tracer descriptor in TR_table.dat.  In the past, the day number
C        was translated here into a NSTEP number (NSTINI), assuming that NSTEP
C        counted timesteps through the model run.  Now, the day number form
C        is retained and the comparison is mafe with IDOR360 to determine
C        if a given tracer should be initialized.

C         NSTINI=NYDO*ISW(6)     ! TIME-STEP ON WHICH TO INIT. TRACER
C         IF(NYDO.EQ.0) NSTINI=1

         IF (NYDO .EQ. 0) NYDO = 1

C         IF (NSTEP.EQ.NSTINI) THEN 

         IF (IDOR360.EQ.NYDO) THEN

         YC0  =TRTAB(2,N)
         ZC0  =TRTAB(3,N)
         RY   =TRTAB(4,N)
         RZ   =TRTAB(5,N)
         RYY  =TRTAB(6,N)
         RZZ  =TRTAB(7,N)
         RYZ  =TRTAB(8,N)
         RVALC=TRTAB(9,N)

            DO K=1,MP$
            DO J=1,NP$

               YCOOR= YPP(J)-YC0
               ZCOOR= ZPP(K)-ZC0
 
               IF ( ABS(RY) .LE. 999.0 ) RYX =1.00/(RY)
               IF ( ABS(RZ) .LE. 999.0 ) RZX =1.00/(RZ)
               IF ( ABS(RYY).LE. 999.0 ) RYYX=1.00/(ABS(RYY)*RYY)
               IF ( ABS(RZZ).LE. 999.0 ) RZZX=1.00/(ABS(RZZ)*RZZ)
               IF ( ABS(RYZ).LE. 999.0 ) RYZX=1.00/(ABS(RYZ)*RYZ)

               X0(J,K,N)= RVALC + RZX*ZCOOR + RYX*YCOOR
     >                  + RYYX*YCOOR*YCOOR + RZZX*ZCOOR*ZCOOR
     >                  + RYZX*YCOOR*ZCOOR

               IF ( X0(J,K,N) .LT. 0.000 ) X0(J,K,N)=0.000          

               XC(J,K,N)=X0(J,K,N)
               X(J,K,N) =X0(J,K,N)*RHO0(J,K)
               X0(J,K,N)=X(J,K,N)

               XX(J,K,N)=0.00   ! SET HIGHER ORDER 
               XY(J,K,N)=0.00   ! MOMENTS TO ZERO 
              XXX(J,K,N)=0.00   ! AS INITIAL COND.
              XYY(J,K,N)=0.00
              XXY(J,K,N)=0.00


            END DO
            END DO
    

         NT=180+N
c         CALL OUTA(NT,-1 ,NDAY,NYEAR,YPP,ZPP ,XC(1,1,N),NP$,MP$,1)
         END IF
            
       END DO

       RETURN
       END

       SUBROUTINE GETMIX(ll,NMAX,ASOR)
 
C     COMPUTE THE DIFFUSION BY KYY AND KZZ FOR THE CONSTITUENT MIXING FIELDS
C     SIMPLE CHEMISTRY AS WELL
C     ALSO COMPUTES DRAG,HEATING AND COOLING FOR TH AND UB
cjer fixed up interpolation of kzz (per J. Bacmeister)

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'
      
ccelf      include 'test_value.inc'
ccelf      logical nantest


C  NCEP/GEOS4 TEMP, UBAR at lowest extended grid point (1.015 km) - UT1KM(NP$,2) - 1=TEMP; 2=UBAR (m/sec)
C        and NCEP UBAR at ZPP(2-3) (3-5 km) - U3KM(NP$,2)
C
       COMMON/CTNCEP/TNCEP(N$,M$), TNCEPSFC(NP$,2), TEMPG42D(NP$,3), 
     >               UT1KM(NP$,2), U3KM(NP$,2)


C   NCEP TEMP, BV freq, on extended coupled model grid;  dtedz is dtheta/dz, XN2E is in 1/sec/day

       COMMON/CTNCEPE/TNCEPE(NP$,MP$), XN2E(NP$,MP$), dtedz(NP$,MP$)


C  common for model compute eddy heating/EP flux correction: 
C      EHFC is in K/sec (theta heating); EPZZ is in cm/sec^2

       COMMON/CEHFC/ ehfc(NP$,MP$), epzz(NP$,MP$)


C  COMMON for WACCM convective mass flux for Kzz, interpolated to current day/dynamics grid
C
ccc46       COMMON/CWCONKZD/cmfd(NP$,M$)



C  COMMON for WACCM combined eddy heating rate (theta heating) in K/sec
C
cccwconv48       COMMON/CQWTD/qwttd(NP$,MP$)



C  COMMON for WACCM heating corrections for current day - QWTD, QVTD(NP$,MP$) are K/sec (theta heating)
C          NOT USED
cccc      COMMON/CWCOND/epwd(NP$,MP$), qwtd(NP$,MP$), qvtd(NP$,MP$)



      REAL KYYMAX,KYYF(NCON),KZZF(NCON)
      REAL        IYYF(NCON),IZZF(NCON)

      REAL XKYY(N$,MP$),XKZZ(NP$,M$),XKYZ(N$,MP$),
     1 XKZY(NP$,M$),XZ(NP$,MP$),ZKYY(N$,MP$),WKYY(N$,MP$),
     2 ALPHA(NP$,MP$),XTMP(NP$,MP$),WKZZ(NP$,M$),ZKZZ(NP$,M$)

      REAL hconv(NP$,MP$), ykzmax(NP$,M$)


      CHARACTER*10 CKP$,CKI$

        CKP$='PLANETARY-'
        CKI$='INERTIAL--'
 
  
        IF(NMAX.EQ.1) then
           CALL GET_KSS  ! get KZZ,KYY from misc. parameterizations
        END IF

        IF (NMAX.EQ.0) RETURN
	IF((ISW(4).EQ.0).AND.(ISW(8).EQ.0).AND.(ISW(9).EQ.0)) RETURN


      IB=0

C     FIRST INTERPOLATE KYY, KZZ, KYZ

cjer code added to fix kzz interpolation

      do k=1,m$
      do j=1,np$
      
           XKZZ(J,K) = 0.5*( KZZ(J,K+1) + KZZ(J,K) ) 

           WKZZ(J,K) = 0.5*( KZZ(J,K+1) + KZZ(J,K) )

           ZKZZ(J,K) = 0.5*( ZIKZZ(J,K+1) + ZIKZZ(J,K) ) 
           
      end do
      end do


      do 10 k=1,m$
      
        xkzy(1,k) = kyz(1,k)
        xkzy(np$,k) = kyz(n$,k)
        
        do 10 j=2,n$
        
          xkzy(j,k) = 0.5*( kyz(j-1,k) + kyz(j,k) )
          
   10 continue


	DO 20 J=1,N$

	   WKYY(J,1)  =KYY(J,1)*0.5
	   WKYY(J,MP$)=KYY(J,M$)*0.5

	   ZKYY(J,1)  =ZIKYY(J,1)*0.5
	   ZKYY(J,MP$)=ZIKYY(J,M$)*0.5

	   XKYZ(J,1)  =KYZ(J,1)*0.5
	   XKYZ(J,MP$)=KYZ(J,M$)*0.5

        DO 20 K=2,M$

           WKYY(J,K)=(KYY(J,K-1)+KYY(J,K))*0.5

           ZKYY(J,K)=(ZIKYY(J,K-1)+ZIKYY(J,K))*0.5

           XKYZ(J,K)=(KYZ(J,K-1)+KYZ(J,K))*0.5

20      CONTINUE
 

	DY22=1./(DY*DY)
	DZ22=1./(DZ*DZ)
        XPZ=EXP(-DZ*0.5/H)
        XMZ=EXP(DZ*0.5/H)



c        WRITE(6,102)
c        WRITE( 6, *   ) ' --------- ------ ----- ------ --- '
c        DO J=NMAX,NMAX
c        WRITE( 6, 101 ) j , KYYPR(J), IYYPR(J), KZZPR(J), IZZPR(J)
c     >                    , DMOLPR(J)
c        WRITE( 6, *   ) ' --------- ------ ----- ------ --- '
c        END DO
101   FORMAT(I2,2X, F10.8,5X, F10.8,5X, F10.8,5X, F10.8,5X,F10.8)
102   FORMAT('T#   PWAVE_KYY     INERT._KYY     BACKG._KZZ',
     >   '     INERT._KZZ     MOL._DIFF. ' )


C  define Kzz max in troposphere - ykzmax(NP$,M$)

      do k=1,m$
      do j=1,np$
      ykzmax(j,k)=(.5+3.*(COSD(ypp(j))**4))*(1.-(zp(k)-zp(1))/15.)*1.e4
ccccccc      ykzmax(j,k)=(1.+11.*(COSD(ypp(j))**4))*(1.-(zp(k)-zp(1))/15.)*1.e4
      end do
      end do



C  LOOP USED TO BE OVER NUMBER OF TRANSPORTED SPECIES
C  NOW IS DONE ONCE PER CALL TO GETMIX.  TRACER INDEX
C  IS DRIVEN IN TWOD.F  (JTB 2/3/94)

      DO 5000 N=NMAX,NMAX         ! 1,NMAX 


C     XC IS THE MIXING RATIO, X IS THE MASS

C     STORE OFF THE NON-UPDATED MIXING RATIO  FOR USE BELOW

        DO 505 K=1,MP$
        DO 505 J=1,NP$
           XTMP(J,K) =XC(J,K,N)
           XC(J,K,N) =X(J,K,N)/RHO0(J,K)

c       nantest=test_nan( xc(j,k,n) )
c       if(nantest) then
c              PRINT*, "NAN in GETMIX:"
c              PRINT*, "      J = ", J
c              PRINT*, "      K = ", K
c              PRINT*, "      N = ", N
c              PRINT*, "      X = ", X(J, K, N)
c              PRINT*, "   RHO0 = ", RHO0(J, K)
c              STOP
c              ENDIF

 505    CONTINUE

        

C     GET ANY SOURCES AND SINKS AND LIFETIME (1/ALPHA)

       CALL GETSOR2(N,ALPHA,SOR(1,1,N),IA,IS,IB)  !FOR THETA TR.
cjer  SOR is smoothed net heating       


C   ADD INERTIAL INST. KYY AND P-WAVE KYY 
C   - turn on SMALL Kyy for theta ONLY - below 15 km - horiz eddy heating term in eq 3.5.2e (EF)
C                          blend in between 15km - 10 km
        KYYPR(1) = 1.

        DO K=1,MP$
           DO J=1,N$

              zwt = (16. - zpp(k))/15.
              if (zwt .le. 0.) zwt = 0.
              if (zwt .ge. .3) zwt = .3

              XKYY(J,K) = KYYPR(N)*WKYY(J,K)*zwt + IYYPR(N)*ZKYY(J,K)
              END DO
           END DO


C   Kzz update:
C
C   inertial instability Kzz (IZZPR, ZKZZ) is 0 for theta, 1 for ang mom - keep these as is
C   adjust WKZZ(NP$,M$) (background Kzz) here vs. latitude below 25 km 
C   do FOR THETA ONLY, ie w'th' is different than w'u'
C   also set a minimum of .1 m2/sec below 16.5 km, and define latitudinal max in troposphere
C                                                      ykzmax(NP$,M$)
        IF (N .eq. 1) THEN
           do 577 k=1,m$
           do 577 j=1,np$
              if (zp(k) .lt. 40.) then
                  zwt = (zp(k) - 13.)/2.
                  if (zwt .le. 1.) zwt = 1.
                  if (zwt .ge. 4.) zwt = 4.

                  ywt = (.25 + .5*(COSD(ypp(j))**2))*zwt
cccccc                  ywt = (.25 + .45*(COSD(ypp(j))**2))*zwt
                  if (ywt .le. .25) ywt = .25
                  if (ywt .ge. 1.) ywt = 1.

                  wkzz(j,k) = wkzz(j,k)*ywt

C  BELOW 30 km, define wkzz by dtedz(NP$,MP$) (dtheta/dz), or BV freq, then set low and high limits
C                                                            ! XN2E(NP$,MP$)

cccc       xn2kz=EXP(ALOG(.3/3.)/ALOG(30./4.)*ALOG(dtedz(j,k)/4.) +ALOG(3.))

       xn2kz = EXP(ALOG(.3/5.)/ALOG(40./7.)*ALOG(XN2E(j,k)/7.)+ALOG(5.))
             wkzz(j,k) = xn2kz*1.e4

c  low latitude lower trop:
              if (ABS(ypp(j)) .le. 40.  .and.  zp(k) .le. 10.) then
                  zkx = 5.e4*(COSD(ypp(j))**2)
                  if (wkzz(j,k) .le. zkx) wkzz(j,k) = zkx
              endif

c  high latitude trop:
              if (ABS(ypp(j)) .gt. 40.  .and.  zp(k) .le. 10.) then
                if (wkzz(j,k) .ge. .7e4) wkzz(j,k)= .7e4
              endif

c
ccc              if (ABS(ypp(j)) .gt. 40.  .and.  zp(k) .le. 15.) then
ccc                if (wkzz(j,k) .ge. ykzmax(j,k)) wkzz(j,k)= ykzmax(j,k)
ccc              endif

C  set MAX above 15 km everywhere:

ccc              if (zp(k) .gt. 15.) then
ccc                  zkx = (.2 + .3*(COSD(ypp(j))**2))*1.e4
ccc                  if (wkzz(j,k) .ge. zkx) wkzz(j,k) = zkx
ccc              endif

       if (zp(k) .lt. 16.5 .and. wkzz(j,k) .le. 1000.) wkzz(j,k) = 1000.


CC WCONV238 - enhance THETA Kzz in tropics at 14-23 km

               if (ABS(ypp(j)) .lt. 45.) then
               if (zp(k) .ge. 15.  .and.  zp(k) .le. 23.) then

                  xzfac0 = 2.
                  if (zp(k) .lt. 17.  .or.  zp(k) .ge. 21.) xzfac0 = 1.5

                  xzfac = xzfac0*COSD(ypp(j))
                  if (xzfac .le. 1.) xzfac = 1.

                  wkzz(j,k) = wkzz(j,k)*xzfac
               endif
               endif

           endif

 577       CONTINUE
        ENDIF


        DO K=1,M$
        DO J=1,NP$
C
C       Modified 4/18/95 by PEM, to avoid irrelevant NAN in ZKZZ.
C           XKZZ(J,K)=KZZPR(N)*WKZZ(J,K)+IZZPR(N)*ZKZZ(J,K)
C     >              +DMOLPR(N)*DMOL(K)         
           XKZZ(J,K) = DMOLPR(N)*DMOL(K)
           IF (KZZPR(N).NE.0.) XKZZ(J,K)=XKZZ(J,K)+KZZPR(N)*WKZZ(J,K)
           IF (IZZPR(N).NE.0.) XKZZ(J,K)=XKZZ(J,K)+IZZPR(N)*ZKZZ(J,K)
C       End modified 4/18/95 by PEM.
C
        END DO
        END DO
           
        IF ( N .EQ. 1 )
     >  CALL OUTA(16,NSTEP,NDAY,NYEAR,YPP,ZP,XKZZ,NP$,M$,0)



C      CALCULATE THE MIXING TERM (1/COS(LAT))*(d(COS(LAT)*KYY*(dX/dy))/dy)
C      and 1/rho*(d(rho*kzz*(dx/dz))/dz)

C---COMPUTE THE LARGEST POSSIBLE TIME STEP FOR THE TRACER DIFFUSION
C---CALCULATION.

	TMIN=1.5*DT

	DO 900 K=1,MP$
        DO 900 J=1,N$
	TMINYY=0.5/(DY22*(XKYY(J,K)+1.))
        TMIN=AMIN1(TMIN,TMINYY)
900	CONTINUE

	DO 901 K=1,M$
        DO 901 J=1,NP$
        TMINZZ=0.5/(DZ22*(XKZZ(J,K)+1.))
        TMIN=AMIN1(TMIN,TMINZZ)
C
C       Added 4/14/95 by PEM for debugging.
C
        npem = DT/TMIN+1
        if (npem.gt.100) then
           print*, "Program Is Runaway at:"
           print*, J, K, NYEAR, MONTH, NDAY, npem
           print*, "DZ22 = ", DZ22, "    XKZZ = ", XKZZ(J,K)
           print*, "WKZZ = ", WKZZ(J,K), "    IZZPR = ", IZZPR(N)
           print*, "ZKZZ = ", ZKZZ(J,K), "    KZZPR = ", KZZPR(N)
           print*, "ZIKZZ = ", ZIKZZ(J,K)
cjer add more printout
           print*,'KZZ(J,K),KZZ(J,K+1) = ',KZZ(J,K),KZZ(J,K+1)	   
	   
           print*, "DMOL = ", DMOL(K), "    DMOLPR = ", DMOLPR(N)
           endif
cjer added for debugging
c      if(j.eq.1.and.k.eq.45) then
c           print*, J, K, NYEAR, MONTH, NDAY, npem
c           print*, "DZ22 = ", DZ22, "    XKZZ = ", XKZZ(J,K)
c           print*, "WKZZ = ", WKZZ(J,K), "    IZZPR = ", IZZPR(N)
c           print*, "ZKZZ = ", ZKZZ(J,K), "    KZZPR = ", KZZPR(N)
c           print*, "ZIKZZ = ", ZIKZZ(J,K)
c           print*, "DMOL = ", DMOL(K), "    DMOLPR = ", DMOLPR(N)
c      endif     
C
C       End Added 4/14/95 by PEM.
901	CONTINUE

650	CONTINUE
	DO 902 K=1,MP$
        DO 902 J=1,NP$
	DUM1X(J,K)=0.
	DUM2X(J,K)=0.
	DUM3X(J,K)=0.
	DUM4X(J,K)=0.
902	CONTINUE


	NSTDIF=DT/TMIN+1   ! RATIO+1 OF DT AND LARGEST ALLOWABLE 
C                          ! TIME STEP

C---NEW TIME STEP NEEDED IF ANY OF THE K'S IS TOO LARGE; 
C---TO PREVENT NEGATIVE MIXING RATIOS

	DTN=DT/NSTDIF 

	if (nstdif.gt.1) PRINT 3401,N,NSTDIF,DTN
3401	FORMAT(' T-STEP CUT FOR CONST.=',I2,'; BY ',I3,
     2         '; NEW DT = ',G14.4)


        IF (nstdif.gt.100) then               !Added 4/19/95 by PEM.
           print*, "NSTDIF triggers STOP..."
           stop
           endif


	DO 4000 NST1=1,NSTDIF   ! LOOP NSTDIF TIMES USING 
C                               ! REDUCED TIME STEP DTN TO COME
C                               ! A FULL DT LATER.

C       KYY

        IF (ISW(8).GE.1) THEN

        KYYMAX=-10.000
        DO K=1,MP$
        DO J=1,N$
           KYYMAX=AMAX1(KYYMAX,XKYY(J,K) )
        END DO
        END DO
  
c        IF (NST1.eq.1) write(6,3402) N , KYYMAX
3402    FORMAT(' ... DIFFUSING TRACER# ',i3,' KYY_max ',E20.4)

	DO 500 K=1,MP$
	DO 510 J=2,N$
        DUM1X(J,K)= (C(J)*XKYY(J,K)*(XC(J+1,K,N)-XC(J,K,N))-
     1    C(J-1)*XKYY(J-1,K)*(XC(J,K,N)-XC(J-1,K,N)))*DY22/CST(J)
510     CONTINUE
        DUM1X(1,K)= C(1)*XKYY(1,K)*(XC(2,K,N)-XC(1,K,N))*DY22/CST(1)
        DUM1X(NP$,K)=-C(N$)*XKYY(N$,K)*(XC(NP$,K,N)-XC(N$,K,N))
     1     *DY22/CST(NP$)
500     CONTINUE

        ENDIF


c***********************************************************
C       KZZ
C 
C  for Theta: bottom level is reset below, so Kzz term reduces dth/dz slope, so it's always cooling
C
       IF (ISW(4).GE.1) THEN

       DO 600  J=1,NP$
       DO 610  K=2,M$
       DUM2X(J,K)=(XPZ*XKZZ(J,K)*(XC(J,K+1,N)-XC(J,K,N))-
     1XMZ*XKZZ(J,K-1)*(XC(J,K,N)-XC(J,K-1,N)))*DZ22
610    CONTINUE
       DUM2X(J,1)=(XPZ*XKZZ(J,1)*(XC(J,2,N)-XC(J,1,N)))*DZ22
       DUM2X(J,MP$)=(-XMZ*XKZZ(J,M$)*(XC(J,MP$,N)-XC(J,M$,N)))*DZ22
600    CONTINUE
   
       ENDIF


C      KYZ 

	IF (ISW(9).GE.1) THEN
	DO 720 J=1,NP$
	DO 710 K=2,M$
        XZ(J,K)=(XC(J,K+1,N)-XC(J,K-1,N))/(DZ*2.)
710     CONTINUE
        XZ(J,1)=0.
        XZ(J,MP$)=0.
720	CONTINUE

	DO 700 K=1,MP$
	DO 730 J=2,N$
	DUM3X(J,K)= 0.5*(C(J)*XKYZ(J,K)*(XZ(J+1,K)+XZ(J,K))-
     1C(J-1)*XKYZ(J-1,K)*(XZ(J,K)+XZ(J-1,K)))/(DY*CST(J))
730	CONTINUE
	DUM3X(1,K)=0.5*(C(1)*XKYZ(1,K)*(XZ(2,K)+XZ(1,K)))/(DY*CST(J))
	DUM3X(NP$,K)=-0.5*(C(N$)*XKYZ(N$,K)*(XZ(N$,K)+XZ(N1$,K)))
     1/(DY*CST(J))
C	CHKYY=0.
C	DO J=1,NP$
C        CHKYY=CHKYY+RHO0(J,K)*DUM3X(J,K)
C	ENDDO
C        PRINT 457,K,CHKYY
C457	FORMAT(' K,CHKYZ=',I5,G14.4)
700	CONTINUE

C     KZY

	DO 820 K=1,MP$
	DO 810 J=2,N$
810     XZ(J,K)=(XC(J+1,K,N)-XC(J-1,K,N))/(2.*DZ)
	XZ(1,K)=0.
820	XZ(NP$,K)=0.

	DO 800 J=1,NP$
	DO 830 K=2,M$
830	DUM4X(J,K)=0.5*(XPZ*XKZY(J,K)*(XZ(J,K+1)+XZ(J,K))
     1 -XMZ*XKZY(J,K-1)*(XZ(J,K)-XZ(J,K-1)))/DZ
	DUM4X(J,1)=0.5*(XPZ*XKZY(J,1)*(XZ(J,2)+XZ(J,1)))/DZ
	DUM4X(J,K)=0.5*(-XMZ*XKZY(J,K-1)*(XZ(J,K)-XZ(J,K-1)))/DZ
C	CHKZZ=0.
C	DO K=1,MP$
C        CHKZZ=CHKZZ+RHO0(J,K)*DUM4X(J,K)
C	ENDDO
C        PRINT 458,J,CHKZZ
C458	FORMAT(' J,CHKZY=',I5,G14.4)
800     CONTINUE
	ENDIF

C
C  EF - add in 20% WACCM eddy heating for current day - QWTTD(NP$,MP$) is in K/sec (theta heating) 
C                                  - this is just for 0-25 km
       DO K=1,MP$
       DO J=1,NP$
         if (N .eq. 1) hconv(J,K) = 0.0      ! QWTTD(J,K)
         if (N .eq. 2) hconv(J,K) = 0.0
       END DO
       END DO



C      DO THE ADVECTION STEP HERE

       DO K=1,MP$
       DO J=1,NP$
          XMIX(J,K,N)=DUM1X(J,K)+DUM2X(J,K)  !ADD UP DIFFUSIVE TERMS 
       END DO                                !FOR EACH X
       END DO
          
       DO 50 K=1,MP$
       DO 50 J=1,NP$
           IF (XC(J, K, N) .GE. 1.0E+14) THEN
              PRINT*, "XC blowup in GETMIX(#2):"
              PRINT*, "XC in GETMIX"
              PRINT*, "      J = ", J
              PRINT*, "      K = ", K
              PRINT*, "      N = ", N
              PRINT*, "   ASOR = ", ASOR
              PRINT*, "    SOR = ", SOR(J, K, N)
              PRINT*, "  DUM1X = ", DUM1X(J, K)
              PRINT*, "  DUM2X = ", DUM2X(J, K)
              PRINT*, "     XC = ", XC(J, K, N)
              STOP
              ENDIF
cjer asor is coming in as 1
c       if(n.eq.1 .and. j.eq.19 .and. k.eq.9) then
c         write(6,*) ' ll = ',ll
c         write(6,*) 'xc before = ',xc(j,k,1)
c	 write(6,*) 'sor,d1x,d2x = ',sor(j,k,1),dum1x(j,k),dum2x(j,k)
c       endif	 

cjer **********here is where xc is updated *******************

       XC(J,K,N)= DTN*( ASOR*SOR(J,K,N)
     > + DUM1X(J,K) + DUM2X(J,K) + hconv(J,K) )
     > + XC(J,K,N)



c       if(n.eq.1 .and. j.eq.19 .and. k.eq.9) then
c         write(6,*) 'xc after = ',xc(j,k,1)
c       endif	 

C
C  NEW (EF, May 2009)
C
C  overwrite lower boundary (MP$=1) to NCEP/GESO4 here, since heat(N$,M$)/cool(N$,M$) arrays 
C  have been interpolated from original (NP$,MP$/R$) grid, so lowest level M$=1 is 2km (not so good)
C    so resetting to NCEP here is like assuming net heating correction
C    Don't worry about advective change or inconsistency w/ net heating used in STRMFUNCTION, since
C    that grid uses the interpolated heat(N$,M$)/cool(N$,M$) arrays
C    w/ MP$ level 1 (1 km) = level 1 of M$ (2km)
C      sfc boundary was also reset after advection in TWODS     
C      use UT1KM(NP$,2),  1=TEMP;  2=UBAR (m/sec) for current day;  TTOTHX(MP$)
C
C     also reset sfc ANG MOM from NCEP UBAR for current day - XC00(NP$,MP$,NCON=2); CST(NP$)
C     RESET UBAR for WCONV75, (turned off for PW40), and set at zpp(2) - U3KM(NP$,2) for WCONV82
C        turned off again for WCONV85
C
       if (K .EQ. 1) then
          if (N .EQ. 1) XC(J,1,1) = ut1km(J,1)*TTOTHX(1)
ccwconv85          if (N .EQ. 2) XC(J,1,2) = 100.*ut1km(J,2)*CST(J) + xc00(J,1,2)
       endif

ccwconv85       if (K .EQ. 2) then
ccwconv85          if (N .EQ. 2) XC(J,2,2) = 100.*u3km(J,1)*CST(J) + xc00(J,2,2)
ccwconv85       endif

ccwconv82       if (K .EQ. 3) then
ccwconv82          if (N .EQ. 2) XC(J,3,2) = 100.*u3km(J,2)*CST(J) + xc00(J,3,2)
ccwconv82       endif



cc      if(n.eq.2 .and. j.eq.21 .and. k.eq.1) then
cc        uuuu = (XC(J,1,2) - xc00(J,1,2))/cst(j)/100.
cc       write(6,*) 'xc after = ', xc(j,1,2), uuuu, ut1km(J,2)
cc      endif	 

           IF (XC(J, K, N) .GE. 1.0E+14) THEN
              PRINT*, "XC blowup in GETMIX(#3):"
              PRINT*, "XC in GETMIX"
              PRINT*, "      J = ", J
              PRINT*, "      K = ", K
              PRINT*, "      N = ", N
              PRINT*, "   ASOR = ", ASOR
              PRINT*, "    SOR = ", SOR(J, K, N)
              PRINT*, "  DUM1X = ", DUM1X(J, K)
              PRINT*, "  DUM2X = ", DUM2X(J, K)
              PRINT*, "     XC = ", XC(J, K, N)
              STOP
              ENDIF

cc     >  /NSTDIF + ( XTMP(J,K) )
cc     >  *(1.-1/NSTDIF))/(1.0+ALPHA(J,K)*DTN)

50    CONTINUE

4000  CONTINUE



C  XKZZ(NP$,M$), WKZZ(NP$,M$), ZKZZ(NP$,M$), CMFD(NP$,M$), TNCEPE(NP$,MP$), XN2E(NP$,MP$)
                                                          ! dtedz(NP$,MP$) is dtheta/dz
cc          if (LL .eq. 12) then
cc             write (445) N$, M$, np$, mp$, N
cc             write (445) yp, ypp, zp, zpp
cc             write (445) xkzz, wkzz, zkzz, ykzmax     ! , cmfd, rtkzz
cc             write (445) xc, ttothx
cc             write (445) tncepe, xn2e, dtedz
cc          endif




C    COMPUTE THE TOTAL MASS AND CONVERT MIXING RATIO TO MASS

	TMASS=0.
	DO 90 K=1,MP$
	DO 90 J=1,NP$
          X(J,K,N)=XC(J,K,N)*RHO0(J,K)
	  TMASS=TMASS+X(J,K,N)
90	CONTINUE

c	PRINT 709,N,TMASS
709	FORMAT('   --->*>*>*> TOTAL MASS FOR CONSTITUENT ',I2,' = ',G14.6)

      NR=18*10+N
c      CALL OUTA(NR,NSTEP,NDAY,NYEAR,YPP,ZPP,XC(1,1,N),NP$,MP$,0)

5000  CONTINUE   !END LOOP OVER SPECIES
c	PRINT 7002
7002    FORMAT('                                  ----END GETMIX --' )

       RETURN
       END




      SUBROUTINE GETSOR2(N,ALPHA,STEMP,IA,IS,IB)

C     GET THE SOURCE MINUS SINK TERM FOR CONSTITUENT N
C     STEMP IS THE SOURCE (MIXING RATIO/SEC)
C     alpha is 1/lifetime of species (1/SEC)
C     IS = 0 IF STEMP=0, IA=0 IF ALPHA=0

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMONT.INC'
      INCLUDE 'COMMOND.INC'
c                                                     ! common of NCEP UBAR/DIFFS in m/sec
      COMMON/CUNCEP/UNCEP(NP$,MP$), UDIFFDAY(NP$,MP$)

C
C  Common of NCEP d(v'T')/dy, d(u'v')/dy from baroclinic waves ONLY
C     interpolated to coupled model grid for current day in TWODS
C     VTC is in K/sec;  UVC is in m/sec^2
C
      COMMON/CVTUV/ VTC(NP$,MP$), UVC(NP$,MP$)



      DIMENSION ALPHA(NP$,MP$),STEMP(NP$,MP$)

C**
C  ZERO OUT SOURCES AND SINKS 
      DO K=1,MP$
      DO J=1,NP$
         STEMP(J,K)=0.00
         ALPHA(J,K)=0.00
      END DO
      END DO
      IS=0
      IA=0
C**
C -- -- -- -- -- -- --



      IF (N.EQ.1) THEN

C     HEATING AND COOLING FOR THETA DEG/SEC
      

      DO K=1,M$
      DO J=1,N$
         GHEATX(J,K)=HEAT(J,K)-COOL(J,K)
      END DO
      END DO
c      print *,'in trcrlib: gheatx(9,1)=',gheatx(9,1)
c      print *,'in trcrlib: heat,cool=',heat(9,1),cool(9,1) 

      DO 55 K=2,M$
      DO 55 J=2,N$
         STEMP(J,K)=( GHEATX(J-1,K-1) + GHEATX(J,K-1)
     1           + GHEATX(J,K) + GHEATX(J-1,K) )/4.0
55    CONTINUE

      DO J=2,N$
         STEMP(J,1)=(GHEATX(J-1,1)+GHEATX(J,1))/2.
         STEMP(J,MP$)=(GHEATX(J-1,M$)+GHEATX(J,M$))/2.
      END DO
      DO K=2,M$
         STEMP(1,K)=(GHEATX(1,K-1)+GHEATX(1,K))/2.
         STEMP(NP$,K)=(GHEATX(N$,K-1)+GHEATX(N$,K))/2.
      END DO
      
      
      STEMP(1,1)=GHEATX(1,1)
      STEMP(1,MP$)=GHEATX(1,M$)
      STEMP(NP$,1)=GHEATX(N$,1)
      STEMP(NP$,MP$)=GHEATX(N$,M$)

C
Cpw45  add in -d(v'T')/dy from NCEP baroclinic waves here
C                     VTC(NP$,MP$) is in K/sec, convert to THETA heating, TTOTHX(MP$)
Cpw46      DO 57 K=1,MP$
Cpw46      DO 57 J=1,NP$
Cpw46        STEMP(J,K) = STEMP(J,K)   ! - VTC(J,K)*TTOTHX(K) - turn off for PW46
Cpw4657    CONTINUE

      DO 56 K=1,MP$
      DO 56 J=1,NP$
      ALPHA(J,K)=0.
56    CONTINUE
      IS=1
      IA=0
      RETURN
      ENDIF




      IF (N.EQ.2) THEN

C     MOM. FORCING FOR UB
      

      NSRFD=ISW(29)*( 1-ISW(35) )
      NSRFD=1


4321  FORMAT( ' SRF.DRAG TIME(DAYS) (just U) ',I3 )


Cpw45  add in -d(u'v')/dy from NCEP baroclinic waves here
C   UVC(NP$,MP$) is in m/sec^2, convert to cm/sec^2, multiply by CST(NP$) as for DRAGX in GETMOM

      DO K=1,MP$
      DO J=1,NP$
         STEMP(J,K)=DRAGX(J,K)    ! - UVC(J,K)*100.*CST(J)    - turn off for PW46      ! M
         ALPHA(J,K)=0.
      END DO
      END DO


      IF(NSRFD.GT.0) THEN 

c      WRITE(6,4321) NSRFD  

Ccelf         IF(ISW(68).EQ.0) THEN
Ccelf         DO J=1,NP$
Ccelf            STEMP(J,1)=-( UBX(J,1)*CST(J) )/( NSRFD*DAYL)
Ccelf            STEMP(J,2)=-( UBX(J,2)*CST(J) )/( NSRFD*DAYL*2)
Ccelf         END DO
Ccelf         END IF


Ccelf   change drag at surface here to time scale of 12 hours (changed from 1 day)
Ccelf   this gives a bit more westerly mom at surface (ALSO done in GETMOM)
Ccelf
Ccelf   relax to NCEP UBAR at surface:  UBX(NP$,MP$), UNCEP(NP$,MP$) is in m/sec - NOT USED
Cc      also load into dragcoup(3,J,1) for output (NO COSINE),  dragcoup(7,NP$,MP$) in COMMOND.INC
CC
ccelf                                           currently, ISW(68)=2 and NSRFD=1



ccpw107 - THIS IS DONE IN GETMOM, and DRAGX loaded in STEMP ABOVE - NO NEED TO REPEAT HERE.....


ccpw107         IF(ISW(68).EQ.2) THEN
ccpw107         DO J=1,NP$
ccpw107            STEMP(J,1) = -( UBX(J,1)*CST(J) )/( DAYL/2.)
ccpw107
ccelf           STEMP(J,1)=-( UBX(J,1)*CST(J) )/( NSRFD*DAYL)
ccelf          STEMP(J,1)= -( (UBX(J,1) - UNCEP(J,1)*100.)*CST(J) )/(DAYL/4.)
ccpw107
ccpw107            dragcoup(3,J,1) = STEMP(J,1)/CST(J)
ccpw107         END DO
ccpw107         END IF


         IF(ISW(68).EQ.1) THEN
           CALL PBLPAR(STEMP)
         END IF

C                         drag at poles set by relaxing to THERMAL WIND- turn OFF!! - EF, 5/09
ccpw44         DO K=1,MP$
ccpw44            J=1
ccpw44            UTH=UTHRMX(j,K)
ccpw44            STEMP(J,K)=-( (UBX(J,K)-UTH)*CST(J) )/( NSRFD*DAYL*5)
ccpw44            J=NP$
ccpw44            UTH=UTHRMX(j,K)
ccpw44            STEMP(J,K)=-( (UBX(J,K)-UTH)*CST(J) )/( NSRFD*DAYL*5)
ccpw44         END DO

ccpw45 - for output:
ccpw46
ccpw46      DO K=1,MP$
ccpw46      DO J=1,NP$
ccpw46         dragcoup(6,J,K) = 0.  ! - UVC(J,K)*100.
ccpw46      END DO
ccpw46      END DO


      END IF



      IS=1
      IA=0
      IB=0
      RETURN
      ENDIF

      IF (N.EQ.30) THEN

C     SOURCES AND SINKS

      CALL PVSOR(STEMP)


      IS=1
      IA=0
      RETURN
      ENDIF


      IF (N.EQ.3) THEN

C     N2O PHOTOLYSIS RATES

      DO 45 K=1,MP$
      DO 45 J=1,NP$
      ALPHA(J,K)=XN2OPHX(J,K)
      STEMP(J,K)= -XC(J,K,N) * ALPHA(J,K)
45    CONTINUE
c      CALL OUTA(23,NSTEP,NDAY,NYEAR,YPP,ZPP,STEMP,Np$,Mp$,0)
      IS=0
      IA=1
      RETURN
      ENDIF

      IF (N.EQ.40) THEN
C     OZONE CHEMISTRY
      F1=0.
      IF (N.EQ.5) F1=1.
      DO 1 K=1,MP$
      KS=K
      IF (K.EQ.MP$) KS=M$
      DO 1 J=1,NP$
      JS=J
      SNA=COS(PI*(YPP(J)-DECL)/180.)
      IF(SNA.LE.0.) SNA=0.
      F2=F1*SNA
      IF (J.EQ.NP$) JS=N$
      ALPHA(J,K)=CHGAMA(K)*F2
      STEMP(J,K)=CHGAMA(K)*X0(J,K,N)/RHO0(J,K)
c     1-CHTHET(K)*(T(JS,KS)-TG(KS))
      STEMP(J,K)=STEMP(J,K)*F2
1     CONTINUE
      IA=1
      IS=1
      RETURN
      ENDIF
      


      RETURN
      END

      SUBROUTINE PVSOR(STEMP)

C     GET THE SOURCE MINUS SINK TERM FOR CONSTITUENT N
C     STEMP IS THE SOURCE (MIXING RATIO/SEC)
C     alpha is 1/lifetime of species (1/SEC)
C     IS = 0 IF STEMP=0, IA=0 IF ALPHA=0

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMONT.INC'
      INCLUDE 'COMMOND.INC'
      DIMENSION STEMP(NP$,MP$)


     
      CALL DERV4B(2,GHEATX,DUM1 ,DZ,N$ ,M$, 2)  ! D(Q)/Dz
      CALL DERV4B(2,DRAGX ,DUM2X,DY,NP$,MP$,1)  ! D(X)/Dy
      CALL DERV4B(2,XC(1,1,1),DUM3X,DZ,NP$,MP$,2)  ! D(TH)/Dz
      CALL DERV4B(2,XC(1,1,2),DUM4X,DY,NP$,MP$,1)  ! D(M)/Dy

      CALL MRS2RBR(DUM1,DUM1X)
   

      DO K=1,MP$
      DO J=1,NP$
         STEMP(J,K)= 1.*(
     >      DUM1X(J,K)*DUM4X(J,K) + DUM2X(J,K)*DUM3X(J,K)
     >                 )*EXP( ZPP(K)/(H*1.0E-5) )
     >                 / CST(J)
      END DO
      END DO

      CALL DERV4B(2,GHEATX,DUM1 ,DY,N$ ,M$, 1)  ! D(Q)/Dy
      CALL DERV4B(2,DRAGX ,DUM2X,DZ,NP$,MP$,2)  ! D(X)/Dz
      CALL DERV4B(2,XC(1,1,1),DUM3X,DY,NP$,MP$,1)  ! D(TH)/Dy
      CALL DERV4B(2,XC(1,1,2),DUM4X,DZ,NP$,MP$,2)  ! D(M)/Dz

      CALL MRS2RBR(DUM1,DUM1X)
   

      DO K=1,MP$
      DO J=1,NP$
         STEMP(J,K)= STEMP(J,K) -1.*(
     >      DUM1X(J,K)*DUM4X(J,K) + DUM2X(J,K)*DUM3X(J,K)
     >                 )*EXP( ZPP(K)/(H*1.0E-5) )
     >                 / CST(J)
      END DO
      END DO

cc      CALL OUTA(20,NSTEP,NDAY,NYEAR,YPP,ZPP,STEMP,Np$,Mp$,0)


      RETURN
      END
      


      SUBROUTINE GETTV
 
C	THIS ROUTINE COMPUTES THE TRANSPORT VELOCITIES FOR THE
C       PRATHER SCHEME. FIRST THE W FIELDS ARE INTERPOLATED ONTO
C       A STAGGERED GRID, WRBR. THEN THE VRBR FIELDS ARE COMPUTED.
C	THE REASON THE ORIGINAL VS FIELDS ARE NOT CHOSEN IS THAT CONTINUITY
C       IN A GRID SENSE HAS TO BE ASSURED. - MRS	



      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'
      INCLUDE 'COMMONR.INC'
      
ccelf      include 'test_value.inc'
ccelf      logical nantest

ccelf      logical lisvalid
      
      COMMON/SLAKWK2/ RSPRIME(NP$,MP$),RSAVE(NP$,MP$)
     1,CRX(NPP$),WRX(NPP$)
	INTEGER IDAYS(13)
	DATA IDAYS/1,32,60,91,121,152,182,213,244,274,305,335,365/

c     some hadley circulation vertical velocity is added here

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
c
c      DO 113 J=1,N$
c      DO 113 K=1,10
c113   DUM1(J,K)=TEHE(J,K,M1) + 
c     1           (TEHE(J,K,M2)-TEHE(J,K,M1))*DX
c
C	DUM1 IS THE JACKMAN HEATING RATES
c
c
c      JTB ZEROED OUT DUM1 7/23/92
c
C      CALL GETGLM(DUM1,GHET,C,N$,M$)
	dr=1./86400.



C     INTERPOLATE WS FIELDS

      csum=0.
      do 899 j=1,np$
      csum=csum+cst(j)
899   continue
      csum=1./csum

      DO 601 K=1,Mp$
      DO 601 J=1,Np$
      VRBR(J,K)=0.
      RSPRIME(J,K)=RHO0(J,K)
      RSAVE(J,K)=0.
      SM(J,K)=RHO0(J,K)
601   CONTINUE

      DO 630 J=1,N$
630   CRX(J+1)=C(J)
      CRX(1)=0.
      CRX(NPP$)=0.

c     Hadley circulation is added in with the dum1 term

      ZP00=1.5*H
      DO 620 K=1,M$

      WGT=1.0
      IF (ZZ(K).GT.ZP00) WGT=EXP(-(ZZ(K)-ZP00)/H)

      DO 100 J=1,N$

c       nantest=test_nan( ws(j,k) )
c       if(nantest) then
c         PRINT*, "NAN in GETTV:"
c         PRINT*, "   J = ", J
c         PRINT*, "   K = ", K
c         STOP
c         ENDIF

      WRX(J+1)=C(J)*(WS(J,K) + dum1(j,k)*dr/dthdz(k))

C   w/ Linux INTEL FORTRAN, use ISNAN function  -  if X is NAN, then ISNAN is true
      
ccelf      lisvalid=test_valid( wrx(j+1) )

      if ( ISNAN ( wrx(j+1) ) ) then
        print *,' wrx not valid in gettv'
	print *,'j,k,wrx,c,ws,dum1,dr,dthdz = ',j,k,wrx(j+1),c(j),
     *    ws(j,k),dr,dthdz(k)
        stop	
      endif
100   continue
      WRX(1)=0.
      WRX(NPP$)=0.
      DO 621 J=1,NP$
      WRBR(J,K)=(WRX(J)+WRX(J+1))*.5/CST(J)

ccelf      lisvalid=test_valid( wrbr(j,k) )
ccelf      if(.not. lisvalid) then

      if ( ISNAN ( wrbr(j,k) ) ) then
        print *,' wrbr(j,k) not valid in gettv'
	print *,'j,k,wrbr,wrx(j),wrx(j+1),cst = ',j,k,wrbr(j,k),wrx(j),
     *    wrx(j+1),cst(j)
        stop  	
      endif
621   continue

c     now correct for mass problems (2 passes)

      do 8899 ikx=1,2
      wsum=0.
      do 623 j=1,np$     
      wsum=wrbr(j,k)*cst(j)+wsum      
623   continue
      werr=wsum*csum
      do 624 j=1,np$
      wrbr(j,k)=wrbr(j,k)-werr
624   continue
8899  continue     

620   CONTINUE


      DO 632 J=1,NP$
632   WRBR(J,MP$)=0.


c     fix up continuity 

         XPDZ=EXP(+DZ*.5/H)
         XMDZ=EXP(-DZ*.5/H)
         dyz=dy/dz
         DTZ=DT/DZ
         DYT=DY/DT

c first do south polar boundary

	VRBR(1,1)=-WRBR(1,1)*XMDZ*CST(1)*DYz/CRX(2)
	VRBR(1,MP$)=WRBR(1,M$)*XPDZ*CST(1)*DYZ/CRX(2)
	do ik=2,M$ 
		VRBR(1,ik)=-dyz*CST(1)*(WRBR(1,ik)*XMDZ
     c			-WRBR(1,ik-1)*XPDZ)/CRX(2)
	enddo
c now do the rest of the v's
	
	
	do ij=2,N$
c lower boundary
		VRBR(ij,1)=-WRBR(ij,1)*XMDZ*CST(ij)*DYz
     c			/CRX(ij+1)+VRBR(ij-1,1)*CRX(ij)
     c			/CRX(ij+1)

c upper boundary
		VRBR(ij,MP$)=+WRBR(ij,M$)*xpDZ*CST(ij)*dyz
     c			/CRX(ij+1)+VRBR(ij-1,MP$)
     c			*CRX(ij)/CRX(ij+1)
c any level
	do ik=2,M$
	VRBR(ij,ik)=-(WRBR(ij,ik)*xMDZ-WRBR(ij,ik-1)*XPDZ
     c		)*CST(ij)/CRX(ij+1)*DYz +
     c		VRBR(ij-1,ik)*CRX(ij)/CRX(ij+1)
     c	
	enddo
	enddo

c correct these to give continuity using the 
c Prather scheme for a 

1125    continue    
	do il=1,NP$
		if(WRBR(il,1) .ge. 0.) then
			RSPRIME(il,1)=RHO0(il,1)-WRBR(il,1)*dtZ
     c			*RHO0(il,1)*XMDZ
			RSAVE(il,2)=RSAVE(il,2)+WRBR(il,1)*dtZ
     c			*RHO0(il,1)*XMDZ
			else
			RSAVE(il,1)=RSAVE(il,1)-WRBR(il,1)*dtZ
     c			*RHO0(il,2)*XPDZ
			RSPRIME(il,2)=RHO0(il,2)+WRBR(il,1)*dtZ
     c			*RHO0(il,2)*XPDZ
		endif
	enddo
	do ik=2,M$
		do il=1,NP$
			if(WRBR(il,ik) .ge. 0.) then
			if(WRBR(il,ik-1) .ge. 0.) then
		RSPRIME(il,ik)=RHO0(il,ik)-WRBR(il,ik)*dtZ
     c			*RHO0(il,ik)*xMDZ
		RSAVE(il,ik+1)=RSAVE(il,ik+1)+WRBR(il,ik)*dtZ
     c			*RHO0(il,ik)*XMDZ
			else
		RSAVE(il,ik+1)=RSAVE(il,ik+1)+WRBR(il,ik)*dtZ
     c			*RSPRIME(il,ik)*XMDZ
		RSPRIME(il,ik)=RSPRIME(il,ik)-WRBR(il,ik)*dtZ
     c			*RSPRIME(il,ik)*XMDZ
			endif
			else
		RSAVE(il,ik)=RSAVE(il,ik)-WRBR(il,ik)*dtZ
     c			*RHO0(il,ik+1)*XPDZ
		RSPRIME(il,ik+1)=RHO0(il,ik+1)+WRBR(il,ik)*dtZ
     c			*RHO0(il,ik+1)*XPDZ
			endif
	enddo
	enddo

	do ik=1,MP$
		do il=1,NP$
		RSPRIME(il,ik)=RSPRIME(il,ik)+RSAVE(il,ik)
		enddo
	enddo

	do ik=1,MP$
		if(VRBR(1,ik) .ge. 0) then
		VRBR(1,ik)=(RSPRIME(1,ik)-RHO0(1,ik))*CST(1)*
     c			DYT/(RSPRIME(1,ik)*CRX(2))
		else
		VRBR(1,ik)=(RSPRIME(1,ik)-RHO0(1,ik))*CST(2)*
     c			DYT/(RSPRIME(2,ik)*CRX(2))
		endif
		do ij=2,N$
c both positive

	if(VRBR(ij-1,ik) .ge. 0. .and. VRBR(ij,ik) .ge. 0.) then
			VRBR(ij,ik)=(RSPRIME(ij,ik)-RHO0(ij,ik)
     c			+VRBR(ij-1,ik)*RSPRIME(ij-1,ik)*CRX(ij)*dt
     c			/CST(ij-1)/DY)*DY*CST(ij)/
     c			(RSPRIME(ij,ik)*CRX(ij+1)*dt)
		endif
c both negative

	if(VRBR(ij-1,ik) .lt. 0. .and. VRBR(ij,ik) .lt. 0.) then
			VRBR(ij,ik)=(RSPRIME(ij,ik)-RHO0(ij,ik)
     c			+VRBR(ij-1,ik)*RSPRIME(ij,ik)*CRX(ij)*dt
     c			/CST(ij)/DY)*CST(ij+1)*DY/
     c			(RSPRIME(ij+1,ik)*CRX(ij+1)*dt)
		endif
c <- and -> (this one is hardest)

	if(VRBR(ij-1,ik) .lt. 0. .and. VRBR(ij,ik) .ge. 0.) then
			RSDP=RSPRIME(ij,ik)+VRBR(ij-1,ik)*dt*CRX(ij)
     c			/CST(ij)*RSPRIME(ij,ik)/DY
			VRBR(ij,ik)=(RSDP-RHO0(ij,ik))/RSDP*DY
     c			*CST(ij)/dt/CRX(ij+1)
		endif
c -> and <-

	if(VRBR(ij-1,ik) .ge. 0. .and. VRBR(ij,ik) .lt. 0.) then
			VRBR(ij,ik)=(RSPRIME(ij,ik)-RHO0(ij,ik)
     c			+VRBR(ij-1,ik)*RSPRIME(ij-1,ik)*CRX(ij)*dt/
     c			(DY*CST(ij-1)))*DY*CST(ij+1)/
     c			(RSPRIME(ij+1,ik)*CRX(ij+1)*dt)
		endif
		enddo
	enddo

C***
96    continue

      DO J=1,MP$
      DO I=1,NP$
         VRB0(I,J)=VRBR(I,J)
         WRB0(I,J)=WRBR(I,J)
      END DO
      END DO

c      CALL W_AERO( 1.0, DUM1X )



      RETURN
      END
