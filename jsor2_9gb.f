C
C    for 9GB run, use X* BBC from NCEP reanalysis, interpolated to coupled model dynamics grid  - BBC1(N$,4)
C                                                                         NOTE: N=N$,  M=M$  here
C
      SUBROUTINE JSOR2(AX,AY,B,CX,CY,F,PSI,N,M,ERROR,ERRMAX,ERRET,NINIT)
      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      
ccelf      include 'test_value.inc'
ccelf     logical nantest

      
      DIMENSION AX(N,M), AY(N,M), B(N,M), CX(N,M), CY(N,M)
      DIMENSION PSI(N,M), F(N,M)
      DIMENSION PSI0(N$,M$)
      DIMENSION B1(N$,M$)
      COMMON /XTIME/ NSTEP,NDAY,NYEAR,DECL,IT1,IT2,IT3,MONTH,IDAYX,
     2               imodul

C                                   BBC1(N$,4) is in m2/sec ,  PSI(N$,M$) is in cm2/sec
      COMMON/CBBC1/BBC1(N$,4)


C     Added 9/18/96 by PEM to flag initialization of PSI, without using
C     NINIT (NSTEP), which has changed function under new timing (timer.f.

      LOGICAL FLAGPSI
      DATA    FLAGPSI /.FALSE./

C     End Added 9/18/96.


C
C----------------------------------------------------------------
C
C  JSOR - A STUPID OVERRELAXATION SOLUTION OF A FORCED, ELLITPIC
C           EQUATION
C
C
C  USES 5-POINT OPERATOR:
C                                + (AY)
C                               
c
C                       + (CX)   x (B, F)  + (AX)
C                             
c
C                                + (CY)                                
C                                                             
C                                                             
C  ASSUMES PSI=0 ALONG ALL BOUNDARIES (+_ DZ,DY FROM J=1,M I=1,N).

C  Note 9/18/96 by PEM.  NINIT is no longer used with the new timing (timer.f).
C  The function of NINIT is now performed by PSI_INIT.

C  IF NINIT \= 1  => USES INITIAL GUESS PROVIDED BY PSI INPUT.     
C  IF NINIT  = 1  => USES INITIAL GUESS PSI=F   
C                           
C1234sTTTTT



      OMEGA_OPT= 0.5 - (PI / (2.0*SQRT(2.0) ) )
     >          * SQRT( (1./N)**2 + (1./M)**2 )

      ALPH_OPT=4.*OMEGA_OPT


      DO 10 J=1,M
      DO 10 I=1,N
         B1(I,J)=1./B(I,J)
cc         IF ( ABS(F(I,J)) .LE. 1.E-15) THEN 
cc            F(I,J)=0.
cc         ENDIF

C     Modified 9/18/96 by PEM to use PSI_INIT instead of NINIT, consistent
C     with new timing (timer.f).

C         IF ( NINIT .EQ. 1) THEN 
C            PSI(I,J)=F(I,J)
C         ENDIF

          IF (.NOT.FLAGPSI) THEN
C             PRINT*, "Initializing PSI ..."
             PSI(I, J) = F(I, J)
             ENDIF

C     End Modified 9/18/96.

 10   CONTINUE


C     Added 9/20/96 by PEM.

          IF (.NOT.FLAGPSI) FLAGPSI = .TRUE.

C     End Added 9/20/96.


      CALL PRMINMAX(B ,'-B-JSOR --',N$, M$,NDBUG) 
      CALL PRMINMAX(F ,'-F-JSOR --',N$, M$,NDBUG) 
      CALL PRMINMAX(CY,'CY-JSOR --',N$, M$,NDBUG) 
      CALL PRMINMAX(AY,'AY-JSOR --',N$, M$,NDBUG) 
      CALL PRMINMAX(CX,'CX-JSOR --',N$, M$,NDBUG) 
      CALL PRMINMAX(AX,'AX-JSOR --',N$, M$,NDBUG) 

c      CALL OUTA(20,NSTEP,NDAY,NYEAR,YPP,ZPP,PSI,N,M,0)
c      CALL OUTA(20,NSTEP,NDAY,NYEAR,YPP,ZPP,Ax,N,M,0)
c      CALL OUTA(20,NSTEP,NDAY,NYEAR,YPP,ZPP,Cx,N,M,0)
c      CALL OUTA(20,NSTEP,NDAY,NYEAR,YPP,ZPP,Ay,N,M,0)
c      CALL OUTA(20,NSTEP,NDAY,NYEAR,YPP,ZPP,Cy,N,M,0)
c      CALL OUTA(20,NSTEP,NDAY,NYEAR,YPP,ZPP,B,N,M,0)
c      CALL OUTA(20,NSTEP,NDAY,NYEAR,YPP,ZPP,F,N,M,0)


      NITER=0
      
      
1000  CONTINUE ! ITERATION LOOP

      DO 90 J=1,M
      DO 90 I=1,N
         PSI0(I,J)=PSI(I,J)
 90   CONTINUE


CCelf
ccelf  set BBC here for initiailization on PSI0(N$=N,M$=M);  -  just do bottom level (2 km)
Ccelf                                        BBC1(N$,4) is for current day, in m2/sec, convert to cm2/sec
CCELF                        
CCELF                           BBC1 here DOES NOT WORK IN THE COUPLED MODEL!!!!
CCELF       ibbct = 1
CCELF
CCELF       do 95 J=1,ibbct
CCELF       do 95 I=1,N
CCELF         PSI0(I,J) = BBC1(I,J)*1.e4
CCELF 95   CONTINUE
CCelf


      IF(NITER.LE.5) ALPH=1.0000  
      IF(NITER.GT.5) ALPH=ALPH_OPT  !OVER-RELAX. COEFF.


C - - RELAX INTERIOR POINTS 

      DO 100 J=2,M-1
      DO 100 I=2,N-1

         PSI(I,J) = 
     >       PSI(I,J) - ALPH*
     >     ( AY(I,J)*PSI(I,J+1) +AX(I,J)*PSI(I+1,J) 
     >     + CY(I,J)*PSI(I,J-1) +CX(I,J)*PSI(I-1,J) 
     >     + B(I,J) *PSI(I,J)   -F(I,J)     )*B1(I,J)

c       nantest=test_nan( psi(i,j) )
c       if(nantest) then
c            PRINT*, "PSI NAN in JSOR2: 1"
c            PRINT*, "    I = ", I
c            PRINT*, "    J = ", J
c            PRINT*, "   AY = ", AY(I,J)
c            PRINT*, "   AX = ", AX(I,J)
c            PRINT*, "   CY = ", CY(I, J)
c            PRINT*, "   CX = ", CX(I, J)
c            PRINT*, "    B = ", B(I, J)
c            PRINT*, "    F = ", F(I, J)
c            PRINT*, "   B1 = ", B1(I, J)
c            STOP
c            ENDIF

 100  CONTINUE

C - - ZERO BC'S AROUND EDGES

      DO I=2,N-1
         J=1
         PSI(I,J) = 
     >       PSI(I,J) - ALPH*
     >     ( AY(I,J)*PSI(I,J+1) +AX(I,J)*PSI(I+1,J) 
     >     + CY(I,J)*0.0000     +CX(I,J)*PSI(I-1,J) 
     >     + B(I,J) *PSI(I,J)   -F(I,J)     )*B1(I,J)

c         IF (PSI(I,J) .NE. PSI(I,J)) THEN
c       nantest=test_nan( psi(i,j) )
c       if(nantest) then
c            PRINT*, "PSI NAN in JSOR2: 2"
c            PRINT*, "    I = ", I
c            PRINT*, "    J = ", J
c            PRINT*, "   AY = ", AY(I,J)
c            PRINT*, "   AX = ", AX(I,J)
c            PRINT*, "   CY = ", CY(I, J)
c            PRINT*, "   CX = ", CX(I, J)
c            PRINT*, "    B = ", B(I, J)
c            PRINT*, "    F = ", F(I, J)
c            PRINT*, "   B1 = ", B1(I, J)
c            STOP
c            ENDIF

         J=M
         PSI(I,J) = 
     >       PSI(I,J) - ALPH*
     >     ( AY(I,J)*0.0000     +AX(I,J)*PSI(I+1,J) 
     >     + CY(I,J)*PSI(I,J-1) +CX(I,J)*PSI(I-1,J) 
     >     + B(I,J) *PSI(I,J)   -F(I,J)     )*B1(I,J)

c       nantest=test_nan( psi(i,j) )
c       if(nantest) then
c            PRINT*, "PSI NAN in JSOR2: 3"
c            PRINT*, "    I = ", I
c            PRINT*, "    J = ", J
c            PRINT*, "   AY = ", AY(I,J)
c            PRINT*, "   AX = ", AX(I,J)
c            PRINT*, "   CY = ", CY(I, J)
c            PRINT*, "   CX = ", CX(I, J)
c            PRINT*, "    B = ", B(I, J)
c            PRINT*, "    F = ", F(I, J)
c            PRINT*, "   B1 = ", B1(I, J)
c            STOP
c            ENDIF

      END DO
      DO J=2,M-1
         I=1
         PSI(I,J) = 
     >       PSI(I,J) - ALPH*
     >     ( AY(I,J)*PSI(I,J+1) +AX(I,J)*PSI(I+1,J) 
     >     + CY(I,J)*PSI(I,J-1) +CX(I,J)*0.0000
     >     + B(I,J) *PSI(I,J)   -F(I,J)     )*B1(I,J)

c       nantest=test_nan( psi(i,j) )
c       if(nantest) then
c            PRINT*, "PSI NAN in JSOR2: 4"
c            PRINT*, "    I = ", I
c            PRINT*, "    J = ", J
c            PRINT*, "   AY = ", AY(I,J)
c            PRINT*, "   AX = ", AX(I,J)
c            PRINT*, "   CY = ", CY(I, J)
c            PRINT*, "   CX = ", CX(I, J)
c            PRINT*, "    B = ", B(I, J)
c            PRINT*, "    F = ", F(I, J)
c            PRINT*, "   B1 = ", B1(I, J)
c            STOP
c            ENDIF

         I=N
         PSI(I,J) = 
     >       PSI(I,J) - ALPH*
     >     ( AY(I,J)*PSI(I,J+1) +AX(I,J)*0.0000
     >     + CY(I,J)*PSI(I,J-1) +CX(I,J)*PSI(I-1,J) 
     >     + B(I,J) *PSI(I,J)   -F(I,J)     )*B1(I,J)

c       nantest=test_nan( psi(i,j) )
c       if(nantest) then
c            PRINT*, "PSI NAN in JSOR2: 5"
c            PRINT*, "    I = ", I
c            PRINT*, "    J = ", J
c            PRINT*, "   AY = ", AY(I,J)
c            PRINT*, "   AX = ", AX(I,J)
c            PRINT*, "   CY = ", CY(I, J)
c            PRINT*, "   CX = ", CX(I, J)
c            PRINT*, "    B = ", B(I, J)
c            PRINT*, "    F = ", F(I, J)
c            PRINT*, "   B1 = ", B1(I, J)
c            STOP
c            ENDIF

      END DO
         I=1
         J=1
         PSI(I,J) = 
     >       PSI(I,J) - ALPH*
     >     ( AY(I,J)*PSI(I,J+1) +AX(I,J)*PSI(I+1,J) 
     >     + CY(I,J)*0.0000     +CX(I,J)*0.0000
     >     + B(I,J) *PSI(I,J)   -F(I,J)     )*B1(I,J)

c       nantest=test_nan( psi(i,j) )
c       if(nantest) then
c            PRINT*, "PSI NAN in JSOR2: 6"
c            PRINT*, "    I = ", I
c            PRINT*, "    J = ", J
c            PRINT*, "   AY = ", AY(I,J)
c            PRINT*, "   AX = ", AX(I,J)
c            PRINT*, "   CY = ", CY(I, J)
c            PRINT*, "   CX = ", CX(I, J)
c            PRINT*, "    B = ", B(I, J)
c            PRINT*, "    F = ", F(I, J)
c            PRINT*, "   B1 = ", B1(I, J)
c            STOP
c            ENDIF

         I=N
         J=1
         PSI(I,J) = 
     >       PSI(I,J) - ALPH*
     >     ( AY(I,J)*PSI(I,J+1) +AX(I,J)*0.0000 
     >     + CY(I,J)*0.0000     +CX(I,J)*PSI(I-1,J) 
     >     + B(I,J) *PSI(I,J)   -F(I,J)     )*B1(I,J)

c        nantest=test_nan( psi(i,j) )
c        if(nantest) then
c            PRINT*, "PSI NAN in JSOR2: 7"
c            PRINT*, "    I = ", I
c            PRINT*, "    J = ", J
c            PRINT*, "   AY = ", AY(I,J)
c            PRINT*, "   AX = ", AX(I,J)
c            PRINT*, "   CY = ", CY(I, J)
c            PRINT*, "   CX = ", CX(I, J)
c            PRINT*, "    B = ", B(I, J)
c            PRINT*, "    F = ", F(I, J)
c            PRINT*, "   B1 = ", B1(I, J)
c            STOP
c            ENDIF

         I=1
         J=M
         PSI(I,J) = 
     >       PSI(I,J) - ALPH*
     >     ( AY(I,J)*0.0000     +AX(I,J)*PSI(I+1,J) 
     >     + CY(I,J)*PSI(I,J-1) +CX(I,J)*0.0000 
     >     + B(I,J) *PSI(I,J)   -F(I,J)     )*B1(I,J)

c       nantest=test_nan( psi(i,j) )
c       if(nantest) then
c            PRINT*, "PSI NAN in JSOR2: 8"
c            PRINT*, "    I = ", I
c            PRINT*, "    J = ", J
c            PRINT*, "   AY = ", AY(I,J)
c            PRINT*, "   AX = ", AX(I,J)
c            PRINT*, "   CY = ", CY(I, J)
c            PRINT*, "   CX = ", CX(I, J)
c            PRINT*, "    B = ", B(I, J)
c            PRINT*, "    F = ", F(I, J)
c            PRINT*, "   B1 = ", B1(I, J)
c            STOP
c            ENDIF

         I=N
         J=M
         PSI(I,J) = 
     >       PSI(I,J) - ALPH*
     >     ( AY(I,J)*0.0000     +AX(I,J)*0.0000 
     >     + CY(I,J)*PSI(I,J-1) +CX(I,J)*PSI(I-1,J) 
     >     + B(I,J) *PSI(I,J)   -F(I,J)     )*B1(I,J)
     

c       nantest=test_nan( psi(i,j) )
c       if(nantest) then
c            PRINT*, "PSI NAN in JSOR2: 9"
c            PRINT*, "    I = ", I
c            PRINT*, "    J = ", J
c            PRINT*, "   AY = ", AY(I,J)
c            PRINT*, "   AX = ", AX(I,J)
c            PRINT*, "   CY = ", CY(I, J)
c            PRINT*, "   CX = ", CX(I, J)
c            PRINT*, "    B = ", B(I, J)
c            PRINT*, "    F = ", F(I, J)
c            PRINT*, "   B1 = ", B1(I, J)
c            STOP
c            ENDIF

C- - - - -


ccelf
CCelf
ccelf   reset BBC here  BBC1(N$,4) is in m2/sec, convert to cm2/sec; just OVERWRITE PSI(N$=N,M$=M)
CCELF                           - just do bottom level (2 km)  - THIS DOES NOT WORK IN THE COUPLED MODEL!!!!
CCELF       do 396 J=1,ibbct
CCELF       do 396 I=1,N
CCELF         PSI(I,J) = BBC1(I,J)*1.e4
CCELF 396   CONTINUE
CCelf


      AVGPSI=0.
      TNORM=0.
      DIFF=0.

      DO 101 J=1,M
      DO 101 I=1,N
         AVGPSI= 
     >        AVGPSI + 
     >        PSI0(I,J)
 101  CONTINUE
      AVGPSI=AVGPSI/(1.*N*M)
  
      DO 102 J=1,M
      DO 102 I=1,N
         DIFF= 
     >        DIFF + 
     >        (PSI(I,J)-PSI0(I,J))**2  
         TNORM= 
     >        TNORM + 
     >        (PSI0(I,J)-AVGPSI)**2
 102  CONTINUE

c      DO I = 1, N
c         DO J = 1, M
c        nantest=test_nan( psi(i,j) )
c        if(nantest) then
c               PRINT*, "PSI NAN in JSOR2: 10"
c               PRINT*, "   I = ", I
c               PRINT*, "   J = ", J
c               STOP
c               ENDIF
c            END DO
c         END DO
 
      ERROR1=DIFF/TNORM

      NITER=NITER+1
c      IF( MOD(NITER,20).EQ.0) WRITE(*,*) 'IN JSOR ',NITER,ERROR1
      IF (NITER .GT. 10000) GO TO 1001   
      IF( (NITER .LT. 10) .OR. (ERROR1 .GT. 1.E-8)) GO TO 1000
    

 1001 CONTINUE

      ERROR=ERROR1
      ERRMAX=TNORM
      ERRET=DIFF    


c       write(6,*) ' IN JSOR ' , error,errmax,erret,niter

 
c      CALL OUTA(20,NSTEP,NDAY,NYEAR,YPP,ZPP,PSI,N,M,0)


      RETURN
      END

