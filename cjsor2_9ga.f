      SUBROUTINE CJSOR2(AX,AY,B,CX,CY,F,PSI,N,M,ERROR,ERRMAX,
     2 ERRET,NINIT)
      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      COMPLEX AX(N,M), AY(N,M), B(N,M), CX(N,M), CY(N,M)
      COMPLEX PSI(N,M), F(N,M)
      COMPLEX PSI0(N$,M$)
      REAL    PSIR(N$,M$)
      COMPLEX B1(N$,M$)   ,AVGPSI, ALPH
      COMMON /XTIME/ NSTEP,NDAY,NYEAR,DECL,IT1,IT2,IT3,MONTH,IDAYX,
     2               imodul
      COMMON /DUMP_1/ IDUMP,IFILE,isw2(isw$)
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
C  IF NINIT \= 1  => USES INITIAL GUESS PROVIDED BY PSI INPUT.     
C  IF NINIT  = 1  => USES INITIAL GUESS PSI=F   
C                           
C1234sTTTTT


      
c      IDUMP=10

c      NDAY=1
c      NYEAR=1



      D_LAT=YP(2)-YP(1)

c      write(6,*) '  Delta Lat. ', d_lat

      DO K=1,M$
      DO J=1,N$
         PSIR(J,K)=REAL( F(J,K) )
      END DO
      END DO


      OMEGA_OPT= 0.5 - (PI / (2.0*SQRT(2.0) ) )
     >          * SQRT( (1./N)**2 + (1./M)**2 )

      ALPH_OPT=4.*OMEGA_OPT


c      WRITE(6,*) ' OMEGA_OPT ',OMEGA_OPT
c      WRITE(6,*) 'PI ',PI


      DO 10 J=1,M
      DO 10 I=1,N
         B1(I,J)=1./B(I,J)
cc         IF ( ABS(F(I,J)) .LE. 1.E-15) THEN 
cc            F(I,J)=0.
cc         ENDIF
         IF ( NINIT .EQ. 1) THEN 
            PSI(I,J)=F(I,J)
         ENDIF
 10   CONTINUE

         DO J=1,M
         DO I=1,N
            IF ( ABS( YP(I) ) .LE. 1.5*D_LAT ) THEN
               PSI(I,J)=CMPLX( 0.00, 0.00 )
                 F(I,J)=CMPLX( 0.00, 0.00 )
                AX(I,J)=CMPLX( 0.00, 0.00 )
                CX(I,J)=CMPLX( 0.00, 0.00 )
                AY(I,J)=CMPLX( 0.00, 0.00 )
                CY(I,J)=CMPLX( 0.00, 0.00 )
            END IF
         END DO
         END DO
C- - - - -


      DO K=1,M$
      DO J=1,N$
         PSIR(J,K)=REAL( PSI(J,K) )
      END DO
      END DO

      DO K=1,M$
      DO J=1,N$
         PSIR(J,K)=AIMAG( PSI(J,K) )
      END DO
      END DO


      NITER=0
1000  CONTINUE ! ITERATION LOOP

      DO 90 J=1,M
      DO 90 I=1,N
         PSI0(I,J)=PSI(I,J)
 90   CONTINUE


      IF(NITER.LE.5) ALPH=CMPLX( 1.000, 0.000 )  
      IF(NITER.GT.5) ALPH=CMPLX( ALPH_OPT, 0.000 )  !OVER-RELAX. COEFF.


c      WRITE( 6, *)'     -----> ', ALPH, NITER

C - - RELAX INTERIOR POINTS 

      DO 100 J=2,M-1
      DO 100 I=2,N-1

         PSI(I,J) = 
     >       PSI(I,J) - ALPH*
     >     ( AY(I,J)*PSI(I,J+1) +AX(I,J)*PSI(I+1,J) 
     >     + CY(I,J)*PSI(I,J-1) +CX(I,J)*PSI(I-1,J) 
     >     + B(I,J) *PSI(I,J)   -F(I,J)     )*B1(I,J)
 100  CONTINUE

C - - ZERO BC'S AROUND EDGES

      DO I=2,N-1
         J=1
         PSI(I,J) = 
     >       PSI(I,J) - ALPH*
     >     ( AY(I,J)*PSI(I,J+1) +AX(I,J)*PSI(I+1,J) 
     >     + CY(I,J)*0.0000     +CX(I,J)*PSI(I-1,J) 
     >     + B(I,J) *PSI(I,J)   -F(I,J)     )*B1(I,J)
         J=M
         PSI(I,J) = 
     >       PSI(I,J) - ALPH*
     >     ( AY(I,J)*0.0000     +AX(I,J)*PSI(I+1,J) 
     >     + CY(I,J)*PSI(I,J-1) +CX(I,J)*PSI(I-1,J) 
     >     + B(I,J) *PSI(I,J)   -F(I,J)     )*B1(I,J)
      END DO
      DO J=2,M-1
         I=1
         PSI(I,J) = 
     >       PSI(I,J) - ALPH*
     >     ( AY(I,J)*PSI(I,J+1) +AX(I,J)*PSI(I+1,J) 
     >     + CY(I,J)*PSI(I,J-1) +CX(I,J)*0.0000
     >     + B(I,J) *PSI(I,J)   -F(I,J)     )*B1(I,J)
         I=N
         PSI(I,J) = 
     >       PSI(I,J) - ALPH*
     >     ( AY(I,J)*PSI(I,J+1) +AX(I,J)*0.0000
     >     + CY(I,J)*PSI(I,J-1) +CX(I,J)*PSI(I-1,J) 
     >     + B(I,J) *PSI(I,J)   -F(I,J)     )*B1(I,J)
      END DO
         I=1
         J=1
         PSI(I,J) = 
     >       PSI(I,J) - ALPH*
     >     ( AY(I,J)*PSI(I,J+1) +AX(I,J)*PSI(I+1,J) 
     >     + CY(I,J)*0.0000     +CX(I,J)*0.0000
     >     + B(I,J) *PSI(I,J)   -F(I,J)     )*B1(I,J)
         I=N
         J=1
         PSI(I,J) = 
     >       PSI(I,J) - ALPH*
     >     ( AY(I,J)*PSI(I,J+1) +AX(I,J)*0.0000 
     >     + CY(I,J)*0.0000     +CX(I,J)*PSI(I-1,J) 
     >     + B(I,J) *PSI(I,J)   -F(I,J)     )*B1(I,J)
         I=1
         J=M
         PSI(I,J) = 
     >       PSI(I,J) - ALPH*
     >     ( AY(I,J)*0.0000     +AX(I,J)*PSI(I+1,J) 
     >     + CY(I,J)*PSI(I,J-1) +CX(I,J)*0.0000 
     >     + B(I,J) *PSI(I,J)   -F(I,J)     )*B1(I,J)
         I=N
         J=M
         PSI(I,J) = 
     >       PSI(I,J) - ALPH*
     >     ( AY(I,J)*0.0000     +AX(I,J)*0.0000 
     >     + CY(I,J)*PSI(I,J-1) +CX(I,J)*PSI(I-1,J) 
     >     + B(I,J) *PSI(I,J)   -F(I,J)     )*B1(I,J)
C- - - - -
C- - - - -
C  ZERO OUT PLANETARY WAVE STREAMFUNCTION IN EQUATORIAL ZONE
         

         

         DO J=1,M
         DO I=1,N
            IF ( ABS( YP(I) ) .LE. 1.5*D_LAT ) THEN
               PSI(I,J)=CMPLX( 0.00, 0.00 )
            END IF
         END DO
         END DO
C- - - - -
  

      AVGPSI=CMPLX(0.,0.)
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
     >        (CABS(PSI(I,J)-PSI0(I,J)))**2  
         TNORM= 
     >        TNORM + 
     >        (CABS(PSI0(I,J)-AVGPSI))**2
 102  CONTINUE
 
      ERROR1=DIFF/TNORM

      NITER=NITER+1
c      NDAY= INT(NITER/10)

      DO K=1,M$
      DO J=1,N$
         PSIR(J,K)=REAL( PSI(J,K) )
      END DO
      END DO

      DO K=1,M$
      DO J=1,N$
         PSIR(J,K)=AIMAG( PSI(J,K) )
      END DO
      END DO


c      IF( MOD(NITER,20).EQ.0) WRITE(*,*) 'IN CJSOR ',NITER,ERROR1
      IF (NITER .GT. 10000) GO TO 1001   
      IF( (NITER .LT. 10) .OR. (ERROR1 .GT. 1.E-8)) GO TO 1000
    

 1001 CONTINUE

      ERROR=ERROR1
      ERRMAX=TNORM
      ERRET=DIFF    


c       write(6,*) ' IN CJSOR ' , error,errmax,erret,niter

 


      RETURN
      END

