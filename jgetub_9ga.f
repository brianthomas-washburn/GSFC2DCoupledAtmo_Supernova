      SUBROUTINE JGETUB

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'

ccelf      include 'test_value.inc'
ccelf      logical nantest


C     Added 9/18/96 by PEM for new timing (timer.f).

      include "timer.inc"

C     End Added 9/18/96 by PEM.


      DIMENSION UBEQ(M$),TYY(M$),STABP(N$,M$),RCOEFP(NP$)


      A2=1.00   !
      A4=1.00-A2
      XKAPPA=0.288

      DO K=1,M$
      DO J=1,N$
         DUM2(J,K)=(
     >              THX(J,K)  +  THX(J+1,K)
     >            + THX(J,K+1)+  THX(J+1,K+1) )/4.00
      END DO
      END DO

      DO K=1,M$
      DO J=1,N$
         DUM1(J,K)=DUM2(J,K)*EXP( -ZZ(K)*XKAPPA/H )
      END DO
      END DO


      RSP = R/H    

      DO K=1,M$
      DO J=2,N$-1
         DUM2(J,K)= 
     >     -RSP*( DUM1(J+1,K)-DUM1(J-1,K) )
     >     / (2*DY)   

      END DO
      END DO

      DO K=1,M$
      DO J=3,N$-2
         DUM2(J,K)= 
     >    A2*( -RSP*( DUM1(J+1,K)-DUM1(J-1,K) )
     >      / (2*DY)  )  
     >  + A4*( -RSP*( DUM1(J+2,K)-DUM1(J-2,K) )
     >     / (4*DY)  ) 
      END DO
      END DO
         
      DO K=1,M$
         J=1
         DUM2(J,K)= 
     >     -RSP*( DUM1(J+1,K)-DUM1(J,K) )
     >     / (2*DY)   
         J=N$
         DUM2(J,K)= 
     >     -RSP*( DUM1(J,K)-DUM1(J-1,K) )
     >     / (2*DY)   
      END DO



      DO J=1,N$
         K=1
         CFF=CF(J)
cc         UTHRM(J,1)= 0.
         UTHRM(J,K)= (0.5*DZ/CFF)*
     >               (DUM2(J,K))
      END DO

      DO K=2,M$
      DO J=1,N$

         CFF=CF(J)+2*S(J)*UTHRM(J,K-1)/(A*C(J))

         UTHRM(J,K)= UTHRM(J,K-1) + (DZ/CFF)*
     >               (DUM2(J,K))
      END DO
      END DO

C      DO K=1,M$
C      DO J=1,N$
C         IF(ABS(YP(J)) .LT. 15.) THEN
C           UTHRM(J,K)= 0.
C         ENDIF
C         IF(ABS(YP(J)) .GT. 185.) THEN
C           UTHRM(J,K)= 0.
C         ENDIF
C      END DO
C      END DO


cc      BDYLAT=20.
      BDYLAT=1.0*ISW(34)


      DO K=1,M$
      DO J=1,N$
         IF(ABS(YP(J)) .LT. 15.) THEN
           UTHRM(J,K)=0.00
          END IF
      END DO
      END DO
            


      DO K=1,M$
      DO J=1,N$
         IF(ABS(YP(J)) .LT. BDYLAT) THEN
           DUM1(J,K)=
     >   (  UTHRM(J+2,K) + 2*UTHRM(J+1,K)
     >      + 4*UTHRM(J,K)
     >     + 2*UTHRM(J-1,K) + UTHRM(J-2,K)  )/10.
          END IF
      END DO
      END DO
            

      DO K=1,M$
      DO J=1,N$
         IF(ABS(YP(J)) .LT. BDYLAT) THEN
           UTHRM(J,K)=DUM1(J,K)
          END IF
      END DO
      END DO
            
c      DO K=1,M$
c      DO J=1,N$
c         IF(ABS(YP(J)) .LT. BDYLAT) THEN
c           WT1=ABS(YP(J))/BDYLAT
c           WT2=1.-WT1
c           UTHRM(J,K)=WT1*UTHRM(J,K)   !+WT2*UB(J,K,2)
c          END IF
c      END DO
c      END DO
            
     
      IF(ISW(37).eq.1) THEN
      WRITE(*,*) ' SMOOTHING THERMAL WIND '

      DO K=1,M$
      DO J=2,N$-1
           DUM1(J,K)=
     >   (  UTHRM(J+2,K) + 2*UTHRM(J+1,K)
     >      + 4*UTHRM(J,K)
     >     + 2*UTHRM(J-1,K) + UTHRM(J-2,K)  )/10.
      END DO
      END DO

      DO K=1,M$
         J=N$
           DUM1(J,K)=
     >     (
     >      UTHRM(J,K)
     >     +UTHRM(J-1,K) )/2. 
         J=N$-1
           DUM1(J,K)=
     >     (UTHRM(J+1,K)
     >      + 2*UTHRM(J,K)
     >     +UTHRM(J-1,K) )/4. 
         J=2
           DUM1(J,K)=
     >     (UTHRM(J+1,K)
     >      + 2*UTHRM(J,K)
     >     +UTHRM(J-1,K) )/4. 
         J=1
           DUM1(J,K)=
     >     (UTHRM(J+1,K)
     >      + UTHRM(J,K)
     >     )/2. 
      END DO

      DO K=1,M$
      DO J=1,N$
         UTHRM(J,K)=DUM1(J,K)
      END DO
      END DO

      END IF


      CALL MRS2RBR(UTHRM,UTHRMX)


      IF(ISW(15) .EQ. 0) THEN
      WRITE(*,*) '_____ USING THERMAL WIND AS UBX  '

      CALL MRS2RBR(UTHRM,DUM1X)

      DO J=1,NP$
         IF(ABS(YPP(J)).GE.30) THEN
         RCOEFP(J)=
     >      1.000
         END IF
         IF(ABS(YPP(J)).LT.30) THEN
         RCOEFP(J)=
     >      1.00 - ABS( ABS( YPP(J) ) - 30. )/30.
         END IF


C     Modified 9/18/96 by PEM to use new timing (timer.f).

C         IF(NSTEP.EQ.1) WRITE(*,*) '___THERM WIND WT ',YPP(J),
C     >                             RCOEFP(J)

         IF(IDOR360.EQ.1) WRITE(*,*) '___THERM WIND WT ',YPP(J),
     >                             RCOEFP(J)

C     End Modified 9/18/96.


      END DO

      DO K=1,MP$
      DO J=1,NP$

         UPREV = UBX(J, K)

         UBX(J,K)=
     >       DUM1X(J,K)*RCOEFP(J)+
     >  ( 1.00-RCOEFP(J) )
     >                    *UBX(J,K)

c      IF (UBX(J, K) .NE. UBX(J, K)) THEN
c       nantest=test_nan(ubx(j,k))
c       if(nantest) then
c         PRINT*, "NAN in JGETUB:"
c         PRINT*, "   J = ", J
c         PRINT*, "   K = ", K
c         PRINT*, "   DUM1X = ", DUM1X(J, K)
c         PRINT*, "   UPREV = ", UPREV
c         STOP
c         ENDIF

      END DO
      END DO


      END IF





      return
      end
