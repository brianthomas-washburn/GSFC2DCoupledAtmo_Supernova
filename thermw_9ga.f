      SUBROUTINE THERMW

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'


C     Added 9/18/96 by PEM for new timing (timer.f).

      include "timer.inc"

C     End Added 9/18/96.


C     Modified 9/18/96 by PEM for new timing (timer.f).

C      IF (NSTEP .EQ. 1 ) THEN
C            WRITE(*,*) 'THERMAL WIND RELAXER ...'
C      ENDIF

c      IF (IDOR360 .EQ. 1 ) THEN
c            WRITE(*,*) 'THERMAL WIND RELAXER ...'
c            ENDIF

C     End Modified 9/18/96.


        MLIM=M$
        DO K=1,MLIM
        DO J=1,N$
           RCOEF2=ABS( CF(J) )
           RCOEF2=AMIN1( RCOEF2, ABS(CF(4)) )  !MAX OUT AT ABOUT 60
           THRMW(J,K) = -RCOEF2*
     >     ( UB(J,K,IT2)-UTHRM(J,K) )
        END DO
        END DO


        IF(ISW(31).GT.0) THEN
c        WRITE(*,*) ' CUTTING THRMW OFF AT LAT=',ISW(31)
        DO K=1,M$
        DO J=1,N$
           IF ( ABS(YP(J)) .LT. 1.0*ISW(31) ) THEN
              THRMW(J,K) = 0.000
           END IF
        END DO
        END DO
        END IF


        IF(ISW(36).EQ.1) THEN
c        WRITE(*,*) ' SMOOTHING THRMW ... '

        DO K=2,M$-1
        DO J=1,N$

           DUM1(J,K)=( THRMW(J,K+1) + THRMW(J,K-1) 
     >            + 2.*THRMW(J,K) )/4.00

        END DO
        END DO  

        DO J=1,N$
           K=M$
           DUM1(J,K)=(             THRMW(J,K-1) 
     >            + 2.*THRMW(J,K) )/3.00
           K=1
           DUM1(J,K)=( THRMW(J,K+1)  
     >            + 2.*THRMW(J,K) )/3.00
        END DO  


        DO K=1,M$

           J=1
           DUM1(J,K)=( THRMW(J,K)*0.5+THRMW(J+1,K)*0.5  )
           J=N$
           DUM1(J,K)=( THRMW(J,K)*0.5+THRMW(J-1,K)*0.5  )

        END DO


        DO K=1,M$
        DO J=1,N$

           THRMW(J,K) = DUM1(J,K)

        END DO
        END DO  

        END IF
      

      RETURN
      END

     





