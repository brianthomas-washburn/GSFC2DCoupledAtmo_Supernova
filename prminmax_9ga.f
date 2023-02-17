      SUBROUTINE PRMINMAX( A, A_LAB, N, M, NDBUG)
      DIMENSION A(N,M)
      CHARACTER*10 A_LAB


      RETURN
      IF( NDBUG .EQ. 0 ) RETURN


      I1=0
      J1=0
      I0=0
      J0=0

      AMIN = 10.E10
      AMAX = -10.E10

      DO 1 J=1,M
      DO 1 I=1,N
         IF(A(I,J) .LT. AMIN) THEN
            I0=I
            J0=J
         ENDIF
         AMIN=AMIN1( A(I,J), AMIN)
         IF(A(I,J) .GT. AMAX) THEN
            I1=I
            J1=J
         ENDIF
         AMAX=AMAX1( A(I,J), AMAX)
 1    CONTINUE

      WRITE(*,*) A_LAB, 'MIN ', AMIN, I0, J0
      WRITE(*,*) A_LAB, 'MAX ', AMAX, I1, J1
      WRITE(*,*)' *********** '

      RETURN  
      END
     



