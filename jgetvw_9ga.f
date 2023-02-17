      SUBROUTINE JGETVW(PSI,VST,WST,N,M)

C     CHANGE PSI TO V AND W

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'

ccelf      INCLUDE 'test_value.inc'
ccelf      logical nantest
      
      DIMENSION PSI(N,M) , VST(N,M) , WST(N,M) ,DM1(Np$,Mp$)

C w/ Linux INTEL FORTRAN, use ISNAN function  -  if PSI(ij,ik) in NAN, then ISNAN is true

      DO K = 1, M
         DO J = 1, N
            IF (ISNAN( PSI(J, K) ) ) THEN
               PRINT*, "Input PSI NAN in JGETVW:"
               PRINT*, "   J = ", J
               PRINT*, "   K = ", K
               STOP
               ENDIF
            END DO
         END DO

      DO 100 K=1,M  
      DO 100 J=2,N-1

      WST(J,K)=
     >  (PSI(J+1,K)-PSI(J-1,K))
     >       / ( C(J)*2.*DY )
100   CONTINUE

      DO 102 K=2,M-1
      DO 102 J=1,N

      VST(J,K)=
     >   -(PSI(J,K+1)-PSI(J,K-1))
     >       / ( C(J)*2.*DZ) 
     >    + PSI(J,K)/( H*C(J) )
102   CONTINUE


      DO 2 K=1,M
      J=1
      WST(J,K)=
     >   (PSI(J+1,K))
     >       / ( C(J)*2*DY )
      J=N
      WST(J,K)=
     >   ( -PSI(J-1,K))
     >       / ( C(J)*2*DY )
2     CONTINUE

      DO 3 J=1,N
      K=1
      VST(J,K)=
     >    -(PSI(J,K+1))
     >       / ( C(J)*2*DZ) 
     >    + PSI(J,K)/( H*C(J) )

      K=M
      VST(J,K)=
     >    -(-PSI(J,K-1))
     >       / ( C(J)*2*DZ) 
     >    + PSI(J,K)/( H*C(J) )

3     CONTINUE


C---KLUGE!!!-KLUGE!!!-KLUGE-------------------
c !SMOOTH W FIELD IN HORIZONTAL
c1234st
      DO K=1,M
      DO J=2,N-1
         DM1(J,K) = 0.25*WST(J-1,K) 
     >            + 0.50*WST(J,K) 
     >            + 0.25*WST(J+1,K)
      END DO
      END DO
     
      DO K=1,M
         DM1(1,K) = 
     >            + 0.66666*WST(1,K) 
     >            + 0.33333*WST(2,K)
         DM1(N,K) = 
     >            + 0.66666*WST(N,K) 
     >            + 0.33333*WST(N-1,K)
      END DO

      DO K=1,M
      DO J=2,N-1
         WST(J,K) = DM1(J,K) 
      END DO
      END DO

      DO J = 1, N
         DO K = 1, M
            IF (ISNAN( WST(J,K) ) ) THEN
               PRINT*, "Output WST NAN in JGETVW:"
               PRINT*, "   J = ", J
               PRINT*, "   K = ", K
               STOP
               ENDIF
            END DO
         END DO

C-------------------------------------------------


      RETURN
      END
