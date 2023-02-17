      SUBROUTINE XGWDRG( UB , TH , Y , ZPZ , N , M, DRG, GKZZ )

      REAL UB(N,M), TH(N,M) , DRG(N,M), Y(N), ZPZ(M), P(100)
      REAL GKZZ(N,M)
      REAL U(100) , B(100) , BV(100,100)  , D(100) , PHI(100)
      REAL F(100) , Z(100)

      INCLUDE 'GWRAY.INC'

      COMMON /CONST2_1/ PI,R,H,CP,TOMEG,A,XKAPPA,IPRINT,NDBUG

ccelf      include 'test_value.inc'
ccelf      logical nantest1, nantest2
    
      DLAT2 = ( Y(2) - Y(1) ) / 2.00

C      WRITE(6,'("---- WILL DO ",I4," GW RAY CALCULATIONS " )') NRAYS

      GRAV  = 980.   

      DO K=1,M
         Z(K) = ZPZ(K) *1.E5
         P(K) = EXP( -ZPZ(K) / 7.00 )
      END DO

      DO K=1,100
      DO J=1,100
         BV(J,K) = 0.00
      END DO
      END DO

        
 
      DO K=2,M-1
      DO J=1,N
         BV(J,K) = ( TH(J,K+1) - TH(J,K-1) ) /
     >             (  Z(K+1) -  Z(K-1) )

      END DO
      END DO

      DO J=1,N
         BV(J,1) = BV(J,2)
         BV(J,M) = BV(J,M-1) 
      END DO

      DO K=1,M
      DO J=1,N
         BV(J,K)   = GRAV * BV(J,K) /  TH(J,K) 
         DRG(J,K)  = 0.00
         GKZZ(J,K) = 0.00
      END DO
      END DO

      DO K=1,M
      DO J=1,N
         IF ( BV(J,K) .GE. 0 ) BV(J,K)  = SQRT( BV(J,K) )
         IF ( BV(J,K) .LT. 0 ) BV(J,K)  = 0.00 

C         IF (BV(J,K).EQ.0.0) THEN
C            PRINT*, "BV ZERO in XGWDRG:"
C            PRINT*, "   J = ", J, ", K = ", K
C            PRINT*, "   TH(J,K+1) = ", TH(J,K+1)
C            PRINT*, "   TH(J,K-1) = ", TH(J,K-1)
C            PRINT*, "      Z(K+1) = ", Z(K+1)
C            PRINT*, "      Z(K-1) = ", Z(K-1)
C            STOP
C            ENDIF
      END DO
      END DO


      DO  NN = 1, NRAYS


          RLATG = RAYINI( 1, NN )
 
          DO J=1,N
             IF ( ABS( RLATG - Y(J) ) .LE. DLAT2 ) JRY = J
          END DO

          A_CIRCUM = A   ! * COS( Y(JRY) * PI / 180. )          

          DO K=1,M 

             U(K) = UB(JRY, K) 
             B(K) = BV(JRY, K)

c              nantest1=test_nan(u(k))
c              nantest2=test_nan(b(k))
c              if(nantest1 .or. nantest2) then
c                PRINT*, "NAN in XGWDRG:"
c                PRINT*, "   K = ", K
c                PRINT*, "   U = ", U(K)
c                PRINT*, "   B = ", B(K)
c                STOP
c                ENDIF

          END DO

          A0 = RAYINI( 3, NN )
          C0 = RAYINI( 4, NN )
          RN = RAYINI( 6, NN )


          CALL GWRAY( U(1), B(1), P(1), D(1), A0 , C0, PHI(1), M )
          CALL GWACC( PHI(1), P(1), Z(1), C0,  F(1), M )

C         FSCALE factor added 10/10/95 by PEM to eliminate problem with
C         excessive gravity-wave induced Kyy at the poles

          FSCALE = 1.0
C          ABSLAT = ABS(Y(JRY))
C          IF (ABSLAT .GE. 85.) THEN
C                FSCALE = 0.0
C             ELSEIF (ABSLAT .GE. 75.) THEN
C                FSCALE = EXP(-1.0*(((ABSLAT - 75.0)/5.0)**2.0))
CC                FSCALE = ((ABSLAT - 90.)/15.)**2.0
CC                FSCALE = ABS((ABSLAT - 90.0)/15.0)
CC                FSCALE = 1 - ((ABSLAT - 75.)/15.)**2.0
C             ELSE
C                FSCALE = 1.0
C             ENDIF

          DO K=1,M
             DRG(JRY,K) = DRG(JRY,K) - ( RN*F(K) / A_CIRCUM )*FSCALE
c              nantest1=test_nan( drg(jry,k) )
c              if(nantest1) then
c                PRINT*, "NAN for DRG in XGWDRG:"
c                PRINT*, "        JRY = ", JRY, ", K = ", K
c                PRINT*, "         RN = ", RN
c                PRINT*, "          F = ", F(K)
c                PRINT*, "   A_CIRCUM = ", A_CIRCUM
c                STOP
c                ENDIF
          END DO

          DO K=1,M
             IF (B(K).GT.0.0) THEN
                CK  = ABS( U(K) - C0 ) / ( B(K)**2 )
             GKZZ(JRY,K) = GKZZ(JRY,K) + (CK*RN*F(K) / A_CIRCUM)*FSCALE
                ENDIF
                
c              nantest1=test_nan( gkzz(jry,k) )
c              if(nantest1) then
c                PRINT*, "NAN for GKZZ in XGWDRG:"
c                PRINT*, "        JRY = ", JRY, ", K = ", K
c                PRINT*, "         RN = ", RN
c                PRINT*, "          F = ", F(K)
c                PRINT*, "   A_CIRCUM = ", A_CIRCUM
c                PRINT*, "         CK = ", CK
c                PRINT*, "          U = ", U(K)
c                PRINT*, "         C0 = ", C0
c                PRINT*, "          B = ", B(K)
c                PRINT*, "     FSCALE = ", FSCALE
c                STOP
c                ENDIF
          END DO

          
CCC          write(6,*) '--- RAY#',nn,rlatg,jry,y(jry),a0,c0,rn

       END DO



       RETURN
       END
