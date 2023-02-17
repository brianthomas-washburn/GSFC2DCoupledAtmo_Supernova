      SUBROUTINE GWRAY( U, B, P, D, A0, C0, PHI, M )
      
ccelf      include 'test_value.inc'
ccelf     logical nantest
 
      REAL U(M) , B(M) , D(M), P(M), PHI(M)

      REAL UW(100)

      IF ( B(2) .LT. 0.001 ) THEN
         DO K = 1,M
            PHI(K) =0.00      
         END DO
         RETURN
      END IF     

      USGN = SIGN( 1.000, (U(2)-C0)  )

      DO K=2,M

         UW(K) = ( U(K) - C0 )*USGN
c      if(uw(k) .le. 0.0) print *,'k,u,c0,usgn,uw=',k,u(k),c0,usgn,uw(k)
         UW(K) = AMAX1( UW(K) , 0.00 )
c	if(uw(k) .le. 0.0) print *,'uw = ',uw(k)

c       nantest=test_nan( uw(k) )
c       if(nantest) then
c            PRINT*, "NAN in GWRAY:"
c            PRINT*, "   K = ", K
c            PRINT*, "   U = ", U(K)
c            PRINT*, "  C0 = ", C0
c            STOP
c            ENDIF

      END DO

      UW(1) = UW(2)

    
      SAT0 = UW(2) / B(2)

      A0   = AMIN1( A0 , SAT0 )

      PHI(2) = P(2)* B(2) * UW(2) * (A0**2) 

      PHI(1) = PHI(2)
      D(1)   = A0
      D(2)   = A0

      DO K=3,M 

         IF ( ABS( U(K) - C0 ) .GT. 0.00 ) THEN

             AK  = SQRT( PHI(K-1) / ( P(K)* B(K) * UW(K) ) )
cjer  I've added this - ok per Julio
             if(phi(k-1).eq.0.0 .and. uw(k).eq.0.0) ak = 0.000
cjer end of fix	      	     
c	     nantest=test_nan(ak)
c	     if(nantest) then
c	       print *,'nan1 in ak in gwray: '
c	       print *,'k,phi,p,b,uw = ',k,phi(k-1),p(k),b(k),uw(k)
c	     endif

         ELSE

             AK = 0.000

         END IF

         IF ( B(K) .GT. 0.001 ) THEN

             SATK = UW(K) / B(K)

         ELSE

             SATK = 0.00

         ENDIF

         AK = AMIN1( AK, SATK )
c	 nantest=test_nan(ak)
c	 if(nantest) then
c	   print *,'nan2 in ak in gwray'
c	   print *,'k,ak,satk = ',k,ak,satk
c	 endif
     
         D(K)   = AK

         PHI(K) = P(K)* B(K) * UW(K) * (AK**2) 

      END DO   
      DO K = 1,M
         PHI(K) = PHI(K) * USGN

c         IF (PHI(K) .NE. PHI(K)) THEN
c         nantest=test_nan( phi(k) )
c         if(nantest) then
c            PRINT*, "NAN in GWRAY:"
c            PRINT*, "   USGN = ", USGN
c            PRINT*, "   P    = ", P(K)
c            PRINT*, "   B    = ", B(K)
c            PRINT*, "   UW   = ", UW(K)
c            PRINT*, "   AK   = ", AK
c            STOP
c            ENDIF

      END DO

      RETURN
      END
      


 
      SUBROUTINE GWACC( PHI, P, Z, C0,  F, M )
      
ccelf      include 'test_value.inc'
ccelf      logical nantest
 
      REAL F(M), Z(M), P(M), PHI(M)

      F(1)=0.00
      F(M)=0.00

      DO K=2,M-1 

         F(K) = ( PHI(K+1) - PHI(K) ) / ( Z(K+1) - Z(K) )

      END DO
  
      k=M
c      F(K) = ( - PHI(K) ) / ( Z(K) - Z(K-1) )
    
      DO K=1,M 

         F(K) = F(K) / P(K)

c         IF (F(K) .NE. F(K)) THEN
c       nantest=test_nan( f(k) )
c       if(nantest) then
c            PRINT*, "NAN in GWACC:"
c            PRINT*, "   K = ", K
c            PRINT*, "   P = ", P(K)
c            PRINT*, "   PHI = ", PHI(K)
c            PRINT*, "   Z = ", Z(K)
c            STOP
c            ENDIF

      END DO

      RETURN
      END
    
       

	
        
 
