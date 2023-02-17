      SUBROUTINE SMTHPSI

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'

C1234.T

c***************************************** DISABLED 9/27/93
c       ICN  = ISW(24)
c       IREM = ICN
c
c       IPSI = INT( IREM/1000 )
c       IREM = IREM-IPSI*1000
c
c       IMER = INT( IREM/100 )
c       IREM = IREM-IMER*100
c
c       IPAS = INT( IREM/10 )
c       IREM = IREM-IPAS*10
c
c       IUBR = INT( IREM/1 )
c       IREM = IREM-IUBR*1
c
c       WRITE( *,1000 ) IPSI ,IMER ,IPAS ,IUBR
c1000   FORMAT( ' IPSI',i2,'; IMER' ,i2,'; IPAS' ,i2,'; IUBR',i2)
c*************************************************************
 
       IPSI=2    !  IN PLACE OF ABOVE 9/27/93
       IMER=0    !

       IF (IPSI.EQ.1) THEN 
          CALL ZISMTH( PSI, N$, M$)
c          WRITE(6,*) '  ZI- ZI- ZI- SMOOTHED PSI FIELD --- --- ----'
       END IF

       IF (IMER.EQ.1) THEN 
          CALL ZISMTH(  VS, N$, M$)
          CALL ZISMTH(  WS, N$, M$)
c          WRITE(6,*) '  ZI- ZI- ZI- SMOOTHED VSTAR AND WSTAR FIELDS'
       END IF


       IF (IPSI.EQ.2) THEN 
          CALL EQSMTH( PSI, N$, M$)
c          WRITE(6,*) '  EQ- EQ- EQ- SMOOTHED PSI FIELD --- --- ----'
       END IF

       IF (IMER.EQ.2) THEN 
          CALL EQSMTH(  VS, N$, M$)
          CALL EQSMTH(  WS, N$, M$)
c          WRITE(6,*) '  EQ- EQ- EQ- SMOOTHED VSTAR AND WSTAR FIELDS'
       END IF

       RETURN
       END

      SUBROUTINE ZISMTH(AA,N,M)

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'
     
      REAL  AA(N,M)
 
C1234.T
      DO K=1,M
      DO J=2,N-1

           DUM1(J,K)=
     >     AA(J,K)  + 
     >         (ZIKYY(J,K)*DT/(DY*DY) ) *
     >     ( AA(J+1,K) - 2*AA(J,K) + AA(J-1,K) )
      
      END DO      
      END DO

      DO K=1,M

           J=N
           DUM1(J,K)=
     >     AA(J,K)  + 
     >         (ZIKYY(J,K)*DT/(DY*DY) ) *
     >     ( AA(J-1,K) - 2*AA(J,K) + AA(J-1,K) )

           J=1
           DUM1(J,K)=
     >     AA(J,K)  + 
     >         (ZIKYY(J,K)*DT/(DY*DY) ) *
     >     ( AA(J+1,K) - 2*AA(J,K) + AA(J+1,K) )
      
      END DO

      DO K=1,M
      DO J=1,N
         AA(J,K)=DUM1(J,K)
      END DO      
      END DO

     
      DO K=2,M-1
      DO J=1,N

           DUM1(J,K)=
     >     AA(J,K)  + 
     >         (ZIKZZ(J,K)*DT/(DZ*DZ) ) *
     >     ( AA(J,K+1) - 2*AA(J,K) + AA(J,K-1) )
      
      END DO      
      END DO

      DO J=1,N

           K=1
           DUM1(J,K)=
     >     AA(J,K)  + 
     >         (ZIKZZ(J,K)*DT/(DZ*DZ) ) *
     >     ( AA(J,K+1) - 2*AA(J,K) + AA(J,K+1) )

           K=M
           DUM1(J,K)=
     >     AA(J,K)  + 
     >         (ZIKZZ(J,K)*DT/(DZ*DZ) ) *
     >     ( AA(J,K-1) - 2*AA(J,K) + AA(J,K-1) )
      
      END DO      


      DO K=1,M
      DO J=1,N
         AA(J,K)=DUM1(J,K)
      END DO      
      END DO

      RETURN
      END


      SUBROUTINE EQSMTH(AA,N,M)

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'
     
      REAL  AA(N,M)
 
C1234.T

      DO K=1,M
      DO J=1,N
         DUM1(J,K)=AA(J,K)
      END DO      
      END DO

      DO K=2,M-1
      DO J=1,N

         IF ( ABS( YP(J) ) .LT. 15. ) THEN
            DUM1(J,K)=
     >      ( AA(J,K+1) + 2*AA(J,K) + AA(J,K-1) )/4.0
         END IF
      
      END DO      
      END DO


      DO J=1,N

           IF ( ABS( YP(J) ) .LT. 15. ) THEN
           K=1
           DUM1(J,K)=
     >     ( AA(J,K+1) + 2*AA(J,K) + AA(J,K+1) )/4.

           K=M
           DUM1(J,K)=
     >     ( AA(J,K-1) + 2*AA(J,K) + AA(J,K-1) )/4.
           END IF
       
      END DO      


      DO K=1,M
      DO J=1,N
         AA(J,K)=DUM1(J,K)
      END DO      
      END DO

      RETURN
      END








