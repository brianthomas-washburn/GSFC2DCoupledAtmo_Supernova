
      SUBROUTINE GWDRAG(II)


      
      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONW.INC'
      INCLUDE 'COMMONT.INC'
      INCLUDE 'GWRAY.INC'

C     THETAG added 3/2/96 by PEM to do JTB's suggested modified call to
C     XGWDRG.

      REAL THETAG(NP$, MP$)


      RR   = 1.
      NRAY = 1
      LUN  = 114

C       write(6,*)' HERE WE ARE '

      IF ( II .EQ. 0 ) THEN         ! READ IN RAY INFO AND CUT OUT
         
         DO WHILE ( RR .GT. -9999. ) 
         READ(LUN,1000) RAY0

         RR = RAY0(1)  
         IF ( RAY0(1) .GT. 9999. ) THEN

            DO J = 1, NP$ 

               RAYINI(1,NRAY) = YPP(J)
               DO L = 2,6 
                  RAYINI(L,NRAY) = RAY0(L) 
                  END DO
               NRAY=NRAY+1

               END DO
         ELSE 

            DO L = 1,6
               RAYINI(L,NRAY) = RAY0(L) 
               END DO
            NRAY=NRAY+1
            END IF
         END DO

         CLOSE(UNIT = LUN)

         NRAYS = NRAY-2

         do n = 1,nrays
            DO L = 1,6
                 RAY0(L) =RAYINI(L,N) 
            end do
c            write(6,1001) ray0
         enddo

         NRAYS = NRAY-2
 
         RETURN

      END IF


C     CALL to XGWDRG modified 3/2/96 by PEM, as per JTB.

C      CALL XGWDRG( UBX , XC(1,1,1) , YPP, ZPP, NP$, MP$, DUM1X, DUM2X )

      DO NDEX = 1, NP$
         THETAG(NDEX, 1) = THG(1)
         THETAG(NDEX, MP$) = THG(M$)
         DO MDEX = 2, (MP$ - 1)
            THETAG(NDEX, MDEX) = 0.5*(THG(MDEX) + THG(MDEX - 1))
            END DO
         END DO

      CALL XGWDRG( UBX, THETAG, YPP, ZPP, NP$, MP$, DUM1X, DUM2X)




 1000 FORMAT(6(F10.0))
 1001 FORMAT("--- STARTING VALUES FOR GW RAYS ", 6(F10.0))

      RETURN
      END
        

