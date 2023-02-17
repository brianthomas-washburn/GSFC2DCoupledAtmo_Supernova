      SUBROUTINE WED_RAD(LL)

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'
      INCLUDE 'COMMONR.INC'  
      INCLUDE 'COMMONW.INC'


      
        WRITE(6,*) '  --- COUPLING JOAN AND XUN (RAD. CODES)  '
       

        CALL GETRAD

        DO K=1,M$
        DO J=1,N$
           DUM1(J,K)=HEAT(J,K)
           DUM2(J,K)=COOL(J,K)
        END DO
        END DO

        CALL RADNLTE(LL)


        DO K=1,M$
        DO J=1,N$

           IF (ZP(K).LE.10) THEN

              HEAT(J,K) = DUM1(J,K)
              COOL(J,K) = DUM2(J,K)

           ENDIF


           IF ((ZP(K).GT.10).AND.(ZP(K).LE.20)) THEN

              WTXZ = ( ZP(K)-10. )/10.
              WTJR = 1.00 - WTXZ
 
              HEAT(J,K) = WTJR*DUM1(J,K) + WTXZ*HEAT(J,K)
              COOL(J,K) = WTJR*DUM2(J,K) + WTXZ*COOL(J,K)

           ENDIF


        END DO
        END DO


        RETURN
        END

