
      SUBROUTINE ELLIPQG

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'
      
ccelf      include 'test_value.inc'
ccelf      logical nantest

      COMMON /SLAKO/ AY(N$,M$),CY(N$,M$),B(N$,M$)
     1,AZ(N$,M$),CZ(N$,M$)
     2,F(N$,M$)

      DIMENSION COYY(N$,M$),COZZ(N$,M$),COZ(N$,M$)
     2, COY(N$,M$)     
      

      DO K=1,M$
      DO J=1,N$
         DUM1X(J,K)= 0.5*( XC00(J,K+1,1)-XC00(J,K,1) )/DZ
     >             + 0.5*( XC00(J+1,K+1,1)-XC00(J+1,K,1) )/DZ
      END DO
      END DO 

       
      DO K=1,M$
      DO J=1,N$ 
          DUM2X(J,K)=CF(J)
      END DO
      END DO

      DO K=1,M$
      DO J=1,N$ 
          COZZ(J,K) = CF(J) *DUM2X(J,K)
          COZ (J,K) =-CF(J) *DUM2X(J,K)/H
          COYY(J,K) = RSTAR(K)*DUM1X(J,K)
          COY(J,K) =  RSTAR(K)*DUM1X(J,K)*S(J)/( C(J)*A )
      END DO
      END DO


C -- NON QG TERMS APPENDED TO FORCING -- -- --
C
C                             ** -W*D(U,T)/DZ
      DO K=1,M$
      DO J=1,N$
         DUM1X(J,K)= 0.5*( THX(J,K+1)-THX(J,K) )/DZ
     >             + 0.5*( THX(J+1,K+1)-THX(J+1,K) )/DZ
         DUM2X(J,K)= 0.5*( UBX(J,K+1)-UBX(J,K) )/DZ
     >             + 0.5*( UBX(J+1,K+1)-UBX(J+1,K) )/DZ
         DUM1(J,K)= GHEATX(J,K) 
     >           -WS(J,K)*DUM1X(J,K)
         DUM2(J,K)= GWMOM(J,K) 
     >           -WS(J,K)*DUM2X(J,K)

c       nantest=test_nan( dum1(j,k) )
c       if(nantest) then
c            PRINT*, "DUM1 NAN in ELLIPQG:"
c            PRINT*, "       J = ", J
c            PRINT*, "       K = ", K
c            PRINT*, "  GHEATX = ", GHEATX(J, K)
c            PRINT*, "      WS = ", WS(J, K)
c            PRINT*, "   DUM1X = ", DUM1X(J, K)
c            PRINT*, "   THX+_ = ", THX(J+1, K)
c            PRINT*, "   THX__ = ", THX(J, K)
c            PRINT*, "   THX++ = ", THX(J+1, K+1)
c            PRINT*, "   THX_+ = ", THX(J, K+1)
c            STOP
c            ENDIF

      END DO
      END DO 

C
C                             ** -V*D(U,T)/DY

      DO K=1,M$
      DO J=1,N$
         DUM1X(J,K)= 0.5*( THX(J+1,K)-THX(J,K) )/DY
     >             + 0.5*( THX(J+1,K+1)-THX(J,K+1) )/DY
         DUM2X(J,K)= 0.5*( UBX(J+1,K)-UBX(J,K) )/DY
     >             + 0.5*( UBX(J+1,K+1)-UBX(J,K+1) )/DY
         DUM1(J,K)= DUM1(J,K) 
     >           -VS(J,K)*DUM1X(J,K)
         DUM2(J,K)= DUM2(J,K) 
     >           -VS(J,K)*DUM2X(J,K)

c       nantest=test_nan( dum1(j,k) )
c       if(nantest) then
c            PRINT*, "DUM1 NAN in ELLIPQG:"
c            PRINT*, "       J = ", J
c            PRINT*, "       K = ", K
c            PRINT*, "      VS = ", VS(J, K)
c            PRINT*, "   DUM1X = ", DUM1X(J, K)
c            PRINT*, "   THX+_ = ", THX(J+1, K)
c            PRINT*, "   THX__ = ", THX(J, K)
c            PRINT*, "   THX++ = ", THX(J+1, K+1)
c            PRINT*, "   THX_+ = ", THX(J, K+1)
c            STOP
c            ENDIF

      END DO
      END DO 
C                   
C                             ** +V*U*tan(LAT)
C                                     ! TANGENT TERM
      DO K=1,M$
      DO J=1,N$
         DUM2X(J,K)= 0.25*( UBX(J+1,K) + UBX(J,K) 
     >                    + UBX(J+1,K+1)+UBX(J,K+1) )
         DUM2(J,K)= DUM2(J,K) 
     >           +VS(J,K)*DUM2X(J,K)*S(J)/(A*C(J))
      END DO
      END DO 


C -- -- -- -- -- 


cjer looks like DUM3 is d(DUM1)/d(DY)
      CALL DERV4B(2,DUM1,DUM3,DY,N$,M$,1)   
      CALL DERV4B(2,DUM2,DUM4,DZ,N$,M$,2)    

      DO K=1,M$
         DUM3(1,K)  = ( DUM1(2,K)-DUM1(1,K))*DY1
         DUM3(N$,K) = ( DUM1(N$,K)-DUM1(N1$,K))*DY1
      END DO


      DO K=1,M$
      DO J=1,N$

         F(J,K)= C(J)*
     >          (  RSTAR(K)*DUM3(J,K) +
     >              CF(J)*DUM4(J,K) )
     
c       nantest=test_nan( f(j,k) )
c       if(nantest) then
c            PRINT*, "F(J,K) NAN in ELLIPQG:"
c            PRINT*, "   J = ", J
c            PRINT*, "   K = ", K
c            PRINT*, "   RSTAR = ", RSTAR(K)
c            PRINT*, "   DUM3 = ", DUM3(J, K)
c            PRINT*, "   CF = ", CF(J)
c            PRINT*, "   DUM4 = ", DUM4(J, K)
c            STOP
c            ENDIF

         AZ(J,K)= COZZ(J,K)*DZX + COZ(J,K)*DZ1
         CZ(J,K)= COZZ(J,K)*DZX - COZ(J,K)*DZ1

         AY(J,K)= COYY(J,K)*DYX + COY(J,K)*DY1
         CY(J,K)= COYY(J,K)*DYX - COY(J,K)*DY1
  
         B(J,K)= -2*( COZZ(J,K)*DZX + COYY(J,K)*DYX )


      END DO
      END DO


      DO K=1,M$
      DO J=1,N$
         DUM1(J,K)= C(J)*
     >              RSTAR(K)*DUM3(J,K) 
         DUM2(J,K)= C(J)*
     >              CF(J)*DUM4(J,K) 
      END DO
      END DO




      RETURN
      END 


      SUBROUTINE XSTRM

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'

      COMMON /SLAKO/ AY(N$,M$),CY(N$,M$),B(N$,M$)
     1,AZ(N$,M$),CZ(N$,M$)
     2,F(N$,M$)

      NPASS=3

      DO K=1,M$
      DO J=1,N$
cjer heat and cool here are in pot temp      
         GHEATX(J,K)=HEAT(J,K)-COOL(J,K) 
     >     + ISW(21) *
     >               ( XMIX(J,K,1)  + XMIX(J+1,K,1)
     >               + XMIX(J,K+1,1)+ XMIX(J+1,K+1,1) )/4.

         DUM1(J,K)=GHEATX(J,K)/TTOTH(K)
      END DO
      END DO 


      CALL OUTA(4,NSTEP,NDAY,NYEAR,YP,ZP,DUM1,N$,M$,0)


      DO NN=1,NPASS

      CALL ELLIPQG

      CALL PRMINMAX(B ,'-B-XSTRM--',N$, M$,NDBUG) 
      CALL PRMINMAX(F ,'-F-XSTRM--',N$, M$,NDBUG) 
      CALL PRMINMAX(CZ,'CZ-XSTRM--',N$, M$,NDBUG) 
      CALL PRMINMAX(AZ,'AZ-XSTRM--',N$, M$,NDBUG) 
      CALL PRMINMAX(CY,'CY-XSTRM--',N$, M$,NDBUG) 
      CALL PRMINMAX(AY,'AY-XSTRM--',N$, M$,NDBUG) 

  
      CALL JSOR2(AY,AZ,B,CY,CZ,F,PSI,N$,M$,ERROR,ERRMAX,ERRET,NSTEP)
      CALL PRMINMAX(PSI,'PSI-XSTRM-',N$, M$,NDBUG) 

      CALL SMTHPSI

      CALL JGETVW(PSI,VS,WS,N$,M$)

      CALL PRMINMAX(Vrbr,'VS-XSTRM--',N$, M$,NDBUG) 
      CALL PRMINMAX(Wrbr,'WS-XSTRM--',N$, M$,NDBUG) 


c           write(*,*) ' Ell. pass - -  ', nn

      END DO


      RETURN
      END
      


      
      SUBROUTINE ALGWS

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'
      DIMENSION ZIA(NP$,MP$)

C
C  XC(*,*,1) -- POTENTIAL TEMP.
C  XC(*,*,2) -- ANGULAR MOM.
C
C-----------------------------------------------
C SOLVES ALGEBRAIC EXPRESSION FOR WSTAR:
C
C    { 1 - (Mz/My)(THy/THz) }WS =
C          (Q - THt)/THz - (X - Ut)(THy/THz)/My
C
C-----------------------------------------------
c   




C--------------------------------------------------------------
     
      CALL DERV4B(1,XC(1,1,1),DUM1X,DY,NP$,MP$,1)  ! d(TH)/dy
      CALL DERV4B(1,XC(1,1,2),DUM2X,DY,NP$,MP$,1)  ! d(M)/dy

      DO K=1,MP$

         DUM1X(1, K) =(XC(2,K,1)-XC(1,K,1))/(2*DY)
         DUM1X(NP$,K)=(XC(NP$,K,1)-XC(N$,K,1))/(2*DY)

         DUM2X(1, K) =(XC(2,K,2))/(2*DY)
         DUM2X(NP$,K)=(-XC(N$,K,2))/(2*DY)

      END DO
 
c>>>       
C            DUM1X == d(TH)/dy
c
C            DUM2X ==  d(M)/dy
c
c
c               BC's: 
c                       d(TH)/dy = 0 at sides 
c              
c                          M = 0     at sides              
c
c>>>




C----------------------------------------------------------------
     
      CALL DERV4B(1,XC(1,1,1),DUM3X,DZ,NP$,MP$,2)  ! d(TH)/dz
      CALL DERV4B(1,XC(1,1,2),DUM4X,DZ,NP$,MP$,2)  ! d(M)/dz

      DO J=1,NP$

         DUM4X(J, 1) =( XC(J,2,2) )/(1.25*DZ)
         DUM4X(J,MP$)=(XC(J,MP$,2)-XC(J,M$,2) )/(2.*DZ)

         DUM3X(J, 1) =(XC(J,2,1)-XC(J,1,1) )/(DZ)
         DUM3X(J,MP$)=(XC(J,MP$,1)-XC(J,M$,1) )/(DZ)
 
      END DO

c>>>       
C            DUM3X == d(TH)/dZ
c
C            DUM4X ==  d(M)/dZ
c
c
c               BC's: 
c                       d(TH)/dz(0+) = d(TH)/dz(0-) at top & bottom
c              
c                          M = 0     at bottom  
c                       d(M)/dz = 0  at top            
c
c>>>



C--------------- NOW CALCULATE FORCING TERMS -------------------


   
      CALL MRS2RBR( HEAT , HEATX )
      CALL MRS2RBR( GWMOM, DRAGX )


      DO K=1,MP$
      DO J=1,NP$
         HEATX(J,K)=
     >      (  HEATX(J,K)  
     >        - ( XC(J,K,1)-XCT(J,K,1) )/DT 
     >           ) / DUM3X(J,K)
      END DO
      END DO



C                         ----   +_ 10 DEG FROM EQ. 
C                         ----   USE W ~ (Q-THt)/THz


      DO K=1,MP$
      DO J=1,NP$

         IF ( (YPP(J) .GT. -10.) .AND. (YPP(J) .LT. 10.)) THEN
            DRAGX(J,K)=0.
            ZIA(J,K)=1.0
         ENDIF

         IF ( (YPP(J) .LE. -10.) .OR. (YPP(J) .GE. 10.)) THEN
            DRAGX(J,K)=
     >          ( DRAGX(J,K) - 
     >          ( XC(J,K,2)-XCT(J,K,2) )/DT )*
     >          ( DUM1X(J,K) / DUM3X(J,K) ) /
     >            DUM2X(J,K)
            ZIA(J,K)=
     >         1. - 
     >         ( DUM1X(J,K)/DUM3X(J,K) ) *
     >         ( DUM4X(J,K)/DUM2X(J,K) )
          ENDIF
     
      END DO
      END DO

      DO K=1,MP$
      DO J=1,NP$
         WRBR(J,K)=
     >       ( HEATX(J,K) - DRAGX(J,K) )/ ZIA(J,K)
      END DO
      END DO


C  Force BC's at top and bottom


      DO J=1,NP$
         WRBR(J,1)=0.
         WRBR(J,MP$)=0.
      END DO



cc      CALL OUTA(7,NSTEP,NDAY,NYEAR,YPP,ZPP,HEATX ,NP$,MP$,0)
cc      CALL OUTA(7,NSTEP,NDAY,NYEAR,YPP,ZPP,DRAGX ,NP$,MP$,0)
cc      CALL OUTA(7,NSTEP,NDAY,NYEAR,YPP,ZPP,WRBR ,NP$,MP$,0)

      CALL PRMINMAX(ZIA , 'ZIA-------',NP$, MP$,NDBUG)   
      CALL PRMINMAX(DUM1X,'DTH/DY----',NP$, MP$,NDBUG)   
      CALL PRMINMAX(DUM2X,'DM/DY-----',NP$, MP$,NDBUG)   
      CALL PRMINMAX(DUM3X,'DTH/DZ----',NP$, MP$,NDBUG)   
      CALL PRMINMAX(DUM4X,'DM/DZ-----',NP$, MP$,NDBUG)  
      CALL PRMINMAX(XCT(1,1,1),'OLD TH----',NP$, MP$,NDBUG)  
      CALL PRMINMAX(XCT(1,1,2),'OLD M-----',NP$, MP$,NDBUG)  
 
      RETURN
      END
