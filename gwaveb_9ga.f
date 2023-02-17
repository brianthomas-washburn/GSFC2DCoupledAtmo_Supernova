
C***********************************************************
C   GRAVITY WAVE PARAMETERIZATION FROM MRS 
C***********************************************************
C

C  NOTE:  ROUTINE GWAVEB is NEVER CALLED -EF, Oct. 2007
C

      SUBROUTINE GWAVEB
      PARAMETER (NPHAZ=5)
      PARAMETER (NPHAZ3=NPHAZ*2-1)

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'

C      SIMPLE GRAVITY WAVE PARAMETERIZATION
      DIMENSION GWMOMW(NP$,MP$,NPHAZ3) , PHASSP(NPHAZ3),D4(MP$) ,
     &          GWMOMO(NP$,mP$)    ! 

      DIMENSION PHSP(NPHAZ)
c      DATA PHSP/0.,5.,10.,15.,20.,25.,30.,35./
c      DATA PHSP/0.,5.,15.,25.,35.,45./
      DATA PHSP/0.,10.,20.,30.,40./
     1,CON1/0.0/,ZB0max/70./,zb0min/55./,ZSTART/15./,clevel/16./

C      THE GRAVITY WAVE DECELERATION IS GIVEN BY K/NH *(UB-PHSP)^3
C      ZSTART IS THE STARTING LEVEL FOR CRITICAL LEVEL CHECK
C      ZB0MAX IS THE MAXIMUM WAVE BREAKING LEVEL IN KM
C      ZB0MIN IS THE MINIMUM WAVE BREAKING LEVEL IN KM
C      THE HORIZONTAL WAVENUMBER IS XKH IN CM-1
C      THE BREAKING LEVEL IS GIVEN BY ZB0+CON1*ABS(PHSP(M/S))
C      PHSP IS THE PHASE SPEEDS CONSIDERED

C************************************
C    JTB-  CHANGED GRID FROM 
C    FROM MRS TO RBR  (10/21/92)
C************************************



      ITG=IT3
      NPHAZ2=NPHAZ
c	nphaz2=3
	XKH=2.*PI/(1.0E10)

C     GET STARTING LEVEL INDEX, KL 
C     ZSTART IS THE STARTING LEVEL IN KM

      DO 50 K=1,MP$
      KL=K
      IF (ZSTART.LT.ZPP(K)) GO TO 55
50    CONTINUE
55    CONTINUE      

      xns=sqrt(xn2(1))  ! BRUNT_VAISAALAAAA FReq 



      do 1 j=1,np$

c     loop over latitudes

C	FIRST ZERO THE FLUX AND DETERMINE THE Breaking Level
C       FOR THE C=0 WAVE

      do 18 k=1,mp$
      gwmomo(j,k)=DRAGX(j,k)  !STORE OFF OLD GWAVE DRAG
      DRAGX(J,K)=0.
      xn=sqrt(xn2(k))
      strx=sqrt(abs(ubx(j,k)*ubx(j,1)/(xn*xns)))
      strx=strx*clevel/gvint(j)
      ztest=2.*H*1.0e-5*alog(strx)
      if (ztest.lt.zpp(k)) goto 555
      kb=k
555   DO 18 II=1,NPHAZ3
      GWMOMW(J,K,II)=0.
18    CONTINUE

c     zb0 is the breaking level for the c=0 wave

      zb0=zpp(kb)
      zb0=amin1(zb0,zb0max)     
      zb0=amax1(zb0,zb0min)

      I3=0
      DO 2 IP=1,NPHAZ2 !number of phase speeds
      I77=2
      IF (IP.EQ.1) I77=1

C     ZB is the breaking level
      
      ZB=ZB0+CON1*ABS(PHSP(IP))
      
C	LOOP OVER THE PHASE SPEEDS

      DO 33 II=1,I77
      I3=I3+1
C     SGN IS THE  SIGN OF THE PHASE SPEED
C     C=0 IS TREATED AS A SPECIAL CASE
      SGN=1-(II-1)*2
      PSP=PHSP(IP)*SGN
      PHASSP(I3)=PSP*1.
      K=KL-1
      UT2=((UBX(J,K))-PSP*100.)

      DO 3 K=KL,MP$
      UT1=((UBX(J,K))-PSP*100.)
C                               CRITICAL LEVEL CHECK
      IF(UT1*UT2.gt.0.) GOTO 331 
      go to 33

331   IF (ZPP(K).GT.ZB) GOTO 10 
C                               BREAKING LEVEL CHECK
       UT2=UT1
         GOTO 3
10      CONTINUE

C	 WAVE IS BREAKING

        UK=XKH/(H*SQRT(XN2(K)))

C	XLR IS INTERMITTENCY FACTOR (XLR<=1)

c        rss=amax1(1.,abs(psp))
c	XLR=.2/(rss)

        XLR=1.0
        PRN=.1
	bf=-xlr*uk*ut1**3
        GWMOMW(J,K,I3)=BF
        DRAGX(J,K)=DRAGX(J,K)+bf
	KZZ(J,K)=KZZ(J,K)+PRN*XLR*UK*UT1**4/XN2(K)
3      CONTINUE ! height loop
33     CONTINUE ! phase speed sign loop
2      CONTINUE ! phase speed loop


C      SMOOTH THE FIELDS Vertically

       DO 30 III=1,5
       DO 20 K=2,M$
       D4(K)=0.5*KZZ(J,K)+0.25*(KZZ(J,K+1)+KZZ(J,K-1))
20     D3(K)=0.5*DRAGX(J,K)+0.25*(DRAGX(J,K+1)+DRAGX(J,K-1))

       DO 21 K=2,M$
       KZZ(J,K)=D4(K)
21     DRAGX(J,K)=D3(K)

30     CONTINUE       

1      CONTINUE

c      smooth horizontally with more smoothing at higher resolution

       nsm=1
       if(n$.gt.20) nsm=6

       DO 71 III=1,nsm
       DO 224 K=1,MP$
       DO 221 J=1,NP$
       IF(J.EQ.1) D3(J)=0.5*DRAGX(J,K)+0.25*(DRAGX(J+1,K))
       IF(J.EQ.N$) D3(J)=0.5*DRAGX(J,K)+0.25*(DRAGX(J-1,K))

       IF(J.GT.1.AND.J.LT.N$) D3(J)=0.5*DRAGX(J,K)+
     1                        0.25*(DRAGX(J+1,K)+DRAGX(J-1,K))
221    CONTINUE
       DO 222 J=1,NP$
222    DRAGX(J,K)=D3(J)
224    CONTINUE
71     CONTINUE

c      smooth gravity wave stress in time

cc       do 75 k=1,m$
cc       do 75 j=1,n$
cc       DRAGX(j,k)=0.5*(DRAGX(j,k)+DRAGX(j,k))
cc75     continue


       DO 355 II=1,NPHAZ3
       nn=400+ifix(phassp(ii))
cc       CALL OUTA(3,NSTEP,NDAY,NYEAR,YP,ZP,GWMOMW(1,1),N$,M$,0) 
355    CONTINUE

cc       CALL OUTA(21,NSTEP,NDAY,NYEAR,YPP,ZPP,DRAGX(1,1),Np$,Mp$,0) 

      RETURN
      END            


       SUBROUTINE JMTNWV( MTNDRG )
 
C
C  SIMPLE CUBIC LAW MOUNTAIN WAVE DRAG 
c  AND CONSISTENT KZZ - ADDED 11/28/94 JTB
C

	INCLUDE 'PARAM.INC'
	INCLUDE 'COMMONC.INC'
	INCLUDE 'COMMOND.INC'
	INCLUDE 'COMMONW.INC'
	INCLUDE 'COMMONT.INC'
 
        REAL MTNDRG(NP$,MP$)

        CUBDRG=1./( 1.5*DAYL*6000.*6000.)   ! DRAG TIME FOR
C                                           ! U=60M/S  --> 1.5 DAYS


        DO K=1,MP$
        DO J=1,NP$
           MTNDRG(J,K) =
     2          CUBDRG*UBX(J,K)*UBX(J,K)
     3                                  *UBX(J,K)
        END DO
        END DO

        DO K=1,MP$      !  ZERO OUT MOUNTAIN WAVE DRAG BELOW MWBRK, 
        DO J=1,NP$      !  AND RAMP-UP TO FULL OVER MWRMP KM
           IF(ZPP(K) .LT. MWBRK(J) ) 
     2              MTNDRG(J,K) = 0.000
           IF( (ZPP(K) .GE. MWBRK(J) ) .AND.
     1         (ZPP(K) .LT. MWBRK(J)+MWRMP ) )  
     2              MTNDRG(J,K) = (1./MWRMP)*(ZPP(K)-MWBRK(J))*
     3                            MTNDRG(J,K) 
        END DO
        END DO

        

        DO K=1,MP$
        DO J=1,NP$
           GWKZZ(J,K) = 
     >          ( UBX(J,K) / XN2ST(K) ) 
     >          *   MTNDRG(J,K)
        END DO
        END DO



        RETURN
        END


