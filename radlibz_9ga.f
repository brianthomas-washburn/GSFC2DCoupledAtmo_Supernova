
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C######################### RADLIBZ.FOR #########################
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C --------------- JANUARY, 1992  -----------------
C --------------- NEW C-K METHOD -----------------

      SUBROUTINE INPUTZ

      PARAMETER (IMS=40)

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'

      COMMON /CURTZ2/ Z0(IMS),P0(IMS),T0(IMS),GDP(IMS),CURT0(IMS,IMS)
     &,A10(IMS,IMS),A01(IMS,IMS),A20(IMS,IMS),A02(IMS,IMS),A11(IMS,IMS)
      COMMON /CURTZ3/ R0(IMS),ZH00(IMS),ZGDP(IMS)
     &  ,ZCURT0(IMS,IMS),ZA01(IMS,IMS),ZA02(IMS,IMS)

      COMMON /ZHU2/ ZM(IMS),TZX(IMS),PRE(IMS),O2(IMS),O3(IMS),
     _              BO2(IMS),BO3(IMS),OXY(IMS),HO23(IMS),
     _              QO3(IMS),QCO2(IMS),QTZX(IMS),QNET(N$,M$)

C --- data files for z in meter, pressure in Pascal, [O3] in m^-3,
C --- [O] in m^-3.

      OPEN(514,FILE='T2DRAD.IN1',STATUS='OLD')
      DO 70 I=1,IMS
  70  READ(514,14) ISER,ZM(I),TZX(I),PRE(I),XCO2,O3(I),OXY(I)
  14  FORMAT(I3,E12.4,F10.2,4E13.4)
      CLOSE(UNIT=514)

C ---  cooling matrices:

      OPEN(514,FILE='CURT2.DAT',STATUS='OLD')
      READ(514,55) (Z0(K),T0(K),P0(K),GDP(K),K=1,IMS)
      READ(514,220) ((CURT0(L1,L2),L1=1,IMS),L2=1,IMS)
      READ(514,220) ((A10(L1,L2),L1=1,IMS),L2=1,IMS)
      READ(514,220) ((A01(L1,L2),L1=1,IMS),L2=1,IMS)
      READ(514,220) ((A20(L1,L2),L1=1,IMS),L2=1,IMS)
      READ(514,220) ((A02(L1,L2),L1=1,IMS),L2=1,IMS)
      READ(514,220) ((A11(L1,L2),L1=1,IMS),L2=1,IMS)
      CLOSE(UNIT=514)
      OPEN(514,FILE='CURT3.DAT',STATUS='OLD')
  55  FORMAT(7E14.5)
 220  FORMAT(6E14.5)
      READ(514,55) (Z0(K),T0(K),P0(K),R0(K),ZGDP(K),ZH00(K),K=1,IMS)
      READ(514,220) ((ZCURT0(L1,L2),L1=1,IMS),L2=1,IMS)
      READ(514,220) ((ZA01(L1,L2),L1=1,IMS),L2=1,IMS)
      READ(514,220) ((ZA02(L1,L2),L1=1,IMS),L2=1,IMS)
      CLOSE(UNIT=514)

c      do k=1,ims
c      write(6,55) Z0(K),T0(K),P0(K),R0(K),ZGDP(K),ZH00(K)
c      end do

      RETURN
      END


      SUBROUTINE RADNLTE(LL)

      INCLUDE 'PARAM.INC'
      PARAMETER (IMS=40,IM1=15,IDM=5,IM2=IMS-IM1+IDM)
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONR.INC'
      INCLUDE 'COMTEMP.INC'
      
C      INCLUDE 'COMMONJT.INC'

      COMMON /ZHU2/ ZM(IMS),TZX(IMS),PRE(IMS),O2(IMS),O3(IMS),
     _              BO2(IMS),BO3(IMS),OXY(IMS),HO23(IMS),
     _              QO3(IMS),QCO2(IMS),QTZX(IMS),QNET(N$,M$)

      COMMON /NLTE/ NCALL, NLTE$
      CHARACTER*3 NLTE$



cxz add begin =====================================
	INTEGER IDAYS(13)
	DATA IDAYS/1,32,60,91,121,152,182,213,244,274,305,335,365/
cxz add end =====================================


      NCALL=0
      READ(102,10) NLTE$
10    FORMAT( A3 )
      REWIND( UNIT=102 )
      WRITE(*,*) 'WHICH CO2-O QUENCHING -- ',NLTE$


      IF(LL.EQ.1) CALL INPUTZ


      CALL NETHTZ
      
ccc      IF( (NSTEP .GT. 12).AND.(ISW(28).EQ.1)) go to 333  ! SKIP TROPO H&C
      IF(ISW(28).EQ.1) go to 333  ! SKIP TROPO H&C





C**********************************************************
C   BEGIN ADDING IN KLUGE-Y TROPO. HEATING AND COOLING    *
C**********************************************************
      
cxz add begin ================= from GETRAD2.FOR ==========

C     ADD TROPOSPHERIC HEAT SOURCE FROM CUNNOLD ET AL.
980      ZP00=1.3*H

      DX=real(IDAYX)/(IDAYS(MONTH+1)-IDAYS(MONTH))
      DO 171 J=1,N$
      DO 171 K=1,M$
  171 DUM1(J,K)=TFROP(J,K,MONTH) + 
     1           (TFROP(J,K,MONTH+1)-TFROP(J,K,MONTH))*DX


c               ------ DUM1 NOW CONTAINS THE INTERPOLATED
C               ------ (FROM MONTHLY MEANS) NMC CAC TEMPERATURE
C               ------ FOR THE CURRENT MODEL DAY

  
      ENHFAC=ISW(39)/100.
      ENHLAT=ISW(40)*1.00
      ENHHGT=ISW(41)*1.00

      DO K=1,M$
      DO J=1,N$

         EQENH=(1.-ABS(YP(J))/ENHLAT)*ENHFAC
         EQENH=AMAX1(EQENH,0.000 )

         EQENH=(1. -   ZP(K) /ENHHGT)*EQENH
         EQENH=AMAX1(EQENH,0.000 )

         DUM2(J,K)=(1.00+EQENH)*DUM1(J,K)
      END DO
      END DO

c               ------ DUM2 NOW CONTAINS EQUATORIALLY ENHANCED
C               ------ NMC CAC TEMPERATURES FOR THE CURRENT 
C               ------ MODEL DAY




      HCFAC=ISW(38)/100.

      F1=1.0*HCFAC
      F2=1.0*HCFAC

      IF(ISW(28).EQ.0) THEN
c      WRITE(6,*) '  USING CUNNOLDIAN TROPO IN XZ RAD  (FACTOR)=',
c     2              ISW(38)

      DO 210 K=1,M$
      IF (ZZ(K).GT.ZP00+1.5*H) GO TO 210
      WGT=1.0
      IF (ZZ(K).GT.ZP00) WGT=EXP( -(ZZ(K)-ZP00)/(2.*H) )
      VNC2=VNC(K)
      IF (WGT.EQ.1.) VNC2=3.0E-6
      DO 200 J=1,N$
      HEAT(J,K)=WGT*VNC2*(DUM2(J,K)-TG(K))*F2 + 
     1        (1.-WGT)*HEAT(J,K)
      COOL(J,K)=(1.-WGT)*COOL(J,K) +
     1        F1*WGT*VNC2*(T(J,K)-THG(K)/TTOTH(K))
200   CONTINUE
210   CONTINUE
      END IF

      IF(ISW(28).EQ.2) THEN
      WRITE(6,*) '  USING JTB-CUNNOLDIAN TROPO IN XZ RAD  (FACTOR)=',
     2              ISW(38)

      DO K=1,M$
      IF (ZZ(K).LE.ZP00+1.5*H) THEN
      WGT=1.0
      IF (ZZ(K).GT.ZP00) WGT=EXP( -(ZZ(K)-ZP00)/(2.*H) )
      VNC2=VNC(K)
      IF (WGT.EQ.1.) VNC2=3.0E-6
      DO J=1,N$
      HEAT(J,K)=WGT*VNC2*( DUM2(J,K)-TG(K) )*F2 + 
     1        (1.-WGT)*HEAT(J,K)
      COOL(J,K)=(1.-WGT)*COOL(J,K) +
     1        F1*WGT*VNC2*( T(J,K)-DUM1(J,K) )
      END DO
      END IF
      END DO
      END IF

      IF(ISW(28).EQ.3) THEN
      WRITE(6,*) '  USING NEWTONIAN COOL TO NMC TEMPS (FACTOR) =',
     2              ISW(38)

      DO K=1,M$
      IF (ZZ(K).LE.ZP00+1.5*H) THEN
      WGT=1.0
      IF (ZZ(K).GT.ZP00) WGT=EXP( -(ZZ(K)-ZP00)/(2.*H) )
      VNC2=VNC(K)
      IF (WGT.EQ.1.) VNC2=3.0E-6
      DO J=1,N$
      COOL(J,K)=(1.-WGT)*COOL(J,K) +
     1        F1*WGT*VNC2*( T(J,K)-DUM2(J,K) )
      END DO
      END IF
      END DO
      END IF



c      CALL OUTA(51,NSTEP,NDAY,NYEAR,YP,ZP,DUM1,N$,M$,0)
c
c********************************************************

C********
C      WRITE (6,40)
C40    FORMAT (2x,'TOTAL HEATING: PRODUC ',/)
C      DO 42 I=1,M$
C      IZ = M$ - I + 1
C      ZKM=ZM(IZ)/1000.0
C      WRITE (6,41) I,ZKM,(HEAT(J,I)*DAYL,J=1,N$)
C41    FORMAT (2X,I4,2X,F8.2,2X,18F6.2)
C42    CONTINUE

C********

C      WRITE (6,50)
C50    FORMAT (2x,'TOTAL COOLING: ZHU ',/)
C      DO 52 I=1,M$
C      IZ = M$ - I + 1
C      ZKM=ZM(IZ)/1000.0
C      WRITE (6,51) I,ZKM,(COOL(J,I)*DAYL,J=1,N$)
51    FORMAT (2X,I4,2X,F8.2,2X,18F6.2)
C52    CONTINUE
C********

c      WRITE (6,95)
c95    FORMAT (2x,'NET RADIATIVE DRIVE: COOLZ ',/)
c      DO 97 I=1,M$
c      IZ = M$ - I + 1
c      ZKM=ZM(IZ)/1000.0
c      WRITE (6,51) I,ZKM,(QNET(J,I)*DAYL,J=1,N$)
c97    CONTINUE

C********

C      CALL OUTA(1,NSTEP,NDAY,NYEAR,YP,ZP,HEAT,N$,M$,0)
C      CALL OUTA(2,NSTEP,NDAY,NYEAR,YP,ZP,COOL,N$,M$,0)

C      CALL GETGLM(HEAT,GHET,C,N$,M$)
C      CALL GETGLM(COOL,GCOL,C,N$,M$)
 
C      CALL OUTA(1,NSTEP,NDAY,NYEAR,YP,ZP,HEAT,N$,M$)
C      CALL OUTA(2,NSTEP,NDAY,NYEAR,YP,ZP,COOL,N$,M$)

 333  continue

      DO 2 K=1,M$
      DO 2 J=1,N$
      HEAT(J,K)=HEAT(J,K)*TTOTH(K)
      COOL(J,K)=COOL(J,K)*TTOTH(K)
2     CONTINUE

      RETURN
      END



      SUBROUTINE NETHTZ
      PARAMETER (IMS=40,IM1=15,IDM=5,IM2=IMS-IM1+IDM)
      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'

      COMMON /ZHU2/ ZM(IMS),TZX(IMS),PRE(IMS),O2(IMS),O3(IMS),
     _              BO2(IMS),BO3(IMS),OXY(IMS),HO23(IMS),
     _              QO3(IMS),QCO2(IMS),QTZX(IMS),QNET(N$,M$)

      DIMENSION SMT(M$)

      DATA RAIR0/287.0D0/,RO2/259.2D0/,RO3/172.8D0/

      BC = 1.3804E-23                  !Boltzman's constant in MKS
      O2MIX = 0.20


C -- loop over latitude:
      DO 5 J=1,N$
      DEC = DECL
      XPHI = YP(J)

C -- loop over altitude:
      DO 3 K=1,M$
      IK = M$ - K + 1
      TOTN = (PRE(K)/(BC*T(J,IK)))        !density in mks
      O3MIX = O3M(J,IK)/1.0E+6            !O3 MIXING RATIO
      O3(K) = O3MIX*TOTN                  !O3 density
      O2(K) = O2MIX*TOTN                  !O2 density
      TZX(K) = T(J,IK)                    !temperature
3     CONTINUE

C ----------------------------------------------------
C -- SMOOTH TEMPERATURE FIELDS IN Z:

C      MM0 = 2
C      MM1 = M$-1

C      DO 400 I=1,4     !NUMBER OF SMOOTH PASSES
C      DO 200 K=MM0,MM1
C200   SMT(K)=0.25*(TZX(K-1)+TZX(K+1))+0.5*TZX(K)
C      DO 300 K=MM0,MM1
C      IK=M$-K+1
C      T(J,K)=SMT(IK)
C300   TZX(K)=SMT(K)
C400   CONTINUE
C--  TOP TEMPERATURE CORRECTION:
C      DTDZ = (T(J,M$-2)-T(J,M$-1))/(ZM(M$-2)-ZM(M$-1))
C      T(J,M$) = T(J,M$-1)+DTDZ*(ZM(1)-ZM(2))
C      TZX(M$)=T(J,M$)
C ---------------------------------------------------------

C   Very crude estimate of the column density of O2 and O3.

      BO2(1)=O2(1)*5.2E3
      BO3(1)=O3(1)*3.1E3
      DO 20 I=2,IMS
      BO2(I)=BO2(I-1)+0.5*(O2(I)+O2(I-1))*(ZM(I-1)-ZM(I))
      BO3(I)=BO3(I-1)+0.5*(O3(I)+O3(I-1))*(ZM(I-1)-ZM(I))
  20  CONTINUE
  
  
      CALL NETHT1(PRE,TZX,ZM,O3,O2,OXY,BO3,BO2
     & ,IMS,QCO2,QO3,HO23,QTZX,DEC,XPHI)


      DO 4 K=1,M$
      IK = M$ - K + 1
      HEAT(J,K) = HO23(IK)/DAYL
      COOL(J,K) =-(QCO2(IK)+QO3(IK))/DAYL
4     QNET(J,K) = QTZX(IK)/DAYL            ! net heating in K/s

5     CONTINUE

      RETURN
      END
      

      
CC ====== The original subroutines for the net radiative heating =======
      
      
      SUBROUTINE NETHT1(PRE,TEM,ZM,O3,O2,OXY,BO3,BO2
     & ,IM,QCO2,QO3,HO2O3,QNET,DEC,PHI)

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
     
C  Net heating rate in K/day  

      PARAMETER(IMS=40,IM1=15,IDM=5,IM2=IMS-IM1+IDM)

      COMMON /CURTZ2/ Z0(IMS),P0(IMS),T0(IMS),GDP(IMS),CURT0(IMS,IMS)
     &,A10(IMS,IMS),A01(IMS,IMS),A20(IMS,IMS),A02(IMS,IMS),A11(IMS,IMS)

      COMMON /CURTZ3/ R0(IMS),ZH00(IMS),ZGDP(IMS)
     &  ,ZCURT0(IMS,IMS),ZA01(IMS,IMS),ZA02(IMS,IMS)

      DIMENSION PRE(IM),TEM(IM),ZM(IM),O3(IM),O2(IM)
     &  ,OXY(IM),BO3(IM),BO2(IM),QCO2(IM),QO3(IM)
     &  ,HO2O3(IM),QNET(IM)


      DIMENSION WK1(IMS),WK2(IMS),WK3(IMS)

      jsep0=isw(19)

      DO 5 J=1,IMS
      DAIR=PRE(J)/(1.381E-23*TEM(J))              ! air density in m^-3
      WK1(J)=1.661*O3(J)/DAIR                     ! ozone mixing ratio
   5  CONTINUE
      CALL COLCO2(TEM,BOXY,QCO2,IMS,1,jsep0)  !changed jsep=1 to 1 JTB 4/15
c                                         !jcon=2 to 1 JTB 4/25/92
      CALL COLO3(TEM,WK1,QO3,IMS,2)

c      SUBROUTINE HEAT1(QH,DEC,PHI,ZM,O2,O3,BO2,BO3,PRE,TEM,
c     *  IM,IEFF2,IEFF3,FCLEAR,RCLOUD,RSURF,ITOP)
c      SUBROUTINE HEAT2(QH,DEC,PHI,ZM,O2,O3,BO2,BO3,PRE,TEM,IM,ITOP)


c      CALL HEAT1(HO2O3,DEC,PHI,ZM,O2,O3,BO2,BO3,PRE,TEM,
c     *  IMS,0,1,0.65,0.7,0.1,1)
      CALL HEAT2(HO2O3,DEC,PHI,ZM,O2,O3,BO2,BO3,PRE,TEM,
     *  IMS,1)


      DO 10 I=1,IMS
      QNET(I)=QCO2(I)+QO3(I)+HO2O3(I)      ! net heating rate in K/day
  10  CONTINUE
      RETURN
      END
      
C  =============== CO2 cooling : ===================================   

      SUBROUTINE COLCO2(T,BOXY,QCO2,IM,JCON,JSEP)
C  This subroutine calculates the cooling rate (QCO2) by Curtis matrix CURT
C  for given p, T, [O] and [CO2]. CURTZ are the off-line calculated Curtis
C  matrices. The units of [O] and p are  M^-3 and Pascal, respectively.. KMS=IM
C  JCON=1 --> linear interpolation; JCON=2 --> quadrature interpolation.
C  JSEP=0 --> LTE; JSEP=1 --> nonLTE in K1 region; JSEP>1 --> complete nonLTE
      PARAMETER (IMS=40,KM=IMS,KMP=KM+1,K1=14,K2=KM-K1,K1P=K1+1)               !  KM=IM
      COMMON /CURTZ2/ Z0(IMS),P0(IMS),T0(IMS),GDP(IMS),CURT0(IMS,IMS)
     &,A10(IMS,IMS),A01(IMS,IMS),A20(IMS,IMS),A02(IMS,IMS),A11(IMS,IMS)
      DIMENSION CURT(KM,KM),THE1(KM),EE(KM,KM),AINV(KM,KMP)
      DIMENSION BINV(K1,K1P),WOK(K2,K1),BAIR(KM),BOXY(KM)
C  EE(KM,KM) is used both for diagonal E matrix and non-LTE Curtix matrix
      DIMENSION T(IM),QCO2(IM)

      DO 10 J=1,IM
      QCO2(J)=0.0
      THE1(J)=THEZ(T(J))
      DO 10 K=1,IM
      EE(J,K)=0.0
  10  CURT(J,K)=CURT0(J,K)

      DO 20 J=1,IM
      DTJ=(T(J)-T0(J))
      DO 20 K=1,IM        ! Get the Curtis matrix by interpolation
      IF(ABS(CURT0(J,K)).LE.0.005) GO TO 20
      DELTA=DETIJ(T0,T,GDP,IM,J,K)
      CURT(J,K)=CURT0(J,K)+A10(J,K)*DTJ+A01(J,K)*DELTA
      IF(JCON.EQ.2) CURT(J,K)=CURT(J,K)
     & +0.5*(A20(J,K)*DTJ**2+A02(J,K)*DELTA**2)+A11(J,K)*DELTA*DTJ
  20  CONTINUE

      IF(JSEP.EQ.0) THEN
      DO 26 J=1,IM
      DO 26 K=1,IM
  26  EE(J,K)=CURT(J,K)
      GO TO 300         ! Skip non-LTE altogether
      ENDIF

      AA21=1.51
      AA32=1.94
      DO 28 J=1,IM
      BAIR(J)=P0(J)/(1.381E-23*T(J))              ! air density in m^-3
  28  CONTINUE
      DO 30 J=1,IM
      PHI1=AVT21(J,BAIR,BOXY,T,IM)/AA21
      PHI2=AVT32(J,BAIR,BOXY,T,IM)/AA32
      PHI=((SPT(T(J),1)+SPT(T(J),3))*PHI1
     &         +SPT(T(J),2)*PHI2)/SPT(T(J),4)
cc      PHI=0.247*P0(J)       ! Simple param. by Houghton(1986, p74)
  30  EE(J,J)=1.96E-3/PHI                     ! 1.96E-3=1/(2*255)

      IF(JSEP.EQ.1) GO TO 100

      DO 40 J=1,IM
      DO 40 K=1,IM
      XX=0.0
      DO 35 LL=1,IM
  35  XX=XX+CURT(J,LL)*EE(LL,K)
      AINV(J,K)=-XX
      IF(J.EQ.K) AINV(J,K)=1.0-XX
  40  CONTINUE
      CALL INVERT(AINV,KM,KMP)
      GO TO 200

 100  CONTINUE

      DO 140 J=1,K1
      DO 140 K=1,K1
      XX=0.0
      DO 135 LL=1,K1
 135  XX=XX+CURT(J,LL)*EE(LL,K)
      BINV(J,K)=-XX
      IF(J.EQ.K) BINV(J,K)=1.0-XX
 140  CONTINUE
      CALL INVERT(BINV,K1,K1P)
      DO 150 J=1,IM
      DO 150 K=K1P,IM
      AINV(J,K)=0.0
      IF(J.EQ.K) AINV(J,K)=1.0
 150  CONTINUE
      DO 160 J=1,K1
      DO 160 K=1,K1
 160  AINV(J,K)=BINV(J,K)
      DO 170 J=1,K2
      DO 170 K=1,K1
      XX=0.0
      DO 165 LL=1,K1
 165  XX=XX+CURT(K1+J,LL)*EE(LL,K)
      WOK(J,K)=XX
 170  CONTINUE
      DO 180 J=1,K2
      DO 180 K=1,K1
      XX=0.0
      DO 175 LL=1,K1
 175  XX=XX+WOK(J,LL)*BINV(LL,K)
      AINV(K1+J,K)=XX
 180  CONTINUE
 200  CONTINUE
      DO 50 J=1,IM
      DO 50 K=1,IM
      XX=0.0
      DO 45 LL=1,IM
  45  XX=XX+AINV(J,LL)*CURT(LL,K)
  50  EE(J,K)=XX
 300  CONTINUE
      DO 60 J=1,IM
      DO 55 K=1,IM
  55  QCO2(J)=QCO2(J)+EE(J,K)*THE1(K)
  60  CONTINUE
      RETURN
      END
  
      FUNCTION DETIJ(T0,T,GDP,IM,II,JJ)
CC  To calculate the DELTA(i,j) used for interpolations (FS81)
      DIMENSION T0(IM),T(IM),GDP(IM)
      DIMENSION DT(151)
      K1=MIN0(II,JJ)
      K2=MAX0(II,JJ)
      XX=0.0
      YY=0.0
      DO 10 K=K1,K2
      XX=XX+GDP(K)
  10  YY=YY+(T(K)-T0(K))*GDP(K)
      DETIJ=YY/XX
      RETURN
      END
  
      FUNCTION THEZ(T)
C  ratio of the Planck black-body function at T to 250
C  and at the wavenumber 6.75E4 m^-1.  47.612466=EXP(970.97/250)-1.0
      THEZ=47.612466/(EXP(970.97/T)-1.0)
      RETURN
      END

      SUBROUTINE INVERT(D,N,M)
CC  
CC  To invert the matrix D,  M=N+1. Mth row can be any values. Single precision
CC
      DIMENSION D(N,M)
      KK=0
      JJ=0
      DO 10 K=1,N
      DO 11 J=1,N
  11  D(J,M)=0.0
      D(K,M)=1.0
      JJ=KK+1
      LL=JJ
      KK=KK+1
  20  IF(ABS(D(JJ,KK)-1.0E-30)) 21,21,22
  21  JJ=JJ+1
      IF(JJ-N) 20,20,99
  99  WRITE(*,98)
  98  FORMAT('ERRORNEOUS INPUT')
      RETURN
  22  IF(LL-JJ) 23,24,23
  23  DO 25 MM=1,M
      DTEMP=D(JJ,MM)
      D(LL,MM)=D(JJ,MM)
  25  D(JJ,MM)=DTEMP
  24  DIV=D(K,K)
      DO 30 LJ=1,M
  30  D(K,LJ)=D(K,LJ)/DIV
      DO 12 I=1,N
      FAC=D(I,K)
      IF(I-K) 15,12,15
  15  DO 31 LJ=1,M
  31  D(I,LJ)=D(I,LJ)-FAC*D(K,LJ)
  12  CONTINUE
      DO 40 J=1,N
  40  D(J,K)=D(J,M)
  10  CONTINUE
      RETURN
      END

      FUNCTION VT1M(T)
      IF(T.LT.200.0) VT1M=2.5E-15/1.0E6
      IF(T.GE.200.0) VT1M=2.5E-15/1.0E6*(1.0+0.03*(T-200.0))
      RETURN
      END
      
      FUNCTION VT2M(T)
      TRA=T/273.3
      VT2M=1.24E-14*TRA*TRA/1.0E6
      RETURN
      END

      FUNCTION VT1N(T)
      COMMON /NLTE/ NCALL, NLTE$
      CHARACTER*3 NLTE$



      IF( NLTE$ .EQ. 'OLD') THEN

        T2=SQRT(T)
        T3=-76.75/T**0.3333333333333
        VT1N=2.32E-9*EXP(T3)+1.0E-14*T2
        VT1N=VT1N/1.0E6

        IF( NCALL .EQ. 0) WRITE(*,*) 'USED OLD CO2-O Q. '

      ENDIF




      IF( NLTE$ .EQ. 'NEW') THEN

        VT1N=5.5E-18    ! New increased CO2-O quenching

        IF( NCALL .EQ. 0) WRITE(*,*) 'USED NEW CO2-O Q. '

      ENDIF


      NCALL=NCALL+1

      RETURN
      END

c--------------------
c  Next 2 func.s changed or removed as a result of new rate 
c  coeff.s for CO2-O quenching (XZ 7/28/92; JTB 8/7/92)
c      


        
c----- OLD ( VT1N = fun[T] )-----

      FUNCTION AVT21(I,BAIR,BOXY,TEM,IM)
      DIMENSION BAIR(IM),BOXY(IM),TEM(IM)
      AVT21=BAIR(I)*VT1M(TEM(I))+BOXY(I)*VT1N(TEM(I))
      RETURN
      END

      FUNCTION AVT32(I,BAIR,BOXY,TEM,IM)
      DIMENSION BAIR(IM),BOXY(IM),TEM(IM)
      AVT32=BAIR(I)*VT2M(TEM(I))+BOXY(I)*VT1N(TEM(I))
      RETURN
      END


c----- NEW  -----

CC      FUNCTION AVT21(I,BAIR,BOXY,TEM,IM)
CC      DIMENSION BAIR(IM),BOXY(IM),TEM(IM)
CC      AVT21=BAIR(I)*VT1M(TEM(I))+BOXY(I)*5.5E-18
CC      RETURN
CC      END
      
CC      FUNCTION AVT32(I,BAIR,BOXY,TEM,IM)
CC      DIMENSION BAIR(IM),BOXY(IM),TEM(IM)
CC      AVT32=BAIR(I)*VT2M(TEM(I))+BOXY(I)*5.5E-18
CC      RETURN
CC      END


c------------



      FUNCTION SPT(T,N)
CC  
CC  To calculate the temperature dependence of band intensity S
CC  in Houghton's appedices. N=1,2,3 and 4 correspond to 626 
CC  fundamental band, 626 first hot band, isotopes fundamental 
CC  band and total 15 micron band, respectively.
CC
      X=T-250.0
      GO TO (10,20,30,40),N
  10  Y=5.0522-5.4693E-4*X-5.2873E-7*X*X
      GO TO 100
  20  Y=3.8067+7.0636E-3*X-2.1024E-5*X*X+1.964E-8*X*X*X
      GO TO 100
  30  Y=3.2338-5.5612E-4*X-5.3153E-7*X*X
      GO TO 100
  40  Y=5.0961+8.8837E-6*X+3.1662E-8*X*X
 100  SPT=10.0**Y
      RETURN
      END


C ======================== Ozone cooling : ===============================


      SUBROUTINE COLO3(T,RR,QO3,IM,JCON)
C  To calculate the O3 cooling rate by Curtis matrix CURT.
C  JCON=1 --> linear interpolation; JCON=2 --> quadrature interpolation.
      PARAMETER (IMS=40,KM=IMS,KMP=KM+1,K1=14,K2=KM-K1,K1P=K1+1)            !  IMS=IM
      COMMON /CURTZ3/ R0(IMS),ZH00(IMS),ZGDP(IMS)
     &  ,ZCURT0(IMS,IMS),ZA01(IMS,IMS),ZA02(IMS,IMS)
      DIMENSION CURT(KM,KM),THE1(KM),EE(KM,KM),AINV(KM,KMP)
      DIMENSION BINV(K1,K1P),WOK(K2,K1),BAIR(KM),BOXY(KM)
C  EE(KM,KM) is used both for diagonal E matrix and non-LTE Curtix matrix
      DIMENSION T(IM),RR(IM),QO3(IM)
      DO 10 J=1,IM
      QO3(J)=0.0
      THE1(J)=THEZ(T(J))
      DO 10 K=1,IM
  10  CURT(J,K)=ZCURT0(J,K)
      DO 20 J=1,IM
      DO 20 K=1,IM        ! Get the Curtis matrix by interpolation
      IF(ABS(ZCURT0(J,K)).LE.0.005) GO TO 20
      DELTA=DETIJR(R0,RR,ZGDP,IM,J,K)
      CURT(J,K)=ZCURT0(J,K)+ZA01(J,K)*DELTA
      IF(JCON.EQ.2) CURT(J,K)=CURT(J,K)+ZA02(J,K)*DELTA**2
  20  CONTINUE
      DO 26 J=1,IM
      FACH=RR(J)/R0(J)
      DO 26 K=1,IM
  26  EE(J,K)=CURT(J,K)*(FACH*ZH00(J))
      DO 60 J=K1P,IM
      DO 55 K=1,IM
  55  QO3(J)=QO3(J)+EE(J,K)*THE1(K)
  60  CONTINUE
      RETURN
      END
  
      FUNCTION DETIJR(R0,RR,GDP,IM,II,JJ)
CC  To calculate the DELTA(i,j) used for interpolations (ZHU, 1992)
      DIMENSION R0(IM),RR(IM),GDP(IM)
      DIMENSION DR(151)
      DELR0=0.15
      K1=MIN0(II,JJ)
      K2=MAX0(II,JJ)
      XX=0.0
      YY=0.0
      DO 10 K=K1,K2
      XX=XX+DELR0*R0(K)*GDP(K)
  10  YY=YY+(RR(K)-R0(K))*GDP(K)
      DETIJR=YY/XX
      RETURN
      END
  
C  ==================== O2 and O3 heating ================================


      SUBROUTINE HEAT1(QH,DEC,PHI,ZM,O2,O3,BO2,BO3,PRE,TEM,
     *  IM,IEFF2,IEFF3,FCLEAR,RCLOUD,RSURF,ITOP)
CCCCC
C       To calculate the daily averaged heating rate (QH: K/day) for O3 
C       and O2, following the parameterization scheme by Strobel (1978).
C       DEC=declination (degree), PHI=latitude (degree).
C       O2,O3=number density (m^-3), BO2,BO3=column density (m^-2).
C       Six bands: Chappius (7500-4500),Huggins (3550-2825),
C       Hadley (2725-2450),Herzberg continuum (2400-2060), 
C       Schumann-Runge continuum (1740-1260),Schumann-Runge (2025-1750).
C       ZM=altitude (m), PRE=pressure (Pa), TEM=temperature. IEFF2,   
C       IEFF3=1 (epsron=1) or 0 are the efficient factors for O2 and O3.
C       FCLEAR=Fraction of clear sky, RCLOUD,RSURF=reflectivities of 
C       clouds and surface for the Chappius band.
C       ITOP=1: The input vectors start from the top boundary. 
CCCCC
      PARAMETER (KM=100)
      DIMENSION QH(IM),ZM(IM),O2(IM),O3(IM)
     *   ,BO2(IM),BO3(IM),PRE(IM),TEM(IM)
      DIMENSION QHR(KM),ZMR(KM),O2R(KM),O3R(KM)
     *   ,BO2R(KM),BO3R(KM),PRER(KM),TEMR(KM)
      DATA OMEGA/0.2617993/,RAIR0/2.87E2/,CP/1.004E3/
      DATA SIGCH/2.85E-25/,FCH/3.7E2/,XCH/6.0/
      DATA SIGHA/8.8E-22/,FHA/5.013/,XHA/2.6/
      DATA SIGHU/1.3E-6/,FHU1/8.5E-2/FHU2/4.9E-2/,EPHU/0.74/
     *     ,SLM/1.27E-2/,SLONG/3.055E3/,SHORT/2.805E3/
      DATA SIGHE2/6.6E-28/,SIGHE3/4.9E-22/,FHE/1.2/,XHE/2.29/
      DATA SIGSRC/1.0E-21/,FSRC/1.1E-3/,EPSRC/0.41/
      IF(ITOP.EQ.1) GO TO 222
      DO 210 I=1,IM
      O2R(I)=O2(I)
      O3R(I)=O3(I)
      BO2R(I)=BO2(I)
      BO3R(I)=BO3(I)
      PRER(I)=PRE(I)
      TEMR(I)=TEM(I)
 210  CONTINUE
      DO 220 I=1,IM
      K=IM+1-I
      O2(I)=O2R(K)
      O3(I)=O3R(K)
      BO2(I)=BO2R(K)
      BO3(I)=BO3R(K)
      PRE(I)=PRER(K)
      TEM(I)=TEMR(K)
 220  CONTINUE
 222  CONTINUE
      RDEC=DEC/57.2958
      RPHI=PHI/57.2958
      DO 10 I=1,IM
      QH(I)=0.0
  10  CONTINUE
CCCCCCCCCC  Diurnal average quantities  CCCCCCCCCCCCCCCCCCCCC
C....T0: half day length (Cogley & Borucki, 1976)
C....UBAR: average direction cosine of sun zenith angle
C....U2BAR: average of squared direction cosines of sun zenith angle
C
      ACB=SIN(RDEC)*SIN(RPHI)
      BCB=COS(RDEC)*COS(RPHI)
      RATIO=ACB/BCB
      IF(RATIO.LE.-1.0) THEN
C                           The sun does not rise
      RETURN
      ENDIF 
      IF(RATIO.GE.1.0) THEN
C                             The sun does not set
      T0=12.0
      GO TO 100
      ENDIF
      T0=ACOS(-ACB/BCB)/OMEGA
 100  U1BAR=2.*(ACB*T0+BCB*SIN(OMEGA*T0)/OMEGA)
      U2BAR=2.*(T0*ACB*ACB+2.*ACB*BCB*SIN(OMEGA*T0)/OMEGA
     * +BCB*BCB*(0.5*T0+SIN(2.*OMEGA*T0)/(4.*OMEGA)))          
      B4=SQRT(2.0*T0/U2BAR)
      A4=0.5*U1BAR*B4
      FDAY=U1BAR*B4/24.0
C  To calculate the effective albedo ALBT for the Chappius band.
C  ALBC=Clear sky albedo and then Cloud sky albedo
      ALBC=0.219/(1.0+0.816/B4)
      ALBC=ALBC*FCLEAR+RCLOUD*(1.0-FCLEAR)
      ALBT=ALBC+(1.0-ALBC)*0.856*RSURF/(1.0-0.144*RSURF)
      DO 50 I=1,IM
      SNO3=BO3(I)*B4
      SNO2=BO2(I)*B4
C   Chappius Band
      ATT1=EXP(-AMIN1(70.0,SIGCH*SNO3))
      ATT2=EXP(-AMIN1(70.0,BO3(IM)*B4+(BO3(IM)-BO3(I))*1.9))
      QIND=O3(I)*FCH*SIGCH*(ATT1+2.*ALBT/B4*ATT2)
      IF(IEFF3.EQ.0) THEN
      EPS=1.0D0-0.08471*XCH
      QIND=QIND*EPS
      ENDIF
      QH(I)=QH(I)+QIND
C   Hartley Band
      ATT=EXP(-AMIN1(70.0,SIGHA*SNO3))
      QIND=O3(I)*FHA*SIGHA*ATT
      IF(IEFF3.EQ.0) THEN
      EPS=1.0D0-0.08471*XHA
      QIND=QIND*EPS
      ENDIF
      QH(I)=QH(I)+QIND
C   Huggins Bands
      ATT1=EXP(-AMIN1(70.0,SIGHU*SNO3*EXP(-SLM*SLONG)))
      ATT2=EXP(-AMIN1(70.0,SIGHU*SNO3*EXP(-SLM*SHORT)))
      QIND=(O3(I)/(SLM*SNO3))*
     &     (FHU1+(FHU2-FHU1)*ATT1-FHU2*ATT2)
      IF(IEFF3.EQ.0) THEN
      QIND=QIND*EPHU
      ENDIF
      QH(I)=QH(I)+QIND
C...OXYGEN ABSORPTION
C   Herzberg Continuum (Ozone and Oxygen)
      ATT=EXP(-AMIN1(70.0,SIGHE2*SNO2+SIGHE3*SNO3))
      QIND2=O2(I)*FHE*SIGHE2*ATT
      QIND3=O3(I)*FHE*SIGHE3*ATT
      IF(IEFF2.EQ.0) THEN
      EPS3=1.0D0-0.08471*XHE
      QIND2=QIND2*EPS2
      ENDIF
      IF(IEFF3.EQ.0) THEN
      EPS3=1.0D0-0.41305*XHE
      QIND3=QIND3*EPS3
      ENDIF
      QH(I)=QH(I)+QIND2+QIND3
C   Schumann-Runge Continuum
      ATT=EXP(-AMIN1(70.0,SIGSCR*SNO2))
      QIND=O2(I)*FSCR*SIGSCR*ATT
      IF(IEFF2.EQ.0) THEN
      QIND=QIND*EPSRC
      ENDIF
      QH(I)=QH(I)+QIND
        FIL=3.43E-3
	FIS=1.35E-3
	SGL=2.9E-23
	SGM=1.54E-22
	SGS=1.1E-21
        IF(IEFF2.EQ.0) THEN
        FIL=9.8E-4
	FIS=4.3E-4
	SGM=1.7E-22
	SGS=1.15E-21
	ENDIF
      ATT1=EXP(-AMIN1(70.0,SGL*SNO2))
      ATT2=EXP(-AMIN1(70.0,SGM*SNO2))
      ATT3=EXP(-AMIN1(70.0,SGS*SNO2))
      QIND=O2(I)*(FIL*ATT1+(FIS-FIL)*ATT2-FIS*ATT3)/SNO2
      QH(I)=QH(I)+QIND
C   Schumann-Runge Band
      ALIT=1.43E-5
      BLIT=9.64E6
      IF(IEFF2.EQ.0) THEN
      ALIT=6.7E-5
      BLIT=3.44E7
      ENDIF
      QIND=1.0E-7*O2(I)/(ALIT*SNO2+BLIT*SQRT(SNO2))
      IF(SNO2.LT.1.0E22) THEN
      QIND=9.03E-26*O2(I)
      IF(IEFF2.EQ.0) QIND=2.43E-26*O2(I)
      ENDIF
      QH(I)=QH(I)+QIND
      QH(I)=FDAY*QH(I)
  50  CONTINUE
      DO 60 I=1,IM
      RAIR=RAIR0
      IF(ZM(I).GT.8.0E4) THEN
      X=(ZM(I)-8.E4)/1.E3
      RAIR=RAIR0/(0.99935D0+1.5266D-3*X-
     1    1.5955D-4*X*X+1.5640D-6*X**3)
      ENDIF
      RHOI=PRE(I)/(RAIR*TEM(I))
      QH(I)=86400.0*QH(I)/(CP*RHOI)
  60  CONTINUE
      IF(ITOP.EQ.1) RETURN
      DO 230 I=1,IM
      QHR(I)=QH(I)
      O2(I)=O2R(I)
      O3(I)=O3R(I)
      BO2(I)=BO2R(I)
      BO3(I)=BO3R(I)
      PRE(I)=PRER(I)
      TEM(I)=TEMR(I)
 230  CONTINUE
      DO 240 I=1,IM
 240  QH(I)=QHR(IM+1-I)
      RETURN
      END

C  ==================== O2 and O3 heating [new] ========================

      SUBROUTINE HEAT2(QH,DEC,PHI,ZM,O2,O3,BO2,BO3,PRE,TEM,IM,ITOP)
CCCCC
C       To calculate the daily averaged heating rate (QH: K/day) for O3 
C       and O2, following the parameterization scheme by Strobel (1978).
C       DEC=declination (degree), PHI=latitude (degree).
C       O2,O3=number density (m^-3), BO2,BO3=column density (m^-2).
C       ZM=altitude (m), PRE=pressure (Pa), TEM=temperature.
C       FCLEAR=Fraction of clear sky, RCLOUD,RSURF=reflectivities of 
C       clouds and surface for the Chappius band.
C       ITOP=1 => Z(1)=TOP BOUNDARY;  ITOP=0 => Z(IM)=TOP BOUNDARY
CCCCC
      DIMENSION QH(IM),ZM(IM),O2(IM),O3(IM)
     &  ,BO2(IM),BO3(IM),PRE(IM),TEM(IM)
      RDEC=DEC/57.2958
      RPHI=PHI/57.2958
      ACB=SIN(RDEC)*SIN(RPHI)
      BCB=COS(RDEC)*COS(RPHI)
      BCB2=BCB*BCB
      DO 10 I=1,IM
  10  QH(I)=0.0
      XA=0.0
      XB=ACB+BCB
      XC=ACB-BCB
      IF(XB.LE.0.0) RETURN          !    The sun does not rise
      IF(XC.GT.0.0) XA=XC           !    The sun does not set
      XM=0.5*(XB+XA)
      XR=0.5*(XB-XA)
      DDX=0.5773502*XR
      BMU1=XM+DDX
      BMU2=XM-DDX
      IF(ITOP.EQ.1) TAUST=BO3(IM)
      IF(ITOP.EQ.0) TAUST=BO3(1)
      DO 30 I=1,IM
      FAC1=HTMU(BMU1,O2(I),O3(I),BO2(I),BO3(I),0.25,TAUST)
     &          /SQRT(BCB2-(BMU1-ACB)**2)
      FAC2=HTMU(BMU2,O2(I),O3(I),BO2(I),BO3(I),0.25,TAUST)
     &          /SQRT(BCB2-(BMU2-ACB)**2)
      QH(I)=XR*(FAC1+FAC2)/3.1415926
      RHOI=PRE(I)/(287.0*TEM(I))
      QH(I)=86400.0*QH(I)/(1004.0*RHOI)
  30  CONTINUE
      RETURN
      END

      FUNCTION HTMU(BMU,O2,O3,BO2,BO3,OME0,BO3T)
CCCCC
C       Heating rate (MKS units) for O3 and O2, following the
C       parameterization scheme by Strobel (1978). BMU=COS(THETA)
C       O2,O3=number density (m^-3), BO2,BO3=column density (m^-2). 
C       Six bands: Chappius (7500-4500),Huggins(3550-2825),
C       Hadley(2725-2450),Herzberg continuum(2400-2060), 
C       Schumann-Runge continuum(1740-1260),Schumann-Runge(2025-1750).
CCCCC
      SNO3=BO3/BMU
      SNO2=BO2/BMU
C...OZONE ABSORPTION
      CHAP=HT12(370.0,2.85E-25,SNO3,1.0,1)              !   Chappius band
      HART=HT12(5.013,8.80E-22,SNO3,1.0,1)              !   Hartley band
      HUGG2=HT12(-3.6E-2,1.8364E-23,SNO3,0.0127,2)      !   Huggins band
      HUGG3=HT12(-4.9E-2,4.3939E-22,SNO3,0.0127,2)      !   Huggins band
      HUBB=6.693/SNO3+HUGG2+HUGG3
      HTSUM1=O3*(CHAP+HART+HUBB)
      ATT=EXP(-AMIN1(70.0,6.6E-28*SNO2+4.9E-22*SNO3))
      HTSUM2=1.2*(6.6E-28*O2+4.9E-22*O3)*ATT            !   Herzberg band
C...OXYGEN ABSORPTION
      SRC1=HT12(1.1E-3,1.0E-21,SNO2,1.0,1)              !   Schumann-Runge continuum 1
      SRC21=HT12(9.8E-4,2.9E-23,SNO2,1.0,2)             !   Schumann-Runge continuum 1
      SRC22=HT12(-5.5E-4,1.7E-22,SNO2,1.0,2)            !   Schumann-Runge continuum 1
      SRC23=HT12(-4.3E-4,1.15E-21,SNO2,1.0,2)           !   Schumann-Runge continuum 1
      HTSRC=O2*(SRC1+SRC21+SRC22+SRC23)
      HTSRB=1.0E-7*O2/(6.7E-5*SNO2+3.44E7*SQRT(SNO2))
      IF(SNO2.LT.1.0E22) HTSRB=2.43E-26*O2
C  Heating rate due to diffuse scatted solar radiation
      TAU=2.85E-25*BO3
      TAUS=2.85E-25*BO3T
      FAC1=TAUS/BMU+1.9*(TAUS-TAU)
      HTSCAT=2.109E-22*O3*OME0*BMU*EXP(-AMIN1(70.0,FAC1))       !   2.11E-22= 2*370*2.85E-25
      HTMU=HTSUM1+HTSUM2+HTSRC+HTSRB+HTSCAT
      RETURN
      END

      
      FUNCTION HT12(F1,SIG,SBN,SM,ICON)
      IF(ICON.EQ.1) HT12=F1*SIG*EXP(-AMIN1(70.0,SIG*SBN))
      IF(ICON.EQ.2) HT12=F1*EXP(-AMIN1(70.0,SIG*SBN))/(SBN*SM)
      RETURN
      END
      


      
