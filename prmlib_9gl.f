C***************************************************************
C* ------------  p r m l i b . f     ------------------------  *
C*                                                             *
C*         CONTAINS TRULY, STUPID PARAMETERIZATIONS            *
C*         AS OPPOSED TO ONLY MODERATELY STUPID (I.E.          *
C*         ``SELF-CONSISTENT'') ONES.                          *
C*                                                             *
C***************************************************************


      SUBROUTINE TROPHT

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'

C GLHDAYC(N$,M$) from LHGCM, WLHDAYC, HACKDAYC from WACCM latent heating, all in K/day - EF 10/07
C G4LHDAYC(N$,M$) is from GEOS 4 - EF 10/08 - THIS IS NOW DONE in RADIATE (5/2009)

      COMMON/CLATH/GLHDAYC(N$,M$), WLHDAYC(N$,M$), HACKDAYC(N$,M$),
     >             G4LHDAYC(N$,M$), LTNHT(N$,M$),  TDIFFDAY(N$,M$),
     >             G4LHDAYCX(NP$,MP$)


cjer added
      real lh_hgt, LTNHT, lhfac

c      real latent(12,n$,5)
      dimension yps(n$),dum5(n$,m$) 

      INTEGER IDAYS(13)
      DATA IDAYS/1,32,60,91,121,152,182,213,244,274,305,335,365/


cjer  use eric flemings data set (I interpolated to dynamics grid)
c which is in K/day
c      open (unit=529,file='latent_fleming.dat',form='unformatted') 
c      read(529) latent
c      close (unit=529)
c      
c      do 50 izz=1,5
c        do 50 ill=1,N$
c           ltnht(ill,izz) = latent(month,ill,izz)/dayl
c   50 continue
c
cjer use Julio's latent heating

      Q_LH  = LH_SW(1)*.01
      HL_Y  = LH_SW(2)*(-1.)
      HL_Z  = LH_SW(3)*(-1.)
      HL_YY = LH_SW(4)*LH_SW(4)*(-1.)
      HL_ZZ = LH_SW(5)*LH_SW(5)*(-1.)
      HL_YZ = LH_SW(6)*LH_SW(7)*(-1.)
      LH_HGT = LH_SW(8)*(1.0)

      IF ( LH_SW(11) .EQ. 0 ) IEXTLH=0
      IF ( LH_SW(11) .EQ. 1 ) IEXTLH=1


C     WRITE(6,*)'  USING JTB LATENT HEAT FOR TROPO '

      DO J=1,N$
         YPS(J)=YP(J)-0.33*DECL   ! NOT QUITE SUBSOLAR POINT
      END DO

      IF (LH_SW(12) .EQ. 1) THEN 
         DO J=1,N$
            YPS(J)=YP(J)          ! REMOVE TIME VARIATION
         END DO
      END IF

      DO K=1,M$
      DO J=1,N$ 
         LTNHT(J,K)=0.
      END DO
      END DO

      DO K=1,M$
      DO J=1,N$ 
         IF((ZP(K) .LT. LH_HGT ) .AND. (ZP(K) .GE. 0.0)) THEN

            LTNHT(J,K) = 
     >            1. + ABS( YPS(J) )/HL_Y + ZP(K)/HL_Z
     >               + (YPS(J)**2)/HL_YY + (ZP(K)**2)/HL_ZZ
     >               + (ABS(YPS(J))*ZP(K))/HL_YZ

            LTNHT(J,K)=AMAX1( LTNHT(J,K), 0.000 )

         ENDIF
      END DO
      END DO

cjer dum5 is the mid-latitude lobe

      DO K=1,M$
      DO J=1,N$
         DUM5(J,K)=0. 
         IF ((k.EQ.2) .OR. (k .EQ.2)) THEN 
            DUM5(J,K) = EXP( -ABS( 45.- ABS(YP(J)) )/30.)
         END IF
 
         DUM5(J,K)=AMAX1( DUM5(J,K), 0.000 )
      END DO
      END DO


      DO K=1,M$
      DO J=1,N$ 
         LTNHT(J,K)=LTNHT(J,K) + DUM5(J,K)*IEXTLH
      END DO
      END DO

      DO K=1,M$
      DO J=1,N$ 
         LTNHT(J,K)=(Q_LH/DAYL)*LTNHT(J,K)
cjer add to make sure no latent heat above top of zone
         IF(ZP(K) .GE. LH_HGT ) LTNHT(J,K)=0. 
      END DO
      END DO

cjer  looks like LTNHT is now in K/sec
c
cjer  end of choices for ltnht


C     Calculate the tropospheric temperature field for current day (DUM3)
C     by interpolating NMC fields.

      DAYM=IDAYX
      DX=1./(IDAYS(MONTH+1)-IDAYS(MONTH))
      DX2=real(IDAYX)/(IDAYS(MONTH+1)-IDAYS(MONTH))

      DO 2111 K=1,M$
         DO 2111 J=1,N$
            DUM3(J,K) = TFROP(J,K,MONTH) + 
     1                  (TFROP(J,K,MONTH+1)-TFROP(J,K,MONTH))*DX2
 2111    CONTINUE


C     Write out the tropospheric temperature field (to file fort.51).
    
c      CALL OUTA(51,NSTEP,NDAY,NYEAR,YP,ZP,DUM3,N$,M$,0)
 

C     Perform 3-point smoothing on latent heating field before adding it
C     to the overall heating.  LTNHT is converted from K/day to K/sec before
C     adding to HEAT.
 
cccccc   CALL SMTHR3(LTNHT,N$,M$)


CCelf  RELAXATION to NMC SURFACE TEMPS NOW  DONE IN TWODS.f and RADIATE.f


cjer  now convert to deg. theta / sec, first, overwrite  LTNHT(N$,M$) w/ GLHDAYC(N$,M$) from LHGCM (K/day)
C
C  use WACCM both latent heating here - WLHDAYC(N$,M$), HACKDAYC(N$,M$) are in K/day, convert to K/sec
C  to account for no clouds in radiation, add 20% increase in LH  LE  6 km, ramp down to NO INCREASE at 10 km
C  this increase warms the troposphere, maximum in tropics of 2.5K at 8 km, decreases to 1K at 4 km and 15 km
C                 and 1-1.5K warming at mid-high latitude troposphere - NO, DON'T USE
C
C  use TD Latent heating from WACCM (seasonal cycle) + GEOS 4 trend/CO2 sensitivity - G4LHDAYC(N$,M$) is in K/day
C                                               and set to 0.0 above 25 km  - EF 10/08
      DO K=1,M$
      DO J=1,N$ 
          
ccwaccm            lhfac = 1.2 - (ZP(K) - 6.)/20.
ccwaccm            if (lhfac .le. 1.) lhfac = 1.
ccwaccm            if (lhfac .ge. 1.2) lhfac = 1.2
ccwaccm            LTNHT(J,K) = (WLHDAYC(J,K) + HACKDAYC(J,K))/86400.    ! *lhfac

            LTNHT(J,K) = G4LHDAYC(J,K)/86400.
            if (ZP(K) .GE. 25.) LTNHT(J,K) = 0.0

ccpw32
ccpw32 - THIS IS NOW DONE in RADIATE (5/2009) IE,  
CCPW32        THIS IS NO LONGER USED HERE, LTNHT is loaded in here, but NOT USED
ccpw32
ccpw32            HEAT(J,K) = HEAT(J,K) + LTNHT(J,K)*TTOTH(K)
ccpw32

      END DO
      END DO


cjer add writeout of latent heat
c      if(daym.eq.15.) then
c        write(68) daym
c        write(68) ltnht
c      endif

c                 
cc       write (772) daym, N$, M$
cc       write (772) ltnht
cc       write (772) glhdayc


      RETURN
      END




      SUBROUTINE GET_KSS

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'

      DIMENSION YPS(NP$)
      REAL KYYMAX,KYYMIN,MESKZZ,MESOPS,MLLAT


C  fixed model Kzzs interpolated to coupled model dynamics grid for current day:  XKZZFX(NP$,MP$) in cm2/sec
C                                                   interpolated in TWODS
      COMMON/CKZZ1/XKZZFX(NP$,MP$)

ccelf      REAL*4 KZZTROPC(NP$,MP$)


cjer      LOGICAL TRINIT
c      REAL    TROPPR(91, 12), LATPR(91)
c      REAL    TROPHT(NP$, 12)
c      DATA    TRINIT /.FALSE./

C ****  SEARCH FOR SUB-SOLAR PT.    
   

      DO J=1,NP$
         YPS(J)=YPP(J)-0.5*DECL
      END DO


C **** DEFINE STRAT. TROP. BOUNDARIES FOR KZZ

      PR_PBLR=1.*KZZ_SW(3)/1000.
      PR_TROP=1.*KZZ_SW(2)/1000.
      PR_STRT=1.*KZZ_SW(1)/1000.
      PR_MESO=1.*KZZ_SW(4)/1000.
      PR_CONV=1.*KZZ_SW(5)/1000.

      EQLAT = 1.*KZZ_SW(7)   ! ABS(LAT) LIMIT OF `EQUATORIAL ZONE'
      MLLAT = 1.*KZZ_SW(8)   ! ABS(LAT) LIMIT OF `MID-LAT. ZONE'

      TRPPEQ=KZZ_SW(9)       ! EQ. TROPOPAUSE HT. (KM)
      TRPPML=KZZ_SW(10)      ! MID-LAT. TROPOPAUSE HT. (KM)
      BDYLYR=KZZ_SW(6)       ! ``BOUNDARY LAYER'' HT. (KM)   
      MESOPS=KZZ_SW(11)      ! MESOPAUSE HT. (KM)
      PGWKZZ=KZZ_SW(12)/100. ! GRAVITY WAVE MIXING FACTOR - looks like KZZ_SW(12)=50  (in KZZ_PARAM.DAT) (EF)
      
      PBLKZZ=10.**PR_PBLR    ! LOWER TROP. KZZ VAL. CM^2/S
      TRPKZZ=10.**PR_TROP    ! UPPER TROP. KZZ VAL. CM^2/S
      STRKZZ=10.**PR_STRT    ! STRAT. KZZ VAL. CM^2/S   
      MESKZZ=10.**PR_MESO    ! MESO.-THERM. KZZ VAL. CM^2/S   
      CNVKZZ=10.**PR_CONV    ! CUMULUS CONVECTION KZZ CM^2/S   




C            ****************
C            * *** KZZ ***  *
C            ****************
C
C
C  PREVIOUS Stuff commented out here, replaced with TROPKZZ methodology as in FIXED model, 
C     based on COUPLED model TBAR, UBAR, and adapted to COUPLED model dynamics grid  - OCTOBER 2007 - EF
C
C   just use lookup table from FIXED model 9fr run:
C
                                               !  KZZTROPC(NP$,MP$)
ccccc           CALL TROPKZZ_COUP(KZZTROPC)


ccelf      DO K=1,MP$
ccelf         DO J=1,NP$ 
ccelf            KZZ(j,k)=STRKZZ   ! 2.0E3
ccelf            END DO
ccelf         END DO

cjer added from Julio's code
      RMAXKZZ=-10000. 


C  for 9GC, DON'T use Julio's  Gravity Wave Kzz ANYWHERE - get negatives and it's VERY BLOTCHY
C      just use FIXED model lookup table  - XKZZFX(NP$,MP$)
C                    
C  for 9GL, ramp down from 16 km - 5 km as upper limit, set max of 10.e4 in lower trop, 100.e4 elsewhere
C     also ensure minimum of .001*1.e4 here for DYNAMICS - this also done in XCOUPLED for the CHEMISTRY
C                                   first find indicies closest to 5 km and 15 km of current dynamics grid
C
C     THIS IS ALL NOW DONE IN XCOUPLED, to be CONSISTENT w/ CHEMISTRY KZZ, so here just recheck limits
C

cc        i5k = 3
cc        i15k = 8
cc
cc        do 877 ik=1,MP$
cc             if (zpp(ik) .lt. 5.4)   i5k = ik
cc             if (zpp(ik) .lt. 15.4) i15k = ik
cc 877    CONTINUE


      DO K=1,MP$
      DO J=1,NP$ 

          KZZ(J,K) = XKZZFX(J,K)

          xkzmx = 100.e4

ccccccc          if (zpp(k) .le. 4.5) xkzmx = 10.e4
ccccccc
ccccccc          if (ABS(ypp(j)) .le. 60.) then
ccccccc            if (zpp(k) .gt. 4.5  .and.  zpp(k) .le. 16.) then 
ccccccc              xkzz5 = 10.e4
ccccccc              if (XKZZFX(J,i5k) .lt. 10.e4) xkzz5 = XKZZFX(J,i5k)
ccccccc
ccccccc            xkzmx = exp( alog(XKZZFX(J,i15k)/xkzz5)/(zpp(i15k)-zpp(i5k))
ccccccc     >            *(zpp(k)-zpp(i5k)) + alog(xkzz5) )
ccccccc            endif
ccccccc          endif

          if (KZZ(J,K) .ge. xkzmx) KZZ(J,K) = xkzmx

          xkzmn = .001*1.e4
          if (KZZ(J,K) .le. xkzmn) KZZ(J,K) = xkzmn


ccelf         IF(ABS(YPS(J)) .LE. EQLAT) THEN
ccelf               IF (ZPP(K) .LT. TRPPEQ)       KZZ(J, K) = TRPKZZ
ccelf               IF (ZPP(K) .LT. BDYLYR)       KZZ(J, K) = PBLKZZ
ccelf               IF (ZPP(K) .LT. (0.4*BDYLYR)) KZZ(J, K) = CNVKZZ
ccelf            ELSE
ccelf               IF (ZPP(K) .LT. TRPPEQ)       KZZ(J, K) = TRPKZZ
ccelf               IF (ZPP(K) .LT. TRPML)        KZZ(J, K) = TRPKZZ
ccelf               IF (ZPP(K) .LT. (0.2*BDYLYR)) KZZ(J, K) = PBLKZZ
ccelf           ENDIF
            
C           Add in Gravity Wave Kzz - looks like PGWKZZ = 0.5

ccelf         GWKZZ0 = PGWKZZ * GWKZZ(J,K)
ccelf         IF (ZPP(K) .GE. TRPPEQ) THEN
ccelf           KZZ(J,K) = AMAX1(STRKZZ,GWKZZ0)
ccelf         ENDIF

         END DO
         END DO


c                                    REAL KZZ(NP$,MP$), YPP(NP$),ZPP(MP$)
ccelf       write (773) NP$, MP$
ccelf       write (773) YPP, ZPP
ccelf       write (773) kzz


ccelf        print *,  '    '
ccelf        print *,  'in GET Ksss:    '
ccelf        print *, EQLAT, STRKZZ, TRPPEQ, TRPML, PGWKZZ
cc        print *,  '    '





ccelf      CALL TRPCONV(CNVKZZ)   !  KLUGY PARAM. OF TROPICAL/SUMMERTIME - DON'T USE - EF 10/07

C                                                  ! kzz_sw(13) appears to be 0, so no smoothing is done here

      if(kzz_sw(13) .gt. 0.) then
        DO NSM=1,KZZ_SW(13)
           CALL SMTHR3( KZZ, NP$, MP$ )
        END DO
      endif


      CALL INRINST  ! CHECK FOR INERTIAL INSTABILITY


      NSM_ZI=INT( N$/24 )+1
    
cjer Julio wasn't sure whether zikzz was actually used
      DO N=1,NSM_ZI
         CALL SMTHR3( ZIKZZ , NP$, MP$ )
         CALL SMTHR3( ZIKYY, N$, M$ )
      END DO

      CALL OUTA(12,NSTEP,NDAY,NYEAR,YPP,ZPP,KZZ,  NP$,MP$,0)
      CALL OUTA(10,NSTEP,NDAY,NYEAR,YP, ZP, KYY,  N$, M$ ,0)
c      CALL OUTA(15,NSTEP,NDAY,NYEAR,YP, ZP, ZIKYY,N$, M$ ,0)
c      CALL OUTA(14,NSTEP,NDAY,NYEAR,YPP,ZPP,ZIKZZ,NP$,MP$,0)


C            ***************************
C            * *** MOL. DIFFUSION ***  *
C            ***************************

      
      DO K=1,M$
         DMOL(K) = 1./RHO(K)
         END DO


C            *************************
C            * *** KYY OVERRIDE ***  *
C            *************************

      IF (ISW(8) .GT. 1) THEN
         WRITE(6,  * ) ' OVER-RIDING KYY PARAMETRIZATION '
         WRITE(6,  * ) ' WITH UNIFORM KYY OF : | | | '
         WRITE(6,  * ) ' ...                   V V V '

         KYYMAX=-1.0E12
         KYYMIN= 1.0E12


         DO K=1,M$
         DO J=1,N$
            KYY(J,K)=10.**( ISW(8) )
         END DO
         END DO
         DO K=1,M$
         DO J=1,N$
            KYYMAX=AMAX1(KYYMAX, KYY(J,K) )
         END DO
         END DO
         DO K=1,M$
         DO J=1,N$
            KYYMIN=AMIN1(KYYMIN, KYY(J,K) )
         END DO
         END DO

         WRITE(6,5000) KYYMIN,KYYMAX

 5000    FORMAT( '  KYYMIN= ',E10.4,' KYYMAX= ',E10.4 )

       ENDIF

      RETURN
      END



      SUBROUTINE TRPCONV(CNVKZZ)
      
cjer  this is still using Paul Meade's latent heat form

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'
      REAL LTNHT1(NP$,MP$),YPS(NP$),ZPS(MP$),FHEAT,QCALC1,QCALC2,QCALC3
      REAL HTROP(NP$), ALTMAX(NP$), DELALT(NP$)


      HL_Y  = LH_SW(2)*(-1.)
      HL_Z  = LH_SW(3)*(-1.)
      HL_YY = LH_SW(4)*LH_SW(4)*(-1.)
      HL_ZZ = LH_SW(5)*LH_SW(5)*(-1.)
      HL_YZ = LH_SW(6)*LH_SW(7)*(-1.)

      LH_HGT= LH_SW(8)*(1.0)




      DO J=1,NP$
          YPS(J) = YPP(J)
          HTROP(J) = 25.0 - (YPS(J)/30.0)**2.0
          ALTMAX(J) = 0.5*(2.0 + HTROP(J))
          DELALT(J) = 0.5*(HTROP(J) - 2.0)
      END DO

      DO K=1,MP$
         ZPS(K) = ZPP(K)
         END DO

      DO K=1,MP$
      DO J=1,NP$ 
         LTNHT1(J,K)=0.
      END DO
      END DO

      FHEAT = 1.0

      DO K=1,MP$
      DO J=1,NP$ 
         IF((ZPP(K) .LT. LH_HGT ) .AND. (ZPP(K) .GE. 0.0)) THEN

            QCALC1 = 0.0

            IF (YPS(J) .LE. 10.0) THEN
                  QCALC2 = 0.425*(1.0 - (DECL/23.5))*
     2                          (1.0 - ((ZPS(K) - 8.0)/8.0)**2.0 - 
     3                                 (ABS(YPS(J) - 10.0)/140.0))
               ELSE
                  QCALC2 = 0.425*(1.0 - (DECL/23.5))*
     2                          (1.0 - ((ZPS(K) - 8.0)/8.0)**2.0 - 
     3                                 ((YPS(J) - 10.0)/25.0)**2.0)
               ENDIF

            IF (YPS(J) .GE. -10.0) THEN
                  QCALC3 = 0.225*(1.0 + (DECL/23.5))*
     2                          (1.0 - ((ZPS(K) - 8.0)/8.0)**2.0 - 
     3                                 (ABS(YPS(J) + 10.0)/140.0))
               ELSE
                  QCALC3 = 0.225*(1.0 + (DECL/23.5))*
     2                          (1.0 - ((ZPS(K) - 8.0)/8.0)**2.0 - 
     3                                 ((YPS(J) + 10.0)/25.0)**2.0)
               ENDIF

            IF (QCALC1 .GT. 0.0) THEN
               LTNHT1(J,K) = LTNHT1(J,K) + QCALC1
               ENDIF

            IF (QCALC2 .GT. 0.0) THEN
               LTNHT1(J,K) = LTNHT1(J,K) + QCALC2
               ENDIF

            IF (QCALC3 .GT. 0.0) THEN
               LTNHT1(J,K) = LTNHT1(J,K) + QCALC3
               ENDIF

            LTNHT1(J,K)=(CNVKZZ)*FHEAT*AMAX1( LTNHT1(J,K), 0.000 )

         ENDIF

         KZZ(J,K) = AMAX1( KZZ(J,K), LTNHT1(J,K) )     

      END DO
      END DO


      RETURN
      END


      SUBROUTINE INRINST

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'
        
cjer        REAL      FSCALE(N$,M$)
c        LOGICAL   FSINIT
c        DATA      FSINIT /.FALSE./


c        IF (.NOT.FSINIT) THEN
c           DO K=1, M$
c              DO J=1, N$
c                 ABSLAT = ABS(YP(J))
c                 IF (ABSLAT .LE. 80.0) THEN
c                       FSCALE(J,K)= 1.0
c                    ELSE
C                       FSCALE(J,K)= EXP(-1.0*(((ABSLAT-80.0)/5.0)**2.0))
C                       FSCALE (J, K) = 0.0
c                       FSCALE (J, K) = 1.0
c                    ENDIF
c                 END DO
c              END DO
c           FSINIT = .TRUE.
c           ENDIF


      ZINLM= 0.001*TOMEG*TOMEG
      ZINKYY=((N$/24.)**2)*5.0E10
      ZINKZZ=1.0E6

      DO K=1,M$
      DO J=1,N$
         DUM1(J,K) = ( UBX(J+1,K)-UBX(J,K) )/ (DY)
         DUM2(J,K) = CF(J)*( CF(J) - DUM1(J,K) )
         DUM3(J,K) = 0.0
      END DO
      END DO

      DO K=1,MP$
      DO J=1,NP$
         DUM1X(J,K)=0.
         DUM2X(J,K)=0.
      END DO
      END DO

      DO K=1,M$
      DO J=1,N$
         IF(DUM2(J,K).LE.ZINLM) THEN
           ZFAC=ZINLM-DUM2(J,K)
           ZFAC=ZFAC/(ZINLM*10.)
           ZFAC=1.-EXP(-ZFAC)
cjer           DUM1X(J,K)=ZINKYY*ZFAC*FSCALE(J,K)
c           DUM3(J,K)=ZINKZZ*ZFAC*FSCALE(J,K)
           DUM1X(J,K)=ZINKYY*ZFAC
           DUM3(J,K)=ZINKZZ*ZFAC
         ENDIF
         ZIKYY(J,K)=DUM1X(J,K)

      END DO
      END DO


      CALL MRS2RBR( DUM3, DUM2X )
      
      DO K=1,MP$
      DO J=1,NP$
           ZIKZZ(J,K)=DUM2X(J,K)
      END DO
      END DO

      RETURN
      END

      SUBROUTINE PBLPAR(PBLDRG)

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'
      REAL     PBLDRG(NP$,MP$),RHST(MP$)

      KTOPPBL=3
      PBLDRG0=1./(DAYL*2.0)

      DO K=1,MP$
         RHST(K)=1.000
      END DO

      DO J=1,NP$
         K=1
         PBLDRG(J,K) = ( (RHST(K+1)*UBX(J,K+1)-RHST(K)*UBX(J,K))
     >                 - (RHST(K)*UBX(J,K)) )*PBLDRG0*CST(J)
      END DO

      DO K=2,KTOPPBL-1
      DO J=1,NP$
         PBLDRG(J,K) =  PBLDRG0*CST(J)*(
     >     ( RHST(K+1)*UBX(J,K+1)-RHST(K)*UBX(J,K) )
     >   - ( RHST(K)*UBX(J,K)-RHST(K-1)*UBX(J,K-1) )
     >                       )
      END DO
      END DO
 
      RETURN
      END

      SUBROUTINE W_AERO(J00, PAR_RAD, WAERO ,NTR)

C**********************************************************
C    CALCULATES AEROSOL TERMINAL FALL SPEEDS              *
C    USING STOKES LAW DRAG DOWN LOW AND KINETIC           *
C    LAW DRAG UP HIGH WITH MATCH IN BETWEEN               *
C                                                         *
C    -- KINETIC DRAG:                                     *
C          (Reid, 1975, JAS, Vol. 32, No. 3 ;             *
C           Jensen & Thomas, 198[8,9] JGR Special Issue   *
C           on Noctilucent Clouds )                       *
C**********************************************************

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'
  
      DIMENSION WAERO(NP$,MP$),V1(MP$),V2(MP$)


      JM1 = MAX0( J00-1, 1  )
      JP1 = MIN0( J00+1, N$ )
      J0  = MIN0( J00  , N$ )

     
      IF(PAR_RAD .LE. 0.000) THEN

        DO K=1,M$
           WAERO(J00,K)=0.000
        END DO

        IF(J00.EQ.1) WRITE(6,101) NTR 


        RETURN
      END IF


      IF(J00.EQ.1) WRITE(6,100) NTR, PAR_RAD
 
 100  FORMAT( '   (#',I2,') -------  THIS IS AN AEROSOL OF ',
     1 F6.2,' uM PARTICLES' )
 101  FORMAT( '   (#',I2,')    ----  NOT AN AEROSOL ' )




      BOLTZ    =  1.381*(1.0E-16)       ! BOLTZMANN CONST. (ERGS/K)
      WMOL_AIR =  1.673*(1.0E-24)*28.9  ! MEAN MOL. WGT. (GRAMS)
      RADCM  =    PAR_RAD*(1.0E-4)      ! PARTICLE RADIUS (CM)


C -------------------   STOKES LAW DRAG 

      DO K=1,M$
      DO J=JM1,J0

         VISCO  = 0.00018 * SQRT( T(J,K)/293. )

         WFALL1 = (2./9.)*( GZ*(RADCM**2)/VISCO )

         DUM1(J,K)=WFALL1

      END DO
      END DO


C -------------------   KINETIC LAW DRAG 

      DO K=1,M$
      DO J=JM1,J0


         PRDYN  = 1000.*1000.*EXP(-ZZ(K)/H )

         V_LIKE = SQRT( T(J,K)*BOLTZ / WMOL_AIR )


         WFALL2 = SQRT(PI/2.0) * GZ*RADCM*V_LIKE / (2.0*PRDYN)
         

         DUM2(J,K)=WFALL2

      END DO
      END DO


C -----------------  MATCH THE TWO DIFFERENT FALL
C -----------------  SPEEDS ACROSS SOME LAYER



      IF(ISW(27).EQ.1) ZMATCH=30.0

      IF(ISW(27).EQ.2) THEN
        ZMATCH=0.0

        DO K=1,M$
          FRE_PTH=1.E-5*EXP( ZZ(K)/H )  ! MEAN MOL. FREE PATH (CM)
          IF(RADCM.GT.FRE_PTH) ZMATCH=ZP(K) 
        END DO          

        IF(J00.EQ.1) WRITE(6,102) ZMATCH, RADCM
        
 102    FORMAT('   ---- KINETIC DRAG LAW STARTS AT Z=' , F6.2
     2        , 'KM   FOR PARTICLE RADIUS=' , G14.6 , ' CM' )
      END IF
  

      DO K=1,M$
      DO J=JM1,J0
 

         WT1=1.00

         IF( ZP(K) .GT. ZMATCH ) THEN 
             WT1=EXP( (ZMATCH-ZP(K))/8.0 )
         ENDIF

         WT2=1.00-WT1

         DUM3(J,K)=WT1*DUM1(J,K)+WT2*DUM2(J,K)
 
      END DO
      END DO

      
      DO K=1,M$
         WAERO(J00,K)=( DUM3(J0,K)+DUM3(JM1,K) )/2.0
      END DO

      DO K=1,M$
         WAERO(J00,K)=AMIN1( WAERO(J00,K), 10.00 )
      END DO

      
      if( J00 .eq. INT(N$/2) ) then
     
      DO K=1,MP$
         V1(K)=WAERO(J00,K)
         V2(K)= DUM2(J00,K)
      END DO

c      CALL OUTA(48,NSTEP,NDAY,NYEAR,YPP(J00),ZP,V1,1,M$,0)
c      CALL OUTA(49,NSTEP,NDAY,NYEAR,YP(J00) ,ZP,V2,1,M$,0)
      end if

      
      RETURN
      END

C
C  SETMASS was taken from trnsp_9ga.f (Prather scheme), no longer used
C  Lin-Rood now used for theta, ang mom transport - EF, Mar 2008
C
      SUBROUTINE SETMASS

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMONT.INC'

C      WRITE(6,'("  RESETTING SM ")' )

      DO J=1,MP$
      DO I=1,NP$
        SM(I,J)=RHO0(I,J)
      END DO
      END DO

      RETURN
      END
