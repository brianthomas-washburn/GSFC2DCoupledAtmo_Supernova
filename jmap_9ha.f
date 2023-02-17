c
C   routine to load diurnally varying Js from the GSFC array JZ(18,PH$,L$,Z$) 
C     into the array JEX for diurnal averaging  
C     (it's then reloaded into the GSFC array JA in UNJMAP)
C     also load into the AER array JMOM(NJR)
C
C   FOR 9CP run use insolation corrected JZ(18,PH$,L$,Z$) array from PDFINSOL routine
C
C      THIS IS UPDATED for 9HA run
C

	SUBROUTINE JMAP(IJC,IKC,ICLOCK)


        include "com2d.h"
        include "comphot.h"

	INTEGER NJR,NKR,NRP,NWLI,NW1,NW2

        PARAMETER (NJR=121,NKR=382,NRP=5,NWLI=77,NW1= 1,NW2=77)
        DOUBLE PRECISION JRMOM
	REAL*8 JZ
        COMMON/JRATE/JRMOM(NJR)

        COMMON/JCORRECT/JZ(18,PH$,L$,Z$)

	REAL*8 JEX
        COMMON/JEXAER/JEX(3)


C
C   common for CGCM-COMBO T1R8 sfc deposition and wet scavenging, 
C    interpolated to L$,Z$ grid (TEMPIN)
C    interpolated to current day and combined in BOUNDC;
C
C    SCAVDAY(20,L$,Z$) is the combined surface deposition + scavenging in 1/sec
C
C   SCAVDAY:  1=ch2o  2=hno2     3=hno3   4=ho2no2   5=ho2   6=h2o2  7=mp    8=no2
C             9=no3  10=n2o5    11=o3    12=brono2  13=hbr  14=br   15=brcl 16=hobr 
C            17=clo  18=clono2  19=hcl  20=hocl
C

        COMMON/CDEPR8/DEPR8(10,14,L$), SCAVR8(20,14,L$,Z$),
     >                SCAVDAY(20,L$,Z$)



C  common for GMI SURFACE EMISSIONS, interpolated to L$ latitude grid
C    for CURRENT DAY --  EMISDAY in #/cm3/sec   (from TEMPIN):
C
C   EMISSIONS : 1=CH2O    2=CO    3=NO (NOx)
C
        COMMON/CDEP/tdep(14), DEP0(14,L$,7), EMIS0(14,L$,3), 
     >              DEPDAY(L$,7), EMISDAY(L$,3)



C  WCONV276 - extra tropospheric ozone production

        COMMON/COZPR/ozpr(360,L$,10)



        SAVE


C     SOLAR PHOTODISSOCIATION RATES FOR CHEMISTRY MODULE IN POINT CALL FORMAT.

	   JRMOM(1)=JZ(ICLOCK,1,IJC,IKC)
	   JRMOM(2)=JZ(ICLOCK,3,IJC,IKC)
	   JRMOM(3)=JZ(ICLOCK,2,IJC,IKC)
	   JRMOM(4)=0.0E0                              ! ODP
	   JRMOM(5)=JZ(ICLOCK,16,IJC,IKC)
	   JRMOM(6)=JZ(ICLOCK,7,IJC,IKC)
	   JRMOM(7)=JZ(ICLOCK,5,IJC,IKC)
	   JRMOM(8)=JZ(ICLOCK,6,IJC,IKC)
	   JRMOM(9)=JZ(ICLOCK,9,IJC,IKC)
	   JRMOM(10)=JZ(ICLOCK,15,IJC,IKC)
	   JRMOM(11)=JZ(ICLOCK,8,IJC,IKC)
	   JRMOM(12)=JZ(ICLOCK,4,IJC,IKC)
	   JRMOM(13)=JZ(ICLOCK,14,IJC,IKC)
	   JRMOM(14)=JZ(ICLOCK,18,IJC,IKC)
	   JRMOM(15)=JZ(ICLOCK,23,IJC,IKC)
	   JRMOM(16)=JZ(ICLOCK,21,IJC,IKC)
	   JRMOM(17)=JZ(ICLOCK,22,IJC,IKC)
	   JRMOM(18)=JZ(ICLOCK,19,IJC,IKC)
	   JRMOM(19)=JZ(ICLOCK,20,IJC,IKC)
	   JRMOM(20)=JZ(ICLOCK,44,IJC,IKC)
	   JRMOM(21)=JZ(ICLOCK,50,IJC,IKC)
	   JRMOM(22)=JZ(ICLOCK,62,IJC,IKC)
	   JRMOM(23)=JZ(ICLOCK,28,IJC,IKC)
	   JRMOM(24)=JZ(ICLOCK,26,IJC,IKC)

C  per Debra, HNO3 washout SHOULD NOT be in JACOB  (March 2011); 
C  it should ONLY be in SOLVER for NOy loss - set J 25 = 0. here
C
	   JRMOM(25) = 0.    ! HNO3W(ijc,ikc)

	   JRMOM(26)=SCAVDAY(6,IJC,IKC)        ! H2O2
	   JRMOM(27)=c(59,IJC,IKC)             ! HCl - not used in JACOB
	   JRMOM(28)=c(59,IJC,IKC)             ! HBr - not used in JACOB
	   JRMOM(29)=c(59,IJC,IKC)             ! HF  - not used in JACOB
	   JRMOM(30)=JZ(ICLOCK,10,IJC,IKC)
	   JRMOM(31)=JZ(ICLOCK,11,IJC,IKC)
	   JRMOM(32)=JZ(ICLOCK,42,IJC,IKC)
	   JRMOM(33)=JZ(ICLOCK,24,IJC,IKC)
	   JRMOM(34)=JZ(ICLOCK,13,IJC,IKC)
	   JRMOM(35)=SCAVDAY(7,IJC,IKC)        ! CH3OOH
	   JRMOM(36)=SCAVDAY(1,IJC,IKC)        ! CH2O
	   JRMOM(37)=JZ(ICLOCK,27,IJC,IKC)
	   JRMOM(38)=JZ(ICLOCK,38,IJC,IKC)
	   JRMOM(39)=JZ(ICLOCK,57,IJC,IKC)
	   JRMOM(40)=JZ(ICLOCK,37,IJC,IKC)
	   JRMOM(41)=JZ(ICLOCK,36,IJC,IKC)
	   JRMOM(42)=JZ(ICLOCK,29,IJC,IKC)
	   JRMOM(43)=JZ(ICLOCK,32,IJC,IKC)
	   JRMOM(44)=JZ(ICLOCK,33,IJC,IKC)
	   JRMOM(45)=JZ(ICLOCK,31,IJC,IKC)
  	   JRMOM(46)=JZ(ICLOCK,30,IJC,IKC)
	   JRMOM(47)=JZ(ICLOCK,39,IJC,IKC)
	   JRMOM(48)=0.0E0                      ! hydrocarb.
	   JRMOM(49)=0.0E0                      ! hydrocarb.
	   JRMOM(50)=0.0E0                      ! hydrocarb.
	   JRMOM(51)=0.0E0                      ! hydrocarb.
	   JRMOM(52)=0.0E0                      ! hydrocarb.
	   JRMOM(53)=0.0E0                      ! hydrocarb.
	   JRMOM(54)=0.0E0                      ! hydrocarb.
	   JRMOM(55)=0.0E0                      ! PAN
           JRMOM(56)=c(59,IJC,IKC)              ! CF2O  - not used in JACOB
	   JRMOM(57)=JZ(ICLOCK,40,IJC,IKC)
	   JRMOM(58)=JZ(ICLOCK,48,IJC,IKC)
	   JRMOM(59)=JZ(ICLOCK,17,IJC,IKC)
	   JRMOM(60)=JZ(ICLOCK,58,IJC,IKC)
	   JRMOM(61)=JZ(ICLOCK,49,IJC,IKC)
	   JRMOM(62)=JZ(ICLOCK,64,IJC,IKC)
	   JRMOM(63)=JZ(ICLOCK,34,IJC,IKC)
	   JRMOM(64)=JZ(ICLOCK,35,IJC,IKC)
	   JRMOM(65)=JZ(ICLOCK,65,IJC,IKC)
	   JRMOM(66)=JZ(ICLOCK,66,IJC,IKC)
	   JRMOM(67)=JZ(ICLOCK,67,IJC,IKC)
	   JRMOM(68)=0.0E0                      ! sulfur
	   JRMOM(69)=0.0E0                      ! sulfur
	   JRMOM(70)=0.0E0                      ! sulfur
	   JRMOM(71)=0.0E0                      ! sulfur
	   JRMOM(72)=0.0E0                      ! sulfur
	   JRMOM(73)=0.0E0                      ! sulfur
	   JRMOM(74)=0.0E0                      ! sulfur
	   JRMOM(75)=0.0E0                      ! sulfur
	   JRMOM(76)=JZ(ICLOCK,54,IJC,IKC)
	   JRMOM(77)=JZ(ICLOCK,68,IJC,IKC)

	   JRMOM(78)=JZ(ICLOCK,47,IJC,IKC)
	   JRMOM(79)=JZ(ICLOCK,25,IJC,IKC)
	   JRMOM(80)=JZ(ICLOCK,12,IJC,IKC)
	   JRMOM(81)=JZ(ICLOCK,41,IJC,IKC)
c                                           ! load JCH4 into 1 array for AER chem
           JRMOM(82)=JZ(ICLOCK,59,IJC,IKC) 
     >             + JZ(ICLOCK,60,IJC,IKC)
     >             + JZ(ICLOCK,61,IJC,IKC)
c                                           ! load JCH4 in JEX array to keep separate
	   JEX(1) = JZ(ICLOCK,59,IJC,IKC)
	   JEX(2) = JZ(ICLOCK,60,IJC,IKC)
	   JEX(3) = JZ(ICLOCK,61,IJC,IKC)


	   JRMOM(83)=0.0E0                      ! iodine
	   JRMOM(84)=0.0E0                      ! iodine
	   JRMOM(85)=0.0E0                      ! iodine
	   JRMOM(86)=0.0E0                      ! iodine
	   JRMOM(87)=0.0E0                      ! iodine
	   JRMOM(88)=0.0E0                      ! iodine
	   JRMOM(89)=0.0E0                      ! iodine
	   JRMOM(90)=0.0E0                      ! iodine
	   JRMOM(91)=0.0E0                      ! iodine
	   JRMOM(92)=0.0E0                      ! iodine
	   JRMOM(93)=0.0E0                      ! iodine
	   JRMOM(94)=0.0E0                      ! CH3CN
	   JRMOM(95)=0.0E0                      ! HCN


	   JRMOM(96)=JZ(ICLOCK,69,IJC,IKC)
	   JRMOM(97)=JZ(ICLOCK,70,IJC,IKC)
	   JRMOM(98)=JZ(ICLOCK,43,IJC,IKC)
	   JRMOM(99)=JZ(ICLOCK,51,IJC,IKC)
	   JRMOM(100)=JZ(ICLOCK,52,IJC,IKC)
	   JRMOM(101)=JZ(ICLOCK,53,IJC,IKC)
	   JRMOM(102)=JZ(ICLOCK,55,IJC,IKC)
	   JRMOM(103)=JZ(ICLOCK,56,IJC,IKC)
	   JRMOM(104)=JZ(ICLOCK,45,IJC,IKC)

C
C   updated for JPL10/9ha runs - Sept 2011 (J[O2]-> O + O1D now separate)
C
           JRMOM(105)=JZ(ICLOCK,79,IJC,IKC)
           JRMOM(106)=JZ(ICLOCK,80,IJC,IKC)
           JRMOM(107)=JZ(ICLOCK,81,IJC,IKC)
           JRMOM(108)=0.                        ! HOONO - as with HNO3, this = 0 for JACOB
       	   JRMOM(109)=SCAVDAY(5,IJC,IKC)        ! HO2
       	   JRMOM(110)=SCAVDAY(11,IJC,IKC)       ! O3

           JRMOM(111)=JZ(ICLOCK,71,IJC,IKC)     ! HFC-134a
           JRMOM(112)=JZ(ICLOCK,72,IJC,IKC)     ! HFC-143a
           JRMOM(113)=JZ(ICLOCK,73,IJC,IKC)     ! HFC-23
           JRMOM(114)=JZ(ICLOCK,74,IJC,IKC)     ! HFC-32
           JRMOM(115)=JZ(ICLOCK,75,IJC,IKC)     ! HFC-125
           JRMOM(116)=JZ(ICLOCK,76,IJC,IKC)     ! HFC-152a
           JRMOM(117)=JZ(ICLOCK,77,IJC,IKC)     ! HFC-227a
           JRMOM(118)=JZ(ICLOCK,78,IJC,IKC)     ! HFC-245fa

           JRMOM(119)=JZ(ICLOCK,46,IJC,IKC)     ! O2 -> O + O1D (now separate)


C  load in CH2O surface emissions for current day (L$ latitude grid)
C    into JRMOM(120) ; EMISDAY(L$,3) is in #/cm3/sec;   1=CH2O
C
       	   JRMOM(120)=0.0E0
           if (ikc .eq. 1) JRMOM(120) = EMISDAY(IJC,1)


C  add tropospheric ozone production into JRMOM(121)
C     ozpr(360,L$,10) in COMMON (computed in FILLA) in #/cm3/second
C
       	   JRMOM(121)=0.0E0

           IF (zalt(ikc) .le. 6.) then
            JRMOM(121) = ozpr(iday360,ijc,ikc)
            if (zalt(ikc) .ge. 4.) JRMOM(121) = .4*ozpr(iday360,ijc,ikc)
            if (zalt(ikc) .ge. 5.) JRMOM(121) = .1*ozpr(iday360,ijc,ikc)
           ENDIF


	RETURN
	END
