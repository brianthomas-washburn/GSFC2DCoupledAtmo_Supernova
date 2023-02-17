
       SUBROUTINE UNDUMP


       include "com2d.h"


       PARAMETER (IISP=54)

       REAL*4 ZZOUT(L$,Z$)

       REAL, allocatable :: latrr(:), zzrr(:)
       REAL, allocatable :: ccinhr(:,:,:), CCIN2(:,:)


C  indicies for mapping AER FAST species indicies into GSFC indicies
C    for the AER hydrocarbons and Iodine family that GSFC doesn't have, just set index to 0
C
      INTEGER ISPMAP0(IISP), ISPADJ(28)
      DATA ISPMAP0/2, 1, 41,  4, 12, 13, 14, 16, 22, 24, 23,  0,  0,  0, 
     >             0, 0,  0,  0,  0,  0, 27, 28, 29, 25, 30, 26, 62, 65,
     >            61,64, 63,  9,  5,  6,  7,  8, 38, 42, 10, 89, 45, 44,
     >            46,68, 47, 79, 43, 60,  0,  0,  0,  0,  0,  0/

      DATA ISPADJ/34, 35, 36, 37, 40, 49, 50, 51, 52, 53, 54, 55, 
     >            69, 70, 71, 72, 75, 81, 82, 83, 84, 85, 86, 87, 
     >            88, 18, 11, 20/



C   9HB - initialize C(S$,L$,Z$) array

            do 767 ik=1,Z$
            do 767 ij=1,L$
            do 767 js=1,S$
 767           c(js,ij,ik) = 1.e-12


c  1  O(3p)	26  CH3OOH		51  CF2ClBr
c  2  O(1d)	27  Cl			52  CHClF2
c  3  O2   	28  ClO			53  C2Cl3F3
c  4  O3   	29  HCl			54  C2Cl F4
c  5  NO   	30  ClONO2		55  C2ClF5
c  6  NO2	31  Odd N		56  HF
c  7  NO3	32  HO2NO2		57  CClFO
c  8  N2O5	33  Cly			58  CF2O
c  9  N		34  CFCl3		59  H2O*
c  10  HNO3	35  CF2Cl2		60  BrCl
c  11  N2O	36  CCl4  		61  Cl2O2
c  12  H   	37  CH3Cl		62  ClOx
c  13  OH	38  Wat			63  ClONO
c  14  HO2	39  Ox			64  Cl2
c  15  H2O	40  CH3CCl3		65  N2
c  16  H2O2	41  O2(1delta)
c  17  H2 	42  N(2D)
c  18  CH4	43  CO2
c  19  CO	44  BrO
c  20  CH3	45  Br
c  21  CH3O2	46  HBr
c  22  CH3O	47  BrONO2
c  23  H2CO	48  Brx
c  24  HCO	49  CH3Br
c  25  HOCl	50  CF3Br

c  Find out which of these numpers are really necessary
c  Get the proper common blocks
c  Create the "right" fort.4 file from an old one (i.e., write a translater)
c  The idea is to read how many species are in the file, and then to load 
c  the arrays based on the location numbers that are given above.
c  Any fort.4 file may be used to start a run, regardless of the number of
c  species.  If there are less, the species will be ignored.  If there are more,
c  the species not present in the fort.4 file will be initialized with 1.E-12 
c  number density, since zeroes may give underflows.
C
C  NEW for run12 (9HB) model - read in latitudes/altitudes with initial profiles
C    and interpolate to current model grid as necessary, so it's automated
C    in going to/from different resolutions - use ALLOCATE function
C

	READ(4,1) ISP, IJI, IKI
        print *,' isp=',isp, iji, iki

        IF (ISP .gt. S$) then
           print *, 'PROBLEM in UNDUMP!!  ISP GT S$ '
           stop
        ENDIF


        allocate(latrr(iji))
        allocate(zzrr(iki))
        allocate(ccinhr(isp,iji,iki))
        allocate(ccin2(iji,iki))


	do 300 ik=1,iki
	do 300 ij=1,iji
	do 300 js=1,isp
          if (js .ne. 59) then
		ccinhr(js,ij,ik) = 1.D-12
          endif
300	continue


 99	format(1x)
   1	FORMAT(4I4,34I2,/,30I2)

   2	FORMAT(4e15.8,I4,e15.8)           ! note change in FORMAT here from I3 to I4
   3	FORMAT(6D13.6)
   4	FORMAT(I5)
   55   FORMAT(10F7.0)
   77   FORMAT(10F7.2)


	read(4,99)
	READ(4,2)TX,TSS,DECD,DAY360,IYR,DAYST
        print *,' day360=',DAY360,' iyr=',IYR,' dayst=',DAYST
	DAY0=DAY360

	read(4,99)
        READ(4,55) (latrr(ij), ij=1,iji)

	read(4,99)
        READ(4,77) (zzrr(ik), ik=1,iki)

	read(4,99)


	DO 200 I=1,ISP

		read(4,4)JS
			DO 150 IJ=1,IJI

                           READ(4,3) (CCINHR(JS,IJ,IK),IK=1,IKI)
C                   print *,'----- in undump: ',i,ij,js,CCINHR(JS,IJ,1)

150			CONTINUE

200	CONTINUE


c 
C  IF current L$,Z$ resolution is IJIxIKI, then just load into C(S$,L$,Z$) array
C    IF NOT, then INTERPOLATE INPUT ARRAY to current grid:
C    CCINHR(ISP,IJI,IKI) -> C(S$,L$,Z$)
C    CCIN2(IJI,IKI), ZZOUT(L$,Z$), LAT4(L$), ZALT90(Z$) are both REAL*4
C

        IF (L$ .eq. IJI  .and.  Z$ .eq. IKI) THEN

            do 757 ik=1,Z$
            do 757 ij=1,L$
            do 757 js=1,ISP
 757           c(js,ij,ik) = ccinhr(js,ij,ik)

        ELSE


          DO 350 js=1,ISP
c                                ensure fields are non-zero before log interp
             DO 850 ik=1,IKI
             DO 850 ij=1,IJI
               ccin2(ij,ik) = ccinhr(js,ij,ik)
               if (ccin2(ij,ik) .lt. 1.D-12) ccin2(ij,ik) = 1.D-12
 850         CONTINUE
C                                ! do LOG INTERPOLATION in both Lat and pres/alt

          CALL BINTERP(latrr, IJI, zzrr, IKI, ccin2, 
     >                 LAT4, L$, ZALT90, Z$, 1, 1, ZZOUT)

             DO 875 ik=1,Z$
             DO 875 ij=1,L$
 875            c(js,ij,ik) = zzout(ij,ik)

350       CONTINUE

        ENDIF


        deallocate(latrr)
        deallocate(zzrr)
        deallocate(ccinhr)
        deallocate(ccin2)



C  initialize transported Cly, Bry families for extended cn58 array in MAIN
C  include BrCl in Bry (not in Cly)
C
C  also  define  M  array  here for inital call to CONTIN       1/17/95
C                          and initialize new N2(L$,Z$) array for N2

         DO 775 ik=1,Z$
         DO 775 ij=1,L$
             C(32,ij,ik) = c(25,ij,ik) + c(26,ij,ik) + c(27,ij,ik)  
     >                   + c(28,ij,ik) + c(29,ij,ik) + c(30,ij,ik) 
     >                   + 2.*c(61,ij,ik) + c(62,ij,ik) + c(63,ij,ik)
     >                   + 2.*c(64,ij,ik) + c(65,ij,ik) 

             C(58,ij,ik) = 2.*c(43,ij,ik) + c(44,ij,ik) + c(45,ij,ik)  
     >                   + c(46,ij,ik) + c(47,ij,ik) + c(60,ij,ik)
     >                   + c(68,ij,ik) + c(79,ij,ik)

             M(ij,ik)  = C(3,ij,ik)/0.21d0
             N2(ij,ik) = 0.79d0*M(ij,ik)
 775     CONTINUE


C  IYR updated Jan. 2012 - DON't use IYR from f8* initial profiles read in above
C    INSTEAD, use IYEARSS from CONTROL.dat/.f, adjust here:
C       0=1935;   15=1950;   65=2000;    2100=165
C
C    and remember that the year (IYR) from the initial profiles is from
C    the year BEFORE IYEARSS, with DAY360 = 360.5, so subtract 1 from IYEARSS

         IYR = IYEARSS - 1935 - 1


C
C  ADJUST AGE:
C
C  YEAR0 - corresponding to surface year of intial profiles - use lat index 9 for all grids
C  YEAR1 - corresponds to IYR/DAY360 read in above for current run - adjust age to match these 
C
          year0 = c(78,9,1)/M(9,1)
          year1 = IYR + 1935. + (DAY360 - 1.5)/360.

          DO 577 ik=1,Z$
          DO 577 ij=1,L$
 577         c(78,ij,ik) = c(78,ij,ik) + (year1 - year0)*M(IJ,IK)



C  Adjust initial source gases from INPUT profiles (BCC1) to value in BC table, 
C   convert to index in BC table for proper year of current run (BCC2):  
C   1=1935,  16=1950,  36=1970,   66=2000,   166=2100;  current range is 1950-2100 (16-166)
C   IYEAR1 (NINT) gets nearest whole number year;    BCTDINPUT(S$,300), ISPADJ(28)


          iyear1 = NINT(year1) - 1934
          if (iyear1 .le. 16)  iyear1 = 16
          if (iyear1 .ge. 166) iyear1 = 166


          DO 777 iss=1,28

            ispp = ISPADJ(iss)

            fac1 = 1.e9
            if (ispp .eq. 20) fac1 = 1.e6

            bcc1 = C(ispp,9,1)/M(9,1)*fac1
            bcc2 = BCTDINPUT(ispp, iyear1)
C                                           set limits where BC = 0
            bcca = 1.
            if (bcc1 .ge. 1.e-7) bcca = bcc2/bcc1


            DO 877 ik=1,Z$
            DO 877 ij=1,L$

               c(ispp,ij,ik) = c(ispp,ij,ik)*bcca
               if (c(ispp,ij,ik) .LT. 1.D-12) c(ispp,ij,ik) = 1.D-12

C   set CH2Br2, CHBr3 = 0 (dont need)

cccc          c(76,ij,ik) = 1.e-12
cccc          c(77,ij,ik) = 1.e-12

 877        CONTINUE

 777      CONTINUE


C
C
c    reset rainout to 0.0 above level 9, rainout only for levels 1-9 (for 46 level model), don't need this 

cchr      DO 752 ik=10,Z$
cchr      DO 752 ij=1,L$
cchr 752	 c(59,IJ,IK) = 0.0



C   Now initialize families for AER Model chemistry

cchr      DO 750 ik=1,Z$
cchr      DO 750 ij=1,L$
C             for new constituents (w/ no initial conditions), just set to 1.D-12 everywhere
C             before using in model: ClOO, BrNO2, Halon-1202, CH2Br2, Bromoform
cchr        c(65,ij,ik) = 1.D-12
cchr        c(75,ij,ik) = 1.D-12
cchr        c(76,ij,ik) = 1.D-12
cchr        c(77,ij,ik) = 1.D-12
cchr        c(79,ij,ik) = 1.D-12


c   initialize total NOy field, CN(31);  set c(32) = 0.0  (not used)

cchr        c(31,ij,ik) = c(9,ij,ik) + c(5,ij,ik) + c(6,ij,ik) + c(7,ij,ik)
cchr     >        + 2.*c(8,ij,ik) + c(38,ij,ik) + c(42,ij,ik)  + c(10,ij,ik)
cchr     >           + c(30,ij,ik) + c(63,ij,ik) + c(47,ij,ik) + c(79,ij,ik) 
cchr     >           + c(66,ij,ik)
cchr
cchr        c(32,ij,ik) = 0.0


C  and initialize CHx, c(73) = CH3O2+CH2O+CH3OOH, and HOx c(74) = H + OH + HO2 + 2*H2O2

cchr        c(73,ij,ik) = c(22,ij,ik) + c(23,ij,ik) + c(24,ij,ik) 
cchr
cchr        c(74,ij,ik) = c(12,ij,ik) + c(13,ij,ik) + c(14,ij,ik)
cchr     >              + 2.*c(16,ij,ik) 


C  Ox        

cchr        c(39,ij,ik) = c(1,ij,ik) + c(2,ij,ik) + c(4,ij,ik)



C  Bry -  the individual members got messed up in the f8 file, so need to rescale them here
C           based on the initial Bry, which is kept

cchr        brsum = 2.*c(43,ij,ik) + c(44,ij,ik) + c(45,ij,ik) + c(79,ij,ik)
cchr     >   + c(46,ij,ik) + c(47,ij,ik) + c(60,ij,ik) + c(68,ij,ik)

cchr        c(43,ij,ik) = c(43,ij,ik)/brsum*c(48,ij,ik)
cchr        c(44,ij,ik) = c(44,ij,ik)/brsum*c(48,ij,ik)
cchr        c(45,ij,ik) = c(45,ij,ik)/brsum*c(48,ij,ik)
cchr        c(46,ij,ik) = c(46,ij,ik)/brsum*c(48,ij,ik)
cchr        c(47,ij,ik) = c(47,ij,ik)/brsum*c(48,ij,ik)
cchr        c(68,ij,ik) = c(68,ij,ik)/brsum*c(48,ij,ik)
cchr        c(60,ij,ik) = c(60,ij,ik)/brsum*c(48,ij,ik)
cchr        c(79,ij,ik) = c(79,ij,ik)/brsum*c(48,ij,ik)
cchr


C  Cly - total Chlorine - C(33)

cchr        c(33,ij,ik) = c(25,ij,ik) + c(26,ij,ik) + c(27,ij,ik)
cchr     >              + c(28,ij,ik) + c(29,ij,ik) + c(30,ij,ik)
cchr     >              + c(60,ij,ik) + 2.*c(61,ij,ik) + c(62,ij,ik)
cchr     >              + c(63,ij,ik) + 2.*c(64,ij,ik) + c(65,ij,ik)

cchr 750	CONTINUE


C
C  And initialize the CN array here, also initialize 
C       NOONRAT(S$,L$,Z$)=1.: noontime radical/24-hr avg radical
C
C       also initialize CDN(2,S$,L$,Z$) - output array of DAY/NITE avgs
C
	do 450 ik=1,Z$
	do 450 ij=1,L$
	do 450 is=1,S$
		cn(is,ij,ik) = c(is,ij,ik)
                noonrat(is,ij,ik) = 1.d0
                cdn(1,is,ij,ik) = 0.0
                cdn(2,is,ij,ik) = 0.0
450	continue



C  for diurnally varying columns, need to initialize CNDC(S$,ntime+3,L$,Z$) (in COMMON)
C     just load in initial (diurnally avged) profiles

        DO 100 IS=1,IISP
          if (ISPMAP0(is) .ne. 0) then
            DO 105 IK=1,Z$
            DO 105 IJ=1,L$
            DO 105 itt=1,18+3
 105           CNDC(ISPMAP0(is), itt, ij,ik) = CN(ISPMAP0(is),ij,ik)
          endif
 100   CONTINUE


Cef
C  load in ratio of starting radical concentrations to the total family 
C     (which has no diurnal cycle)
C      NOONRAT(S$,L$,Z$), this will just be converted back to the initial value in FILLA
C        (except for O2(1D), which is just started from the initial value)
Cef
Cef  - NO! No longer do this when transporting radicals (everything is transported), 
Cef      for initial time step, just set NOONRAT = 1 for everything (set above)
Cef
C
Cef         DO 700 ik=1,Z$
Cef         DO 700 ij=1,L$
C  Ox
Cef           NOONRAT(1,ij,ik) = cn(1,ij,ik)/cn(39,ij,ik)
Cef           NOONRAT(2,ij,ik) = cn(2,ij,ik)/cn(39,ij,ik)

Cef           NOONRAT(4,ij,ik) = cn(4,ij,ik)/cn(39,ij,ik) 
C      !  need Ozone for the mesosphere in OX
C
Cef           NOONRAT(41,ij,ik) = 1.     !  to start, just set Noon ratio = 1. for O2(1D)
Cef
C  HOx
Cef           NOONRAT(12,ij,ik) = cn(12,ij,ik)/cn(74,ij,ik)
Cef           NOONRAT(13,ij,ik) = cn(13,ij,ik)/cn(74,ij,ik)
Cef           NOONRAT(14,ij,ik) = cn(14,ij,ik)/cn(74,ij,ik)
Cef           NOONRAT(16,ij,ik) = cn(16,ij,ik)/cn(74,ij,ik)
Cef
C  CHx
Cef           NOONRAT(22,ij,ik) = cn(22,ij,ik)/cn(73,ij,ik)
Cef           NOONRAT(23,ij,ik) = cn(23,ij,ik)/cn(73,ij,ik)
Cef           NOONRAT(24,ij,ik) = cn(24,ij,ik)/cn(73,ij,ik)
Cef
C  NOz
Cef           NOONRAT(5,ij,ik) = cn(5,ij,ik)/cn(32,ij,ik)
Cef           NOONRAT(6,ij,ik) = cn(6,ij,ik)/cn(32,ij,ik)
Cef           NOONRAT(7,ij,ik) = cn(7,ij,ik)/cn(32,ij,ik)
Cef           NOONRAT(8,ij,ik) = cn(8,ij,ik)/cn(32,ij,ik)
Cef           NOONRAT(9,ij,ik) = cn(9,ij,ik)/cn(32,ij,ik)
Cef           NOONRAT(38,ij,ik) = cn(38,ij,ik)/cn(32,ij,ik)
Cef           NOONRAT(42,ij,ik) = cn(42,ij,ik)/cn(32,ij,ik)
Cef
C  Clx
Cef           NOONRAT(25,ij,ik) = cn(25,ij,ik)/cn(33,ij,ik)
Cef           NOONRAT(26,ij,ik) = cn(26,ij,ik)/cn(33,ij,ik)
Cef           NOONRAT(27,ij,ik) = cn(27,ij,ik)/cn(33,ij,ik)
Cef           NOONRAT(28,ij,ik) = cn(28,ij,ik)/cn(33,ij,ik)
Cef           NOONRAT(29,ij,ik) = cn(29,ij,ik)/cn(33,ij,ik)
Cef           NOONRAT(30,ij,ik) = cn(30,ij,ik)/cn(33,ij,ik)
Cef           NOONRAT(61,ij,ik) = cn(61,ij,ik)/cn(33,ij,ik)
Cef           NOONRAT(62,ij,ik) = cn(62,ij,ik)/cn(33,ij,ik)
Cef           NOONRAT(63,ij,ik) = cn(63,ij,ik)/cn(33,ij,ik)
Cef           NOONRAT(64,ij,ik) = cn(64,ij,ik)/cn(33,ij,ik)
Cef
C  Brx
Cef           NOONRAT(43,ij,ik) = cn(43,ij,ik)/cn(48,ij,ik)
Cef           NOONRAT(44,ij,ik) = cn(44,ij,ik)/cn(48,ij,ik)
Cef           NOONRAT(45,ij,ik) = cn(45,ij,ik)/cn(48,ij,ik)
Cef           NOONRAT(46,ij,ik) = cn(46,ij,ik)/cn(48,ij,ik)
Cef           NOONRAT(47,ij,ik) = cn(47,ij,ik)/cn(48,ij,ik)
Cef           NOONRAT(60,ij,ik) = cn(60,ij,ik)/cn(48,ij,ik)
Cef           NOONRAT(68,ij,ik) = cn(68,ij,ik)/cn(48,ij,ik)
Cef
Cef 700     CONTINUE
Cef

	RETURN
	END

