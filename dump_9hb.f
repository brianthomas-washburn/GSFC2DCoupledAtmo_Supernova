C
	SUBROUTINE DUMP

C
C   THIS now writes out to a unit number based on IYR (for long runs), so we don't have 
C     to start from the very beginning if something crashes.....
C
C
C   THIS is UPDATED TO WRITE OUT L$, Z$, LAT(L$), and ZALT90(Z$) arrays
C

        include "com2d.h"


	DIMENSION CNS(S$,L$,Z$)

        SAVE

c  find out if there are any other needed commons
c  are there any variables here we don't need?
c
	DO 100 JS=1,S$
	   DO 100 IJC=1,L$
	      DO 100 IKC=1,Z$
		 CNS(JS,IJC,IKC)=0.0E0
 100		 CONTINUE

	DO 200 JS=1,S$
	   DO 200 IJC=1,L$
	      DO 200 IKC=1,Z$
		 CNS(JS,IJC,IKC) = CN(JS,IJC,IKC)
c                                                    ! set to min of 1.d-12, other than rainout, to avoid 0.0
          if (JS .ne. 58  .and.  JS .ne. 59) then
             IF (CNS(JS,IJC,IKC) .LT. 1.D-12) CNS(JS,IJC,IKC) = 1.D-12
          endif
 200		 CONTINUE


c   set ISP = S$ (now for 80 constituents w/ AER fast chemistry) to dump out profiles
C
        ISP = S$
  
        print *,' isp in dump=',isp

C                                                               ! need to use OPEN and CLOSE statements here
        OPEN (iyr+1935, form='formatted', status='unknown')


           WRITE(iyr+1935, 1) ISP, L$, Z$
           WRITE(iyr+1935, 99)
 99           FORMAT(1X)
           tss=0.0e0

           WRITE(iyr+1935, 2) T,TSS,DECD,DAY360,IYR,DYT360
           WRITE(iyr+1935, 99)

           WRITE(iyr+1935, 55) LAT4
           WRITE(iyr+1935, 99)

           WRITE(iyr+1935, 77) ZALT90
           WRITE(iyr+1935, 99)


c  ISP is the number of constituents in the run

	DO 10 JS=1,ISP

           WRITE(iyr+1935, 4) JS

		DO 10 IJ=1,L$

			WRITE(iyr+1935, 3) (CNS(JS,IJ,IK),IK=1,Z$)

10    CONTINUE

        CLOSE (iyr+1935)


    1	FORMAT(4I4,34I2,/,30I2)
    2	FORMAT(4D15.8,I4,D15.8)
    3	FORMAT(6D13.6)
    4	FORMAT(I5,10X,A8)
   55   FORMAT(10F7.0)
   77   FORMAT(10F7.2)

	WRITE(6,7171)DYT360
7171       FORMAT(' WRITE OUT DUMP AT DAY TOTAL=',1Pe15.8)
	RETURN
	END
