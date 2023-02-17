C
       SUBROUTINE DYNOUT(im10, iiyr, im324)


       include "com2d.h"

C
       CHARACTER*45 CCF
       CHARACTER*4 CCYR(33)

       DATA CCYR/'junk', '1979', '1980', '1981', '1982', '1983', '1984', 
     >           '1985', '1986', '1987', '1988', '1989', '1990', '1991', 
     >           '1992', '1993', '1994', '1995', '1996', '1997', '1998', 
     >           '1999', '2000', '2001', '2002', '2003', '2004', '2005',
     >           '2006', '2007', '2008', '2009', '2010'/


c
        SAVE

C
C  now store 1st and last time periods in extra arrays:  
C     W74(L$,Z$X1,2), KZZ74(L$,Z$X1,2), KYY74(L$+1,Z$X,2), TEMP74(L$,Z$X,2)  are in COMMON
C     W11(L$,Z$X1,2), KZZ11(L$,Z$X1,2), KYY11(L$+1,Z$X),2, TEMP11(L$,Z$X,2)  are in COMMON
C     WALL(L$,Z$X1,74), KZZTOT(L$,Z$X1,74), KYYALL(L$+1,Z$X,74), TEMPALL(L$,Z$X,74) are in COMMON
C
C   now at the end of the year (im10=72) store in w74, kzz74, etc....
C        ie, load previous years' im10=72 into w74(1), and load in current years' im10=72 into w74(2)
C        also do special cases (initiailize) for 1st and last years

      IF (IM10 .EQ. 72) THEN

         if (iiyr .eq. 1) then
            do 601 ik=1,Z$X1
            do 601 ij=1,L$
              w74(ij,ik,2) = wall(ij,ik,2)
 601	      kzz74(ij,ik,2) = kzztot(ij,ik,2)

            do 612 ik=1,Z$X
            do 612 ij=1,L$+1
 612	       kyy74(ij,ik,2) = kyyall(ij,ik,2)

            do 613 ik=1,Z$X
            do 613 ij=1,L$
 613	       temp74(ij,ik,2) = tempall(ij,ik,2)
         endif


         do 401 ik=1,Z$X1
         do 401 ij=1,L$
           w74(ij,ik,1) = w74(ij,ik,2)
           kzz74(ij,ik,1) = kzz74(ij,ik,2)

           w74(ij,ik,2) = wall(ij,ik,73)
           kzz74(ij,ik,2) = kzztot(ij,ik,73)
  401    CONTINUE

         do 412 ik=1,Z$X
         do 412 ij=1,L$+1
           kyy74(ij,ik,1) = kyy74(ij,ik,2)
           kyy74(ij,ik,2) = kyyall(ij,ik,73)
 412	 CONTINUE

         do 413 ik=1,Z$X
         do 413 ij=1,L$
           temp74(ij,ik,1) = temp74(ij,ik,2)
           temp74(ij,ik,2) = tempall(ij,ik,73)
 413	 CONTINUE
       ENDIF


C    and do a similar thing for im10=1 in w11(1,2), also load in wall(1) from period 73 of 2 years ago w74(1)
c      also do special cases (initiailize) for 1st and last years
C      then write out array for 74 time periods
                   
      IF (IM10 .EQ. 1) THEN

         if (iiyr .eq. 1) then
            do 1601 ik=1,Z$X1
            do 1601 ij=1,L$
              w11(ij,ik,1) = wall(ij,ik,2)
              w74(ij,ik,1) = wall(ij,ik,2)

     	      kzz11(ij,ik,1) = kzztot(ij,ik,2)
     	      kzz74(ij,ik,1) = kzztot(ij,ik,2)
 1601	    CONTINUE


            do 1612 ik=1,Z$X
            do 1612 ij=1,L$+1
               kyy11(ij,ik,1) = kyyall(ij,ik,2)
 1612	       kyy74(ij,ik,1) = kyyall(ij,ik,2)

            do 1613 ik=1,Z$X
            do 1613 ij=1,L$
               temp11(ij,ik,1) = tempall(ij,ik,2)
 1613	       temp74(ij,ik,1) = tempall(ij,ik,2)


C  also initialize clim arrays first time through here
C    TEMPCL(L$,Z$X,74), KYYCL(L$+1,Z$X,74), WCL(L$,Z$X1,74), KZZCL(L$,Z$X1,74)  are all in COMMON

            do 1900 imm=1,74
               do 1901 ik=1,Z$X
               do 1901 ij=1,L$
 1901             tempcl(ij,ik,imm) = 0.0

               do 1902 ik=1,Z$X
               do 1902 ij=1,L$+1
 1902             kyycl(ij,ik,imm) = 0.0

               do 1903 ik=1,Z$X1
               do 1903 ij=1,L$
                  wcl(ij,ik,imm) = 0.0
 1903             kzzcl(ij,ik,imm) = 0.0
 1900  CONTINUE

         endif


         do 201 ik=1,Z$X1
         do 201 ij=1,L$
           w11(ij,ik,2) = wall(ij,ik,2)

           wall(ij,ik,2) = w11(ij,ik,1)
           wall(ij,ik,74) = w11(ij,ik,2)
           wall(ij,ik,1) = w74(ij,ik,1)


           kzz11(ij,ik,2) = kzztot(ij,ik,2)

           kzztot(ij,ik,2) = kzz11(ij,ik,1)
           kzztot(ij,ik,74) = kzz11(ij,ik,2)
           kzztot(ij,ik,1) = kzz74(ij,ik,1)
  201    CONTINUE

         do 212 ik=1,Z$X
         do 212 ij=1,L$+1
           kyy11(ij,ik,2) = kyyall(ij,ik,2)

           kyyall(ij,ik,2) = kyy11(ij,ik,1)
           kyyall(ij,ik,74) = kyy11(ij,ik,2)
           kyyall(ij,ik,1) = kyy74(ij,ik,1)
 212	 CONTINUE

         do 213 ik=1,Z$X
         do 213 ij=1,L$
           temp11(ij,ik,2) = tempall(ij,ik,2)

           tempall(ij,ik,2) = temp11(ij,ik,1)
           tempall(ij,ik,74) = temp11(ij,ik,2)
           tempall(ij,ik,1) = temp74(ij,ik,1)
 213	 CONTINUE



C   WRITE OUT the transport arrays for each year (it's actually for the previous year)
C
      if (iiyr .ge. 2) then

       ccf = 'dyn_9gtp_'//ccyr(iiyr)//'.xdr'

       OPEN (463, FILE = ccf, convert="big_endian",
     >   form = 'unformatted', access='sequential', status = 'new')
            write (463) wall
            write (463) kyyall
            write (463) kzztot
            write (463) tempall
       CLOSE (463)



C  ALSO, load in for climatologies here (then write out on last year) - arrays were initialized above
C   WALL(L$,Z$X1,74), KYYALL(L$+1,Z$X,74), KZZTOT(L$,Z$X1,74), TEMPALL(L$,Z$X,74)
C    TEMPCL(L$,Z$X,74), KYYCL(L$+1,Z$X,74), WCL(L$,Z$X1,74), KZZCL(L$,Z$X1,74)  are all in COMMON
C
      do 2900 imm=1,74
         do 2940 ik=1,Z$X1
         do 2940 ij=1,L$
           wcl(ij,ik,imm) = wcl(ij,ik,imm) + WALL(ij,ik,imm)/32.
 2940	   kzzcl(ij,ik,imm) = kzzcl(ij,ik,imm) + KZZTOT(ij,ik,imm)/32.

         do 2930 ik=1,Z$X
         do 2930 ij=1,L$+1
 2930      kyycl(ij,ik,imm) = kyycl(ij,ik,imm) + KYYALL(ij,ik,imm)/32.

         do 2920 ik=1,Z$X
         do 2920 ij=1,L$
 2920     tempcl(ij,ik,imm) = tempcl(ij,ik,imm) + TEMPALL(ij,ik,imm)/32.

 2900 CONTINUE

      endif


C  after writing out, need to restore/update to current year:

         do 301 ik=1,Z$X1
         do 301 ij=1,L$
            wall(ij,ik,2) = w11(ij,ik,2)
            w11(ij,ik,1) = w11(ij,ik,2)

            kzztot(ij,ik,2) = kzz11(ij,ik,2)
            kzz11(ij,ik,1) = kzz11(ij,ik,2)
  301    CONTINUE

         do 312 ik=1,Z$X
         do 312 ij=1,L$+1
           kyyall(ij,ik,2) = kyy11(ij,ik,2) 
 312       kyy11(ij,ik,1) = kyy11(ij,ik,2) 

         do 313 ik=1,Z$X
         do 313 ij=1,L$
           tempall(ij,ik,2) = temp11(ij,ik,2)
 313       temp11(ij,ik,1) = temp11(ij,ik,2)

       ENDIF


C  finally, do special case for final year:
C      set period 1 = Dec of previous year as usual,  and set period 74 = 73

       IF (iiyr .eq. 32  .and.  im10 .eq. 72) THEN
ccccctest9gt      IF (iiyr .eq. 2  .and.  im10 .eq. 72) THEN

         do 1401 ik=1,Z$X1
         do 1401 ij=1,L$
           wall(ij,ik,1) = w74(ij,ik,1)
           wall(ij,ik,74) = wall(ij,ik,73)

           kzztot(ij,ik,1) = kzz74(ij,ik,1)
           kzztot(ij,ik,74) = kzztot(ij,ik,73)
 1401    CONTINUE

         do 1412 ik=1,Z$X
         do 1412 ij=1,L$+1
           kyyall(ij,ik,1) = kyy74(ij,ik,1)
 1412	   kyyall(ij,ik,74) = kyyall(ij,ik,73)

         do 1413 ik=1,Z$X
         do 1413 ij=1,L$
           tempall(ij,ik,1) = temp74(ij,ik,1)
 1413	   tempall(ij,ik,74) = tempall(ij,ik,73)


      ccf= 'dyn_9gtp_'//ccyr(iiyr+1)//'.xdr'

       OPEN (463, FILE = ccf, convert="big_endian",
     >   form = 'unformatted', access='sequential', status = 'new')
            write (463) wall
            write (463) kyyall
            write (463) kzztot
            write (463) tempall
       CLOSE (463)



C   also, load in final year for climatology, and write out
C      TEMPCL(L$,Z$X,74), KYYCL(L$+1,Z$X,74), WCL(L$,Z$X1,74), KZZCL(L$,Z$X1,74)  are all in COMMON
C
      do 3900 imm=1,74
         do 3940 ik=1,Z$X1
         do 3940 ij=1,L$
           wcl(ij,ik,imm) = wcl(ij,ik,imm) + WALL(ij,ik,imm)/32.
 3940	   kzzcl(ij,ik,imm) = kzzcl(ij,ik,imm) + KZZTOT(ij,ik,imm)/32.

         do 3930 ik=1,Z$X
         do 3930 ij=1,L$+1
 3930      kyycl(ij,ik,imm) = kyycl(ij,ik,imm) + KYYALL(ij,ik,imm)/32.

         do 3920 ik=1,Z$X
         do 3920 ij=1,L$
 3920     tempcl(ij,ik,imm) = tempcl(ij,ik,imm) + TEMPALL(ij,ik,imm)/32.
 3900 CONTINUE


C  for endpoints (im=1 and 74), it's best to just reload assuming wrap-around, 
C     otherwise there's a slight averaging difference for endpoints
C
         do 4940 ik=1,Z$X1
         do 4940 ij=1,L$
           wcl(ij,ik,1) = wcl(ij,ik,73)
     	   kzzcl(ij,ik,1) = kzzcl(ij,ik,73)

           wcl(ij,ik,74) = wcl(ij,ik,2)
     	   kzzcl(ij,ik,74) = kzzcl(ij,ik,2)
 4940	 CONTINUE

         do 4930 ik=1,Z$X
         do 4930 ij=1,L$+1
           kyycl(ij,ik,1) = kyycl(ij,ik,73)
 4930      kyycl(ij,ik,74) = kyycl(ij,ik,2)

         do 4920 ik=1,Z$X
         do 4920 ij=1,L$
          tempcl(ij,ik,1) = tempcl(ij,ik,73)
 4920     tempcl(ij,ik,74) = tempcl(ij,ik,2)



      OPEN (463, FILE = 'dyn_9gtp_clim.xdr',
     >    convert="big_endian", form='unformatted', 
     >    access='sequential', status = 'new')
            write (463) wcl
            write (463) kyycl
            write (463) kzzcl
            write (463) tempcl
       CLOSE (463)

C                            end iiyr=32, im10=72 loop
      ENDIF


      RETURN
      END

