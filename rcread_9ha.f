c***********************************************************************
c                          subroutine reacread                         *
c this subroutine reads the files reac_base.dat, reac_work.dat, and    *
c reac_het.dat.  The files are assumed to be copied into fort.29,      *
c fort.30, and fort.31                                                 *
c***********************************************************************

      SUBROUTINE RCREAD


      include "com2d.h"

      save 



C   first read AER-GSFC gas phase reaction rate mapping - IAGMAP(382) is in COMMON
C

        OPEN (unit=802, 
     >  NAME='/misc/kah07/fleming/9h0/aer-gsfc_reac_9ha.lis',
     >		type='OLD', form='FORMATTED')

           READ (802,103)
           READ (802,103)
           READ (802,103)
           READ (802,103)
           READ (802,103)
           READ (802,103)

           do 200 idum=1,382 
              READ (802, 204) iak, igk
    	      IAGMAP(iak) = igk
 200	   CONTINUE

 204  format (4x, I3, 28x, I4) 

       CLOSE(unit=802)



c	OPEN(UNIT=29, NAME='reac_base.dat',TYPE='OLD',
c     c		FORM='FORMATTED')
c	OPEN(UNIT=30, NAME='reac_work.dat',TYPE='OLD',
c     c		FORM='FORMATTED')
c	OPEN(UNIT=31, NAME='reac_het.dat',TYPE='OLD',
c     c		FORM='FORMATTED')

c read reaction rate date from file:
c skip title and ruler:
      
      read(29,103)
      read(29,103)

c read the base set reaction rate values:

      do 501 i=1,rb$

         read(29,101) ibdy(i),k0(i),e(i),khi(i),ehi(i)

c hno3 gets special treatment: assume hno3 is in place 37 (6 parameters):

         if (i.eq.37) then
            khno3(1)=k0(i)
            ehno3(1)=e(i)
            read(29,102) khno3(2),ehno3(2)
            read(29,102) khno3(3),ehno3(3)
         endif
 501  continue

c read the working set reaction rate values:

      read(30,103)
      read(30,103)

      do 502 i=1,rw$

         read(30,101) ibdyw(i),k0w(i),ew(i),khiw(i),ehiw(i)

 502  continue


c read the heterogeneous reaction rate gammas: fort.31

      read(31,103)
      read(31,103)
      read(31,103)
      read(31,103)
      read(31,103)
      read(31,103)
      read(31,103)

      do 503 i=1,rh$
         read(31,104) gs(i)
 503  continue

      read(31,103)

      do 504 i=1,rh$
         read(31,104) gn(i)
 504  continue

      read(31,103)

      do 505 i=1,rh$
         read(31,104) gi(i)
 505  continue


 101  format (i1,34x,e9.3,f9.1,4x,e7.1,f9.1)
 102  format (37x,e7.1,f9.1)
 103  format (80x)
 104  format (35x,e7.1)

c	CLOSE(UNIT=29)
c	CLOSE(UNIT=30)
c	CLOSE(UNIT=31)

      return
      end
