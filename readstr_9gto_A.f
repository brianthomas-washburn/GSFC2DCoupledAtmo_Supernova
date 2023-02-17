
        SUBROUTINE READSTR(IIYR) 


C  read in yearly files of heating rates and temperatures/Kyys for use in STREAMFUNCTION


       include "com2d.h"


       CHARACTER*47 CCF
       CHARACTER*51 CHH
       CHARACTER*42 CCFEGW

       CHARACTER*4 CCYR(32)

       DATA CCYR/'1979', '1980', '1981', '1982', '1983', '1984', '1985',
     >           '1986', '1987', '1988', '1989', '1990', '1991', '1992',
     >           '1993', '1994', '1995', '1996', '1997', '1998', '1999',
     >           '2000', '2001', '2002', '2003', '2004', '2005', '2006',
     >           '2007', '2008', '2009', '2010'/

	SAVE



C  
C  *********    Read  in  yearly lookup  tables from MERRA/MLS  for 1979 - 2010  *********
C
C   KYYIN, EPIN, TEMP5IN, UBARIN, EHFIN, EHFZIN, QBARY(L$S=91,Z$S=117,72) in COMMON 
C

       ccf='kyy_merra_'//ccyr(iiyr)//'.xdr'

        OPEN (75, FILE = ccf, convert="big_endian",
     >    form = 'unformatted', access = 'stream', status = 'old')
             READ (75) kyyin
             READ (75) epin
             READ (75) temp5in
             READ (75) ubarin
             READ (75) ehfin
             READ (75) ehfzin
             READ (75) qbary
        CLOSE (75)


c
c  read in MERRA Heating rates,  HEATIN(L$S=91,Z$S=117,72) (in K/day) in COMMON
C      

        chh = 
     > 'qnet_merra_sc_'//ccyr(iiyr)//'.xdr'

        OPEN (76, FILE = chh, convert="big_endian",
     >    form = 'unformatted', access='stream', status = 'old')
            READ (76) HEATIN
        CLOSE (76)


c
c  read in equatorial gravity wave momentum source (in m/sec^2)
c    FEGWIN(L$S=91,Z$S=117,72) is in COMMON - read in  for each year:
C

       ccfegw = 'fegw_'//ccyr(iiyr)//'.xdr'

        OPEN (76, FILE = ccfegw, convert="little_endian",
     >    form = 'unformatted', access='stream', status = 'old')
            READ (76) FEGWIN
        CLOSE (76)


       RETURN
       END
