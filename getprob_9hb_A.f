
	SUBROUTINE GETPROB


C     THIS routine reads in TPROB45/LPROB45 arrays for proper YEAR. 
C     This routine is called ON FIRST DAY OF THE YEAR ONLY
C

       include "com2d.h"


       CHARACTER*40 CCF
       CHARACTER*4 CCYR

ccccccc       DATA CCYR/'clim', '1958', '1959', '1960', '1961', '1962', '1963', 
ccccccc     >           '1964', '1965', '1966', '1967', '1968', '1969', '1970', 
ccccccc     >           '1971', '1972', '1973', '1974', '1975', '1976', '1977', 
ccccccc     >           '1978', '1979', '1980', '1981', '1982', '1983', '1984', 
ccccccc     >           '1985', '1986', '1987', '1988', '1989', '1990', '1991', 
ccccccc     >           '1992', '1993', '1994', '1995', '1996', '1997', '1998', 
ccccccc     >           '1999', '2000', '2001', '2002', '2003', '2004', '2005',
ccccccc     >           '2006', '2007', '2008', '2009', '2010'/


	SAVE

C
C   read in and load parcel model temperature and latitude PDF arrays
C     all are on the 45x76 grid from MERRA/MLS data for 19790-2010
c     TPROB45(211,L$,Z$,72), for 120-330K is IN COMMON,  and LPROB45(181,L$,Z$,72) for 90S-90N 
C
C
C   for 9JG run, use Interannual for 1958-2010, climatology otherwise
C     now use internal write/read to convert INTEGER year to CHAR (Jan 2012)

         ccyr = 'clim'

         iyx = IYR + 1935
         IF (IYR .ge. 23  .and.  IYR .le. 75) write(ccyr, '(I4.4)' ) iyx

                                         ! set here if climatology is set for ALL YEARS in MAIN
         if (iclim .eq. 1) ccyr = 'clim'


C  read in TPROB45(125,45,76,72), LPROB45(70,45,76,72) for proper year, and the starting/ending indicies:
C      itr1(45,76,72), itr2(45,76,72), ilr1(45,76,72), ilr2(45,76,72)

        ccf = 'tprobm_45x76_'//ccyr//'.xdr'


        OPEN (75, FILE = ccf, convert="big_endian",
     >    form = 'unformatted', access = 'stream', status = 'old')
              READ (75) TPROB45
              READ (75) ITR1
              READ (75) ITR2

              READ (75) LPROB45
              READ (75) ILR1
              READ (75) ILR2
        CLOSE (75)


	RETURN
	END
