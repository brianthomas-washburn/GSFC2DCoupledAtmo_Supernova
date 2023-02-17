       SUBROUTINE PHOTIN

C   Received from Randy Kawa 1/24/95 - adapted from ptest_2d.f
C	use this with photcodes7.f subroutines
C
C  - updated Jan. 2002 to have 100 pressure levels for new resolution model

C   rlam -- array of wavelength values (nlam) used to generate s table
C   sza_tab -- array of sza values (nsza) used to generate s table
C   o3_tabin -- array of col o3 values (no3,npr) used to generate s table
C   stabin -- array of s data (nsza,no3s,nprs,nwavelengths)
C
C   July 2011 - Randy's arrays with 100 pressure levels are now interpolated in LOG
C         (in the vertical) to the 2D model pressure levels of the current Z$ grid
C

        include "com2d.h"
        include "comphot.h"

        INTEGER lun3,lun4

        REAL o3_tabin(no3,npr), stabin(nsza,no3,npr,nlam)
        REAL o2jdatin(nsza,no3,npr), ztab(npr), tt7(npr), uu7(nz)


        write (6,*) ' enter photin program ' 

c        Read in cross-section, quantum yield, and solar flux data. 

         lun3=41
         lun4=42 
       
	OPEN (UNIT=lun3,FILE='stab100_2d69.xdr', convert="big_endian", 
     1		FORM='unformatted',access='stream')
	OPEN (UNIT=lun4,FILE='indx100_2d69.dat',
     1		FORM='formatted')

C PURPOSE: read in radiation source function data table for subsequent
C		interpolation.  Also read table indexing arrays.
C CATEGORY:
C INPUTS:
C   lun4 -- unit number of file containing table index values, ie., pressures,
C		wavelengths, sza, and o3 cols (no3,npr) used to generate stabin
C		Values were written at same time as stabin.
C   lun3 -- unit number of file containing s data table.  
C		Data must be on wavelength grid consistent with x-sections.
C		Table is written by wrt_stab.pro
C   nlam -- number of wavelength intervals used (checked against read-in)
C   nsza -- number of zenith angles in table
C   no3s -- number of column o3 values in table
C OPTIONAL INPUT PARAMETERS:
C KEYWORD PARAMETERS:
C OUTPUTS:
C   rlam -- array of wavelength values (nlam) used to generate s table
C   sza_tab -- array of sza values (nsza) used to generate s table
C   o3_tabin -- array of col o3 values (no3,npr) used to generate s table
C   stabin -- array of s data (nzens,no3s,nprs,nwavelengths)
C   o2jdatin -- array (nzens,no3s,nprs) of J(O2) values from srcfun output
C COMMON BLOCKS:
C SIDE EFFECTS:
C PROCEDURE:
C   Straightforward.
C RESTRICTIONS:
C   
C REQUIRED ROUTINES: Open (FORM='system') and Close units in calling
C			 program
C MODIFICATION HISTORY: 
C   Created 930824 - SR Kawa
C   Modified 950404 to read in J(O2) values from srcfun appended to s table
C   $Header: /home/kawa/3d/rdstab.f $
C-
C   Modified 020115 to have 100 pressure levels for new high res model
C
C40


C Read in dimensions of index arrays and check for consistency
	READ (lun4,*) npr_in,nlam_in,nsza_in,no3_in
	PRINT *,'s data array dimensions:',nsza_in,no3_in,npr_in,nlam_in
	write (6,*) 
        IF((npr_in.NE.npr).OR.(nlam_in.NE.nlam).OR.(nsza_in.NE.nsza).OR.
     1		(no3_in.NE.no3)) THEN 
	   PRINT *,'Array dimensions of table indices do not match',
     1		' those expected',npr_in,npr,nlam_in,nlam,nsza_in,
     2		nsza,no3_in,no3
	   RETURN
	ENDIF

C  Read in indexing values
	READ (lun4,*) pr_tab
        print *,' right before printout of pr_tab'
        print *,pr_tab
        print *,' right after printout of pr_tab'
	READ (lun4,*) rlam
        print *,rlam
        print *,' right after printout of rlam'
	READ (lun4,*) sza_tab
        print *,sza_tab
        print *,' right after printout of sza_tab'
	READ (lun4,*) o3_tabin
c        print *,o3_tabin

C Read in s function data in binary format.  
	READ (lun3) stabin
	READ (lun3) o2jdatin


        do 151 ik=1,npr
 151       ztab(ik) = 7.*ALOG(1013./pr_tab(ik))


C
C  interpolate arrays to current 2D model pressure grid (nz=Z$) for use in RFLUXINTERP
C    use ZALT90(Z$), ztab(npr), INTERPOLATE LOGARITHMICALLY (no 0s in the input files)
C       tt7(npr), uu7(nz)
C
C       o3_tabin(no3,npr)         -> o3_tab(no3,nz)
C       stabin(nsza,no3,npr,nlam) -> stab(nsza,no3,nz,nlam)
C       o2jdatin(nsza,no3,npr)    -> o2jdat(nsza,no3,nz)
C

        do 100 ij=1,no3

           do 101 ik=1,npr
 101          tt7(ik) = o3_tabin(ij,ik)

           CALL LINTERP(ztab, npr, tt7, zalt90, nz, 1, uu7)

           do 102 ik=1,nz
 102          o3_tab(ij,ik) = uu7(ik)

 100    CONTINUE


        do 200 im7=1,nlam
        do 200 ij=1,no3
        do 200 iii=1,nsza

           do 201 ik=1,npr
 201          tt7(ik) = stabin(iii,ij,ik,im7)

           CALL LINTERP(ztab, npr, tt7, zalt90, nz, 1, uu7)

           do 202 ik=1,nz
 202          stab(iii,ij,ik,im7) = uu7(ik)

 200    CONTINUE


        do 300 ij=1,no3
        do 300 iii=1,nsza

           do 301 ik=1,npr
 301          tt7(ik) = o2jdatin(iii,ij,ik)

           CALL LINTERP(ztab, npr, tt7, zalt90, nz, 1, uu7)

           do 302 ik=1,nz
 302          o2jdat(iii,ij,ik) = uu7(ik)

 300    CONTINUE


	RETURN
	END
