C  PUT IN NOX SOURCE from other ionization (in addition to GCRs)
C  AYelland 12September2020

      subroutine ionSourceIn

      include "com2d.h"
      INTEGER, parameter :: ZION = 46
      REAL*4 :: SNCRions_in(1:ZION), SNCRions_out(1:Z$)
C     REAL*4 :: ionSource(1:L$,1:Z$)    ***DECLARED IN "com2d_9hb_ions.h"***

C  Printing out subrouting readin data

C       print *, 'L$ = ', L$
C       print *, 'LAT4 = '
C          do IJ = 1,L$
C              print *, IJ, LAT4(IJ)
C          end do
C       print *,

C       print *, 'ZION = ', ZION
C       print *, 'zz46 = '
C          do IK = 1,ZION
C              print *, IK, zz46(IK)
C          end do
C       print *,

C       print *, 'Z$ = ', Z$
C       print *, 'ZALT90 = '
C          do IK = 1,Z$
C              print *, IK, ZALT90(IK)
C          end do
C       print *,

C     Intitalizing Arrays with zeros.

      do IK=1,ZION
         SNCRions_in(IK)=0.d0
      end do

      do IK=1,Z$
         SNCRions_out(IK) = 0.d0
      end do

      do IK=1,Z$
         do IJ=1,L$
            ionSource(IJ,IK)=0.d0
         end do
      end do

C  Read-in of Ionization values from Flux & Ionization octave program ("SNCR_Flux_and_Ionization_($parsec distance)pc.m")

      OPEN(UNIT=121, FILE='SNCRionization.dat') !These values are in ((ions)/(cm^2*s))
      READ(121,'(E16.8)') (SNCRions_in(IK), IK=1,ZION)
      CLOSE(UNIT=121)

C       print *, "------------------- SNCRions_in ---------------------"
C       print *,
C         do IK = 1,ZION
C             print *, IK, SNCRions_in(IK)
C         end do
C       print *,


C  Interpolate ionization from old model grid [SNCRions_in(1:ZION)] to new model grid [SNCRions_out(1:Z$)]

C       print *, "------------ SNCRions_out, before interp ------------"
C       print *,
C         do IK = 1,Z$
C             print *, IK, SNCRions_out(IK)
C         end do
C       print *,

      CALL LINTERP(zz46, ZION, SNCRions_in, ZALT90, Z$, 1, 
     >               SNCRions_out)
     
C       print *, "------------ SNCRions_out, after interp ---------------"
C       print *,
C         do IK = 1,Z$
C             print *, IK, SNCRions_out(IK)
C         end do
C       print *,


      OPEN(122, FILE="SNCRionization_afterInterp.dat")
      do IK = 1,Z$
          WRITE(122,'(1X,E16.8)') SNCRions_out(IK)
      end do
      CLOSE(122)

 
C  Since SNCR input is isotropic, copy single location ionization across all latitudes
C     Convert from ((ions)/(cm^2*s)) to ((ions)/(cm^3*s)) using “deltaz”, the height of each bin converted from km to cm:

      OPEN(123, FILE="ionSource_cm2.dat")
      OPEN(124, FILE="ionSource_cm3.dat")

      do IK = 1,Z$
         do IJ = 1,L$
            ionSource(IJ,IK)=SNCRions_out(IK)

            !Write out ionSource in ((ions)/(cm^2*s)) for debugging:
               WRITE(123,'(1X,E16.8,$)') ionSource(IJ,IK)

            !Write out ionSource in ((ions)/(cm^3*s)) for debugging:
               ionSource(IJ,IK) = ionSource(IJ,IK)/(deltaz(IJ,IK)*1.0E5)
               WRITE(124,'(1X,E16.8,$)') ionSource(IJ,IK)

         end do
         WRITE(123,*)
         WRITE(124,*)
      end do

      CLOSE(123)
      CLOSE(124)

      RETURN
      END