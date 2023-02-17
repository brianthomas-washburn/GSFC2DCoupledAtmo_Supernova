C  PUT IN NOX SOURCE from other ionization (in addition to GCRs)
C  AYelland 3October2020

      subroutine GCRionsIn

      include "com2d.h"
      INTEGER, parameter :: ZION = 46
      REAL*4 :: GCRions_in(1:ZION), GCRions_out(1:Z$)
C     REAL*4 :: GCRions(1:L$,1:Z$)    ***DECLARED IN "com2d_9hb_ions.h"***

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
         GCRions_in(IK)=0.d0
      end do

      do IK=1,Z$
         GCRions_out(IK) = 0.d0
      end do

      do IJ=1,L$
         do IK=1,Z$
            GCRions(IJ,IK)=0.d0
         end do
      end do

C  Read-in of Ionization values from Flux & Ionization octave program ("SNCR_Flux_and_Ionization_($parsec distance)pc.m")

      OPEN(UNIT=121, FILE='GCRions46.dat') !Units = ((ions)/(cm^2*s))
      READ(121,'(E16.8)') (GCRions_in(IK), IK=1,ZION)
      CLOSE(UNIT=121)

C       print *, "------------------- GCRions_in ---------------------"
C       print *,
C         do IK = 1,ZION
C             print *, IK, GCRions_in(IK)
C         end do
C       print *,


C  Interpolate ionization from old model grid [GCRions_in(1:ZION)] to new model grid [GCRions_out(1:Z$)]

C       print *, "------------ GCRions_out, before interp ------------"
C       print *,
C         do IK = 1,Z$
C             print *, IK, GCRions_out(IK)
C         end do
C       print *,

      CALL LINTERP(zz46, ZION, GCRions_in, ZALT90, Z$, 1, GCRions_out)
     
C       print *, "------------ GCRions_out, after interp ---------------"
C       print *,
C         do IK = 1,Z$
C             print *, IK, GCRions_out(IK)
C         end do
C       print *,

      OPEN(122, FILE="GCRions76.dat")
      do IK = 1,Z$
          WRITE(122,'(1X,E16.8)') GCRions_out(IK)
      end do
      CLOSE(122)


C  Since GCR input is isotropic, copy single location ionization across all latitudes

      OPEN(123, FILE="GCRions_cm2.dat")
      OPEN(124, FILE="GCRions_cm3.dat")

      do IK = 1,Z$
        do IJ = 1,L$
          GCRions(IJ,IK)=GCRions_out(IK)

          !Write out GCRions in ((ions)/(cm^2*s)) for debugging:
            WRITE(123,'(1X,E16.8,$)') GCRions(IJ,IK)

          !Write out GCRions in ((ions)/(cm^3*s)) for debugging:
            GCRions(IJ,IK) = GCRions(IJ,IK)/(deltaz(IJ,IK)*1.0E5) !deltaz is in (km)
            WRITE(124,'(1X,E16.8,$)') GCRions(IJ,IK)

        end do
        WRITE(123,*)
        WRITE(124,*)
      end do

      CLOSE(123)
      CLOSE(124)
      
      RETURN
      END
