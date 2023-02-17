C  PUT IN NOX SOURCE from other ionization (in addition to GCRs)
C  AYelland 3October2020

      subroutine ionsRatio_calc

      include "com2d.h"
C     REAL*4 :: ionsRatio(1:L$,1:Z$), ionSource(1:L$,1:Z$), GCR(1:L$,1:Z$)
C           ***DECLARED IN "com2d_9hb_ions.h"***

C Intitalizing Array with zeros.

      do IJ=1,L$
         do IK=1,Z$
            ionsRatio(IJ,IK)=0.d0
         end do
      end do

C Calculating the ratio of ionization between "ionSource" (max SNCR ions) and "GCRions" (background ions)
C     >>> This is used as the multiplier for the lightning effect calculation in "solv_9hb.f"

      do IJ = 1,L$
         do IK = 1,Z$
            ionsRatio(IJ,IK) = ionSource(IJ,IK) / GCRions(IJ,IK)
         end do
      end do

C Saving the ionsRatio/lightning multiplier to output file
C       OPEN(122, FILE="gcr_model.dat")
C       do IK = 1,Z$
C          do IJ = 1,L$
C             WRITE(122,'(1X,E16.8,$)') GCR(IJ,IK)
C          end do
C          WRITE(122,*)
C       end do
C       CLOSE(122)

      OPEN(123, FILE="ionsRatio.dat")
      do IK = 1,Z$
         do IJ = 1,L$
            WRITE(123,'(1X,E16.8,$)') ionsRatio(IJ,IK)
         end do
         WRITE(123,*)
      end do
      CLOSE(123)
      
      RETURN
      END
