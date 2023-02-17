C
C
        SUBROUTINE XCOUPLED_IN(IYR0, IIDAY360)

C
C   this sets up the initializations for the coupled model
C

        include "comcfg.h"
        include "timer.inc"


        PARAMETER (NSWTCHS=70)
        COMMON /SWITC0/ ISWD(NSWTCHS), LH_SW(12), KZZ_SW(15)

        INTEGER IYR0, IIDAY360, L$C, Z$C


cjer common for dynamics writeout

        common /dump_2/ iwrite,yrjer
        data iwrite /0/


        read(60,*) iswo3
        read(60,*) iswh2o
        read(60,*) iswhetchem
        read(60,*) iswvolc
        read(60,*) iswfixed
        read(60,*) iswcircout
        read(60,*) iswsolcyc

        close(unit=60)

        write(6,*) 'iswo3 = ',iswo3
        write(6,*) 'iswh2o = ',iswh2o
        write(6,*) 'iswhetchem = ',iswhetchem
        write(6,*) 'iswvolc = ',iswvolc
        write(6,*) 'iswfixed = ',iswfixed
        write(6,*) 'iswcircout = ',iswcircout
        write(6,*) 'iswsolcyc = ',iswsolcyc



C   Added 8/28/96, by PEM - Initialize 2D Dynamics Calculation
c      print *,'before twods'

        CALL TWODS

c      print *,'after twods'
C   End Added 8/28/96.


C     Initialize the new timing routines (timer.f).  This needs to be after
C     CALL TWODS so that the ISW array (here ISWD) in COMMON /SWITC0/ will be
C     initialized.

        IDOY0  = IIDAY360
        ICHSTPD0 = 1
        IDYSTPD0 = ISWD(6)

        CALL TIME_INIT(IYR0, IDOY0, ICHSTPD0, IDYSTPD0, -1)


        RETURN    
        END
