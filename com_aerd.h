C-*-Fortran-*-  -  this is used for emacs/Fortran mode capability, or just do M-x fortran-mode
C        from AER (diurnal.i) now com_aerd.h   (7/8/02)

      INTEGER NTIME, NOON, MIDNIT, RIS90, SET90, SUNRIS, SUNSET
      INTEGER DAYAVG, NITAVG, DIURAVG
      INTEGER NDIVDAY, NDIVTWI, NDIVNIT

C                                                     ! NTIME = 18
      PARAMETER (NDIVDAY=5, NDIVTWI=1, NDIVNIT=5,
     &     NTIME=2*(NDIVDAY+NDIVTWI) + NDIVNIT + 1)

C                                                     ! NTIME = 33
Cef      PARAMETER (NDIVDAY=10, NDIVTWI=1, NDIVNIT=10,
Cef     &     NTIME=2*(NDIVDAY+NDIVTWI) + NDIVNIT + 1)

      PARAMETER (NOON=1, MIDNIT=(NTIME+1)/2,
     &     DAYAVG=NTIME+1, NITAVG=NTIME+2, DIURAVG=NTIME+3)

      COMMON/DIURTIME/RIS90,SET90,SUNRIS,SUNSET
