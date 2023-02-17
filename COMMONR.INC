C-*-Fortran-*-  -  this is used for emacs/Fortran mode capability, or just do M-x fortran-mode
C
      REAL LAT

      COMMON/RADINI0/MZ,MZ1,MZTOT,ILBEDO(N$)
      COMMON/RADINI1/PTOT(91),RHOTOT(91),TTOT(N$,91),LAT(37),ZPREV(M$),
     +              ZJRREV(91),ZJRMID(45),ZJR(91)
      COMMON/RADCAL0/HEATAV(N$,45),COOLJR(N$,45),VNCJR(N$,45)
      INCLUDE 'COMMONR1.INC'
