C-*-Fortran-*-  -  this is used for emacs/Fortran mode capability, or just do M-x fortran-mode
C
C --------------------------------------------------------
C      2-D MODEL PARAMETER STATEMENTS
C

      INTEGER S$C

c dynamics module grid points in horizontal 
      PARAMETER (N$=36)
c dynamics module grid points in vertical
      PARAMETER (M$=46)
      
      PARAMETER (N1$=N$-1)
      PARAMETER (M1$=M$-1)

c  number of constituents = NCON
      PARAMETER (NCON=2)

      PARAMETER (MP$=M$+1)
      PARAMETER (NP$=N$+1)
      PARAMETER (NPP$=NP$+1)
      PARAMETER (MPP$=MP$+1)

C    NUMBER OF ZONAL WAVENUMBERS IN PW SCHEME
      PARAMETER (MMAX$ = 5 )

      PARAMETER (NSWTCHS=70)
C     parameter (isw$=26+mmax$+ncon)
      parameter (isw$=70)
      
C     CHEMISTRY MODULE PARAMETERS - remember to CHANGE THIS w/ COM2D changes!!!
      PARAMETER (S$C=80+10)
      PARAMETER (NS$=45)
      PARAMETER (MS$=76)      
