C-*-Fortran-*-  -  this is used for emacs/Fortran mode capability, or just do M-x fortran-mode
C        from AER (grid.i) now com_aerg.h   (7/8/02)
C
      INTEGER NLT, NHT, NLT1, NHT1, NLT2, NHT2
      INTEGER NTN, NLV, NPOL, NDIST

      PARAMETER (NLT=45,NHT=76,DZETA=0.5/3.,Psurf=1000.)
      PARAMETER (NLT1=NLT+1,NHT1=NHT+1,NLT2=2*NLT+1,NHT2=2*NHT+1)
      PARAMETER (NTN=14,NLV=21,NPOL=8,NDIST=131,TD1=170.)
