C-*-Fortran-*-  -  this is used for emacs/Fortran mode capability, or just do M-x fortran-mode
C           %W%  %G%    
C           $Id: species.i,v 1.42 2007/06/25 19:08:27 dkweis Exp $

      INTEGER NFSP, NRSP, NTSP, NPLQ ,NFAM, NFP, NSSP, NCSP
      INTEGER NJR, NKR, NRP, NWLI, NW1, NW2
      INTEGER NWSRB, NSR, NODF, NSPINT
      INTEGER NO3PL, NO3JR, NZJR, NZKR

      PARAMETER (NFSP=54,NRSP=56,NTSP=116,NPLQ=115,NFAM=10,NFP=10)
      PARAMETER (NSSP=NTSP-NFSP,NCSP=NFSP+NRSP)
      PARAMETER (NJR=121,NKR=382,NRP=5,NWLI=77,NW1= 1,NW2=77)
      PARAMETER (NWSRB=0,NSR=15,NODF= 6,NSPINT=12)
      PARAMETER (NO3PL=40,NO3JR=5,NZJR=1,NZKR= 8)
      PARAMETER (IO3     =  4,     INOX    = 91,     ICLX    = 92)
      PARAMETER (IBRX    = 93,     INOZ    =  0,     IHNO3   = 39)
      PARAMETER (IPAN    = 96,     IHCL    = 23,     IHOCL   = 24)
      PARAMETER (IHOBR   = 44,     IIN2O5   = 36,    IH2O    = 98)
      PARAMETER (IO2     =114,     IM      =116,     INOY    =112)
      PARAMETER (ICLY    =113,     ISHNO3  = 97,     ICLNO3  = 25)
      PARAMETER (ISH2O   = 99,     IISO2    =106,    IH2SO4  =108)
      PARAMETER (IOX=90, IBRNO3= 45)
      PARAMETER (INOXPR  = 74,     INOZPR  =  0,     IHNO3PR =  0)
      PARAMETER (IPANPR  = 85,     INOZLO  =  0,     IHNO3LO =  0)
      PARAMETER (IPANLO  = 86,     ICLYPR  =  0,     ICLXPR  = 77)
      PARAMETER (IHCLLO  =  0,     IHCLPR  =  0)
      PARAMETER (IQ1O2   =  1,     IQ1O3   =  2,     IQ2O3   =  3)
      PARAMETER (IQ1NO   =  5,     IQ2HNO3 = 25)
