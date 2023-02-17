C-*-Fortran-*-  -  this is used for emacs/Fortran mode capability, or just do M-x fortran-mode
C---jv_std5.i---COMMON BLOCKS for UCI photolysis code, version 5.0 (1Nov93)
C----includes new SR & NO plus allowance for Mie scattering (see jv_stdm.cmn)
C----set for 78/79 wavelengths (NWLI) ,41/42 altitudes(NHT1), 48+4=52 J-values(NJR)
C
C             @(#)jv_std5.i	1.4  07/15/96
C
C        UCI parameters
      CHARACTER*80 TITLE0,TIFILE(8),TITRAJ(8),TITLEV
      COMMON /TITLES/ TITLE0,TITLEJ,TITRAJ,TIFILE,TITLEV
      COMMON /ATMOS/T(NHT1),REFO3(NHT1),P(NHT1),DM(NHT1),DO3(NHT1),
     $     Z(NHT1),AER(NHT1),WTAU(NHT1,NHT1),RAD,RFLECT,SZA,U0,TANHT,
     $     ZZHT,NLBATM,NC
      COMMON /CCSCAT/ TAU(NHT1+1),PIRAY(NHT1+1),PIAER(NHT1+1),
     .        FLTAU(NHT1+1),TTFBOT,NTT
      COMMON /CCWVL/ WBIN(NWLI+1),WL(NWLI),FL(NWLI)
      COMMON /CCSRB/ODF(NODF,NSR),FNO(NSR),QNO(NODF,NSR),
     $     O2X(NODF,NSR,3),QO2(NWLI,3),QO3(NWLI,3),TQQ(3,3),
     $     QMI1M,QMIWVL,ISR(NSR)
      COMMON /CCJVL/ VNO(NHT1),VO2(NHT1),FFF(NWLI,NHT1)
      COMMON /CCMIE/ FFFMIE(7,NHT1),ASTRAT,WEDG(8),KEDG(8),MSTRAT
