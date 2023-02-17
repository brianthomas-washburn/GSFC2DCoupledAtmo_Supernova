C-*-Fortran-*-  -  this is used for emacs/Fortran mode capability, or just do M-x fortran-mode
C
C --------------------------------------------------------
C        2-D MODEL COMMON BLOCKS
C            FOR DYNAMICS

      COMMON /PSI/ PSI(N$,M$)

C     COMMON /SLAKO/ AY(N$,M$),CY(N$,M$),B(N$,M$),AZ(N$,M$),CZ(N$,M$)
C    1,F(N$,M$)

      COMMON /MAIND/ UB(N$,M$,3),TH(N$,M$,3)

      COMMON /DERIV/ UBY(N$,M$),UBZ(N$,M$),DTHY(N$,M$),DTHZ(N$,M$)

      COMMON /SMOOTH/ FLX3(N$,M$)

      COMMON /XTIME/ NSTEP,NDAY,NYEAR,DECL,IT1,IT2,IT3,MONTH,IDAYX,
     c               imodul

      COMMON /COOLX/ GCOL(M$),T(N$,M$),COOL(N$,M$),VNC(M$),olr(N$)

      COMMON /HEAT/ GHET(M$),HEAT(N$,M$),ALT(N$,M$),GHEATX(N$,M$)

      COMMON /MOM/ PWMOM(N$,M$),GWMOM(N$,M$),RAL(M$),
     1GWMOMX(N$,M$),UTHRM(N$,M$),THRMW(N$,M$),
     2UTHRMX(NP$,MP$),QUBX(NP$,MP$),THRMWX(NP$,MP$)

      REAL KYY2,KYY,KYZ,KZZ
      COMMON /KS/ KYY2(N$,M$),KYY(N$,M$),KYZ(N$,M$),KZZ(NP$,MP$),
     1PWKYY(NP$,MP$),       ZIKYY(N$,M$),         ZIKZZ(NP$,MP$),
     2DMOL(MP$), GWKZZ(NP$,MP$)

      COMMON /FLUXES/ VERF(N$,M$),HORF(N$,M$),EPFLX(NP$,MP$)

      COMMON /TWODWK/ DUM1(N$,M$),DUM2(N$,M$),DUM3(N$,M$),D3(MP$),
     1DUM4(N$,M$),DUMV(NPP$,MPP$),YPV(NPP$),DUM1X(NP$,MP$),
     2DUM2X(NP$,MP$),D3X(N$+10),ZPW(MPP$),DUM3X(NP$,MP$),
     3DUM4X(NP$,MP$)

      COMMON /TROP / TFROP(N$,M$,13),TEHE(N$,10,12)

      COMMON /DUMP_1/ IDUMP,IFILE,isw2(isw$)

      COMMON /CHEM_1/ CHGAMA(MP$),CHTHET(MP$)

      COMMON /UMAX/ UMAX,UMIN,TMAXx,TMINx,VMAX,VMIN,WMAX,WMIN

      COMMON /DIFFU/ DIFZU(N$,M$),DIFZT(N$,M$),TMIX(N$,M$),
     1               UMIX(N$,M$)         

      COMMON/O3/CO2M(N$,M$),O3M(N$,M$),H2OM(N$,M$),WLAT(18),ZWAT(23),
     1WATZ(18,23,12),XLSBUV(18),ZSBUV(17),SBUV(18,17,12)

      COMMON/JLS/QBARY(Np$,MP$),QBAR(NP$,MP$),ubw(n$,m$)

      common/gforc/ glat(181),gavg(181),gvar(181),gvint(np$)

      COMMON/CDRAG/ dragcoup(7,NP$,MP$), tempbx(NP$,MP$), 
     >              HEATG5X(NP$,MP$), COOLG5X(NP$,MP$)

C -------------------- END COMMON BLOCKS ---------------------------- 
