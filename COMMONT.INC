C-*-Fortran-*-  -  this is used for emacs/Fortran mode capability, or just do M-x fortran-mode
C
C --------------------------------------------------------
C        2-D MODEL COMMON BLOCKS
C            FOR TRANSPORT


      COMMON /RESID/ WS(N$,M$),VS(N$,M$)
      COMMON /VELRBR/ WRBR(NP$,MP$),VRBR(NP$,MP$),
     1                WRB0(NP$,MP$),VRB0(NP$,MP$),
     2                GZPW(NP$,MP$),GYPW(NP$,MP$)

      COMMON /WNDJTB/ RALP(MP$),UBX(NP$,MP$),VFRC(NP$,MP$),
     1THX(NP$,MP$),HEATX(NP$,MP$),DRAGX(NP$,MP$),
     2PSIX(NP$,MP$)


      COMMON /TRACER_D/ X_D(NP$,MP$,NCON),SOR_D(NP$,MP$,NCON),
     1X0_D(NP$,MP$,NCON),RHO0_D(NP$,MP$),RHOX_D(NP$,MP$),RHOY_D(NP$,MP$)
     2,XX_D(NP$,MP$,NCON),XY_D(NP$,MP$,NCON),XXX_D(NP$,MP$,NCON),
     3XXY_D(NP$,MP$,NCON),XYY_D(NP$,MP$,NCON),XC_D(NP$,MP$,NCON),
     4SM_D(NP$,MP$),
     5XC00_D(NP$,MP$,NCON),XCT_D(NP$,MP$,NCON),
     5XMIX_D(NP$,MP$,NCON),
     6RHO00_D(NP$,MP$),RHOX0_D(NP$,MP$),RHOY0_D(NP$,MP$)


      COMMON /N2O/ XN2OPHX(NP$,MP$)

      dimension X(NP$,MP$,NCON),SOR(NP$,MP$,NCON),
     1X0(NP$,MP$,NCON),RHO0(Np$,Mp$),RHOX(Np$,Mp$),RHOY(Np$,Mp$),
     2XX(NP$,MP$,NCON),XY(NP$,MP$,NCON),XXX(NP$,MP$,NCON),
     3XXY(NP$,MP$,NCON),XYY(NP$,MP$,NCON),XC(NP$,MP$,NCON),
     4SM(NP$,MP$),
     5XC00(NP$,MP$,NCON),XCT(NP$,MP$,NCON),
     5XMIX(NP$,MP$,NCON),
     6RHO00(Np$,Mp$),RHOX0(Np$,Mp$),RHOY0(Np$,Mp$)


	equivalence (X_D(1,1,1),X(1,1,1))
	equivalence (SOR_D(1,1,1),SOR(1,1,1))
	equivalence (X0_D(1,1,1),X0(1,1,1))
	equivalence (RHO0_D(1,1),RHO0(1,1))
	equivalence (RHOX_D(1,1),RHOX(1,1))
	equivalence (RHOY_D(1,1),RHOY(1,1))
	equivalence (XX_D(1,1,1),XX(1,1,1))
	equivalence (XY_D(1,1,1),XY(1,1,1))
	equivalence (XXX_D(1,1,1),XXX(1,1,1))
	equivalence (XXY_D(1,1,1),XXY(1,1,1))
	equivalence (XYY_D(1,1,1),XYY(1,1,1))
	equivalence (XC_D(1,1,1),XC(1,1,1))
	equivalence (SM_D(1,1),SM(1,1))
	equivalence (XC00_D(1,1,1),XC00(1,1,1))
	equivalence (XCT_D(1,1,1),XCT(1,1,1))
	equivalence (XMIX_D(1,1,1),XMIX(1,1,1))
	equivalence (RHO00_D(1,1),RHO00(1,1))
	equivalence (RHOX0_D(1,1),RHOX0(1,1))
	equivalence (RHOY0_D(1,1),RHOY0(1,1))

C -------------------- END COMMON BLOCKS ---------------------------- 












