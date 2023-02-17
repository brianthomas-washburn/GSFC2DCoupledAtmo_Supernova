C-*-Fortran-*-  -  this is used for emacs/Fortran mode capability, or just do M-x fortran-mode
C
	COMPLEX HEIGHT,psic,psa,psb,qprmy,qprime
c
	common/waveq/psic(n$,m$,mmax$),psa(n$,m$,mmax$),
     1  psb(n$,m$,mmax$),QBRYWX(NP$,MP$),MWBRK(NP$),MWRMP,
     2  QPRMY(N$,M$),QPRIME(N$,M$),QBRYW(N$,M$),MWVNS(4,MMAX$),
     3  GDISS(N$,M$,MMAX$),EQDMP(N$,M$)
	COMMON/PWHITE/HEIGHT(n$,mmax$)
        REAL MWBRK,MWRMP
