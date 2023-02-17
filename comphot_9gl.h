C-*-Fortran-*-  -  this is used for emacs/Fortran mode capability, or just do M-x fortran-mode
C
        INTEGER nlat,nlon,nz,nlam,nsza,no3,nx,nts,npr

ccccc	PARAMETER (nlat=18,nlon=1,nz=46)	

	PARAMETER (nlat=L$,nlon=1,nz=Z$)	
	PARAMETER (nlam=69,nsza=20,no3=12)
	PARAMETER (nx=20,nts=200)
	PARAMETER (npr=100)	       
c                       npr is now 100, for the 100 pressure levels of all the look-up tables
c   Model inputs
        COMMON/PHTAB1/sza2d(L$), ovho3(L$,Z$),
     *  sflux(nlat,nz,nlam)

c   Table arrays
        COMMON/PHTAB2/xtab(nlam,nts,nx), stab(nsza,no3,nz,nlam),
     *  pr_tab(npr), rlam(nlam), sza_tab(nsza), o3_tab(no3,nz)

c   J(O2) table array

        COMMON/PHTAB3/o2jdat(nsza,no3,nz), o2jint(nlat,nz)

