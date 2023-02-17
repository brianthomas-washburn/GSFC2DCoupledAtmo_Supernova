C-*-Fortran-*-  -  this is used for emacs/Fortran mode capability, or just do M-x fortran-mode
C           @(#)psc.i	1.3  05/18/00
      REAL CBOLTZ,Av,pi,rad1,mu1,den1,rad2,mu2,den2,nfact1,nfact2,den
      PARAMETER (CBOLTZ=1.38E-19)
      PARAMETER (Av = 6.02e23)   ! Avogadro's number
      PARAMETER (pi = 3.14159)  
      PARAMETER (rad1 = 0.5E-4)  ! radius of NAT particles (cm)
      PARAMETER (mu1 = 117.)     ! molecular weight of NAT (grams/mole)
      PARAMETER (den1 = 1.62)    ! density of NAT (grams/cm**3)
      PARAMETER (rad2 = 7.0E-4)  ! radius of ICE particles (cm)
      PARAMETER (mu2 = 18.)      ! molecular weight of ice (grams/mole)
      PARAMETER (den2 = 0.928)   ! density of ice (grams/cm**3)
C        FACTOR TO CONVERT SOLID HNO3 NUMBER DENSITY TO NUMBER OF PARTICLES
      PARAMETER (nfact1 = (3.*mu1)/(Av*den1*4.*pi*rad1**3))
C        FACTOR TO CONVERT SOLID H2O NUMBER DENSITY TO NUMBER OF PARTICLES
      PARAMETER (nfact2 = (3.*mu2)/(Av*den2*4.*pi*(rad2**3-rad1**3)))
c         overall density of PSC2 particles, ICE and NAT
      PARAMETER (deni = (den1*rad1**3 + den2*(rad2**3-rad1**3))/rad2**3)
C         supersaturation factor for NAT and ICE
      PARAMETER (SSnat=1.0)
      PARAMETER (SSice=1.0)
