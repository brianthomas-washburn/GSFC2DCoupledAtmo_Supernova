c

      SUBROUTINE HETCHEM(RAER,nkr,xrhcl,xrclono2,ij,ik)

chettest             SUBROUTINE HETCHEM(RAER,nkr,xhcl,xclono2,ij,ik,iang)



C  *************************************************************************************************
C
C   routine to ONLY COMPUTE THE HET REACTION RATES WITHIN THE DIURNAL LOOP
C       the NAT and ICE aerosol PSCs are now computed separately once per day in subroutine NATICE
C       which was taken from David's JPL-2000 HETCHEM routine
C
C      adapted from David's hetchem_9aa_jpl00.f  - EF 3/12/03
C
C  *************************************************************************************************


c     David B. Considine, 9/29/93

c     This is a long-overdue subroutine to automate the calculation of 
c     the amount of HNO3 that is in solid phase, and calculate the 
c     surface area density of Type 1 NAT aerosols.  It uses probability
c     density distributions that are calculated offline and read in 
c     in the initialization section of the subroutine.

c     VERSION 2 created 10/5/93:  This version includes Type 2 PSCs.
c     The type 2 calculation occurs before the type 1 calculation.

c     VERSION 3 created 11/1/93:  This version transports dehydration.
c     changes to previous version allowed keeping account of
c      dehydration,
c     but did not allow transport.  This version also uses diffusion
c     terms for solid HNO3, H2O, and dehydration.
c    
c     VERSION 4 created 11/19/93:  This version calculates the 
c     heterogeneous
c     reaction rates here rather than in react.  Reaction rates are
c     calculated every day.
c
c     VERSION 4a created 11/30/93:  This version sediments the nat
c     aerosols by generating an average fall velocity based on the
c     fall velocity of both the NAT and the ICE aerosols, using f(ice)
c     and f(nat)
c
c     VERSION 4b created 12/9/93: This version is created for a model
c     where H2O is transported, so the ice aerosol stuff works just like
c     the nat aerosols: c(15) is transported h2o, c(61) doesn't exist
c     anymore (it used to be dehydration when c(15) was a constant
c     field taken from LIMS data), c(67) is solid H2O.
c
c     VERSION 4c created 12/9/93: This version calculates the fall 
c     velocities
c     of the nat and ice aerosols using the mode radius and a formula
c     from Kasten, JAM 7 944-947 (1968).
c
c     VERSION 4d created 12/15/93: This version takes allows for a 
c     correction due to supersaturation.
c
c     VERSION 4e created 12/23/93: This version calculates the rates
c     for the sulfate ClONO2+H2O and ClONO2+HCl using the zonal mean
c     temp, as suggested by Ravishankara.
c
c     VERSION 4f created 12/28/93: This version allows you to turn off
c     the PSC calculation and retain the calculation of sulfate gammas
c     and reaction rates.
c
c     VERSION 4g created 12/29/93: This version allows you to 
c     independently
c     turn off type 1 or type 2 PSCs.  It also allows independent
c     specifications of supersaturation corrections for Type 1 and 
c     Type 2 PSCs
c
c     VERSION 4h: 3/29/84 This version provides a switch to turn 
c     sedimentation
c     on or off.  It also sediments implicitly and uses a correction 
c     term to account for effective lognormal particle size 
c     distribution fluxes
c
c     VERSION 4i: 4/15/94: This version adjusts the sedimentation
c     time step to allow for fall velocities greater that 1 gridbox/
c     timestep.
c
c     VERSION 4j: 6/17/94: This version corrects the surface area 
c     density
c     calculation mistake that was discovered by Dave Usinski, where
c     in the conversion to surface area density I was 
c     using ln^2(sigma)/4 instead of 5ln^2(sigma)/2.
c
c     VERSION 4k: 9/30/94: This version uses a new probability 
c     distribution
c     calculated from 15 years of nmc data (1979-1993, inclusive).
c     the file is named "probdist2.dat"
c
c     VERSION 4l: 10/03/94: This version changes the probability 
c     distribution
c     daily rather than monthly, interpolating between nearest two 
c     months.
c    
c     VERSION 4m: 11/28/94: This version calculates a sticking 
c     coefficient
c     and pseudo-bimolecular reaction rate for HOCl+HCl -> Cl2+H2O. In
c     addition, the Hanson and Ravishankara formulation for ClONO2 + HCl
c     and ClONO2 + H2O presented in J. Phys Chem. 98, 5728 (1994) is 
c     used.
c
c     VERSION 4m_newmod: 01/08/95: This version is changed to be 
c     compatible with the new model structure.
c
c     VERSION hetchem_95_a.f: 02/16/95: This version includes the
c     heterogeneous reaction BrONO2+H2O -> HOBr+HNO3
c----------------------------------------------------------------------c


      include 'com2d.h'

C                                          KHA(RH$,L$,Z$) array is just internal to this SUBROUTINE
      INTEGER NKR

      REAL*8 RAER(NKR)                        !!, kha(7,L$,Z$)
      REAL*8 xhcl, xclono2, xrhcl, xrclono2

      REAL*8 pph2o, logpph2o, pphcl, ppclono2, droprad
      REAL*8 tt, wtpct, molsulf, aw, tmin, dtemp 
      REAL*8 vn2o5, vclono2, vhocl, vhobr, vbrono2, gamhydsulf
      REAL*8 khgs, khgn, khgi


      SAVE                         ! save all internal variables



CPRINT
CPRINT                  if (ij .eq. 9 .and. ik .eq. 15  .or.  
CPRINT     >                ij .eq. 1 .and. ik .eq. 8) then
ccc                  write(59,*) '  '
CPRINT         write(59,113) iday360,ij,ik, nkr, raer(6), raer(7), 
CPRINT     >       raer(90), raer(92), raer(164), raer(165), raer(169)
CPRINT 113              format('HETCHEM',4I4, 5X, 1P7E10.2)
CPRINT                     write(59,*) ' '
CPRINT                   end if



C   load in xhcl, xclono2, set lower limit (1.e-30 is good) - avoids them going to 0. which gives NaNs

            xhcl = xrhcl
            xclono2 = xrclono2

            if (xhcl    .lt. 1.e-30)    xhcl = 1.e-30
            if (xclono2 .lt. 1.e-30) xclono2 = 1.e-30



C  IK40 = index for 40 km - defined in MAIN, set KHs = 0.0 above IK40  AND BELOW 6 km and RETURN - ZALT(Z$X)
C     KHG(rht$=42,L$,Z$) in COMMON, and array RAER for the AER fast chemistry
C       ie, set KHG=0 here, then go to 2000 below which loads KHG array into the RAER array

       IF (ik .gt. IK40  .or.  ZALT(ik) .lt. 6) then

            do 456 ihr=1,rht$
 456           khg(ihr,ij,ik) = 0.0

ccef            do 555 ihr=164,170
ccef 555           RAER(ihr) = 0.0

          go to 2000
       ENDIF


c     Section to calculate heterogeneous reaction rates:



c Partial pressure of H2O in mb:

                pph2o=(cn(15,ij,ik))
     c                /m(ij,ik)*press(ik)

                logpph2o=DLOG(pph2o)

c Partial pressure of HCl and ClONO2 in mb: - HCl and ClONO2 now within DIURNAL LOOP

                pphcl = xhcl/m(ij,ik)*press(ik)
                ppclono2 = xclono2/m(ij,ik)*press(ik)


c assume a monodisperse aerosol with a droplet radius of .1 micron

         droprad=0.1d-4


C  now loop through temp. PDF, first initialize all summing variables, TPROB(211,L$,Z$)
C                                  loop through from 120-330 K, use ITP1(L$,Z$), ITP2(L$,Z$)
       do 1004 igf = 1,rht$
 1004     khg(igf,ij,ik) = 0.0


       do 4000 itemp=itp1(ij,ik)+119, itp2(ij,ik)+119
            dtemp = DBLE(itemp) ! get real*8 value
                                                        !if using zonal mean temps, use TEMP(L$,Z$X) - REAL*8
            if (IZONAVGT .eq. 1) dtemp = TEMP(ij,ik)

c get wt pct of sulfate aerosols by calling subroutine wtpctsulf:
c subroutine also calculates molality of sulfate and water activity,
c needed below:

      call wtpctsulf(dtemp,tt,pph2o,wtpct,molsulf,aw,tmin)


c get gamma for BrONO2 + H2O on sulfate:

      gs(6) = 1.d0/(1.d0/0.805+1./(0.114d0+DEXP(29.24-0.396*wtpct)))


c get gamma for N2O5 + H2O on sulfate - replace gs(3) which was
c read in from reac_het:

      gs(3) = gamhydsulf(tt,wtpct)


c get gammas for ClONO2 + HCl, ClONO2 + H2O, and HOCl + HCl by
c calling subroutine gamclonitsulf:

       call gamclonitsulf(ij,ik,tt,molsulf,wtpct,pphcl,ppclono2,
     c               aw,droprad,gs(1),gs(2),gs(5),gs(7),iday360)


c      write(6,*)ij,ik,dtemp,molsulf,wtpct,pphcl,ppclono2,aw,
c     c          droprad,gs(1),gs(2),gs(5)


c thermal velocity of N2O5, ClONO2, HOCl, BrONO2, NO3, and Cl2:

      vn2o5=(8.d0*1.38D-16*dtemp/108./1.67D-24/3.141593)**.5
      vhocl=(8.d0*1.38D-16*dtemp/52.46/1.67D-24/3.141593)**.5
      vhobr=(8.d0*1.38D-16*dtemp/96.92/1.67D-24/3.141593)**.5
      vclono2=(8.d0*1.38D-16*dtemp/97.46/1.67D-24/3.141593)**.5
      vbrono2=(8.d0*1.38D-16*dtemp/141.91/1.67D-24/3.141593)**.5

cccc      vno3=(8.d0*1.38D-16*dtemp/62.005/1.67D-24/3.141593)**.5
cccc      vcl2=(8.d0*1.38D-16*dtemp/70.91/1.67D-24/3.141593)**.5


C  KHG(rht$,L$,Z$) - now varies diurnally, and load into AER R array 
C     for H2O reactions, multiply by H2O to be compatible with AER fast chemistry code
C
c  ClONO2 + HCl -> Cl2 + HNO3

c      khg(1,ij,ik)=vclono2/(4.0*(xhcl + 1.0))*aerosol(ij,ik)*gs(1) 
c      khg(1+7,ij,ik)=vclono2/(4.0*(xhcl + 1.0))*nataer(ij,ik)*gn(1)
c      khg(1+14,ij,ik)=vclono2/(4.0*(xhcl + 1.0))*aerice(ij,ik)*gi(1)

      khgs = vclono2/(4.d0)*aerosol(ij,ik)*gs(1) 
      khgn = vclono2/(4.d0)*nataer(ij,ik)*gn(1)
      khgi = vclono2/(4.d0)*aerice(ij,ik)*gi(1)

      khg(1,ij,ik)    = khg(1,ij,ik)    + khgs*TPROB(itemp-119,ij,ik)
      khg(1+14,ij,ik) = khg(1+14,ij,ik) + khgn*TPROB(itemp-119,ij,ik)
      khg(1+28,ij,ik) = khg(1+28,ij,ik) + khgi*TPROB(itemp-119,ij,ik)


cc  ggsfc(rht$,L$,Z$)  array is for output checking only

      ggsfc(1,ij,ik) = gs(1) 
      ggsfc(1+14,ij,ik) = gn(1) 
      ggsfc(1+28,ij,ik) = gi(1) 
                                                    

cc      if (ij .eq. 9  .and.  ik .eq. 10) then
cc         write(28,1051)
cc      write(28,1050)ij,ik,khg(1,ij,ik),gs(1), aerosol(ij,ik),vclono2, 
cc     >     aw, wtpct, xhcl, m(ij,ik), dtemp
cc 1050    format(2I4, 4x, '91 ', 4x, 1P9E11.3)
cc         write(28,1051)
cc 1051    format(4x)
cc      endif
cc

c ClONO2 + H2O -> HOCl + HNO3

c      khg(2,ij,ik)=vclono2/(4.*(cn(15,ij,ik)+1.))*aerosol(ij,ik)*gs(2)
c      khg(2+7,ij,ik)=vclono2/(4.*(cn(15,ij,ik)+1.))*nataer(ij,ik)*gn(2)
c      khg(2+14,ij,ik)=vclono2/(4.*(cn(15,ij,ik)+1.))*aerice(ij,ik)*gi(2)

      khgs = vclono2/4.d0*aerosol(ij,ik)*gs(2)
      khgn = vclono2/4.d0*nataer(ij,ik)*gn(2)
      khgi = vclono2/4.d0*aerice(ij,ik)*gi(2)

      khg(2,ij,ik)    = khg(2,ij,ik)    + khgs*TPROB(itemp-119,ij,ik)
      khg(2+14,ij,ik) = khg(2+14,ij,ik) + khgn*TPROB(itemp-119,ij,ik)
      khg(2+28,ij,ik) = khg(2+28,ij,ik) + khgi*TPROB(itemp-119,ij,ik)


      ggsfc(2,ij,ik) = gs(2) 
      ggsfc(2+14,ij,ik) = gn(2) 
      ggsfc(2+28,ij,ik) = gi(2) 


cc      if (ij .eq. 9  .and.  ik .eq. 10) then
cc         write(28,1051)
cc      write(28,2050)ij,ik,khg(2,ij,ik),gs(2),aerosol(ij,ik),vclono2,
cc     >    aw, wtpct, xhcl, m(ij,ik), dtemp
cc 2050    format(2I4, 4x, '93 ', 4x, 1P9E11.3)
cc         write(28,1051)
cc      endif



c N2O5 + H2O -> 2HNO3

c      khg(3,ij,ik)=vn2o5/(4.0*(cn(15,ij,ik)+1.))*aerosol(ij,ik)*gs(3)
c      khg(3+7,ij,ik)=vn2o5/(4.0*(cn(15,ij,ik)+1.))*nataer(ij,ik)*gn(3)
c      khg(3+14,ij,ik)=vn2o5/(4.0*(cn(15,ij,ik)+1.))*aerice(ij,ik)*gi(3)

      khgs = vn2o5/4.d0*aerosol(ij,ik)*gs(3)
      khgn = vn2o5/4.d0*nataer(ij,ik)*gn(3)
      khgi = vn2o5/4.d0*aerice(ij,ik)*gi(3)

      khg(3,ij,ik)    = khg(3,ij,ik)    + khgs*TPROB(itemp-119,ij,ik)
      khg(3+14,ij,ik) = khg(3+14,ij,ik) + khgn*TPROB(itemp-119,ij,ik)
      khg(3+28,ij,ik) = khg(3+28,ij,ik) + khgi*TPROB(itemp-119,ij,ik)


      ggsfc(3,ij,ik)    = gs(3) 
      ggsfc(3+14,ij,ik) = gn(3) 
      ggsfc(3+28,ij,ik) = gi(3) 


chetest
cc      if (ij .eq. 9  .and.  ik .eq. 10) then
cc         write(28,1051)
cc      write(28,3050)ij,ik,khg(3,ij,ik),gs(3),aerosol(ij,ik), vn2o5, 
cc     >    aw, wtpct, xhcl, m(ij,ik), dtemp
cc 3050    format(2I4, 4x, '92 ', 4x, 1P9E11.3)
cc         write(28,1051)
cc      endif
chetest



c N2O5 + HCl -> ClONO + HNO3

c      khg(4,ij,ik)=vn2o5/(4.0*(xhcl + 1.0))*aerosol(ij,ik)*gs(4)
c      khg(4+7,ij,ik)=vn2o5/(4.0*(xhcl + 1.0))*nataer(ij,ik)*gn(4)
c      khg(4+14,ij,ik)=vn2o5/(4.0*(xhcl + 1.0))*aerice(ij,ik)*gi(4)

      khgs = vn2o5/(4.d0)*aerosol(ij,ik)*gs(4)
      khgn = vn2o5/(4.d0)*nataer(ij,ik)*gn(4)
      khgi = vn2o5/(4.d0)*aerice(ij,ik)*gi(4)

      khg(4,ij,ik)    = khg(4,ij,ik)    + khgs*TPROB(itemp-119,ij,ik)
      khg(4+14,ij,ik) = khg(4+14,ij,ik) + khgn*TPROB(itemp-119,ij,ik)
      khg(4+28,ij,ik) = khg(4+28,ij,ik) + khgi*TPROB(itemp-119,ij,ik)


      ggsfc(4,ij,ik) = gs(4) 
      ggsfc(4+14,ij,ik) = gn(4) 
      ggsfc(4+28,ij,ik) = gi(4) 


cc      if (ij .eq. 9  .and.  ik .eq. 10) then
cc         write(28,1051)
cc      write(28,4050)ij,ik,khg(4,ij,ik),gs(4),aerosol(ij,ik),vn2o5, 
cc     >    aw, wtpct, xhcl, m(ij,ik), dtemp
cc 4050    format(2I4, 4x, 'xx ', 4x, 1P9E11.3)
cc         write(28,1051)
cc      endif




c HOCl + HCl -> Cl2 + H2O

c      khg(5,ij,ik)=vhocl/(4.0*(xhcl + 1.0))*aerosol(ij,ik)*gs(5)
c      khg(5+7,ij,ik)=vhocl/(4.0*(xhcl + 1.0))*nataer(ij,ik)*gn(5)
c      khg(5+14,ij,ik)=vhocl/(4.0*(xhcl + 1.0))*aerice(ij,ik)*gi(5)

      khgs = vhocl/(4.d0)*aerosol(ij,ik)*gs(5)
      khgn = vhocl/(4.d0)*nataer(ij,ik)*gn(5)
      khgi = vhocl/(4.d0)*aerice(ij,ik)*gi(5)

      khg(5,ij,ik)    = khg(5,ij,ik)    + khgs*TPROB(itemp-119,ij,ik)
      khg(5+14,ij,ik) = khg(5+14,ij,ik) + khgn*TPROB(itemp-119,ij,ik)
      khg(5+28,ij,ik) = khg(5+28,ij,ik) + khgi*TPROB(itemp-119,ij,ik)


      ggsfc(5,ij,ik) = gs(5) 
      ggsfc(5+14,ij,ik) = gn(5) 
      ggsfc(5+28,ij,ik) = gi(5) 


cc      if (ij .eq. 9  .and.  ik .eq. 10) then
cc         write(28,1051)
cc      write(28,5050)ij,ik,khg(5,ij,ik),gs(5),aerosol(ij,ik),vhocl,  
cc     >    aw, wtpct, xhcl, m(ij,ik), dtemp
cc 5050    format(2I4, 4x, '90 ', 4x, 1P9E11.3)
cc         write(28,1051)
cc      endif


c                                                    ! KH5 sometimes goes negative above 27 km, reset here
cccc      if (khg(5,ij,ik) .le. 0.) khg(5,ij,ik) = 0.0



c BrONO2 + H2O -> HOBr + HNO3

c      khg(6,ij,ik)=vbrono2/(4.0*(cn(15,ij,ik)+1.))*aerosol(ij,ik)*gs(6) 
c      khg(6+7,ij,ik)=vbrono2/(4.0*(cn(15,ij,ik)+1.))*nataer(ij,ik)*gn(6)
c      khg(6+14,ij,ik)=vbrono2/(4.*(cn(15,ij,ik)+1.))*aerice(ij,ik)*gi(6)

      khgs = vbrono2/4.d0*aerosol(ij,ik)*gs(6) 
      khgn = vbrono2/4.d0*nataer(ij,ik)*gn(6)
      khgi = vbrono2/4.d0*aerice(ij,ik)*gi(6)

      khg(6,ij,ik)    = khg(6,ij,ik)    + khgs*TPROB(itemp-119,ij,ik)
      khg(6+14,ij,ik) = khg(6+14,ij,ik) + khgn*TPROB(itemp-119,ij,ik)
      khg(6+28,ij,ik) = khg(6+28,ij,ik) + khgi*TPROB(itemp-119,ij,ik)


      ggsfc(6,ij,ik) = gs(6) 
      ggsfc(6+14,ij,ik) = gn(6) 
      ggsfc(6+28,ij,ik) = gi(6) 


cc      if (ij .eq. 9  .and.  ik .eq. 10) then
cc         write(28,1051)
cc      write(28,6050)ij,ik,khg(6,ij,ik),gs(6),aerosol(ij,ik),vbrono2, 
cc     >    aw, wtpct, xhcl, m(ij,ik), dtemp
cc 6050    format(2I4, 4x, '123', 4x, 1P9E11.3)
cc         write(28,1051)
cc      endif




c HOBr + HCl -> Cl2 + H2O
                             ! set gs(7) = 0.0:  following JPL-2002/JPL-2006
      gs(7) = 0.d0

c      khg(7,ij,ik)=vhobr/(4.0*(xhcl + 1.0))*aerosol(ij,ik)*gs(7)
c      khg(7+7,ij,ik)=vhobr/(4.0*(xhcl + 1.0))*nataer(ij,ik)*gn(7) 
c      khg(7+14,ij,ik)=vhobr/(4.0*(xhcl + 1.0))*aerice(ij,ik)*gi(7)

      khgs = vhobr/(4.d0)*aerosol(ij,ik)*gs(7)
      khgn = vhobr/(4.d0)*nataer(ij,ik)*gn(7) 
      khgi = vhobr/(4.d0)*aerice(ij,ik)*gi(7)

      khg(7,ij,ik)    = khg(7,ij,ik)    + khgs*TPROB(itemp-119,ij,ik)
      khg(7+14,ij,ik) = khg(7+14,ij,ik) + khgn*TPROB(itemp-119,ij,ik)
      khg(7+28,ij,ik) = khg(7+28,ij,ik) + khgi*TPROB(itemp-119,ij,ik)


      ggsfc(7,ij,ik) = gs(7) 
      ggsfc(7+14,ij,ik) = gn(7) 
      ggsfc(7+28,ij,ik) = gi(7) 


C
C   Now set up 7 new HET reactions for JPL-2006 - fort 9fk run, set all = 0.0 below
C
C
C   WCONV261 - these reactions (8-14) are ZERO'd below (not used)
C      so comment out here to save CPU
C
C
c  NO3 + H2O -> HNO3 + OH
C
cccc      khgs = vno3/4.d0*aerosol(ij,ik)*gs(8)
cccc      khgn = vno3/4.d0*nataer(ij,ik)*gn(8)
cccc      khgi = vno3/4.d0*aerice(ij,ik)*gi(8)
cccc
cccc      khg(8,ij,ik)    = khg(8,ij,ik)    + khgs*TPROB(itemp-119,ij,ik)
cccc      khg(8+14,ij,ik) = khg(8+14,ij,ik) + khgn*TPROB(itemp-119,ij,ik)
cccc      khg(8+28,ij,ik) = khg(8+28,ij,ik) + khgi*TPROB(itemp-119,ij,ik)


cccc      ggsfc(8,ij,ik)    = gs(8)
cccc      ggsfc(8+14,ij,ik) = gn(8)
cccc      ggsfc(8+28,ij,ik) = gi(8)


C
c  Cl2 + HBr -> BrCl + HCl
C
cccc      khgs = vcl2/4.d0*aerosol(ij,ik)*gs(9)
cccc      khgn = vcl2/4.d0*nataer(ij,ik)*gn(9)
cccc      khgi = vcl2/4.d0*aerice(ij,ik)*gi(9)

cccc      khg(9,ij,ik)    = khg(9,ij,ik)    + khgs*TPROB(itemp-119,ij,ik)
cccc      khg(9+14,ij,ik) = khg(9+14,ij,ik) + khgn*TPROB(itemp-119,ij,ik)
cccc      khg(9+28,ij,ik) = khg(9+28,ij,ik) + khgi*TPROB(itemp-119,ij,ik)


cccc      ggsfc(9,ij,ik)    = gs(9)
cccc      ggsfc(9+14,ij,ik) = gn(9)
cccc      ggsfc(9+28,ij,ik) = gi(9)



c  N2O5 + HBr -> HNO3 + BrONO/BrNO2

cccc      khgs = vn2o5/(4.d0)*aerosol(ij,ik)*gs(10)
cccc      khgn = vn2o5/(4.d0)*nataer(ij,ik)*gn(10)
cccc      khgi = vn2o5/(4.d0)*aerice(ij,ik)*gi(10)

cccc      khg(10,ij,ik)    = khg(10,ij,ik)    + khgs*TPROB(itemp-119,ij,ik)
cccc      khg(10+14,ij,ik) = khg(10+14,ij,ik) + khgn*TPROB(itemp-119,ij,ik)
cccc      khg(10+28,ij,ik) = khg(10+28,ij,ik) + khgi*TPROB(itemp-119,ij,ik)


cccc      ggsfc(10,ij,ik) = gs(10) 
cccc      ggsfc(10+14,ij,ik) = gn(10) 
cccc      ggsfc(10+28,ij,ik) = gi(10) 


C
c  HOCl + HBr -> BrCl + H2O

cccc      khgs = vhocl/(4.d0)*aerosol(ij,ik)*gs(11)
cccc      khgn = vhocl/(4.d0)*nataer(ij,ik)*gn(11)
cccc      khgi = vhocl/(4.d0)*aerice(ij,ik)*gi(11)

cccc      khg(11,ij,ik)    = khg(11,ij,ik)    + khgs*TPROB(itemp-119,ij,ik)
cccc      khg(11+14,ij,ik) = khg(11+14,ij,ik) + khgn*TPROB(itemp-119,ij,ik)
cccc      khg(11+28,ij,ik) = khg(11+28,ij,ik) + khgi*TPROB(itemp-119,ij,ik)


cccc      ggsfc(11,ij,ik) = gs(11) 
cccc      ggsfc(11+14,ij,ik) = gn(11) 
cccc      ggsfc(11+28,ij,ik) = gi(11) 



c  ClONO2 + HBr -> BrCl + HNO3

cccc      khgs = vclono2/(4.d0)*aerosol(ij,ik)*gs(12) 
cccc      khgn = vclono2/(4.d0)*nataer(ij,ik)*gn(12)
cccc      khgi = vclono2/(4.d0)*aerice(ij,ik)*gi(12)

cccc      khg(12,ij,ik)    = khg(12,ij,ik)    + khgs*TPROB(itemp-119,ij,ik)
cccc      khg(12+14,ij,ik) = khg(12+14,ij,ik) + khgn*TPROB(itemp-119,ij,ik)
cccc      khg(12+28,ij,ik) = khg(12+28,ij,ik) + khgi*TPROB(itemp-119,ij,ik)


cccc      ggsfc(12,ij,ik) = gs(12)
cccc      ggsfc(12+14,ij,ik) = gn(12)
cccc      ggsfc(12+28,ij,ik) = gi(12)
                                                    


c  HOBr + HBr -> Br2 + H2O

cccc      khgs = vhobr/(4.d0)*aerosol(ij,ik)*gs(13)
cccc      khgn = vhobr/(4.d0)*nataer(ij,ik)*gn(13) 
cccc      khgi = vhobr/(4.d0)*aerice(ij,ik)*gi(13)

cccc      khg(13,ij,ik)    = khg(13,ij,ik)    + khgs*TPROB(itemp-119,ij,ik)
cccc      khg(13+14,ij,ik) = khg(13+14,ij,ik) + khgn*TPROB(itemp-119,ij,ik)
cccc      khg(13+28,ij,ik) = khg(13+28,ij,ik) + khgi*TPROB(itemp-119,ij,ik)


cccc      ggsfc(13,ij,ik) = gs(13) 
cccc      ggsfc(13+14,ij,ik) = gn(13) 
cccc      ggsfc(13+28,ij,ik) = gi(13) 



c  BrONO2 + HCl -> BrCl + HNO3

cccc      khgs = vbrono2/(4.d0)*aerosol(ij,ik)*gs(14) 
cccc      khgn = vbrono2/(4.d0)*nataer(ij,ik)*gn(14)
cccc      khgi = vbrono2/(4.d0)*aerice(ij,ik)*gi(14)

cccc      khg(14,ij,ik)    = khg(14,ij,ik)    + khgs*TPROB(itemp-119,ij,ik)
cccc      khg(14+14,ij,ik) = khg(14+14,ij,ik) + khgn*TPROB(itemp-119,ij,ik)
cccc      khg(14+28,ij,ik) = khg(14+28,ij,ik) + khgi*TPROB(itemp-119,ij,ik)


cccc      ggsfc(14,ij,ik) = gs(14)
cccc      ggsfc(14+14,ij,ik) = gn(14)
cccc      ggsfc(14+28,ij,ik) = gi(14)
   


 4000   CONTINUE



cc      if (ij .eq. 9  .and.  ik .eq. 10) then
cc         write(28,1051)
cc      write(28,7050)ij,ik,khg(7,ij,ik),gs(7),aerosol(ij,ik),vhobr, 
cc     >    aw, wtpct, xhcl, m(ij,ik), dtemp
cc 7050    format(2I4, 4x, '125', 4x, 1P9E11.3)
cc         write(28,1051)
cc      endif


cccc      if (iyr .eq. 60  .or.  iyr .eq. 65) then
ccc      if (iday360 .eq. 10) then
ccc      if (ik .eq. 11) then
ccc        write(18,103) 
ccc        write(18,102) ij, ik, khg(1,ij,ik), khg(2,ij,ik), khg(3,ij,ik), 
ccc     >    khg(4,ij,ik), khg(5,ij,ik), khg(6,ij,ik), khg(7,ij,ik)
cc
ccc        write(18,104) vclono2, c(29,ij,ik), c(15,ij,ik), 
ccc     >        aerosol(ij,ik), nataer(ij,ik), aerice(ij,ik)
cc
ccc        write(18,104) gs(1), gs(2), gs(3), gs(4), gs(5), gs(6), gs(7)
cc
cc        write(16) iday360, ij, ik, gs(1), gs(2), gs(5), gs(6), 
cc     >            c(15,ij,ik), c(29,ij,ik), m(ij,ik), dtemp
cc
cc       write(18,105) iday360, ij, ik, wtpct, gs(1), gs(2), gs(5), gs(6), 
cc     >                pph2o, pphcl, dtemp, tmin, tt
cc
cc 102    format(2I3, 4x, 1P7E10.2)
cc 104    format(10X, 1P7E10.2)
cc 103    format(5x)
cc
cc 105    format(3I3, 2x, 7E12.4, 3F8.2)
cc      ENDIF
ccc      ENDIF
ccc     ENDIF
cc
c
c      if(ij.eq.1.and.ik.eq.10)write(6,*)'in hetchem:'
c      if(ij.eq.1.and.ik.eq.10)write(6,*)khg(1,ij,ik),
c     ckhg(2,ij,ik),khg(3,ij,ik),khg(4,ij,ik),
c     ckhg(5,ij,ik),dtemp,ah2o,
c     cwtpct,gn(1),gn(2),gn(5)



C  set het reactions = 0.0  at bottom level

       if (ik .eq. 1) then 
          do 2004 igf = 1,rht$
 2004        khg(igf,ij,1) = 0.0
       endif



 2000   CONTINUE



cc  for 9FO run, use 9FN and TURN ON  NEW HET reactions that have only small impact
CC     ie, turn OFF REACTIONS that are significant (BOTH SULFATES;  ClONO2+HBr on NAT, ICE;  Cl2+HBr (Ice)
C
C   for 9GL run, turn OFF ALL new reactions from JPL-2006
C

       do 404 iihh=8,14 
           khg(iihh,ij,ik) = 0.E0
           khg(iihh+14,ij,ik) = 0.E0
 404       khg(iihh+28,ij,ik) = 0.E0

ccc           khg(12+14,ij,ik) = 0.E0

ccc           khg( 9+28,ij,ik) = 0.E0
ccc           khg(12+28,ij,ik) = 0.E0




c  load TOTAL GSFC het reaction array into the AER array, after looping through temperature PDF


       RAER(91)  = khg(1,ij,ik) 
       RAER(161) = khg(1+14,ij,ik) 
       RAER(166) = khg(1+28,ij,ik) 

ccc       RAER(166) = khg(2,ij,ik)*cn(15,ij,ik)

       RAER(93)  = khg(2,ij,ik)
       RAER(160) = khg(2+14,ij,ik)
       RAER(165) = khg(2+28,ij,ik)

ccc       RAER(164) = khg(3,ij,ik)*cn(15,ij,ik) 

       RAER(92)  = khg(3,ij,ik)
       RAER(158) = khg(3+14,ij,ik)
       RAER(163) = khg(3+28,ij,ik)

       RAER(159) = khg(4+14,ij,ik) 
       RAER(164) = khg(4+28,ij,ik) 

       RAER(90)  = khg(5,ij,ik) 
       RAER(162) = khg(5+14,ij,ik) 
       RAER(167) = khg(5+28,ij,ik) 

ccc       RAER(170) = khg(6,ij,ik)*cn(15,ij,ik)

       RAER(122) = khg(6,ij,ik)
       RAER(169) = khg(6+28,ij,ik)


       RAER(124) = khg(7,ij,ik)
       RAER(168) = khg(7+28,ij,ik)

       RAER(285) = khg(8+28,ij,ik)

       RAER(275) = khg(9+28,ij,ik)

       RAER(277) = khg(10+14,ij,ik) 
       RAER(276) = khg(10+28,ij,ik) 

       RAER(278) = khg(11+28,ij,ik) 

       RAER(279) = khg(12+14,ij,ik) 
       RAER(280) = khg(12+28,ij,ik) 

       RAER(282) = khg(13,ij,ik)
       RAER(281) = khg(13+28,ij,ik)

       RAER(284) = khg(14,ij,ik)
       RAER(283) = khg(14+28,ij,ik)


      RETURN
      END




      subroutine wtpctsulf(t,tt,pph2o,wtpct,molsulf,aw,tmin)

c subroutine to calculate weight % h2so4 from JPL 00, Table A-1, p. 58
c input variables:
c     t - temperature
c     pph2o - partial pressure of h2o, in mbar
c     
c output variables:
c     wtpct - wt % sulfate.
c     molsulf - molality of sulfate
c     aw - water activity

c internal variables:
c     svph2o - saturation vapor pressure of h2o, in mbar
c     y1,y2 -intermediate polynomials
c     a1,b1,c1,d1,a2,b2,c2,d2 - arrays storing values for aw calc
c     sc - array storing constants for h2o svp calc.

c delare and initialize variables:

      integer i
      real*8 t
      real*8 a1(3),b1(3),c1(3),d1(3),a2(3),b2(3),c2(3),d2(3)
      real*8 sc(4),aw,y1,y2,molsulf,wtpct,tt,tmin,pph2o

      data sc /18.452406985,-3505.1578807,-330918.55082,12725068.262/

      data a1 /12.37208932,11.820654354,-180.06541028/
      data b1 /-0.16125516114,-0.20786404244,-0.38601102592/
      data c1 /-30.490657554,-4.807306373,-93.317846778/
      data d1 /-2.1133114241,-5.1727540348,273.88132245/
      data a2 /13.455394705,12.891938068,-176.95814097/
      data b2 /-0.1921312255,-0.23233847708,-0.36257048154/
      data c2 /-34.285174607,-6.4261237757,-90.469744201/
      data d2 /-1.7620073078,-4.9005471319,267.45509988/

c if t is lower than minimum acceptable for given h2o
c partial pressure, constrain.  Note the 3rd order polynomial
c function being used here was determined by first calculating the
c wt pct h2so4 for 170 < t < 230 and 5.5e-5 < pph2o < 5.5e-4 mbar
c (1 - 10 ppmv at 50 mbar), and finding the minimum temperature
c as a function of pph2o before the wt pct calculation went negative.
c I added the 1. as a kluge - with my first try, there were still
c a few negative wt pcts showing up.

      tmin=1.+172.369141+7.66720e4*pph2o-1.51781376e8*pph2o**2
     c     +1.2052752e11*pph2o**3
   

      tmin = DMIN1(195.d0,tmin)
      tt = DMAX1(t,tmin)


c calculate saturation vapor pressure of h2o

      svph2o=DEXP(sc(1)+sc(2)/tt+sc(3)/tt**2+sc(4)/tt**3)  

c calculate water activity

      aw=pph2o/svph2o
      aw=DMIN1(1.0d0,aw)
      aw=DMAX1(1.d-3,aw)

c calculate sulfate molality

      if(aw.le.0.05)i=1
      if(aw.gt.0.05.and.aw.lt.0.85)i=2
      if(aw.ge.0.85)i=3

      y1=a1(i)*aw**b1(i)+c1(i)*aw+d1(i)
      y2=a2(i)*aw**b2(i)+c2(i)*aw+d2(i)
 
      molsulf=y1+((tt-190.d0)*(y2-y1))/70.d0

      wtpct=9800.d0*molsulf/(98.d0*molsulf+1000.d0)
      return
      end


      function gamhydsulf(t,wtpct)

c function to calculate sticking coefficient for N2O5+H2O on
c sulfate aerosol reaction - according to Note 17 of JPL 00   -- IDENTICAL to JPL-06 (note 27)

c input variables:
c     t - temperature
c     wtpct - wt pct of h2so4, as from function wtpctsulf
c internal variables:
c     kc0,kc1,kc2 - arrays of coefficients for calculating k0,k1,k2
c     k0,k1,k2 - functions of wt pct for calculating gamma
c output variable:
c     gamhydsulf- sticking coeff. for n2o5+h2o on sulfate

c declare and initialize variables:
      
      real*8 kc0(4),kc1(4),kc2(4),k0,k1,k2,t,wtpct,gamhydsulf

      data kc0 /-25.5265,-0.133188,0.0093084,-9.0194d-5/
      data kc1 /9283.76,115.345,-5.19258,0.0483464/
      data kc2 /-851801.,-22191.2,766.916,-6.85427/

c calculate k0,k1,k2

      k0=kc0(1)+kc0(2)*wtpct+kc0(3)*wtpct**2+kc0(4)*wtpct**3
      k1=kc1(1)+kc1(2)*wtpct+kc1(3)*wtpct**2+kc1(4)*wtpct**3
      k2=kc2(1)+kc2(2)*wtpct+kc2(3)*wtpct**2+kc2(4)*wtpct**3

c calculate gamma:

      gamhydsulf = DEXP(k0+k1/t+k2/t**2)

      return
      end


     
      subroutine gamclonitsulf(ij,ik,t,molsulf,wtpct,pphcl,ppclono2,
     c                         aw,droprad,gamclono2hcl,gamclono2h2o,
     c                         gamhoclhcl,gHOBr_HCl,iday360)

c subroutine to calculate sticking coefficients for ClONO2+H2O,
c ClONO2 + HCl, and HOCl+HCl on sulfate. TABLES A-2, A-3, A-4 in JPL-00 (pp 59-62)

c input variables
c     t - temperature
c     molsulf - molality of sulfate
c     wtpct - weight percent of sulfate
c     pphcl - partial pressure of hcl (mbar)
c     ppclono2 - partial pressure of clono2 (mbar)
c     aw - water activity
c     droprad - droplet radius (cm)

c output variables
c     gamclono2hcl
c     gamclono2h2o
c     gamhoclhcl

c internal variables
c     setch - Setchenow coefficient (M^-1)
c     z1,z2,z3 - to calculate rhoh2so4
c     rhosulf - solution density (g/cm^3)
c     molarsulf - molarity of sufate
c     henhcl - Henry's law coeff for HCl
c     molarhcl - molarity of hcl
c     to -to calculate viscosity
c     visc - viscosity of sulfate solution
c     diffclono2 - diffusion coeff. for clono2
c     ah - acid activity
c     kh2o,kh,khyd - surface rxn rates?
c     henclono2 - Henry's law coeff for ClONO2
c     thermclono2 - thermal speed of clono2
c     fclono2 -
c     gamrclono2
c     gambhcl
c     gams
c     fhcl - 
c     denom - intermediate to fhcl - no phys. interp
c     gamsp
c     gambhclp
c     gamb

c declare and initialize variables:


      real*8 t,molsulf,wtpct,pphcl,ppclono2,aw,droprad,gamclono2hcl,
     c       gamclono2h2o,gamhoclhcl,
     c       setchhocl,diffhocl,khocl,henhocl,thermhocl,
     c       gamrhocl,lhocl,fhocl

      real*8 z1,z2,z3,rhosulf,kh2o,kh,khyd,henclono2,thermclono2,
     >       molarsulf,molfrac,setch,henhcl,molarhcl,to,visc,diffclono2,
     >       ah,khcl,gambh2o,lclono2,fclono2,gamrclono2,
     >       gambhcl,gams, fhcl,denom,gamsp,gambhclp,gamb,gamclono2
  
      REAL*8 c_HOBr, SHOBr, HHOBr, DHOBr, kHOBr_HCl, GHOBrrxn, 
     >       lHOBr, fHOBr, gHOBr_HCl


c h2so4 solution density:

      z1=0.12364-5.6d-7*t**2
      z2=-0.02954+1.814d-7*t**2
      z3=2.343d-3-1.487d-6*t-1.324d-8*t**2
      rhosulf=1.d0+z1*molsulf+z2*molsulf**1.5+z3*molsulf**2

c molarity of sulfate:

      molarsulf=rhosulf*wtpct/9.8d0

c mole fraction:

      molfrac=wtpct/(wtpct+(100.d0-wtpct)*98.d0/18.d0)

c Setchenow coefficient:

      setch= 0.306+24.d0/t

c Henry's law coefficient for HCl:

      henhcl=(0.094-0.61*molfrac+1.2*molfrac**2)*
     c       DEXP(-8.68+(8515.-10718.*molfrac**0.7)/t)

c Molarity of HCl (convert hcl part pressure to atm from mb)

      molarhcl=henhcl*pphcl/1013.d0


c viscosity of sulfate:

      to=144.11+0.166*wtpct-0.015*wtpct**2+2.18d-4*wtpct**3

      visc=(169.5+5.18*wtpct-0.0825*wtpct**2+3.27d-3*wtpct**3)*
     c     t**(-1.43)*DEXP(448./(t-to))

c diffusion coefficient for clono2

      diffclono2=5.d-8*t/visc


c acid activity

      ah=DEXP(60.51-0.095*wtpct+0.0077*wtpct**2-1.61d-5*wtpct**3
     c       -(1.76+2.52d-4*wtpct**2)*t**0.5
     c       +(253.05*wtpct**0.076-805.89)/t**0.5)

c surface rxn rates:

      kh2o=1.95d10*DEXP(-2800./t)
      kh=1.22d12*DEXP(-6200./t)
      khyd=kh2o*aw+kh*ah*aw
      khcl=7.9d11*ah*diffclono2*molarhcl

c Henry's law coefficients for ClONO2

      henclono2=1.6d-6*DEXP(4710./t)*DEXP(-setch*molarsulf)

c Thermal speed of ClONO2:

      thermclono2=1474.d0*t**0.5
c 

      gambh2o=4.d0*0.082*henclono2*t*(diffclono2*khyd)**0.5
     c        /thermclono2


c reacto-diffusive length:

      lclono2=(diffclono2/(khyd+khcl))**0.5
c 

      fclono2=1.d0/DTANH(droprad/lclono2)-lclono2/droprad

c

      gamrclono2=fclono2*gambh2o*(1.d0+khcl/khyd)**0.5

      gambhcl=gamrclono2*khcl/(khcl+khyd)

      gams=66.12*henclono2*molarhcl*DEXP(-1374./t)

      denom=1.d0+0.612*(gams+gambhcl)*ppclono2/pphcl
      fhcl=1.d0/denom

      gamsp=fhcl*gams
      gambhclp=fhcl*gambhcl

      gamb=gambhclp+gamrclono2*khyd/(khcl+khyd)


c Finally, the sticking coefficients:

      gamclono2=1.d0/(1.d0+1./(gamsp+gamb))
      gamclono2hcl=gamclono2*(gamsp+gambhclp)/(gamsp+gamb)
      gamclono2h2o=gamclono2-gamclono2hcl

cccc      if (ij .eq. 1 .or. ij .eq. 4 .or. ij .eq. 9 .or.  ij .eq. 13) then
cc      if (iday360 .eq. 10) then
cc      if (ij .eq. 9) then
cc        write(17,103) 
cc        write(17,102) ij, ik, t,molsulf,wtpct,pphcl,ppclono2,aw,droprad
cc        write(17,104) z1,z2,z3,rhosulf,molarsulf,molfrac,setch, 
cc     >               henhcl,molarhcl,to,visc,diffclono2
cc        write(17,104) ah,kh2o,kh,khyd,khcl,henclono2,
cc     >      thermclono2,gambh2o,lclono2,fclono2,gamrclono2
cc        write(17,104) gambhcl,gams,denom,fhcl,gamsp,gambhclp,gamb,
cc     >               gamclono2, gamclono2hcl, gamclono2h2o
cc 102    format(2I3, 4x, 1P11E10.2)
cc 104    format(1P12E10.2)
cc 103    format(5x)
cc      ENDIF
cc      ENDIF
cc

c OK, now do HOCl + HCL:

      setchhocl=0.0776d0+59.18/t
      diffhocl=6.4d-8*t/visc
      khocl=1.25d9*ah*diffhocl*molarhcl
      henhocl=1.91d-6*DEXP(5862.4/t)*DEXP(-setchhocl*molarsulf)
      thermhocl=2009.d0*t**0.5

      gamrhocl=4.d0*0.082*henhocl*t*(diffhocl*khocl)**0.5/thermhocl
      lhocl=(diffhocl/khocl)**0.5

      fhocl=1./DTANH(droprad/lhocl)-lhocl/droprad

      denom=1.d0+1.d0/(fhocl*gamrhocl*fhcl)
      gamhoclhcl=1.d0/denom

c
c  //Table 5.  HOBr + HCl - taken from AER code - but currently set gs(7) = 0.0 above
C  JPL-02, JPL-06 do not have definitive recommendations - further data is needed to model this reaction

      c_HOBr    =  1478.025*SQRT(t)
      SHOBr     =  0.0776+59.18/T  
      HHOBr     =  30.
      DHOBr     =  1.E-8
      kHOBr_HCl =  1.4E5*henhcl*pphcl/1013.25                            ! HCl partial press in atm
      GHOBrrxn  =  4.*HHOBr*0.082*T*SQRT(DHOBr*kHOBr_HCl)/c_HOBr     
      lHOBr	=  SQRT(DHOBr/kHOBr_HCl)
      fHOBr     =  1./TANH(droprad/lHOBr)- lHOBr/droprad         
      gHOBr_HCl =  1./(1.+1./(fHOBr*GHOBrrxn*FHCl)) 

      return
      end
