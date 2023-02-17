c******************************************************************************
c                             subroutine reaction                             *
c                                                                             *
c  the subroutine uses the data read in reacread to calculate the reaction    *
c  rates in each case.                                                        *
c                                                                             *
c******************************************************************************
C
C
C  THIS IS for JPL-2006 - updated June 2007 (EF)
C

       SUBROUTINE REACTION


       include "com2d.h"

c
c   K-rates now calculated in REAL*8 w/ AER chemistry, (EF, Feb. 2003)
C

       REAL*8 ts, xx, yyz, k01, k2, k3, kir, k00, kinf


       SAVE


C  calculate reaction rates for base reaction set.
C
c   ibdy(i) = 0  -> calculates two-body rates with NO temperature dependence
C
c   ibdy(i) = 1  -> calculates two-body rates WITH temperature dependence
C
C   ibdy = 3 is three-body reaction
C
c   if a special calculation is required, then ibdy(i) = 4 or 5
C
C   for BASE9BD - include longitudinal temperature PDF 
C    - loop through from 120-330 K, TPROB(211,L$,Z$), ITP1(L$,Z$), ITP2(L$,Z$) 
C
C      REAL*8 K(RB$=284,L$,Z$)
C


      do 100  ik=1,Z$
         do 110  ij=1,L$


         do 200 irz=1,rb$
 200            k(irz,ij,ik) = 0.0



      do 1000 ir=1,rb$


C           2-body reaction, NO temperature dependence

        if (ibdy(ir) .eq. 0) k(ir,ij,ik) = k0(ir)


C                                 2-body reaction with temperature dependence
        if (ibdy(ir) .eq. 1) then
           do 4000 itemp=itp1(ij,ik)+119, itp2(ij,ik)+119
             ts = DBLE(itemp)                               ! get real*8 value
                                             !if using zonal mean temps, use TEMP(L$,Z$X)-REAL*8
             if (IZONAVGT .eq. 1) ts = TEMP(ij,ik)

             kir = k0(ir)*DEXP(-e(ir)/ts)
             k(ir,ij,ik) = k(ir,ij,ik) + kir*TPROB(itemp-119,ij,ik)
 4000      CONTINUE
        endif


c   Note: for 3-body reactions, now include M in rate to be consistent with AER

        if(ibdy(ir) .eq. 3) then
           do 4002 itemp=itp1(ij,ik)+119, itp2(ij,ik)+119
             ts = DBLE(itemp)   ! get real*8 value

             if (IZONAVGT .eq. 1) ts = TEMP(ij,ik)

             xx=k0(ir)*(ts/300.D0)**(-e(ir))
             yyz=khi(ir)*(ts/300.D0)**(-ehi(ir))

             kir = m(ij,ik)*(xx/(1.D0 + xx*m(ij,ik)/yyz)*
     >             .6D0**(1.D0/(1.D0 + (DLOG10(xx*m(ij,ik)/yyz))**2)))

             k(ir,ij,ik) = k(ir,ij,ik) + kir*TPROB(itemp-119,ij,ik)
 4002      CONTINUE
         ENDIF

 1000       continue



c special cases in base reaction set:

c old reaction k(8): o(1d)+m->o(3p)+m - now done the usual way, consistent with AER model in K18, K19


c reaction 1: O+O2+M -> O3+M  -  now include M in the reaction rate to be consistent w/ AER

      do 4006 itemp=itp1(ij,ik)+119, itp2(ij,ik)+119
        ts = DBLE(itemp)                             ! get real*8 value
        if (IZONAVGT .eq. 1) ts = TEMP(ij,ik)

        xx = k0(1)*(ts/300.D0)**(-e(1))
       k(1,ij,ik) = k(1,ij,ik) + xx*m(ij,ik)*TPROB(itemp-119,ij,ik)
 4006  CONTINUE



c reaction 82: O(1D)+N2+M->N2O+M 
C  - seems like we DO include M here in reaction rate to be consistent w/ AER JACOB (similar to O+O2+M above)
C      (ie, there's NO M in JACOB for this reaction - R212)

        do 4012 itemp=itp1(ij,ik)+119, itp2(ij,ik)+119
          ts = DBLE(itemp)                   ! get real*8 value
          if (IZONAVGT .eq. 1) ts = TEMP(ij,ik)

          xx = k0(82)*(ts/300.D0)**(-e(82))
          k(82,ij,ik) = k(82,ij,ik) + xx*m(ij,ik)*TPROB(itemp-119,ij,ik)
 4012   CONTINUE




C  reaction 14:  OH+CH4

      do 4001 itemp=itp1(ij,ik)+119, itp2(ij,ik)+119
        ts = DBLE(itemp)                            ! get real*8 value
        if (IZONAVGT .eq. 1) ts = TEMP(ij,ik)

        kir = k0(14)*DEXP(-e(14)/ts)*ts**(0.667D0)
        k(14,ij,ik) = k(14,ij,ik) + kir*TPROB(itemp-119,ij,ik)
 4001 CONTINUE




c the 6-parameter rate constant (rxn #37)

      do 4003 itemp=itp1(ij,ik)+119, itp2(ij,ik)+119
        ts = DBLE(itemp)                            ! get real*8 value
        if (IZONAVGT .eq. 1) ts = TEMP(ij,ik)

  	k01 = khno3(1)*DEXP(ehno3(1)/ts)
	k2 = khno3(2)*DEXP(ehno3(2)/ts)
	k3 = khno3(3)*DEXP(ehno3(3)/ts)

        kir = k01 + k3*m(ij,ik)/(1.D0 + k3*m(ij,ik)/k2) 
        k(37,ij,ik) = k(37,ij,ik) + kir*TPROB(itemp-119,ij,ik)
 4003  CONTINUE

ccc  this is the same 
ccc            k(37,ij,ik)=khno3(1)*exp(ehno3(1)/ts)+
ccc     c           khno3(3)*exp(ehno3(3)/ts)*m(ij,ik)/
ccc     c           (1.+khno3(3)*exp(ehno3(3)/ts)*m(ij,ik)/
ccc    c           (khno3(2)*exp(ehno3(2)/ts)))


c now, do the reactions for which we have forward rate constants 
c and equilibrium parameters 
C
C  NOTE: 
C    the forward rate constants (k46, k34, k102, k190, K271) are 3-body reactions, 
C    which now are multiplied by M (above) to be consistent w/ AER model
C    so k31, k70, k106, k191, K272 here also include M and are consistent w/ AER
  

        do 4008 itemp=itp1(ij,ik)+119, itp2(ij,ik)+119
          ts = DBLE(itemp)                             ! get real*8 value
          if (IZONAVGT .eq. 1) ts = TEMP(ij,ik)

          kir = k(46,ij,ik)/(k0(31)*DEXP(e(31)/ts))
          k(31,ij,ik) = k(31,ij,ik) + kir*TPROB(itemp-119,ij,ik)
 4008   CONTINUE


        do 4009 itemp=itp1(ij,ik)+119, itp2(ij,ik)+119
          ts = DBLE(itemp)                             ! get real*8 value
          if (IZONAVGT .eq. 1) ts = TEMP(ij,ik)

          kir = k(34,ij,ik)/(k0(70)*DEXP(e(70)/ts))
          k(70,ij,ik) = k(70,ij,ik) + kir*TPROB(itemp-119,ij,ik)
 4009   CONTINUE


        do 4010 itemp=itp1(ij,ik)+119, itp2(ij,ik)+119
          ts = DBLE(itemp)                             ! get real*8 value
          if (IZONAVGT .eq. 1) ts = TEMP(ij,ik)

          kir = k(102,ij,ik)/(k0(106)*DEXP(e(106)/ts))
          k(106,ij,ik) = k(106,ij,ik) + kir*TPROB(itemp-119,ij,ik)
 4010   CONTINUE


        do 4017 itemp=itp1(ij,ik)+119, itp2(ij,ik)+119
          ts = DBLE(itemp)                             ! get real*8 value
          if (IZONAVGT .eq. 1) ts = TEMP(ij,ik)

          kir = k(190,ij,ik)/(k0(191)*DEXP(e(191)/ts))
          k(191,ij,ik) = k(191,ij,ik) + kir*TPROB(itemp-119,ij,ik)
 4017   CONTINUE


        do 4027 itemp=itp1(ij,ik)+119, itp2(ij,ik)+119
          ts = DBLE(itemp)                             ! get real*8 value
          if (IZONAVGT .eq. 1) ts = TEMP(ij,ik)

          kir = k(271,ij,ik)/(k0(272)*DEXP(e(272)/ts))
          k(272,ij,ik) = k(272,ij,ik) + kir*TPROB(itemp-119,ij,ik)
 4027   CONTINUE



C reaction 64 and 65:  HO2+HO2 and HO2+HO2+M = H2O2+O2   
C    these are combined into K(64) to be consistent w/ AER model (AER rate #17, type 7)

        do 4011 itemp=itp1(ij,ik)+119, itp2(ij,ik)+119
          ts = DBLE(itemp)                                 ! get real*8 value
          if (IZONAVGT .eq. 1) ts = TEMP(ij,ik)

          kir = k0(64)*DEXP(-e(64)/ts) + 
     >         (k0(65)*DEXP(-e(65)/ts))*m(ij,ik)

          k(64,ij,ik) = k(64,ij,ik) + kir*TPROB(itemp-119,ij,ik)
 4011   CONTINUE




c  special case for OH+CO-> H+CO2 channel, K(36) - see JPL-2006, pp. 2-2, 2-4, 2-10 - THIS IS READY TO GO!!!!

           do 4022 itemp=itp1(ij,ik)+119, itp2(ij,ik)+119
             ts = DBLE(itemp)                             ! get real*8 value
             if (IZONAVGT .eq. 1) ts = TEMP(ij,ik)

             k00 =  k0(36)*(ts/300.D0)**(-e(36))
             kinf = khi(36)*(ts/300.D0)**(-ehi(36))

             kir = (k00/(1.D0 + (k00/(kinf/m(ij,ik)) ) ) )*0.6D0**
     >            ( 1.D0/(1.D0 + (DLOG10( k00/(kinf/m(ij,ik)) ) )**2 ) )

             k(36,ij,ik) = k(36,ij,ik) + kir*TPROB(itemp-119,ij,ik)
 4022      CONTINUE




C  SPECIAL CASES FOR REACTION TYPE 5 - radical prod/loss NOW included in AER model for MESOSPHERE - May 03


c reaction 90: O+O+M->O2+M  -  DON'T include M here in reaction rate to be constitent w/ AER JACOB

        do 4014 itemp=itp1(ij,ik)+119, itp2(ij,ik)+119
          ts = DBLE(itemp)                              ! get real*8 value
          if (IZONAVGT .eq. 1) ts = TEMP(ij,ik)

          kir = k0(90)*DEXP(-e(90)/ts)
          k(90,ij,ik) = k(90,ij,ik) + kir*TPROB(itemp-119,ij,ik)
 4014   CONTINUE



c reaction 111:  CO+O+M -> CO2+M  - DO NOT include M here in reaction rate to be constitent w/ AER JACOB

            k(111,ij,ik) = 2.0D-37


C Note: for K(140) ClO+O2(1D)+M->ClO3+M,   M is NOT included in the reaction rate as per the AER code (JACOB)



c now do the working reaction rate set: - DON'T DO SINCE IT'S NOT USED ANYWHERE
cef
cef            do 1001 ir=1,rw$
cef
cef               if(ibdyw(ir).eq.0)then
cef
cef                  kw(ir,ij,ik)=k0w(ir)*exp(-ew(ir)/ts)
cef
cef               endif
cef
cef               if(ibdyw(ir).eq.1)then
cef
cef                  xx=k0w(ir)*(ts/300.)**(-ew(ir))
cef                  yyz=khiw(ir)*(ts/300.)**(-ehiw(ir))
cef                  kw(ir,ij,ik)=xx/(1.+xx*m(ij,ik)/yyz)*
cef     c                 .6**(1./(1.+(log10(xx*m(ij,ik)/yyz))**2))
cef
cef               endif
cef
cef 1001       continue
cef
c special cases for the working reaction - rate set:


C
C  for 9HB (base JPL10), redo O1D reaction/quenching reactions here, 
C                                      based on Jim's spreadsheet

C  CFCl3
            k(83,ij,ik)  = .88*23.D-11
            k(206,ij,ik) = .12*23.D-11
C  CF2Cl2
            k(80,ij,ik)  = .86*14.D-11
            k(207,ij,ik) = .14*14.D-11
C  CCl4
            k(54,ij,ik)  = .9*33.D-11
            k(226,ij,ik) = .1*33.D-11
C  CHClF2
            k(125,ij,ik) = .7*10.D-11
            k(205,ij,ik) = .3*10.D-11
C  CBrClF2
            k(112,ij,ik) = .65*15.D-11
            k(208,ij,ik) = .35*15.D-11
C  CBrF3
            k(113,ij,ik) = .4*10.D-11
            k(209,ij,ik) = .6*10.D-11
C  CH2FCF3
            k(215,ij,ik) = .35*4.9D-11
            k(216,ij,ik) = .65*4.9D-11
C  CF3CH3
            k(224,ij,ik) = .8*4.4D-11
            k(225,ij,ik) = .2*4.4D-11
C  CHF3
            k(203,ij,ik) = .15*.91D-11
            k(204,ij,ik) = .85*.91D-11
C  C2Cl2F4
            k(100,ij,ik) = .75*13.D-11
            k(221,ij,ik) = .25*13.D-11
C  CH3CCL2F
            k(114,ij,ik) = .7*26.D-11
            k(213,ij,ik) = .3*26.D-11
C  CH3CCLF2
            k(117,ij,ik) = .75*22.D-11
            k(214,ij,ik) = .25*22.D-11
C  CBr2F2
            k(174,ij,ik) = .45*22.D-11
            k(210,ij,ik) = .55*22.D-11
C  C2Br2F4
            k(123,ij,ik) = .75*16.D-11
            k(223,ij,ik) = .25*16.D-11
C  CH2F2
            k(201,ij,ik) = .3*5.1D-11
            k(202,ij,ik) = .7*5.1D-11
C  CHF2CF3
            k(218,ij,ik) = .75*12.D-11
            k(219,ij,ik) = .25*12.D-11
C  CH3CHF2
            k(211,ij,ik) = .55*17.5D-11
            k(212,ij,ik) = .45*17.5D-11


 110     continue
 100  continue


c for testing, write out CO+OH reaction stuff:

ccccccc           if (iday360 .le. 2) then
ccccccc              write(557,1051)
ccccccc              write(557,1057) iday360, k0(36), e(36), khi(36), ehi(36)
ccccccc              write(557,1058) iday360, k0(57), e(57), khi(57), ehi(57)
ccccccc              write(557,1051)
ccccccc              write(557,1059) k(36,22,5), k(57,22,5), k(36,22,35), 
ccccccc     >                        k(57,22,35), k(36,22,70), k(57,22,70)  
ccccccc              write(557,1051)
ccccccc              write(557,1051)
ccccccc          endif 
ccccccc
ccccccc 1051    format(4x)
ccccccc 1057  format(1x, 'iday360 =', I4, 4x,'K36 INput : ', 1P4E11.3)
ccccccc 1058  format(1x, 'iday360 =', I4, 4x,'K57 INput : ', 1P4E11.3)
ccccccc 1059  format(1x, 'K36, K57 = ', 1P6E10.2)

      return
      end
