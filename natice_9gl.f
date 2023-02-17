c

      SUBROUTINE NATICE(idor)


C  *************************************************************************************************
C
C   routine to compute the NAT and ICE aerosols, which was the first part of David's HETCHEM routine
C      this is CALLED ONCE PER DAY - the HET REACTIONS ARE NOW COMPUTED WITHIN THE DIURNAL LOOP
C
C      adapted from  hetchem_9aa_jpl00.f  - EF 3/12/03
C
C  *************************************************************************************************
C

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

c     common blocks:

      include 'com2d.h'

c     dimension internal arrays, define variables
c     note: this routine will not work without modification if the
c     latitude resolution of the model is changed from 10 degree bins.

      integer ijmax,ikmax,pscnaton,psciceon,
     c        nosedice,nosednat,imaxiceht(L$),imaxnatht(L$),
     c        iminiceht(L$),iminnatht(L$),init, NKR

c these store temp prob distributions
      real hno3sol(L$,Z$),hno3solnew(L$,Z$),hno3gas(L$,Z$)
      real h2osol(L$,Z$),h2osolnew(L$,Z$),h2ogas(L$,Z$)
      real fice(L$,Z$),fnat(L$,Z$),vfallnat(L$,Z$),
     >     vfallice(L$,Z$),fnucnat(L$,Z$),fnucice(L$,Z$)
      real h2ost(L$,Z$), tprobnow(L$,Z$,211), totna(L$,Z$)

cccc      real tprob(12,L$,Z$,131), tprobr(12,18,30,131)


c     initialization:

      data init /0/

      SAVE                         ! save all internal variables


C  First store in initial H2O, which is used at the end of HETCHEM 
C      to compute H2O loss due to ice fallout for output - EF 11/97

         do 225 ik=1,Z$
         do 225 ij=1,L$
 225        h2ost(ij,ik) = cn(15,ij,ik)


C  ONLY go up to 40 km NOW w/ high res model - ik40 in COMMON - EF, 7/07

c     read in temperature probability density distributions:
c     also set switch to turn pscs on or off:

      if(init.eq.0)then
         init=1
         ijmax=L$
         ikmax=IK40
c to turn on PSCs, set pscnaton and psciceon = 1
         pscnaton=1
         psciceon=1
c to turn off sedimentation, set nosednat and nosedice = 1
         nosednat=0
         nosedice=0

         if(pscnaton.eq.1.or.psciceon.eq.1)then
cc            open(87,file='probdist3.dat',
cc     c                                  form='unformatted')
cc            read(87)tprobr
cc            close(87)
cc
C  tprobr(12,18,30,131) - January 85S doesn't look right, just interpolate between Dec-Feb.
C
cc      do 887 ittt=1,131
cc      do 887 ik=1,30
cc 887  tprobr(1,1,ik,ittt)= (tprobr(12,1,ik,ittt)+tprobr(2,1,ik,ittt))/2.


C now load in probability distribs. into new grid (bi-linear interpolation didn't work)
C   LATIN(18), ZZ30(30), LAT(L$), zalt(Z$X), tprobr(12,18,30,131), tprob(12,L$,Z$,131)
C

cc        do 115 ik=1,ik40
cc           do 117 ik1=1,30
cc             if (zalt(ik) .ge. zz30(ik1)-1.9908/2.  .and.  
cc     >           zalt(ik) .le. zz30(ik1)+1.9908/2.) iko=ik1
cc 117       CONTINUE
cc
cc        do 120 ij=1,L$
cc           do 122 ij1=1,18
cc             if (lat(ij) .ge. latin(ij1)-5.  .and.  
cc     >           lat(ij) .le. latin(ij1)+5.) ijo=ij1
cc 122       CONTINUE
cccc
cc           do 110 ittt=1,131
cc           do 110 imm=1,12
cc 110         tprob(imm,ij,ik,ittt) = tprobr(imm,ijo,iko,ittt)
cc
cc 120     CONTINUE
cc 115     CONTINUE

         endif
      endif


      if(pscnaton.ne.1.and.psciceon.ne.1)go to 1500

c get correct probability distribution for day:

cc      imon=(iday360-1)/30+1
cc      imonp=imon+1
cc      if(imonp.eq.13)imonp=1
c idom = integer day of month
cc      idom=mod(iday360-1,30)
cc      fract=idom*1./30.

c  load in proper PDF for day:  tprob_psc(211,L$,Z$), tprobnow(L$,Z$,211) 
C   just start at 170K, as David has done, since there are no points colder than 170K below 60 km
C                                      ADJUSTED to start at 120K for coupled model adjusted PDFs
      do 33 ij=1,ijmax
         do 33 ik=1,ikmax
            do 33 l=1,211

cc               tprobnow(ij,ik,l)=tprobr(imon,ij,ik,l)+fract*
cc     c         (tprobr(imonp,ij,ik,l)-tprobr(imon,ij,ik,l))

cc             tprobnow(ij,ik,l)=tprob(imon,ij,ik,l)+fract*
cc   c         (tprob(imonp,ij,ik,l)-tprob(imon,ij,ik,l))

cccc              tprobnow(ij,ik,l) = tprob_psc(l+50,ij,ik)

               tprobnow(ij,ik,l) = tprob_psc(l,ij,ik)
 33   continue


      fluxcorr=1.0        ! allows to adjust for lognormal mass flux

      delgrid=delzc*100.  ! vertical grid resolution - use delzc (meters), convert to centimeters 
cccccccccef      delgrid=2.e5        ! vertical grid resolution, centimeters 

      natsatrat=10.       ! saturation ratio for nat aerosols
      satratice=1.4       ! saturation ratio for ice aerosols
c      natsatrat=1.
c      satratice=1.
      rnat=1.e-6          ! mode radius of nat, meters 
      rice=1.e-5          ! mode radius of ice, in meters


c Code to calculate nat and ice fall velocity (rnat and rice are the
c radiuses in meters)

c constants:

      a=1.249
      b=0.42
      cc=0.87
      bet=1.458e-6
      s=110.4
      rhonat=1.6e3
      rhoice=1.e3
      sigsq=1.3323e-19
      grav=9.8
      
      vfnatmax=0.       ! maximum fall velocity for nats
      vficemax=0.       ! maximum fall velocity for ice


      do 2000 ij=1,ijmax

         imaxnatht(ij)=0.
         iminnatht(ij)=ikmax
         imaxiceht(ij)=0.
         iminiceht(ij)=ikmax

         do 2000 ik=1,ikmax

c mean free path:

            mfp=.22508/(sigsq*m(ij,ik)*1.e6)

c dynamic viscosity:

            dynvis=bet*temp(ij,ik)**1.5/(temp(ij,ik)+s)

c ice fall velocity (in cm/s):

            vfallice(ij,ik)=0.2222*rhoice*rice*rice*grav/dynvis*100.*
     c                    (1.+ mfp/rice*(a+b*exp(-cc*rice/mfp)))

            if(vfallice(ij,ik).gt.vficemax.and.fice(ij,ik).gt.0.)then
               vficemax=vfallice(ij,ik)
               if(ik.gt.imaxiceht(ij))imaxiceht(ij)=ik
               if(ik.lt.iminiceht(ij))iminiceht(ij)=ik
            endif

            if(nosedice.eq.1)vfallice(ij,ik)=0.

c nat fall velocity (in cm/s):

            vfallnat(ij,ik)=0.2222*rhonat*rnat*rnat*grav/dynvis*100.*
     c                    (1.+ mfp/rnat*(a+b*exp(-cc*rnat/mfp)))

            if(fnat(ij,ik).gt.0.)then
 
               if(fice(ij,ik).eq.0..and.vfallnat(ij,ik).gt.vfnatmax)
     c            vfnatmax=vfallnat(ij,ik)

               if(fice(ij,ik).gt.0..and.fice(ij,ik).le.fnat(ij,ik).and.
     c         (fice(ij,ik)/fnat(ij,ik)*vfallice(ij,ik)
     c         +(1.-fice(ij,ik)/fnat(ij,ik))*vfallnat(ij,ik)).gt.
     c         vfnatmax)
     c         vfnatmax= fice(ij,ik)/fnat(ij,ik)*vfallice(ij,ik)
     c         +(1.-fice(ij,ik)/fnat(ij,ik))*vfallnat(ij,ik)

               if(fice(ij,ik).gt.fnat(ij,ik).and.vfallice(ij,ik).gt.
     c         vfnatmax)vfnatmax=vfallice(ij,ik)

               if(ik.gt.imaxnatht(ij))imaxnatht(ij)=ik
               if(ik.lt.iminnatht(ij))iminnatht(ij)=ik

            endif

            if(nosednat.eq.1)vfallnat(ij,ik)=0.

c            if(j.eq.9)write(6,*)k,vfallnat(ij,ik),vfallice(ij,ik)

 2000 continue

c      write(6,*)vficemax,vfnatmax
c      write(6,2128)(imaxnatht(ii),ii=1,L$)
c      write(6,2128)(iminnatht(ii),ii=1,L$)
c      write(6,2128)(imaxiceht(ii),ii=1,L$)
c      write(6,2128)(iminiceht(ii),ii=1,L$)

 2128 format(18i3)

c     ICE particle size particulars:  The variable convice converts
c     the number density of solid phase H2O into a surface area
c     density, assuming:

      rgice=1.e2*rice
      sigice=1.8
      gpmolice=2.991e-23
      densice=1.

      convice=gpmolice/densice*3./rgice*
     c                exp(-5.*alog(sigice)*alog(sigice)/2.)

c      write(6,*)convice

c     Sediment solid phase H2O:  
c     H2O is sedimented in same way as NAT aerosol (see below),  h2osol(L$,Z$)

         do 100 ij=1,ijmax
               h2osol(ij,IK40+1)=0.           ! add an extra level (=0) to accomodate code below, line 360
            do 100 ik=1,ikmax
               h2osol(ij,ik)=cn(67,ij,ik)
               if(h2osol(ij,ik).le.0.)h2osol(ij,ik)=0.
               h2ogas(ij,ik)=cn(15,ij,ik)
               if(h2ogas(ij,ik).le.0.)then
                  h2ogas(ij,ik)=1.e-9*m(ij,ik)
                  write(6,*)'h2o gas less than 0 at:', ij,ik
               endif
 100     continue

      if(psciceon.eq.1)then

c max dt for sedimentation (depends on max fall velocity):

         dtsedice=delgrid/(fluxcorr*vficemax)

         if(dt/dtsedice.le.1.)then
            dtsedice=dt             ! no need to change timestep
            idtsedice=1
         else
            idtsedice=int(dt/dtsedice)+1
            dtsedice=dt/(idtsedice*1.)   
         endif

         do 200 ii=1,idtsedice
            do 200 ij=1,ijmax
               if(iminiceht(ij).le.imaxiceht(ij))then
                  imin=max0(1,iminiceht(ij)-ii)
                  do 201 ik=imin,imaxiceht(ij)

c frac gridbox volume moving to lower box

                     fracice=fluxcorr*vfallice(ij,ik)*dtsedice/delgrid
                     ph2osol=fracice*h2osol(ij,ik+1)
                     lh2osol=fracice*h2osol(ij,ik)
                     h2osol(ij,ik)=h2osol(ij,ik)+ph2osol-lh2osol
 201              continue
               endif
 200     continue
      endif
c      write(6,*)dtsedice,idtsedice
c     NAT particle size particulars:  the following variable 
c     conv converts the
c     number density of HNO3 molecules in solid phase into a surface
c     area density, assuming that the size distribution is lognormal,
c     as in Dye et al., J. Geophys. Res 97, 8015-8034, (1992).

      rg=1.e2*rnat              ! mode radius of NAT particles
      sigma=1.8             ! "stand dev"
      natgpmol=1.943e-22    ! grams/molecule of HNO3(3H2O)
      dens=1.6              ! grams/cm^3 of NAT

      conv=natgpmol/dens*3./rg*exp(-5.*alog(sigma)*alog(sigma)/2.)

c     Sediment solid phase HNO3:  HNO3 is divided into a gas-phase and a
c     solid-phase component.  Each of these fields contains the number
c     density (molecules/cm^3) of HNO3 in that particular phase.  The
c     solid-phase HNO3 is being advected as a passive tracer by the
c     model (though not in this subroutine).  Here, we account for
c     sedimentation assuming an average NAT fall velocity, vfallnat,
c     and a vertical grid point separation of 2 kilometers.

c     Actually, assuming that one fall
c     velocity (for all particles) is appropriate is a
c     a kluge - the larger particles will fall out faster.  Later
c     versions of this subroutine can deal with this.

c     this following loop may be replaceable with an equivalence
c     statement (and remove last loop in subroutine, too)

c     make subroutine understandable,   hno3sol(L$,Z$)

      if(pscnaton.eq.1)then
         do 1 ij=1,ijmax
               hno3sol(ij,IK40+1)=0.           ! add an extra level (=0) to accomodate code below, line 448
            do 1 ik=1,ikmax
               hno3sol(ij,ik)=cn(66,ij,ik) 
               hno3gas(ij,ik)=cn(10,ij,ik)
                                                              !  save initial total HNO3
               totna(ij,ik) = cn(10,ij,ik) + cn(66,ij,ik)
 1       continue


c max dt for sedimentation (depends on max fall velocity):

         dtsednat=delgrid/(fluxcorr*vfnatmax)

         if(dt/dtsednat.le.1.)then
            dtsednat=dt             ! no need to change timestep
            idtsednat=1
         else
            idtsednat=int(dt/dtsednat)+1
            dtsednat=dt/(idtsednat*1.)   
         endif
c         write(6,*)dtsednat,idtsednat
         
         do 2 ii=1,idtsednat
            do 2 ij=1,ijmax
               if(iminnatht(ij).le.imaxnatht(ij))then
                  imin=max0(1,iminnatht(ij)-ii)
                  do 21 ik=imin,imaxnatht(ij)

c fice and fnat are the fractions of the grid box that are saturated 
c with respect to ice and nat aerosols, respectively.  If ice is 
c present we assume that nat is also present, and the nat was used 
c as the condensation nuclei for the ice.  Thus, in the fraction 
c of the atmosphere saturated with respect to ice, nat will fall out 
c with the velocity of ice.  In the fraction not saturated, nat will 
c fall out with its regular sedimentation velocity.

                     if(fnat(ij,ik).gt.0.)then
                        if(fice(ij,ik)/fnat(ij,ik).lt.1.)then
                           frac=fluxcorr*dtsednat/delgrid*
     c                     (vfallnat(ij,ik)*(1.-fice(ij,ik)/fnat(ij,ik))
     c                     +vfallice(ij,ik)*fice(ij,ik)/fnat(ij,ik))
                        else
                           frac=fluxcorr*dtsednat/delgrid*
     c                          vfallice(ij,ik)
                        endif
                     else
                        frac=fluxcorr*dtsednat/delgrid*vfallnat(ij,ik)
                     endif

                     phno3sol=frac*hno3sol(ij,ik+1)
                     lhno3sol=frac*hno3sol(ij,ik)
                     hno3sol(ij,ik)=hno3sol(ij,ik)+phno3sol-lhno3sol
 21               continue
               endif
 2       continue
      endif


c     Main loop: most of the work is done here:

      do 3 ij=1,ijmax
         do 3 ik=1,ikmax

            hno3amb=hno3sol(ij,ik)+hno3gas(ij,ik) ! ambient hno3 conc.
            h2oamb=h2osol(ij,ik)+h2ogas(ij,ik) ! ambient h2o conc.

c     partial pressures of H2O, HNO3, in Torr:

            pphno3=760./1013.*press(ik)*hno3amb/m(ij,ik)

c     to calculate nucleation temp, tnuc:

            pphno3nuc=pphno3/natsatrat

c     the definition of ambient h2o is different because ice aerosols 
c     are forming,  and includes solid phase H2O.
c     For NAT formation, you only want to consider gas phase H2O.

            pph2o=760./1013.*press(ik)*h2ogas(ij,ik)/m(ij,ik)

c     coefficients for Hanson + Mauersberger eq. saturation temps:

            b=(-2.7836*alog10(pph2o)-alog10(pphno3)+38.9855)/
     c        (.009179-.00088*alog10(pph2o))


            bnuc=(-2.7836*alog10(pph2o)-alog10(pphno3nuc)+38.9855)/
     c        (.009179-.00088*alog10(pph2o))

            cc=-11397./(.009179-.00088*alog10(pph2o))

c     HNO3 saturation temperature:

            tsat=(-b+sqrt(b*b-4*cc))/2.
            tnuc=(-bnuc+sqrt(bnuc*bnuc-4*cc))/2.

c            if(j.eq.1.and.k.eq.11)write(6,*)tsat

c     get H2O saturation temperature, based on Marti and Mauersberger,
c     GRL 20 363-366 (1993).

            ah2o=-2663.5
            bh2o=12.537

c     partial pressure of H2O, in Pascal:

            parth2o=100.*press(ik)*h2oamb/m(ij,ik)
            parth2onuc=parth2o/satratice

c     H2O saturation temperature:

            tsath2o=ah2o/(alog10(parth2o)-bh2o)

c     H2O nucleation temperature:

            tnuch2o=ah2o/(alog10(parth2onuc)-bh2o)

c            if(j.eq.1.and.k.eq.11)write(6,*) tsat,tsath2o
c     figure out how much HNO3 comes out of gas phase:

            if(pscnaton.eq.1)then
               if(tsat.gt.205.)then

                  hno3min=0.
                  hno3max=0.
                  fnat(ij,ik)=0.
                  fnucnat(ij,ik)=0.

               else

c     get total amt of HNO3 removed at gridpoint:  I'm assuming here    
c     that the probability density distributions give the prob of 
c     finding a temp T between 170 and 260 K:  - changed to 120K-330K  (EF, Aug 2008)

                  hno3min=0.
                  hno3max=0.
                  fnat(ij,ik)=0.
                  fnucnat(ij,ik)=0.

                  do 4 itemp=120,330

                     rtemp=float(itemp) ! get real value
                     logpphno3=-1.*(2.7836+.00088*rtemp)*alog10(pph2o)
     c                   +38.9855-11397./rtemp+.009179*rtemp

c     calc. partial press of hno3 in eq at rtemp:, if logpphno3 too small, set to zero

                     if (logpphno3 .ge. -35.) hno3eq = 10.**logpphno3
                     if (logpphno3 .lt. -35.) hno3eq = 0.


c     convert to number density

                     hno3eq=m(ij,ik)*hno3eq/press(ik)*(1013./760.)

                     if(rtemp.lt.tnuc)then

                        hno3min=hno3min+(hno3amb-hno3eq)
     c                          *tprobnow(ij,ik,itemp-119)
                        fnucnat(ij,ik)=fnucnat(ij,ik)+
     c                              tprobnow(ij,ik,itemp-119)
                     endif

                     if(rtemp.lt.tsat)then

                        hno3max=hno3max+(hno3amb-hno3eq)
     c                          *tprobnow(ij,ik,itemp-119)
                        fnat(ij,ik)=fnat(ij,ik)+
     c                              tprobnow(ij,ik,itemp-119)

                     endif

 4                continue

               endif
            endif

c     Now do the same for H2O:

            if(psciceon.eq.1)then
               if(tsath2o.gt.200.)then
                  h2omin=0.
                  h2omax=0.
                  fice(ij,ik)=0.
                  fnucice(ij,ik)=0.

               else

c     get total amt of H2O removed at gridpoint:  I'm assuming here
c     that the probability density distributions give the prob of 
c     finding a temp T between 170 and 260 K:  - changed to 120K-330K  (EF, Aug 2008)

                  h2omin=0.
                  h2omax=0.
                  fice(ij,ik)=0.
                  fnucice(ij,ik)=0.


                  do 400 itemp=120,330

                     rtemp=float(itemp) ! get real value
                     logpph2o=ah2o/rtemp+bh2o

c     calc. partial press (in pascal) of h2o in eq at rtemp:

                     h2oeq=10.**logpph2o

c     convert to number density

                     h2oeq=m(ij,ik)*h2oeq/press(ik)*(1./100.)

                     if(rtemp.lt.tnuch2o)then

                        h2omin=h2omin+(h2oamb-h2oeq)
     c                          *tprobnow(ij,ik,itemp-119)
                        fnucice(ij,ik)=fnucice(ij,ik)+
     c                              tprobnow(ij,ik,itemp-119)

                     endif

                     if(rtemp.lt.tsath2o)then

                        h2omax=h2omax+(h2oamb-h2oeq)
     c                          *tprobnow(ij,ik,itemp-119)
                        fice(ij,ik)=fice(ij,ik)+
     c                              tprobnow(ij,ik,itemp-119)

                     endif

 400              continue

               endif
            endif


c     store new calculations of solid phase HNO3 and H2O:

            if(pscnaton.eq.1)then

               if(hno3sol(ij,ik).lt.hno3min)then
                  hno3sol(ij,ik)=hno3min
                  hno3gas(ij,ik)=hno3amb-hno3min
               endif

               if(hno3sol(ij,ik).gt.hno3max)then
                  hno3sol(ij,ik)=hno3max
                  hno3gas(ij,ik)=hno3amb-hno3max
               endif


c     convert solid-phase HNO3 to surface area density using conv. 
c     calculated above:

               nataer(ij,ik)=conv*hno3sol(ij,ik)
c               if(j.eq.1.and.k.eq.11)write(6,*)hno3gas(ij,ik),hno3rm,
c     c         hno3sol(ij,ik),dhno3gp


            endif

            if(psciceon.eq.1)then

               if(h2osol(ij,ik).lt.h2omin)then
                  h2osol(ij,ik)=h2omin
                  h2ogas(ij,ik)=h2oamb-h2omin
               endif

               if(h2osol(ij,ik).gt.h2omax)then
                  h2osol(ij,ik)=h2omax
                  h2ogas(ij,ik)=h2oamb-h2omax
               endif

               aerice(ij,ik)=convice*h2osol(ij,ik)

            endif


 3    continue
c                                 ! UPdate CN arrays here, ALSO NEED TO UPDATE THE C array here for TRANSPORT
      do 5 ij=1,ijmax
         do 5 ik=1,ikmax

            if(pscnaton.eq.1)then
               cn(66,ij,ik) = hno3sol(ij,ik) 
               if (cn(66,ij,ik).le.0.) cn(66,ij,ik) = 1.E-12
               c(66,ij,ik) = cn(66,ij,ik)

               cn(10,ij,ik) = hno3gas(ij,ik) 
c              (diffusion term figured in elsewhere)
               if (cn(10,ij,ik).le.0.) cn(10,ij,ik) = 1.E-12
               c(10,ij,ik) = cn(10,ij,ik)


cnoy   adjust NOy to reflect change in TOTAL HNO3 (GAS phase + SOLID) which occurs via sedimentation
C      so:   NOy = NOy + NEW total HNO3 - OLD total HNO3;    totna(L$,Z$)

               cn(31,ij,ik) = cn(31,ij,ik) + 
     >                 (cn(10,ij,ik) + cn(66,ij,ik)) - totna(ij,ik)

               if (cn(31,ij,ik) .le. 0.) then
                  write(6,*)'noy less than 0 at:' ,ij,ik, iday360
                  cn(31,ij,ik)=1.e-12
               endif

               c(31,ij,ik) = cn(31,ij,ik)
            endif


            if(psciceon.eq.1)then
               cn(67,ij,ik) = h2osol(ij,ik)
               if(cn(67,ij,ik).le.0.) cn(67,ij,ik) = 1.E-12
               c(67,ij,ik) = cn(67,ij,ik)

               cn(15,ij,ik) = h2ogas(ij,ik)
               c(15,ij,ik) = cn(15,ij,ik)
            endif

 5    continue

cc          if (idor .ge. indays-359) then
cc              if (mod(idor-15,IOUTP) .eq. 0.0) then
cc               write(88)nataer,aerice
cc              end if
cc          end if


 10   format(i10)
 11   format(6e11.3)

 1500 continue


c  Now compute how much H2O lost from ice fallout, put into MIXING RATIO/SECOND
C        load into ICE1(L$,Z$) for output (in COMMON)

         do 227 ik=1,Z$
         do 227 ij=1,L$
 227        ICE1(ij,ik) = (cn(15,ij,ik) - h2ost(ij,ik))/m(ij,ik)/DT

      RETURN
      END
