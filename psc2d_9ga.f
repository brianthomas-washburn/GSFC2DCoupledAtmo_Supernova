      subroutine psc(jlat,psclat,lprint, YIN, ZIN)
      
c solid NAT and ICE optical thicknesses
c this has been modified in order to estimate the heating perturbation due to
c tropical ice particles predicted by the coupled 2d model
c
      include 'PARAM.INC'
      include 'COMMONR1.INC'
      
      common /solid1/ ssh2o(N$,ms$),sshno3(N$,ms$)
         
      common /pscld/ ipsc(45),pscabs(45),pscext(45),psch(45),
     * psc96(45),psc15(45)
     

      !integer NS$, MS$

      logical lprint 
      real los, yin(ns$),zin(ms$)
      real m,natsat
      dimension pll(45)
      dimension p2d(ms$),pl2d(ms$),h2opsc(45),ptemp(ms$),htemp(ms$)
      dimension hno3psc(45),h1temp(ms$)
c
c can have heating due to ice or nat or both
      logical nat,ice
      data nat/.true./, ice/.true./
c
      data los/2.69e19/, pi/3.1415927/
      data denice/0.917/
c following are grams/molecule
      data gpmolice/2.991e-23/, gpmolnat/1.943e-22/
      
      N2d=ms$
      NZ=45
c
      if(lprint) write(6,2)
    2 format('  printout from psc')
     
      do 15 n=1,NZ
      pll(n)=alog(pl(n))
   15 continue
c
c define p2d array
      do 20 n=1,N2d
           p2d(n) = 1013.*EXP(-zin(n)/7.)
ccelf      p2d(n) = 1013.*exp(-.2844*((n-1)+0.5))
   20 continue
      
      if(.not. lprint) go to 860
c      if(jlat .ne. 3 .and. jlat .ne. 20) go to 860
      write(6,*) 'jlat = ',jlat 
      write(6,*) 'solid h2o (mol/cc) before vert interp'
      do 801 i=1,15
c      write(6,*) i,p2d(i),ssh2o(jlat,i)
      write(6,*) i,p2d(i),sshno3(jlat,i)
  801 continue
  860 continue
  
c reverse psc arrays to get them top down
c      write(6,*) 'solid h2o, reversed in vertical'
      do 25 n=1,N2d
      nn=N2d-n+1
      ptemp(n)=p2d(nn)
      htemp(n)=ssh2o(jlat,nn)
      h1temp(n)=sshno3(jlat,nn)
      if(htemp(n).lt.1.) htemp(n)=0.0
      if(h1temp(n).lt.1.) h1temp(n)=0.0
c      write(6,*) 'n,ptemp,htemp: ',n,ptemp(n),htemp(n)
c      write(6,*) 'n,ptemp,h1temp: ',n,ptemp(n),h1temp(n)
   25 continue
   
      do 26 n=1,N2d
   26 pl2d(n)=alog(ptemp(n))
           
c  interpolate htemp profile to model vertical grid
      call inter(h2opsc,pll,NZ,htemp,pl2d,N2d)
      call inter(hno3psc,pll,NZ,h1temp,pl2d,N2d)
      do 27 n=1,NZ
      if(h2opsc(n).lt.0.0) h2opsc(n)=0.0
      if(hno3psc(n).lt.0.0) hno3psc(n)=0.0
   27 continue
      
c      if(lprint) write(6,800)
c  800 format('  solid h2o after interpolation')
c      do 830 n=1,NZ
c      if(lprint) write(6,805) n,pl(n),hi(n),h2opsc(n) 
c  805 format(i4,2f12.4,e12.4)
c  830 continue
c  870 continue

      if(lprint) write(6,800)
  800 format('  solid hno3 after interpolation')
      do 830 n=1,NZ
      if(lprint) write(6,805) n,pl(n),hi(n),hno3psc(n) 
  805 format(i4,2f12.4,e12.4)
  830 continue
  870 continue

c
c  main loop
c
      do 500 n=10,NZ
c
      ipsc(n)=0
      pscabs(n)=0.
      pscext(n)=0.
      psch(n)=0.
      psc96(n)=0.
      psc15(n)=0.
c
c don't allow these type of clouds for pl greater than 250 mb
c or for altitudes greater than 25 km
c
      if(pl(n) .gt. 250.) go to 500
      if(hi(n) .gt. 25.) go to 500
       
      if(.not. ice) go to 300
c

      if(h2opsc(n).le.0.0) go to 300

c
c forming ice psc
c
cjer
c  ice polar stratospheric cloud cross sections
      if(abs(psclat) .gt. 30.) then

c         ice  (rmod = 10.0)
c
c          rmod=10.0*1.e-4
          r3mean=2095.*1.e-12
c         solar
          solabs=1.3e-7
          solext=9.2e-6
          onemic=9.2e-6
c         ir
          sigh = 4.1e-6
          sig96= 3.9e-6
          sig15= 5.5e-6

      else
c
c         ice  (rmod = 2.0) a log-normal size dist.
c
c          rmod=2.0*1.e-4
c          r3mean=16.76*1.e-12
c         solar
c          solabs=3.0e-9
c          solext=4.2e-7
c          onemic=4.1e-7
c         ir
c          sigh = 6.6e-8
c          sig96= 4.9e-8
c          sig15= 1.3e-7
c
c          ice  (rmod = 6.0)
c
          rmod=6.0e-4
          r3mean=452.6*1.e-12
c          solar
          solabs=3.7e-8
          solext=3.4e-6
          onemic=3.4e-6
c          ir
           sigh = 1.2e-6
           sig96= 1.1e-6
           sig15= 1.9e-6
c
c          ice  (rmod = 8.0)  a log-normal size dist. (see my psc paper)
c
c          rmod=8.0e-4
c          r3mean=1072.9*1.e-12
c          solar
c          solabs=7.4e-8
c          solext=6.0e-6
c          onemic=6.0e-6
c          ir
c          sigh = 2.3e-6
c          sig96= 2.2e-6
c          sig15= 3.5e-6
c
c         ice  (rmod = 10.0)
c
c          rmod=10.0*1.e-4
c          r3mean=2095.*1.e-12
c         solar
c          solabs=1.3e-7
c          solext=9.2e-6
c          onemic=9.2e-6
c         ir
c          sigh = 4.1e-6
c          sig96= 3.9e-6
c          sig15= 5.5e-6
c
c          ice  (rmod = 15.0)
c
c          rmod=15.0*1.e-4
c          solar
c          solabs=2.1e-7
c          solext=1.5e-5
c          ir
c          sigh = 6.8e-6
c          sig96= 6.8e-6
c          sig15= 9.2e-6

      end if

c
      if(lprint) write(6,*)
      if(lprint) write(6,*) ' n, pressure = ',n,pl(n)
      if(lprint) write(6,199) rmod
  199 format(/,'  ice psc,  rmod (cm) = ',e12.2,/)
c
      ipsc(n)=1
      cloud(n)=1
      
c  mass of ice in cc of air
      amt = gpmolice * h2opsc(n)
c      if(lprint) write(6,203) h2opsc(n),amt
  203 format(' molec., grams of ice in cc of air = ',2e12.3)
  
c  volume of ice in cc of air
      volt = amt/denice
c  average particle volume (cm)
c      volp=4./3.*pi*rmod**3
      volp=4./3.*pi*r3mean
c      if(lprint) write(6,205) volt,volp
  205 format(' total volume of ice = ',e12.3,/,
     * ' volume of one particle = ',e12.3)
     
c  number of ice particles in cc of air
      den=volt/volp
c      if(lprint) write(6,207) den
  207 format(' number of ice particles in cc of air = ',e12.3)
  
c  column density in grid box
      htkm=hie(n) - hie(n+1)
      colden=den*( htkm )*1.e5
c      if(lprint) write(6,209) colden
  209 format(' column density of ice in grid box (cm-2) =',e12.3)
c
      go to 400
c
  300 continue
c
c  forming nat psc's
c
      if(.not. nat) go to 400
c
c      if(lprint) write(6,13)
   13 format(/,' checking for nat')
c
c  nat polar stratospheric cloud cross sections
c
c taesler et al give nat density of 1.62 
c
c  hno3 (70%)  (rmod = 0.8)   
c      rmod=0.8*1.e-4
c      dennat=1.62
c      solabs=1.0e-8
c      solext=3.8e-8
c      sigh = 1.2e-8
c      sig96= 1.1e-8
c      sig15= 8.0e-9
c  hno3 (40%)  (rmod = 0.8)
      rmod=0.8e-4
      dennat=1.62
      solabs=1.6e-8
      solext=3.7e-8
      onemic=8.5e-8
      sigh = 1.3e-8
      sig96= 1.1e-8
      sig15= 9.9e-9
c
c (begin aase comments)
c  hno3 (40%)  (rmod = 0.5)
c      rmod=0.5e-4
c      dennat=1.62
c      solabs=4.3e-9
c      solext=7.9e-9
c      onemic=3.1e-8
c      sigh = 2.9e-9
c      sig96= 2.6e-9
c      sig15= 2.3e-9
c
c      rmod=1.2e-4
c      onemic=1.6e-7
c (end aase comments)

      if(lprint) write(6,10) rmod
   10 format(/,'  nat psc,   rmod (cm) = ',e12.2,/)

      if(lprint) write(6,301) n,pl(n),tl(n)
  301 format(' n, pl,tl:  ',i4,2f10.2)
c
c  mass of nat in cc of air
      amt=gpmolnat * hno3psc(n)
      if(amt.lt.0.0) go to 500
      ipsc(n)=1
      cloud(n)=1

c  volume of nat in cc of air
      volt=amt/dennat
c  volume of one particle
      volp=4./3.*pi*rmod**3
      if(lprint) write(6,305) volt,volp
  305 format(' total volume of nat in cc of air = ',e12.3,/,
     * ' volume of one particle = ',e12.3)

c  number of nat particles in cc of air
      den=volt/volp
      if(lprint) write(6,307) den
  307 format(' number of nat particles in cc of air = ',e12.3)

c  column density in grid box
      htkm=hie(n) - hie(n+1)
      colden=den*( htkm )*1.e5
      if(lprint) write(6,309) colden
  309 format(' column density of nat in grid box=',e12.3)

  400 continue
c
c  optical depths
c
c  for near ir solar
      pscabs(n)=solabs*colden
      pscext(n)=solext*colden
c  for ir
      opsch=sigh*colden
      opsc96=sig96*colden
      opsc15=sig15*colden

c
      if(.not.lprint) go to 121
      if(ipsc(n) .eq. 0) go to 500
      write(6,110)
  110 format(/,'  ir optical depths',/)
      write(6,116)
  116 format('   opsch    opsc96    opsc15')
      write(6,117) opsch,opsc96,opsc15
  117 format(3f10.3)
  121 continue

      psch(n)=1.0 - exp(-1.66*opsch)
      psc96(n)=1.0 - exp(-1.66*opsc96)
      psc15(n)=1.0 - exp(-1.66*opsc15)

      if(.not.lprint) go to 131
      write(6,120) 
  120 format(/,'  ir emissivities',/)
      write(6,118)
  118 format('    psch   psc96   psc15')
      write(6,117) psch(n),psc96(n),psc15(n)
      write(6,125)
  125 format(/,'  solar near infrared optical depths',/)
      write(6,130)
  130 format('  pscext   pscabs')
      write(6,117) pscext(n),pscabs(n)
  131 continue
c
c  set psc=0 if optical depths all less than 1.e-4
      eps=1.e-4
      if( (pscext(n).lt. eps) .and. (psch(n).lt. eps) .and.
     *  (psc96(n).lt. eps) .and. (psc15(n).lt. eps) ) then
       
         ipsc(n)=0
         cloud(n)=0.
         
      endif
c
  500 continue
c
c  do fcloud, fclear again for pscs
c
c  cloud top, total and fractional cloudiness parameters
c  this code is taken from solar1 of the glas gcm
c    cloud fraction is the maximum of the fractions of all cloud types
      nz1=NZ+1
      fcloud = 0.0
      fclear = 1.0
      ntopt = nz1
      ntopf = nz1
      icld=1
               do 540 n=icld,NZ
      xx = cloud(n)
               if (xx.lt.0.01)  go to 540
               if (xx.gt.0.99)  go to 530
      fc = amax1(xx,fcloud)
      fcloud = fc
      fclear = 1.0 - fcloud
               if (ntopf.lt.NZ)  go to 540
      ntopf = n
               go to 540
  530          continue
               if (ntopt.lt.NZ)  go to 540
      ntopt = n
      fclear = 0.0
  540          continue
               if (fclear.gt.0.99)  go to 600
      if (fcloud.lt.0.01)  fcloud=1.00
  600 continue
      ntop=min0(ntopt,ntopf)
c
c
      return
      end
