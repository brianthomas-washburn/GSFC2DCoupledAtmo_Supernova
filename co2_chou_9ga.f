	subroutine trco2v(nz,pa,ta,pu,lprint)
c	subroutine trco2v(nz,pa,ta,pu,lprint,flag,flag1)
c
c  co2 15 micron transmissions from parameterization of chou and kouvaris 
c  (1988) which allows for variable co2 mixing ratios.
c
c        include 'param.inc'
        logical flag,flag1
	real delp(48),delt(48)
	real pu(48), ta(49), pa(48)
	logical datread,lprint,co2var
	common /co2/ tauco2(48,48),taustd(49,49),fsn(19600)
	common /coefs/ absorb(31,28,3),datread
	common /co2mr/ co2vmr(45),co2lw(48),co2var
c	data npt /nv3$/
c
c  data is read in on first call only
c
	if (.not. datread) call co2coef

      if(lprint) then
        print *,'co2lw in trco2v'
        do i=1,npt
          print *,i,pa(i),co2lw(i)
        end do
      endif
c
c  initialize variables
c
      npt=nz+3
	tauco2(1,1) = 1.0
	do 10 i=2,npt
	    delp(i) = 794.1*co2lw(i)*(pu(i)-pu(i-1))
	    tauco2(i,i) = 1.0
c       if(flag .and. flag1) print *,'i,delp,co2lw,pu(i),pu(i-1) = ',
c     *  i,delp(i),co2lw(i),pu(i),pu(i-1)
10	continue
c       if(flag .and. flag1) print *,'at 1 in trco2v'
c
c  loop over pressure levels
c
	do 20 i=1,npt-1
	    wt = 0.0
	    pdw = 0.0
	    tdw = 0.0
	    do 30 j=i+1,npt
c       if(flag .and. flag1) print *,'i,j = ',i,j
c       if(flag .and. flag1) print *,'pa,ta,delp = ',pa(j),ta(j),delp(j)
c
c  compute effective temp, pressure and absorber amount
c
		wt = wt + delp(j)
		pdw = pdw + pa(j)*delp(j)
		tdw = tdw + ta(j)*delp(j)
		pefflog = alog10(pdw/wt)
		teff = tdw/wt
		wtlog = alog10(wt)
c       if(flag .and. flag1) print *,'pefflog,teff,wtlog = ',
c     *   pefflog,teff,wtlog 	
c
c  interpolate to find the absorption
c
		call interp2d(pefflog,wtlog,a,1)
c	 if(flag .and. flag1) print *,'a = ',a	
		call interp2d(pefflog,wtlog,alpha,2)
c	 if(flag .and. flag1) print *,'alpha = ',alpha	
		call interp2d(pefflog,wtlog,beta,3)
c	 if(flag .and. flag1) print *,'beta = ',beta	
		t250 = teff - 250.0
		a = 10.0**(-a)
c
c  find transmission at given temperature
c
		tauco2(i,j) = 1.0-a*(1.0+(alpha+beta*t250)*t250)
c      if(flag .and. flag1) print *,'a,tau = ', a,tauco2(i,j)
		tauco2(j,i) = tauco2(i,j)
 30	    continue
 20	continue
	return
	end

	subroutine interp2d(plog,wlog,extrap,index)
c
c  this subroutine uses the tables from chou and kouvaris (1988)
c  to interpolate the absorption and quadratic fit coeficients in
c  log10 pressure (mb) and log10 co2 amount (atm-cm).  if index=1
c  then a is returned, index=2 then alpha is returned and index=3
c  then beta is returned.
c
        logical datread
	common /coefs/ absorb(31,28,3),datread
	real lev1(31),lev2(28)
	data scale1 /0.2/, scale2 /0.3/, pmin /-3.0/, wmin /-5.1/
	data lev1 /-3.0,-2.8,-2.6,-2.4,-2.2,-2.0,-1.8,-1.6,-1.4,-1.2,
     *		   -1.0,-0.8,-0.6,-0.4,-0.2, 0.0, 0.2, 0.4, 0.6, 0.8,
     *		    1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8,
     *		    3.0/
	data lev2 /-5.1,-4.8,-4.5,-4.2,-3.9,-3.6,-3.3,-3.0,-2.7,-2.4,
     *		   -2.1,-1.8,-1.5,-1.2,-0.9,-0.6,-0.3, 0.0, 0.3, 0.6,
     *		    0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0/
c
	i = int((plog - pmin)/scale1) + 1
	j = int((wlog - wmin)/scale2) + 1
	delp = (plog - lev1(i))/scale1
	delw = (wlog - lev2(j))/scale2
	h1 = delp*(absorb(i+1,j,index)-absorb(i,j,index)) + 
     *   absorb(i,j,index)
	h2 = delp*(absorb(i+1,j+1,index)-absorb(i,j+1,index)) +
     *   absorb(i,j+1,index)
	extrap = delw*(h2 - h1) + h1
	return
	end


	subroutine co2coef
c
c  USED IN THE CHOU CO2 PARAMETERIZATION
c  this subroutine reads in the coeficients for interpolating
c  the absorption coef.  see chou and kouvaris (1988).
c
	logical datread,co2var
c        include 'param.inc'
	data datread /.false./
        common /co2mr/ co2vmr(45),co2lw(48),co2var
	common /coefs/ coef(31,28,3),datread

c	open(33,file='$RADLIB/co2coef.dat',status='old',
c     *						form='formatted')

	do 10 k=1,3
	    do 20 j=1,28
		read(113,*) (coef(i,j,k),i=1,31)
 20	    continue
 10 	continue
	do 30 i=1,31
	    do 40 j=1,28
		coef(i,j,1) = coef(i,j,1)*1.0e-3 
		coef(i,j,2) = coef(i,j,2)*1.0e-4 
		coef(i,j,3) = coef(i,j,3)*1.0e-6
 40	    continue
 30	continue 
c
c  set datread so that this routine is called only once
c
	datread = .true.

	close (113)

	return
	end
