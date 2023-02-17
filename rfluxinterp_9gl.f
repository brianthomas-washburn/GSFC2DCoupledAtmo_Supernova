       SUBROUTINE RFLUXINTERP

C   Received from Randy Kawa 1/24/95 - adapted from photcodes7.f 

        include "com2d.h"
        include "comphot.h"

cccccc      	real tempuse(nlat,nz,nx,nlam),tuse(nlat)
      	real tuse(nlat)
        integer locsza(L$)
 
c       interpolate radiative flux function values to model conds. 

c       For each input solar zenith angle, find the first element of
c       tabled sza_tab values that is greater than it.  Use this
c       table element and previous table element to determined 
c       interpolated value. 
        do ij=1,L$
           do is=1,nsza
              if (sza_tab(is).gt.sza2d(ij)) go to 333 
           end do 
 333       locsza(ij) = is 
        end do
            
        do ij=1,nlat
           if (locsza(ij).eq.(nsza+1)) then 
c             Point is in darkness.  Set s's to zero. 
              do ik=1,nz 
                 o2jint(ij,ik) = 0.
              do il=1,nlam
                 sflux(ij,ik,il) = 0. 
              end do 
              end do 
           else  
              ijj = locsza(ij)                       
              tuse(ij)=(sza2d(ij)-sza_tab(ijj-1))/
     1         (sza_tab(ijj)-sza_tab(ijj-1)) 
c             For each input overhead o3 column find the first element
c              of tabled o3_tab values that is > than it.  Use this
c              table element and previous table element to determine
c              interpolated value. 
C
C   WCONV277 - lookup table arrays are interpolated to current model pressure grid in PHOTIN
C              stab(nsza,no3,nz,nlam); o2jdat(nsza,no3,nz); o3_tab(no3,nz),  nz=Z$
C
C              so just define ikt=ik here; bilinear interpolation is in SZA and ozone column
C

              do ik=1,nz

                 ikt = ik
C
c  this finds the pr_tab pressure level that's closest to 2D model pressure level in question
C
cc277                 if (press(ik) .gt. pr_tab(1)) then
cc277                    ikt=1
cc277                    goto 8888
cc277                 endif
cc277
cc277                 if (press(ik) .lt. pr_tab(100)) then 
cc277                    ikt=100
cc277                    goto 8888
cc277                 endif

cc277                 do 777 iklu=1,99
cc277                   if (press(ik) .lt. pr_tab(iklu)  .and.  
cc277     >                 press(ik) .gt. pr_tab(iklu+1)) then

cc277                     if (abs(alog(press(ik)/pr_tab(iklu))) .le. 
cc277     >                   abs(alog(press(ik)/pr_tab(iklu+1)))) ikt=iklu
cc277                     if (abs(alog(press(ik)/pr_tab(iklu))) .gt. 
cc277     >                   abs(alog(press(ik)/pr_tab(iklu+1)))) ikt=iklu+1
cc277                     goto 8888
cc277                   endif
cc277 777             continue


 8888            do is=1,no3 
                    if (o3_tab(is,ikt) .gt. ovho3(ij,ik)) go to 334 
                 end do 

 334             ikk = is 
                 ikkm = ikk-1 
                 if ((ikk.gt.1).and.(ikk.le.no3)) then 
                    u=(ovho3(ij,ik)-o3_tab(ikkm,ikt))/
     1               (o3_tab(ikk,ikt)-o3_tab(ikkm,ikt))          

c                   do bilinear interpolation at ik for each wavelength
c                   from numerical recipes, p.96

	            do il=1,nlam       
	               sflux(ij,ik,il)=(1.-tuse(ij))*(1.-u)*
     1                  stab(ijj-1,ikkm,ikt,il)+tuse(ij)*(1.-u)*
     1		        stab(ijj,ikkm,ikt,il)+tuse(ij)*u*
     1                  stab(ijj,ikk,ikt,il)+
     1	         	(1.-tuse(ij))*u*stab(ijj-1,ikk,ikt,il)
                    end do

	    o2jint(ij,ik)=(1.-tuse(ij))*(1.-u)*o2jdat(ijj-1,ikkm,ikt)
     1			+tuse(ij)*(1.-u)*o2jdat(ijj,ikkm,ikt)
     2			+tuse(ij)*u*o2jdat(ijj,ikk,ikt)
     3	         	+(1.-tuse(ij))*u*o2jdat(ijj-1,ikk,ikt)
                                
                 else if (ikk.eq.1) then 
c                   write (6,33) ij,ik,ovho3(ij,ik),o3_tab(1,ikt)
 33                 format(' Ovhd o3 col(',i3,i3,') of ',e10.3,
     1               ' too thin! Lowest tabled val=',e10.3)  
      
                    do il=1,nlam
	             sflux(ij,ik,il)=(1.-tuse(ij))*stab(ijj-1,1,ikt,il)+
     1                  tuse(ij)*stab(ijj,1,ikt,il)
	            end do        
	        o2jint(ij,ik)=(1.-tuse(ij))*o2jdat(ijj-1,1,ikt)+
     1                  tuse(ij)*o2jdat(ijj,1,ikt)

                 else 
c                    write (6,34) ij,ik,ovho3(ij,ik),o3_tab(no3,ikt)
  34                 format(' Ovhd o3 col(',i3,i3,') of ',e10.3,
     1               ' too thick! Highest tabled val=',e10.3) 
                    do il=1,nlam
	          sflux(ij,ik,il)=(1.-tuse(ij))*stab(ijj-1,no3,ikt,il)+
     1                  tuse(ij)*stab(ijj,no3,ikt,il)
	            end do 
	          o2jint(ij,ik)=(1.-tuse(ij))*o2jdat(ijj-1,no3,ikt)+
     1                  tuse(ij)*o2jdat(ijj,no3,ikt)

                 end if  
              end do 
	    end if
           end do 


        return
        end 
