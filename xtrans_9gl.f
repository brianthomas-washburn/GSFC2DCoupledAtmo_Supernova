C
C
	SUBROUTINE XTRANSPORT

C
C  routine to do all transport (ADVECTION + DIFFUSION)                        - EF 2/28/03
C     code is just been taken out of MAIN and put into this separate routine
C
C    ALL ARRAYS HERE ARE IN COMMON
C


        include "com2d.h"


        REAL*8 mass1(2), mass2(2), mass6(2)
        REAL*8 massyy(2), massyz(2), masszz(2)
        REAL*8 sumyy(2), sumyz(2), sumzz(2)


	SAVE


cef
cef        if (iday360.ge.158) write(27) c, m, w, v


c
C    ************************    ADVECTION    -   PPM  TRANSPORT ROUTINE  ************************
c
C
c  TPCORE2D uses mixing ratio as input xmr(L$,Z$X,T$), and converts to density in routine
c  Note, cn58 and c58 are just transported species, cn58, c58 (T$,L$,Z$X), XMR(L$,Z$X,T$)  in COMMON

	DO 5700 IT1=1,itrans
	DO 5700 IK=1,Z$
	DO 5700 IJ=1,L$
     	    xmr(IJ,IK,IT1) = c(inttra(IT1),IJ,IK)/M(IJ,IK)
5700	CONTINUE

c
c   Now extend transported species to 58 levels, use new c58 array 
c
	do 5701 IT1=1,itrans
	do 5701 IK=Z$+1,Z$X
	do 5701 IJ=1,L$
	    xmr(IJ,IK,IT1) = c58(IT1,IJ,IK)/M(IJ,IK)
5701	CONTINUE
c
c 
C  determine the number of time steps per day for advection
C    compute DT for a maximum w of +- 5 cm/sec, and for a maximum v of +-2000 cm/sec
C    this should handle almost all if not all extremes. To ensure that the Courant 
C    condition is not violated for the given latitude and altitude grid, use the smallest DT of the two
C    the w limit of +- 5 cm/sec is set in VWMASS,  DELYC and DELZC (constituent grid) are in meters, REAL*8 
c
      idtw = IDNINT(87600./(delzc*100.)*5.) 
      idtv = IDNINT(87600./(delyc*100.)*2000.) 
      iadt = MAX(idtv, idtw) + 1
      if (L$ .eq. 18  .and.  Z$ .eq. 46) iadt = 2

c  need to decrease the time step here at the end of 2002 so that NOy doesn't blow up at ~84 km, 65S-75S
      if (L$ .eq. 18  .and.  Z$ .eq. 46  .and.  
     >    iyr .eq. 67   .and.  iday360 .ge. 340) iadt = 3


        DT8 = DT8/DBLE(IADT)

        do 8846 ie=1,IADT
             CALL TPCORE2D
 8846	continue

	DT8 = DT8*DBLE(IADT)

c convert back to density - load into the C array

	do 5900 IK=1,z$
	do 5900 IJ=1,l$
	do 5900 IT1=1,itrans
           c(inttra(IT1),IJ,IK) = xmr(IJ,IK,IT1)*m(IJ,IK)
cccccccc           cn(inttra(IT1),IJ,IK) = xmr(IJ,IK,IT1)*m(IJ,IK)
5900	CONTINUE

	do 5901 IK=1,Z$X
	do 5901 IJ=1,L$
	do 5901 IT1=1,itrans
           c58(IT1,IJ,IK) = xmr(IJ,IK,IT1)*m(IJ,IK)
cccccccc           cn58(IT1,IJ,IK) = xmr(IJ,IK,IT1)*m(IJ,IK)
5901	CONTINUE


cef
cef        if (iday360.ge.158) write(27) c, m, xmr, w, v



ccm
ccm       do 116 it1=1,itrans
ccm         mass1(it1) = 0.d0
ccm         mass2(it1) = 0.d0
ccm         mass6(it1) = 0.d0
ccm
ccm         massyy(it1) = 0.d0
ccm         massyz(it1) = 0.d0
ccm         masszz(it1) = 0.d0
ccm
ccm         sumyy(it1) = 0.d0
ccm         sumyz(it1) = 0.d0
ccm         sumzz(it1) = 0.d0
ccm 116   CONTINUE
c             
ccm                        
ccm       do 117 it1=1,itrans
ccm       do 117 ik=1,Z$X
ccm       do 117 ij=1,L$     
ccm         mass1(it1) = mass1(it1) + c58(it1,ij,ik)*cosc(ij)*deltaz(ij,ik)
ccm 117   CONTINUE
ccm

cc       write (27) L$, Z$X, T$, IYDT 
cc       write (27) c58

c
c     *********************     DIFFUSION   ***********************************
C
C    Diffusion is now done together with the time step needed for KZZ
C
c
C  NEWDIFY (array DIFNY), use a 1/IYDT of a day time step, to accomodate new FINER GRID RESOLUTION, 
C   calculate diffusion terms for each transported constituent, NOW DONE FOR Z$X LEVELS *******
C                         ST(T$,L$,Z$X) in COMMON
C
C   IYDT = number of time steps per day, based on the number of latitudes (L$) 
C
C   check mass conservation, to ensure that the total number of molecules within the model 
C       domain are the same before and after diffusion, COSC(L$), TCOS are REAL*8 in COMMON
C        also, the total area and conversion of deltaz to cm are not included here, since they are just 
C        constants and cancel when taking the ratio of mass1/mass2
C
c
cc       IYDT = INT((L$/18)**2) + 1
cc       if (L$ .eq. 18) IYDT = 1
cc   
cc       DT8 = DT8/DBLE(IYDT)
cc
cc 
cc        DO 1554 ied=1,IYDT
cc 
cc 	  DO 6102 IK=1,Z$
cc 	  DO 6102 IJ=1,L$
cc 	  DO 6102 IT1=1,ITRANS
cc 6102        ST(IT1,IJ,IK) = c(inttra(it1),IJ,IK)/m(IJ,IK)
C                                                              Extend diffusion array up to extended levels
cc          DO 6103 IK=Z$+1,Z$X
cc 	  DO 6103 IJ=1,L$
cc 	  DO 6103 IT1=1,ITRANS
cc 6103        ST(IT1,IJ,IK) = c58(it1,IJ,IK)/m(IJ,IK)
cc    
cc   
cc                CALL NEWDIFY
cc  
cc  
cc        do 7702 IT1=1,itrans
cc            do 7715 IK=1,z$
cc            do 7715 IJ=1,l$
cc 7715  c(inttra(IT1),IJ,IK) = c(inttra(IT1),IJ,IK)- DIFNY(IT1,IJ,IK)*DT8
cc     
cc      	   do 7754 IK=1,Z$X
cc     	   do 7754 IJ=1,L$
cc 7754   c58(IT1,IJ,IK) = c58(IT1,IJ,IK) - DIFNY(IT1,IJ,IK)*DT8
cc 7702   CONTINUE
cc 
cc 1554   CONTINUE
cc 
cc     	DT8 = DT8*DBLE(IYDT)
cc
cc
c
c ************************     KYZ DIFFUSION   ***********************************
C
c
C  NEWDIFYZ (array DIFNYZ), use a 1/IYDT of a day time step, to accomodate new FINER GRID RESOLUTION, 
C   calculate diffusion terms for each transported constituent, NOW DONE FOR Z$X LEVELS *******
C                         ST(T$,L$,Z$X) in COMMON
C
C   IYDT = number of time steps per day, based on the number of latitudes (L$) 
C
C   check mass conservation, to ensure that the total number of molecules within the model 
C       domain are the same before and after diffusion, COSC(L$), TCOS are REAL*8 in COMMON
C        also, the total area and conversion of deltaz to cm are not included here, since they are just 
C        constants and cancel when taking the ratio of mass1/mass2
C
c
cc         IYZDT = 1      !  IYDT*10
cc         DT8 = DT8/DBLE(IYZDT)
cc 
cc         DO 1557 ied=1,IYZDT
cc 
cc 	  DO 7102 IK=1,Z$
cc 	  DO 7102 IJ=1,L$
cc 	  DO 7102 IT1=1,ITRANS
cc 7102        ST(IT1,IJ,IK) = c(inttra(it1),IJ,IK)/m(IJ,IK)
C                                                              Extend diffusion array up to extended levels
cc          DO 7103 IK=Z$+1,Z$X
cc 	  DO 7103 IJ=1,L$
cc 	  DO 7103 IT1=1,ITRANS
cc 7103        ST(IT1,IJ,IK) = c58(it1,IJ,IK)/m(IJ,IK)
cc    
cc   
cc           CALL NEWDIFYZ
cc  
cc  
cc        do 8803 IT1=1,itrans
cc            do 7725 IK=1,Z$
cc            do 7725 IJ=1,L$
cc7725  c(inttra(IT1),IJ,IK) = c(inttra(IT1),IJ,IK)- DIFNYZ(IT1,IJ,IK)*DT8
cc
cc     
cc      	   do 7774 IK=1,Z$X
cc     	   do 7774 IJ=1,L$
cc 7774   c58(IT1,IJ,IK) = c58(IT1,IJ,IK) - DIFNYZ(IT1,IJ,IK)*DT8
cc8803   CONTINUE
cc 
cc 1557   CONTINUE
cc 
cc     	DT8 = DT8*DBLE(IYZDT)
cc
cc

c *********************   END   KYZ DIFFUSION   ***********************************


c
c     *********************     KZZ DIFFUSION   ***********************************
C
C
c  NEWDIFZ (array DIFNZ), with a 1/IZDT of a day time step, to accomodate up to 100 m2/sec Kzz in mesosphere
C    IZDT = number of time steps per day; 
C    calculate diffusion terms for each transported constituent, NOW DONE FOR Z$X LEVELS *******
C    ST(T$,L$,Z$X) in COMMON,    key on DELZC in meters, REAL*8
C  

ccccc        IZDT = 48 and 54 are no good,  60 and 72 work w/ 0.5 km for 8-30 km, although 72 is better
cccc        IZDT = 72                   16 works for 1 km for 0-60 km, w/ 2 km above 60 km where Kzz is large
cccc                                    12 works for 2 km
cccc        IZDT = 30
cccc        IZDT = 16
cccc        IZDT = 24

        if (delzc .ge. 1900.) IZDT = 12
        if (delzc .ge. 900.  .and.  delzc .lt. 1900.) IZDT = 24        ! 16 - better to use 24 here......
        if (delzc .ge. 400.  .and.  delzc .lt. 900.) IZDT = 60         ! 72
        if (delzc .ge. 200.  .and.  delzc .lt. 400.) IZDT = 120 

        if (delzc .ge. 95.  .and.  delzc .lt. 200.) IZDT = 240         ! 6 minute time step


cccccccccc        IZDT = 12


        if (IYR .eq. 35  .or.  IYR .eq. 55) then 
           if (iday360 .le. 2)  print *, 'in XTRANS:   IZDT = ', izdt
        endif


        DT8 = DT8/DBLE(IZDT)

      DO 1551 ied=1,IZDT

	DO 6150 IK=1,Z$
	DO 6150 IJ=1,L$
	DO 6150 IT1=1,ITRANS
 6150      ST(IT1,IJ,IK) = c(inttra(it1),IJ,IK)/m(IJ,IK)
C                                                              Extend diffusion array up to Z$X levels
	DO 6151 IK=Z$+1,Z$X
	DO 6151 IJ=1,L$
	DO 6151 IT1=1,ITRANS
 6151      ST(IT1,IJ,IK) = c58(it1,IJ,IK)/m(IJ,IK)


           CALL NEWDIFY
           CALL NEWDIFYZ
           CALL NEWDIFZ

 
         do 7800 IT1=1,itrans

             do 7810 ik=1,Z$
             do 7810 ij=1,L$
               c(inttra(it1),ij,ik) = c(inttra(it1),ij,ik) - 
     >     (DIFNY(it1,ij,ik) + DIFNYZ(it1,ij,ik) + DIFNZ(it1,ij,ik))*DT8
c                                                                             ! reset neg values to small pos
       IF (c(inttra(it1),ij,ik) .LT. 1.D-12) c(inttra(it1),ij,ik)=1.D-12
 7810	  CONTINUE

             do 7850 ik=1,Z$X
             do 7850 ij=1,L$
               c58(it1,ij,ik) = c58(it1,ij,ik) - 
     >     (DIFNY(it1,ij,ik) + DIFNYZ(it1,ij,ik) + DIFNZ(it1,ij,ik))*DT8
c                                                                             ! reset neg values to small pos
              IF (c58(it1,ij,ik) .LT. 1.D-12) c58(it1,ij,ik) = 1.D-12
 7850	   CONTINUE
 
 7800  CONTINUE

 1551	continue

    	DT8 = DT8*DBLE(IZDT)


cef
cef        if (iday360.ge.158) write(27) c, m, st, DIFNY, DIFNYZ, DIFNz

C
c *********************   END  DIFFUSION   ***********************************


cc       write (27) ekyy, ekyz, ekzz
cc       write (27) difny, difnyz, difnz, cosc, deltaz, m, c58 



c  check mass conservation
c      REAL*8 difny(T$,L$,Z$X), difnz(T$,L$,Z$X), difnyz(T$,L$,Z$X) , massyy(2), massyz(2), masszz(2) 
c
ccm
ccm       do 217 it1=1,itrans
ccm       do 217 ik=1,Z$X
ccm       do 217 ij=1,L$     
ccm         mass2(it1) = mass2(it1) + c58(it1,ij,ik)*cosc(ij)*deltaz(ij,ik)
ccm
ccm      massyy(it1)= massyy(it1) + difny(it1,ij,ik)*cosc(ij)*deltaz(ij,ik)
ccm       sumyy(it1) = sumyy(it1) + difny(it1,ij,ik)
ccm
ccm      massyz(it1)=massyz(it1) + difnyz(it1,ij,ik)*cosc(ij)*deltaz(ij,ik)
ccm       sumyz(it1) = sumyz(it1) + difnyz(it1,ij,ik)
ccm
ccm      masszz(it1)= masszz(it1) + difnz(it1,ij,ik)*cosc(ij)*deltaz(ij,ik)
ccm       sumzz(it1) = sumzz(it1) + difnz(it1,ij,ik)
ccm 217   CONTINUE
ccm
ccm
ccm       do 216 it1=1,itrans
ccm  	  mass6(it1) = (mass2(it1) - mass1(it1))/mass1(it1)
ccm
ccm          massyy(it1) = massyy(it1)/sumyy(it1)
ccm          massyz(it1) = massyz(it1)/sumyz(it1)
ccm          masszz(it1) = masszz(it1)/sumzz(it1)
ccm 216   CONTINUE
ccm
C  write out
C
cc        if (iday360 .eq. 300) then
cc          if (it .eq. 1) write(26,504)
cc          if (it .eq. 1) write(26,504)
ccccc          write(26,504)

cc        write(26,505) iday360, mass1(1), mass6(1), mass1(2), mass6(2), 
cc     >  massyy(1), massyz(1), masszz(1), massyy(2), massyz(2), masszz(2)

cc         sumyy(1), sumyz(1), sumzz(1), sumzz(2), 


cc     write(26,505) iday360, (difnyz(1,9,ik), ik=5,46,12), 
cc   >     (difny(1,9,ik), ik=10,46,15), (difnz(1,9,ik), ik=10,46,15)

cc        endif
cc            write(26,505) iday360, iyr, IT, inttra(IT), 
cc     >                    mass1, mass2, mass6, mass6a
cc 505  format(I3, 2x, 1P4D11.3, 3x, 1P3D10.2, 3x, 1P3D10.2)
cc 504  format(1X)



c
C  Now convert DIFNY, DIFNZ, DIFNYZ to mixing ratio/second for output 
C     - DIFFY(T$,L$,Z$X), DIFFZ(T$,L$,Z$X), DIFFYZ(T$,L$,Z$X), 
C
	   do 4950 IK=1,Z$X
	   do 4950 IJ=1,L$
           do 4950 IT1=1,T$ 
              diffy(IT1,IJ,IK) = -DIFNY(IT1,IJ,IK)/M(IJ,IK)
cccccc             diffyz(IT1,IJ,IK) = -DIFNYZ(IT1,IJ,IK)/M(IJ,IK)
 4950         diffz(IT1,IJ,IK) = -DIFNZ(IT1,IJ,IK)/M(IJ,IK)


        RETURN
        END
