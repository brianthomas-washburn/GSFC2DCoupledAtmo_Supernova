      SUBROUTINE OUTA(N,NT,NDAY,NYEAR,X,Y,A,N3$,M3$,isw6)

cjer use of isw6 switch is no longer used
cjer  frequency of output determined in main by iwrite

      INCLUDE 'PARAM.INC'
c      COMMON /DUMP_1/ IDUMP,IFILE,isw2(isw$)
      common /dump_2/ iwrite,yrjer
      DIMENSION A(N3$,M3$),X(N3$),Y(M3$)

      INCLUDE 'timer.inc'

      INTEGER IFREQ

C     GENERALIZED OUTPUT ROUTINE FOR ARRAYS 
C
C     N = THE LABEL (A NUMBER), NT = THE STEP
C     NDAY,NYEAR THE DAY AND YEAR add 1 to N in these routines

      R=FLOAT(N)
      DAY=NDAY
      YEAR=NYEAR
      STEP=NT
      R1=N3$
      R2=M3$


c nt is the number of dynamics time steps in a chemical time step

      IF ( iwrite .eq. 1 .and. nt .eq. 1) THEN

        WRITE(10) R,DAY,yrjer,STEP,R1,R2
        WRITE(10) X,Y,A
 

      ENDIF

      RETURN



      END


      SUBROUTINE OUTP(N,NT,NDAY,NYEAR,X,Y,A,N3$,M3$,K3$,isw6)

      INCLUDE 'PARAM.INC'
      COMMON /DUMP_1/ IDUMP,IFILE,isw2(isw$)
      common /dump_2/ iwrite,yrjer
      DIMENSION A(N3$,M3$,K3$),X(N3$),Y(M3$)

      INCLUDE 'timer.inc'

      INTEGER IFREQ


C     GENERALIZED OUTPUT ROUTINE FOR ARRAYS 
C
C     N = THE LABEL (A NUMBER), NT = THE STEP
C     NDAY,NYEAR THE DAY AND YEAR add 1 to N in these routines
c     isw is a switch array read in from isw.dat
c     if isw = 0 then values are output based on IDUMP
c     if isw=-1 then the value is not output
c     if isw >0 then the value is output modulo with isw instead of
c       idump
c      numbers for variables (add 1 to index)
c     0)='HEATING RATE K/DAY 
c     1)='COOLING RATE K/DAY   
c     2)='GW DECELER. M/S/DAY  
c     3)='PW DECELER. M/S/DAY  
c     4)='PSI CM^2/S           
c     5)='VS CM/S              
c     6)='WS CM/S              
c     7)='U-BAR M/S            
c     8)='TEMP. K              
c     9)='KYY CM^2/S           
c     10)='Q BAR Y 1/CM/S      
c     11)='KZZ CM^2/S          
c     12)='RALEIGH FRICT. 1/S  
c     13)='NEWTONIAN COOL 1/S  
c     14)='GLOBAL MEAN T K     
c     15)='PW AMP CM2/S2       
c     16)='DTHDZ 1/CM/S        
c     17)='SOURCE   X/S        
c     18)='2-D Newtonian Cooling'    
c     19)='VERT FLUX X*CM/S    
c     20)='HOR FLUX CONVG X/S  
c     21)='VERT FLUX CONV X/S  
c     22)='CONSTITUENT         
c     23)='VELOCITY VECTORS    
C     24)=' OZONE MIXING RATIO USED IN RADIATION
C     25)=' H2O USED IN RADIATION
C
C   FOR CONSTITUENTS, CONSTITUENT NUMBER IS M (10 constituents allowed)
C
C  N = 130+M ='Source X/t       ' 
C  N = 140+M ='lifetime 1/t   ' 
C  N = 150+M ='unused   ' 
C  N = 180+M ='CONSTITUENT        ' 

C   FOR PLANETARY WAVES, WAVE NUMBER IS MW (10 wavenumbers allowed)

c  N = 200+mw ='Wave Amplitude     
c  N = 210+MW ='WAVE PHASE         '
C  N = 220+MW ='HORIZ EP FLUX     
C  N = 230+MW ='VERTICAL EP FLUX
C  N = 240+MW ='PV FLUX '
C     
C  For Gravity Waves,phase speed +400

C  N = 99 IS END OF TIME LOOP marker

C   Variable output scheme:
C     Isw is an array which controls how often the variable
c is output  using the modulo function with the step. ISW is
c dimensioned 26+m+mw  and  the field is output based on the
c isw value. If isw=1 the field is output  every step, -1 -
c never output, 10 every 10 steps. If isw is 0 then isw6 or
c idump is used to determine the output frequency.  Isw6 takes 
c priority over idump. The Isw array is read  in
c during the program start. The array switches.dat contains a
c second  modulo function for output, idump, this is the global
c function which gives the  minimum modulo when isw6 and 
c isw(j) = 0. Ex/ if idump = 10, isw6=0 and isw(j)=1 then 1 is
c used. If idump=10 and isw(j)=-1 then no dump occurs. If 
c isw(j)=0 and idump =10 but isw6=1 then 1 is used.
c  To summarize, the priority is
c    isw6     highest  0 means check the next priority
c    isw(j)   second   0 means check next priority
c    idump    lowest
        n1=n
        if((n1.le.25).or.(n1.ge.50)) n2=n1

        if(n1.gt.100.and.n1.lt.200) then 
c       constituent
        ix5=n1*.1        
        ncn=n1-ix5*10
c       ncn is the second and third digit of n1
        n2=26+ncn
        endif

        if((n1.gt.200).and.(n1.lt.300)) then 
c       planetary waves
        ix5=n1*.1
        ncp=n1-ix5*10
c       ncp is the second and third digit of n1
        n2=30 + ncp
        endif

        if(n1.gt.300) then
c       gravity waves, switch on the zonal wind
        n2=3
        endif

        is2=idump
        if (n2.le.isw$) is2=isw2(n2)
        if (isw6.gt.0) is2=min0(is2,isw6)
c       pick the min modulo
        if ((isw6.eq.0).and.(is2.eq.0)) is2=idump
c         print 888,is2,isw6,n,idump

 888	format(' check from outa, is2,isw6,n,idump',4i5)
c        no output of variable

c	IF (N.EQ.99) GO TO 77
	if (is2.lt.0) return 
      If (is2.eq.0) return

C     IFREQ added 3/26/97 by PEM to use new timing.  DYSTP initialized
C     to ISW(6) by MAIN call to TIME_INIT.

      IFREQ = IS2/ DYSTPD

C      IF (MOD(NT,Is2).NE.0) RETURN
cjer      IF ((MOD(NDAY, IFREQ).NE.0) .OR. (NT .GT. 1)) RETURN
      IF ((MOD(idoy360, IFREQ).NE.0) .OR. (NT .GT. 1)) RETURN

C         print 889,n,is2,nt
C 889     format(' outputting file ',i5,' is2,nt =',2i5)
 889     format(' outputting file ',i3,' DAY=',i4,' YEAR=',I2 )
77      IF (N.NE.99) GO TO 45 
      IFILE=IFILE+1
45    R=FLOAT(N)+.1
      DAY=NDAY
      YEAR=NYEAR
      STEP=NT
      R1=N3$
      R2=M3$
      R3=K3$

C     If conditional on WRITE block added 6/28/94 by P. E. Meade ...
      IF (iwrite .eq. 1) THEN

         WRITE(110) R,DAY,YEAR,STEP,R1,R2
         WRITE(110) R3
         WRITE(110) X,Y,A

      ENDIF

      RETURN
      END
