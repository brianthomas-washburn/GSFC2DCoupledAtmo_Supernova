C **********************************************************************
C
      SUBROUTINE NWT2DD(FSP,JR,DELT,ESTEP,EPER,ITSTEP,ITDAYS,LAT,LEV,
     $ IERR,LOLA,gfaer, gfnat, gfice, TEMPAER,PRAER,ADAER, iday360, iyr)
C
C **********************************************************************
C             $Id: nwt2dd.f,v 1.30 2007/04/09 15:38:03 dkweis Exp $

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include  'com_aers.h'
      include  'com_aerd.h'
      include  'com_aerg.h'

c      include 'species.i'
c      include 'diurnal.i'
c      include 'grid.i'

      REAL gfaer, gfnat, gfice, TEMPAER, PRAER, ADAER

      DOUBLE PRECISION LO,JR,JRMOM
      DIMENSION FSP(NFSP,NTIME),JR(NJR,NTIME),DELT(NTIME)
      DIMENSION ABA(NTSP),B(NFSP)
      INTEGER LTPW(NRSP),INDBC(NRSP),NSTPD
      REAL TSTEP,TSCMIN,PSC1,PSC2,SPSC1,SPSC2,VSED1,VSED2
      REAL*8 xhcl, xclono2
      REAL SURFACE,AMASS,AH2SO4
      LOGICAL LHET
      CHARACTER*24 SPID,PLQID

      COMMON/CLOSURE/IFAM(NFAM),ICF(NFAM-1,NFP),ISF(NFSP,NFAM-1)
      COMMON /SID/SPID(NTSP),PLQID(NPLQ)
      COMMON /TSTPD/NSTPD,TSTEP,TSCMIN
      COMMON /TREAT/LTPW,INDBC
      COMMON /SULFATE/SURFACE(NLT,NHT),AMASS(NLT,NHT),AH2SO4(NLT,NHT)
      COMMON /PSC/PSC1(NLT,NHT),SPSC1(NLT,NHT),VSED1(NLT,NHT1),
     $     PSC2(NLT,NHT),SPSC2(NLT,NHT),VSED2(NLT,NHT1)
C        THE FOLLOWING 3 COMMON BLOCKS ARE USED AS TEMPORARY VARIABLES
C        TO PASS SP, KRATE, AND JRATE TO JACOB AT A GIVEN TIME OF DAY
      COMMON/SP/SPNOON(NTSP)
      COMMON/JRATE/JRMOM(NJR)
      COMMON/KRATE/R(NKR)
      COMMON/RATE/RA(NRP,NKR)
C        THE FOLLOWING 2 COMMON BLOCKS ARE USED TO PASS DATA BACK FROM JACOB
      COMMON /PROD/PR(NCSP)/LOSS/LO(NCSP)/QLOS/QL(NCSP)
      COMMON /MATRIX/A(NFSP,NFSP)
      COMMON/FLAGS/IFL(90)

C
C         PRINT OUT THE SUBROUTINE ID
      IF(IFL(42).EQ.0) THEN
         WRITE(16,*) '$Id: nwt2dd.f,v 1.30 2007/04/09 15:38:03 dkweis $'
         IFL(42)=1
      END IF
C
C         SET UP THE ERROR PARAMETER AND INITIALIZE ABA
      IERR=0
      DO I=1,NTSP
         ABA(I)=0.
      END DO


ccc        LOLA = 10
        ijp = 14
        ikp = 1
ccc
c        print *, ' '
c        print *, 'NWT2DD-1:    '
c        print *, R
c        print *, ' '
c        print *, spnoon, ihcl, ihocl, ICLNO3, IHOBR, IIN2O5,
c     >     lat,lev,iang, it, nf, ntime, nfsp
c        print *, ' '
c        print *, ' '


C    ITSTEP IS MAX NUMBER OF NEWTON ITERATIONS ALLOWED - set to 200 for ALL YEARS

         itstep = 200      ! 10
ccccccc         if (iyr .ge. 115) itstep = 200

ccccccc       print *,' itstep=',itstep


CHET
C
C        LHET IMPLIES HETEROGENEOUS REACTIONS NEED TO BE CALCULATED EVERY STEP
C
C      lhet=.false.
C      if(surface(lat,lev).gt.0. .or. spsc1(lat,lev).gt.0.
C     $      .or. spsc2(lat,lev).gt.0.) lhet=.true.

c   If HCl small (ie, <5. #den), then loop through all reactions, 
C                   if HET reaction w/ HCl (ie, type=24,25,28,32,33,42), then Rate=0.0

        lhet=.true.
        HCLBEFORE=SPNOON(IHCL)
        IF(SPNOON(IHCL).LT.5.) THEN
           do ikr=1,nkr
              itst=nint(ra(1,ikr))
              if(itst.eq.24 .or.itst.eq.25 .or. itst.eq.28 .or. 
     $      itst.eq.29 .or. itst.eq.32 .or. itst.eq.33 .or. itst.eq.42) 
     >              R(ikr)=0.
           end do
        END IF
CHET

C
*
*   ITERATE ITTT DAYS UNTIL DIURNAL PERIODIC BOUNDARY CONDITION
*   IS SATISFIED
*
c      print *,' itdays=',itdays
      ITDAYS=1
c      print *,' itdays=',itdays

      DO 170 ITTT=1,ITDAYS
cc      if(lola.gt.0) write(16,*) 'iteration day ',ittt
C       WRITE(16,*) 'INITIAL FROM ALLSP'
C       WRITE(16,*) (SPNOON(I),I=1,NFSP)
C       WRITE(16,*) 'REACTION RATE'
C       WRITE(16,*) (R(IRR),IRR=1,NKR)
*
*   CALCULATE FAST SPECIES NUMBER DENSITIES AT EACH K-TH TIME
*   DURING THE DAY, UNDER THE DIURNAL MODE
*
c      IF(LOLA.GT.1) THEN
c         WRITE(16,*) 'INITIAL TIME STEP'
c         write(16,702) 'CLNO3=',SPNOON(ICLNO3)
c         WRITE(16,702) 'HOCL=',SPNOON(IHOCL)
c         WRITE(16,702) ' HCL=',SPNOON(IHCL)
c      END IF
C

      DO 150 K=2,NTIME
c        IF(LOLA.GT.0) WRITE(16,*) 'TIME INTERVAL ',K
      IANG=K-1


      DO I=1,NFSP
         FSP(I,IANG)=SPNOON(I)

cpp         if(LAT .eq. 9)then
cpp            if(LEV .eq. 11 .or. LEV .eq. 21 .or. LEV .eq. 31)then
cpp               IF(I.eq.1)THEN
cpp         print *,' fsp=',fsp(i,iang),' iang=',iang
cpp         ENDIF
cpp            endif
cpp            endif
      END DO

c         print *,' delt=',delt(iang)

      IF(DELT(IANG).EQ.0.) GO TO 151
ccccccc      IF(DELT(IANG).EQ.0.) GO TO 150


C           GET J-RATES

      DO 7 IJR=1,NJR
          JRMOM(IJR)=JR(IJR,K)
c         print *,' ijr=(do 7)',ijr,' jrmom=',jrmom(ijr)
 7       CONTINUE

C
C       WRITE(16,*) 'J RATE FOR TIME STEP',K
C       WRITE(16,*) (JRMOM(IJR),IJR=1,NJR)
*
*   CALCULATE NEW NUMBER DENSITIES FOR ALL FAST SPECIES,
*   PROCEEDING FAMILY BY FAMILY
*
C
C   ITSTEP IS MAXIMUM NUMBER OF NEWTON ITERATIONS ALLOWED - defined above


      DO 80 IT=1,ITSTEP

         IF(LOLA.GT.1 .and. lat .eq. ijp .and. lev .eq. ikp) then
              write(16,*) '   '
              write(16,*) '   '
             WRITE(16,*) 'NEWTON ITERATION, IANG, LOLA  ',IT, IANG, LOLA
         END IF
cpp         if(LAT .eq. 9)then
cpp            if(LEV .eq. 11 .or. LEV .eq. 21 .or. LEV .eq. 31)then
cpp            print *,' iang(nwt2dd,timestep)=',iang,' it=',it
cpp            endif
cpp            endif

      DO 50 NF=2,NFAM
*
*   INITIALIZE FOR JACOBIAN MATRIX CREATED USING JUST ONE FAST FAMILY
*
       N0=IFAM(NF-1)+1
       N1=IFAM(NF)
       N=N1-N0+1
       NFM=NF-1
*
*   ZERO OUT THE JACOBIAN MATRIX TO INITIALIZE
*
      DO 17 I=N0,N1
      DO 17 II=N0,N1
   17 A(I,II)=0.
c
C         SKIP THIS FAMILY IF ALL SPECIES ARE ZERO
      DO I=N0,N1
         IF(SPNOON(I).GT.0.) GO TO 18
      END DO
      GO TO 50
   18 CONTINUE
*
*   CALCULATE THE PR(PRODUCTION), LO(LOSS) AND A(JACOBIAN) TERMS
*
CHET
C           UPDATE KRATES DUE TO CHANGES IN HCL


      IF(LHET .AND. ICF(NFM,2).EQ.ICLX) THEN
         IF(DABS((SPNOON(IHCL)-HCLBEFORE)/HCLBEFORE).GT.1.E-4) THEN
C            if(lola.gt.1)
C     $           WRITE(16,*) 'CALLING KRCALC WITHIN NWT2DD FOR ',
C     $           'GRID POINT ',LAT,LEV,' TIME STEP ',K,' ITERATION ',IT

ccccccccccccccc            CALL KRCALCH(R,LAT,LEV)

cchet           CALL KRHET(R, ADAER, TEMPAER, PRAER, gfaer,
cchet        >             WTD, gfnat, gfice, LAT, LEV, iday360)


                xhcl = SPNOON(IHCL)
                xclono2 = SPNOON(ICLNO3) 

                CALL HETCHEM(R,nkr,xhcl,xclono2,lat,lev)


C  define HCl and ClONO2 for call to HETCHEM , (LAT,LEV) need to call HETCHEM for each grid point separately

c##            write(16,*) 'From nwt2dd, called krcalch for day ',ittt,
c##     $           ' step ',k
c##            WRITE(16,*) 'REACTION RATES'
c##            WRITE(16,*) (R(IRR),IRR=1,NKR)


            HCLBEFORE=SPNOON(IHCL)

             IF(SPNOON(IHCL).LT.5.) THEN
                do ikr=1,nkr
                   itst=nint(ra(1,ikr))
                   if(itst.eq.24 .or. itst.eq.25 .or. itst.eq.28 .or.
     $     itst.eq.29 .or. itst.eq.32 .or. itst.eq.33 .or. itst.eq.42)
     $                 R(IKR)=0.
                end do
             END IF
         END IF

         hclhetloss=(R(90)*SPNOON(IHOCL)
     $        + R(91)*SPNOON(ICLNO3) + R(124)*SPNOON(IHOBR)
     $        + R(159)*SPNOON(IIN2O5) + R(161)*SPNOON(ICLNO3)
     $        + R(162)*SPNOON(IHOCL) + R(164)*SPNOON(IIN2O5)
     $        + R(166)*SPNOON(ICLNO3) + R(167)*SPNOON(IHOCL)
     $        + R(168)*SPNOON(IHOBR) + R(283)*SPNOON(IBRNO3)
     $        + R(284)*SPNOON(IBRNO3)) *DELT(IANG)             ! DELT(IANG) is just DTIMEA(1->NTIME-1)
         IF(HCLHETLOSS.GT.SPNOON(IHCL)) THEN
            FACT=SPNOON(IHCL)/HCLHETLOSS
            do ikr=1,nkr
               itst=nint(ra(1,ikr))
               if(itst.eq.24 .or. itst.eq.25 .or. itst.eq.28 .or.
     $           itst.eq.29 .or. itst.eq.32 .or. itst.eq.33 .or.
     $              itst.eq.42) R(ikr)=R(ikr)*fact
            end do

ccc      if (lat .eq. ijp  .and.  lev .eq. ikp) then
ccc         write(49,1051)
ccc         write(49,1050) iday360, lat, lev, iang, SPNOON(IHCL), 
ccc     >                  hclhetloss, fact, R(90), R(91), R(125)
ccc 1050    format(4I4, 3x, 1P6E12.3)
ccc         write(49,1051)
ccc 1051    format(4x)
ccc      endif

         END IF
      END IF
CHET


      CALL JACOB(NFM)


c     IF(LHET .AND. ICF(NFM,2).EQ.ICLX) THEN
c       IF(LOLA.GT.1) THEN
c        write(16,*) '!! TIME STEP ',K,'  ITERATION ',IT
c        WRITE(16,*) 'R162=',R(162),'  CLNO3+HCL=',
c    $        R(161)*SPNOON(ICLNO3)*DELT(IANG)
c        WRITE(16,*) 'R163=',R(163),'  HOCL+HCL=',
c    $        R(162)*SPNOON(IHOCL)*DELT(IANG)
c        write(16,702) 'CLNO3=',SPNOON(ICLNO3),'PR(CLNO3)=',
c    $        PR(ICLNO3),  !  *DELT(IANG),
c    $        '  LO(CLNO3)=',LO(ICLNO3),  ! *spnoon(iclno3)*DELT(IANG),
c    $     '  DELTA=',(pr(iclno3)-lo(iclno3)*spnoon(iclno3))*DELT(IANG)
c        WRITE(16,702) 'HOCL=',SPNOON(IHOCL),'PR(HOCL)=',
c    $        PR(IHOCL)*DELT(IANG),
c    $        '  LO(HOCL)=',LO(IHOCL)*spnoon(ihocl)*DELT(IANG),
c    $        '  DELTA=',(pr(iHOCL)-lo(ihocl)*spnoon(ihocl))*DELT(IANG)
c        WRITE(16,702) ' HCL=',SPNOON(IHCL),'PR(HCL)=',
c    $        PR(IHCL)*DELT(IANG),
c    $        '  LO(HCL)=',LO(IHCL)*spnoon(ihcl)*DELT(IANG),
c    $        '  DELTA=',(pr(ihcl)-lo(ihcl)*spnoon(ihcl))*DELT(IANG)
c 702    FORMAT(A6,1PE12.4,3X,A10,1PE12.4,2X,A14,1PE12.4,A8,1PE12.4)
c       END IF
c     END IF
*
*   DEFINE LINEAR SYSTEM MATRIX ELEMENTS(A)
*
c      if(nfm.eq.1)print *,' before do i=n0,n1(nwt2dd)',n0,n1

      DO I=N0,N1
         DO II=N0,N1
            A(I,II)=-DELT(IANG)*A(I,II)
         END DO
      END DO
      DO 10 I=N0,N1
         A(I,I)=A(I,I)+1.0
         AUX=LO(I)*SPNOON(I)-PR(I)
         B(I)=-DELT(IANG)*AUX-SPNOON(I)+FSP(I,IANG)
   10 CONTINUE
c      if(nfm.eq.1)print *,' after 10 continue(nwt2dd)'
      IF(LOLA.GT.1 .and. lat.eq.ijp .and. lev.eq.ikp .and. nf.eq.2) THEN
         WRITE(16,*) '    SPECIES DENSITIES FOR FAMILY ',NF-1
         WRITE(16,6001) (SPID(JF),SPNOON(JF),JF=N0,N1)
         WRITE(16,*) '    PRODUCTION TERMS FOR FAMILY ',NF-1
         WRITE(16,6001) (SPID(JF),PR(JF),JF=N0,N1)
         WRITE(16,*) '    LOSS TERMS FOR FAMILY ',NF-1
         WRITE(16,6001) (SPID(JF),LO(JF),JF=N0,N1)
         WRITE(16,*) '    B AFTER JACOB FOR FAMILY ',NF-1
         WRITE(16,6001) (SPID(JF),B(JF),JF=N0,N1)
 6001    FORMAT(3(2X,A8,'=',E20.12))
         write(16,*) '    A Matrix for family ',nf-1
         do ii=n0,n1
            write(16,7007) (a(i,ii),i=n0,n1)
         end do
 7007    format(10E11.4)
      END IF
*
*   TEST THE JACOBIAN FOR SINGULARITY AND RETURN THE CHANGES IN NUMBER
*   DENSITIES(DELTA SP), PLACING THEM IN THE LAST COLUMN OF ARRAY A
*   AS SOLUTIONS TO THE SYSTEM  A!* DELTA SP!= -F!.
*
C          lineqn subroutine inserted here
      SUMA=0.
      DO 20 J=N0,N1
      DO 20 I=N0,N1
   20 SUMA=SUMA+DABS(A(I,J))
      TOLER=(SUMA/(N*N))*1.E-12
c      if(nfm.eq.1)print *,' before do 300(nwt2dd)'
      DO 300 M=N0,N1
      MP1=M+1
      TEMP=0.
      ITEMP=M
      DO 100 I=M,N1
      IF(DABS(A(I,M)).GT.DABS(TEMP)) THEN
         TEMP=A(I,M)
         ITEMP=I
      END IF
  100 CONTINUE
      IF(DABS(TEMP).LE.TOLER) THEN
*!!!!!!!!!!     ERROR:  SINGULAR JACOBIAN
         IERR=3
cc         WRITE(16,*) 'SINGULAR MATRIX FOR FAMILY',NF-1,'  STEP',IT,
cc     .        '  TIME',K,'  DAY',ITTT

         IF (IT .eq. 1  .or.  IT .eq. 10) then
           WRITE(16,505) NF-1, IT, K, lat, lev, iday360
 505  format('SINGULAR MATRIX FOR FAMILY', I3,'    STEP',I3, 
     >  '       TIME', I3, '          LAT,  Lev,  iday360 =', 3I4)
ctol         if (lev .ge. 5) RETURN
         END IF
ccc         RETURN
      END IF
      II=ITEMP-M
c      if(nfm.eq.1)print *,' before if(ii..(nwt2dd)',ii
      IF(II.EQ.0) THEN
         DO 130 I=M,N1
  130    A(M,I)=A(M,I)/TEMP
         B(M)=B(M)/TEMP
      ELSE
         DO 110 I=M,N1
            ATEMP=A(M,I)
            A(M,I)=A(M+II,I)
            A(M+II,I)=ATEMP
  110    A(M,I)=A(M,I)/TEMP
         BTEMP=B(ITEMP)
         B(ITEMP)=B(M)
         B(M)=BTEMP/TEMP
      END IF
      IF(M.NE.N1) THEN
         DO 210 I=MP1,N1
         DO 210 J=MP1,N1
  210    A(I,J)=A(I,J)-A(I,M)*A(M,J)
         DO 215 I=MP1,N1
  215    B(I)=B(I)-A(I,M)*B(M)
      END IF
  300 CONTINUE
      DO 400 M=N1-1,N0,-1
      DO 400 L=N1,M+1,-1
  400 B(M)=B(M)-A(M,L)*B(L)
c            end of subroutine lineqn
*
*   UPDATE NUMBER DENSITIES IN WORKING ARRAY SPNOON
*
c      if(nfm.eq.5)print *,' before do 30(nwt2dd)',n0,n1
      DO 30 I=N0,N1

C  MAKE SURE THAT NEW NUMBER DENSITIES WILL BE POSITIVE

cpp         if(LAT .eq. 9)then
cpp            if(LEV .eq. 11 .or. LEV .eq. 21 .or. LEV .eq. 31)then
c               IF(NFM.EQ.5 .and. I.eq.23)THEN
c         print *,' before spnoon(i)=',spnoon(i),' i=',i
cpp               IF(NFM.EQ.1 .and. I.eq.1)THEN
cpp         print *,' before spnoon(i)=',spnoon(i),' i=',i,' pr=',
cpp     *  pr(1),' lo=',lo(1)
cpp         ENDIF
cpp            endif
cpp            endif

      IF(B(I)+SPNOON(I).LE.0.) B(I)=-SPNOON(I)/2.0
      SPNOON(I)=B(I)+SPNOON(I)

cpp         if(LAT .eq. 9)then
cpp            if(LEV .eq. 11 .or. LEV .eq. 21 .or. LEV .eq. 31)then
cpp               IF(NFM.EQ.1 .and. I.eq.1)THEN
cpp         print *,' spnoon(i)=',spnoon(i),' i=',i,' b(i)=',b(i)
c               IF(NFM.EQ.5 .and. I.eq.23)THEN
c         print *,' spnoon(i)=',spnoon(i),' i=',i,' b(i)=',b(i)
cpp               ENDIF
cpp            endif
cpp            endif

 30   CONTINUE
c   30 SPNOON(I)=B(I)+SPNOON(I)

      IF(LOLA.GT.1 .and. lat.eq.ijp .and. lev.eq.ikp .and. nf.eq.2) THEN
         WRITE(16,*) '    DELTA SPECIES DENSITIES FOR FAMILY',NF-1
         WRITE(16,6001) (SPID(JF),B(JF),JF=N0,N1)
         WRITE(16,*) '    NEW SPECIES DENSITIES FOR FAMILY',NF-1
         WRITE(16,6001) (SPID(JF),SPNOON(JF),JF=N0,N1)
      END IF
c$$$*
c$$$*   CHECK CLOSURE SPECIES NUMBER DENSITIES FOR MASS CONSERVATION
c$$$*
c$$$      DO 199 I=NFM,NFM1
c$$$      ICFI=ICF(NFM,1)
c$$$      IF(ICFI.NE.0) THEN
c$$$C        THIS IS A CLOSURE FAMILY
c$$$C          IF DOING ALL FAMILIES AT ONCE, USE N2,N3 FOR LOOP 192
c$$$C        N2=IFAM(I)+1
c$$$C        N3=IFAM(I+1)
c$$$         SSUM=0.0
c$$$         DO 192 II=N0,N1
c$$$  192    SSUM=SSUM+ISF(II,NFM)*SPNOON(II)
c$$$         IF(ICFI.GT.1) THEN
c$$$C          THERE ARE EXTERNAL SPECIES
c$$$            DO 194 II=2,ICFI
c$$$            IEXSP=ICF(NFM,II+1)
c$$$  194       SSUM=SSUM+ISF(IEXSP,NFM)*SPNOON(IEXSP)
c$$$         END IF
c$$$         SPFAM=SPNOON(ICF(NFM,2))
c$$$  198    IF(DABS(SSUM-SPFAM)/SPFAM .GT. 2.0E-2) then
c$$$            write(16,*) 'mass conservation violated: family',nf-1,
c$$$     .           ' time step',k,' day',ittt
c$$$            write(16,*) 'ssum=',ssum,'  spfam=;',spfam
c$$$         end if
c$$$      END IF
c$$$  199 CONTINUE

   50 CONTINUE
*        END OF LOOP OVER FAMILIES
C
*
*   CHECK RELATIVE ERRORS ABA WITH CONVERGENCE TOLERANCE ESTEP
*
      DO 51 I=1,NFAM-1
         CALL JACOB(I)
   51 CONTINUE
      DO 52 I=1,NFSP
      IF(SPNOON(I).LT.1.0E-30) THEN
         ABA(I)=1.0E-10
      ELSE
         AUX=LO(I)*SPNOON(I)-PR(I)
         B(I)=-DELT(IANG)*AUX-SPNOON(I)+FSP(I,IANG)
         ABA(I)=DABS(B(I)/(SPNOON(I)+LO(I)*SPNOON(I)*DELT(IANG)))
      END IF
   52 CONTINUE
      IF(LOLA.GT.1 .and. lat.eq.ijp .and. lev.eq.ikp .and. nf.eq.2) THEN
        write(16,6670) eSTEP,(aba(i),i=1,nfsp)
 6670   format('End of NEWTON iterations, checking for convergence,',
     .       ' error cutoff=',f10.5/'aba(1-nfsp)=',11f10.5/
     .       2x,12f10.5/2x,12f10.5/2x,12f10.5/2x,12f10.5)
      END IF
      DO 53 I=1,NFSP
      IF(ABA(I).GE.ESTEP)GO TO 80
   53 CONTINUE
C         THE ACCURACY IS O.K., QUIT NEWTON ITERATIONS
      IF(LOLA.GT.0 .and. lat .eq. ijp .and. lev .eq. ikp .and. nf.eq.2)  
     >    WRITE(16,6162) ITTT,K,IT
 6162 FORMAT(' DAY',I5,' TIME INTERVAL',I3,': CONVERGENCE OBTAINED',
     .     ' IN',I5,' ITERATIONS')
      GO TO 85
   80 CONTINUE
*        END OF LOOP OVER MAXIMUM NUMBER OF NEWTON ITERATIONS
*!!!!!!!!!!!!   ERROR:  TOO MANY NEWTON ITERATIONS
      IERR=1
      IF(LOLA.NE.0 .and. lat .eq. ijp .and. lev .eq. ikp .and. nf.eq.2) 
     >    WRITE(16,6161) ITTT,K,IT-1
 6161 FORMAT(' DAY',I5,' TIME INTERVAL',I3,': CONVERGENCE NOT OBTAINED',
     .     ' IN',I5,' ITERATIONS')
*
*   UPDATE NUMBER DENSITIES SP WITH THE WORKING ARRAY VALUES IN SPNOON
*
   85 CONTINUE
*
*
*   RESCALE NOX FAMILY - NO LONGER NEEDED w/ HNO3 as a fast species,  and this is NOW DONE in UNFILLA
*
ccnox      DO 199 NFM=1,NFAM-1
ccnox      ICFI=ICF(NFM,1)
ccnox      ICF2=ICF(NFM,2)
ccnox      IF(ICF2.NE.INOX) GO TO 199
C        THIS IS A CLOSURE FAMILY
ccnox      N2=IFAM(NFM)+1
ccnox      N3=IFAM(NFM+1)
ccnox      SSUM=0.0
ccnox      DO 192 II=N2,N3
ccnox  192 SSUM=SSUM+ISF(II,NFM)*SPNOON(II)
ccnox      SPFAM=SPNOON(ICF2)
ccnox      IF(ICFI.GT.1) THEN
ccnox         SPEXT=0.
ccnox         DO 194 II=2,ICFI
ccnox            IEXSP=ICF(NFM,II+1)
ccnox  194       SPEXT=SPEXT+ISF(IEXSP,NFM)*SPNOON(IEXSP)
ccnox      END IF
ccnox      if(ssum.lt.1.e-20) then
ccnox         write(16,*) 'ssum=0. for inox, day ',ittt,
ccnox     .     ' time step ',k
ccnox         stop 'ssum=0'
ccnox      end if
ccnox      if(SPFAM.LT.SPEXT) then
ccnox         SPNOON(ICF2)=1.01*SPEXT
ccnox         SPFAM=0.01*SPEXT
ccnox      ELSE
ccnox         SPFAM=SPFAM-SPEXT
ccnox      end if
ccnox      DO II=N2,N3
ccnox
cpp         if(LAT .eq. 9)then
cpp            if(LEV .eq. 11 .or. LEV .eq. 21 .or. LEV .eq. 31)then
cpp               IF(NFM.EQ.1 .and. I.eq.1)THEN
cpp         print *,' spnoon(ii)=',spnoon(ii),' pr=',pr(1),' lo=',lo(1)
c               IF(NFM.EQ.5 .and. I.eq.23)THEN
c         print *,' spnoon(ii)=',spnoon(ii)
cpp               ENDIF
cpp            endif
cpp            endif
ccnox
ccnox         SPNOON(II)=SPNOON(II)*SPFAM/SSUM
ccnox
cpp         if(LAT .eq. 9)then
cpp            if(LEV .eq. 11 .or. LEV .eq. 21 .or. LEV .eq. 31)then
cpp               IF(NFM.EQ.1 .and. I.eq.1)THEN
cpp         print *,' spnoon(ii)=',spnoon(ii),' spfam=',spfam,' ssum=',
cpp     *  ssum
c               IF(NFM.EQ.5 .and. I.eq.23)THEN
c         print *,' spnoon(ii)=',spnoon(ii),' spfam=',spfam,' ssum=',
c     *  ssum
cpp               ENDIF
cpp            endif
cpp            endif
ccnox
ccnox      END DO
ccnox  199 CONTINUE
C
 151  CONTINUE


C  now load in AER HET REACTIONS into GSFC KH array for use in SOLVER/OX

          CALL HETMAP(R,nkr,lat,lev,iang,gfaer, gfnat, gfice)

  150 CONTINUE
C                   END OF LOOP OVER TIME STEPS WITHIN THE DAY

      DO I=1,NFSP
         FSP(I,NTIME)=SPNOON(I)
      END DO


CHET
C  need to call HETMAP once more to update the KH array for the NOONTIME value - just use the 
C     current values in the R array; otherwise the KH array will be 0.0 at diurnal time step, iang=18
C     (this is just for the KH array used in SOLVER/OX - the R array is not used at all past this point)

         iang0 = 18
         CALL HETMAP(R,nkr,lat,lev,iang0,gfaer, gfnat, gfice)
CHET

ccc      if (ij .eq. ijp  .and.  ik .eq. ikp) write(49,1051)
ccc      if (ij .eq. ijp  .and.  ik .eq. ikp) write(49,1051)

C
c      print *,' lola(before call plq2dd)=',lola

cpp         if(LAT .eq. 9)then
cpp            if(LEV .eq. 11 .or. LEV .eq. 21 .or. LEV .eq. 31)then
cpp               IF(NFM.EQ.1 .and. I.eq.1)THEN
cpp         print *,' before call plq2dd spnoon(i)=',spnoon(i)
c               IF(NFM.EQ.5 .and. I.eq.23)THEN
c         print *,' before call plq2dd spnoon(i)=',spnoon(i)
cpp               ENDIF
cpp            endif
cpp            endif

C
ccplq       CALL PLQ2DD(FSP,JR,DELT,LTPW,LAT,LEV,LHET,LOLA)

c##      CALL DO3PLQ(FSP,JR,DELT,LOLA,0)
C         UPDATE SLOW SPECIES CONCENTRATIONS IF IN EQUILIBRIUM
ccplq      ISP=IHNO3
ccplq      IF(LO(ISP).GT.TSCMIN .AND. INDBC(ISP-NFSP).GT.0) THEN
ccplq         SP1=PR(ISP)/LO(ISP)
ccplq         IF(SPNOON(ISP).LT.1.0E-30) THEN
ccplq            ABA(ISP)=1.0E-10
ccplq         ELSE
ccplq            ABA(ISP)=(SP1-SPNOON(ISP))/SPNOON(ISP)
ccplq         END IF
ccplq         SPNOON(ISP)=SP1
ccplq         SPNOON(INOZ)=SPNOON(INOX)-SPNOON(IHNO3)
ccplq         CALL RSCNCL(SPNOON)
ccplq         IF(LOLA.NE.0) write(16,668) spid(isp),spnoon(isp),ittt
ccplq  668    format('LIN. EQUILIBRIUM FOR ',A8,E10.4,
ccplq     .        ' ITERATION DAY ',I4)
ccplq      END IF
ccplq      ISP=IO3
ccplq      IF((LO(ISP)+SPNOON(ISP)*2.*QL(ISP)).GT.TSCMIN .AND.
ccplq     .     INDBC(ISP-NFSP).GT.0) THEN
ccplq         D=LO(ISP)*LO(ISP) + 4.*PR(ISP)*QL(ISP)
ccplq         IF(D.LT.0.D0) STOP 'QUADRATIC EQUILIBRIUM EQUATION'
ccplq         SP1=(-LO(ISP)+DSQRT(D))/(2.*QL(ISP))
ccplq         IF(SPNOON(ISP).LT.1.0E-30) THEN
ccplq            ABA(ISP)=1.0E-10
ccplq         ELSE
ccplq            ABA(ISP)=(SP1-SPNOON(ISP))/SPNOON(ISP)
ccplq         END IF
ccplq         SPNOON(ISP)=SP1
ccplq         IF(LOLA.NE.0) write(16,669) spid(isp),spnoon(isp),ittt
ccplq  669    format('QUAD. EQUILIBRIUM FOR ',A8,E10.4,
ccplq     .        ' ITERATION DAY ',I4)
ccplq      END IF
C         NEED TO UPDATE FSP IF DID RSCNCL AFTER SLOW SPECIES EQUILIBRIUM
ccplq      DO I=1,NFSP
ccplq         FSP(I,NTIME)=SPNOON(I)
ccplq      END DO
C
cpp         if(LAT .eq. 9)then
cpp            if(LEV .eq. 11 .or. LEV .eq. 21 .or. LEV .eq. 31)then
cpp               IF(NFM.EQ.1 .and. I.eq.1)THEN
cpp         print *,' before if(lola..spnoon(i)=',spnoon(i)
c               IF(NFM.EQ.5 .and. I.eq.23)THEN
c         print *,' before if(lola..spnoon(i)=',spnoon(i)
cpp               ENDIF
cpp            endif
cpp            endif

      IF(LOLA.NE.0 .and. lat.eq.ijp .and. lev.eq.ikp .and. nf.eq.2) THEN
         WRITE(16,*) 'NOONTIME DENSITIES FOR ITERATION DAY ',ITTT
         WRITE(16,6101) (SPID(ISP),FSP(ISP,NTIME),ISP=1,NFSP)
 6101    FORMAT(5(2X,A8,'=',E12.4))
      END IF
      DO 153 I=NFSP+1,NCSP
      IF(ABA(I).GE.ESTEP)GO TO 170
  153 CONTINUE
C
C   CHECK PRESENT DAY'S NOON-TIME DENSITIES WITH PREVIOUS DAY'S
C   FOR PERIODICITY
C
      DO 155 I=1,NFSP
      SP1=FSP(I,1)
      IF(SP1.LT.1.D-2) GO TO 155
      IF(DABS(SPNOON(I)-SP1)/SP1.GT.EPER) GO TO 170
  155 CONTINUE
      IF(LOLA.NE.0 .and. lat .eq. ijp .and. lev .eq. ikp .and. nf.eq.2) 
     >       WRITE(16,*) 'PERIODIC SOLUTION IN',ITTT,' DAYS'
      RETURN
C
  170 CONTINUE
*          END OF LOOP OVER ITDAYS:  NO PERIODIC SOLUTION FOUND
*!!!!!!!!!!!!!!    ERROR:  NO PERIODIC SOLUTION
      IERR=2
      IF(LOLA.NE.0 .and. lat .eq. ijp .and. lev .eq. ikp .and. nf.eq.2) 
     >      WRITE(16,*) 'NO PERIODIC SOLUTION IN',ITDAYS,' DAYS'
      RETURN
C
      END
