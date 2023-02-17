C **********************************************************************
C
      SUBROUTINE TIMEAV(SP,NSP,DTIMEA,I,J, TFD, L$)
c      SUBROUTINE TIMEAV(SP,NSP,DTIME,DLITE)
C
C **********************************************************************
C             @(#)timeav.f	1.5  06/29/99

      include  'com_aerd.h'
c      include 'diurnal.i'


      DOUBLE PRECISION SP(NSP,NTIME+3),DTIMEA(NTIME),SPDAY,SPNITE
c      DOUBLE PRECISION SP(NSP,NTIME+3),DTIME(NTIME),SPDAY,SPNITE

      INTEGER L$
      REAL TFD(L$)

      COMMON/FLAGS/IFL(90)

C
C         WRITE THE SUBROUTINE ID THE FIRST TIME AROUND
      IF(IFL(44).EQ.0) THEN
         WRITE(16,*) '@(#)timeav.f	1.5  06/29/99'
         IFL(44)=1
      END IF
C

      HDAY=TFD(I)*86400.
      HNITE=(1.-TFD(I))*86400.

ccp         if(I .eq. 9)then
ccp            if(J .eq. 11 .or. J .eq. 21 .or. J .eq. 31)then
ccp               print *,' tfd(timeav)=',tfd(i),' hday=',hday,
ccp     *  ' hnite=',hnite
ccp            endif
ccp            endif
      

C      HDAY=DLITE*86400.
C      HNITE=(1.-DLITE)*86400.
cvd$ novector
c      print *,' sp(1,21,23)=',sp(1,1),sp(21,1),sp(23,1)
c      print *,' sp(2,23)=',sp(2,1),sp(23,1)

ccp         if(I .eq. 9)then
ccp            if(J .eq. 11 .or. J .eq. 21 .or. J .eq. 31)then
ccp               print *,' sp(1)=',sp(1,1)
c               print *,' sp(23)=',sp(23,1)
ccp            endif
ccp            endif

      DO 150 ISP=1,NSP
      SPDAY=0.
      SPNITE=0.
cvd$ novector
      DO IT=1,SUNSET-1

ccp         if(I .eq. 9)then
ccp            if(J .eq. 11 .or. J .eq. 21 .or. J .eq. 31)then
ccp         if(isp.eq.1)print *,' spday=',spday
c         if(isp.eq.23)print *,' spday=',spday
ccp            endif
ccp            endif

         SPDAY=SPDAY+(SP(ISP,IT)+SP(ISP,IT+1))*0.5*DTIMEA(IT)

c         SPDAY=SPDAY+(SP(ISP,IT)+SP(ISP,IT+1))*0.5*DTIME(IT)


ccp         if(I .eq. 9)then
ccp            if(J .eq. 11 .or. J .eq. 21 .or. J .eq. 31)then
ccp         if(isp.eq.1)then
ccp            print *,' spday=',spday,' sp(isp,it)=',sp(isp,it),
ccp     *  ' sp(isp,it+1)=',sp(isp,it+1),' dtimea=',dtimea(it),' it=',it
c         if(isp.eq.23)then
c            print *,' spday=',spday,' sp(isp,it)=',sp(isp,it),
c     *  ' sp(isp,it+1)=',sp(isp,it+1),' dtimea=',dtimea(it),' it=',it
ccp            endif
ccp            endif
ccp            endif


      END DO
      DO IT=SUNSET,SUNRIS-1

ccp         if(I .eq. 9)then
ccp            if(J .eq. 11 .or. J .eq. 21 .or. J .eq. 31)then
ccp         if(isp.eq.1)print *,' spnite=',spnite
c         if(isp.eq.23)print *,' spnite=',spnite
ccp            endif
ccp            endif


         SPNITE=SPNITE+(SP(ISP,IT)+SP(ISP,IT+1))*0.5*DTIMEA(IT)
c         SPNITE=SPNITE+(SP(ISP,IT)+SP(ISP,IT+1))*0.5*DTIME(IT)

ccp         if(I .eq. 9)then
ccp            if(J .eq. 11 .or. J .eq. 21 .or. J .eq. 31)then
ccp        if(isp.eq.1)then
ccp            print *,' spnite=',spnite,' sp(isp,it)=',sp(isp,it),
ccp     *  ' sp(isp,it+1)=',sp(isp,it+1),' dtimea=',dtimea(it),' it=',it
c        if(isp.eq.23)then
c            print *,' spnite=',spnite,' sp(isp,it)=',sp(isp,it),
c     *  ' sp(isp,it+1)=',sp(isp,it+1),' dtimea=',dtimea(it),' it=',it
ccp            endif
ccp            endif
ccp            endif


      END DO
      DO IT=SUNRIS,NTIME-1

ccp         if(I .eq. 9)then
ccp            if(J .eq. 11 .or. J .eq. 21 .or. J .eq. 31)then
ccp         if(isp.eq.1)print *,' spday=',spday
c         if(isp.eq.23)print *,' spday=',spday
ccp            endif
ccp            endif


         SPDAY=SPDAY+(SP(ISP,IT)+SP(ISP,IT+1))*0.5*DTIMEA(IT)
c         SPDAY=SPDAY+(SP(ISP,IT)+SP(ISP,IT+1))*0.5*DTIME(IT)

ccp         if(I .eq. 9)then
ccp            if(J .eq. 11 .or. J .eq. 21 .or. J .eq. 31)then
ccp         if(isp.eq.1)then
ccp            print *,' spday=',spday,' sp(isp,it)=',sp(isp,it),
ccp     *  ' sp(isp,it+1)=',sp(isp,it+1),' dtimea=',dtimea(it),' it=',it
c         if(isp.eq.23)then
c            print *,' spday=',spday,' sp(isp,it)=',sp(isp,it),
c     *  ' sp(isp,it+1)=',sp(isp,it+1),' dtimea=',dtimea(it),' it=',it
ccp            endif
ccp            endif
ccp            endif


      END DO
  100 CONTINUE
      IF(HDAY.GT.0.) THEN
         SP(ISP,DAYAVG)=SPDAY/HDAY
      ELSE
         SP(ISP,DAYAVG)=SP(ISP,1)
      END IF

ccp         if(I .eq. 9)then
ccp            if(J .eq. 11 .or. J .eq. 21 .or. J .eq. 31)then
ccp         if(isp.eq.1)then
ccp            print *,' sp(isp,dayavg)=',sp(isp,dayavg),
ccp     *  ' spday=',spday,' hday=',hday,' sp(isp,1)=',sp(isp,1),
ccp     *  ' dayavg=',dayavg
c         if(isp.eq.23)then
c            print *,' sp(isp,dayavg)=',sp(isp,dayavg),
c     *  ' spday=',spday,' hday=',hday,' sp(isp,1)=',sp(isp,1),
c     *  ' dayavg=',dayavg
ccp            endif
ccp            endif
ccp            endif


      IF(HNITE.GT.0.) THEN
         SP(ISP,NITAVG)=SPNITE/HNITE
      ELSE
         SP(ISP,NITAVG)=SP(ISP,MIDNIT)
      END IF
      SP(ISP,DIURAVG)=(SPDAY+SPNITE)/86400.


ccp         if(I .eq. 9)then
ccp            if(J .eq. 11 .or. J .eq. 21 .or. J .eq. 31)then
ccp        if(isp.eq.1)then
ccp            print *,' sp(isp,nitavg)=',sp(isp,nitavg),
ccp     *  ' spnite=',spnite,' hnite=',hnite,' sp(isp,diurnag)=',
ccp     *  sp(isp,diuravg),' nitavg=',nitavg,' diuravg=',diuravg,
ccp     *  ' ntime=',ntime
c         if(isp.eq.23)then
c            print *,' sp(isp,nitavg)=',sp(isp,nitavg),
c     *  ' spnite=',spnite,' hnite=',hnite,' sp(isp,diurnag)=',
c     *  sp(isp,diuravg),' nitavg=',nitavg,' diuravg=',diuravg
ccp            endif
ccp            endif
ccp            endif

  150 CONTINUE
      RETURN
      END

