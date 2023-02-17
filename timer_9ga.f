
C     >>>>>>>>>>>>>>>  timer.f  <<<<<<<<<<<<<<<

C     timer.f is a new timing routine for the Goddard 2D coupled model
C     which handles the internal timing for both the dynamics and
C     chemistry parts of the model.  SUBROUTINE TIME_INIT is called to
C     initialize the timing, while SUBROUTINE TIME_TICK is called to update
C     the timing.

C     Time is kept in terms of an idealized 360-day year, with the effective
C     day length increased so that the total number of seconds in a year
C     is the same as for a real 365-day year.  Leap years are not included.
C     Corresponding day numbers, etc., for a real 365-day year are calculated
C     from the 360-day calendar values.

C     Fundamental timekeeping is done in terms of a decimal day number, so
C     so that fundamental timesteps corresponding to a fraction of a day
C     can be correctly treated.

C     timer.f also contains FUNCTION CALC_DECL, which calculates the
C     solar declination (in degrees) as a function of the day of the year.
C     (365-day real calendar).  This function is included here since it
C     depends only on the internal timekeeping of the model.  SUBROUTINE
C     TIME_INIT and SUBROUTINE TIME_TICK call FUNCTION CALC_DECL to update
C     the solar declination variable SOLDEC, which is contained in
C     COMMON /SOLAR/ in timer.inc.

C     The timer.f routines modify the timing variables in the COMMON blocks
C     in timer.inc.  See that file for explanation of COMMON block variables.



C     SUBROUTINE TIME_INIT:  
C        IYR0:      Initial model year number (INTEGER).
C        IDOY0:     Starting day of the year, 360-day calender (INTEGER).
C        CHSTPD0:   Number of fundamental (chemistry) timesteps per day
C                   (INTEGER).
C        DYSTPD0:   Number of dynamics timesteps per day (INTEGER).


      SUBROUTINE TIME_INIT(IYR0, IDOY0, CHSTPD0, DYSTPD0, IDOYST0)

      INTEGER IYR0, IDOY0, CHSTPD0, DYSTPD0, IDOYST0

      include "timer.inc"

      REAL    DAYL
      DATA    DAYL /86400./

      INTEGER MNTHDOY(13)
      DATA    MNTHDOY/1, 32, 60, 91, 121, 152, 182,
     2                213, 244, 274, 305, 335, 366/

      REAL    DEGTORAD
      DATA    DEGTORAD /0.017453293/     !     PI/180.0


C     If the number of timesteps per day is less than 1, reset it to 1 and
C     assign to CHSTPD.  Calculate the fundamental timestep in terms of a
C     fraction of a day (DELT_DAY) for the 360-day idealized calender, and in
C     terms of seconds (DELT_SEC).  Note that DELT_SEC is scaled so that 360
C     days of length DELT_SEC are equal to 365 days of normal length (86400
C     seconds).

      CHSTPD = CHSTPD0
      IF (CHSTPD .LT. 1) CHSTPD = 1
      DELT_DAY = 1.0/FLOAT(CHSTPD)
      DELT_SEC = DAYL*(365./360.)/FLOAT(CHSTPD)


C     Initialize the number of dynamics timesteps per day (DYSTPD).  Adjust
C     DYSTPD to be at least 1.  Calculate the number of dynamical timesteps
C     per chemistry timestep (DYPRCH) and adjust it to be an integer value
C     of at least 1.  Recalculate DYSTPD to be consistent with CHSTPD and
C     DYPRCH.

      DYSTPD = DYSTPD0
      IF (DYSTPD .LT. 1) DYSTPD = 1
      DYPRCH = INT(DYSTPD/CHSTPD)
      DYPRCH = MAX0(1, DYPRCH)
      DYSTPD = DYPRCH*CHSTPD


C     Initialize the step-of-day to 1.

      CHSTEP = 1


C     Set the year number (IYEARC) to the initial year number, and initialize
C     the year-of-run counter (IYOR) to 1.

      IYEARC  = IYR0
      IYOR   = 1


C     Initialize the integer day-of-year (IDOY360) to IDOY0.  Make sure that
C     IDOY360 is at least 1, but less than 361.  Otherwise, adjust.
C     Note: Year numbers (IYEARC, IYOR) are NOT changed.

      IDOY360 = IDOY0

 100  IF (IDOY360 .LT. 1) THEN
         IDOY360 = IDOY360 + 360
         GO TO 100
         END IF

 200  IF (IDOY360 .GE. 361) THEN
         IDOY360 = IDOY360 - 360
         GO TO 100
         END IF


C     Initialize the decimal day-of-year (DOY360) and day-of-run counters
C     (DOR360) for the 360-day idealized calendar.

C      DOY360  = FLOAT(IDOY0) + 0.5
      DOY360  = FLOAT(IDOY0)
      DOR360  = 1.
      IDOR360 = INT(DOR360)


C     The month-of-year for the 360-day idealized calendar (MOY360/IMOY360)
C     is calculated from the decimal day number.  The formula used returns
C     MOY360 = 1 for DOY360 = 1.0, and MOY360 = 12 for DOY360 = 360.99...
C     (within machine precision).  Day-of-the-month (IDOM360) is calculated
C     from IDOY360 and MOY360.

      MOY360  = (DOY360 - 1.0)/30. + 1.0
      IMOY360 = INT(MOY360)
      IDOM360 = IDOY360 - 30*(IMOY360 - 1)


C     The day-of-year (DOY365/IDOY365) and effective day-of-run (DOR365) for
C     the 365-day real calendar are calculated from the corresponding
C     quantities in the 360-day idealized calender.  The formula used
C     correctly returns a value for DOY365 that is greater than or equal to
C     1.0 (for DOY360 = 1.0), and less than 366.0 (for DOY360 = 360.99...).

      DOY365  = 1.0 + (365./360.)*(DOY360 - 1.0)
      IDOY365 = INT(DOY365)
      DOR365  = 1.0 + (365./360.)*(DOR360 - 1.0)
      IDOR365 = INT(DOR365)


C     Determine the month-of-year for the 365-day real calendar (MOY365/
C     IMOY365) using the first-day-of-month lookup table (MNTHDOY).

      MDEX = 1
 300  IF (IDOY365 .GE. MNTHDOY(MDEX)) THEN
         MDEX = MDEX + 1
         GOTO 300
         END IF
      IMOY365 = MDEX - 1
      MOY365 = FLOAT(IMOY365) +
     2         SNGL(DOY365 - MNTHDOY(IMOY365))/
     3         FLOAT(MNTHDOY(IMOY365 + 1) - MNTHDOY(IMOY365))


C     Calculate the day-of-month for the 365-day real calendar using
C     IDOY365 and the first-day-of-month lookup table (MNTHDOY).

      IDOM365 = IDOY365 - (MNTHDOY(IMOY365) - 1)


C     Initialize the static day number if this is a static model run.
C     If IDOYST is 1 or greater, model will stop changing the day number
C     when IDOYST is reached, allowing for permanent equinox, etc., model
C     runs.  If IDOYST is initially greater than 360, it is adjusted to be
C     between 1 and 360.

      IDOYST = IDOYST0
 400  IF (IDOYST .GE. 361) THEN
         IDOYST = IDOYST - 360
         GO TO 400
         END IF


C     CALL SUBROUTINE CALC_SOLAR to initialize the solar declination (SOLDEC,
C     SINDEC, COSDEC) and Sun-Earth distance (ORBRAD).

      CALL CALC_SOLAR

      RETURN
      END




C     SUBROUTINE TIME_TICK:


      SUBROUTINE TIME_TICK

      include "timer.inc"

      INTEGER MNTHDOY(13)
      DATA    MNTHDOY/1, 32, 60, 91, 121, 152, 182,
     2                213, 244, 274, 305, 335, 366/

      REAL    DEGTORAD
      DATA    DEGTORAD /0.017453293/     !     PI/180.0


C     The step-of-day is incremented by 1.  If it exceeds the number of
C     timesteps per day, then it must be a new day and CHSTEP is reset to 1.

      CHSTEP = CHSTEP + 1
      IF (CHSTEP .GT. CHSTPD) CHSTEP = 1


C     The day-of-year (DOY360/IDOY360) and day-of-run (DOR360) for the
C     360-day idealized calendar are incremented by the day-fraction
C     timestep (DELT_DAY).

      DOY360  = DOY360 + DELT_DAY
      IDOY360 = INT(DOY360)
      DOR360  = DOR360 + DELT_DAY
      IDOR360 = INT(DOR360)


C     If this is a static run (IDOYST >= 1), check to see if the day number
C     has reached the static day number (IDOYST).  If so, adjust the day
C     number.

      IF ((IDOYST .GE. 1) .AND. (IDOY360 .GT. IDOYST)) THEN
         IDOY360 = IDOYST
         DOY360  = FLOAT(IDOY360)
         ENDIF


C     If IDOY360 indicates that the change to a new year has occured
C     (IDOY360 greater than or equal to 361), DAY360 is decremented by
C     360.0 so that it will be at least 1.0, but less than 361.0, and
C     IDOY360 is recalculated.  The year number (IYEARC) and year-of-run
C     counter (IYOR) are incremented by 1.

      IF (IDOY360 .GE. 361) THEN
         DOY360  = DOY360 - 360.0
         IDOY360 = INT(DOY360)
         IYEARC   = IYEARC + 1
         IYOR  = IYOR + 1
         END IF


C     The month-of-year for the 360-day idealized calendar (MOY360/IMOY360)
C     is calculated from the decimal day number.  The formula used returns
C     MOY360 = 1 for DOY360 = 1.0, and MOY360 = 12 for DOY360 = 360.99...
C     (within machine precision).  Day-of-the-month (IDOM360) is calculated
C     from IDOY360 and MOY360.

      MOY360  = (DOY360 - 1.0)/30. + 1.0
      IMOY360 = INT(MOY360)
      IDOM360 = IDOY360 - 30*(IMOY360 - 1)


C     The day-of-year (DOY365/IDOY365) and effective day-of-run (DOR365) for
C     the 365-day real calendar are calculated from the corresponding
C     quantities in the 360-day idealized calender.  The formula used
C     correctly returns a value for DOY365 that is greater than or equal to
C     1.0 (for DOY360 = 1.0), and less than 366.0 (for DOY360 = 360.99...).

      DOY365  = 1.0 + (365./360.)*(DOY360 - 1.0)
      IDOY365 = INT(DOY365)
      DOR365  = 1.0 + (365./360.)*(DOR360 - 1.0)
      IDOR365 = INT(DOR365)


C     Determine the month-of-year for the 365-day real calendar (MOY365/
C     IMOY365) using the first-day-of-month lookup table (MNTHDOY).

      MDEX = 1
 100  IF (IDOY365 .GE. MNTHDOY(MDEX)) THEN
         MDEX = MDEX + 1
         GOTO 100
         END IF
      IMOY365 = MDEX - 1
      MOY365 = FLOAT(IMOY365) +
     2         SNGL(DOY365 - MNTHDOY(IMOY365))/
     3         FLOAT(MNTHDOY(IMOY365 + 1) - MNTHDOY(IMOY365))


C     Calculate the day-of-month for the 365-day real calendar using
C     IDOY365 and the first-day-of-month lookup table (MNTHDOY).

      IDOM365 = IDOY365 - (MNTHDOY(IMOY365) - 1)


C     CALL FUNCTION CALC_DECL to update the solar declination (SOLDEC).

      CALL CALC_SOLAR

      RETURN
      END




C     SUBROUTINE CALC_SOLAR:


      SUBROUTINE CALC_SOLAR

      include "timer.inc"

      REAL    DAY, DECL


C     Declination constants from 1995 Astronomical Almanac, p.C24.
C     These are technically appropriate only for 1995, but variations are
C     small on the scale of decades.

      REAL    CDEC1,  CDEC2,  CDEC3,  CDEC4,  CDEC5,  CDEC6,  CDEC7
      REAL    CDEC1A, CDEC2A, CDEC3A, CDEC4A, CDEC5A, CDEC6A, CDEC7A
      REAL    RAU0,  RAU1,  RAU2

      DATA    CDEC1 /280.466/, CDEC2 /0.9856474/,
     2        CDEC3 /356.967/, CDEC4 /0.9856003/,
     3        CDEC5 /  1.916/, CDEC6 /0.020/,
     4        CDEC7 / 23.440/

      DATA    RAU0 /1.00014/, RAU1 /-0.01671/, RAU2 /-0.00014/


C     Constant DEGTORAD converts from degrees to radians.

      REAL    DEGTORAD
      DATA    DEGTORAD /0.017453293/     !     PI/180.0


C     FUNCINIT used to initialize constants on first call ONLY.

      LOGICAL FUNCINIT
      DATA    FUNCINIT /.FALSE./


C     MEANL is the mean solar longitude (abberation-corrected), MEANA is
C     the mean anomaly, ECLIPL is the ecliptic longtiude, and OBLIQU is the
C     approximate obliquity of the ecliptic.

      REAL    MEANL, MEANA, ECLIPL, OBLIQU, SINEL


C     If this is the first call to this function, initialize the DECxA
C     varibles, which are the declination formula constants converted from
C     degrees to radians.  Since we compile with the -static switch, these
C     values are retained from one call to the next.

ccelf      IF (.NOT. FUNCINIT) THEN
         FUNCINIT = .TRUE.
         CDEC1A = CDEC1*DEGTORAD
         CDEC2A = CDEC2*DEGTORAD
         CDEC3A = CDEC3*DEGTORAD
         CDEC4A = CDEC4*DEGTORAD
         CDEC5A = CDEC5*DEGTORAD
         CDEC6A = CDEC6*DEGTORAD
         CDEC7A = CDEC7*DEGTORAD
ccelf         ENDIF


C     Formulae for the mean solar longitude (MEANL), mean anomaly (MEANA),
C     ecliptic longitude (ECLIPL), and obliquity of ecliptic (OBLIQU) are from
C     the 1995 Astronomical Almanac, p. C24.  Simplification to the form for
C     the year 2000 introduces negligible error.

C     Formulae necessarily use DOY365, the day-of-year for the 365-day
C     real calendar.

      MEANL  = CDEC1A + CDEC2A*DOY365
      MEANA  = CDEC3A + CDEC4A*DOY365
      ECLIPL = MEANL + CDEC5A*SIN(MEANA) + CDEC6A*SIN(2.0*MEANA)
      OBLIQU = CDEC7A
      SOLDEC = ASIN(SIN(OBLIQU)*SIN(ECLIPL))/DEGTORAD
      SINDEC = SIN(SOLDEC*DEGTORAD)
      COSDEC = COS(SOLDEC*DEGTORAD)

      ORBRAD = RAU0 + RAU1*COS(MEANA) + RAU2*COS(2.0*MEANA)

      RETURN
      END
