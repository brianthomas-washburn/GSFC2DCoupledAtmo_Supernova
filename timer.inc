C-*-Fortran-*-  -  this is used for emacs/Fortran mode capability, or just do M-x fortran-mode
C
C     >>>>>>>>>>>>>>>>>>>  timer.inc  <<<<<<<<<<<<<<<<<<<


      INTEGER IYEARC, IYOR
      INTEGER IDOY360, IDOY365, IMOY360, IMOY365, IDOM360, IDOM365
      INTEGER IDOR360, IDOR365, IDOYST
      REAL    DOY360,  DOY365,  MOY360,  MOY365,  DOR360,  DOR365
      REAL    DELT_SEC, DELT_DAY
      INTEGER CHSTPD, CHSTEP, DYSTPD, DYPRCH
      REAL    SOLDEC, SINDEC, COSDEC, ORBRAD

      COMMON  /YEARS/ IYEARC, IYOR
      COMMON  /YR360/ DOY360, IDOY360, DOR360, IDOR360, MOY360, IMOY360,
     2                IDOM360, IDOYST
      COMMON  /YR365/ DOY365, IDOY365, DOR365, IDOR365, MOY365, IMOY365,
     2                IDOM365
      COMMON  /TIMER/ DELT_SEC, DELT_DAY, CHSTPD, CHSTEP, DYSTPD, DYPRCH
      COMMON  /SOLAR/ SOLDEC, SINDEC, COSDEC, ORBRAD


C     COMMON /YEARS/:

C     IYEARC:     Year number.  Initialized to the value of IYR0, a parameter
C                 passed to SUBROUTINE TIME_INIT.
C     IYOR:       Current year-of-run.  Keeps track of the number of years for
C                 which the model has executed.  Initialized to 1.  Note that
C                 IYOR is incremented on January 1 of each year, so if a
C                 simulation was started sometime during the year, less than
C                 a year will have passed when IYOR is incremented to 2.


C     COMMON /YR360/:

C     DOY360:     Day-of-year (decimal) for the 360-day idealized calendar.
C                 Initialized by the value IDOY0 passed to SUBROUTINE
C                 TIME_INIT and reset to 1.0 at the start of each new year.
C     IDOY360:    Day-of-year (integer) for the 360-day idealized calendar.
C     DOR360:     Day-of-run (decimal) for the 360-day idealized calendar.
C                 Keeps track of the absolute day number of the model run.
C                 Initialized to 1.0 but not reset at the start of a new
C                 model year.
C     IDOR360:    Day-of-run (integer) for the 360-day idealized calendar.
C     MOY360:     Month-of-year (decimal) for the 360-day idealized calendar.
C                 Initialized to 1.0 and reset to 1.0 at the start of each
C                 new year.
C     IMOY360:    Month-of-year (integer) for the 360-day idealized calendar.
C     IDOM360:    Day-of-month (integer) for the 360-day idealized calendar.
C     IDOYST:     Day-of-year for a static model run.  IDOYST less than 1
C                 indicates a normal, time-dependent run.  If IDOYST is
C                 greater than 1, it indicates the day-of-year for the 360-day
C                 idealized calendar at which the model becomes static.
C                 DOY360, IDOY360, MOY360, IMOY360, IDOM360, DOY365, IDOY365,
C                 MOY365, IMOY365, IDOM365, IYEARC, and IYOR stop increasing.
C                 DOR360 and DOR365 are unaffected.  Some useful values for
C                 IDOYST include:
C                      1:  First Day of the Year.
C                     79:  NH Spring Equinox / SH Fall Equinox (March 21).
C                    169:  NH Summer Equinox / SH Winter Equinox (June 21).
C                    260:  NH Fall Equinox / SH Spring Equniox (Sept. 21).
C                    350:  NH Winter Equinox / SH Summer Equinox (Dec. 21).
C                 Note that the model will run normally up until IDOYST is
C                 reached.  If it is desired that the model start from
C                 IDOYST, IDOY0 should be the same as IDOYST0 in CALL to
C                 SUBROUTINE TIME_INIT.


C     COMMON /YR365/:

C     DOY365:     Day-of-year (decimal) for the 365-day real calendar.
C     IDOY365:    Day-of-year (integer) for the 365-day real calendar.
C     DOR365:     Day-of-run (decimal) for the 365-day real calendar.
C     IDOR365:    Day-of-run (integer) for the 365-day real calendar.
C     MOY365:     Month-of-year (decimal) for the 365-day calendar.
C     IMOY365:    Month-of-year (integer) for the 360-day idealized calendar.
C     IDOM365:    Day-of-month (integer) for the 360-day idealized calendar.

C     Note: Since this module counts days (or day fractions) in terms of
C     the 360-day idealized calendar, the 365-day real calendar values are
C     actually effective values calculated from the 360-day calendar values,
C     using a linear scaling between the two years.  Since the 360-day year
C     consists of 12 30-day months, and the lengths of months in the real
C     calendar vary slightly, the month-of-year and day-of-month values
C     on the two calenders may not always agree.


C     COMMON /TIMER/:

C     DELT_SEC:   The length (in seconds) of the fundamental model timestep.
C                 Determined by the length of a day and the number of
C                 timesteps per day (specified by parameter CHSPD0 passed to
C                 SUBROUTINE TIME_INIT), adjusted so that 360 model days
C                 (with CHSPD0 timesteps per day) have the correct number of
C                 seconds for a real 365-day year.
C     DELT_DAY:   The length (in days, 360-day calendar) of the fundamental
C                 model timestep.
C     CHSTPD:     The number of chemistry timesteps per day for the 360-day
C                 idealized calendar.
C     CHSTEP:     Current chemistry step-of-day.  Keeps track of the step
C                 during the day when CHSPD is greater than 1.
C     DYSTPD:     The number of dynamics timesteps per day for the 360-day
C                 idealized calendar.  Adjusted so that there will be an
C                 integer number of dynamics timesteps in each chemistry
C                 timestep.
C     DYPRCH:     The number of dynamics timesteps per chemistry timesteps.
C                 Adjusted to be an integer number.


C     COMMON /SOLAR/:

C     SOLDEC:     Solar declination angle (determined from date).
C     SINDEC:     Sine of solar declination angle.
C     COSDEC:     Cosine of solar declination angle.
C     ORBRAD:     Current Sun-Earth distance (includes orbital ellipticity),
C                 in Astronomical Units.
