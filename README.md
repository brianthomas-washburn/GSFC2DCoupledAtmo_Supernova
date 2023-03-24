GSFC2DCoupledAtmo_Supernova

This file includes some basic information on the NASA Goddard 
Space Flight Center (GSFC) 2D (latitude/altitude) atmospheric chemistry/dynamnics model.
This refers to the newer verison of the model that includes 
temperature feedback and has user-definable latitude and altitude grids.

This version is further modified by Brian Thomas and Alexander Yelland,
primarily for the purpose of studying effect of nearby supernovae.

DISCLAIMER: This version of the model is NOT endorsed or supported by the orginal
NASA GSFC development team.  All modifications, and any errors therein, are solely
the responsibility of the authors noted here.

This README file written by Brian Thomas and Alexander Yelland, December 2020.
Revised for github repository, February 2023.

This software comes with NO WARRANTY or guarantee of fitness for any purpose.
The model has been extensively tested and used for certain scientific 
simulation purposes, but you must validate your implementation of the model.

Shared under Creative Commons license:
Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
   https://creativecommons.org/licenses/by-nc-sa/4.0/

Questions about this README can be sent to Brian Thomas, brian.thomas@washburn.edu; 
however BT is NOT a primary developer of the model and may not be able to answer 
all questions about the model!

This README contains basic information on using the model. For more extensive
scientific details see these papers:
* A model study of the impact of source gas changes on the stratosphere for 1850–2100, 
E. L. Fleming, C. H. Jackman, R. S. Stolarski, and A. R. Douglass, Atmos. Chem. Phys., 
11, 8515–8541, 2011, doi:10.5194/acp-11-8515-2011
See in particular the appendices.
* Jackman, C. H., and E. L. Fleming (2014), Stratospheric ozone response to a solar 
irradiance reduction in a quadrupled CO2 environment, Earth’s Future, 2, 331–340, 
doi:10.1002/2014EF000244.
* Fleming, E. L., et al. (2015), The impact of current CH4 and N2O atmospheric loss 
process uncertainties on calculated ozone abundances and trends,
J. Geophys. Res. Atmos., 120, 5267–5293, doi:10.1002/2014JD022067.


A list of the chemical constituents included in the model can be found
in a file named "2Dmodel_ConstituentList_constit15_9ua.pdf"

The model grid (altitude and latitude) can be adjusted by the user as needed.
It is currently configured with: 
* 76 altitude bins running from the ground to approximately 92 km, on pressure levels,
  with a 1 km resolution at the lower 60 bins and 2 km resolution above that
* 45 latitude bins running from 88deg South to 88deg North

Input files such as ionization profiles can be configured on a different grid
and then interpolated to the model grid after read in.

A common usage of this model by Brian Thomas and others has been to simulate
effects of increased atmospheric ionization due to various causes (solar energetic
particles, supernova photons and cosmic rays, etc.).
Ionization rates can be input in units of ions per square centimeter per second.
A conversion is done to volume rate (ions per cubic centimeter per second),
using the altitude bin sizes in the model.
The values are then used as a source of NOy with a conversion factor of 1.25.

We have also recently implemented an enhancement of the model lightning rate
based on the enhancement of ionization due to supernova cosmic rays
compared to the background galactic cosmic ray ionization.
It is important to note that, while there is good evidence that lightning
is connected to atmospheric ionization, the connection is not quantitatively
established.  Here, we assume a simplistic 1-to-1 correlation in enhancement.

The code package includes scripts written in NCL that post-process
some of the model output into netCDF files for easier storage and analysis.
To run these scripts you must install NCL:
https://www.ncl.ucar.edu/

System requirements to use the model:
* The gfortran compiler
* For post-processing using the NCL scripts, a working version of NCL

Basic model structure:
* Code files written in fortran (ending in ".f" or ".f90").
  This model has been developed over several decades and so there 
  are sections in older style fortran and some areas in newer style.
* A "main" file that setups up the model, calls subroutines, 
  runs over a time loop, and outputs results.
* Various subroutines for read-ins, interpolation, doing physics/chemistry, etc.
* A number of data files including boundary conditions for some constituents and
  climatologies for certain parameters that are not actively simulated.
  Data files generally end in ".xdr" or ".dat"
* Header or common files with declarations of model variables.
  These files end in ".h" or ".INC"

*** NOTE for github users:
    For tracking on github, xdr and dat input files had to be compressed.
    Therefore, you must gunzip those files before attempting to run the model.
    Runscript has a line to do this in case you forget.

The model is setup and run using a run script "Run_2D_Coupled.sh"
You must modify this script to set locations of some input data,
as well as where output data should be put.

Basic run procedure:
* Edit the Run script:
 - Define locations of input data and where output data goes
 - The Run script runs the model in 2-year intervals.  You can define how many intervals
  will be run.  Each interval produces a "restart" file which is used as the initial
  condition set for the next interval.  
 - The model can be run with "paleo" (pre-industrial) boundary conditions.
  If you wish to do a run using this, you can use initial conditions 
  produced by a "paleobase" run, which was a long spin-up run.
 - You can choose which year to start.  These years refer to the climatology used
  for certain forcings.  The earliest year is 1935.  For a run using the "paleobase"
  initial conditions use 1975.
 - Pre-industrial runs still use the 20th century forcings, 
  so you should not treat pre-industrial runs as simulating any particular time period.  
  The intent is to provide a background free of anthropogenic constituents, 
  espeically ozone depleting compounds.
  Pre-industrial runs are best used comparatively - a perturbed run with, e.g. ionization 
  input, compared against a control run.

The run script also contains flags that should be set to turn on/off 
run options such as ionization input and lightning enhancement implementation. In the same 
run script, you must specify what ionization profiles are being read-in from the 
input data location mentioned above. From there, there are two main logical variables 
that act as flags to determine whether you would like to do a basic run (without 
ionization or lightning enhancement factored in), a run with ionization, or a run with both 
ionization and lightning enhancement. NOTE: Since the lightning enhancement rate is based on
a ratio dependent on the galactic cosmic ray background (GCR) and the ionization profile
being read-in, it is not possible to do a run involving lightning enhancement without ionization.


* Run the Run script from the command line.  It is recommended to use this form, which 
will allow the process to run even if the user logs out, and redirects diagnostic
output to a log file:

> nohup ./Run_2D_Coupled.sh > RunLog.log &

* The Run script automatically places output files in the designated location,
runs the NCL scripts that convert model output into netCDF files with defined metadata,
and finally packages and compresses the raw output for storage if needed later.

You can use whatever software you prefer to access the netCDF files for analysis.
Included here (see analysisdir) are some NCL code files that can be used for analysis
and visualization of model output.

This model is large, complex, and has been developed over decades.  It is NOT easy to use,
so take care and validate as much as possible, especially if you make any changes.

Good luck and happy computing!

