#!/bin/bash
#
# This script sets up parameters and input files, compiles code,
#   runs, and deals with output of the GSFC coupled 2D atmospheric chem-climate model
# 
# Parameters below can be set to do multiple 2-year runs in order to do long runs (e.g. 20 yrs)
# Output is moved to a storage directory and processing code is called to convert to netCDF files
#
#  BASE9HB = 2012 base model (run12 - January 2012) - 2011/2012 updates merged with MERRA IAV
#

echo 'Start Time: '
date
echo 'Start Time:' > start.txt
date > startTime.txt

#------------------------------ EDIT THESE VALUES ------------------------#

case="HinputRun_12Jan2023"
#"ions_400kJpm2-6months_AllAbove70km_23Nov22"
#"SNCR20pc_100yr_70kmUpOnly_8Nov22"
#"SNCR100pc_3000yr_lightning_20Sep2021"

#ensure that this putdir actually exists (code checks below)
putdir="/home/bthomas/storage/atmoRuns/${case}"

#analysisdir is where NCL code for processing output lives:
analysisdir="/home/bthomas/storage/atmoRuns/analysis"

#number of 2-year long runs to do (so, 10 here is a 20 year run)
rnum=5

#year to start at:
startYr=1975 #1935  #2010
#For a run starting after the 'paleobaseRun', want to start at 1975; or can set whatever year for steady state

#which control file to use?
controlFile=control_coupled_steady-state.dat #control file for steady-state run
#control_coupled_steady-state-1935.dat  #this control file for paleobase run starting in 1935 with stead-state boundary conditions
#control_ne_coupled.dat
#control_coupled-1935start.dat #this control file for paleobase run starting in 1935

#which initial conditions file (f8->f4) to use? 
# (the year here, eg. 1974, is the ending year of the previous run)
f8file=f8-1974_paleobaseRun_6July2021.dat #this is the version to use for a "paleo" run (without anthropogenic influence) starting in 1975
#f4_ICs-paleobase.dat #this one has CFCs and other stuff 0'ed to make "clean" conditions - used to start paleobase run
#f8_2010_9hb.dat - modern day conditions

#read-in ionization data? (0=no, 1=yes)
ionsflag=0
  
#ionization data file due to supernova cosmic rays
# Note: This does not have to be ionization from supernova cosmic rays,
#       could be any ionization source you like.
#       We've just named it this since we are currently using
#       this code in that context.
SNCRionization=./ion_input/ions_400kJpm2-6months_AllAbove70km.dat
#SNCR20pc-ions_cm2-100.dat

#how long to wait (in days) until the ionization data is read in? (1 year = 360 days)
#ionStart must be less than 720!!!
ionStart=360

#lightning? (0=no, 1=yes) >>> (changes calculation within solv_9hb.f)
lightnflag=0
#lightnStart=ionStart
if [ $ionsflag -eq 0 ] ; then
    lightnflag=0
    #lighting can only included with an ionSource; see ionRatio_calc.f & solv_9hb.f
fi


#additional H input? (0=no, 1=yes)
Hinflag=1
#how long to wait (in days) until H input starts? (1 year = 360 days)
#must be less than 720!!!
HinStart=180
#when to stop?
HinStop=540
#must be less than 720!!!  

#paleo BCs (which tempin subroutine and BCs file to use)? 0=no, 1=yes
paleo=1

#--- use these flags for debugging: ---#

#compile model? 0=no 1=yes
compile=1
#run the model? 0=no 1=yes 
runmodel=1
#store data?  0=no 1=yes
store=1
#compress data?  0=no 1=yes
compress=1

#--- end debug flags section ---#

#------------------------------ end EDIT section -------------------------#

#Check that the putdir directory exists:
if [ ! -d ${putdir} ]; then
    echo "  ERROR: " ${putdir} " directory does not exist"
    echo "   Aborting "
    exit
fi

#Check that ion input file exists:
if [ $ionsflag -eq 1 ] ; then
    if [ ! -f ${SNCRionization} ]; then
	echo "  ERROR: " ${SNCRionization} "  does not exist"
	echo "   Aborting "
	exit
    fi
fi

cp ${f8file} fort.4
# the fort.4 file is a 'restart' or initial conditions file;
# the f8 file is the initial conditions file
# in this code version, an 'f8' file is output each year (fort.year),
# which can be used as a fort.4 for restart

cp ${controlFile} fort.7
# control input code reads in this startYear.txt file separately
echo ${startYr} > startYear.txt

#F68H2OMR.DAT is read in from unit 12 in rdata.f - no longer used w/ compute H2O
#F40GCRMM.DAT is read in from unit 96 in gcrs.f
#aerosol.dat is read in from unit 21 in aerosols_xe.f
#stab_2d69.xdr is read in from unit 41 in photin_xb.f
#indx_2d69.dat is read in from unit 42 in photin_xb.f
#xsecttd10_jpl10.dat is read in from unit 23 in crossin_9ha_jpl10.f
#latheat_base8d.xdr is read in from unit 78 in tempin_9hb.f
#trop2d_base8e.xdr is read in from unit 74 in tempin_9hb.f 
#troph2oday_2d.dat is read in from unit 60 in tempin_9hb.f 
#gwbc_base9ba.xdr is read in from unit 63 in tempin_9hb.f 
#chibbc2_base8xf.xdr is read in from unit 75 in tempin_9hb.f 
#co2bc.xdr is read in from unit 79 in tempin_9hb.f
#h2obc.xdr is read in from unit 79 in tempin_9hb.f

cp bc_wmo10_9ha.dat          fort.20
cp reac_9ha_jpl10.dat        fort.29
cp reac_9fk_het_jpl06.dat    fort.31
cp reac_work_me.dat          fort.30
cp wflUmean.dat              fort.92
cp xsect10_jpl10.dat         fort.93
cp davg.dat                  fort.94
cp xsect_no_9dg.dat          fort.95

cp com2d_9hb_Hinput.h          com2d.h
cp comphot_9gl.h             comphot.h
cp com_aerg_9da.h            com_aerg.h
cp com_aers_9ha.h            com_aers.h
cp PARAM_9ha.INC             PARAM.INC


#  copy files for coupled model:

cp f13_start.dat             fort.13
#dynamics files - these are all formatted ASCII files (binary won't work this way in LINUX)
#
cp ir2.dat                   fort.1
cp switches_9gl.dat          fort.9
cp co2fix.dat                fort.15
cp teheweak.dat              fort.516
cp latent_heat.dat           fort.40
cp kzz_param.dat             fort.44
cp K_table.dat               fort.45
cp TR_table.dat              fort.43
cp topav.dat                 fort.47
cp f60_switches.dat          fort.60
cp isw.dat                   fort.81
#cp idealheightsc.dat         fort.89
cp mrsheightsa.dat           fort.89
cp nmcta.dat                 fort.90
cp ozonea.dat                fort.91
cp n2oa.dat                  fort.85
cp water.dat                 fort.595
#######elfcp nmc_ubar.dat              fort.97
#cp mwbrk.dat                fort.98
cp mwbrk1.dat                fort.98
cp P_wave_number_9gl.dat     fort.99

cp rayleigh.dat              fort.101
cp co2var.dat                fort.113
cp gwray_ini.dat             fort.114


#introducing ionization data
echo  ${ionsflag} > ionsflag.txt
echo  ${ionStart} > ionStart.txt
cp    ${SNCRionization} SNCRionization.dat 

#Hinput settings
echo  ${Hinflag} > Hinflag.txt
echo  ${HinStart} > HinStart.txt
echo  ${HinStop} > HinStop.txt

#introducing lightning calculation
echo ${lightnflag} > lightnflag.txt

####################################################
# ------- Compile Fortran source code files -------#
####################################################

###ifort -mcmodel=small -fast \
#######ifort -f77rtl -fast -parallel -O3 -tpp7 -ip -static \
#######ifort -f77rtl -fast -O3 -tpp7 -ip -static \        - on ASSESS "-tpp7" does not work
#ifort -f77rtl -fast -O3 -ip -static \

if [ $compile -eq 1 ] ; then
    rm go.exmod

  if [ $paleo -eq 1 ] ; then
	  cp tempin_9hb_paleo.f       tempin_9hb_A.f
	  cp bc_wmo10_9ha-paleo.dat   fort.20
  fi
  if [ $paleo -eq 0 ] ; then
  	cp tempin_9hb_Orig.f   tempin_9hb_A.f
  fi

#linterp48_9aa.f -- 24March2021: this was in the list, but does not exist as a file?
  
  gfortran -fdollar-ok -fno-align-commons -std=legacy -cpp -fno-automatic -fmax-errors=1 -Wno-argument-mismatch \
    main_9hb.f\
    aerchem_9hb.f\
    aerosol_mc2_9hb_A.f\
    aerosols_9fe.f\
    bcin_9ha_A.f\
    bd22d_9ha.f\
    binterp_9aa.f\
    binterp8_9fe.f\
    boundc_9hb.f\
    chinit_aer_E.f\
    cjsor2_9ga.f\
    clock_aer_G_9fe.f\
    co2_chou_9ga.f\
    colden_9gl.f\
    coldiur_9gl.f\
    control_9hb.f\
    crosec_9ha.f\
    crossin_9ha_jpl10_A.f\
    csin_no_9fe.f\
    dailyint_9hb.f\
    dump_9hb.f\
    dynout_9gtp_A.f\
    epd10_9gn.f\
    fastdiur_9hb.f\
    fill_col.f\
    fill.f\
    filla_9hb.f\
    gcrs_9fe.f\
    GCRionsIn.f\
    getdyn_9hb_A.f\
    getkyz_9gl.f\
    getmomx_9gg.f\
    getprob_9hb_A.f\
    getrad_9ga.f\
    gwave_9fe.f\
    gwaveb_9ga.f\
    gwdrag_9ga.f\
    gwray_9ga.f\
    hetchem_9gl.f\
    hetmap_9fk.f\
    input_9hb.f\
    ionSourceIn.f\
    ionsRatio_calc.f\
    irrad_9gj.f\
    jacob_9ha.f\
    jgetub_9ga.f\
    jgetvw_9ga.f\
    jmap_9ha.f\
    jsor2_9gb.f\
    kmap_9ha.f\
    kwave_9fe.f\
    kwave_coup_9gl.f\
    kyybackot_9gl.f\
    kzzadj_9gl.f\
    life_9hb.f\
    linterp_9aa.f\
    linterp8_9aa.f\
    mult_9gl.f\
    natice_9gl.f\
    newdify_9fe.f\
    newdifyz_9fe.f\
    newdifz_9fe.f\
    no_photdis_9gl.f\
    nwt2dd_9gl.f\
    outputs_9ga.f\
    ox_9ga.f\
    pdfinsol_9hb.f\
    pdfint_9gl.f\
    photheat_9gl.f\
    photheatin_9gl.f\
    photin_9gl.f\
    prminmax_9ga.f\
    prmlib_9gl.f\
    psc2d_9ga.f\
    radiate_9hb.f\
    radlib_9gk_nosad.f\
    radlibz_9ga.f\
    radt_9fe.f\
    rain_9hb.f\
    rcread_9ha_A.f\
    rdaer_9ga.f\
    reac_9hb.f\
    readstr_9gto_A.f\
    rfluxinterp_9gl.f\
    setdaily_9hb.f\
    setupa_9hb_G.f\
    smthpsi_9ga.f\
    solflin_9fe.f\
    solv_9hb_Hinput.f\
    solv58_9ga.f\
    sorad_9hb.f\
    spdr_9hb.f\
    srband_9fe.f\
    streamf_9gtp_A.f\
    stridag_9ga.f\
    strmlib_9ga.f\
    tempin_9hb_A.f\
    templon_9gn.f\
    thermw_9ga.f\
    timeav_9af.f\
    timeavj_9hb.f\
    timer_9ga.f\
    tpcore_9fe.f\
    tpcore_dyn_9gl.f\
    trcrlib_9gl.f\
    tropkzz_9ga.f\
    twodlib_9gl_A.f\
    twods_9hb.f\
    ubardiag_9gt_A.f\
    undump_nb.f\
    unfilla_9hb.f\
    unjmap_9ha.f\
    vwmass_9fe.f\
    wed_rad_9ga.f\
    xcoupled_9hb.f\
    xcoupled_in_9ga.f\
    xgwdrg_9ga.f\
    xtrans_9gl.f -o go.exmod

  rm *.u
  rm *.o
  #./go.exmod

  echo " "
  echo "COMPILED."
  echo " "

fi #compile

  if [ $ionsflag -eq 1 ] ; then
      echo " "
      echo "SNCRionization= " ${SNCRionization}
      echo " "
  fi

  if [ $Hinflag -eq 1 ] ; then
      echo " "
      echo "H input run."
      echo "  starting at " $HinStart
      echo "  stopping at " $HinStop
      echo " "
  fi
  
for (( r=1; r <= rnum; r++)) #loop to do multiple sequential runs:
do
  
  if [ $r -gt 1 ] ; then
  	# Set the start year for the next run:
  	startYr=$[${startYr}+2] #2 for 2-year run
  	# Same control file will be used, except we need to update the start year,
  	#   control input code reads in this startYear.txt file separately
  	echo $startYr > startYear.txt
  fi #updating and startYr

  echo " "
  echo "run number: " ${r}
  echo "${r}" > runNum.txt
  echo "startYr: " ${startYr}

  if [ $r -gt 1 ] ; then 
    echo "fort.4 file: " fort.${f8yr} 
  fi

  if [ $runmodel -eq 1 ] ; then
  	echo "running model... "
  	#run the model:
  	./go.exmod > ${case}_r${r}_startYr${startYr}.log
  fi
    
  if [ $store -eq 1 ] ; then
	
  	# files fort.14 and fort.33 are original GSFC developers output files.
  	#   I (BThomas 2018) am keeping them here,
  	#   but have created my own that I will use for my analysis
  	mv fort.14 ${putdir}/f14_A_2yr_${case}_r${r}_startYr${startYr}.dat
  	# fort.14 is output file containing column O3 data
  	#   It is the developers output file for O3; 
  	mv fort.33 ${putdir}/f33_A_2yr_${case}_r${r}_startYr${startYr}.dat
  	# fort.33 has various outputs from 'main', including lat, alt, pressure arrays,
  	#   some dynamics stuff (velocities W and V),
  	#   and cnout, which is specified chemical constituents plus some other data
  	#   Basically, it is the developers main 'output for analysis' file.
  	#
  	# f8 (restart files):
  	cp fort.${f8yr} ${putdir}/
  	#
  	# Thomas output files:
  	mv BToutput_alt-lat_daily.out ${putdir}/specout_daily-2yr_${case}_r${r}_startYr${startYr}.out
  	mv BToutput_alts-O3-NO2_TUV_daily.out ${putdir}/alts-O3-NO2_TUV_daily-2yr_${case}_r${r}_startYr${startYr}.out
  	mv BToutput_coldenhno3_daily.out ${putdir}/coldenhno3_daily-2yr_${case}_r${r}_startYr${startYr}.out
  	mv BToutput_coldens_daily.out ${putdir}/colden_daily-2yr_${case}_r${r}_startYr${startYr}.out

  	# run makeNCfile code to convert the text output files to netCDF:
  	echo 'Running NCL code on output... '

    echo $putdir > nclInput_colden.txt
    	echo colden_daily-2yr_${case}_r${r}_startYr${startYr}.out >> nclInput_colden.txt
    	echo coldenhno3_daily-2yr_${case}_r${r}_startYr${startYr}.out >> nclInput_colden.txt
    	cp nclInput_colden.txt $analysisdir
  	ncl -Q ${analysisdir}/makeNCfile_coupled2D-output_coldens.ncl

  	echo $putdir > nclInput_alt-lat.txt
  	  echo specout_daily-2yr_${case}_r${r}_startYr${startYr}.out >> nclInput_alt-lat.txt
      cp nclInput_alt-lat.txt $analysisdir
    ncl -Q ${analysisdir}/makeNCfile_coupled2D-output_alt-lat.ncl

  	# also store the logfiles:
  	cp ${case}_r${r}_startYr${startYr}.log ${putdir}/
	
  fi #store data

  # need the initial conditions file for next run
  # the last fort.year file from this run is copied into f4 for next run:    
    f8yr=$[$startYr+1]
    cp fort.${f8yr} fort.4 
  # and copy out this last fort.year file to storage:
    # cp fort.${f8yr} ${putdir}/

done #loop over multiple runs


if [ $compress -eq 1 ] ; then

  echo " "
  echo 'Cleaning up output (zipping text files)...'
  echo " " 

  tar -czf ${putdir}/output_textfiles_colden-${case}.tar.gz ${putdir}/colden*.out
    rm ${putdir}/colden*.out
  tar -czf ${putdir}/output_textfiles_specout-${case}.tar.gz ${putdir}/specout*.out
    rm ${putdir}/specout*.out
  tar -czf ${putdir}/output_textfiles_f8-f14-f33_${case}.tar.gz ${putdir}/f*.dat ${putdir}/fort.*
    rm ${putdir}/f*.dat ${putdir}/fort.*
  tar -czf ${putdir}/alts-O3-NO2_TUV_daily-2yr_${case}.tar.gz ${putdir}/alts*start*.out
    rm ${putdir}/alts*start*.out
    
fi #compress data


echo 'End time:' > doneTime.txt
date | cat >> doneTime.txt
cat start.txt startTime.txt doneTime.txt > doneEmail.txt
mail -s Coupled2DModel_FinishedNotification brian.thomas@washburn.edu < doneEmail.txt

echo 'End time:'
date

rm BToutput_*.out

cp ${case}.log ${putdir}/

# end of script
