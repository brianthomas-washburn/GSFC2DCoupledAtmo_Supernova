;This NCL code is used to do analysis (mainly plotting) of "specout" files
; from the NASA GSFC coupled 2D atmo chem/dynaimcs model.
; By Brian Thomas, v1 2019, current version 2023

;You can specify a single file to look at the data,
; or specify two files to compare using eithr a ratio or % difference.
;
;There is an option to convert values to "mixing ratio",
; which is computed using O2 values read in from the file

;You can run over a range of input files from 2-year duration runs,
; to do that you need to uncomment the time loop -- just search for "time loop"
; to find where to comment/uncomment

begin

;--------------------------- Paths and Parameters --------------------------

  ;specify case data:
;  NOTE: case1 is used as the Control case for ratio or perdiff comparisons
  case1 = "controlRun_19July2021"
  ;"paleobaseRun_6July2021"
  indir1   = "path/"+case1+"/"

  ;specify another case data: (eg. to compare to a control run)
  case2 = "case"
  indir2   = "path/"+case2+"/"
  
  ;starting year for this set of runs:
  stYr = 1975 ;1935

  ;number of runs in the set of runs (2 years per run)
  itmax = 5 ;10 ;20 
  
;If NOT using the time loop, specify the specific run year here:
  ;Can specify it directly, or using the run number (it)
  it = 10
  yr = stYr+2*(it-1)
  
  ;What data to extract and to be plotted? (this goes in output filenames)
  ; (This variable is read from both files when doing a comparison)
  datPlot = "O3" ;"HOx" ;"NOy"  ;"HNO3_rainout"    "Temp"  

  ;Look at data, ratio, or percent difference?
  analyze = "perdiff" ;"data" ; "ratio" 

  pname = datPlot+analyze
  
  ;Convert concentrations to mixing ratios? [Only applies if looking at data, not ratio/perdiff]
  convMR = "no"

  ;any scale factor to apply? (only applies when looking at data, not ratio/perdiff)
  scale = 1/1.0E12
  datAunits = "O3 concentration 1E12 molecules/cm3" ;specify appropriate units here as needed ;"K" ;
  
  ;specify plot file name:
  plotFileName ="altTimePlot_"+case2+"_"+pname+"_startYr"
  ;"altTimePlot_O3ratio_startYr"
  ;"altTimePlot_O3perdiff_startYr"
  ;"altTimePlot_O3concentration_contour_startYr"
  ; "temperature" ;

  ;title for plot
  plotTitle = pname+" "+case2+" vs "+case1
  ;"NOy perdiff"+" "+case2+" vs "+case1
  ;"O3 ratio"+" "+case2+" vs "+case1
  ;"O3 concentration 1E12 molecules/cm3"
  ;"Temperature "+" "+case2+" vs "+case1 

  ; You can plot contours OR line profile (not both at the same time): 
  plotAltTimeCont = "yes" ;Plot contours?                 
  plotlineProfile = "no" ;Plot line profile? 

  ;for line profile plot only, choose time point to plot:
  tchoice = 90

  ;choose latitude to plot
  latchoice = 22
  ;to loop over latitudes set latstart and latstop to be different
  ;  and set maximum number of lats and step 
  ;to just do one lat set start/stop equal and latloopnum=1
  latstart = 1
  latstop = 22
  latloopnum = 23
  latstep = 2
  
  ;(may need to change other plot parameters below)

  
;---- end of parameter setting ----;
;
;..... shouldn't need to edit anything below except perhaps in plot settings...
;

;loop over latitudes:
  do il=1,latloopnum,latstep
    latchoice = il

;Now we loop over each 2-year run input file:
;time loop:
  do it=1,itmax    
    yr = stYr+2*(it-1)
    
    infile1 = "specout_daily-2yr_"+case1+"_r"+it+"_startYr"+yr+".out.nc"
  ;"specout_daily-2yr_atmoRun_control2020_r1_startYr1975.out.nc"

    infile2 = "specout_daily-2yr_"+case2+"_r"+it+"_startYr"+yr+".out.nc"
  ;"specout_daily-2yr_SNCR10pc_5Dec2020_r1_startYr1975.out.nc"

  ;Read file:
    f1 = addfile(indir1+infile1,"r")
 
  ; plot file name:
    plotFile = plotFileName+yr

;-------------------------- Get data and other info --------------------------

   ;get data from file 1:
    dat1 := f1->$datPlot$   
  
  ;Also get the latitude and altitude (level) arrays:
    lat = f1->lat4arr
    lev = f1->press
    
  ;Set the maximum values for the arrays:
    sizes = dimsizes(dat1)      
    tmax = sizes(0)
    latmax = sizes(1)
    levmax = sizes(2)  
  
    printVarSummary(dat1)
    printVarSummary(lat)
    printVarSummary(lev)

    print("------------------------------------"+yr) 
    print(dat1(0,0,18))
    
    plotFile := plotFile+"_lat"+lat(latchoice)+"deg"

;------------------------------ Analysis ------------------------------

    if(analyze.eq."data") then
      datA := dat1
      datA = datA*scale ;scale factor (for nicer values on plot, or whatever)
      datA@units=datAunits
      printVarSummary(datA)
    
    ;Convert the molecules/cm3 concentration to mixing ratio:
      if(convMR.eq."yes") then
        print("------------- dat1 max: "+max(dat1)+" min: "+min(dat1))
        datO2 = f1->O2
        print("------------- O2 max: "+max(datO2)+" min: "+min(datO2))
        ; printVarSummary(datO2)
        datO2 = datO2/0.21
        dat1mr = dat1/datO2 
        dat1mr = dat1mr*(1.e6)
        dat1mr@units="ppm"
        dat1mr!0 = "time"
        dat1mr!1 = "lat"
        dat1mr!2 = "lev"
        dat1mr&lat = lat
        dat1mr&lev = lev
        print("------------- mixing ratio in ppm")
        printVarSummary(dat1mr)
        print("------------- dat1 max mr: "+max(dat1mr)+" min mr: "+min(dat1mr))
        datA := dat1mr
        printVarSummary(datA)
      end if
    end if ;data


    if(analyze.eq."ratio") then
    ;Load data from control file:
      f2 = addfile(indir2+infile2,"r")
      dat2 = f2->$datPlot$
      printVarSummary(dat2)
    ;Compute ratio:
      ratio = new((/tmax,latmax,levmax/),float)
      ratio@_FillValue = -9999.99     ;change fill value to this
      ratio@units="Ratio: case/Control"
      ratio!0 = "time"
      ratio!1 = "lat"
      ratio!2 = "lev"
      ratio&lat = lat
      ratio&lev = lev
      ratio = dat2/dat1
      printVarSummary(ratio)
      printMinMax(ratio,False)
      datA := ratio
    end if ;ratio


    if(analyze.eq."perdiff") then
    ;Load data from control file:
      f2 = addfile(indir2+infile2,"r")
      dat2 = f2->$datPlot$
      printVarSummary(dat2)
    ;Compute percent difference
      perdiff = new((/tmax,latmax,levmax/),float)
      perdiff@_FillValue = -9999.99     ;change fill value to this                                 
      perdiff@units="Percent difference: (case-Control)/Control"
      perdiff!0 = "time"
      perdiff!1 = "lat"
      perdiff!2 = "lev"
      perdiff&lat = lat
      perdiff&lev = lev
      perdiff = 100*(dat2-dat1)/dat1
      printVarSummary(perdiff)
      printMinMax(perdiff,False)
      datA := perdiff
    end if ;perdiff

;------------------------------ Plotting ------------------------------

    wks = gsn_open_wks("png",plotFile)

    if(plotAltTimeCont.eq."yes") then
    
      datAs := datA(time|:,lev|:,lat|:)
      printVarSummary(datAs)
      res = True
      res@tiMainString = plotTitle
      res@cnFillOn = True
      res@cnLinesOn = False
    ;res@lbOrientation = "Vertical"
  
      colors = "gsltod"
      res@cnFillPalette = colors
    
    ;--- Here, I want a wider range in level (altitude) than the data covers
    ; in order for that to work in plotting, I have to extend the data range,
    ; so that it starts at 1000 mb instead of 850
     
    ;-Create the new level array with 1 extra element
      lev_new           = new(levmax+1,typeof(lev))
      lev_new(1:levmax) = lev
      lev_new(0)        = 1000 
       
    ;-Create the new data array extended into new lat/alt area (with *empty* values there)   
      datp = new((/tmax,levmax+1,latmax/),typeof(datAs))
      datp(:,1:levmax,:) = datAs
      datp!0 = "time"
      datp!1 = "lev"
      datp!2 = "lat"
      datp&lev = lev_new
      datp&lat = lat

      res@tmYLMode        = "Explicit"                     ; Force labels where we want them
      res@tmYLValues      = (/0.001,0.01,0.1,1,10,100,1000/)
      res@tmYLLabels      = "" + res@tmYLValues
      ;res@trYMinF = 100
      res@tiYAxisString   = "Pressure (hPa)"
      res@tiXAxisString   = "Time (days)"

    ;These can be specified, or made blank, or commented out to let plot routine fill them:     
      res@gsnLeftString = " "
      res@gsnRightString = " "
       
    ;don't draw yet because contour label bar will be modified below   
      res@gsnFrame             = False       ; don't draw yet
      res@gsnDraw              = False       ; don't advance frame yet

    ;Here you can specify specific scale bar values:     
    ;   res@cnLevelSelectionMode = "ManualLevels" ; manual levels
    ;   res@cnMinLevelValF       =  1.0E6 ;0.0 ;datmin
    ;   res@cnMaxLevelValF       =  9.01E6 ;datmax
    ;   res@cnLevelSpacingF      =  0.5E6

    ;Other scale bar settings:
    ;   res@lbOrientation        = "Vertical"     
      res@pmLabelBarOrthogonalPosF = 0.2
      res@vpHeightF = 0.5

    ;moves location of x axis title:
      res@tiXAxisOffsetYF = 0.25

    ;declare the plot:   
      datp_reorder = datp(lev|:,time|:,lat|:)
      plot = gsn_csm_pres_hgt(wks,datp_reorder(:,:,latchoice),res)


      ;add text to specify latitude:
      txres = True
      txres@txFontHeightF = 0.03
      dum = gsn_add_text(wks,plot,"lat "+lat(latchoice)+"deg",100,1.0,txres)

;    ; Retrieve contour levels, for use in formatting numbers on label bar
;      getvalues plot@contour
;      "cnLevels" : levels
;      end getvalues
;      res@lbLabelStrings = sprintf("%3.1f",levels)   ; Format the labels    
;
;    ;have to re-declare plot since changed a resource:
;      plot = gsn_csm_pres_hgt(wks,datp(tchoice,:,:),res)
      
    ;actually draw the plot:   
      draw(plot)
      frame(wks)
    end if ;plotAltTimeCont


    if(plotlineProfile.eq."yes") then
    
      datp := datA(time|:,lev|:,lat|:)
      printVarSummary(datp)
    
      res = True
      res@trYReverse = True ;reverse Y-axis
      res@xyYStyle = "Log"  ;make Y-axis logarithmic
      res@xyLineThicknessF = 4.0
      res@trYMinF = 1.E-3
      res@tmXBLabelFontHeightF = 0.015
      
      res@tiXAxisString = plotTitle
      
      plot = gsn_csm_xy(wks,datp(tchoice,:,latchoice),lev,res)
    end if ;plotlineProfile


  end do ;end of time loop
  end do ;end of lat loop
  end ;end of script

;This stuff is older, left here for reference (commented out)
  /;
    res = True
    res@cnFillOn = True
    res@cnLinesOn = False
    res@lbOrientation = "Vertical"
    
    res@gsnYAxisIrregular2Log = True
    res@xyYStyle = "Log"  ;make Y-axis logarithmic
          
    colors = "gsltod"
    res@cnFillPalette = colors

    ;    res@tiMainString = "yr "+yr+" day "+tday
              
    ;  plot = gsn_csm_contour(wks,datp,res)
    ;    printVarSummary(datp(tday,:,:))
    ;          print(datp(tday,:,:))
    print("make plot")
    plot = gsn_csm_pres_hgt(wks,datp(tday,:,:),res)
    print("plot made")
  ;/