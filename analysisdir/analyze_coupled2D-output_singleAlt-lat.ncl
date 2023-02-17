;This NCL code is used to do analysis (mainly plotting) of a single level (altitude) from
; "specout" files
; from the NASA GSFC coupled 2D atmo chem/dynaimcs model.

; based on analyze_coupled2D-output_coldens.ncl

;v1.0 July2021 by Brian Thomas

load "./gsnColorRange.ncl"

begin
 
;---- Variables for data to read in, what to plot, etc: ---------;
;  
;2 cases, for comparing, or just look at one
; If doing percent difference, then case2 is the "new", case1 is control
  case1 =  "case1"
  case2 =  "case2"
  ;"SNCR20pc_100yr_lightning_8Sep2021"
  ;"SNCR50pc_1000yr_26Aug2021"
  ;"SNCR50pc_1000yr_lightning_15Sep2021"
  
  in_dir1 = "path/"+case1+"/"
  in_dir2 = "path/"+case2+"/"
; 
;Analyze/Plot data, or percent difference between data sets?
  opt = "perdiff"
;if "data", which data file to look at?
  opt_df = 2
;
;What level (altitude) to look at? (0 = lowest level; main arrays are in pressure levels)
  levChoice = 0
; 
;How many 2-year runs?  
  itmax = 10 ;20
;Number of time points in each run (NOTE: below assumes daily data):
  rtmax = 2*360
;Set the time axis label:
  timeLabel = "Days"
;setup time variables array for plotting  
  tmax = itmax*rtmax
  times = ispan(0,tmax-1,1)
;  printVarSummary(times)
;  t1 = 0 ;start time for plot
;  t2 = tmax ;end time for plot
;Year of start of run sequence:
  stYr = 1975 
;  
;How many latitudes?
  latmax = 45
;  
;What data to be plotted? (this goes in output filenames)
; (This variable is read from both files)
  datPlot =  "O3" ;"Temp" ; "NO2"; ;"NOy" ;"HNO3_rainout"
;

;other scale factor? Set to 0 for no scale factor to be used
  scale = 0 ;1.02962e-15 ;For HNO3 rainout, converts to mg of N per m^2 per s
  scaleUnits = "mg of N per m^2 per s"
;
;Average over latitude? (global average) This also turns on line plot for the result.  
  latav = "no"
;
;For contour plot, use red-blue color scheme?
;  This should be used for O3 perdiff lat v time plots
  plot_redBlue = "yes" ;"yes" 
  colors = "BlRe" ;color table to use
  
;------------- END of setup variables -----------------;
  
;Arrays to store all of the data (from multiple run files):
  data1 = new((/itmax*rtmax,latmax/),float)
  data2 = new((/itmax*rtmax,latmax/),float)
;Loop over runs:  
  dst = 0 ;array time position to start in, for filling total data array
  dsp = 0 ; (stop position)
  do it=1,itmax
;    print(it)
;File 1:
    in_file1a = "specout_daily-2yr_"+case1
    in_file1b = "_r"+it
    yr = stYr+it*2-2
    in_file1c = "_startYr"+yr
    in_file1 = in_file1a+in_file1b+in_file1c+".out.nc"
    f1 = addfile(in_dir1+in_file1,"r")
;    
    dat1 := f1->$datPlot$
;
;Get the pressure levels and altitudes:
;    pressureArray = f1->press
;    altArray = f1->zalt90arr
    levChoicePress = f1->press(levChoice) ;pressureArray(levChoice)
    print(levChoicePress)
    levChoiceAlt = f1->zalt90arr(levChoice) ;altArray(levChoice)
    print(levChoiceAlt)
;    
;File 2:    
    in_file2a = "specout_daily-2yr_"+case2
    in_file2b = "_r"+it                                  
    yr = stYr+it*2-2                         
    in_file2c = "_startYr"+yr          
    in_file2 = in_file2a+in_file2b+in_file2c+".out.nc"
    f2 = addfile(in_dir2+in_file2,"r")  
;
    dat2 := f2->$datPlot$
;
;    printVarSummary(dat1)
;    printVarSummary(dat2)
;
    dsp = (dst+rtmax)-1
;    print(dst+"  "+dsp)
    data1(dst:dsp,:) = dat1(:,:,levChoice)
    data2(dst:dsp,:) = dat2(:,:,levChoice)
;
;    printVarSummary(data1)
;    printVarSummary(data2)
;
    dst = dsp+1
;    
  end do ;it (loop on multiple runs)   
;
;Read in the latitude values (all files are the same so it doesn't matter which one here):    
  lat = f1->lat4arr
;  print(lat)
    
;Plot data (just one file, not comparing the two):
  if(opt.eq."data") then
    if(opt_df.eq.1) then
      datp = data1
      cname = case1
    end if
    if(opt_df.eq.2) then
      datp = data2
      cname = case2
    end if
;other scale factor?
    if (scale.ne.0) then
      datp = datp*scale
      datp@units = scaleUnits 
    end if;1.0
    ;printVarSummary(datp)
  end if ;plot data or compare

;Compare the two data sets (percent difference):
  if(opt.eq."perdiff") then
    cname = case2
    datp = (data2-data1)/data1
    datp = 100*datp
    copy_VarMeta(data2,datp)
    datp@units = "percent difference"
    printVarSummary(datp)
  end if ;compare (perdiff)    

;Latitude average (global average)
; (see https://www.ncl.ucar.edu/Document/Functions/Built-in/wgt_areaave.shtml):
  if(latav.eq."yes") then
    rad = 4.0*atan(1.0)/180.0
    clat = cos(lat*rad)
    datp_lAvg = dim_avg_wgt_n_Wrap(datp,clat,1,1)
    datp_lAvg@_FillValue = -9999.99
    printVarSummary(datp_lAvg)
;    print(datp_lAvg)
;     
;Plotting:
    if(opt.eq."perdiff") then
      axisLabel="% change "+datPlot+" (glob avg)"+datp@units+" at "+levChoiceAlt+" km alt"
    end if
    if(opt.eq."data") then
      axisLabel=datPlot+" (glob avg) "+datp@units+" at "+levChoiceAlt+" km alt"
    end if
    plotFile := datPlot+"_"+opt+"_globalAvg-lev"+levChoice+"-"+cname
    wks   := gsn_open_wks ("png",plotFile )  ; send graphics to PNG file                            
    res                  = True                ; plot mods desired                                  
    res@xyLineThicknessF = 4.0                                                                      
    res@tiXAxisString = timeLabel ;"Days"            
    res@tiYAxisString = axisLabel 
    res@tiMainString  = " " ;plotFile                    
    res@tiMainFontHeightF = 0.020 ;changes font size       
;    res@trXMinF = t1    
;    res@trXMaxF = t2

    res@tiMainString     = cname 

    ;make the time axis labels years instead of days:        
;    res@tmXBMode        = "Explicit"         
;    res@tmXBValues      = ispan(0,5760,360)      
;    res@tmXBLabels      = (/-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/) ;"" + res@tmXBValues
;    res@tmXBMinorValues = ispan(0,5760,90)
;    res@tiXAxisSTring   = "Years"
    
    plot = gsn_csm_xy(wks,times,datp_lAvg,res)

  end if ;lat av yes
;
;    
;If not lat av, then do contour plots (lat vs time)
    if(latav.eq."no") then
;reorder for plotting:                                                                          
      datp := datp(lat|:,time|:)
      printVarSummary(datp)    
      
;Plotting:  
      plotFile = datPlot+"_"+opt+"_lat-time-"+cname ;
      
      wks   = gsn_open_wks ("png",plotFile )  ; send graphics to PNG file 
      res   = True                ; plot mods desired

;For contours, define the min/max values and step between values:
      datmin = min(datp)
      datmax = max(datp)    
      datstep = 0.1

;For blue/red style:
      if(plot_redBlue.eq."yes")then
        grid_center = 0
        cmapsize = 97 ;dependent on color map
        cmapcenter = 49 ;dependent on color map
        
        gsn_range = gsnColorRange(datmin, datmax, datstep, grid_center, cmapsize, cmapcenter)
        res@gsnMaximize = True
        colors = "BlRe"
        cmap = read_colormap_file(colors) ;these lines are needed to flip color map around
        cmap = cmap(::-1,:)               ; so that red is negative values
        res@cnFillPalette = cmap(gsn_range@ColorStart-2:gsn_range@ColorEnd-2,:) ;colors
      end if
;--
;For other style:
      if(plot_redBlue.eq."no")then
        res@cnFillPalette = colors
      end if
      
      res@cnFillOn         = True                ; color on
;  res@lbLabelStride    = 4                   ; every other label
      res@lbOrientation    = "Vertical"          ; vertical label bar
;  res@pmLabelBarOrthogonalPosF = -0.04       ; move label bar closer to axis
      res@cnLinesOn        = False               ; turn off contour lines
      
      res@tiMainString     = cname+"~C~   level "+levChoice+" - alt "+levChoiceAlt+" km" ;yr ;-1935

      res@tiXAxisString = timeLabel ;"Days"                             
      res@tiYAxisString = "Latitude" 
      
;  res@vpXF             = 0.12                ; default is 0.2 (aspect ratio)
;  res@vpYF             = 0.8                 ; default is 0.8
;  res@vpHeightF        = 0.4                 ; default is 0.6
;  res@vpWidthF         = 0.75                ; default is 0.6
      
;  res@cnLevelSelectionMode = "ManualLevels" ; manual levels
;  res@cnMinLevelValF       =  140
;  res@cnMaxLevelValF       =  420
;  res@cnLevelSpacingF      =  20 
      
      plot = gsn_csm_lat_time(wks,datp(:,:),res) 
      
    end if ;lat av no
    
end