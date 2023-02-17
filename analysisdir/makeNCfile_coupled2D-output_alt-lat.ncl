;This code takes multi-dimensional output from GSFC 2D coupled model;
;  that output includes chemical constituent concentrations,
;  as well as temperature and some dynamics variables,
;  all as function of altitude and latitude and time (in days)                       
;
; This code assumes input file has 2 years of daily data (with 360 day years);
;
;The code here reads in a single-column text file, divides the data
;  into individual variables (separate ones for each constituent),
;  and then creates an netCDF file with this data and associated attributes.
;The output file is more usable (contains metadata, etc.), and smaller.
;
;BCThomas, June 2018

begin
  
;  in_dir      = "/storage/atmoRuns/test11June2018"
;  in_file1    = "BToutput_alt-lat_daily_test11June2018.out"

  files = new((/2/),string)
  files = asciiread("nclInput_alt-lat.txt",(/2,1/),"string")
  in_dir = files(0)
  in_file1 = files(1)

  print(in_dir+" , "+in_file1)

;some parameter values:                                
  timemax = 720   ;number of time writes
  latmax  = 45    ;number of latitude bins
  altmax  = 76    ;number of altitude bins
  chemmax = 80 ;65    ;number of chemical constituents output
  varmax  = 6     ;number of other variables (e.g. temp, winds)

  prewrites = altmax+latmax+altmax+altmax*latmax ;first we have the press(alt) array, 
                          ;then the lat4(lat) array, 
                          ;then the zalt90(alt) array, then the zkm(lat,alt) array

  ncols       = 1	    ;number of columns - all data in this file is in one column
  nrows       = prewrites+timemax*2+timemax*altmax+timemax*altmax*latmax+timemax*altmax*latmax*(chemmax+varmax)   ;number of rows
  ;here, number of rows is array pre-writes, 
  ; and then 2 time index values at each time write, 
  ; and then an alt value at each time, 
  ; and then a lat value at each alt, at each time, 
  ; and then the list of chem constits and other variables at each lat, alt, and time.

  in_data1    = asciiread(in_dir+"/"+in_file1,(/nrows/),"float")	;array for reading in data
;for this data, don't need fill value, all values are valid.
;  in_data1@_FillValue = 0.000000E+00 ;input file has this as fill
;  in_data1@_FillValue = -9999.99     ;change fill value to this

;Now, I want to put various sets of data into their own arrays:
;array for pressure level values:
  pressarr = new((/altmax/),float)
;Fill this array from data:
  in1 := altmax-1
  pressarr = in_data1(0:in1)
;
  pressarr@long_name = "pressure"
  pressarr@units = "mb"
  pressarr@_FillValue = -9999.99
  pressarr@missing_value = -9999.99
  pressarr!0 = "lev"
  ;printVarSummary(pressarr)
;
;array for latitude values:
  lat4arr = new((/latmax/),float)
;Fill this array from data:
  in2 := in1+1
  in3 := in2+latmax-1
  lat4arr = in_data1(in2:in3)
;
  lat4arr@long_name = "latitude"
  lat4arr@units     = "degrees_north"
  lat4arr@_FillValue = -9999.99
  lat4arr@missing_value = -9999.99
  lat4arr!0 = "lat"
  ;printVarSummary(lat4arr)
;
;array for fixed-grid altitude levels:
  zalt90arr = new((/altmax/),float)
;Fill this array from data:
  in4 := in3+1
  in5 := in4+altmax-1
  zalt90arr = in_data1(in4:in5)
;
  zalt90arr@long_name = "altitude (fixed grid)"
  zalt90arr@units     = "km"
  zalt90arr@_FillValue = -9999.99
  zalt90arr@missing_value = -9999.99
  zalt90arr!0 = "alt_fixed"
  ;printVarSummary(zalt90arr)
  ;print(zalt90arr)
;
;array for dynamic real altitude (determined from space-time dependent scale height/temperature):
  zkmarr = new((/latmax,altmax/),float)
;Fill this array from data:
  n := in5+1
  do k=0,altmax-1
    do j=0,latmax-1
      zkmarr(j,k) = in_data1(n) ;last index count here should be 3540+76=3616
      n := n+1 ;last value here should be 3541+76=3617
    end do
  end do
;
  zkmarr@long_name = "altitude (dynamic, determined from pressure)"
  zkmarr@units     = "km"
  zkmarr@_FillValue = -9999.99
  zkmarr@missing_value = -9999.99
  zkmarr!0 = "alt_dyn"
  ;printVarSummary(zkmarr)
  ;print(zkmarr)
;
;
;array for days of year (360-day year)
  iday360 = new((/timemax/),float)
;array for chemical constituents (e.g. O3, NO2...)
  chemarr = new((/timemax,chemmax,latmax,altmax/),float)
;array for other variables (e.g. temperature)
  vararr  = new((/timemax,varmax,latmax,altmax/),float)
;
;Now, fill these arrays:
  ;n here should be 3541+76=3617
  ;print("this value should be 3617:  "+n)
  tc = 0
  do t=0,timemax-1
    ;there's a time index value we don't need, then the iday360 value (day of year)
    ;print("t "+t+" "+n+" "+" ival "+in_data1(n))
    n := n+1 ;move to next item
    ;print(" ")
    ;print("t "+t+" "+n+" "+" iday360val "+in_data1(n))
    iday360(t) = in_data1(n)
    n := n+1 ;and move to next item
    do k=0,altmax-1
      ;if(t.eq.tc) then
;      if(k.eq.0) then
;        print("k "+t+" "+k+" "+n+" "+in_data1(n))
;      end if
      n := n+1 ;skip the altitude index value
;      ;print(n+" "+in_data1(n))
      do j=0,latmax-1
        ;if(t.eq.tc) then 
;        if(k.eq.0.and.j.eq.0) then 
;          print("j "+t+" "+k+" "+j+"  "+n+" "+in_data1(n))
;        end if
        n := n+1 ;skip the latitude index value
        do s=0,chemmax-1
          ;print("s "+t+" "+k+" "+j+" "+s+"  "+n+" "+in_data1(n))
          chemarr(t,s,j,k) = in_data1(n) 
          n := n+1 ;move to next item
;          ;print("chemarrr: "+chemarr(t,s,j,k))
        end do ;chemicals, index s
        do v=0,varmax-1
          ;print("v "+t+" "+k+" "+j+" "+s+" "+v+"  "+n+" "+in_data1(n))
          vararr(t,v,j,k)  = in_data1(n)          
          n := n+1 ;move to next item
        end do ;other variables, index v
      end do ;lat
    end do ;alt
;    n := n+1
  end do ;time
  print("last n value: "+n)
;
  iday360@long_name = "day of year (360-day year)"
  iday360@units = "days"
  iday360@_FillValue = -9999.99
  iday360@missing_value = -9999.99
  iday360!0 = "day_of_year"
  ;printVarSummary(iday360)
;
;END of read-in and array-filling.

;Now, we'll take the chemical constituent array and make variables for each constit:
  O_3p = new((/timemax,latmax,altmax/),float)
  O_3p = chemarr(:,0,:,:)
  O_3p@long_name = "O(3p) concentration"
  O_3p@units     = "molecules cm^-3"
  O_3p@_FillValue = -9999.99
  O_3p@missing_value = -9999.99
  O_3p!0 = "time"
  O_3p!1 = "lat"
  O_3p!2 = "lev"
  O_3p&lat = lat4arr
  O_3p&lev = pressarr
  ;printVarSummary(O_3p)
;  print(O_3p)

  O_1d = new((/timemax,latmax,altmax/),float)
  O_1d = chemarr(:,1,:,:)
  O_1d@long_name = "O(1d) concentration"
  O_1d@units     = "molecules cm^-3"
  O_1d@_FillValue = -9999.99
  O_1d@missing_value = -9999.99
  O_1d!0 = "time"
  O_1d!1 = "lat"
  O_1d!2 = "lev"
  O_1d&lat = lat4arr
  O_1d&lev = pressarr
  ;printVarSummary(O_1d)
;  print(O_1d)

  O2 = new((/timemax,latmax,altmax/),float)
  O2 = chemarr(:,2,:,:)
  O2@long_name = "O2 concentration"
  O2@units     = "molecules cm^-3"
  O2@_FillValue = -9999.99
  O2@missing_value = -9999.99
  O2!0 = "time"
  O2!1 = "lat"
  O2!2 = "lev"
  O2&lat = lat4arr
  O2&lev = pressarr
  ;printVarSummary(O2)
;  print(O2)

  O3 = new((/timemax,latmax,altmax/),float)
  O3 = chemarr(:,3,:,:)
  O3@long_name = "O3 concentration"
  O3@units     = "molecules cm^-3"
  O3@_FillValue = -9999.99
  O3@missing_value = -9999.99
  O3!0 = "time"
  O3!1 = "lat"
  O3!2 = "lev"
  O3&lat = lat4arr
  O3&lev = pressarr
  ;printVarSummary(O3)
;  print(O3)

  NO = new((/timemax,latmax,altmax/),float)
  NO = chemarr(:,4,:,:)
  NO@long_name = "NO concentration"
  NO@units     = "molecules cm^-3"
  NO@_FillValue = -9999.99
  NO@missing_value = -9999.99
  NO!0 = "time"
  NO!1 = "lat"
  NO!2 = "lev"
  NO&lat = lat4arr
  NO&lev = pressarr
  ;printVarSummary(NO)
;  print(NO)

  NO2 = new((/timemax,latmax,altmax/),float)
  NO2 = chemarr(:,5,:,:)
  NO2@long_name = "NO2 concentration"
  NO2@units     = "molecules cm^-3"
  NO2@_FillValue = -9999.99
  NO2@missing_value = -9999.99
  NO2!0 = "time"
  NO2!1 = "lat"
  NO2!2 = "lev"
  NO2&lat = lat4arr
  NO2&lev = pressarr
  ;printVarSummary(NO2)
;  print(NO2)  

  NO3 = new((/timemax,latmax,altmax/),float)
  NO3 = chemarr(:,6,:,:)
  NO3@long_name = "NO3 concentration"
  NO3@units     = "molecules cm^-3"
  NO3@_FillValue = -9999.99
  NO3@missing_value = -9999.99
  NO3!0 = "time"
  NO3!1 = "lat"
  NO3!2 = "lev"
  NO3&lat = lat4arr
  NO3&lev = pressarr
  ;printVarSummary(NO3)
;  print(NO3)

  N2O5 = new((/timemax,latmax,altmax/),float)
  N2O5 = chemarr(:,7,:,:)
  N2O5@long_name = "N2O5 concentration"
  N2O5@units     = "molecules cm^-3"
  N2O5@_FillValue = -9999.99
  N2O5@missing_value = -9999.99
  N2O5!0 = "time"
  N2O5!1 = "lat"
  N2O5!2 = "lev"
  N2O5&lat = lat4arr
  N2O5&lev = pressarr
  ;printVarSummary(N2O5)
;  print(N2O5)

  N = new((/timemax,latmax,altmax/),float)
  N = chemarr(:,8,:,:)
  N@long_name = "N concentration"
  N@units     = "molecules cm^-3"
  N@_FillValue = -9999.99
  N@missing_value = -9999.99
  N!0 = "time"
  N!1 = "lat"
  N!2 = "lev"
  N&lat = lat4arr
  N&lev = pressarr
  ;printVarSummary(N)
;  print(N)

  HNO3 = new((/timemax,latmax,altmax/),float)
  HNO3 = chemarr(:,9,:,:)
  HNO3@long_name = "HNO3 concentration"
  HNO3@units     = "molecules cm^-3"
  HNO3@_FillValue = -9999.99
  HNO3@missing_value = -9999.99
  HNO3!0 = "time"
  HNO3!1 = "lat"
  HNO3!2 = "lev"
  HNO3&lat = lat4arr
  HNO3&lev = pressarr
  ;printVarSummary(HNO3)
;  print(HNO3)

  N2O = new((/timemax,latmax,altmax/),float)
  N2O = chemarr(:,10,:,:)
  N2O@long_name = "N2O concentration"
  N2O@units     = "molecules cm^-3"
  N2O@_FillValue = -9999.99
  N2O@missing_value = -9999.99
  N2O!0 = "time"
  N2O!1 = "lat"
  N2O!2 = "lev"
  N2O&lat = lat4arr
  N2O&lev = pressarr
  ;printVarSummary(N2O)
;  print(N2O)

  H = new((/timemax,latmax,altmax/),float)
  H = chemarr(:,11,:,:)
  H@long_name = "H concentration"
  H@units     = "molecules cm^-3"
  H@_FillValue = -9999.99
  H@missing_value = -9999.99
  H!0 = "time"
  H!1 = "lat"
  H!2 = "lev"
  H&lat = lat4arr
  H&lev = pressarr
  ;printVarSummary(H)
;  print(H)

  OH = new((/timemax,latmax,altmax/),float)
  OH = chemarr(:,12,:,:)
  OH@long_name = "OH concentration"
  OH@units     = "molecules cm^-3"
  OH@_FillValue = -9999.99
  OH@missing_value = -9999.99
  OH!0 = "time"
  OH!1 = "lat"
  OH!2 = "lev"
  OH&lat = lat4arr
  OH&lev = pressarr
  ;printVarSummary(OH)
;  print(OH)
  
  HO2 = new((/timemax,latmax,altmax/),float)
  HO2 = chemarr(:,13,:,:)
  HO2@long_name = "HO2 concentration"
  HO2@units     = "molecules cm^-3"
  HO2@_FillValue = -9999.99
  HO2@missing_value = -9999.99
  HO2!0 = "time"
  HO2!1 = "lat"
  HO2!2 = "lev"
  HO2&lat = lat4arr
  HO2&lev = pressarr
  ;printVarSummary(HO2)
;  print(HO2)

  H2O = new((/timemax,latmax,altmax/),float)
  H2O = chemarr(:,14,:,:)
  H2O@long_name = "H2O concentration"
  H2O@units     = "molecules cm^-3"
  H2O@_FillValue = -9999.99
  H2O@missing_value = -9999.99
  H2O!0 = "time"
  H2O!1 = "lat"
  H2O!2 = "lev"
  H2O&lat = lat4arr
  H2O&lev = pressarr
  ;printVarSummary(H2O)
;  print(H2O)

  H2O2 = new((/timemax,latmax,altmax/),float)
  H2O2 = chemarr(:,15,:,:)
  H2O2@long_name = "H2O2 concentration"
  H2O2@units     = "molecules cm^-3"
  H2O2@_FillValue = -9999.99
  H2O2@missing_value = -9999.99
  H2O2!0 = "time"
  H2O2!1 = "lat"
  H2O2!2 = "lev"
  H2O2&lat = lat4arr
  H2O2&lev = pressarr
  ;printVarSummary(H2O2)
;  print(H2O2)

  H2 = new((/timemax,latmax,altmax/),float)
  H2 = chemarr(:,16,:,:)
  H2@long_name = "H2 concentration"
  H2@units     = "molecules cm^-3"
  H2@_FillValue = -9999.99
  H2@missing_value = -9999.99
  H2!0 = "time"
  H2!1 = "lat"
  H2!2 = "lev"
  H2&lat = lat4arr
  H2&lev = pressarr
  ;printVarSummary(H2)
;  print(H2)

  CH4 = new((/timemax,latmax,altmax/),float)
  CH4 = chemarr(:,17,:,:)
  CH4@long_name = "CH4 concentration"
  CH4@units     = "molecules cm^-3"
  CH4@_FillValue = -9999.99
  CH4@missing_value = -9999.99
  CH4!0 = "time"
  CH4!1 = "lat"
  CH4!2 = "lev"
  CH4&lat = lat4arr
  CH4&lev = pressarr
  ;printVarSummary(CH4)
;  print(CH4)

  CO = new((/timemax,latmax,altmax/),float)
  CO = chemarr(:,18,:,:)
  CO@long_name = "CO concentration"
  CO@units     = "molecules cm^-3"
  CO@_FillValue = -9999.99
  CO@missing_value = -9999.99
  CO!0 = "time"
  CO!1 = "lat"
  CO!2 = "lev"
  CO&lat = lat4arr
  CO&lev = pressarr
  ;printVarSummary(CO)
;  print(CO)

  CO2 = new((/timemax,latmax,altmax/),float)
  CO2 = chemarr(:,19,:,:)
  CO2@long_name = "CO2 concentration"
  CO2@units     = "molecules cm^-3"
  CO2@_FillValue = -9999.99
  CO2@missing_value = -9999.99
  CO2!0 = "time"
  CO2!1 = "lat"
  CO2!2 = "lev"
  CO2&lat = lat4arr
  CO2&lev = pressarr
  ;printVarSummary(CO2)
;  print(CO2)   

  ClOx = new((/timemax,latmax,altmax/),float)
  ClOx = chemarr(:,20,:,:)
  ClOx@long_name = "ClOx concentration"
  ClOx@units     = "molecules cm^-3"
  ClOx@_FillValue = -9999.99
  ClOx@missing_value = -9999.99
  ClOx!0 = "time"
  ClOx!1 = "lat"
  ClOx!2 = "lev"
  ClOx&lat = lat4arr
  ClOx&lev = pressarr
  ;printVarSummary(ClOx)
;  print(ClOx)   
     
  CH3O2 = new((/timemax,latmax,altmax/),float)
  CH3O2 = chemarr(:,21,:,:)
  CH3O2@long_name = "CH3O2 concentration"
  CH3O2@units     = "molecules cm^-3"
  CH3O2@_FillValue = -9999.99
  CH3O2@missing_value = -9999.99
  CH3O2!0 = "time"
  CH3O2!1 = "lat"
  CH3O2!2 = "lev"
  CH3O2&lat = lat4arr
  CH3O2&lev = pressarr
  ;printVarSummary(CH3O2)
;  print(CH3O2)

  CH2O = new((/timemax,latmax,altmax/),float)
  CH2O = chemarr(:,22,:,:)
  CH2O@long_name = "CH2O concentration"
  CH2O@units     = "molecules cm^-3"
  CH2O@_FillValue = -9999.99
  CH2O@missing_value = -9999.99
  CH2O!0 = "time"
  CH2O!1 = "lat"
  CH2O!2 = "lev"
  CH2O&lat = lat4arr
  CH2O&lev = pressarr
  ;printVarSummary(CH2O)
;  print(CH2O)

  CH3OOH = new((/timemax,latmax,altmax/),float)
  CH3OOH = chemarr(:,23,:,:)
  CH3OOH@long_name = "CH3OOH concentration"
  CH3OOH@units     = "molecules cm^-3"
  CH3OOH@_FillValue = -9999.99
  CH3OOH@missing_value = -9999.99
  CH3OOH!0 = "time"
  CH3OOH!1 = "lat"
  CH3OOH!2 = "lev"
  CH3OOH&lat = lat4arr
  CH3OOH&lev = pressarr
  ;printVarSummary(CH3OOH)
;  print(CH3OOH)  

  HOCl = new((/timemax,latmax,altmax/),float)
  HOCl = chemarr(:,24,:,:)
  HOCl@long_name = "HOCl concentration"
  HOCl@units     = "molecules cm^-3"
  HOCl@_FillValue = -9999.99
  HOCl@missing_value = -9999.99
  HOCl!0 = "time"
  HOCl!1 = "lat"
  HOCl!2 = "lev"
  HOCl&lat = lat4arr
  HOCl&lev = pressarr
  ;printVarSummary(HOCl)
;  print(HOCl)

  ClO3 = new((/timemax,latmax,altmax/),float)
  ClO3 = chemarr(:,25,:,:)
  ClO3@long_name = "ClO3 concentration"
  ClO3@units     = "molecules cm^-3"
  ClO3@_FillValue = -9999.99
  ClO3@missing_value = -9999.99
  ClO3!0 = "time"
  ClO3!1 = "lat"
  ClO3!2 = "lev"
  ClO3&lat = lat4arr
  ClO3&lev = pressarr
  ;printVarSummary(ClO3)
;  print(ClO3)

  Cl = new((/timemax,latmax,altmax/),float)
  Cl = chemarr(:,26,:,:)
  Cl@long_name = "Cl concentration"
  Cl@units     = "molecules cm^-3"
  Cl@_FillValue = -9999.99
  Cl@missing_value = -9999.99
  Cl!0 = "time"
  Cl!1 = "lat"
  Cl!2 = "lev"
  Cl&lat = lat4arr
  Cl&lev = pressarr
  ;printVarSummary(Cl)
;  print(Cl)  

  ClO = new((/timemax,latmax,altmax/),float)
  ClO = chemarr(:,27,:,:)
  ClO@long_name = "ClO concentration"
  ClO@units     = "molecules cm^-3"
  ClO@_FillValue = -9999.99
  ClO@missing_value = -9999.99
  ClO!0 = "time"
  ClO!1 = "lat"
  ClO!2 = "lev"
  ClO&lat = lat4arr
  ClO&lev = pressarr
  ;printVarSummary(ClO)
;  print(ClO)  

  HCl = new((/timemax,latmax,altmax/),float)
  HCl = chemarr(:,28,:,:)
  HCl@long_name = "HCl concentration"
  HCl@units     = "molecules cm^-3"
  HCl@_FillValue = -9999.99
  HCl@missing_value = -9999.99
  HCl!0 = "time"
  HCl!1 = "lat"
  HCl!2 = "lev"
  HCl&lat = lat4arr
  HCl&lev = pressarr
  ;printVarSummary(HCl)
;  print(HCl)  

  ClONO2 = new((/timemax,latmax,altmax/),float)
  ClONO2 = chemarr(:,29,:,:)
  ClONO2@long_name = "ClONO2 concentration"
  ClONO2@units     = "molecules cm^-3"
  ClONO2@_FillValue = -9999.99
  ClONO2@missing_value = -9999.99
  ClONO2!0 = "time"
  ClONO2!1 = "lat"
  ClONO2!2 = "lev"
  ClONO2&lat = lat4arr
  ClONO2&lev = pressarr
  ;printVarSummary(ClONO2)
;  print(ClONO2)  

  NOy = new((/timemax,latmax,altmax/),float)
  NOy = chemarr(:,30,:,:)
  NOy@long_name = "NOy concentration"
  NOy@units     = "molecules cm^-3"
  NOy@_FillValue = -9999.99
  NOy@missing_value = -9999.99
  NOy!0 = "time"
  NOy!1 = "lat"
  NOy!2 = "lev"
  NOy&lat = lat4arr
  NOy&lev = pressarr
  ;printVarSummary(NOy)
;  print(NOy)  

  NOx = new((/timemax,latmax,altmax/),float)
  NOx = chemarr(:,31,:,:)
  NOx@long_name = "NOx concentration"
  NOx@units     = "molecules cm^-3"
  NOx@_FillValue = -9999.99
  NOx@missing_value = -9999.99
  NOx!0 = "time"
  NOx!1 = "lat"
  NOx!2 = "lev"
  NOx&lat = lat4arr
  NOx&lev = pressarr
  ;printVarSummary(NOx)
;  print(NOx)  

  Cly = new((/timemax,latmax,altmax/),float)
  Cly = chemarr(:,32,:,:)
  Cly@long_name = "Cly concentration"
  Cly@units     = "molecules cm^-3"
  Cly@_FillValue = -9999.99
  Cly@missing_value = -9999.99
  Cly!0 = "time"
  Cly!1 = "lat"
  Cly!2 = "lev"
  Cly&lat = lat4arr
  Cly&lev = pressarr
  ;printVarSummary(Cly)
;  print(Cly)  

  CFCl3 = new((/timemax,latmax,altmax/),float)
  CFCl3 = chemarr(:,33,:,:)
  CFCl3@long_name = "CFCl3 concentration"
  CFCl3@units     = "molecules cm^-3"
  CFCl3@_FillValue = -9999.99
  CFCl3@missing_value = -9999.99
  CFCl3!0 = "time"
  CFCl3!1 = "lat"
  CFCl3!2 = "lev"
  CFCl3&lat = lat4arr
  CFCl3&lev = pressarr
  ;printVarSummary(CFCl3)
;  print(CFCl3)  

  CF2Cl2 = new((/timemax,latmax,altmax/),float)
  CF2Cl2 = chemarr(:,34,:,:)
  CF2Cl2@long_name = "CF2Cl2 concentration"
  CF2Cl2@units     = "molecules cm^-3"
  CF2Cl2@_FillValue = -9999.99
  CF2Cl2@missing_value = -9999.99
  CF2Cl2!0 = "time"
  CF2Cl2!1 = "lat"
  CF2Cl2!2 = "lev"
  CF2Cl2&lat = lat4arr
  CF2Cl2&lev = pressarr
  ;printVarSummary(CF2Cl2)
;  print(CF2Cl2)  

  CCl4 = new((/timemax,latmax,altmax/),float)
  CCl4 = chemarr(:,35,:,:)
  CCl4@long_name = "CCl4 concentration"
  CCl4@units     = "molecules cm^-3"
  CCl4@_FillValue = -9999.99
  CCl4@missing_value = -9999.99
  CCl4!0 = "time"
  CCl4!1 = "lat"
  CCl4!2 = "lev"
  CCl4&lat = lat4arr
  CCl4&lev = pressarr
  ;printVarSummary(CCl4)
;  print(CCl4)  

  CH3Cl = new((/timemax,latmax,altmax/),float)
  CH3Cl = chemarr(:,36,:,:)
  CH3Cl@long_name = "CH3Cl concentration"
  CH3Cl@units     = "molecules cm^-3"
  CH3Cl@_FillValue = -9999.99
  CH3Cl@missing_value = -9999.99
  CH3Cl!0 = "time"
  CH3Cl!1 = "lat"
  CH3Cl!2 = "lev"
  CH3Cl&lat = lat4arr
  CH3Cl&lev = pressarr
  ;printVarSummary(CH3Cl)
;  print(CH3Cl)  

  HO2NO2 = new((/timemax,latmax,altmax/),float)
  HO2NO2 = chemarr(:,37,:,:)
  HO2NO2@long_name = "HO2NO2 concentration"
  HO2NO2@units     = "molecules cm^-3"
  HO2NO2@_FillValue = -9999.99
  HO2NO2@missing_value = -9999.99
  HO2NO2!0 = "time"
  HO2NO2!1 = "lat"
  HO2NO2!2 = "lev"
  HO2NO2&lat = lat4arr
  HO2NO2&lev = pressarr
  ;printVarSummary(HO2NO2)
;  print(HO2NO2)  

  Ox = new((/timemax,latmax,altmax/),float)
  Ox = chemarr(:,38,:,:)
  Ox@long_name = "Ox concentration"
  Ox@units     = "molecules cm^-3"
  Ox@_FillValue = -9999.99
  Ox@missing_value = -9999.99
  Ox!0 = "time"
  Ox!1 = "lat"
  Ox!2 = "lev"
  Ox&lat = lat4arr
  Ox&lev = pressarr
  ;printVarSummary(Ox)
;  print(Ox)    

  CH3CCl3 = new((/timemax,latmax,altmax/),float)
  CH3CCl3 = chemarr(:,39,:,:)
  CH3CCl3@long_name = "CH3CCl3 concentration"
  CH3CCl3@units     = "molecules cm^-3"
  CH3CCl3@_FillValue = -9999.99
  CH3CCl3@missing_value = -9999.99
  CH3CCl3!0 = "time"
  CH3CCl3!1 = "lat"
  CH3CCl3!2 = "lev"
  CH3CCl3&lat = lat4arr
  CH3CCl3&lev = pressarr
  ;printVarSummary(CH3CCl3)
;  print(CH3CCl3)  

  O2_1delta = new((/timemax,latmax,altmax/),float)
  O2_1delta = chemarr(:,40,:,:)
  O2_1delta@long_name = "O2(1delta) concentration"
  O2_1delta@units     = "molecules cm^-3"
  O2_1delta@_FillValue = -9999.99
  O2_1delta@missing_value = -9999.99
  O2_1delta!0 = "time"
  O2_1delta!1 = "lat"
  O2_1delta!2 = "lev"
  O2_1delta&lat = lat4arr
  O2_1delta&lev = pressarr
  ;printVarSummary(O2_1delta)
;  print(O2_1delta)  

  HONO = new((/timemax,latmax,altmax/),float)
  HONO = chemarr(:,41,:,:)
  HONO@long_name = "HONO concentration"
  HONO@units     = "molecules cm^-3"
  HONO@_FillValue = -9999.99
  HONO@missing_value = -9999.99
  HONO!0 = "time"
  HONO!1 = "lat"
  HONO!2 = "lev"
  HONO&lat = lat4arr
  HONO&lev = pressarr
  ;printVarSummary(HONO)
;  print(HONO)    

  Br2 = new((/timemax,latmax,altmax/),float)
  Br2 = chemarr(:,42,:,:)
  Br2@long_name = "Br2 concentration"
  Br2@units     = "molecules cm^-3"
  Br2@_FillValue = -9999.99
  Br2@missing_value = -9999.99
  Br2!0 = "time"
  Br2!1 = "lat"
  Br2!2 = "lev"
  Br2&lat = lat4arr
  Br2&lev = pressarr
  ;printVarSummary(Br2)
;  print(Br2)    

  BrO = new((/timemax,latmax,altmax/),float)
  BrO = chemarr(:,43,:,:)
  BrO@long_name = "BrO concentration"
  BrO@units     = "molecules cm^-3"
  BrO@_FillValue = -9999.99
  BrO@missing_value = -9999.99
  BrO!0 = "time"
  BrO!1 = "lat"
  BrO!2 = "lev"
  BrO&lat = lat4arr
  BrO&lev = pressarr
  ;printVarSummary(BrO)
;  print(BrO)   

  Br = new((/timemax,latmax,altmax/),float)
  Br = chemarr(:,44,:,:)
  Br@long_name = "Br concentration"
  Br@units     = "molecules cm^-3"
  Br@_FillValue = -9999.99
  Br@missing_value = -9999.99
  Br!0 = "time"
  Br!1 = "lat"
  Br!2 = "lev"
  Br&lat = lat4arr
  Br&lev = pressarr
  ;printVarSummary(Br)
;  print(Br)   

  HBr = new((/timemax,latmax,altmax/),float)
  HBr = chemarr(:,45,:,:)
  HBr@long_name = "HBr concentration"
  HBr@units     = "molecules cm^-3"
  HBr@_FillValue = -9999.99
  HBr@missing_value = -9999.99
  HBr!0 = "time"
  HBr!1 = "lat"
  HBr!2 = "lev"
  HBr&lat = lat4arr
  HBr&lev = pressarr
  ;printVarSummary(HBr)
;  print(HBr)   

  BrONO2 = new((/timemax,latmax,altmax/),float)
  BrONO2 = chemarr(:,46,:,:)
  BrONO2@long_name = "BrONO2 concentration"
  BrONO2@units     = "molecules cm^-3"
  BrONO2@_FillValue = -9999.99
  BrONO2@missing_value = -9999.99
  BrONO2!0 = "time"
  BrONO2!1 = "lat"
  BrONO2!2 = "lev"
  BrONO2&lat = lat4arr
  BrONO2&lev = pressarr
  ;printVarSummary(BrONO2)
;  print(BrONO2)   

  Bry = new((/timemax,latmax,altmax/),float)
  Bry = chemarr(:,47,:,:)
  Bry@long_name = "Bry concentration"
  Bry@units     = "molecules cm^-3"
  Bry@_FillValue = -9999.99
  Bry@missing_value = -9999.99
  Bry!0 = "time"
  Bry!1 = "lat"
  Bry!2 = "lev"
  Bry&lat = lat4arr
  Bry&lev = pressarr
  ;printVarSummary(Bry)
;  print(Bry)   

  CH3Br = new((/timemax,latmax,altmax/),float)
  CH3Br = chemarr(:,48,:,:)
  CH3Br@long_name = "CH3Br concentration"
  CH3Br@units     = "molecules cm^-3"
  CH3Br@_FillValue = -9999.99
  CH3Br@missing_value = -9999.99
  CH3Br!0 = "time"
  CH3Br!1 = "lat"
  CH3Br!2 = "lev"
  CH3Br&lat = lat4arr
  CH3Br&lev = pressarr
  ;printVarSummary(CH3Br)
;  print(CH3Br)   

  CBrF3 = new((/timemax,latmax,altmax/),float)
  CBrF3 = chemarr(:,49,:,:)
  CBrF3@long_name = "CBrF3 concentration"
  CBrF3@units     = "molecules cm^-3"
  CBrF3@_FillValue = -9999.99
  CBrF3@missing_value = -9999.99
  CBrF3!0 = "time"
  CBrF3!1 = "lat"
  CBrF3!2 = "lev"
  CBrF3&lat = lat4arr
  CBrF3&lev = pressarr
  ;printVarSummary(CBrF3)
;  print(CBrF3)   

  CBrClF2 = new((/timemax,latmax,altmax/),float)
  CBrClF2 = chemarr(:,50,:,:)
  CBrClF2@long_name = "CBrClF2 concentration"
  CBrClF2@units     = "molecules cm^-3"
  CBrClF2@_FillValue = -9999.99
  CBrClF2@missing_value = -9999.99
  CBrClF2!0 = "time"
  CBrClF2!1 = "lat"
  CBrClF2!2 = "lev"
  CBrClF2&lat = lat4arr
  CBrClF2&lev = pressarr
  ;printVarSummary(CBrClF2)
;  print(CBrClF2)   

  CHClF2 = new((/timemax,latmax,altmax/),float)
  CHClF2 = chemarr(:,51,:,:)
  CHClF2@long_name = "CHClF2 concentration"
  CHClF2@units     = "molecules cm^-3"
  CHClF2@_FillValue = -9999.99
  CHClF2@missing_value = -9999.99
  CHClF2!0 = "time"
  CHClF2!1 = "lat"
  CHClF2!2 = "lev"
  CHClF2&lat = lat4arr
  CHClF2&lev = pressarr
  ;printVarSummary(CHClF2)
;  print(CHClF2)   

  C2Cl3F3 = new((/timemax,latmax,altmax/),float)
  C2Cl3F3 = chemarr(:,52,:,:)
  C2Cl3F3@long_name = "C2Cl3F3 concentration"
  C2Cl3F3@units     = "molecules cm^-3"
  C2Cl3F3@_FillValue = -9999.99
  C2Cl3F3@missing_value = -9999.99
  C2Cl3F3!0 = "time"
  C2Cl3F3!1 = "lat"
  C2Cl3F3!2 = "lev"
  C2Cl3F3&lat = lat4arr
  C2Cl3F3&lev = pressarr
  ;printVarSummary(C2Cl3F3)
;  print(C2Cl3F3)   

  C2Cl2F4 = new((/timemax,latmax,altmax/),float)
  C2Cl2F4 = chemarr(:,53,:,:)
  C2Cl2F4@long_name = "C2Cl2F4 concentration"
  C2Cl2F4@units     = "molecules cm^-3"
  C2Cl2F4@_FillValue = -9999.99
  C2Cl2F4@missing_value = -9999.99
  C2Cl2F4!0 = "time"
  C2Cl2F4!1 = "lat"
  C2Cl2F4!2 = "lev"
  C2Cl2F4&lat = lat4arr
  C2Cl2F4&lev = pressarr
  ;printVarSummary(C2Cl2F4)
;  print(C2Cl2F4) 

  C2ClF5 = new((/timemax,latmax,altmax/),float)
  C2ClF5 = chemarr(:,54,:,:)
  C2ClF5@long_name = "C2ClF5 concentration"
  C2ClF5@units     = "molecules cm^-3"
  C2ClF5@_FillValue = -9999.99
  C2ClF5@missing_value = -9999.99
  C2ClF5!0 = "time"
  C2ClF5!1 = "lat"
  C2ClF5!2 = "lev"
  C2ClF5&lat = lat4arr
  C2ClF5&lev = pressarr
  ;printVarSummary(C2ClF5)
;  print(C2ClF5) 

  HF = new((/timemax,latmax,altmax/),float)
  HF = chemarr(:,55,:,:)
  HF@long_name = "HF concentration"
  HF@units     = "molecules cm^-3"
  HF@_FillValue = -9999.99
  HF@missing_value = -9999.99
  HF!0 = "time"
  HF!1 = "lat"
  HF!2 = "lev"
  HF&lat = lat4arr
  HF&lev = pressarr
  ;printVarSummary(HF)
;  print(HF) 

  CF2O = new((/timemax,latmax,altmax/),float)
  CF2O = chemarr(:,56,:,:)
  CF2O@long_name = "CF2O concentration"
  CF2O@units     = "molecules cm^-3"
  CF2O@_FillValue = -9999.99
  CF2O@missing_value = -9999.99
  CF2O!0 = "time"
  CF2O!1 = "lat"
  CF2O!2 = "lev"
  CF2O&lat = lat4arr
  CF2O&lev = pressarr
  ;printVarSummary(CF2O)
;  print(CF2O) 

  BrOx = new((/timemax,latmax,altmax/),float)
  BrOx = chemarr(:,57,:,:)
  BrOx@long_name = "BrOx concentration"
  BrOx@units     = "molecules cm^-3"
  BrOx@_FillValue = -9999.99
  BrOx@missing_value = -9999.99
  BrOx!0 = "time"
  BrOx!1 = "lat"
  BrOx!2 = "lev"
  BrOx&lat = lat4arr
  BrOx&lev = pressarr
  ;printVarSummary(BrOx)
;  print(BrOx) 

  H2O_rainout = new((/timemax,latmax,altmax/),float)
  H2O_rainout = chemarr(:,58,:,:)
  H2O_rainout@long_name = "H2O_rainout"
  H2O_rainout@units     = "molecules cm^-3 s^-1"
  H2O_rainout@_FillValue = -9999.99
  H2O_rainout@missing_value = -9999.99
  H2O_rainout!0 = "time"
  H2O_rainout!1 = "lat"
  H2O_rainout!2 = "lev"
  H2O_rainout&lat = lat4arr
  H2O_rainout&lev = pressarr
  ;printVarSummary(H2O_rainout)
;  print(H2O_rainout) 

  BrCl = new((/timemax,latmax,altmax/),float)
  BrCl = chemarr(:,59,:,:)
  BrCl@long_name = "BrCl concentration"
  BrCl@units     = "molecules cm^-3"
  BrCl@_FillValue = -9999.99
  BrCl@missing_value = -9999.99
  BrCl!0 = "time"
  BrCl!1 = "lat"
  BrCl!2 = "lev"
  BrCl&lat = lat4arr
  BrCl&lev = pressarr
  ;printVarSummary(BrCl)
;  print(BrCl)        

  Cl2O2 = new((/timemax,latmax,altmax/),float)
  Cl2O2 = chemarr(:,60,:,:)
  Cl2O2@long_name = "Cl2O2 concentration"
  Cl2O2@units     = "molecules cm^-3"
  Cl2O2@_FillValue = -9999.99
  Cl2O2@missing_value = -9999.99
  Cl2O2!0 = "time"
  Cl2O2!1 = "lat"
  Cl2O2!2 = "lev"
  Cl2O2&lat = lat4arr
  Cl2O2&lev = pressarr
  ;printVarSummary(Cl2O2)
;  print(Cl2O2)        

  OClO = new((/timemax,latmax,altmax/),float)
  OClO = chemarr(:,61,:,:)
  OClO@long_name = "OClO concentration"
  OClO@units     = "molecules cm^-3"
  OClO@_FillValue = -9999.99
  OClO@missing_value = -9999.99
  OClO!0 = "time"
  OClO!1 = "lat"
  OClO!2 = "lev"
  OClO&lat = lat4arr
  OClO&lev = pressarr
  ;printVarSummary(OClO)
;  print(OClO) 

  ClNO2 = new((/timemax,latmax,altmax/),float)
  ClNO2 = chemarr(:,62,:,:)
  ClNO2@long_name = "ClNO2 concentration"
  ClNO2@units     = "molecules cm^-3"
  ClNO2@_FillValue = -9999.99
  ClNO2@missing_value = -9999.99
  ClNO2!0 = "time"
  ClNO2!1 = "lat"
  ClNO2!2 = "lev"
  ClNO2&lat = lat4arr
  ClNO2&lev = pressarr
  ;printVarSummary(ClNO2)
;  print(ClNO2) 

  Cl2 = new((/timemax,latmax,altmax/),float)
  Cl2 = chemarr(:,63,:,:)
  Cl2@long_name = "Cl2 concentration"
  Cl2@units     = "molecules cm^-3"
  Cl2@_FillValue = -9999.99
  Cl2@missing_value = -9999.99
  Cl2!0 = "time"
  Cl2!1 = "lat"
  Cl2!2 = "lev"
  Cl2&lat = lat4arr
  Cl2&lev = pressarr
  ;printVarSummary(Cl2)
;  print(Cl2)        

  ClOO = new((/timemax,latmax,altmax/),float)
  ClOO = chemarr(:,64,:,:)
  ClOO@long_name = "ClOO concentration"
  ClOO@units     = "molecules cm^-3"
  ClOO@_FillValue = -9999.99
  ClOO@missing_value = -9999.99
  ClOO!0 = "time"
  ClOO!1 = "lat"
  ClOO!2 = "lev"
  ClOO&lat = lat4arr
  ClOO&lev = pressarr
  ;printVarSummary(ClOO)
;  print(ClOO)        

  HNO3s = new((/timemax,latmax,altmax/),float)
  HNO3s = chemarr(:,65,:,:)
  HNO3s@long_name = "HNO3 (solid) concentration"
  HNO3s@units     = "molecules cm^-3"
  HNO3s@_FillValue = -9999.99
  HNO3s@missing_value = -9999.99
  HNO3s!0 = "time"
  HNO3s!1 = "lat"
  HNO3s!2 = "lev"
  HNO3s&lat = lat4arr
  HNO3s&lev = pressarr
  ;printVarSummary(HNO3s)
;  print(HNO3s)

  LiqH2O = new((/timemax,latmax,altmax/),float)
  LiqH2O = chemarr(:,66,:,:)
  LiqH2O@long_name = "H2O (solid/liquid) concentration"
  LiqH2O@units     = "molecules cm^-3"
  LiqH2O@_FillValue = -9999.99
  LiqH2O@missing_value = -9999.99
  LiqH2O!0 = "time"
  LiqH2O!1 = "lat"
  LiqH2O!2 = "lev"
  LiqH2O&lat = lat4arr
  LiqH2O&lev = pressarr
  ;printVarSummary(LiqH2O)
;  print(LiqH2O)  

  HOBr = new((/timemax,latmax,altmax/),float)
  HOBr = chemarr(:,67,:,:)
  HOBr@long_name = "HOBr concentration"
  HOBr@units     = "molecules cm^-3"
  HOBr@_FillValue = -9999.99
  HOBr@missing_value = -9999.99
  HOBr!0 = "time"
  HOBr!1 = "lat"
  HOBr!2 = "lev"
  HOBr&lat = lat4arr
  HOBr&lev = pressarr
  ;printVarSummary(HOBr)
;  print(HOBr)

  CH3CCl2F = new((/timemax,latmax,altmax/),float)
  CH3CCl2F = chemarr(:,68,:,:)
  CH3CCl2F@long_name = "CH3CCl2F concentration"
  CH3CCl2F@units     = "molecules cm^-3"
  CH3CCl2F@_FillValue = -9999.99
  CH3CCl2F@missing_value = -9999.99
  CH3CCl2F!0 = "time"
  CH3CCl2F!1 = "lat"
  CH3CCl2F!2 = "lev"
  CH3CCl2F&lat = lat4arr
  CH3CCl2F&lev = pressarr
  ;printVarSummary(CH3CCl2F)
;  print(CH3CCl2F)

  CH3CClF2 = new((/timemax,latmax,altmax/),float)
  CH3CClF2 = chemarr(:,69,:,:)
  CH3CClF2@long_name = "CH3CClF2 concentration"
  CH3CClF2@units     = "molecules cm^-3"
  CH3CClF2@_FillValue = -9999.99
  CH3CClF2@missing_value = -9999.99
  CH3CClF2!0 = "time"
  CH3CClF2!1 = "lat"
  CH3CClF2!2 = "lev"
  CH3CClF2&lat = lat4arr
  CH3CClF2&lev = pressarr
  ;printVarSummary(CH3CClF2)
;  print(CH3CClF2)

  CHCl2CF3 = new((/timemax,latmax,altmax/),float)
  CHCl2CF3 = chemarr(:,70,:,:)
  CHCl2CF3@long_name = "CHCl2CF3 concentration"
  CHCl2CF3@units     = "molecules cm^-3"
  CHCl2CF3@_FillValue = -9999.99
  CHCl2CF3@missing_value = -9999.99
  CHCl2CF3!0 = "time"
  CHCl2CF3!1 = "lat"
  CHCl2CF3!2 = "lev"
  CHCl2CF3&lat = lat4arr
  CHCl2CF3&lev = pressarr
  ;printVarSummary(CHCl2CF3)
;  print(CHCl2CF3)

  C2Br2F4 = new((/timemax,latmax,altmax/),float)
  C2Br2F4 = chemarr(:,71,:,:)
  C2Br2F4@long_name = "C2Br2F4 concentration"
  C2Br2F4@units     = "molecules cm^-3"
  C2Br2F4@_FillValue = -9999.99
  C2Br2F4@missing_value = -9999.99
  C2Br2F4!0 = "time"
  C2Br2F4!1 = "lat"
  C2Br2F4!2 = "lev"
  C2Br2F4&lat = lat4arr
  C2Br2F4&lev = pressarr
  ;printVarSummary(C2Br2F4)
;  print(C2Br2F4)

  CHx = new((/timemax,latmax,altmax/),float)
  CHx = chemarr(:,72,:,:)
  CHx@long_name = "CHx concentration"
  CHx@units     = "molecules cm^-3"
  CHx@_FillValue = -9999.99
  CHx@missing_value = -9999.99
  CHx!0 = "time"
  CHx!1 = "lat"
  CHx!2 = "lev"
  CHx&lat = lat4arr
  CHx&lev = pressarr
  ;printVarSummary(CHx)
;  print(CHx) 

  HOx = new((/timemax,latmax,altmax/),float)
  HOx = chemarr(:,73,:,:)
  HOx@long_name = "HOx concentration"
  HOx@units     = "molecules cm^-3"
  HOx@_FillValue = -9999.99
  HOx@missing_value = -9999.99
  HOx!0 = "time"
  HOx!1 = "lat"
  HOx!2 = "lev"
  HOx&lat = lat4arr
  HOx&lev = pressarr
  ;printVarSummary(HOx)
;  print(HOx) 

  CBr2F2 = new((/timemax,latmax,altmax/),float)
  CBr2F2 = chemarr(:,74,:,:)
  CBr2F2@long_name = "CBr2F2 concentration"
  CBr2F2@units     = "molecules cm^-3"
  CBr2F2@_FillValue = -9999.99
  CBr2F2@missing_value = -9999.99
  CBr2F2!0 = "time"
  CBr2F2!1 = "lat"
  CBr2F2!2 = "lev"
  CBr2F2&lat = lat4arr
  CBr2F2&lev = pressarr
  ;printVarSummary(CBr2F2)
;  print(CBr2F2) 

  CH2Br2 = new((/timemax,latmax,altmax/),float)
  CH2Br2 = chemarr(:,75,:,:)
  CH2Br2@long_name = "CH2Br2 concentration"
  CH2Br2@units     = "molecules cm^-3"
  CH2Br2@_FillValue = -9999.99
  CH2Br2@missing_value = -9999.99
  CH2Br2!0 = "time"
  CH2Br2!1 = "lat"
  CH2Br2!2 = "lev"
  CH2Br2&lat = lat4arr
  CH2Br2&lev = pressarr
  ;printVarSummary(CH2Br2)
;  print(CH2Br2) 

  CHBr3 = new((/timemax,latmax,altmax/),float)
  CHBr3 = chemarr(:,76,:,:)
  CHBr3@long_name = "CHBr3 concentration"
  CHBr3@units     = "molecules cm^-3"
  CHBr3@_FillValue = -9999.99
  CHBr3@missing_value = -9999.99
  CHBr3!0 = "time"
  CHBr3!1 = "lat"
  CHBr3!2 = "lev"
  CHBr3&lat = lat4arr
  CHBr3&lev = pressarr
  ;printVarSummary(CHBr3)
;  print(CHBr3) 

  Age = new((/timemax,latmax,altmax/),float)
  Age = chemarr(:,77,:,:)
  Age@long_name = "Clock Tracer (Age) concentration?"
  Age@units     = "molecules cm^-3 ?"
  Age@_FillValue = -9999.99
  Age@missing_value = -9999.99
  Age!0 = "time"
  Age!1 = "lat"
  Age!2 = "lev"
  Age&lat = lat4arr
  Age&lev = pressarr
  ;printVarSummary(Age)
;  print(Age)

  BrONO = new((/timemax,latmax,altmax/),float)
  BrONO = chemarr(:,78,:,:)
  BrONO@long_name = "Clock Tracer (BrONO) concentration"
  BrONO@units     = "molecules cm^-3"
  BrONO@_FillValue = -9999.99
  BrONO@missing_value = -9999.99
  BrONO!0 = "time"
  BrONO!1 = "lat"
  BrONO!2 = "lev"
  BrONO&lat = lat4arr
  BrONO&lev = pressarr
  ;printVarSummary(BrONO)
;  print(BrONO)

  OzoneDiag = new((/timemax,latmax,altmax/),float)
  OzoneDiag = chemarr(:,79,:,:)
  OzoneDiag@long_name = "Ozone Diagnostic"
  OzoneDiag@units     = "DU km^-1"
  OzoneDiag@_FillValue = -9999.99
  OzoneDiag@missing_value = -9999.99
  OzoneDiag!0 = "time"
  OzoneDiag!1 = "lat"
  OzoneDiag!2 = "lev"
  OzoneDiag&lat = lat4arr
  OzoneDiag&lev = pressarr
  ;printVarSummary(OzoneDiag)
;  print(OzoneDiag)

;Now, the other model output variables:
  
  Temp = new((/timemax,latmax,altmax/),float)
  Temp = vararr(:,0,:,:)
  Temp@long_name = "Temperature"
  Temp@units     = "K"
  Temp@_FillValue = -9999.99
  Temp@missing_value = -9999.99
  Temp!0 = "time"
  Temp!1 = "lat"
  Temp!2 = "lev"
  Temp&lat = lat4arr
  Temp&lev = pressarr
  ;printVarSummary(Temp)
;  print(Temp)    

  UbarD = new((/timemax,latmax,altmax/),float)
  UbarD = vararr(:,1,:,:)
  UbarD@long_name = "mean zonal wind speed"
  UbarD@units     = "m s^-1"
  UbarD@_FillValue = -9999.99
  UbarD@missing_value = -9999.99
  UbarD!0 = "time"
  UbarD!1 = "lat"
  UbarD!2 = "lev"
  UbarD&lat = lat4arr
  UbarD&lev = pressarr
  ;printVarSummary(UbarD)
;  print(UbarD) 

  HeatOut = new((/timemax,latmax,altmax/),float)
  HeatOut = vararr(:,2,:,:)
  HeatOut@long_name = "heating rate"
  HeatOut@units     = "K/day"
  HeatOut@_FillValue = -9999.99
  HeatOut@missing_value = -9999.99
  HeatOut!0 = "time"
  HeatOut!1 = "lat"
  HeatOut!2 = "lev"
  HeatOut&lat = lat4arr
  HeatOut&lev = pressarr
  ;printVarSummary(HeatOut)
;  print(HeatOut)     

  CoolOut = new((/timemax,latmax,altmax/),float)
  CoolOut = vararr(:,3,:,:)
  CoolOut@long_name = "cooling rate"
  CoolOut@units     = "K/day"
  CoolOut@_FillValue = -9999.99
  CoolOut@missing_value = -9999.99
  CoolOut!0 = "time"
  CoolOut!1 = "lat"
  CoolOut!2 = "lev"
  CoolOut&lat = lat4arr
  CoolOut&lev = pressarr
  ;printVarSummary(CoolOut)
;  print(CoolOut)                                                                                       
  W = new((/timemax,latmax,altmax/),float)
  W = vararr(:,4,:,:)
  W@long_name = "vertical velocity"
  W@units     = "cm s^-1"
  W@_FillValue = -9999.99
  W@missing_value = -9999.99
  W!0 = "time"
  W!1 = "lat"
  W!2 = "lev"
  W&lat = lat4arr
  W&lev = pressarr
  ;printVarSummary(W)
;  print(W)                                                                                       

  V = new((/timemax,latmax,altmax/),float)
  V = vararr(:,5,:,:)
  V@long_name = "horizontal velocity (north-south)"
  V@units     = "cm s^-1"
  V@_FillValue = -9999.99
  V@missing_value = -9999.99
  V!0 = "time"
  V!1 = "lat"
  V!2 = "lev"
  V&lat = lat4arr
  V&lev = pressarr
  ;printVarSummary(V)
;  print(V) 

; This block creates a NetCDF file out of the data that was read in.
;------------------------------------------------------------
;  out_dir = "./"
;  fn = out_dir+"_"+in_file1+".nc"
  out_dir = in_dir
  fn = in_file1+".nc"
  system("rm -f "+out_dir+"/"+fn) ;remove old version

  out_file    = addfile(out_dir+"/"+fn,"c")
;Global attributes:
  out_file@title	= "Output of GSFC 2D atmo model (coupled version) brian.thomas@washburn.edu" ;set title of file
  out_file@creation_date = systemfunc("date")      ;set creation date from system

;output data values:
  out_file->press = pressarr
  out_file->lat4arr = lat4arr
  out_file->zalt90arr = zalt90arr
  out_file->zkmarr = zkmarr
  out_file->iday360 = iday360
  out_file->O_3p = O_3p
  out_file->O_1d = O_1d
  out_file->O2 = O2
  out_file->O3 = O3
  out_file->NO = NO
  out_file->NO2 = NO2
  out_file->NO3 = NO3
  out_file->N2O5 = N2O5
  out_file->N = N
  out_file->HNO3 = HNO3
  out_file->N2O = N2O
  out_file->H = H
  out_file->OH = OH
  out_file->HO2 = HO2
  out_file->H2O = H2O
  out_file->H2O2 = H2O2
  out_file->H2 = H2
  out_file->CH4 = CH4
  out_file->CO = CO
  out_file->CO2 = CO2
  out_file->ClOx = ClOx
  out_file->CH3O2 = CH3O2
  out_file->CH2O = CH2O
  out_file->CH3OOH = CH3OOH
  out_file->HOCl = HOCl
  out_file->ClO3 = ClO3
  out_file->Cl = Cl
  out_file->ClO = ClO
  out_file->HCl = HCl
  out_file->ClONO2 = ClONO2
  out_file->NOy = NOy
  out_file->NOx = NOx
  out_file->Cly = Cly
  out_file->CFCl3 = CFCl3
  out_file->CF2Cl2 = CF2Cl2
  out_file->CCl4 = CCl4
  out_file->CH3Cl = CH3Cl
  out_file->HO2NO2 = HO2NO2
  out_file->Ox = Ox
  out_file->CH3CCl3 = CH3CCl3
  out_file->O2_1delta = O2_1delta
  out_file->HONO = HONO
  out_file->Br2 = Br2
  out_file->BrO = BrO
  out_file->Br = Br
  out_file->HBr = HBr
  out_file->BrONO2 = BrONO2
  out_file->Bry = Bry
  out_file->CH3Br = CH3Br
  out_file->CBrF3 = CBrF3
  out_file->CBrClF2 = CBrClF2
  out_file->CHClF2 = CHClF2
  out_file->C2Cl3F3 = C2Cl3F3
  out_file->C2Cl2F4 = C2Cl2F4
  out_file->C2ClF5 = C2ClF5
  out_file->HF = HF
  out_file->CF2O = CF2O
  out_file->BrOx = BrOx
  out_file->H2O_rainout = H2O_rainout
  out_file->BrCl = BrCl
  out_file->Cl2O2 = Cl2O2
  out_file->OClO = OClO
  out_file->ClNO2 = ClNO2
  out_file->Cl2 = Cl2
  out_file->ClOO = ClOO
  out_file->HNO3s = HNO3s
  out_file->LiqH2O = LiqH2O
  out_file->HOBr = HOBr
  out_file->CH3CCl2F = CH3CCl2F
  out_file->CH3CClF2 = CH3CClF2
  out_file->CHCl2CF3 = CHCl2CF3
  out_file->C2Br2F4 = C2Br2F4
  out_file->CHx = CHx
  out_file->HOx = HOx
  out_file->CBr2F2 = CBr2F2
  out_file->CH2Br2 = CH2Br2
  out_file->CHBr3 = CHBr3
  out_file->Age = Age
  out_file->BrONO = BrONO
  out_file->OzoneDiag = OzoneDiag
  ;
  out_file->Temp = Temp
  out_file->UbarD = UbarD
  out_file->HeatOut = HeatOut
  out_file->CoolOut = CoolOut
  out_file->W = W
  out_file->V = V

end

