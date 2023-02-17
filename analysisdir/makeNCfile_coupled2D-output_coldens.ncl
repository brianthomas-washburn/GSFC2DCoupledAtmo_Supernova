
;This code takes column density output from GSFC 2D coupled model,
;coldens(time,lat)
;
; This code assumes input file has 2 years of daily data (with 360 day years)
;
;Structure of output from model:
;        DATA INSP/3,4,5,6,10,31,59/
;! 3,  4,  5,  6,   10,   31,  59
;! O2, O3, NO, NO2, HNO3, NOy, H2O precip
;
;Also reads in rainout of HNO3, which was output by GSFC model
; to a different file (but, same structure).
;
;write out in model at each time-write is:
;         do ij=1,L$
;            write(5131,'(1x,5e16.8)') 
;     >            (ncolgd(is,ij,1),is=1,nsp)
;---------------

;The code here reads in the text file, divides the data
;  into individual variables (separate ones for each constituent),
;  and then creates a netCDF file with this data and associated attributes.
;The output file is more usable (contains metadata, etc.), and smaller.
;
;BCThomas, June 2018

begin
  
;  in_dir      = "/storage/atmoRuns/test11June2018"
;  in_file1    = "BToutput_coldens_daily_test11June2018.out"
;  in_file2    = "BToutput_coldenhno3_daily_test11June2018.out"

  files = new((/3/),string)
  files = asciiread("nclInput_colden.txt",(/3,1/),"string")
  in_dir = files(0)
  in_file1 = files(1)
  in_file2 = files(2)

  print(in_dir+" , "+in_file1+" , "+in_file2)
  

;some parameter values:                                
  timemax = 720   ;number of time writes
  latmax  = 45    ;number of latitude bins
  chemmax = 7     ;number of chemical constituents output

  ncols = chemmax
  nrows = timemax*latmax

  in_data1    = asciiread(in_dir+"/"+in_file1,(/nrows,ncols/),"float")	;array for reading in data
  ;printVarSummary(in_data1)
  ;print(in_data1)

  ;Now I want to structure data in an array that is (time,lat,data)
  data1 = new((/timemax,latmax,ncols/),float) ;make new array
  data1@_FillValue = -9999.99 ;set fill value

  ;Loop to fill data array from input
  nt=0
  nl=0
  do i=0,(timemax*latmax)-1 ;loops over all the in_data values
    do s=0,chemmax-1 
      data1(nt,nl,s) = in_data1(i,s)
    end do
;    print(i+" , "+nt+" , "+nl+" , "+in_data1(i,1)+" , "+data1(nt,nl,1))
    if(nl.le.(latmax-1)) then ;nl steps over latitude - if we're less than latmax, increment it
      nl=nl+1
    end if
    if(nl.gt.(latmax-1)) then ;but if we reach latmax, then it's time to start at 0 again
      nl=0
      nt=nt+1 ;and, it's time to go to the next time index in data
    end if
  end do
  ;---

;Need to set up latitude array:
  lat4arr = new((/latmax/),float)
  lat4arr = ispan(-88,88,4)
  lat4arr@long_name = "latitude"
  lat4arr@units     = "degrees_north"
  lat4arr@_FillValue = -9999.99
  lat4arr@missing_value = -9999.99
  lat4arr!0 = "lat"
  ;print(lat4arr)
  ;printVarSummary(lat4arr)  

;Define variables for the chemical constituents:
  O2 = new((/timemax,latmax/),float)
  O2 = data1(:,:,0)
  O2@long_name = "O2 column density"
  O2@units     = "molecules cm^-2"
  O2@_FillValue = -9999.99
  O2@missing_value = -9999.99
  O2!0 = "time"
  O2!1 = "lat"
  O2&lat = lat4arr
  ;printVarSummary(O2)
  ;print(O2)

  O3 = new((/timemax,latmax/),float)
  O3 = data1(:,:,1)
  O3@long_name = "O3 column density"
  O3@units     = "molecules cm^-2"
  O3@_FillValue = -9999.99
  O3@missing_value = -9999.99
  O3!0 = "time"
  O3!1 = "lat"
  O3&lat = lat4arr
  ;printVarSummary(O3)                                                                            
  ;print(O3) 

  NO = new((/timemax,latmax/),float)
  NO = data1(:,:,2)
  NO@long_name = "NO column density"
  NO@units     = "molecules cm^-2"
  NO@_FillValue = -9999.99
  NO@missing_value = -9999.99
  NO!0 = "time"
  NO!1 = "lat"
  NO&lat = lat4arr
  ;printVarSummary(NO)
  ;print(NO) 

  NO2 = new((/timemax,latmax/),float)
  NO2 = data1(:,:,3)
  NO2@long_name = "NO2 column density"
  NO2@units     = "molecules cm^-2"
  NO2@_FillValue = -9999.99
  NO2@missing_value = -9999.99
  NO2!0 = "time"
  NO2!1 = "lat"
  NO2&lat = lat4arr
  ;printVarSummary(NO2)

  HNO3 = new((/timemax,latmax/),float)
  HNO3 = data1(:,:,4)
  HNO3@long_name = "HNO3 column density"
  HNO3@units     = "molecules cm^-2"
  HNO3@_FillValue = -9999.99
  HNO3@missing_value = -9999.99
  HNO3!0 = "time"
  HNO3!1 = "lat"
  HNO3&lat = lat4arr
  ;printVarSummary(HNO3) 

  NOy = new((/timemax,latmax/),float)
  NOy = data1(:,:,5)
  NOy@long_name = "NOy column density"
  NOy@units     = "molecules cm^-2"
  NOy@_FillValue = -9999.99
  NOy@missing_value = -9999.99
  NOy!0 = "time"
  NOy!1 = "lat"
  NOy&lat = lat4arr
  ;printVarSummary(NOy) 

  H2O_rainout = new((/timemax,latmax/),float)
  H2O_rainout = data1(:,:,0)
  H2O_rainout@long_name = "H2O rainout"
  H2O_rainout@units     = "molecules cm^-2 s^-1"
  H2O_rainout@_FillValue = -9999.99
  H2O_rainout@missing_value = -9999.99
  H2O_rainout!0 = "time"
  H2O_rainout!1 = "lat"
  H2O_rainout&lat = lat4arr
  ;printVarSummary(H2O_rainout) 

;------ Now, also add in HNO3 rainout, which is output from model in a separate file:

  ncols := 1
  nrows := timemax*latmax

  in_data2    = asciiread(in_dir+"/"+in_file2,(/nrows,ncols/),"float")  ;array for reading in data 
  ;printVarSummary(in_data2)

  ;Now I want to structure data in an array that is (time,lat,data)                              

  data2 = new((/timemax,latmax,ncols/),float) ;make new array                                    
  data2@_FillValue = -9999.99 ;set fill value                                                   

  ;Loop to fill data array from input                                                             
  nt=0
  nl=0
  do i=0,(timemax*latmax)-1 ;loops over all the in_data values                        
    data2(nt,nl,0) = in_data2(i,0)
    ;print(i+" , "+nt+" , "+nl+" , "+in_data1(i,6)+" , "+data1(nt,nl,6))                          
    if(nl.le.(latmax-1)) then ;nl steps over latitude - if we're less than latmax, increment it   
      nl=nl+1
    end if
    if(nl.gt.(latmax-1)) then ;but if we reach latmax, then it's time to start at 0 again          
      nl=0
      nt=nt+1 ;and, it's time to go to the next time index in data                             
    end if
  end do
  ;---  
  HNO3_rainout = new((/timemax,latmax/),float)
  HNO3_rainout = data2(:,:,0)
  HNO3_rainout@long_name = "HNO3 rainout"
  HNO3_rainout@units     = "molecules cm^-2 s^-1"
  HNO3_rainout@_FillValue = -9999.99
  HNO3_rainout@missing_value = -9999.99
  HNO3_rainout!0 = "time"
  HNO3_rainout!1 = "lat"
  HNO3_rainout&lat = lat4arr
  ;printVarSummary(HNO3_rainout)      

; This block creates a NetCDF file out of the data that was read in.
;------------------------------------------------------------
;  out_dir = "./"
  out_dir = in_dir
;  fn = out_dir+"_"+in_file1+".nc"
  fn = in_file1+".nc"
  system("rm -f "+out_dir+"/"+fn) ;remove old version

  out_file    = addfile(out_dir+"/"+fn,"c")
;Global attributes:
  out_file@title	= "Output of GSFC 2D atmo model (coupled version); column densities. brian.thomas@washburn.edu" ;set title of file
  out_file@creation_date = systemfunc("date")      ;set creation date from system

;output data values:
  out_file->lat4arr = lat4arr
  out_file->O2 = O2
  out_file->O3 = O3
  out_file->NO = NO
  out_file->NO2 = NO2
  out_file->HNO3 = HNO3
  out_file->NOy = NOy
  out_file->H2O_rainout = H2O_rainout
  out_file->HNO3_rainout = HNO3_rainout

end
