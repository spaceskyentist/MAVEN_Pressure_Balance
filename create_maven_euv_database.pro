pro create_maven_euv_database, dates=dates, n_days=n_days, start_date=start_date, end_date=end_date, idlsave_pathname=idlsave_pathname

  ; ----------------------------------------------------------
  ; Set up the dates to get data for:
  ; ----------------------------------------------------------
  ; if you want data stored from a particular date, then set it.
  ; Otherwise, get data from first full month of mission data (L2)
  ; https://lasp.colorado.edu/maven/sdc/public/

  if n_elements(dates) eq 0 then begin ; if no array of dates given, then use start and end dates
    if n_elements(start_date) ne 0 then start_date=start_date else start_date='2014-10-15'

    if n_elements(end_date) ne 0 then begin
      end_date=end_date ; use specified end date
    endif else begin
      tmp=strsplit(time_string(systime(/seconds)),'/', /extract) ; current UTC date, seconds since 1970
      end_date = tmp[0] ; current date, without the time
    endelse

    tmp = strsplit(start_date,'-',/extract)
    start_year = uint(tmp[0])
    start_month = uint(tmp[1])
    start_ind = start_month-1 ;start index for month

    tmp = strsplit(end_date,'-',/extract)
    end_year = uint(tmp[0])
    end_month = uint(tmp[1])
    end_ind = end_month-1 ; end index

    n_years = end_year-start_year ; number of years to get data for

    years = indgen(n_years+1)+ start_year  ; set up year array
    years = strtrim(string(years),2) ; turn into strings with no whitespace

    months = indgen(12) + 1 ; set up array of months to go through
    months = strtrim(string(months),2) ; get string of the months
    months[0:8] = '0'+months[0:8] ; add zero to the single digit months for correct formatting


    dates = [start_date] ; initialize the start date
    if n_years eq 0 then begin
      if start_ind-end_ind eq 0 then goto, ltonemonth ; if you want less than a month of data
      dates=[dates, years[0]+'-'+months[start_ind+1:end_ind]+'-'+'01'] ; for only a year of data
    endif else begin
      if start_month ne 12 then dates=[dates, years[0]+'-'+months[start_ind+1:*]+'-'+'01'] ; if the starting month is in december, continue to next year
      for i = 1,n_years do dates=[dates, years[i]+'-'+months[*]+'-'+'01']
      tmp = where(dates eq strtrim(string(end_year)+'-'+months[end_ind]+'-'+'01',2))
      dates = dates[0:tmp] ; remove any months after the end date month
    endelse

    ltonemonth: if end_date ne dates[-1] then dates = [dates, end_date] ; add the end date
  endif



  ; ----------------------------------------------------------
  ; Download the Data for the Specified Dates
  ; ----------------------------------------------------------

  n_dates = n_elements(dates)

  setup_colortable ; use the colortable provided by Scott Thaller


  ; get and save data onto external harddrive
  for i=0.d,n_dates-2 do begin

    if n_elements(n_days) ne 0 then tr = [dates[i], time_string(time_double(dates[i])+n_days*86400.D,prec=-3)] $
    else tr = [dates[i], dates[i+1]] ; range of dates to get data for

    print, ' '
    print, 'Retrieving data for ', tr[0], ' through ', tr[1]
    print, ' '
    wait, 3

    timespan, tr
    
    mvn_euv_l3_load ;get FISM model
    get_data, 'mvn_euv_l3_minute',data=data ; data.v wavelength
    if size(data,/type) eq 0 then begin
      print, 'issue with file ', tr[0], ' through ', tr[1]
    endif
    ieuv = where(data.v gt 0 and data.v lt 91) ; EUV ionizing wavelengths
    z = total(data.y[*,ieuv],2)
    store_Data, 'mvn_euv_l3_euvirr', data={x:data.x, y:z},dlimit={ytitle:'EUV irr (0-91 nm)', ysubtitle:'[W/m2]'}
    
    
  
    

    ; -------------------------------------------------------------
    ; IDL Save Variables
    ; -------------------------------------------------------------


  date_str = strsplit(time_string(data.x[0]), '/', /extract)
  date_str = date_str[0]
  end_date_str =strsplit(time_string(data.x[-1]), '/', /extract)
  end_date_str = end_date_str[0]
  save_fn = 'database_vars_'+date_str+'_to_'+end_date_str+'.sav'
    
    
    ;break up time
    time = data.x
    tmp = strsplit(time_string(time),'/', /extract)
    yr_mon_day = fltarr(n_elements(time), 3)
    tmp3 = strarr(n_elements(time))
  
    for tt=0,n_elements(tmp)-1 do begin
      yr_mon_day[tt,*] = float(strsplit(tmp[tt,0],'-',/extract))
      tmp3[tt] = tmp[tt,0]
    endfor
  
    sec_of_day = time - time_double(tmp3)
  
    database_array = [ [time], [yr_mon_day], [sec_of_day], [z] ]
  
    ;save_fn = 'database_vars_'+date_str+'_to_'+end_date_str+'.sav'
    save, database_array, filename=idlsave_pathname+save_fn
  
    print, 'Saved file:', idlsave_pathname+save_fn
    delvar, date_str,time, tmp, yr_mon_day, tmp3, sec_of_day, database_array
    ;endif
  
    store_data,'*',/delete ; delete the tplot variables from IDL memory
  
    nextdatafile:
    print, 'Continuing to next date range file...'
    
  endfor ; over dates

print, 'Completed data retrieval.'


end