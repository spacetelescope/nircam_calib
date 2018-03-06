pro find_centroids

  files = file_search('mag5/nrca5_TA_timeseries_NRCFLATA5GRTS*uncal.fits')
  for mm = 6,10 do begin
     mfiles = file_search('mag'+strcompress(string(mm),/remove_all)+'/nrca5_TA_timeseries_NRCFLATA5GRTS*uncal.fits')
     files = [files,mfiles]
  endfor
    
  
  ;sort the input files by true peak location, and
  ;then magnitude
  
  ;extract the center locations of all files,
  centerstrings = strarr(n_elements(files))
  for i = 0,n_elements(files)-1 do begin
     cstr1 = strpos(files[i],'center_')
     cstr2 = strpos(files[i],'noise')
     slen = (cstr2-1) - (cstr1+7)
     fname = files[i]
     centerstrings[i] = strmid(fname,cstr1+7,slen)
  endfor
  print,'centerstrings',centerstrings
  
  ;find unique center values
  centers = centerstrings[uniq(centerstrings, sort(centerstrings))]
  centers = centers[sort(centers)]
  print,'centers',centers
  
  ;sort the files with each center value
  allfiles = ['xx']
  for i = 0,n_elements(centers)-1 do begin
     print,'xxx',centers[i]
     centerfiles = files[WHERE(STRMATCH(files, '*'+centers[i]+'*', /FOLD_CASE) EQ 1)]
     print,'1',centerfiles
     centerfiles = centerfiles[sort(centerfiles)]
     print,'2',centerfiles
     allfiles = [allfiles,centerfiles]
  endfor
  allfiles = allfiles[1:*]
  files = allfiles

  ;half width of centroid box. To match the flight software
  ;set to 4
  hw = 4
  ;half width of checkbox. To march flight software, set to 1
  checkboxhw = 1
  
  openw,lun,'centroid_table_full_dataset.tab',/get_lun,width=200
  printf,lun,'File    Noise_Realization   True_Centroid_x     True_Centroid_y     Measured_Centroid_x    Measured_Centroid_y     Initial_Centoid_x    Initial_Centroid_y   Centoid_x_Error   Centroid_y_Error    Centroid_Radial_Error_Pix    Centroid_Radial_Error_Arcsec   Source_Magnitude_in_Simulator'
  
  for i = 0,n_elements(files)-1 do begin
     fits_read,files[i],image
     im1 = image[*,*,0] * 1.
     im2 = image[*,*,1] * 1.
     im3 = image[*,*,2] * 1.
     d32 = im3 - im2
     d21 = im2 - im1
     ta_image = im3 * 0.
     for x = 0,n_elements(im3[*,0])-1 do begin
        for y = 0,n_elements(im3[0,*])-1 do begin
           ta_image[x,y] = min([d32[x,y],d21[x,y]])
        endfor
     endfor

     slash = strpos(files[i],'/')
     fileandpath = files[i]
     dirname = strmid(fileandpath,0,slash+1)
     fname = strmid(fileandpath,slash+1,strlen(fileandpath)-(slash+1))
     taname = dirname + 'ta_image_for_' + fname
     fits_write,taname,ta_image

     init_center = indgen(2)
     peak = mycentr3(ta_image,hw,checkboxhw,10,init_center)

     ; Get the real centroid location
     cstr1 = strpos(files[i],'center_')
     cstr2 = strpos(files[i],'uncal.fits')
     slen = (cstr2-1) - (cstr1+7)
     cstring = strmid(files[i],cstr1+7,slen)
     under = strpos(cstring,'_')
     centerx = float(strmid(cstring,0,under))
     centery = float(strmid(cstring,under+1,strlen(cstring)-(under+1)))
     noisepos = strpos(files[i],'noiserealiz_')
     noise = strmid(fileandpath,noisepos+12,1)
     
     realx = centerx - 0.5
     ;if centerx eq 16 then realx = centerx + 0.5
     realy = centery - 0.5
     if centery eq 16 then realy = centery + 0.5
     if centery eq 15 then realy = centery + 0.5

     deltax = realx - peak[0]
     deltay = realy - peak[1]
     deltar = sqrt(abs(deltax)^2 + abs(deltay)^2)

     pixscale = 0.063           ;arcsec/pixel for A5
     deltararcsec = deltar * pixscale

     magloc = strpos(files[i],'_mag')
     under = strpos(files[i],'_',magloc+1)
     magstr = float(strmid(files[i],magloc+4,under-(magloc+4)))
     
     line = files[i]+'  '+string(noise)+'  '+string(realx)+'  '+string(realy)+'  '+string(peak[0])+'  '+string(peak[1])+'  '+string(init_center[0])+'  '+string(init_center[1])+'  '+string(deltax)+'  '+string(deltay)+'  '+string(deltar)+'  '+string(deltararcsec)+'  '+string(magstr)
     print,files[i],peak,init_center
     printf,lun,line

  endfor

  close,lun
  free_lun,lun
  
end
