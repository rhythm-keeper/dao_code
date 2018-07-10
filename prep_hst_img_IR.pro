pro prep_hst_img_IR,input
; --
; -- short script like code to convert input hst drizzle products into their proper units for daophot --
; 
; -- v0.0 -- RLBeaton, 28 Jan 2015
; -- v0.1 -- RLBeaton, 03 Feb 2015
;       * resolved issue with data formatting such that the image can go directly into daophot from this script
;       * putting information from the SCI extension header into the output header file for daophot so no information is lost
;       * adding a history comment with date/time
; -- v1.0 -- Hatt, 13 May 2015
;	* this code accomodates both the drizzled and individual images
;	* added parameters to output header to keep track of WCS info
;	* drizzled images are adjusted for plate scale
;	* "hot" or NaN pixels from the drizzled images are set to 1e10 so they are ignored in DAOPHOT
; -- 

; -- copied from mkopt by D. Nidever -- 
 ; A list file was input
print,input

 if strmid(input[0],0,1) eq '@' then begin
  inp = strmid(input[0],1)

   ; Loading the files
   readcol,inp,files,format='A',/silent
   nfiles = n_elements(files)
  endif else begin

   ; Probably an array of filenames
   if n_elements(input) gt 1 then begin
     files = input
     nfiles = n_elements(files)
      ; A globbed list
   endif else begin
     files = file_search(input)
     nfiles = n_elements(files)
   endelse
 endelse

  ; -- start loop
 for i=0, nfiles-1 do begin

  ; -- 0 -- test for file
   dum = file_search(files(i))
  if dum eq '' then begin
   print,files(i),' DOES NOT EXIST'
  end

   ; -- 
  thisfile = files(i)
   ; -- strip ".fits"
  file = strtrim(file_basename(thisfile,'.fits'),2)

  ; -- 1 -- read in header in ext0, 
  ;         read in image in SCI,
  ;         read in image header in SCI

   ; -- EXT0 header
  head0 = headfits(thisfile)

   ; -- SCI header
  jnkdata  = mrdfits(thisfile,'SCI',headsci)
  data = readfits(thisfile,jnkhead,EXTEN_NO='1')
  data_corrected_for_flags=data

  flagged_data=readfits(thisfile,jnkhead,EXTEN_NO='3')
  flagged_data_ind=where(flagged_data gt 0)	

 data_corrected_for_flags[flagged_data_ind]=1e19

 data_fudged_for_flags=data
 ; loop through and correct bad points for interpolated values
 for bad_point_el=0L,n_elements(flagged_data_ind)-1 do begin
 
         ind=array_indices(data,flagged_data_ind[bad_point_el])
 
 
         min_range_x = ind[0]-1 > 0
         max_range_x = ind[0]+1 < 1014-1
         min_range_y = ind[1]-1 > 0
         max_range_y = ind[1]+1 < 1014-1
 
         surrounding_array=data_corrected_for_flags[min_range_x:max_range_x,min_range_y:max_range_y]
 
         ; do not include central pixel
         remove_bad_pix=where(surrounding_array lt 1e19,num_good)
         surrounding_array=surrounding_array[remove_bad_pix]
 
         if num_good eq 0 then continue ; cannot correct pixel-no valid neighbors
 
         data_fudged_for_flags[ind[0],ind[1]]=median(surrounding_array,/even)
 
 endfor



  ; HATT MODIFICATION: GET RID OF NaN VALUES
  ind=where(finite(data) eq 0,count)
  print,'# NaN values='+strtrim(count,2)
  data[ind]=1e19

  ; -- 2.1 -- pull easy header parameters 
  ncomb  = sxpar(headsci,'NCOMBINE')
  skymed = sxpar(headsci,'MDRIZSKY')
  totexp = sxpar(head0,'EXPTIME')

  ; HATT MODIFICATION: KEEP WCS INFO
  CRPIX1 = sxpar(headsci,'CRPIX1')
  CRPIX2 = sxpar(headsci,'CRPIX2')
  CRVAL1 = sxpar(headsci,'CRVAL1')
  CRVAL2 = sxpar(headsci,'CRVAL2')
  CTYPE1 = sxpar(headsci,'CTYPE1')
  CTYPE2 = sxpar(headsci,'CTYPE2')
  CD1_1 = sxpar(headsci,'CD1_1')
  CD1_2 = sxpar(headsci,'CD1_2')
  CD2_1 = sxpar(headsci,'CD2_1')
  CD2_2 = sxpar(headsci,'CD2_2')
  CUNIT1 = sxpar(headsci,'CUNIT10')
  CUNIT2 = sxpar(headsci,'CUNIT20')

  ; -- 2.2 -- pull exptimes 
  ; HATT MODIFICATION: IF DRIZZLED IMAGES
  ; USE AVERAGE EXPOSURE TIME OTHERWISE MEAXEXP==TOTEXP 
  ;if files[i] eq 'f1_h_drz.fits' or files[i] eq 'f2_h_drz.fits' then begin
  ;exps = fltarr(ncomb)
  ;for j = 0, ncomb-1 do begin
  ; no = string(j,format='(i03)')
  ; tempkey = 'D'+no+'DEXP'
  ; exps[j] = sxpar(head0,tempkey)
  ;endfor
  ; meanexp = mean(exps)
  ;endif else begin
   meanexp=totexp ; for individual frame
  ;endelse

  ; HATT MODIFICATION: USE MEDIAN DRIZSKY FOR FLT FILES
  ; SINCE DRZ DOES NOT KEEP THAT INFO
  ;if files[i] eq 'f1_h_drz.fits' or files[i] eq 'f2_h_drz.fits' then begin
  ; sky_values=[]
  ; if files[i] eq 'f1_h_drz.fits' then flt_list='f1_h.list'
  ; if files[i] eq 'f2_h_drz.fits' then flt_list='f2_h.list'
  ; readcol,flt_list,flt_file_names,format='(A)'  
  ; for flt_loop=0,n_elements(flt_file_names)-1 do begin
  ; 	jnkdata  = mrdfits(flt_file_names[flt_loop],'SCI',flt_headsci,/silent)
  ;       sky_values=[sky_values,sxpar(flt_headsci,'MDRIZSKY')]
  ; endfor
  ; skymed=median(sky_values)
  ;endif

  ; -- 3 -- image manipulations
  ;sky_cps = skymed/meanexp
  
  ; HATT MODIFICATION: SCALE SKY BACKGROUND BY PIXEL AREA
  ; scale the sky background by the pixel area
  ; ACS/WFC
  ;new_scale=0.03333
  ;plate_scale=0.05
  ; WFC3/IR
  ;if files[i] eq 'f1_h_drz.fits' or files[i] eq 'f2_h_drz.fits' then begin
  ; new_scale=0.0666
  ; plate_scale=0.13
  ; sky_cps=sky_cps*(new_scale/plate_scale)^2.
  ;endif
  ;print,'Sky for ',thisfile,' in cps ',sky_cps
  ;print,'Total Exposure for ',thisfile,' in seconds ',totexp
  ;data = data + sky_cps
  ;data = data * totexp
  data_corrected_for_flags = data_corrected_for_flags * totexp
  data_fudged_for_flags = data_fudged_for_flags * totexp
  


gainA=sxpar(head0,'ATODGNA')
gainB=sxpar(head0,'ATODGNB')
gainC=sxpar(head0,'ATODGNC')
gainD=sxpar(head0,'ATODGND')
readA=sxpar(head0,'READNSEA')
readB=sxpar(head0,'READNSEB')
readC=sxpar(head0,'READNSEC')
readD=sxpar(head0,'READNSED')

avg_gain=avg([gainA,gainB,gainC,gainD])
avg_readnoise=avg([readA,readB,readC,readD])

print,'avg gain, read noise',avg_gain,avg_readnoise


;data = data/avg_gain
data_corrected_for_flags = data_corrected_for_flags/avg_gain
data_fudged_for_flags = data_fudged_for_flags/avg_gain





  ; -- 4 -- write out
  ;mwrfits,data,file+'_dao.fits',head0,/create,/ascii
   ; -- reset header parameter EXTEND to 'false' because there are no extensions
  sxaddpar,head0,'EXTEND','F'
   ; -- add in sci head details that are important
  sxaddpar,head0,'NCOMBINE',ncomb,'number of image sets combined during CR rejecti'
  sxaddpar,head0,'MDRIZSKY',skymed,'Sky value computed by AstroDrizzle'
  ; HATT MODIFICATION: WRITE WCS INFO TO HEADER
  sxaddpar,head0,'CRPIX1',CRPIX1
  sxaddpar,head0,'CRPIX2',CRPIX2
  sxaddpar,head0,'CRVAL1',CRVAL1
  sxaddpar,head0,'CRVAL2',CRVAL2
  sxaddpar,head0,'CTYPE1',CTYPE1
  sxaddpar,head0,'CTYPE2',CTYPE2
  sxaddpar,head0,'CD1_1',CD1_1
  sxaddpar,head0,'CD1_2',CD1_2
  sxaddpar,head0,'CD2_1',CD2_1
  sxaddpar,head0,'CD2_2',CD2_2
  sxaddpar,head0,'CUNIT10',CUNIT1
  sxaddpar,head0,'CUNIT20',CUNIT2
  sxaddpar,head0,'WCSNAME0','OPUS    '
  sxaddpar,head0,'WCSNAME',sxpar(headsci,'WCSNAME')
  ; -- add in calculated sky cps
  ;sxaddpar,head0,'MSKYCPS',sky_cps,'Sky value in CPS added to image'
   ; -- add in history statement
  ;sxaddhist,'Converted from drizzle output into counts with sky for daophot on '+systime(),head0
  ; -- typecast data to 32 bit since this is what DAOPHOT wants
  fits_write,file+'_dao.fits',float(data_corrected_for_flags),head0
  fits_write,file+'_dao_noflag.fits',float(data_fudged_for_flags),head0
 endfor


end
