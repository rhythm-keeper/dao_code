pro prep_hst_acs_img,input

; -- copied from mkopt by D. Nidever --
; A list file was input

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
if dum eq '' then print,files(i),' DOES NOT EXIST'

thisfile = files(i)
file = strtrim(file_basename(thisfile,'.fits'),2)
; -- EXT0 header
head0 = headfits(thisfile)
for ext_loop=1,4,3 do begin


; write CR points coord to file
bad_data=readfits(thisfile,jnkhead,EXTEN_NO=strtrim(ext_loop+2,2))
bad_points=where(bad_data gt 0) ; 4096 = CR
;openw,cr_lun,file+'_'+strtrim(ext_loop,2)+'_daopam.crxy',/get_lun
;for cr_loop=0,n_cr_points-1 do begin
;	printf,cr_lun,string(array_indices(cr_data,cr_points[cr_loop])+1)
;endfor
;free_lun,cr_lun


; -- SCI header
jnkdata = mrdfits(thisfile,ext_loop,headsci)
data = readfits(thisfile,jnkhead,EXTEN_NO=strtrim(ext_loop,2))

data_corrected_for_flags=data
data_corrected_for_flags[bad_points]=1e19

data_fudged_for_flags=data
; loop through and correct bad points for interpolated values
for bad_point_el=0L,n_elements(bad_points)-1 do begin

	ind=array_indices(data,bad_points[bad_point_el])

	;print,bad_points[bad_point_el],ind
	
	min_range_x = ind[0]-1 > 0
	max_range_x = ind[0]+1 < 4095
	min_range_y = ind[1]-1 > 0
	max_range_y = ind[1]+1 < 2047

	surrounding_array=data_corrected_for_flags[min_range_x:max_range_x,min_range_y:max_range_y]
	
	; do not include central pixel
	remove_bad_pix=where(surrounding_array lt 1e19,num_good)
	surrounding_array=surrounding_array[remove_bad_pix]

	if num_good eq 0 then continue ; cannot correct pixel-no valid neighbors

	data_fudged_for_flags[ind[0],ind[1]]=median(surrounding_array,/even)
	
endfor


; -- 2.1 -- pull easy header parameters

ncomb = sxpar(headsci,'NCOMBINE')
; Keep WCS info
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

skymed = sxpar(headsci,'MDRIZSKY')
totexp = sxpar(head0,'EXPTIME')
sky_cps = skymed/totexp

gainA=sxpar(head0,'ATODGNA')
gainB=sxpar(head0,'ATODGNB')
gainC=sxpar(head0,'ATODGNC')
gainD=sxpar(head0,'ATODGND')
readA=sxpar(head0,'READNSEA')
readB=sxpar(head0,'READNSEB')
readC=sxpar(head0,'READNSEC')
readD=sxpar(head0,'READNSED')

avg_gain_c1=(gainA+gainB)/2.
avg_gain_c2=(gainC+gainD)/2.

avg_read_c1=(readA+readB)/2.
avg_read_c2=(readC+readD)/2.


print,avg_gain_c1,avg_gain_c2,avg_read_c1/avg_gain_c1,avg_read_c2/avg_gain_c2

if ext_loop eq '1' then begin
	;data=data/avg_gain_c2
	data_fudged_for_flags=data_fudged_for_flags/avg_gain_c2
	data_corrected_for_flags=data_corrected_for_flags/avg_gain_c2

endif
if ext_loop eq '4' then begin
	;data=data/avg_gain_c1
	data_fudged_for_flags=data_fudged_for_flags/avg_gain_c1
	data_corrected_for_flags=data_corrected_for_flags/avg_gain_c1
endif

print,files(i);,avg_gain,avg_read

; -- add in sci head details that are important
;sxaddpar,head0,'NCOMBINE',ncomb,'number of image sets combined during CR rejecti'
;sxaddpar,head0,'MDRIZSKY',skymed,'Sky value computed by AstroDrizzle'
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
;sxaddpar,head0,'AVG_GAIN',avg_gain
;sxaddpar,head0,'AVG_READ',avg_read
;sxaddpar,head0,'WCSNAME0','OPUS '
sxaddpar,head0,'WCSNAME',sxpar(headsci,'WCSNAME')
; -- add in calculated sky cps
;sxaddpar,head0,'MSKYCPS',sky_cps,'Sky value in CPS added to image'
; -- add in history statement
;sxaddhist,'Converted from drizzle output into counts with sky for daophot on '+systime(),head0
; -- typecast data to 32 bit since this is what DAOPHOT wants
fits_write,file+'_'+strtrim(ext_loop,2)+'_dao.fits',float(data_corrected_for_flags),head0
fits_write,file+'_'+strtrim(ext_loop,2)+'_dao_noflag.fits',float(data_fudged_for_flags),head0



endfor ; ext_loop
endfor ; nfiles

end
