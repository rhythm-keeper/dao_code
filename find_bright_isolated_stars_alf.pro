pro find_bright_isolated_stars_alf 


; loop through alf files

alf_files=file_search('*daopam.alf',count=n_alf_files)


max_mag=16

; read in alf files, find matching lst file
for curr_alf_file=0,n_alf_files-1 do begin

	readcol,alf_files[curr_alf_file],id_alf,x_alf,y_alf,mag_alf,magerr_alf,sky_alf,nit_alf,chi_alf,sharp_alf,skipline=3,/silent

	alf_filename_base=file_basename(alf_files[curr_alf_file],'.alf')

	all_x_alf=x_alf
	all_y_alf=y_alf
	all_mag_alf=mag_alf

	ind=where(mag_alf le max_mag and $
		  magerr_alf lt 0.09+0.015*exp(mag_alf-16.75)  and $
		  sharp_alf lt 0.2+0.015*exp(mag_alf-16.75) and $
        	  sharp_alf gt -0.2-0.015*exp(mag_alf-16.75) )			


			

	id_alf=id_alf[ind]
	x_alf=x_alf[ind]
	y_alf=y_alf[ind]
	mag_alf=mag_alf[ind]

	; make sure bright stars are not neighbors
	;for star=0,n_elements(id_alf)-1 do begin

	;	
	;	ind=where(all_mag_alf le max_mag+1)

	;	 dx=x_alf[star]-all_x_alf[ind] 
        ;         dy=y_alf[star]-all_y_alf[ind]
        ;         distances=sqrt( dx^2. + dy^2. )
	;	 ;ind=where(distances lt 10 and distances ne 0,num)

	;	 ;if num gt 0 then mag_alf[star]=99

	;	; look for stars up to 1 mag fainter than bright star and reject
	;	

	;endfor
	; remove bright stars that are too close to other bright stars
	; this should be rare
	ind=where(mag_alf le max_mag)
        id_alf=id_alf[ind]
        x_alf=x_alf[ind]
        y_alf=y_alf[ind]
        mag_alf=mag_alf[ind]	


	query='head -3 '+alf_files[curr_alf_file]+' > '+alf_filename_base+'.lst'
        spawn,query

	 foreach curr_id,id_alf do begin

                        query='grep '+string(39B)+' '+strtrim(string(curr_id,format='(I)'),2)+' '+string(39B)+' '+alf_files[curr_alf_file]
                        spawn,query,output

                        query='echo '+string(39B)+output+string(39B)+' >> '+alf_filename_base+'.lst'

                        spawn,query

	endforeach


	

endfor ; alf file loop

end
