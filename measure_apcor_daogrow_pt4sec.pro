pro measure_apcor_daogrow_pt4sec

; get parameters from .gro file daogrow

readcol,'name.list',name,format='(A)'

openw,lunout,'all_meas_apcor.dat',/get_lun


for i=0,n_elements(name)-1 do begin


	apcor_meas=[]
        apcor_meas_err=[]
        alf_mags=[]



	; using one of the sets

	; read in alf file for matching purposes
	readcol,name[i]+'.alf',id_alf,x_alf,y_alf,mag_alf,magerr_alf,skipline=3,/silent


	; read in aperture file

	; read in aperture file
        query='wc -l < '+name[i]+'a.ap'
        spawn,query,nlines
        if float(nlines) le 3 then begin
                printf,lunout,name[i]+' -99 -99 -99 -99 -99 -99'
                continue
        endif



	openr,lun,name[i]+'a.ap',/get_lun

	readf,lun,line ; header
	readf,lun,line ; header
		

	while ~eof(lun) do begin

		readf,lun,id,x,y,ap1,ap2,ap3,ap4,ap5,ap6,ap7,ap8,ap9,ap10,ap11,ap12	
		readf,lun,sky,chi,sharp,aperr1,aperr2,aperr3,aperr4,aperr5,aperr6,aperr7,aperr8,aperr9,aperr10,aperr11,aperr12 


		; get min aperture error of largest ap
		all_aps   =[ap1,ap2,ap3,ap4,ap5,ap6,ap7,ap8,ap9,ap10,ap11,ap12]
		all_aperrs=[aperr1,aperr2,aperr3,aperr4,aperr5,aperr6,aperr7,aperr8,aperr9,aperr10,aperr11,aperr12]

		;apind=where(all_aperrs eq min(all_aperrs))
		;best_ap_index=min(apind)

		; ignore stars with bad ap phot
	
		bad_ap=where( all_aps[4] ge 99. or all_aps[4] lt 0 , num_bad_ap)
		if num_bad_ap ge 1 then continue
		bad_appers=where( all_aperrs[4] ge 9. , num_bad_aperr)
		if num_bad_aperr ge 1 then continue
		
		;for offset_loop=10,4,-1 do begin
                ;        all_aps[11] -= growth_offsets[offset_loop]
                ;endfor

		distances=sqrt( (x-x_alf)^2. + (y-y_alf)^2.  )
                distind=where( distances eq min(distances) )
		
		print,min(distances),all_aps[4],mag_alf[distind]
		apcor_meas=[apcor_meas, all_aps[4]-mag_alf[distind] ]
		apcor_meas_err=[apcor_meas_err,sqrt(all_aperrs[4]^2.+magerr_alf[distind]^2.)]

		alf_mags=[alf_mags,mag_alf[distind]]


	endwhile
	ind=where(apcor_meas gt -0.25 and apcor_meas lt 0.25,num)
	;print,name[i],apcor_meas[ind]


	if num gt 0 then print,name[i],median(apcor_meas[ind],/even),wmean(apcor_meas[ind],apcor_meas_err[ind])

	   if n_elements(apcor_meas) gt 0 then begin

                        ;print,name[i],median(apcor_meas,/even),mad(apcor_meas)/sqrt(n_elements(apcor_meas)),wmean(apcor_meas,apcor_meas_err),stdev(apcor_meas)

                        for el=0,n_elements(apcor_meas)-1 do begin
                                printf,lunout,name[i]+' '+strtrim(apcor_meas[el],2)+' '+strtrim(apcor_meas_err[el],2)+' '+strtrim(alf_mags[el],2)
                        endfor
                        ;printf,lunout,name[i]+' '+strtrim(median(apcor_meas,/even),2)+' '+strtrim(mad(apcor_meas)/sqrt(n_elements(apcor_meas)),2)
                endif


	;print,name[i],median(apcor_meas,/even),mad(apcor_meas)/sqrt(n_elements(apcor_meas))

	;printf,lunout,name[i],median(apcor_meas,/even),mad(apcor_meas)/sqrt(n_elements(apcor_meas))

	free_lun,lun


endfor

free_lun,lunout


end
