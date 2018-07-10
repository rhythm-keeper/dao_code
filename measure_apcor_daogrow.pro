pro measure_apcor_daogrow


for ext_loop=1,4,3 do begin

	ext=strtrim(ext_loop,2)
	
	readcol,'all_'+ext+'.list',name,format='(A)'
	
	
	openw,lunout,'all_'+ext+'_meas_apcor.dat',/get_lun
	
	
	for i=0,n_elements(name)-1 do begin
	
	
		apcor_x=[]
		apcor_y=[]
		apcor_meas=[]
		apcor_meas_err=[]
		ap_mags=[]
		alf_mags=[]
	
		; set incremental growth dmags
		;readcol,name[i]+'a.cur',ap,set1,set2,set3,set4,/silent
	
		;if n_elements(ap) eq 0 then continue
	
		;growth_offsets=set2
	
		; using one of the sets
	
		; read in alf file for matching purposes
		readcol,name[i]+'.alf',id_alf,x_alf,y_alf,mag_alf,magerr_alf,skipline=3,/silent
	
		openr,lun,name[i]+'a.ap',/get_lun


		; read in aperture file
                query='wc -l < '+name[i]+'a.ap'
                spawn,query,nlines
                if float(nlines) le 3 then begin
			printf,lunout,name[i]+' -99 -99 -99 -99 -99 -99'	
			continue
		endif
	
		readf,lun,line ; header
		readf,lun,line ; header
			
	
		while ~eof(lun) do begin
	
			readf,lun,id,x,y,ap1,ap2,ap3,ap4,ap5,ap6,ap7,ap8,ap9,ap10,ap11,ap12	
			readf,lun,sky,chi,sharp,aperr1,aperr2,aperr3,aperr4,aperr5,aperr6,aperr7,aperr8,aperr9,aperr10,aperr11,aperr12 
	
	
			; get min aperture error of largest ap
			all_aps   =[ap1,ap2,ap3,ap4,ap5,ap6,ap7,ap8,ap9,ap10,ap11,ap12]
			all_aperrs=[aperr1,aperr2,aperr3,aperr4,aperr5,aperr6,aperr7,aperr8,aperr9,aperr10,aperr11,aperr12]
	
			;apind=where(all_aperrs eq min(all_aperrs))
			;best_ap_index=max(apind)
	
			best_ap_index=11
	
			;if best_ap_index ne 11 then begin ; largest aperture not the smallest uncertainty, so correct to largest ap
	
			;	; ignore stars with bad ap phot
	
		
			;	bad_ap=where( all_aps[best_ap_index:10] ge 99. or all_aps[best_ap_index:10] lt 0 , num_bad_ap)
			;	if num_bad_ap ge 1 then continue
			;	bad_appers=where( all_aperrs[best_ap_index:10] ge 9. , num_bad_aperr)
			;	if num_bad_aperr ge 1 then continue
			;	
	
			;	for offset_loop=best_ap_index,10 do begin
	                ;	        all_aps[best_ap_index] += growth_offsets[offset_loop]
	                ;	endfor
	
			;endif
	
	
			;print,name[i],all_aps[best_ap_index],all_aperrs[best_ap_index]	
			bad_ap=where( all_aps[best_ap_index] ge 99. or all_aps[best_ap_index] lt 0 , num_bad_ap)
			if num_bad_ap ge 1 then continue
			bad_appers=where( all_aperrs[best_ap_index] ge 9. , num_bad_aperr)
			if num_bad_aperr ge 1 then continue
	
	
			; IGNORE BRIGHT STARS CLOSEST TO GALAXY
			;if x lt 250 then continue 
	
	
			distances=sqrt( (x-x_alf)^2. + (y-y_alf)^2.  )
	
	
			min_dist=min(distances)
	
	                distind=where( distances eq min_dist )
			
			apcor_x=[apcor_x,x_alf[distind]]
			apcor_y=[apcor_y,y_alf[distind]]		
			apcor_meas=[apcor_meas, all_aps[best_ap_index]-mag_alf[distind] ]
			;print,all_aps[best_ap_index],mag_alf[distind]
			apcor_meas_err=[apcor_meas_err,sqrt(all_aperrs[best_ap_index]^2.+magerr_alf[distind]^2.)]
	
			ap_mags=[ap_mags,all_aps[best_ap_index]]
			alf_mags=[alf_mags,mag_alf[distind]]
	
		endwhile

		free_lun,lun

		;ind=where(apcor_meas gt -0.15 and apcor_meas lt 0.15)
		;print,name[i],apcor_meas[ind]
		;print,name[i],apcor_meas
		
		
		if n_elements(apcor_meas) gt 0 then begin 
	
			;print,name[i],median(apcor_meas,/even),mad(apcor_meas)/sqrt(n_elements(apcor_meas)),wmean(apcor_meas,apcor_meas_err),stdev(apcor_meas)
	
			;for el=0,n_elements(apcor_meas)-1 do begin
			;	printf,lunout,name[i]+' '+$
			;		strtrim(apcor_x[el],2)+' '+$
			;		strtrim(apcor_y[el],2)+' '+$
			;		strtrim(apcor_meas[el],2)+' '+strtrim(apcor_meas_err[el],2)+' '+$
			;		strtrim(ap_mags[el],2)+' '+$
			;		strtrim(alf_mags[el],2)
			;endfor
			;printf,lunout,name[i]+' '+strtrim(median(apcor_meas,/even),2)+' '+strtrim(mad(apcor_meas)/sqrt(n_elements(apcor_meas)),2)
		endif



		vals=findgen(100)/99.*(0.5-0)-0.5
		resid=999.
		curr_best_val=999

		
		ind=where(apcor_meas gt -0.3 and apcor_meas lt 0.,num_valid)
		foreach val,vals do begin
			;print,val,total((val-apcor[ind])^2.)
			
			tot_weight=total(1./apcor_meas_err[ind]^2.)
			resid_new=0
			for el=0,n_elements(apcor_meas[ind])-1 do begin
				resid_new+=(val-apcor_meas[ind[el]])^2./apcor_meas_err[ind[el]]^2.
			endfor
			resid_new = resid_new / tot_weight
			;resid_new=total((val-apcor_meas[ind])^2./apcor_meas_err[ind]^2.)
			;print,val,resid_new,resid
			;print,apcor_meas
			if resid_new lt resid then begin
				resid=resid_new
				curr_best_val=val
			endif
		endforeach

		print,curr_best_val
		printf,lunout,name[i]+' '+strtrim(curr_best_val,2)+' '+strtrim(num_valid,2)


		cgplot,ap_mags,apcor_meas,psym=cgsymcat(16),yrange=[-0.4,0],ystyle=1
		oploterror,ap_mags,apcor_meas,apcor_meas_err,psym=3
		cgoplot,[10,16],[curr_best_val,curr_best_val],linestyle=1
                pausing=''
                read,pausing


	

	
	
	endfor ; images loop
	
	free_lun,lunout

endfor ; ext_loop


end
