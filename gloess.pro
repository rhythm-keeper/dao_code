function gloess,data,smooth_factor




;==========================
; First bin the data
;==========================
; Set boundaries for LF
max_data_val=max(data)
min_data_val=min(data)
binsize=0.01
n_data_bins=ceil( (max_data_val-min_data_val)/binsize )


binned_array_count=fltarr(n_data_bins)
binned_array_value=fltarr(n_data_bins)
for bin_el=0L,n_data_bins-1 do begin

        ; Get stars within the current bin range
        within_range_index=where( data ge min_data_val+binsize*bin_el and $
                                data lt min_data_val+binsize*(bin_el+1) , number_within_range)

        binned_array_count[bin_el]=number_within_range
        binned_array_value[bin_el]=min_data_val+binsize*(bin_el+0.5)

endfor

;=============================================
; Now loop through an evenly spaced
; grid of points based on the smooth_factor, 
; compute the distance between points, compute
; the weighted contribution based on a Gaussian
; then record the weighted mean of all points
;=============================================

weight=fltarr(n_data_bins) ; for containing a Gaussian weights from center of bin
new_LF=fltarr(n_data_bins) ; to contained smoothed bin

for bin_el=0L,n_data_bins-1 do begin

        ; The position used as the current reference point
        curr_position=min_data_val+binsize*(bin_el+0.5)

        ; Loop through other data points and compute their contribution
        for other_bin_el=0L,n_data_bins-1 do begin

                distance=abs(curr_position - (min_data_val+binsize*(other_bin_el+0.5)) )
                weight[other_bin_el]=exp( -1.0*(distance/smooth_factor)^2. )

        endfor ; done looping through all bins

        ; Compute weighted mean
        dum=0
        total_weight=total(weight)
        for bin_loop=0L,n_data_bins-1 do begin
                dum+=binned_array_count[bin_loop]*weight[bin_loop]/total_weight
        endfor
        new_LF[bin_el]=dum


endfor ; done evaluating function at other bins


; Normalize dataset so that total number of edge response equals
; the input number of stars
norm_factor=total(binned_array_count)/total(new_LF)
new_LF*=norm_factor


; Sobel [-1,0,+1]

; Sobel [-1,0,+1]
new_sobel=[]
sigma = []
xvals=[]
for loop=1,n_data_bins-2 do begin
    xvals     = [xvals,     min_data_val+binsize*(loop+0.5)]
    new_sobel = [new_sobel, (new_LF[loop+1] - new_LF[loop-1])]
    sigma    =  [sigma,     sqrt(new_LF[loop+1]+new_LF[loop-1])]
endfor

weight = abs(new_sobel)/sigma
total_weight = total(weight)

new_sobel = new_sobel * weight/total_weight


ind=where(new_sobel eq max(new_sobel))
peak_loc=xvals[ind]
print,'peak location: '+strtrim(peak_loc,2)

return,peak_loc

end

