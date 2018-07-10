#!/usr/bin/env

import numpy as np
import sys
import os

#import pyfits
import astropy.io.fits as pyfits


def filelines(filename):
	linenum=0
	with open(filename) as f:
		for row in f:
			linenum+=1
	return linenum


def weighted_avg(mags,magerrs):

	tot_weight=0
	dum=0
	for mag_el,magerr_el in zip(mags,magerrs):
		weight=1./magerr_el**2
		tot_weight+=1./magerr_el**2
		dum+=mag_el*weight

	avg_mag=dum/tot_weight

	inv_sqr_mag=[1./x**2 for x in magerrs]
	avg_magerr=np.sqrt( 1./np.sum(inv_sqr_mag) )

	return [avg_mag,avg_magerr]

	

def get_avg_mag_magerr(mag,magerr):
	"""
	Flux averaged magnitudes
	"""

	fluxes=[]
	fluxerrs=[]
	
	for mag_el,magerr_el in zip(mag,magerr):
		flux_el=10**( (mag_el-25.)/(-2.5) )
		fluxerr_el=flux_el*magerr_el/1.0857

		fluxes.append(flux_el)
		fluxerrs.append(fluxerr_el)

	# wavg flux
	tot_weight=0
	dum=0
	for flux_el,fluxerr_el in zip(fluxes,fluxerrs):
		weight=1./fluxerr_el**2
		tot_weight+=1./fluxerr_el**2
		dum+=flux_el*weight	

	avg_flux=dum/tot_weight

	inv_sqr_flux=[1./x**2 for x in fluxerrs]
	avg_fluxerr=np.sqrt( 1./np.sum(inv_sqr_flux) )

	#print fluxes,fluxerrs
	#print avg_flux,avg_fluxerr
	
	avg_mag=-2.5*np.log10(avg_flux)+25
	avg_magerr=avg_fluxerr*1.0857/avg_flux

	return [avg_mag,avg_magerr]



def get_avg_apcors():

	first_band_apcors=list()
	first_band_apcorerrs=list()
	second_band_apcors=list()
	second_band_apcorerrs=list()
	
	# first band
	for image_name in first_band_image_names:
	
		with open(apcor_filename) as f:
	
			# record data values for a given band+image combo 
			# before averging
			temp_apcor=list()
			temp_apcorerr=list()
		
			for row in f:
				parts=row.split()
			
		
				if parts[0] == image_name:
					temp_apcor.append(float(parts[1]))
					temp_apcorerr.append(float(parts[2]))   
		
		avg_mag,avg_magerr=weighted_avg(temp_apcor,temp_apcorerr)
		first_band_apcors.append(avg_mag)
		first_band_apcorerrs.append(avg_magerr)
	

	# second band
	for image_name in second_band_image_names:
	
		with open(apcor_filename) as f:
	
			# record data values for a given band+image combo 
			# before averging
			temp_apcor=list()
			temp_apcorerr=list()
		
			for row in f:
				parts=row.split()
			
		
				if parts[0] == image_name:
					temp_apcor.append(float(parts[1]))
					temp_apcorerr.append(float(parts[2]))   
		
		avg_mag,avg_magerr=weighted_avg(temp_apcor,temp_apcorerr)
		second_band_apcors.append(avg_mag)
		second_band_apcorerrs.append(avg_magerr)


		
	return [first_band_apcors,first_band_apcorerrs,second_band_apcors,second_band_apcorerrs]



def get_average_gain(sample_filename):


	# get average gain
	hdulist = pyfits.open(sample_filename+'.fits')
	GAINA=hdulist[0].header['ATODGNA']
	GAINB=hdulist[0].header['ATODGNB']
	GAINC=hdulist[0].header['ATODGNC']
	GAIND=hdulist[0].header['ATODGND']


	return np.average([GAINA,GAINB,GAINC,GAIND])


def calculate_zps():

	band1_zp=list()
	band2_zp=list()

	# get average again
	# only need one image since using the same detector
	avg_gain = get_average_gain(first_band_image_names[0])

	for image_name,curr_meas_apcor in zip(first_band_image_names,first_band_apcors):

		# get exposure time from FITS image	
		hdulist = pyfits.open(image_name+'.fits')
		EXPTIME=hdulist[0].header['EXPTIME']

		zp_el=(band1_phot_zp-25)+2.5*np.log10(EXPTIME/avg_gain)-band1_inf_apcor+curr_meas_apcor
		#print band1_inf_apcor,curr_meas_apcor,zp_el

		band1_zp.append(zp_el)


	for image_name,curr_meas_apcor in zip(second_band_image_names,second_band_apcors):
											       
		# get exposure time from FITS image	
		hdulist = pyfits.open(image_name+'.fits')
		EXPTIME=hdulist[0].header['EXPTIME']

		zp_el=(band2_phot_zp-25)+2.5*np.log10(EXPTIME/avg_gain)-band2_inf_apcor+curr_meas_apcor
											       
		band2_zp.append(zp_el)



	return [band1_zp,band2_zp]





# Important input
filename=sys.argv[1]
apcor_filename=sys.argv[2]
ext_no=sys.argv[3]

n_images=filelines("all_"+ext_no+".list")
second_band_index=6 # index in list where second filter appears 
apcor_04inf_f814w=0.09763
apcor_04inf_f606w=0.09526
f814w_zp=25.517
f606w_zp=26.405

band1_phot_zp=f814w_zp
band1_inf_apcor=apcor_04inf_f814w
band2_phot_zp=f606w_zp
band2_inf_apcor=apcor_04inf_f606w


# get names of images per band
first_band_image_names=list()
second_band_image_names=list()
row_index=0
with open("all_"+ext_no+".list") as f:
	for row in f:
		parts=row.split()

		if row_index < second_band_index: 
			first_band_image_names.append(parts[0])
		else:
			second_band_image_names.append(parts[0])
		row_index+=1



# average aperture correction per frame
first_band_apcors,first_band_apcorerrs,second_band_apcors,second_band_apcorerrs = get_avg_apcors()

band1_zp,band2_zp = calculate_zps()
print (band1_zp,band2_zp)
exit(0)


num_rows_per_object=int(n_images / 6) + 1


mag_arr=np.zeros(n_images)
magerr_arr=np.zeros(n_images)

# count number of stars in file
nlines=filelines(filename)
nstars=(nlines-3)/num_rows_per_object

with open(filename) as f:
	
	# Header
	head1=f.readline()
	head2=f.readline()
	head3=f.readline()

	outfile1=open("f814w_"+ext_no+".cal","w")
	outfile1.write(head1)
	outfile1.write(head2)
	outfile1.write(head3)
	outfile1.close()


	outfile2=open("f606w_"+ext_no+".cal","w")
	outfile2.write(head1)
	outfile2.write(head2)
	outfile2.write(head3)
	outfile2.close()

	
	star_num=1


	while star_num <= nstars:

		# put all data for a star into a long list
		all_parts=list()
		for i in range(num_rows_per_object):
			all_parts.append(f.readline().split())

		# flatten list
		parts=[item for sublist in all_parts for item in sublist]

		# loop through star values and record them
		part_num=0

		id_el=int(parts[part_num])
		part_num+=1
		x_el=float(parts[part_num])
		part_num+=1
		y_el=float(parts[part_num])
		part_num+=1
		
	
		p_mag_magerr_pairs=int( (len(parts)-3-2)/2. )  # exclude id,x,y and chi,sharp
		for col in range(0, p_mag_magerr_pairs): 
			mag_arr[col]=float(parts[part_num])
			part_num+=1			
			magerr_arr[col]=float(parts[part_num])
			part_num+=1

		# now do chi, sharp
		chi_el=float(parts[part_num])
		part_num+=1	
		sharp_el=float(parts[part_num])
		part_num+=1


		# add zero-points to magnitudes
		row_index=0
		while row_index < n_images:	
			if row_index < second_band_index:
				mag_arr[row_index] += band1_zp[row_index]
			else:
				mag_arr[row_index] += band2_zp[row_index-second_band_index]
			row_index+=1


		# convert to fluxes, get average mag
		avg_mag_band1,avg_magerr_band1 = get_avg_mag_magerr(mag_arr[0:second_band_index],magerr_arr[0:second_band_index])

		avg_mag_band2,avg_magerr_band2 = get_avg_mag_magerr(mag_arr[second_band_index:n_images],magerr_arr[second_band_index:n_images])


		# write to file
		outfile1=open("f814w_"+ext_no+".cal","a")
		id_string=str(id_el)
		x_string =("{0:.3f}").format(x_el)
		y_string =("{0:.3f}").format(y_el)
		mag1_string=("{0:.3f}").format(avg_mag_band1)
		magerr1_string=("{0:.3f}").format(avg_magerr_band1)
		fake_sky="100"
		fake_nit="10"
		chi_string=("{0:.2f}").format(chi_el)
		sharp_string=("{0:.3f}").format(sharp_el)

		outfile1.write("   "+id_string+" "+\
			x_string+" "+y_string+" "+\
			mag1_string+" "+magerr1_string+"   "+\
			fake_sky+"       "+fake_nit+"    "+\
			chi_string+"    "+sharp_string+"\n")

		outfile1.close()


		outfile2=open("f606w_"+ext_no+".cal","a")
		id_string=str(id_el)
		x_string =("{0:.3f}").format(x_el)
		y_string =("{0:.3f}").format(y_el)
		mag2_string=("{0:.3f}").format(avg_mag_band2)
		magerr2_string=("{0:.3f}").format(avg_magerr_band2)
								    
		outfile2.write("   "+id_string+" "+\
			x_string+" "+y_string+" "+\
			mag2_string+" "+magerr2_string+"   "+\
			fake_sky+"       "+fake_nit+"    "+\
			chi_string+"    "+sharp_string+"\n")


								    
		outfile2.close()



		
		star_num+=1

		# clear list array
		del all_parts
		# clear other single instance vars
		#del id_el, x_el, y_el, chi_el, sharp_el


