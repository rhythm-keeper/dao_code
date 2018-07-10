pro reduce_master_srt_file,input


readcol,'id.list',id


for i=0,n_elements(id)-1 do begin

	query="grep -A1 "+string(39B)+" "+strtrim(string(id[i],format='(I)'),2)+" "+string(39B)+" "+input+".srt"
	spawn,query
	print," "


endfor


end
