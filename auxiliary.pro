; auxiliary.pro consists of various auxiliary functions for align.pro and edgeoffkeep_hanning.pro
; To compile:
;	GBTIDL -> .compile auxiliary.pro

; ---------------

function make_filepath,id,i,directory=directory,filetype=filetype,format=format

	; make_filepath
	; Usage: make_filepath(id, i [, directory=directory ] [, filetype=filetype ] [, format=format ])
	;
	; Description: Returns an input filepath, depending on format. Designed to work for raw GBT data.
	;
	; Parameters:
	; id (required, input, string): The ID of the dataset, beginning with "AGBT" and ending before the underscore,
	; 	such as "AGBT11B_051" or "AGBT16A_433"
	; i (required, input, int): The file number of the dataset
	; directory (optional, input, string): The relative directory path where the data can be found. This
	; 	requires the last backslash, such as "./reduced/"
	; filetype (optional, input, string): The filetype of the desired data. This requires the first period,
	; 	such as ".raw.acs.fits"
	; format (optional, input, int): An optional int (1, 2, or 3) to specify the format of the string. 
	;	1 for ".raw.acs.fits" files, where file numbers under 10 are appended with a 0
	;	2 for ".raw.vegas.A.fits" files in their own ".raw.vegas" directory, where file numbers under 10 are appended with a 0
	; 	3 for ".raw.vegas.A.fits" files in their own ".raw.vegas" directory, where file numbers under 4 are appended with 00
	;	This code is built around a few specific datasets, so these formats may be too specific or hardcoded for broad use
	; fpath (output, string): The resulting filepath 
	;
	; Examples:
	; 	GBTIDL -> make_filepath("AGBT11B_051", 1)
	;	"AGBT11B_051_01.raw.acs.fits"
	;
	; 	GBTIDL -> make_filepath("AGBT16A_433", 13, "./reducedData/", ".fits")
	;	"./reducedData/AGBT16A_433_13.raw.vegas/AGBT16A_433_13.raw.vegas.A.fits"

	if id eq "AGBT11B_051" or id eq "AGBT12A_266" or id eq "AGBT13A_312" then format = 1 else $
	if id eq "AGBT16A_433" then format = 2 else format = 3

	; default to format 1
	if not keyword_set(format) then format = 1
	if format ne 2 or format ne 3 then format = 1

	if not keyword_set(directory) then directory = ""

	if format eq 1 then begin
		if not keyword_set(filetype) then filetype = ".raw.acs.fits"
		if i gt 9 then begin
			fpath = strcompress(directory+id+"_"+string(i)+filetype,/remove_all)
		endif else begin
			fpath = strcompress(directory+id+"_0"+string(i)+filetype,/remove_all)
		endelse
	endif else if format eq 2 then begin
		if not keyword_set(filetype) then filetype = ".raw.vegas.A.fits"
		if i gt 9 then begin
			fpath = strcompress(directory+id+"_"+string(i)+".raw.vegas/"+id+"_"+string(i)+filetype,/remove_all)
		endif else begin
			fpath =strcompress(directory+id+"_0"+string(i)+".raw.vegas/"+id+"_0"+string(i)+filetype,/remove_all)
		endelse	
	endif else if format eq 3 then begin
		if not keyword_set(filetype) then filetype = ".raw.vegas.A.fits"
		if i lt 4 then begin
			fpath = strcompress(directory+id+"_00"+string(i)+".raw.vegas/"+id+"_00"+string(i)+filetype,/remove_all)
		endif else begin
			fpath = strcompress(directory+id+"_0"+string(i)+".raw.vegas/"+id+"_0"+string(i)+filetype,/remove_all)
		endelse
	endif

	return, fpath

end

; ---------------

function setdifference, a, b

	; setdifference
	; Usage: setdifference(a, b)
	;
	; Description: Returns the difference between two sets (arrays) a and b
	;	I probably found this code online. Thanks!
	;
	; Parameters:
	; a (required, input, array): Original array
	; b (required, input, array): Second array; the array to be "subtracted" from the original array
	; r (output, array): Difference of the two arrays (b subtracted from a)
	;
	; Examples:
	; 	GBTIDL -> setdifference([1,4,5,6], [5])
	;	[1,4,6]
	;
	;	GBTIDL -> setdifference([1,4,5,6], [7])
	;	[1,4,5,6]

	mina = Min(a,Max=maxa)
	minb = Min(b,Max=maxb)
	if (minb GT maxa) or (maxb lt mina) then return, a
	r = where((histogram(a, Min=mina, Max=maxa) NE 0) AND $
		(histogram(b, Min=mina, Max=maxa) EQ 0), count)
	if count eq 0 then return, -1 else return, r + mina
end

; ---------------

function get_my_nintegrations,scans

	; get_my_nintegrations
	; Usage: get_my_nintegrations(scans)
	;
	; Description: Given a series of scan numbers in a file, this returns an array of all the integration numbers
	;
	; Parameters:
	; scans (required, input, array): An array of scan numbers for which the user wants to obtain the number of
	; 	integrations
	; intarray (output, array): The array with the number of integrations per scan, for each of the scan numbers
	;
	; Examples:
	; 	GBTIDL -> get_my_nintegrations, [9,10,11,12]
	;	[225,225,224,225]
	;	; for a dataset where scans 9, 10, and 12 all have 225 integrations each, and scan 11 has 224 integrations

	; for some reason, GBTIDL requires me to initalize arrays with a value. Not to worry, as this can be removed at the end, but it seems inefficient
	intarray = [0]

	for i=0,n_elements(scans)-1 do begin
			info = scan_info(scans[i])
			nint = info.n_integrations
			intarray = [intarray, nint]
	endfor

	remove,0,intarray
	return,intarray

end

; ---------------

pro findlonggaps,scanid,value=value

	; findlonggaps
	; Usage: findlonggaps, scanid
	;
	; Description: Used to identify any particularly large gaps between integrations or missing integrations.
	;	This will print a list of all the the differences in longitude from integration to integration in a scan.
	;	It will also print the largest, smallest, and mean gap in longitude.
	; 	This needs work; so far, it's hardcoded just to include longitude_axis and I don't know how to change it
	; 	!g.s[0] doesn't take kindly to other string parameters unfortunately
	;	Kind of a useless function unfortunately but it can be interesting supposedly
	;
	; Parameters:
	; scanid (required, input, int); The scan number for which the user wants to see gaps in integrations
	;
	; Examples:
	; 	GBTIDL -> findlonggaps(9)
	;	[0.03	0.035	0.025	...	0.03]
	;	0.035
	; 	0.025
	; 	0.03
	;	; This tells us that the maximum longitude gap between integrations is 0.035, the minimum is 0.025,
	; 	; and the mean is 0.03. Supposedly, this can help locate discrepancies

	freeze

	;find number of integrations
	get, scan=scanid, int=0
	info = scan_info(scanid)
	nint = info.n_integrations
	; Obtain arrays of all longitudes. This can be changed to latitude_axis, or any other parameter of interest
	array1 = [!g.s[0].longitude_axis]
	array2 = [!g.s[0].longitude_axis]

	for i=1,nint-1 do begin
		gettp, scanid, int=i, /quiet
		array1 = [array1, !g.s[0].longitude_axis]
		array2 = [array2, !g.s[0].longitude_axis]
	endfor
	remove,nint-1,array1 ;all values but last
	remove,0,array2 ;all values but first
	differences = array2-array1 ;2-1 of value, and so on

	print, differences
	print, max(differences)
	print, min(differences)
	print, mean(differences)

end

; ---------------

pro getmissingints,id,start_ind,end_ind,ints,proc=proc,sname=sname,missingfiles=missingfiles,directory=directory, filetype=filetype, format=format

	; getmissingints
	; Usage: getmissingints, id, start_ind, end_ind, ints [, proc=proc ] [, sname=sname ] [, missingfiles=missingfiles ] [, directory=directory ] [, fileteype=filetype ] [, format=format ]
	;
	; Description: For a certain set of scans or files, return all the scan numbers with fewer integrations than expected and how many they have
	; 
	; Parameters:
	; id (required, input, string): The ID of the dataset, beginning with "AGBT" and ending before the underscore,
	; 	such as "AGBT11B_051" or "AGBT16A_433"
	; start_ind (required, input, int): The file number with which to start (such as 1 for "AGBT11B_51_01")
	; end_ind (required, input, int): The file number with which to end (such as 104 for "AGBT11B_051_104"). The function will analyze files start_end through
	;	end_ind in order.
	; ints (required, input, int/array): The number of integrations each scan should have. If more than one is acceptable, this can be passed as an array
	; proc (optional, input, string): Only analyze scans with this procedure (such as "Track", "RALongMap", etc.). Compatible with bash scripting.
	; sname (optional, input, string): Only analyze scans pointed at this source (such as "M31", "3C48", etc.). Compatible with bash scripting.
	; missingfiles (optional, input, array): Skip these file numbers when iterating from start_ind and end_ind
	; directory (optional, input, string): The relative directory path where the data can be found. This
	; 	requires the last backslash, such as "./reduced/"
	; filetype (optional, input, string): The filetype of the desired data. This requires the first period,
	; 	such as ".raw.acs.fits"
	; format (optional, input, int): An optional int (1, 2, or 3) to specify the format of the string. 
	;	1 for ".raw.acs.fits" files, where file numbers under 10 are appended with a 0
	;	2 for ".raw.vegas.A.fits" files in their own ".raw.vegas" directory, where file numbers under 10 are appended with a 0
	; 	3 for ".raw.vegas.A.fits" files in their own ".raw.vegas" directory, where file numbers under 4 are appended with 00
	;	This code is built around a few specific datasets, so these formats may be too specific or hardcoded for broad use
	;
	; Examples:
	; 	GBTIDL -> getmissingints, "AGBT11B_051", 1, 3, 225
	; 	File        1 scan           35          224 integrations
	;          	26          17          27
	; 	scans
	;	; Lists all the scans from AGBT11B_051_01 through AGBT11B_051_03 without 225 integrations
	;
	;	GBTIDL -> getmissingints, "AGBT12A_266", 1, 3, 300, proc="Track", sname="M31Ext*"
	;	File        1 scan         4235          294 integrations
	;	File        2 scan            8           60 integrations
	;	File        2 scan           10            3 integrations
	;	File        2 scan           11           10 integrations
	;	File        2 scan           18           29 integrations
	;			8          20          24
	;	scans
	;	; Lists all the scans from AGBT12A_266_01 through AGBT12A_266_03, using the procedure "Track" and any source beginning with "M31Ext"
	;	; without 300 scans

	; set defaults: Map scans
	if not keyword_set(missingfiles) then missingfiles=[0]
	if not keyword_set(proc) then proc="*Map"

	scans = [0]
	for i=start_ind, end_ind do begin
		ind = where(i eq missingfiles,count)
		if count eq 0 then begin
			fpath = make_filepath(id, i, directory=directory, filetype=filetype, format=format)
			filein, fpath
			;summary ;can comment this out, I'm just curious

			if not keyword_set(sname) then sn = get_scan_numbers(procedure=proc) else sn = get_scan_numbers(source=sname, procedure=proc)
			nints = get_my_nintegrations(sn)

			for j = 0,n_elements(nints)-1 do begin
				ind = where (nints[j] eq ints, count)
				if count eq 0 then print,'File ',i,' scan ',sn[j],' ',nints[j],' integrations'
			endfor
			scans = [scans, n_elements(sn)]
		endif
	endfor

	remove,0,scans
	print,scans,' scans'
end

; ---------------

pro summaries,id,start_ind,end_ind,missingfiles=missingfiles,directory=directory,filetype=filetype,format=format

	; summaries
	; Usage: summaries, id, start_ind, end_ind [, missingfiles=missingfiles ] [, directory=directory ] [, filetype=filetype ] [, format=format ]
	;
	; Description: Prints the summary for every file with a given id
	;
	; Parameters:
	; id (required, input, string): The ID of the dataset, beginning with "AGBT" and ending before the underscore,
	; 	such as "AGBT11B_051" or "AGBT16A_433"
	; start_ind (required, input, int): The file number with which to start (such as 1 for "AGBT11B_51_01")
	; end_ind (required, input, int): The file number with which to end (such as 104 for "AGBT11B_051_104"). The function will analyze files start_end through
	;	end_ind in order.
	; missingfiles (optional, input, array): Skip these file numbers when iterating from start_ind and end_ind
	; directory (optional, input, string): The relative directory path where the data can be found. This
	; 	requires the last backslash, such as "./reduced/"
	; filetype (optional, input, string): The filetype of the desired data. This requires the first period,
	; 	such as ".raw.acs.fits"
	; format (optional, input, int): An optional int (1, 2, or 3) to specify the format of the string. 
	;	1 for ".raw.acs.fits" files, where file numbers under 10 are appended with a 0
	;	2 for ".raw.vegas.A.fits" files in their own ".raw.vegas" directory, where file numbers under 10 are appended with a 0
	; 	3 for ".raw.vegas.A.fits" files in their own ".raw.vegas" directory, where file numbers under 4 are appended with 00
	;	This code is built around a few specific datasets, so these formats may be too specific or hardcoded for broad use
	;
	; Examples:
	;	GBTIDL -> summaries,"AGBT11B_051",1,3
	;	Scan           Source      Vel    Proc Seq    RestF nIF nInt nFd     Az    El
	;	-------------------------------------------------------------------------------
	;		5             3C48   -260.0   OffOn   1    1.420   4   20   1  249.7  77.5
	;	...
	; 	    39       M31_M33_01   -260.0 RALongM  31    1.420   1  225   1  306.6  15.3

	if not keyword_set(missingfiles) then missingfiles=[0]

	scans = [0]
	for i=start_ind, end_ind do begin
			ind = where(i eq missingfiles,count) 
		if count eq 0 then begin
			fpath = make_filepath(id, i, directory=directory, filetype=filetype, format=format)
			filein, fpath
			summary
		endif
	endfor

end

; ---------------

pro getrecords,id,file,directory=directory,filetype=filetype,format=format

	; getrecords
	; Usage: getrecords, id, file [, directory=directory ] [, filetype=filetype ] [, format=format ]
	;
	; Description: Displays gettp for every scan and polarization for a given id and file number, for the purpose of quick visual inspection
	;	Note that this only displays gettp and only two polarizations per scan. This function can be generalized to include more than gettp
	;	and more records per scan, such as records for differing IFNUM or FDNUM
	;
	; Parameters:
	; id (required, input, string): The ID of the dataset, beginning with "AGBT" and ending before the underscore,
	; 	such as "AGBT11B_051" or "AGBT16A_433"
	; file (required, input, int): The file number of the file requested
	; directory (optional, input, string): The relative directory path where the data can be found. This
	; 	requires the last backslash, such as "./reduced/"
	; filetype (optional, input, string): The filetype of the desired data. This requires the first period,
	; 	such as ".raw.acs.fits"
	; format (optional, input, int): An optional int (1, 2, or 3) to specify the format of the string. 
	;	1 for ".raw.acs.fits" files, where file numbers under 10 are appended with a 0
	;	2 for ".raw.vegas.A.fits" files in their own ".raw.vegas" directory, where file numbers under 10 are appended with a 0
	; 	3 for ".raw.vegas.A.fits" files in their own ".raw.vegas" directory, where file numbers under 4 are appended with 00
	;	This code is built around a few specific datasets, so these formats may be too specific or hardcoded for broad use
	;
	; Examples:
	;	GBTIDL -> getrecords,"AGTBT11B_51",1
	;	Scan:     5 (IF:0 FD:0 PL:0)    Tsys:  16.70
	;	Scan:     5 (IF:0 FD:0 PL:1)    Tsys:  16.57
	; 	...
	;	Scan:    39 (IF:0 FD:0 PL:1)    Tsys:  20.59
	;	; Note that gettp spectra will be displayed in the GBTIDL Plotter window as the function executes, so that the user can visually inspect scans

	; set defaults: directory "reduced/"
	if not keyword_set(directory) then directory = "reduced/"

	fpath = make_filepath(id,file,directory=directory,filetype=filetype,format=format)
	filein, fpath
	sn = get_scan_numbers()
	for i=0,n_elements(sn)-1 do begin
		for j=0,1 do begin
			gettp,sn[i],plnum=j,int=k
		endfor
	endfor

end

; ---------------

pro finddiscrepantints,scan

	; finddiscrepantints
	; Usage: finddiscrepantints, scan
	;
	; Description: Displays three neighboring integrations plotted atop each other, for each set of integrations on a scan, for the purpose of
	;	quickly visually inspecting any highly discrepant integration. This function requires that a file has already been loaded, and the specified
	; 	scan must be in the loaded file
	; 
	; Parameters:
	; scan (required, input, int): The scan for which the user wants to inspect the integrations
	;
	; Examples
	;	GBTIDL -> finddiscrepantints,39
	;	Scan:    39 (IF:0 FD:0 PL:0)    Tsys:  21.35
	;	int       0
	;	Scan:    39 (IF:0 FD:0 PL:0)    Tsys:  21.35
	;	Scan:    39 (IF:0 FD:0 PL:0)    Tsys:  21.17
	;	Scan:    39 (IF:0 FD:0 PL:0)    Tsys:  21.33
	;	...
	;	Scan:    39 (IF:0 FD:0 PL:0)    Tsys:  20.21
	;	; Note that the integration spectra, overlaid atop each other, will be displayed in the GBTIDL Plotter window as the function executes:
	;	; for example, integrations 0, 1, and 2 together. This is so the user can quickly see if any integrations are far different from their neighbors

	freeze
	gettp,scan,int=0
	info = scan_info(scan)
	nint = info.n_integrations
	for i=0,nint-3 do begin
		freeze
		print, "int", i
		gettp,scan,int=i 
		copy,0,1
		gettp,scan,int=i+1
		copy,0,2
		gettp,scan,int=i+2
		copy,0,3
		unfreeze
		show,1
		oshow,2,color=!green
		oshow,3,color=!blue
		; This shows integrations of all three spectra. Pause for a second so the user can see it.
		wait,1
	endfor

end

; ---------------

pro getcoords,id,start_ind,end_ind,u,writepath,proc=proc,directory=directory,filetype=filetype,format=format

	; getcoords
	; Usage: getcoords, id, start_ind, end_ind, u, writepath [, proc=proc ] [, directory=directory ]
	;
	; Description: Records ID, scan number, polarization number, sequence number, azimuth, elevation, longitude, and latitude for a given ID and range of files
	;	and saves them to a files. This takes awhile to run and might not be entirely worthwhile, but can obtain lots of coordinates rather quickly.
	;
	; Parameters:
	; id (required, input, string): The ID of the dataset, beginning with "AGBT" and ending before the underscore,
	; 	such as "AGBT11B_051" or "AGBT16A_433"
	; start_ind (required, input, int): The file number with which to start (such as 1 for "AGBT11B_51_01")
	; end_ind (required, input, int): The file number with which to end (such as 104 for "AGBT11B_051_104"). The function will analyze files start_end through
	;	end_ind in order.
	; u (required, input, int): File unit number required for IDL to write to a file. This can be obtained by running "GBTIDL -> get_lun, u" and is usually
	;	a small integer
	; writepath (required, input, string): Relative path to which the file should be written
	; proc (optional, input, string): Only return scans with this procudure (such as "Track", "RALongMap", etc.). Compatible with bash scripting.
	; directory (optional, input, string): The relative directory path where the data can be found. This
	; 	requires the last backslash, such as "./reduced/"
	;
	; Examples:
	;	GBTIDL -> getcoords,"AGBT11B_051",3,3,100,"test.txt"
	;	Scan:    13 (IF:0 FD:0 PL:0)    Tsys:  16.56
	;		272.28115       69.216145       15.001065       36.250417
	;	Scan:    13 (IF:0 FD:0 PL:1)    Tsys:  16.50
	;		272.28115       69.216145       15.001065       36.250417
	;	...
	;	Scan:    39 (IF:0 FD:0 PL:1)    Tsys:  20.59
	;		306.58279       15.321734       14.919557       37.812662
	;	; Output of test.txt:
	;	AGBT11B_051_03          13           0           5       272.28115
	;		69.216145       15.001065       36.250417
	;	AGBT11B_051_03          39           1          31       306.58279
	;		15.321734       14.919557       37.812662
	;	; az/el and long/lat of each scan can be found here
  
	; default to map-type scans
	if not keyword_set(proc) then proc="*Map"

	for k=start_ind,end_ind do begin
		fpath = make_filepath(id, k, directory=directory, filetype=filetype, format=format)
		filein, fpath

		sn = get_scan_numbers(procedure=proc)
		sz = size(sn, /n_elements)
		; int = get_my_nintegrations(sn) ; not using these I guess - no need

		for j=0,sz-1 do begin
			scan = sn[j]
			for i=0,1 do begin
				if proc eq "*Map" then begin
					gettp,scan,plnum=i
				endif else begin
					getfs,scan,plnum=i,tcal=1.6738227704817121 ; surely I can get away with one (1) hardcoding
				endelse
				
				data = [!g.s[0].azimuth, !g.s[0].elevation, !g.s[0].longitude_axis, !g.s[0].latitude_axis]
				print,data
				openw,u,writepath,/append
				printf,u,!g.s[0].projid, !g.s[0].scan_number, !g.s[0].polarization_num, !g.s[0].procseqn, data
				close,u
				
				; wait,1 ; why was this here?
			endfor
		endfor
	endfor

end

; ---------------

pro processduplicates,id,start_ind,end_ind,directory=directory,filetype=filetype,format=format

	; processduplicates
	; Usage: processduplicates, id, start_ind, end_ind [, directory=directory ] [, filetype=filetype ] [, format=format ]
	;
	; Description: This function goes through all files with a specified ID and range of file numbers and flags any duplicate scans. The default
	;	directory is "../formap/", default filetype is ".keep.fits", and default format is 1, but this can certainly be changed
	;	Some functions inadvertently duplicate old data if a certain integration is missing, and this function exists to correct for that.
	; 	This particular function, however, was only used for the output of LSRDopplerCorrection.pro
	;	This function ONLY looks at scan number, integration number, and polarization number. In the case of multiple IFNUMs or FDNUMs, it may
	;	wrongly flag differing IFNUMs or FDNUMs as duplicates
	;	Essentially this is an extremely specific function and the specific issue has been fixed so I kind of doubt this will be used ever
	;
	; Parameters:
	; id (required, input, string): The ID of the dataset, beginning with "AGBT" and ending before the underscore,
	; 	such as "AGBT11B_051" or "AGBT16A_433"
	; start_ind (required, input, int): The file number with which to start (such as 1 for "AGBT11B_51_01")
	; end_ind (required, input, int): The file number with which to end (such as 104 for "AGBT11B_051_104"). The function will analyze files start_end through
	;	end_ind in order.
	; directory (optional, input, string): The relative directory path where the data can be found. This
	; 	requires the last backslash, such as "./reduced/"
	; filetype (optional, input, string): The filetype of the desired data. This requires the first period,
	; 	such as ".raw.acs.fits"
	; format (optional, input, int): An optional int (1, 2, or 3) to specify the format of the string. 
	;	1 for ".raw.acs.fits" files, where file numbers under 10 are appended with a 0
	;	2 for ".raw.vegas.A.fits" files in their own ".raw.vegas" directory, where file numbers under 10 are appended with a 0
	; 	3 for ".raw.vegas.A.fits" files in their own ".raw.vegas" directory, where file numbers under 4 are appended with 00
	;	This code is built around a few specific datasets, so these formats may be too specific or hardcoded for broad use
	;
	; Examples:
	; 	GBTIDL -> processduplicates, "AGBT11B_051", 105, 106
	; 	duplicate	105	10515	224	0
	;	; polarization 0 of integration 224, scan 10515, of file AGBT11B_051_105 was duplicated and flagged
	;
	;	GBTIDL -> processduplicates, "AGBT16A_433", 6, 8, directory="./M31-M33-data/", filetype=".raw.acs.fits", format=3
	;	duplicate	7	707		223	1
	;	duplicate	7	707		223	1
	;	; polarization 1 of integration 223, scan 707, of file AGBT16A_433_07 was duplicated twice, and both duplicates were flagged
	;	; note that this function will probably never need to be run on the raw data (with filetype ".raw.acs.fits")

	; default to directory "formap/"
	if not keyword_set(directory) then directory = "formap/"

	for file=start_ind,end_ind do begin
		infile = make_filepath(id,file,directory=directory,filetype=filetype,format=format)
		filein, infile

		nr = nrecords()

		getrec,0
		for i=0,nr-2 do begin
			copy,0,1
			getrec,i+1
			; if the parameters are exactly the same as the previous, flag it as a duplicate
			if (!g.s[0].polarization_num eq !g.s[1].polarization_num) and (!g.s[0].integration eq !g.s[1].integration) and (!g.s[0].scan_number eq !g.s[1].scan_number) then begin
				flagrec,i+1,idstring="duplicate"
				print, "duplicate", file, !g.s[0].scan_number, !g.s[0].integration, !g.s[0].polarization_num
			endif
		endfor
	endfor
end

; ---------------

pro accumdeeppointings,id,start_ind,end_ind,pointing=pointing,stype=stype,tcal=tcal,indir=indir,infiletype=infiletype,outdir=outdir,outfiletype=outfiletype,format=format

	; accumdeeppointings
	; Usage: accumdeeppointings, id, start_end, end_ind [, pointing=pointing ] [, stype=stype ] [, tcal=tcal ] [, indir=indir ] [, infiletype=infiletype ] [, outdir=outdir ] [, outfiletype=outfiletype ] [, format=format ]
	;
	; Description: This function accumulates and averages deep pointing scans. It defaults to source name "M31Ext1", tcal = 1.674, input directory "M31-M33-data/",
	;	input filetype ".raw.acs.fits", output directory "deepPointings/", and output filetype ".keep.fits", but all these can be changed. This is supposed to
	;	gather and average all integrations, scans, and polarizations of a deepPointing source into one keep file.
	;
	; Parameters:
	; id (required, input, string): The ID of the dataset, beginning with "AGBT" and ending before the underscore,
	; 	such as "AGBT11B_051" or "AGBT16A_433"
	; start_ind (required, input, int): The file number with which to start (such as 1 for "AGBT11B_51_01")
	; end_ind (required, input, int): The file number with which to end (such as 104 for "AGBT11B_051_104"). The function will analyze files start_end through
	;	end_ind in order.
	; pointing (optional, input, int): The pointing number to consolidate
	; stype (optional, input, string): The source type for which to consolidate pointings. This should not have a number at the end, for example, "M31Ext"
	; indir (optional, input, string): The relative directory path where the data can be found. This
	; 	requires the last backslash, such as "./reduced/"
	; infiletype (optional, input, string): The filetype of the desired input data. This requires the first period,
	; 	such as ".raw.acs.fits"
	; outdir (optional, input, string): The relative directory path where the data should be output. This
	; 	requires the last backslash, such as "./reduced/"
	; outfiletype (optional, input, string): The filetype of the desired output data. This requires the first period,
	; 	such as ".raw.acs.fits"
	; format (optional, input, int): An optional int (1, 2, or 3) to specify the format of the string. 
	;	1 for ".raw.acs.fits" files, where file numbers under 10 are appended with a 0
	;	2 for ".raw.vegas.A.fits" files in their own ".raw.vegas" directory, where file numbers under 10 are appended with a 0
	; 	3 for ".raw.vegas.A.fits" files in their own ".raw.vegas" directory, where file numbers under 4 are appended with 00
	;	This code is built around a few specific datasets, so these formats may be too specific or hardcoded for broad use
	; 
	; Examples:
	; 	GBTIDL -> accumdeeppointings, "AGBT12A_266", 1, 13, 1, "M31Ext", tcal=1.77
	;	; This compiles all deepPointings M31Ext1, searching from AGBT12A_266_01 to AGBT12A_266_13 and tcal 1.77, into file "deepPointings/AGBT12A_266_01.keep.fits"

	; this seems a bit hard-code-y, but set defaults: M31Ext, input M31-M33-data/[file].raw.acs.fits, and output deepPointings/[file].keep.fits
	if not keyword_set(pointing) then pointing = 1
	if not keyword_set(stype) then stype = "M31Ext"
	if not keyword_set(tcal) then tcal = 1.6738227704817121
	if not keyword_set(indir) then indir = "M31-M33-data/"
	if not keyword_set(infiletype) then infiletype = ".raw.acs.fits"
	if not keyword_set(outdir) then outdir = "deepPointings/"
	if not keyword_set(outfiletype) then outfiletype = ".keep.fits"

	for i=start_ind,end_ind do begin
		infile = make_filepath(id,i,directory=indir,filetype=infiletype,format=format)
		outfile = make_filepath(id,i,directory=outdir,filetype=outfiletype,format=format)
		filein, infile
		fileout, outfile

		pointing = strcompress("M31Ext"+string(pointing), /remove_all)
		sn = get_scan_numbers(source=pointing)
		sz = size(sn, /n_elements)

		; wait all the integrations too? I don't think so
		for j = 0,sz-1 do begin
			for p = 0,1 do begin
				getfs,sn[j],tcal=tcal,plnum=p
				accum
			endfor
		endfor
	endfor
	ave
	keep
end