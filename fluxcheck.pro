pro fluxcheck,id,start_ind,end_ind,u,writepath,directory=directory,filetype=filetype,format=format,proc=proc,leftbound=leftbound,rightbound=rightbound,order=order

	; fluxcheck
	; Usage: fluxcheck, id, start_ind, end_ind, u, writepath [, directory=directory ] [, filetype=filetype ] [, format=format ] [, proc=proc ]
	;
	; Description: This function lets the user manually select a region, fits and subtracts a baseline, and extracts the peak intensity from the spectrum.
	; 	It prints the ID, scan number, polarization number, peak flux, standard deviation, az/el, and RA/Dec into a text file.
	;	This defaults to a 2nd-order baseline and zooms into the range 1.418 - 1.423 MHz (neutral hydrogen observations), but these can be changed
	;	It also defaults to ".raw.acs.fits" files
	;	Since this outputs both elevation and peak intensity, this can be used to find calibration parameters from calibration/OffOn scans. Some work
	; 	to be done to convert the output into a machine-readable CSV file. findTcal.py is designed to conver this txt file into a CSV.
	;	auxiliary.pro must be compiled first
	;
	;   To compile:
	;       GBTIDL -> .com auxiliary.pro
	;       GBTIDL -> .com fluxcheck.pro
	;
	;	Code by Bethany Grimm
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
	; directory (optional, input, string): The relative directory path where the data can be found. This
	; 	requires the last backslash, such as "./reduced/"
	; filetype (optional, input, string): The filetype of the desired data. This requires the first period,
	; 	such as ".raw.acs.fits"
	; format (optional, input, int): An optional int (1, 2, or 3) to specify the format of the string. 
	;	1 for ".raw.acs.fits" files, where file numbers under 10 are appended with a 0
	;	2 for ".raw.vegas.A.fits" files in their own ".raw.vegas" directory, where file numbers under 10 are appended with a 0
	; 	3 for ".raw.vegas.A.fits" files in their own ".raw.vegas" directory, where file numbers under 4 are appended with 00
	;	This code is built around a few specific datasets, so these formats may be too specific or hardcoded for broad use
	; proc (optional, input, string): Only analyze scans with this procedure (such as "Track", "RALongMap", etc.). Compatible with bash scripting.
	;	Defaults to "OffOn"
	; leftbound (optional, input, float): The lower bound of the baseline region, in MHz - the user will still manually set the region within
	;	this bound. Defaults to 1.418 MHz
	; rightbound (optional, input, float): The upper bound of the baseline region, in MHz - the user will still manually set the region within
	;	this bound. Defaults to 1.423 MHz
	; order (optional, input, int): The order of the polynomial used to fit the baseline. Defaults to 2nd order
	;
	; Examples:
	;	GBTIDL -> fluxcheck, "AGBT11B_051", 1, 104, 100, "output.txt"
	;	Projid ScanNum Plnum Sflux StdSflux MeanTcal Az El Ra Dec
	;	AGBT11B_051_01           6           0       27.814917     0.016367706
	;		1.4678938       249.71561       77.464752       24.422085       33.159630
	;	...
	;	AGBT11B_051_104       10407           1       26.778819     0.020003691
	;       1.4552635       67.344756       26.600959       24.421477       33.159755
	;	; Note that user input is required to manually set baselining regions.
  
  for k=start_ind,end_ind do begin
  	; In case make_filepath doesn't work -
	; fpath = strcompress(id+string(k)+".raw.acs.fits", /remove_all)

	; set default filetype to ".raw.acs.fits" and default procedure to "OffOn"
	if not keyword_set(filetype) then filetype=".raw.acs.fits"
	if not keyword_set(proc) then proc="OffOn"
	fpath = make_filepath(id,k,directory=directory,filetype=filetype,format=format)
	filein, fpath

	 sn = get_scan_numbers(procedure=proc)
	 sz = size(sn, /n_elements)

	 for j=1,sz-1,2 do begin
		 scan = sn[j]
		 print,"Scan number: ",scan
		 for i=0,1 do begin
			nfit,2
  			
	  		getps,scan,plnum=i
	  		freq
	  		setxunit,"GHz"
	  		
	  		setx,1.418,1.423
	  		
	  		setregion
	  		showregion
			
	  		bshape,modelbuffer=1
			
	  		a=getdata(1)
	  		x=getxarray()
	  		ok=where(x ge 1.42 and x lt 1.421)
			
	  		data = [mean(a[ok]), stdev(a[ok]), !g.s[0].mean_tcal, !g.s[0].azimuth, !g.s[0].elevation, !g.s[0].longitude_axis, !g.s[0].latitude_axis]
	  		print,data
	  		openw,u,writepath,/append
	  		printf,u,!g.s[0].projid, !g.s[0].scan_number, !g.s[0].polarization_num, data
	  		close,u
	  		
	  		wait,1
		end
	end
  	print,!g.s[0].projid
  end

end