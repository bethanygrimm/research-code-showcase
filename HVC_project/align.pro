pro align, id, sstart, send, reffile, indir=indir, infiletype=infiletype, outdir=outdir, outfiletype=outfiletype, format=format

  ; align
  ; Usage: align, id, sstart, send, reffile [, indir=indir ] [, infiletype=infiletype ] [, outdir=outdir ] [, outfiletype=outfiletype ] [, format=format ]
  ;
  ; Description: Aligns all files from sstart to send to reference file reffile. This is to correct for small errors in velocity and align spectra
  ;   belonging to the same map together. Defaults to align ".keep.fits" files from the directory "newreduced/",
  ;   and output as ".keep.fits" files into the directory "formap/".
  ;   auxiliary.pro must be compiled first
  ;
  ;   for AGBT11B_051 session 01 and 02 are different from the sessions
  ;   03-104.  This means that the velocities for 01 and 02 can not be
  ;   simply aligned to match the remaining data
  ;
  ;   add a check to make sure data array has changed
  ;   if it hasn't changed then do not save the data
  ;
  ;   To compile:
  ;       GBTIDL -> .com auxiliary.pro
  ;       GBTIDL -> .com align.pro
  ;
  ;   Code by Toney Minter 
  ;   Edited by Bethany Grimm
  ;
  ; Parameters:
  ; id (required, input, string): The ID of the dataset, beginning with "AGBT" and ending before the underscore,
  ; 	such as "AGBT11B_051" or "AGBT16A_433"
  ; sstart (required, input, int): The file number with which to start (such as 1 for "AGBT11B_51_01")
  ; send (required, input, int): The file number with which to end (such as 104 for "AGBT11B_051_104"). The function will analyze files sstart through
  ;	send in order.
  ; reffile (required, input, int): The file number with which all files should be aligned (the reference file)
  ; indir (optional, input, string): The relative directory path where the data can be found. This
  ; 	requires the last backslash, such as "./reduced/"
  ; infiletype (optional, input, string): The filetype of the desired input data. This requires the first period,
  ; 	such as ".raw.acs.fits"
  ; outdir (optional, input, string): The relative directory path where the data should be output. This
  ; 	requires the last backslash, such as "./reduced/"
  ; outfiletype (optional, input, string): The filetype of the desired output data. This requires the first period,
  ; 	such as ".raw.acs.fits"
  ; format (optional, input, int): An optional int (1, 2, or 3) to specify the format of the string. 
  ;	  1 for ".raw.acs.fits" files, where file numbers under 10 are appended with a 0
  ;	  2 for ".raw.vegas.A.fits" files in their own ".raw.vegas" directory, where file numbers under 10 are appended with a 0
  ; 	3 for ".raw.vegas.A.fits" files in their own ".raw.vegas" directory, where file numbers under 4 are appended with 00
  ;	  This code is built around a few specific datasets, so these formats may be too specific or hardcoded for broad use
  ;
  ; Examples:
  ;   GBTIDL -> align, "AGBT11B_051", 3, 104, 5
  ;   ; Align all scans in files AGBT11B_051_03 through AGBT11B_051_104 to the first scan in AGBT11B_051_05
  ;
  ;   GBTIDL -> align, "AGBT11B_051", 3, 104, 5, indir="testDir/", outdir="testMap/"
  ;   ; Align all scans in files AGBT11B_051_03 through AGBT11B_051_104 to the first scan in AGBT11B_051_05. This finds files from directory "testDir/"
  ;   ; and outputs them in directory "testMap/"

  ; check=1
  ; get first integration of the first file to use as vel ref.
  ; set defaults: input directory "newreduced/", output directory "formap/", both input and output file types '.keep.fits"
  if not keyword_set(indir) then indir = "newreduced/"
  if not keyword_set(infiletype) then infiletype = ".keep.fits"
  if not keyword_set(outdir) then outdir = "formap/"
  if not keyword_set(outfiletype) then outfiletype = ".keep.fits"

  infile = make_filepath(id,reffile,directory=indir,filetype=infiletype,format=format)
  filein, infile
  Scans = get_scan_numbers()
  get,scan=Scans[0],int=0,plnum=0
  accum
  old = getdata(0)

  for session_num=sstart,send do begin
  
    ; In case the "make_filepath" function doesn't work -
    ; if session_num lt 10 then begin ;Open files corresponding to their session number
      ; ff=['AGBT11B_051_0',strcompress(string(fix(session_num)),/remove_all),'.keep.fits']
      ; ff=strjoin(ff,'')
    ; endif else begin
      ; ff=['AGBT11B_051_',strcompress(string(fix(session_num)),/remove_all),'.keep.fits']
      ; ff=strjoin(ff,'')
    ; endelse
    ; fr=['formaptest/',ff]
    ; fi=['newreduced/',ff]
    ; filein,strjoin(fi)
    ; fileout, strjoin(fr)

    infile = make_filepath(id,session_num,directory=indir,filetype=infiletype,format=format)
    outfile = make_filepath(id,session_num,directory=outdir,filetype=outfiletype,format=format)
    filein, infile
    fileout, outfile
    
    Scans=get_scan_numbers()
    nintArr=get_my_nintegrations(Scans)
    
    for i=0,n_elements(Scans)-1 do begin
      start_int = 0
      for j=start_int,nintArr[i]-1 do begin
        for p=0,1 do begin
          get,scan=Scans[i],int=j,plnum=p
          vchange=vshift()
          print,session_num,scans[i],j,p,vchange
          gshift,vchange 
          new=getdata()
          if ~(max(new eq old)) then begin ; data is not equal
              keep
          endif else begin
              print,"duplicate data",session_num,scans[i],j,p,vchange
          endelse
          old=new
        endfor
      endfor
    endfor
  endfor
end
