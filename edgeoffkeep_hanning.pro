PRO runedgeoffkeep,id,start_ind,end_ind,sname,velbaseline,tcalpol0,tcalpol1,badscans=badscans,order=order,spw=spw,smooth=smooth,rms=rms,indir=indir,infiletype=infiletype,outdir=outdir,outfiletype=outfiletype

    ; runedgeoffkeep
    ; Usage: runedgeoffkeep, id, start_ind, end_ind, sname, velbaseline, tcalpol0, tcalpol1 [, badscans=badscans ] [, order=order ] [, spw=spw ] [, smooth=smooth ] [, rms=rms ] [, indir=indir ] [, infiletype=infiletype ] [, outdir=outdir ] [, outfiletype=outfileype ]
    ;
    ; Description: Run edgeoffkeep for multiple files in a row, starting from file start_end to file end_ind
    ;   This defaults to input directory "M31-M33-data/" and input filetype ".raw.acs.fits", and output directory "reduced/" and output filetype ".keep.fits"
    ;
    ;   To compile:
    ;       GBTIDL -> .com edgeoffkeep_hanning.pro
    ;
    ; Parameters:
    ; id (required, input, string): The ID of the dataset, beginning with "AGBT" and ending before the underscore,
    ; 	such as "AGBT11B_051" or "AGBT16A_433"
    ; start_ind (required, input, int): The file number with which to start (such as 1 for "AGBT11B_51_01")
    ; end_ind (required, input, int): The file number with which to end (such as 104 for "AGBT11B_051_104"). The function will analyze files start_end through
    ;	end_ind in order.
    ; sname (required, input, string): Only analyze scans pointed at this source (such as "M31", "3C48", etc.). Compatible with bash scripting.
    ; velbaseline (required, input, float array): The velocity range for which a baseline should be taken and subtracted. The values must be in increasing order
    ;   and have an even number of n_elements
    ; tcalpol0 (required, input, float): The Tcal for polarization 0
    ; tcalpol1 (required, input, float): The Tcal for polarization 1
    ; spw (optional, input, int): The ifnum to analyze. Defaults to 0
    ; smooth (optional, input, int): If smooth is not set to 0, the data will be Hanning smoothed. Defaults to 1 (defaults to being Hanning smoothed).
    ; rms (optional, input, float): rms in K. Defaults to 0.101
    ; badscans (optional, input, int array): The "bad" scans to skip when running edgeoffkeep. Note that this will count for *all* the files: if badscans is set
    ;   to [1,2], *all* scan numbers 1 or 2 in the given set of files will be ignored
    ; order (optional, input, int): The order of the polynomial used to fit the baseline. Defaults to 2nd order
    ; indir (optional, input, string): The relative directory path where the data can be found. This
    ; 	requires the last backslash, such as "./reduced/"
    ; infiletype (optional, input, string): The filetype of the desired input data. This requires the first period,
    ; 	such as ".raw.acs.fits"
    ; outdir (optional, input, string): The relative directory path where the data should be output. This
    ; 	requires the last backslash, such as "./reduced/"
    ; outfiletype (optional, input, string): The filetype of the desired output data. This requires the first period,
    ; 	such as ".raw.acs.fits"
    ;
    ; Examples:
    ;    GBTIDL -> runedgeoffkeep, "AGBT11B_051_01",1,104,"M31_M33_01",[-600,-300,300,600],1.73,1.70,badscans=[29],order=3,spw=0,smooth=1

    ; set defaults
    if not keyword_set(indir) then indir="M31-M33-data/"
    if not keyword_set(infiletype) then infiletype=".raw.acs.fits"
    if not keyword_set(outdir) then outdir="reduced/"
    if not keyword_set(outfiletype) then outfiletype=".keep.fits"

    for i=start_ind,end_ind do begin
        if i gt 9 then begin
            myfile = strcompress(indir+id+"_"+string(i)+infiletype,/remove_all)
            outfile = strcompress(outdir+id+"_"+string(i)+outfiletype/remove_all)
        endif else begin
            myfile = strcompress(indir+id+"_0"+string(i)+infiletype,/remove_all)
            outfile = strcompress(outdir+id+"_0"+string(i)+outfiletype,/remove_all)
        endelse

        ; badscans will vary for each file. might need to fix this in a way that's not just hardcoding (az/el > 0.0)

        edgeoffkeep,myfile,sname,velbaseline,tcalpol0,tcalpol1,badscans=badscans,order=order,outfile=outfile,spw=spw,smooth=smooth,rms=rms
    endfor

end

PRO edgeoffkeep,myfile,sname,velbaseline,tcalpol0,tcalpol1,badscans=badscans,keepchans=keepchans,order=order,outfile=outfile,spw=spw,smooth=smooth,rms=rms

    ; edgeoffkeep
    ; Usage: runedgeoffkeep, id, start_ind, end_ind, sname, velbaseline, tcalpol0, tcalpol1 [, badscans=badscans ] [, order=order ] [, spw=spw ] [, smooth=smooth ] [, rms=rms ]
    ;
    ; Description:
    ;   ==============================================================================
    ;   Spencer A. Wolfe ----- 1/14/11
    ;   Program to average and smooth spectra by treating the ends of a map as off source positions.
    ;   Inputs a scan range, as well as Tcal for both polarizations for each observing session.
    ;   Calculates Ta = Tcal*{ON - <OFF>}/{<OFF_calon - OFF_caloff>}
    ;   and then corrects for atmosphere --> Ta* = Ta *
    ;   exp{tau0/sin(el)}/eta  ..... Tcal is found from scal.pro for the OffOn data
    ;   ASSUMPTION: tau=0.01
    ;   1/18/11
    ;   also truncates spectra to desired frequency range (HI for now) and smooths ~5km/s
    ;   1/19/11
    ;   RFI STUFF --- For now, this just does a simple rms test that
    ;               compares the rms of a spectral region to that of the
    ;               theoretical value.  if rms > 3*(rms theo), then the
    ;               data are not added to the final fits file. This takes
    ;               care of most of the integrations with "ripply"
    ;               baselines. Will need to deal with the RFI spikes in th
    ;               OH and recomb lines later...
    ;
    ;   ORIGINAL BY Spencer Wolfe (WVU)
    ;   Heavly Edited by Dominic A. Ludovici (WVU) for processing of LGG140 data.
    ;
    ;   Rewritten to incorporate folding, baselining, and boxcaring into same script 
    ;   in order to improve speed and functionality.  DAL 5/17/2011
    ;   changed all 'get'and 'select' commands to 'gettp' commands to get
    ;   idlToSdfits to work.  7/1/2011 DJP
    ;   provided ability to set fileout and IFnum in inputs as well as simplified baseline inputs
    ;   Also commented out "replace" commands and set units to "velo" before "keep" statements
    ;   (returning to "chan" afterwards).
    ;   1/7/2012 DJP
    ;   Added parameters to allow selection of a given IF (spw) and to pick the number of channels
    ;   to smooth the data (smooth)
    ;   3/20/12 DJP 
    ;   baselining changed based on velocities.  keepchans still in
    ;   channels.  Also setting tcal in !g.s[0] headers
    ;   6/27/16 DJP
    ;   6/5/26 Toney Minter replace smooth with hanning function
    ;   6/9/26 Toney Minter comment out with ";;" unecessary steps for total
    ;   power mapping data with reference integrations at edge of map
    ;   replace break with return in if statements
    ;   comment out sig_state for gettp
    ;   velbaseline in km/s so need to multiply by 1000 (not divide) to get m/s
    ;   flag data with high rms above theoretical
    ;   3 seconds per integration, 2.5 km/s resolution, expected noise
    ;   about 101.4 mK
    ;   7/31/26 (hopefully) clarifying edits made by Bethany Grimm
    ;   To compile:
    ;       GBTIDL -> .com auxiliary.pro
    ;       GBTIDL -> .com edgeoffkeep_hanning.pro
    ;   ===============================================================================
    ;
    ; Parameters:
    ; id (required, input, string): The ID of the dataset, beginning with "AGBT" and ending before the underscore,
    ; 	such as "AGBT11B_051" or "AGBT16A_433"
    ; start_ind (required, input, int): The file number with which to start (such as 1 for "AGBT11B_51_01")
    ; end_ind (required, input, int): The file number with which to end (such as 104 for "AGBT11B_051_104"). The function will analyze files start_end through
    ;	end_ind in order.
    ; sname (required, input, string): Only analyze scans pointed at this source (such as "M31", "3C48", etc.). Compatible with bash scripting.
    ; velbaseline (required, input, float array): The velocity range for which a baseline should be taken and subtracted. The values must be in increasing order
    ;   and have an even number of n_elements
    ; tcalpol0 (required, input, float): The Tcal for polarization 0
    ; tcalpol1 (required, input, float): The Tcal for polarization 1
    ; spw (optional, input, int): The ifnum to analyze. Defaults to 0
    ; smooth (optional, input, int): If smooth is not set to 0, the data will be Hanning smoothed. Defaults to 1 (defaults to being Hanning smoothed).
    ; rms (optional, input, float): rms in K. Defaults to 0.101
    ; badscans (optional, input, int array): The "bad" scans to skip when running edgeoffkeep. Note that this will count for *all* the files: if badscans is set
    ;   to [1,2], *all* scan numbers 1 or 2 in the given set of files will be ignored
    ; order (optional, input, int): The order of the polynomial used to fit the baseline. Defaults to 2nd order
    ;
    ; Examples:
    ;    GBTIDL -> edgeoffkeep,'M31-M33-data/AGBT11B_051_01.raw.acs.fits','M31_M33_01',[-600,-300,300,600],1.73,1.70,badscans=[29],order=3,outfile='reduced/AGBT11B_051_01.keep.fits',spw=0,smooth=1

    ;; run LSRDopplerCorrection on data after edgeoffkeep

    ; open the file
    filein,myfile
    
    ;select all scans of a source that are mapping scans
    allscans=get_scan_numbers(source=sname,procedure='*Map')

    ;remove badscans
    if not keyword_set(badscans) then badscans=get_scan_numbers(azimuth=0.0,elevation=0.0)
    goodscans=setdifference(allscans,badscans)

    if not keyword_set(velbaseline) then begin
        print,'Must provide at least 4 velocities to define baseline region'
        return
    endif

    if not keyword_set(tcalpol0) or not keyword_set(tcalpol1) then begin
        print,'Must provide Tcal values for both polarization'
        return
    endif

    ;set 'if' (spectral window) to process
    if not keyword_set(spw) then spw=0

    ;set output filename
    if not keyword_set(outfile) then begin
        fileout,'edgeoffkeep_out.fits'	
    endif else begin
        fileout,outfile
    endelse	

    ;if smooth ne 0 then we will hanning smooth the data
    if not keyword_set(smooth) then smooth=1

    if not keyword_set(rms) then rms=0.101 ; Kelvin


    chan											   
    freeze
    for hh = 0,n_elements(goodscans)-1 do begin
        h=goodscans(hh) ;this is just to avoid editing all of spencers old code
    ;====================================
    ; RA_Long scans total power and not frequency swich using signal
    ;====================================
        print,'Working on scan: ',h
        get, scan=h, int=0, ifnum=spw, plnum=0, sig='T', cal='T' ; call scan to get the number of intergrations
        info=scan_info(h)
        nint=info.n_integrations
        offs=[0,1,2,3,nint-4,nint-3,nint-2,nint-1] ; set limits for the off scans
        print, "Using:", offs
        for e = 0,n_elements(offs)-1 do begin ; loop through limits and polarizations
            for f = 0,1 do begin
                gettp, h, intnum=offs[e], plnum=f, ifnum=spw, cal_state=1;, sig_state=1 ; get calon for the off for sig=T
                copy, 0, 1 
                gettp, h, intnum=offs[e], plnum=f, ifnum=spw, cal_state=0;, sig_state=1 ; get caloff for the off for sig=T
                copy, 0, 2 
                subtract, 1, 2      ; (calon - caoloff) for the off
                accum, f            ; accumulate for each polarization
    ;;            gettp, h, intnum=offs[e], plnum=f, ifnum=spw, cal_state=1, sig_state=0 ; get calon for the off for sig=F
    ;;            copy, 0, 3 
    ;;            gettp, h, intnum=offs[e], plnum=f, ifnum=spw, cal_state=0, sig_state=0 ; get caloff for the off for sig=F
    ;;            copy, 0, 4 
    ;;            subtract, 3, 4      ; (calon - caoloff) for the off
    ;;            accum, f+2          ; accumulate for each polarization
            endfor ; end for f
        endfor ; end for e
        ave, 0,/quiet               ; <calon - caloff> for pol 0
        data_copy,!g.s[0],X1
        ave, 1,/quiet               ; <calon - caloff> for pol 1
        data_copy,!g.s[0],X2        
    ;;    ave, 2,/quiet
    ;;    data_copy,!g.s[0],X3
    ;;    ave, 3,/quiet
    ;;    data_copy,!g.s[0],X4
        sclear, 0                   ; clear accumulators 1 and 2
        sclear, 1
    ;;    sclear, 2
    ;;    sclear, 3
        for f = 0, 1 do begin       ; loop through polarizations
            sclear
            for e=0,n_elements(offs)-1 do begin
                gettp,h,intnum=offs[e],plnum=f,ifnum=spw;,sig_state=1 ; get the off ints for sig=T
                accum
            endfor  ; end for e
            ave                     ; <off> for sig=T
            ; smooth the offs - reduce their noise by about factor of 3
            boxcar,10 ; do not decimate - need same number channels for on-<off>
            copy, 0, 1+f
            sclear
    ;;        for e=0,n_elements(offs)-1 do begin
    ;;            gettp,h,intnum=offs[e],plnum=f,ifnum=spw,sig_state=0 ; get the off ints for sig=F
    ;;            accum
    ;;        endfor
    ;;        ave                     ; <off> for sig=F
    ;;        copy, 0, 3+f
    ;;        sclear
        endfor  ; end for f
        for g=0,nint-1 do begin    ; loop through the ons and polarization
            for f=0,1 do begin
                gettp,h,intnum=g,plnum=f,ifnum=spw;,sig_state=1 ;cal_on+cal_off/2 for sig=T
                ; print,"Int ",g
            copy, 0, 5+f
    ;;            gettp,h,intnum=g,plnum=f,ifnum=spw,sig_state=0 ;cal_on + cal_off /2 for sig=F
    ;;            copy, 0, 7+f
            endfor ; end for f
            subtract, 5, 1, 9       ; on - <off>
            subtract, 6, 2, 10
    ;;        subtract, 7, 3, 11
    ;;        subtract, 8, 4, 12
            set_data_container,X1
            divide, 9, 0         ; (on - <off>)/<calon - caloff> for pol 0
            scale, tcalpol0*exp(0.01/sin(!g.s[0].elevation*3.14159/180.0))/0.99 ; Ta*
    ;;        copy,0,13               ;sig0
    ;;        set_data_container,X3
    ;;        divide, 11, 0        ; (on - <off>)/<calon - caloff> for pol 0
    ;;        scale, tcalpol0*exp(0.01/sin(!g.s[0].elevation*3.14159/180.0))/0.99 ; Ta*
    ;;        copy,0,14               ;ref0
    ;;        fold,13,14
            !g.s[0].mean_tcal=tcalpol0
        bregion=round(veltochan(!g.s[0],velbaseline*1000.))  ;converts velocities to channels for baselining
            nregion,bregion
            if keyword_set(order) then baseline,nfit=order
            x=dcextract(!g.s[0],min(bregion),max(bregion)) ; extracts spectrum over regions of interest
            set_data_container,x
            data_free,x
            if smooth ne 0 then hanning,/decimate
            !g.s[0].units='Ta*'
            bregion=round(veltochan(!g.s[0],velbaseline*1000.)) ;converts velocities to channels for baselining
            stats,bregion[0]-1,bregion[1]-1,ret=mystats,/quiet
            ;;print,"stats-0",mystats.rms,h,g,!g.s[0].elevation
            if mystats.rms gt 3*rms then begin  
            print,"flag this data",mystats.rms,h,f,g
            endif else begin
            velo
            keep
            endelse
            chan
            set_data_container,X2
            divide, 10, 0
            scale, tcalpol1*exp(0.01/sin(!g.s[0].elevation*3.14159/180.0))/0.99 ; pol 1
    ;;        copy,0,13               ;sig1                                 
    ;;        set_data_container,X4
    ;;        divide, 12, 0
    ;;        scale, tcalpol1*exp(0.01/sin(!g.s[0].elevation*3.14159/180.0))/0.99 ; pol 1
    ;;        copy,0,14               ;ref1
    ;;        fold,13,14
            !g.s[0].mean_tcal=tcalpol1
        bregion=round(veltochan(!g.s[0],velbaseline*1000.))  ;converts velocities to channels for baselining
        nregion,bregion
            if keyword_set(order) then baseline,nfit=order
            x=dcextract(!g.s[0],min(bregion),max(bregion))  ; extracts spectrum over regions of interest
            set_data_container,x
            data_free,x
        if smooth ne 0 then  hanning,/decimate
        !g.s[0].units='Ta*'
            bregion=round(veltochan(!g.s[0],velbaseline*1000.))  ;converts velocities to channels for baselining
            stats,bregion[0]-1,bregion[1]-1,ret=mystats,/quiet
            ;;print,"stats-1",mystats.rms,h,g,!g.s[0].elevation
            if mystats.rms gt 3*rms then begin  
            print,"flag this data",mystats.rms,h,f,g
            endif else begin
            velo
            keep
            endelse
        endfor ; end for g
        data_free,X1
        data_free,X2
    ;;    data_free,X3
    ;;    data_free,X4SW
    endfor  ; end for h

    unfreeze

end