'''
This script takes FITS files, assembles moment-0 maps, masks out noise, and iterates Gaussian fits
in order to find the outflow opening and position angles.
Also user input and iteration are required at many points so this works better as a Jupyter notebook
'''

#imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math
from spectral_cube import SpectralCube
import astropy.units as u
from astropy.utils import data
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs import utils
from astropy.coordinates import SkyCoord
import json
from scipy.optimize import curve_fit
from labels import label_update
from labels import correct_angle

#have the full label list ready
jsonpath = "./json/labels.json"

data = json.loads(open(jsonpath).read())
l_list = data["Data"]

#define arrays and functions
OAr = np.zeros((21,36))
OAb = np.zeros((21,36))
X_OAr = np.zeros((21,36))
X_OAb = np.zeros((21,36))

A  = np.zeros((21,2))
Mu = np.zeros((21,2))
Sig= np.zeros((21,2))
C  = np.zeros((21,2))
Mu_err = np.zeros((21,2))
Sig_err = np.zeros((21,2))

def gaussian(x, A, mu, sig, c):
    '''
    A Gaussian function.

    Inputs:
        x (float or numpy array): x of the function
        A (float): amplitude of the function
        mu (float): mean of the function
        sig (float): width of the function
        c (float): y-shift of the function

    Output:
        output (float or numpy array): resultant Gaussian function
    '''
    return A*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))+c

def transform(c, header):
    '''
    This function transforms RA/Dec coordinates to pixel coordinates.

    Inputs:
        c (astropy coordinate object): RA/Dec coordinates
        header (FITS header object): header of FITS file
    
    Outputs:
        cx (float): x-coordinate of c, in pixels
        cy (float): y-coordinate of c, in pixels
    '''
    ra1,dec1 = header['CRVAL1'], header['CRVAL2']
    d_ra,d_dec = header['CDELT1'], header['CDELT2']
    p_ra1,p_dec1 = header['CRPIX1'], header['CRPIX2']
    ra2,dec2 = c.ra.degree, c.dec.degree
    p_ra2 = p_ra1 + ((ra2-ra1)/d_ra)
    p_dec2 = p_dec1 + ((dec2-dec1)/d_dec)
    return(p_ra2,p_dec2)

def lines(c: bool, blrr: bool, fig, vname: str, cx: float, cy: float, r: int, top: int):
    '''
    This function plots lines on top of an outflow contour plot, showing the outflow position and opening angles.

    Inputs:
        c (bool): True for inclination-corrected angles, and False for non-inclination-corrected angles
        blrr (bool): boolean meaning "blue left red right": True if the blueshifted lobe is centered on the left side
            of the plot and redshifted lobe is on the right, and False if vice versa. Unfortunately, not every outflow
            looks neatly like this: a more sophisticated method might be needed (or a little bit of hardcoding)
        fig (matplotlib figure object): figure upon which to plot the lines
        vname (str): name of source - this is also used as the key for which to obtain label information
        cx (float): x-coordinate of protostellar center, in pixels
        cy (float): y-coordinate of protostellar center, in pixels
        r (int): radius of plot, in pixels
        top (int): number of pixels to "cut off" the top, such that the lines do not overlap with the lables
    
    Outputs:
        fig (matplotlib figure object): figure, but with lines plotted on top
    '''
    data = json.loads(open(jsonpath).read())
    l_list = data["Data"]
    
    op_key = "Opening corrected" if c else "Opening angle"
    try:
        xcb, x1b, x2b, xcr, x1r, x2r = [],[],[],[],[],[]
        xarrays = [xcb, x1b, x2b, xcr, x1r, x2r]
        for i in range(len(xarrays)):
            xarrays[i] = np.linspace(cx-r, cx+r, 4*r)

        target = {}
        for i in l_list:
            if i["Name"] == vname:
                target = i

        #Already multiplied by 3.3302 for full width quarter maximum
        cb = target["Position angle"]["Blue"]
        opb = target[op_key]["Blue"]
        b1 = cb - opb/2
        b2 = cb + opb/2
        cr = target["Position angle"]["Red"]
        opr = target[op_key]["Red"]
        r1 = cr - opr/2
        r2 = cr + opr/2

        ycb, y1b, y2b, ycr, y1r, y2r = [],[],[],[],[],[]
        yarrays = [ycb, y1b, y2b, ycr, y1r, y2r]
        angles = [cb, b1, b2, cr, r1, r2]
        print(angles)

        for i in range(len(yarrays)):
            yarrays[i] = (xarrays[i] - cx) * (math.tan((angles[i]+90)*math.pi/180)) + cy
            yarrays[i] = np.where(yarrays[i]>(cy+(r-top)), None, yarrays[i])
            
            if (angles[i] == cb) or (angles[i] == cr):
                ls = "--"
                w = 0.5
            else:
                ls = "-"
                w = 0.8
            
            # what if opening angle changes sign?
            # what if abs(angle) < 1?

            # if left, xarrays[1]>cx becomes 0
            # this is true if blue side & blrr & not changed sign
            # or red side & blrr & changed sign
            # or red side & not blrr & not changed sign
            # or blue side & not blrr & changed sign
            # also in this case, add 1 if abs(angle[i]) < 0
            
            th = cx
            if(angles[i]==0):
                continue
            
            elif((blrr and (i<3) and (angles[i]/angles[0] > 0)) or
               (blrr and (i>2) and (angles[i]/angles[3] < 0)) or
               ((not blrr) and (i>2) and (angles[i]/angles[3] > 0)) or
               ((not blrr) and (i<3) and (angles[i]/angles[0] < 0)) or
               ((not blrr) and (i<3) and (angles[i] > 180))):
                print(str(i) + " left")
                if (abs(angles[i]) < 1):
                    th += 1
                xarrays[i] = np.where(xarrays[i]>cx, None, xarrays[i])
            
            else:
                print(str(i) + " right")
                if (abs(angles[i]) < 1):
                    th -= 1
                xarrays[i] = np.where(xarrays[i]<(cx), None, xarrays[i])
            
            plt.plot(xarrays[i],yarrays[i],linestyle=ls,c="black",lw=w)
            
        # This function has issues plotting vertical lines. The following lines of code plot a vertical line
        # through the center
        # xv = np.linspace(cx, cx, r)
        # yv = np.linspace(cy-r, cy, r)
        # plt.plot(xv, yv, linestyle="-",c="black",lw=0.8)

    except IndexError:
        pass
    return fig

def find_range(PA, a_guess, cutoffs):
    '''
    The angles-from center of every pixel in the masked moment-0 map of the outflows are assembled into a histogram.
    Taking a Gaussian fit identifies the mean outflow position angle (center) and outflow opening angle (FWHM), but
    fitting the Gaussian is difficult if there are few or missing pixels.
    This function is to help find a range for which the histogram Gaussian should be fit.

    Inputs:
        PA (array of floats): the bin values for the histogram
        a_guess (2-element array of floats): initial guess for the range, in degrees
        cutoffs (2 or 0-element array of floats): if the histogram requires two separate Gaussian fits, this specifies
            at what value (in degrees) to end the first range and start the second range. if the histogram only needs 
            one Gaussian fit, this can be passed as an empty array.

    Outputs:
        edge1 (bool): returns True if the range falls on the leftmost edge of the histogram
        edge2 (bool): returns True if the range falls on the rightmost edge of the histogram
        e1 (int): bin index at which range starts
        e2 (int): bin index at which range ends
        c1 (int): if two Gaussian fits, bin index at which first range ends
        c2 (int): if two Gaussian fits, bin index at which second range starts
    '''
    edge1 = False
    edge2 = False
    c1 = 0
    c2 = 0
    e1 = np.argmax(PA > a_guess[0])
    if e1 == 0:
        edge1 = True
    e2 = np.argmax(PA > a_guess[1])
    if e2 == 0:
        edge2 = True
        e2 = len(PA) - 1
    if len(cutoffs) == 2:
        c1 = np.argmax(PA > cutoffs[0])
        c2 = np.argmax(PA > cutoffs[1])
    return (edge1, edge2, e1, e2, c1, c2)

def sample_fit_plots(a_guess,cutoffs,binsize,X,shift,p0_i):
    '''
    This function iteratively fits and plots Gaussians, slightly varying edge parameters each time.

    Inputs:
        a_guess (2-element array of floats): predicted range of angles the outflow spans
        cutoffs (2 or 0-element array of floats): predicted range of angles with shallow/mising outflow
        binsize (int): we iterate over binsizes of 3, 4, and 5
        X (array of ints): histogram of angles, as an array
        shift (float): amount in degrees histogram is to be shifted
        p0_i (array of floats): guess for Gaussian fit - record this for reproducibility
            [amplitude, mu, sigma, y-shift]

    Outputs:
        mus (array of floats): array of mu (center) values obtained from each Gaussian
        sigmas (array of floats): array of sigma (width parameter) values obtained from each Gaussian
        Gaussian plots will be displayed as output as well
    '''
    binrange = np.arange((-92.5+shift),(92.5+shift),binsize)
    PA = np.arange((-90+shift),(90+shift),binsize)

    i_it = []
    mus = []
    sigmas = []
    fig1 = plt.figure()
    (edge1, edge2, e1, e2, c1, c2) = find_range(PA, a_guess, cutoffs)
    lsh = 1
    ls = lsh * 2 + 1
    side1 = []
    side2 = []
    for i in range(ls):
        if edge1 and i <= lsh:
            side1.append(0)
        else:
            side1.append(e1 + i - lsh - 1)
        if edge2 and i > (lsh + 1):
            side2.append(len(PA) - 1)
        else:
            side2.append(e2 + i - lsh - 1)
    i_it.append(side1)
    i_it.append(side2)
    # print(i_it)

    for i in range(ls):
        for j in range(ls):
            if len(cutoffs) == 2:
                ii = np.concatenate((np.arange(i_it[0][i],c1), np.arange(c2,i_it[1][j])))
            else:
                ii = np.arange(i_it[0][i],i_it[1][j])
            ax = plt.subplot(ls,ls,(j*ls+i+1))
            output = ax.hist(X.flatten(),bins=binrange)
            # print(output)
            target_func = gaussian
            popt, pcov = curve_fit(target_func,PA[ii],output[0][ii],p0=p0_i, maxfev=100000)
            x = np.linspace((-90+shift),(90+shift),180)
            mu_i = popt[1]
            if mu_i > 90:
                mu_i1 = mu_i-180
            else:
                mu_i1 = mu_i
            sigma_i = popt[2]
            ax.plot(x, gaussian(x, *popt), 'g--',label=('fit: $\mu$='+str(round(mu_i1,2))+', $\sigma$='+str(round(sigma_i,2))))
            plt.xticks([])
            plt.yticks([])
            mus.append(float(mu_i))
            sigmas.append(float(abs(sigma_i)))
            #manually check the best (default init guess) and use that to derive error
            #then update
            #repeat for red. at the moment this will only work for outflows that can be analyzed with one gaussian

    return [mus, sigmas]

def calc_angles(mus,sigmas,valid):
    '''
    This function calculates position angle, opening angle, and errors for both.

    Inputs:
        mus (array of floats): mu values found from Gaussian fits, output from sample_fit_plots
        sigmas (array of floats): sigma values found from Gaussian fits, output from sample_fit_plots
        valid (array of ints): indices for which the fit is valid (and not wildly inaccurate)

    Outpus:
        pa (float): position angle of outflow, in degrees
        oa (float): opening angle of outflow, in degrees
        epa (float): error of position angle, in degrees
        eoa (float): error of opening angle, in degrees
    '''
    mus = np.ndarray.flatten(np.array(mus))[valid]
    sigmas = np.ndarray.flatten(np.array(sigmas))[valid]
    pa = np.average(mus)
    oa = np.average(sigmas) * 3.3302
    epa = (np.max(mus) - np.min(mus)) * 0.5
    eoa = (np.max(sigmas) - np.min(sigmas)) * 3.3302 * 0.5
    return (pa, oa, epa, eoa)

#get data for both blue and red sides
#blue is higher frequency

#recreate the moment maps and define parameters for blue and red sides
filename = "HOPS-347_12CO_image_taper500k.pbcor.fits" #input file name
bparams = (3.45793e+011, 3.45787e+011) #input frequency range in Hz for blueshifted lobe
rparams = (3.45782e+011, 3.45775e+011) #input frequency range in Hz for redshifted lobe: should be lower than bparams
v = "" #if this is a binary source, replace "" with "-A" or the appropriate suffix
filepath = "./downloads/" #directory where to find the downloaded files
cubepath = "./cubes/" #directory where to save the asembled datacubes
maskpath = "./anglemask/" #directory where to save the masked outflows with angle information
histpath = "./histograms/" #directory where to save angle histograms
contpath = "./contours/" #directory where to save contour plots


#Creating the datacubes
if True:
    pathname = filepath + filename
    name = filename.split("_")[0].upper()
    vname = name + v
    fitsname = vname + ".fits"
    f = fits.open(pathname)
    header = f[0].header

    cube = SpectralCube.read(pathname).with_spectral_unit(u.Hz)

    #blue side
    (b1, b2) = bparams
    bcube = cube.spectral_slab(b1 * u.Hz, b2 * u.Hz)
    bnii_cube = bcube.with_spectral_unit(u.km/u.s,
                                        velocity_convention='radio',
                                        rest_value=345.796* u.GHz)
    bmoment = bnii_cube.moment(order=0)
    # bmoment = bcube.moment(order=0)
    bsavepath = cubepath + "b" + fitsname
    bmoment.write(bsavepath, overwrite = True)

    #red side
    (r1, r2) = rparams
    rcube = cube.spectral_slab(r1 * u.Hz, r2 * u.Hz)
    rnii_cube = rcube.with_spectral_unit(u.km/u.s,
                                        velocity_convention='radio',
                                        rest_value=345.796* u.GHz)
    rmoment = rnii_cube.moment(order=0)
    rsavepath = cubepath + "r" + fitsname
    rmoment.write(rsavepath, overwrite = True)

    image_datab = fits.getdata(bsavepath,ext=0)
    image_datar = fits.getdata(rsavepath,ext=0)

#Plot blueshifted and redshifted moment maps
fig = plt.figure(figsize = (16,4))

ax = plt.subplot(1,3,1)
plt.imshow(image_datab,origin='lower')
plt.clim(vmin=0)
plt.colorbar()

ax = plt.subplot(1,3,2)
plt.imshow(image_datar,origin='lower')
plt.clim(vmin=0)
plt.colorbar()

#This one shows a composite image of both lobes
image_data = image_datar+image_datab
ax = plt.subplot(1,3,3)
plt.imshow(image_data,origin='lower')
plt.clim(vmin=0)
plt.colorbar()

plt.show()

#Set noise value, based on number of integration slices. This may be manually adjusted if needed
noise = (((bparams[0]-bparams[1])+(rparams[0]-rparams[1]))/1000000 + 2) * 0.02
# noise = 0.3 #If auto-calculated noise doesn't do, try this
print(noise)

#First contour plot
nanmax = np.nanmax(image_datab+image_datar)
steps = np.arange(noise, nanmax, noise)
plt.contour(image_datab, steps, colors='b')
plt.contour(image_datar, steps, colors='r')
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')

#Radius of source, in pixels: this can be guessed first, and then adjusted to see what best masks out fringe noise
radb,radr = 100,170
r = max(radb,radr)

#Find center from FITS file
#Find coordinates from Table 6 
no_coord = True
c=""
for i in l_list:
    try:
        #Change vname to actual index if inaccurate
        if i["Name"] == vname:
            c = i["RADEC"]
            no_coord = False
    except KeyError:
        pass
if (c=="" or no_coord):
    #May set coordinates manually, if they're not in the header
    cx,cy = 0,0

else:
    c = SkyCoord(c, unit=(u.hourangle,u.degree))
    cx,cy = transform(c,header)
    # cx,cy = 481.32,469.8 #May set coordinates manually, if they're not in the header

#Plotting continues here
#This shows a preview of the masked moment-0 map
if True:
    pb = patches.Circle((cx, cy), radb, fill=False)
    pr = patches.Circle((cx, cy), radr, fill=False)

    #this masks both plots
    maskb = np.where(image_datab > noise,1.0,0.0 )
    maskr = np.where(image_datar > noise,1.0,0.0 )

    #show the masked plots
    fig = plt.figure(figsize = (12,5))

    ax = plt.subplot(1,2,1)
    plt.imshow(maskb,origin='lower')
    ax = plt.gca()
    ax.add_patch(pb)
    ax.plot(cx, cy, "k.")
    plt.colorbar()

    ax = plt.subplot(1,2,2)
    plt.imshow(maskr,origin='lower')
    ax = plt.gca()
    ax.add_patch(pr)
    ax.plot(cx, cy, "k.")
    plt.colorbar()

    plt.show()

#This segment is only needed if extra areas need to be masked out

extra_masking = True

if extra_masking:
    #note that i corresponds to y, and j with x
    #masks out circular patches on the blueshifted lobe
    exb = [] #y-coord of patch center
    eyb = [] #x-coord of patch center
    erb = [] #radius of patch center
    #masks out circular patches on the redshifted lobe
    exr = [] #y-coord of patch center
    eyr = [] #x-coord of patch center
    err = [] #radius of patch center
    
    for i in range(len(image_datab)):
        for j in range(len(image_datab[0])):
            for k in range(len(exb)):
                rade = math.sqrt(math.pow((i-eyb[k]), 2) + math.pow((j-exb[k]), 2))
                if rade < erb[k]:
                    maskb[i][j] = 0
                    image_datab[i][j] = 0
            for k in range(len(exb)):
                rade = math.sqrt(math.pow((i-eyr[k]), 2) + math.pow((j-exr[k]), 2))
                if rade < err[k]:
                    maskr[i][j] = 0
                    image_datar[i][j] = 0

            #more masking parameters can be manually input here for either red or blue
            #good if an entire region above/below or left/right of a certain line can be scrapped
            # if j < 475 or i < 475:
            #     maskb[i][j] = 0
            #     image_datab[i][j] = 0
            # if i > 550 or j > 570 or j < 400:
            #     maskr[i][j] = 0
            #     image_datar[i][j] = 0

#now, mask out everything outside the noise radius
for i in range(len(image_datab)):
    for j in range(len(image_datab[0])):
        rad = math.sqrt(math.pow((i-cy), 2) + math.pow((j-cx), 2))
        if rad > radb:
            maskb[i][j] = 0
            image_datab[i][j] = 0
        if rad > radr:
            maskr[i][j] = 0
            image_datar[i][j] = 0

#show the masked plots
fig = plt.figure(figsize = (12,5))
ax = plt.subplot(1,2,1)
plt.imshow(maskb,origin='lower')
plt.colorbar()
ax = plt.subplot(1,2,2)
plt.imshow(maskr,origin='lower')
plt.colorbar()
plt.show()

#does the histogram need to be shifted by a certain number of pixels to accurately show the whole shape?
#since the default cutoffs are at -90 and 90 degrees, we don't want a histogram that starts to peak
#around 90, cuts, and continues on the other side of the plot at -90 (can't fit the Gaussian that way)
shift = 0

if True:
    image_header = fits.getheader(pathname, ext=0)

    PA_grid = np.full((image_header['NAXIS1'],image_header['NAXIS2']),0.0)

    #calculate the angles from the center
    for i in range(PA_grid.shape[0]):
        for j in range(PA_grid.shape[1]):
            x = i-cy
            y = j-cx
            if x ==0:
                PA_grid[i][j] = -np.arctan(np.inf)/np.pi*180
            else:
                val = -np.arctan(y/x)/np.pi*180
                if(val <= -90+shift):
                    val += 180
                PA_grid[i][j] = val

    maskb[maskb==0]=np.nan
    maskr[maskr==0]=np.nan

    #Show the angle distribution on the mask as a colorbar, then save this
    fig = plt.figure(figsize = (13,5))
    ax = plt.subplot(1,2,1)
    B = PA_grid*maskb
    plt.imshow(B,origin='lower')
    plt.colorbar(label="Position Angle")
    ax.set_xlim(cx-(r+10), cx+(r+10))
    ax.set_ylim(cy-(r+10), cy+(r+10))
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_yticklabels([])
    ax.set_yticks([])
    ax.set_xlabel("RA")
    ax.set_ylabel("DEC")
    ax = plt.subplot(1,2,2)
    R = PA_grid*maskr
    plt.imshow(R,origin='lower')
    plt.colorbar(label="Position Angle")
    ax.set_xlim(cx-(r+10), cx+(r+10))
    ax.set_ylim(cy-(r+10), cy+(r+10))
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_yticklabels([])
    ax.set_yticks([])
    ax.set_xlabel("RA")
    ax.set_ylabel("DEC")
    fig.suptitle(vname)
    plt.show()
    anglemaskpath = maskpath + vname + ".png"
    fig.savefig(anglemaskpath)

# May define bin sizes manually: 3, 4, and 5 give the best results
binsizes = [3,4,5]

#Let's stick with 5 just to get an idea-
binrange = np.arange((-92.5+shift),(87.5+shift),5)
# brb = binrange
# brr = binrange
# brb = np.arange((85-2.5), (130-2.5), 3)
# brr = np.arange((25-2.5), (90-2.5), 5)

#Show angle histograms for blueshifted and redshifted lobes
fig = plt.figure(figsize = (12,5))
ax = plt.subplot(1,2,1)
ax.grid()
outputb = plt.hist(B.flatten(),bins=binrange)
ax = plt.subplot(1,2,2)
ax.grid()
outputr = plt.hist(R.flatten(),bins=binrange)
plt.show()

#Initial guesses

histpathb = histpath + vname + "_opening_angle_blue_hist.png"
histpathr = histpath + vname + "_opening_angle_red_hist.png"

#Initial guesses for fit parameters: change these to refine the fit
#[amplitude, center, width, y-shift]
p0_b=[100,-65,5,0]
p0_r=[800,-20,5,0]

#If either the blue or red lobe need to be "split", that is two Gaussians computed for two peaks of the
#same histogram, change to True. If one fit will do, change to False
bsplit = True
rsplit = True

#Only if two fits are anticipated: initial guess for the second (rightmost) Gaussian
#[amplitude, center, width, y-shift]
p0_b2 = [50,-35,5,0]
p0_r2 = [200,30,5,0]

#Specify at which angle values to start and end the fit
#Noise in angle ranges not of interest may interfere
ab_guess = [-90,-25]
ar_guess = [-90,90]

#Again, if two fits are anticipated for either histogram, specify where to end the first Gaussian
#range and begin the second
cutoffs_b = [-50,-60]
cutoffs_r = [0,20]

#Iterate through many candidate fits

mus_b = []
sigmas_b = []
mus_r = []
sigmas_r = []

#Case 1: single Gaussian fit
#Iterate fits over the blueshifted lobe, and record values
if not bsplit:
    for b in binsizes:
        [mus_bi, sigmas_bi] = sample_fit_plots(ab_guess,cutoffs_b,b,B,shift,p0_b)
        mus_b.append(mus_bi)
        sigmas_b.append(sigmas_bi)
    #identify best overall to find error

#Repeat, iterating fits over the redshifted lobe, and record values
if not rsplit:
    for b in binsizes:
        [mus_ri, sigmas_ri] = sample_fit_plots(ar_guess,cutoffs_r,b,R,shift,p0_r)
        mus_r.append(mus_ri)
        sigmas_r.append(sigmas_ri)

#Case 2: two Gaussian fits. this gets a little more complicated, with more arrays to track
#But essentially, it's the same as earlier: just do everything twice
mus_b1 = []
sigmas_b1 = []
mus_b2 = []
sigmas_b2 = []
if bsplit:
    ab_guess1 = [ab_guess[0],cutoffs_b[0]]
    ab_guess2 = [cutoffs_b[1],ab_guess[1]]
    for b in binsizes:
        [mus_bi, sigmas_bi] = sample_fit_plots(ab_guess1,[],b,B,shift,p0_b)
        mus_b1.append(mus_bi)
        sigmas_b1.append(sigmas_bi)
    for b in binsizes:
        [mus_bi, sigmas_bi] = sample_fit_plots(ab_guess2,[],b,B,shift,p0_b2)
        mus_b2.append(mus_bi)
        sigmas_b2.append(sigmas_bi)

mus_r1 = []
sigmas_r1 = []
mus_r2 = []
sigmas_r2 = []
if rsplit:
    ar_guess1 = [ar_guess[0],cutoffs_r[0]]
    ar_guess2 = [cutoffs_r[1],ar_guess[1]]
    for b in binsizes:
        [mus_ri, sigmas_ri] = sample_fit_plots(ar_guess1,[],b,R,shift,p0_r)
        mus_r1.append(mus_ri)
        sigmas_r1.append(sigmas_ri)
    for b in binsizes:
        [mus_ri, sigmas_ri] = sample_fit_plots(ar_guess2,[],b,R,shift,p0_r2)
        mus_r2.append(mus_ri)
        sigmas_r2.append(sigmas_ri)

#Here's where to filter out poorly-fitting fits, if necessary
#numpy arrays of all the indices where all the fits were valid (and not wildly inaccurate)
valid_b = np.arange(0,27)
valid_r = np.arange(0,27)
valid_b2 = np.array([3,6,25])
valid_r2 = np.arange(0,27)

#Angle calculations happen here
if True:
    pa_b = 0
    oa_b = 0
    epa_b = 0
    eoa_b = 0
    if not bsplit:
        (pa_b, oa_b, epa_b, eoa_b) = calc_angles(mus_b, sigmas_b, valid_b)
    pa_r = 0
    oa_r = 0
    epa_r = 0
    eoa_r = 0
    if not rsplit:
        (pa_r, oa_r, epa_r, eoa_r) = calc_angles(mus_r, sigmas_r, valid_r)

    pa_b1 = 0
    oa_b1 = 0
    pa_b2 = 0
    oa_b2 = 0
    if bsplit:
        (pa_b1, oa_b1, epa_b1, eoa_b1) = calc_angles(mus_b1, sigmas_b1, valid_b)
        (pa_b2, oa_b2, epa_b2, eoa_b2) = calc_angles(mus_b2, sigmas_b2, valid_b2)
        pa_b = (pa_b1 + pa_b2)/2
        oa_b = abs(pa_b1 - pa_b2)
        epa_b = np.sqrt(epa_b1**2 + epa_b2**2)/2
        eoa_b = np.sqrt(epa_b1**2 + epa_b2**2)
    pa_r1 = 0
    oa_r1 = 0
    pa_r2 = 0
    oa_r2 = 0
    if rsplit:
        (pa_r1, oa_r1, epa_r1, eoa_r1) = calc_angles(mus_r1, sigmas_r1, valid_r)
        (pa_r2, oa_r2, epa_r2, eoa_r2) = calc_angles(mus_r2, sigmas_r2, valid_r2)
        pa_r = (pa_r1 + pa_r2)/2
        oa_r = abs(pa_r1 - pa_r2)
        epa_r = np.sqrt(epa_r1**2 + epa_r2**2)/2
        eoa_r = np.sqrt(epa_r1**2 + epa_r2**2)

    #for unsplit outflows
    #Position angle: average of mus
    #Opening angle: average of sigmas * 3.3302
    #PA error: greatest difference in mus
    #OA error: greatest difference in sigmas * 3.3302

#Update json database with newly found angles
label_update(vname, pa_b, pa_r, epa_b, epa_r, oa_b, oa_r, eoa_b, eoa_r)
#Apply inclination correction too
correct_angle(vname)

#Plot histograms with best found fit
#For reasons I don't remember this requires the fit to be found yet again?
#But it varies very little from the previously found best fit, if at all

binrange = np.arange((-92.5+shift),(87.5+shift),5)
x = np.linspace((-90+shift),(90+shift),180)

(edge1, edge2, e1, e2, c1, c2) = find_range(binrange, ab_guess, cutoffs_b)
if (len(cutoffs_b) != 0):
    iib = np.concatenate((np.arange(e1,c1), np.arange(c2,e2)))
else:
    iib = np.arange(e1, e2)

plt.figure()
plt.hist(B.flatten(),bins=binrange)
s_b = oa_b / 3.3302
if not bsplit:
    p0_b_new = [p0_b[0], pa_b, s_b, p0_b[3]]
    popt, pcov = curve_fit(gaussian,binrange[iib],outputb[0][iib],p0=p0_b_new, bounds=((-np.inf, pa_b-0.001, s_b-0.001, -np.inf), (np.inf, pa_b+0.001, s_b+0.001, np.inf)), maxfev=100000)
    # print(popt)
    # print(pcov)
    plt.plot(x, gaussian(x, popt[0], pa_b, s_b, popt[3]), 'g--',label=('fit: $\mu$='+str(round(pa_b,2))+', $\sigma$='+str(round(s_b,2))))
if bsplit:
    s_b1 = oa_b1 / 3.3302
    s_b2 = oa_b2 / 3.3302
    iib = np.arange(e1,c1)
    p0_b_new = [p0_b[0], pa_b1, s_b1, p0_b[3]]
    popt, pcov = curve_fit(gaussian,binrange[iib],outputb[0][iib],p0=p0_b_new, bounds=((-np.inf, pa_b1-0.001, s_b1-0.001, -np.inf), (np.inf, pa_b1+0.001, s_b1+0.001, np.inf)), maxfev=100000)
    plt.plot(x, gaussian(x, popt[0], pa_b1, s_b1, popt[3]), 'g--',label=('side 1: $\mu$='+str(round(pa_b1,2))))
    iib = np.arange(c2,e2)
    p0_b_new = [p0_b2[0], pa_b2, s_b2, p0_b2[3]]
    popt, pcov = curve_fit(gaussian,binrange[iib],outputb[0][iib],p0=p0_b_new, bounds=((-np.inf, pa_b2-0.001, s_b2-0.001, -np.inf), (np.inf, pa_b2+0.001, s_b2+0.001, np.inf)), maxfev=100000)
    plt.plot(x, gaussian(x, popt[0], pa_b2, s_b2, popt[3]), 'g--',label=('side 2: $\mu$='+str(round(pa_b2,2))))
    plt.plot(np.linspace(pa_b,pa_b,10),np.linspace(0,np.max(outputb[0]),10),'k:',label=('center: $\mu$='+str(round(pa_b,2))+', $\sigma$='+str(round(s_b,2))))
plt.legend()
plt.xlabel("PA Angle (Degrees)")
plt.ylabel("pixel count")
plt.title(vname + " blue")
plt.savefig(histpathb,dpi=300)
plt.show()
plt.close()

(edge1, edge2, e1, e2, c1, c2) = find_range(binrange, ar_guess, cutoffs_r)
if (len(cutoffs_r) != 0):
    iir = np.concatenate((np.arange(e1,c1), np.arange(c2,e2)))
else:
    iir = np.arange(e1, e2)

plt.figure()
plt.hist(R.flatten(),bins=binrange)
s_r = oa_r / 3.3302
if not rsplit:
    p0_r_new = [p0_r[0], pa_r, s_r, p0_r[3]]
    popt, pcov = curve_fit(gaussian,binrange[iir],outputr[0][iir],p0=p0_r_new, bounds=((-np.inf, pa_r-0.001, s_r-0.001, -np.inf), (np.inf, pa_r+0.001, s_r+0.001, np.inf)), maxfev=100000)
    # print(popt)
    # print(pcov)
    plt.plot(x, gaussian(x, popt[0], pa_r, s_r, popt[3]), 'g--',label=('fit: $\mu$='+str(round(pa_r,2))+', $\sigma$='+str(round(s_r,2))))
if rsplit:
    s_r1 = oa_r1 / 3.3302
    s_r2 = oa_r2 / 3.3302
    iir = np.arange(e1,c1)
    p0_r_new = [p0_r[0], pa_r1, s_r1, p0_r[3]]
    popt, pcov = curve_fit(gaussian,binrange[iir],outputr[0][iir],p0=p0_r_new, bounds=((-np.inf, pa_r1-0.001, s_r1-0.001, -np.inf), (np.inf, pa_r1+0.001, s_r1+0.001, np.inf)), maxfev=100000)
    plt.plot(x, gaussian(x, popt[0], pa_r1, s_r1, popt[3]), 'g--',label=('side 1: $\mu$='+str(round(pa_r1,2))))
    iir = np.arange(c2,e2)
    p0_r_new = [p0_r2[0], pa_r2, s_r2, p0_r2[3]]
    popt, pcov = curve_fit(gaussian,binrange[iir],outputr[0][iir],p0=p0_r_new, bounds=((-np.inf, pa_r2-0.001, s_r2-0.001, -np.inf), (np.inf, pa_r2+0.001, s_r2+0.001, np.inf)), maxfev=100000)
    plt.plot(x, gaussian(x, popt[0], pa_r2, s_r2, popt[3]), 'g--',label=('side 2: $\mu$='+str(round(pa_r2,2))))
    plt.plot(np.linspace(pa_r,pa_r,10),np.linspace(0,np.max(outputr[0]),10),'k:',label=('center: $\mu$='+str(round(pa_r,2))+', $\sigma$='+str(round(s_r,2))))
plt.legend()
plt.xlabel("PA Angle (Degrees)")
plt.ylabel("pixel count")
plt.title(vname + " red")
plt.savefig(histpathr,dpi=300)
plt.show()
plt.close()

#Finally, create the contour plot

r = max(radb, radr)
r += 0
space = r/10
top = 1 if r > 130 else 0
sbl = 5 if r > 150 else 0

#"blue left red right": True if the blueshifted lobe is centered on the left side of the plot and 
# redshifted lobe is on the right, and False if vice versa. Unfortunately, not every outflow
# looks neatly like this: a more sophisticated method might be needed (or a little bit of hardcoding)
blrr = False

hdu = fits.open(pathname)[0]
wcs = WCS(hdu.header)
contourpath = contpath + vname + ".png"
correctpath = contpath + vname + "-corrected.png"

def contour():
    fig, ax = plt.subplots(subplot_kw=dict(projection=wcs,slices=('x','y',0,0)))
    plt.contour(image_datab, steps, colors='b')
    plt.contour(image_datar, steps, colors='r')
    plt.xlabel("RA (ICRS)")
    plt.ylabel("DEC (ICRS)")

    #Add labels
    data = json.loads(open(jsonpath).read())
    l_list = data["Data"]
    no_data = False
    bmaj = fits.getval(pathname, 'BMAJ', ext=0)
    bmin = fits.getval(pathname, 'BMIN', ext=0)
    bpa = fits.getval(pathname, 'BPA', ext=0)
    #Extract labels from l_list
    #May have multiple values for dPa and e_Pa (deg)
    indices = []
    for j in range(len(l_list)):
        list_name = l_list[j]['Name']
        ind = list_name.rfind(name)
        if ind == -1:
            continue
        try:
            #some are binary sources and have a dash after the name
            if (list_name[ind+len(name)] == "-"):
                indices.append(j)
        except IndexError:
            #some are single sources
            indices.append(j)
    # print(indices)
    if (len(indices) == 0):
        no_data = True

    if (not no_data):
        #More than one index indicates a multiple source: print dPa and e_Pa for each star in source
        d0 = l_list[indices[0]] #Dictionary corresponding to first instance of source, will have source's name, class, Lbol, and Tbol
        # print(d0)
        label_1 = "Class "
        label_1 += d0['Class']
        try:
            Lbol = d0['Lbol']
            label_1 += ", $L_{bol}="
            label_1 += Lbol
            label_1 += "L_\odot$"
        except KeyError:
            pass
        try:
            Tbol = d0['Tbol']
            label_1 += ", $T_{bol}="
            label_1 += Tbol
            label_1 += "$K"
        except KeyError:
            pass
        ax.annotate(label_1, (cx-(r),cy+(r-(top*10))))

        #In case of multiple sources, update multiple labels accordingly
        for i in range(len(indices)):
            label_i = ""
            di = l_list[indices[i]] #Dictionary with source's data
            # print(di)
            if(len(indices) > 1):
                source = di['Name'][len(name)+1:]
                label_i += "Source "
                label_i += source
                label_i += ": "

            try:
                dPa = di['dPa']
                label_i += "$PA_{disk}="
                label_i += dPa
            except KeyError:
                label_i += "$"
            try:
                error = di['e_pa (deg)']
                label_i += " \pm "
                label_i += error
                label_i += "^\circ$"
            except KeyError:
                label_i += "^\circ$"
            #relics of hardcoding to make the labels look nice
            # if small:
            #     plt.annotate(label_i, loc=(0.5, (0.9-0.03*i)), fontsize='xx-small')
            # else:
            #     plt.legend(label_i, loc=(0.5, (0.9-0.05*i)), fontsize='small')
            ax.annotate(label_i, (cx-(r), cy+(r-(space*(i+1+top)))))

            #Annotate sources
            try:
                c = di['RADEC']
                c = SkyCoord(c, unit=(u.hourangle,u.deg))
                px, py = transform(c,header)
                plt.plot(px,py,'k+')
                #relics of hardcoding to make the labels look nice
                # if(source=="A"):
                #     px -= 5
                #     py -= 5
                ax.annotate(source, (px+5,py-5))
            except ValueError:
                pass
            except KeyError:
                pass
            except UnboundLocalError:
                pass

    #Add beam size ellipse
    scale = abs(header['CDELT1'])
    bmaj_p, bmin_p = bmaj/scale, bmin/scale
    ax.add_patch(patches.Ellipse((cx+(r-5), cy-(r-5)), bmaj_p, bmin_p, angle=bpa, color='black'))

    #Add scalebar
    #1 arcsec = 400 au = 1/3600 degrees
    #if 400 au exceeds 2r-20 pixels, try 200 au instead
    #also this was never a problem ever
    au400 = True
    width = abs(1/(3600*scale))
    h = 2
    if r > 100:
        h = 4
    if r > 140:
        h = 6
    if r > 200:
        h = 10
    if (width > (2*r - 20)):
        width /= 2
        au400 = False
    ax.add_patch(patches.Rectangle((cx-(r-10), cy-(r-10)), width, h, edgecolor="black", facecolor="black"))
    text="400au" if au400 else "200au"
    ax.annotate(text, (cx-(r-10), cy-(r+sbl)), backgroundcolor = '1')

    return fig

t = 20 #how many pixels to take off the "top" of the lines, so that they don't overlap with the labels

#If angle information exists, plot lines where the outflows are
#Dashed for center outflow line, solid for opening angle lines
fig = contour()
f1 = lines(False, blrr, fig, vname, cx, cy, r, t)
#Zoom in
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
ax.set_xlim(cx-(r+10), cx+(r+10))
ax.set_ylim(cy-(r+10), cy+(r+10))
plt.title(vname)
plt.savefig(contourpath,dpi=300)

#If inclination corrected angles exist, plot lines on a new graph
try:
    fig = contour()
    fc = lines(True, blrr, fig, vname, cx, cy, r, t)
    #Zoom in
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(cx-(r+10), cx+(r+10))
    ax.set_ylim(cy-(r+10), cy+(r+10))
    plt.title(vname + " with inclination correction")
    plt.savefig(correctpath,dpi=300)
except KeyError:
    pass