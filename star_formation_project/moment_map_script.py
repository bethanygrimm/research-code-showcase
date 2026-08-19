'''
This script creates moment-0 maps (with a user-defined frequency range) from FITS files and saves them to a directory.
It also displays moment-0 maps as a PNG and saves these to a directory.
Currently accommodates display of four images at a time
'''

#imports
import numpy as np
from spectral_cube import SpectralCube
from astropy import log, units as u
from astropy.io.fits import getval
from labels import label_list
import aplpy
import matplotlib as plt
import json
import csv
log.setLevel('ERROR')

#How many plots to display in one figure? This is kind of hardcoded right now
n_plots = 4

figname = "Figures_1.png" #Output figure name
#Input FITS file names
fitsdir = "./downloads/" #Relative path of directory where FITS files can be found
cubedir = "./moment0s/" #Relative path where moment-0 files are to be saved
figdir = "./figures/" #Relative path where moment-0 images are to be saved
sources=["HH111mms_12CO_image_taper500k.pbcor.fits", 
         "HH212mms_12CO_image_taper500k.pbcor.fits", 
         "HH270mms1_12CO_image_taper500k.pbcor.fits", 
         "HH270mms2_12CO_image_taper500k.pbcor.fits"]
#Manually input frequency range, in Hz, to use for the moment-0 map
params=[(3.45798e+011, 3.45773e+011), 
        (3.45812e+011, 3.45781e+011), 
        (3.45813e+011, 3.45759e+011), 
        (3.45794e+011, 3.45776e+011)]
ellipses=[]
for i in sources:
    readpath = fitsdir + i #Input directory where FITS files can be found
    bmaj = getval(readpath, 'BMAJ', ext=0)
    bmin = getval(readpath, 'BMIN', ext=0)
    bpa = getval(readpath, 'BPA', ext=0)
    ellipses.append((bmaj, bmin, bpa))
    
if (not(len(sources)==n_plots)) or (not(len(params)==n_plots)):
    raise Exception("Incorrect number of sources or parameters")

small = False

names=[]
for i in sources:
    names.append(i.split("_")[0].upper())

fitsnames = []
for i in names:
    fitsnames.append(i + ".fits")

#Assemble datacubes! 
cubes = []
moments = []
for i in range(n_plots):
    # print(sources[i])
    readpath = fitsdir + sources[i]
    cube = (SpectralCube.read(readpath).with_spectral_unit(u.Hz))
    (bound1, bound2) = params[i]
    cube = cube.spectral_slab(bound1 * u.Hz, bound2 * u.Hz)
    nii_cube = cube.with_spectral_unit(u.km/u.s,
                                       velocity_convention='radio',
                                       rest_value=345.796* u.GHz)  
    cubes.append(nii_cube)
    
    moment_0 = nii_cube.moment(order=0)
    # moment_0 = (cube.moment(order=0))
    moments.append(moment_0)
    savepath = cubedir + fitsnames[i]
    moment_0.write(savepath, overwrite = True)
    moment_0.hdu

#now find the label list!
#'.json/labels_orig.json' can be created with the label_list function in labels.py
label_list()

data = json.loads(open('./json/labels_orig.json').read())
l_list = data["Data"]
no_data = False

#plotting the moment-0 maps

fig, ax = plt.pyplot.subplots(2,2)
fig.set_size_inches(9, 9)
# fig.tight_layout()
fig.subplots_adjust(wspace=.6)
#Hide default axes
for a in ax:
    for b in a:
        b.set_visible(False)

for i in range(n_plots):
    f = aplpy.FITSFigure(moments[i].hdu, figure=fig, subplot=(2, 2, (i+1)))
    name = names[i]
    f.set_title(name)

    (bmaj, bmin, bpa) = ellipses[i]

    #Extract labels from l_list
    #May have multiple values for dPa and e_Pa (deg)
    indices = []
    for j in range(len(l_list)):
        list_name = l_list[j]['Name']
        ind = list_name.rfind(name)
        if ind == -1:
            continue
        try:
            if (list_name[ind+len(name)] == "-"):
                indices.append(j)
        except IndexError:
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
        f.add_label(0.5, 0.95, label_1, True, 'white', weight='bold', size='small')

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
            if small:
                f.add_label(0.5, (0.9 - 0.03*i), label_i, True, 'white', weight='bold', size='xx-small')
            else:
                f.add_label(0.5, (0.9 - 0.05*i), label_i, True, 'white', weight='bold', size='small')

    # f.show_ellipses([300, 300], [10, 10], edgecolor='white', facecolor='white')
    # print((f.pixel2world(300, 300))[0])
    c = f.pixel2world(723,300)
    f.show_ellipses(c[0], c[1], bmaj, bmin, bpa, edgecolor='white', facecolor='white')
    
    # print(f.ax.get_xlim())
    # print(f.ax.get_ylim())
    f.ax.set_xlim(255.5, 767.5)
    f.ax.set_ylim(255.5, 767.5)

    f.show_colorscale(vmin = 0, cmap='plasma')
    f.add_colorbar()
    f.colorbar.set_axis_label_text("$Jy/beam\:km/s$")
    f.colorbar.set_axis_label_pad(1)
savepath = figdir + figname
fig.savefig(savepath)