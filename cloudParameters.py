#!/usr/bin/env python3

'''
Read all Gaussian fit CARTA log files from a given directory, and convert them
into a machine-readable CSV file
This is still a little bit hardcoded as it requires the CARTA log fitting
format, and for files to be in certain directories
'''

import os
import glob
import csv
import math

dirpath = "./results/gaussianFits/wholeRegions/"
source = "regionH" #source to analyze: this assumes that there is a directory 
    # with this name
dirpath = "./results/gaussianFits/" + source #this is where all CARTA Gaussian
    # fit logs should be
outputpath = "./results/cloudParameters/" + source + ".csv"

#ideally fit should already be done in km/s (or MHz), but if channel information is known, change these values
#before running, ensure that all log files for a given source are in the same directory, and that directory paths are correct

# nChans = 621
# vel0 = 400.263
# velN = -399.7536
nChans = 1244
vel0 = 400.585
velN = -400.7198

#simple auxiliary functions to convert channels to km/s given a range of channels

def chanToKmS(value,nChans,vel0,velN):
    '''
    Description: convert channel to velocity. Probably doesn't actually have to
        be in km/s, despite the name, but it will help with units later on.
        This function takes into account the channel 0 velocity and is thus
        good for computing center velocities or quantities where the channel
        position actually matters

    Inputs:
        value (int): channel to be converted to velocity
        nChans (int): number of total channels in the datacube
        vel0 (float): velocity corresponding to channel 0
        velN (float): veloctiy corresponding to the maximum channel
    
    Outputs:
        vel (float): velocity corresponding with given channel
    '''
    return vel0 + value * ((velN - vel0)/nChans)

def chanToKmSE(value,nChans,vel0,velN):
    '''
    Description: convert channel to velocity. Probably doesn't actually have to
        be in km/s, despite the name, but it will help with units later on.
        This function does not take into account the channel 0 velocity and is 
        thus good for computing errors, FWHMs, or quantities where the channel
        position doesn't actually matter

    Inputs:
        value (int): channel to be converted to velocity
        nChans (int): number of total channels in the datacube
        vel0 (float): velocity corresponding to channel 0
        velN (float): veloctiy corresponding to the maximum channel
    
    Outputs:
        vel (float): velocity corresponding with given channel
    '''
    return value * abs((velN - vel0)/nChans)

# To-do: improve so this is less hard-coded
# This for-loop extracts and calculates values of interest for every file labeled "*Stats*" in a given directory
# These should be the .tsv files downloaded from CARTA from the Statistics regions
# longLatFloat = True returns the center longitude and latitude of the region as floats
# longLatDeg = True returns the center longitude and latitude of the region in H:M:S/D:M:S format
# scaleFactor = True is kind of useless but essentially returns the number of pixels in the region
# but this is kind of already in the Stats file
# recommended to only set one boolean to True at a time so they don't all print at once
# Note that this returns values from an overall region, and not a Gaussian fit

longLatFloat = False
longLatDeg = False
scaleFactor = False

for filename in glob.glob(os.path.join(dirpath, '*Stats*')):
    #file by file
    with open(os.path.join(os.getcwd(), filename), 'r') as f:
        #line by line
        lines = (f.readlines())
        sum, mean = 0, 0
        for l in lines:
            if longLatFloat:
                if "wcs:FK5" in l:
                    # print(l.split())
                    # print(l.split()[1][22:-1])
                    # print(filename)
                    long = (int(l.split()[1][17:18])*15) + (float(l.split()[1][19:21])/4) + (float(l.split()[1][22:-1])/240)
                    lat = (int(l.split()[2][0:2])) + (float(l.split()[2][3:5])/60) + (float(l.split()[2][6:-2])/3600)
                    print(long,"\t",lat)
            if longLatDeg:
                if "wcs:FK5" in l:
                    long = (l.split()[1][17:26])
                    lat = (l.split()[2][0:10])
                    print(long,"\t",lat)
            if scaleFactor:
                if "Sum" in l:
                    sum = float(l.split()[1])
                elif "Mean" in l:
                    mean = float(l.split()[1])
                    print(sum)
                    print(mean)
                    scaleFactor = sum/mean
                    print(scaleFactor)

#make sure we don't duplicate pixels - floor pixels and check
#convert coordinates to actual numbers

pixels = []
unit = ""
csvData = [["Pixel", "Long (degrees)", "Lat (degrees)", "amp (K)", "amp e (K)", "center (km/s)", "center e (km/s)", "FWHM (km/s)", "FWHM e (km/s)", "integral (K*km/s)", "integral e (K*km/s)"]]

#extract information from log files
#to-do (if time allows): this could be less hard-coded
for filename in glob.glob(os.path.join(dirpath, '*.txt')):
    #file by file
    with open(os.path.join(os.getcwd(), filename), 'r') as f:
        #line by line
        lines = (f.readlines())
        skip = False
        isPoint = False
        for l in lines:
            #don't record duplicate pixels
            if "pixel" in l:
                isPoint = True
                pix1 = math.floor(float(l.split()[3][1:10]))
                pix2 = math.floor(float(l.split()[4][0:9]))
                pCurrent = (pix1,pix2)
                # print(pCurrent)
                for p in pixels:
                    if pCurrent == p:
                        skip = True
                if not skip:
                    pixels.append((pix1,pix2))
            if "wcs:FK5" in l:
                long = (int(l.split()[3][1:2])*15) + (float(l.split()[3][3:5])/4) + (float(l.split()[3][6:-2])/240)
                lat = (int(l.split()[4][0:2])) + (float(l.split()[3][3:5])/60) + (float(l.split()[3][6:-2])/3600)
            if "amp1" in l:
                amp = float(l.split()[2])
                ampE = float(l.split()[5])
            if "center1" in l:
                center = float(l.split()[2])
                centerE = float(l.split()[5])
                if "km/s" in l.split()[3]:
                    unit = "km/s"
                elif "Channel" in l.split()[3]:
                    unit = "Channel"
                    center = chanToKmS(center, nChans, vel0, velN)
                    centerE = chanToKmSE(centerE, nChans, vel0, velN)
            if "fwhm1" in l:
                fwhm = float(l.split()[2])
                fwhmE = float(l.split()[5])
                if unit == "Channel":
                    fwhm = chanToKmSE(fwhm, nChans, vel0, velN)
                    fwhmE = chanToKmSE(fwhmE, nChans, vel0, velN)
            if "integral" in l:
                integral = float(l.split()[4])
                integralE = float(l.split()[9])
                if unit == "Channel":
                    integral = chanToKmSE(integral, nChans, vel0, velN)
                    integralE = chanToKmSE(integralE, nChans, vel0, velN)
        if (not skip) and isPoint:
            # print("appending")
            csvData.append([pCurrent,long,lat,amp,ampE,center,centerE,fwhm,fwhmE,integral,integralE])

#write to CSV
with open(outputpath, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(csvData)

print(outputpath)