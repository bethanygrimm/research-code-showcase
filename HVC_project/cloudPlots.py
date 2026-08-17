#!/usr/bin/env python3

'''
There are a number of plots on this .py file, and thus might be more conducive
for a .ipynb notebook. This script reads the HVC parameters from CSV files, as
defined in cloudParameters.py. It also derives new quantities and puts them in
a new CSV file.
Plots this script makes include
    Spatial scatter plots of the HVCs, either labeled or colorbar-ed with a
    quantity of interest
    Angle from M31 vs. velocity plots
    Dynamic mass vs. measured HI mass
'''

import csv
import numpy as np
import matplotlib.pyplot as plt
import astropy
from astropy import units as u
import astropy.coordinates as coord
import glob
import os
import sys
from velConversion import velConversion
from angleFromM31 import angleFromM31

# defining a color scheme
# thank you to the colorblind-friendly color scheme checker! https://davidmathlogic.com/colorblind/#
c1 = "#000000"
c2 = "#F54E4E"
c3 = "#5656FB"
c4 = "#8ED47B"
c5 = "#947BD4"
c6 = "#F6B911"

# files to use to find HVC parameters and Magellanic Stream HVC parameters
f1 = "HVCData/M31HVC+bridge_table.txt"
f2 = "HVCData/M31M33_deep_cloudlist.txt"
f3 = "HVCData/M33_HVC_table.txt"
f4 = "HVCData/asu.tsv" #table of HVCs
useData = [f1, f3]
deHeij = [237, 304, 307, 391, 402]

path = "./" #relative path
filepath = path + "data/cloudParameters.csv" #path with all the
    # cloud paramters CSV files. this can be obtained from cloudParameters.py
outputpath = path + "data/derivedCloudParameters.csv"
    # output path

def avg(reader, key, errorKey):
    '''
    Description: return the average and average error from list of dicts reader
    where key is the key to be averaged, and errorKey is the error of key
    '''
    sum = 0
    errorSum = 0
    points = 0
    for row in reader:
        # print(row)
        points += 1
        sum += float(row[key])
        errorSum += ((float(row[errorKey]))**2)
    if points == 0:
        return 0, 0
    else:
        sum = sum/points
        errorSum = np.sqrt(errorSum)/points
    return sum, errorSum

def weightedAvg(reader, key, errorKey, weight, errorWeight):
    '''
    Description: return the weighted average and weighted average error from
    list of dicts reader, where key is the key to be averaged, weight is the
    key with the weighting, errorKey is the error of key, and errorWeight is
    the error of the weighting
    '''
    numSum = 0
    denomSum = 0
    errorSum = 0
    for row in reader:
        numSum += float(row[key]) * float(row[weight])
        denomSum += float(row[weight])
    if numSum == 0:
        return 0,0
    else:
        weightedSum = numSum / denomSum
        for row in reader:
            errorSum +=(float(row[errorKey]) * float(row[weight]) / denomSum)**2
            errorSum += (float(row[errorWeight]) * ((denomSum * float(row[key])) - numSum) / denomSum)**2
        errorSum = np.sqrt(errorSum)
    return weightedSum, errorSum

def icrsToGal(long, lat):
    '''
    Description: convert longitude (RA) and latitude (Dec) from ICRS frame to
    the galactic frame; return galactic longitude (l) and latitude (b)
    '''
    c_icrs = coord.SkyCoord(ra=long*u.degree, dec=lat*u.degree, frame='icrs')
    c_gal= c_icrs.galactic
    l = (c_gal.l.deg)
    b = (c_gal.b.deg)
    return l, b    

#read from text files
#this segment collects galactic longitude, latitude, angle from M31, LSR
#velocity, and LGSR velocity from the HVCs and galaxies in the files specified

cloudL = []
cloudB = []
cloudAngleM31 = []
cloudVlsr = []
cloudVlgsr = []
cloudName = []
cloudLists = [cloudL, cloudB, cloudAngleM31, cloudVlsr, cloudVlgsr]
galaxyL = []
galaxyB = []
galaxyAngleM31 = []
galaxyVlsr = []
galaxyVlgsr = []
galaxyName = []
galaxyLists = [galaxyL, galaxyB, galaxyAngleM31, galaxyVlsr, galaxyVlgsr]
msRa = []
msDec = []
msVlsr = []
msVlgsr = []
msAngleM31 = []

for file in useData:
    with open(file,"r") as f:
        lines = f.readlines()
    for i in range(2,len(lines)-1):
        # print(lines[i].split())
        if not ("Map" in lines[i] or "Deep" in lines[i] or "MS" in lines[i]):
            for j in range(len(cloudLists)):
                cloudLists[j].append(float(lines[i].split()[j]))
            cloudName.append(lines[i].split()[5])
        if "M31" in lines[i] or "M33" in lines[i]:
            for j in range(len(galaxyLists)):
                galaxyLists[j].append(float(lines[i].split()[j]))
            galaxyName.append(lines[i].split()[5][1:4])

with open(f4,"r") as f:
    lines = f.readlines()
    for i in lines:
        # print(i.split())
        try:
            if int(i.split()[0]) in deHeij:
                ra = float(i.split()[3])*15 + float(i.split()[4])/4
                dec1 = float(i.split()[5])
                dec2 = float(i.split()[6])/60
                if dec1 < 0:
                    dec2 = dec2 * -1
                dec = dec1 + dec2
                msRa.append(ra)
                msDec.append(dec)
                msVlsr.append(float(i.split()[7]))
                msVlgsr.append(float(i.split()[8]))
                l, b = icrsToGal(ra, dec)
                msAngleM31.append(angleFromM31(l, b))
        except (ValueError, IndexError):
            pass

#read from the whole cloud parameters
#derive (spatial) areas (there should be 3) -> area in cm^2
#derive uncertainty in area (just 1)
#derive (Gaussian) area (specify whether it should be whole cloud, average, or intensity-weighted average)
#derive NHI
#derive cloud mass
#derive VGLSR velocities
#derive angle from M31

# weighting = 1
#0 for whole-cloud fit properties, 1 for arithmetic average, 2 for intensity-weighted average
areaToUse = 0
#0 for Gaussian, 1 for approximation 1
#default to 0, being the most reliable approximation

NHICONST = 1.82e+18 #1/cm^2
SOLARMASS = 2e+30 #kg
DIST = 800 #kpc
KPCTOCM = 3e+21 #centimeters
PIXTOARCMIN = 10.0/3.0 #arcmin, 1 degree = 60 arcmin
HMASS = 0.001008 #kg/mol
AVOGADRO = 6.022e+23

csvData = [["Name", "l (deg)", "b (deg)", "Area (approx.) (arcmin^2)", "Area (Gaussian) (arcmin^2)", "Area e (arcmin^2)", "Amp (K)", "Amp e (K)", "Center (VLSR) (km/s)", "Center (VLGSR) (km/s)", "Center e (km/s)", "FWHM (km/s)", "FWHM e (km/s)", "Integral (K*km/s)", "Integral e (K*km/s)", "NHI (atoms/cm^2)", "NHI e (atoms/cm^2)", "MHI (Msun)", "MHI e (Msun)", "M Dyn (Msun)", "M Dyn e (Msun)", "Angle (degrees)"]]
#so it's all 3 areas plus error
#all gaussian properties plus VGLSR
#NHI, MHI, and angle

regionName = []
regionAmp = []
regionVlsr = []
regionVlgsr = []
regionFWHM = []
regionAngleM31 = []
regionLats = []
regionLongs = []
regionL = []
regionB = []
regionNhi = []
regionMhi = []
regionMdyn = []
regionMratio = []
with open(filepath, encoding="utf-16") as csvfile:
    csvreader = csv.DictReader(csvfile, delimiter = ',')
    for row in csvreader:
        # print(row)
        # tempweighting = 0
        name = row["Name"]
        regionName.append(name)
        area1 = int(row["Area (approx.)"]) *  PIXTOARCMIN**2
        # area2 = int(row["Area 2 (approx.)"]) *  PIXTOARCMIN**2
        # scaleFactor = float(row["Scale (pix)"])
        boundingArea = float(row["Bounding area (px^2)"])
        major = float(row["Major (arcmin)"])
        minor = float(row["Minor (arcmin)"])
        majorE = float(row["Minor e (arcmin)"])
        minorE = float(row["Major e (arcmin)"])
        area3 = major * minor * np.pi * 0.25 #arcmin
        areaE = np.pi * 0.25 * np.sqrt((major**2)*(majorE) + (minor**2)*(minorE))
        scaleFactor = boundingArea / (area3/(PIXTOARCMIN**2))
        rad = np.sqrt(major * minor) * 60.0 / 1000.0 #kpc
        radE = np.sqrt((majorE**2 * 0.25 * (minor/major)) + (minorE**2 * 0.25 * (major/minor))) * 60.0 / 1000.0 
        long = float(row["Long"])
        regionLongs.append(long)
        lat = float(row["Lat"])
        regionLats.append(lat)
        l, b = icrsToGal(long, lat)
        regionL.append(l)
        regionB.append(b)
        # print(l,b)
        # print(name, long, lat, major, minor)
        try:
            centerVLSR = float(row["Center (km/s)"])
            centerVLGSR = velConversion(l,b,centerVLSR)[2]
        except ValueError:
            centerVLSR, centerVLGSR = np.NaN, np.NaN
        try:
            # print(name, scaleFactor)
            amp = float(row["Amp (K)"]) * scaleFactor
            ampE = float(row["Amp e (K)"]) * scaleFactor
            centerE = float(row["Center e (km/s)"])
            FWHM = float(row["FWHM (km/s)"])
            FWHME = float(row["FWHM e (km/s)"])
            integral = float(row["integral (K*km/s)"]) * scaleFactor
            integralE = float(row["integral e (K*km/s)"]) * scaleFactor
        except ValueError:
            amp, ampE, centerE, FWHM, FWHME, integral, integralE = np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN
            # tempWeighting = 1
        regionAmp.append(amp)
        regionFWHM.append(FWHM)
        nhi = NHICONST * integral
        nhiE = NHICONST * integralE
        regionNhi.append(nhi)
        if areaToUse == 0:
            area = area3
        else:
            area = area1
        # else:
        #     area = area2
        areaCM = area * (np.pi**2) * DIST**2 * KPCTOCM**2 / ((60*180)**2)
        natoms = nhi * areaCM
        natomsE = nhiE * areaCM
        mhi = natoms * HMASS / (AVOGADRO * SOLARMASS)
        mhiE = natomsE * HMASS / (AVOGADRO * SOLARMASS)
        regionMhi.append(mhi)
        mDyn = (2.5e+5) * rad * FWHM**2
        mDynE = (2.5e+5) * np.sqrt((radE*(FWHM**2))**2 + (rad*2*FWHM*FWHME)**2)
        regionMdyn.append(mDyn)
        regionMratio.append(mDyn / mhi)
        # print(long)
        # print(lat)
        angle = angleFromM31(l,b)
        csvData.append([name, l, b, area1, area3, areaE, amp, ampE, centerVLSR, centerVLGSR, centerE, FWHM, FWHME, integral, integralE, nhi, nhiE, mhi, mhiE, mDyn, mDynE, angle])
        if centerVLSR != 0:
            regionVlsr.append(centerVLSR)
            regionVlgsr.append(centerVLGSR)
        else:
            regionVlsr.append(np.NaN)
            regionVlgsr.append(np.NaN)
        regionAngleM31.append(angle)

        # print(centerVLSR, centerVLGSR, nhi, areaCM, natoms, mhi, angle, mDyn, mDynE)

with open(outputpath, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(csvData)

print(outputpath)

#read from CSV. we need all files from the cloudParameters
#specifying which regions are M31 HVCs, which regions are ambiguous, and which
#regions are new

regionsM31 = ["A", "H", "J"]
regionsAmbiguous = ["B", "C", "D", "G", "I", "K"]
regionsNew = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O"]
rM31Longs = []
rM31Lats = []
rAmbLongs = []
rAmbLats = []
rNewAngleM31 = []
rNewVlsr = []
rNewVlgsr = []
for i in range(len(regionName)):
    if regionName[i][-1:] in regionsM31:
        rM31Longs.append(regionLongs[i])
        rM31Lats.append(regionLats[i])
    if regionName[i][-1:] in regionsAmbiguous:
        rAmbLongs.append(regionLongs[i])
        rAmbLats.append(regionLats[i])
    if regionName[i][-1:] in regionsNew:
        rNewAngleM31.append(regionAngleM31[i])
        rNewVlsr.append(regionVlsr[i])
        rNewVlgsr.append(regionVlgsr[i])
# print(rM31Lats)
# print(rM31Longs)

galaxyLong = []
galaxyLat = []

#to-do: improve
#this map is kind of hardcoded but there's no way around it
#This segment involves creating spatial scatter plots, with values of interest
#plotted as the colorbar

mapLongs = [20.833,12,12,10.625,10.625,14.75,14.75,14.375,14.375,24.5,24.5,22.5,22.5,21.75,21.75,20.833]
mapLats = [42.15,42.15,40.15,40.15,37.825,37.825,36.159,36.159,34.333,34.333,36.15,36.15,38.15,38.15,39.883,39.833]
cb = "MHI" #label/name of plot
cbUnit = cb + " (solar masses)" #label/name of plot with unit
colorMap = regionMhi #which dataset to use for the colormap

for i in range(len(galaxyL)):
    gc = coord.SkyCoord(l=galaxyL[i]*u.degree, b=galaxyB[i]*u.degree, frame='galactic')
    gc_icrs = gc.transform_to("icrs")
    long = (gc_icrs.ra.deg)
    lat = (gc_icrs.dec.deg)
    galaxyLong.append(long)
    galaxyLat.append(lat)

# didn't end up needing these but they may come in handy
scaledRegionNhi = np.array(regionNhi) / 1e+18
scaledRegionMhi = np.array(regionMhi) / 1e+5

fig, ax = plt.subplots(1,1)
plt.axis("equal")
plt.fill(mapLongs, mapLats, alpha=0.2, label="mapped region")
if cb == "" or cb == "updated":
    ax.scatter(regionLongs,regionLats,marker="+",c=c1,s=15,label="HVCs between M31 and M33")
    ax.scatter(rM31Longs,rM31Lats,marker="o",c=c6,s=20,label="HVCs associated with M31")
    ax.scatter(rAmbLongs,rAmbLats,marker="D",c=c5,s=15,label="Ambiguous regions")
    cbTitle = ", " + cb
else:
    ax.scatter(regionLongs,regionLats,marker="o",c="white",s=10)
    sc = ax.scatter(regionLongs,regionLats,marker="o",c=colorMap,s=15,cmap="viridis")
    clb = plt.colorbar(sc,label=cbUnit)
    cbTitle = ", " + cbUnit

ax.scatter(galaxyLong,galaxyLat,marker="s",c=c2,s=30)
for i in range(len(galaxyName)):
    plt.annotate(galaxyName[i],(galaxyLong[i]+0.25,galaxyLat[i]+0.25),size="small")
for i in range(len(regionName)):
    # plt.annotate(round(regionFWHM[i],1),(regionLongs[i],regionLats[i]+0.1),size="x-small")
    plt.annotate(regionName[i][-1:],(regionLongs[i]-0.1,regionLats[i]+0.1),size="x-small")
ax.annotate(galaxyLong,galaxyLat)
ax.xaxis.set_inverted(True)

plt.xlabel("Longitude (degrees)")
plt.ylabel("Latitude (degrees)")
plt.title("Known clouds and possible regions" + cbTitle)
plt.legend(loc="lower right")

plt.savefig("./figures/cloudsRegions"+cb+".png",dpi=300 )

#this segment plots VLGSR or VLSR as a function of distance from M31
#if annotations = True is set, the Magellanic Stream HVCs will not be plotted for room
#if annotations = False is set, the Magellanic Stream HVCs will be plotted, and a legend
#will be shown instead
#if plotVlgsr = True is set, LGSR velocity will be plotted
#if plotVlgsr = False is set, LSR velocity will be plotted

annotations = True
# depending on how many things we want to plot, annotating every point is going to get cluttered. a legend may be preferable
plotVlgsr = False

cloudV, regionV, galaxyV, msV = [], [], [], []
velType = ""
fileStr = ""
if plotVlgsr:
    cloudV = cloudVlgsr
    regionV = regionVlgsr
    galaxyV = galaxyVlgsr
    newV = rNewVlgsr
    msV = msVlgsr
    velType = "LGSR"
else:
    cloudV = cloudVlsr
    regionV = regionVlsr
    galaxyV = galaxyVlsr
    newV = rNewVlsr
    msV = msVlsr
    velType = "LSR"

f = plt.figure()
# f.set_size_inches(10,5)

#ignore
sample = [0,5,15]
sample2 = [100,-100,-200]
v1 = 4.7
v2 = 10.2
min = -75
max = 5
vx1 = np.linspace(v1,v1,50)
vy1 = np.linspace(min,max,50)
vx2 = np.linspace(v2,v2,50)
vy2 = np.linspace(min,max,50)
plt.scatter(cloudAngleM31,cloudV,c=c2,marker="o",s=15,label="HVCs around M31 and M33")
plt.scatter(regionAngleM31,regionV,c=c1,marker="+",s=15,label="HI cloud bridge")
plt.scatter(galaxyAngleM31,galaxyV,c=c3,marker="s",s=30,label="Galaxies")
plt.scatter(rNewAngleM31,newV,c=c4,marker="+",s=20)
if annotations:
    for i in range(len(regionName)):
        plt.annotate(regionName[i][-1:],(regionAngleM31[i]+0.2,regionV[i]-1),size="x-small")
    for i in range(len(galaxyName)):
        plt.annotate(galaxyName[i],(galaxyAngleM31[i],galaxyV[i]+5),size="small")
else:
    plt.scatter(msAngleM31,msV,c=c6,marker="p",s=15,label="Magellanic Stream clouds")
    plt.legend()
    fileStr = "WithMS"
# plt.plot(vx1,vy1,c="k")
# plt.plot(vx2,vy2,c="k")

plt.xlabel("Angle from M31 (degrees)")
plt.ylabel(velType + " Velocity (km/s)")
plt.title("Angle from M31 vs. " + velType + " Velocity of Various Clouds")

plt.savefig("./figures/angleM31V"+velType+fileStr+".png",dpi=300)

# This plot plots dynamic mass vs. measured HI mass on a log-log scale

maxMdyn = np.nanmax(regionMdyn)
linear = np.linspace(0,maxMdyn,100)
plt.scatter(regionMdyn, regionMhi, marker="+", c=c1, label="Clouds")
for i in range(len(regionName)):
    plt.annotate(regionName[i][-1:],(regionMdyn[i]*1.05,regionMhi[i]*1.1),size="x-small")
plt.loglog(linear,linear,"k--",label="1:1 mass ratio")

plt.xlabel("Dynamic mass (solar masses)")
plt.ylabel("Measured HI mass (solar masses)")
plt.legend()
plt.title("Dynamic mass vs. measured HI mass")
plt.savefig("./figures/Mratio.png",dpi=300)
