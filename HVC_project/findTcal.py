#!/usr/bin/env python3

'''
Given a .txt or .csv file containing peak intensities and elevations of 
calibration scans of the source 3C48, as output by fluxcheck.pro, fit, plot,
and compute the actual calibration temperatures for both polarizations
'''

import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

f = 1.420 #frequency in GHz
JyToK = 1.985 #conversion rate from Jansky to K
sfluxCutoff = 29 #maximum value for sflux. Anything exceeding this will be
    #flagged and removed
tau = 0.01 #zenith opacity, default 0.01

#Does the csv file containing all the peak fluxes and elevations exist? If not,
#set csvExists to False in order to initalize it
#txtFilepath should be the txt file containing this data
#this can be obtained by running fluxcheck.pro
csvExists = True
txtFilepath = "./results/sflux.txt"
csvFilepath = "./results/sflux.csv"
figurePath = "./figures/" #relative path in which to save figures

#Define some colors
c1="#C17BD4"
c2="#8ED47B"
c3="#947BD4"
c4="#BBD47B"
c5="#7B8ED4"

def rmFlag(flagInd, *args):
    '''
    Description: Remove flagged indices from a list or lists

    Inputs:
        flagInd (list of ints): flagged indices to be removed from the dataset

    *Args:
        arg (list): List for which the user wants to remove flags
    
    Outputs:
        newLists (list of lists): Lists with flagged indices removed. This list
        contains as many lists as were passed into *args
    '''
    flagInd.sort(reverse=True)
    newLists = []
    for j in args:
        for i in flagInd:
            del j[i]
        newLists.append(j)
    return newLists

def getInd(ind, *args):
    '''
    Description: Obtain specified indices from a list or lists

    Inputs:
        ind (list of ints): indices to be obtained from the dataset

    *Args:
        arg (list): List for which the user wants to obtain values at certain
        indices
    
    Outputs:
        newLists (list of lists): Lists that consist only of the values at the
        specified indices given by ind. This list contains as many lists as
        were passed into *args
    '''
    newLists = []
    for j in args:
        newList = []
        for i in ind:
            newList.append(j[i])
        newLists.append(newList)
    return newLists

def f1(el,a,b,tau):
    # Test function 1
    return a + b*tau/(np.sin(el*np.pi/180.0))

def f2(el,a,b,tau):
    # Test function 2
    return a*np.exp(b*tau/(np.sin(el*np.pi/180.0)))

def f3(el,a,tau):
    '''
    Description: Equation for flux intensity as a function of elevation

    Inputs:
        el (list of floats): elevation (degrees)
        a (float): maximum intensity (K)
        tau (float): zenith opacity (unitless)
    
    Outputs:
        sflux (list of floats): flux intensity as a function of elevation (K)
    '''
    return a*np.exp(-1*tau/(np.sin(el*np.pi/180.0)))

#coefficients found from Perley & Butler for 3C48
#so essentially this code only works for source 3C48
a0 = 1.3253
a1 = -0.7553
a2 = -0.1914
a3 = 0.0498
T3C48 = JyToK * np.power(10,(a0 + a1*np.log10(f) + a2*np.power(np.log10(f),2) + a3*np.power(np.log10(f),3)))

# No need to run if sflux.csv already exists
# Note that this only works for a specific text format and is kind-of single purpose: improve if time allows

csvExists = True

if not csvExists:
    sftxt = open(txtFilepath, "r+")
    sflines = (sftxt.readlines())

    #sflux.txt needs both lines consolidated to one row

    c=" "
    sfrows = []

    #what if the file is formatted differently? This seems too hardcoded
    #but the .txt file output from sflux.pro is kind of weird so it works
    for i in range(1,len(sflines),2):
        line1 = sflines[i].strip()
        line1 = c.join(line1.split())
        line2 = sflines[i+1].strip()
        line2 = c.join(line2.split())
        sflist = line1.split(c)
        sflist.extend(line2.split(c))
        for j in range(1,3):
            sflist[j] = int(sflist[j])
        for j in range(3,10):
            sflist[j] = float(sflist[j])
        sfrows.append(sflist)

    with open(csvFilepath, "w+") as sfcsv:
        csv_writer = csv.writer(sfcsv, delimiter=",")
        csv_writer.writerow(sflines[0].strip().split())
        for i in sfrows:
            csv_writer.writerow(i)

#Populate lists with previously calculated values
sfrows = []
sflux = []
stdSflux = []
meanTcal = []
az = []
el = []
ra = []
dec = []

with open(csvFilepath, "r") as sfcsv:
    csv_reader = csv.reader(sfcsv)
    for row in csv_reader:
        try:
            sflux.append(float(row[3]))
            stdSflux.append(float(row[4]))
            meanTcal.append(float(row[5]))
            az.append(float(row[6]))
            el.append(float(row[7]))
            ra.append(float(row[8]))
            dec.append(float(row[9]))
            sfrows.append(row)
        except ValueError:
            pass

#Flag and exclude outliers
flag = []
for i in range(len(sfrows)):
    if sflux[i] > sfluxCutoff:
        flag.append(i)

sflux, stdSflux, meanTcal, az, el, ra, dec = rmFlag(flag, sflux, stdSflux, meanTcal, az, el, ra, dec)

#Separate by polarization: even and odd numbers
#This only works if the polarizations are actually listed as every other entry
pol0 = np.arange(0,len(sflux),2)
pol1 = np.arange(1,len(sflux),2)

sflux0, stdSflux0, meanTcal0, az0, el0, ra0, dec0 = getInd(pol0, sflux, stdSflux, meanTcal, az, el, ra, dec)
tcal0 = np.mean(meanTcal0)

sflux1, stdSflux1, meanTcal1, az1, el1, ra1, dec1 = getInd(pol1, sflux, stdSflux, meanTcal, az, el, ra, dec)
tcal1 = np.mean(meanTcal1)

# fig = plt.figure(figsize=(13,5))
fig = plt.figure()
# ax1 = plt.subplot(1,2,1)
ax2 = plt.subplot(1,1,1)

# ax1.errorbar(az0,sflux0,yerr=stdSflux0,c=c1,fmt="None",elinewidth=1)
# ax1.scatter(az0,sflux0,c=c1,s=5,marker="v")
# ax1.errorbar(az1,sflux1,yerr=stdSflux1,c=c2,fmt="None",elinewidth=1)
# ax1.scatter(az1,sflux1,c=c2,s=5,marker="^")
# ax1.set_title("Source Flux vs Azimuth")
# ax1.set_xlabel("Azimuth ($^{\circ}$)")
# ax1.set_ylabel("Source Flux (K)")
# ax1.legend(["YY", "XX"],loc=2)

ax2.errorbar(el0,sflux0,yerr=stdSflux0,c=c1,fmt="None",elinewidth=1)
ax2.scatter(el0,sflux0,c=c1,s=5,marker="v")
ax2.errorbar(el1,sflux1,yerr=stdSflux1,c=c2,fmt="None",elinewidth=1)
ax2.scatter(el1,sflux1,c=c2,s=5,marker="^")
ax2.set_title("Source Flux vs Elevation")
ax2.set_xlabel("Elevation ($^{\circ}$)")
ax2.set_ylabel("Source Flux (K)")

xrange = np.linspace(np.min(el), np.max(el))
e = 1e-5

# popt,pcov = curve_fit(f2,el,sflux,bounds=((0,-np.inf,0.01-e),(np.inf,0,0.01+e)))
# ax2.plot(xrange,f2(xrange,popt[0],popt[1],popt[2]),c=c5)
# a1 = popt[0]
# b1 = popt[1]
# l1 = "$" + str(round(a1,2)) + "\cdot e^{(" + str(round(b1,2)) + "*\\tau/\sin(el))}$"

popt,pcov = curve_fit(f3,el0,sflux0,bounds=((0,0.01-e),(np.inf,0.01+e)))
ax2.plot(xrange,f3(xrange,popt[0],popt[1]),c=c3)
a2 = popt[0]
l2 = "$" + str(round(a2,3)) + "\cdot e^{(-\\tau/\sin(el))}$"

popt,pcov = curve_fit(f3,el1,sflux1,bounds=((0,0.01-e),(np.inf,0.01+e)))
ax2.plot(xrange,f3(xrange,popt[0],popt[1]),c=c4)
a3 = popt[0]
l3 = "$" + str(round(a3,3)) + "\cdot e^{(-\\tau/\sin(el))}$"

ax2.legend(["YY", "XX", l2, l3], loc=2)

plt.savefig(figurePath+"sfluxVsEl.png",dpi=300)

# factor * t / sin(el)
# or was it factor * exp(t/sin(el)) ?

# 1.98 5 K per jansky

tcalpol0 = tcal0 * T3C48 / a2
tcalpol1 = tcal1 * T3C48 / a3

print("Calculated T of 3C48: " + str(round(T3C48,3)) + "K")
print("AYY = " + str(round(a2,3)))
print("AXX = " + str(round(a3,3)))
print("Original Tcal for YY = " + str(round(tcal0,3)) + "K")
print("Original Tcal for XX = " + str(round(tcal1,3)) + "K")
print("Corrected Tcal for YY = " + str(round(tcalpol0,5)) + "K")
print("Corrected Tcal for XX = " + str(round(tcalpol1,5)) + "K")