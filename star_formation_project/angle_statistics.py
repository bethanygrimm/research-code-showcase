#imports
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import json
import math
from scipy import stats
from scipy.optimize import curve_fit
from scipy.stats import rice
import scipy.special as scimath
from numpy.polynomial.polynomial import polyfit
import emcee
import corner
import csv

#function definitions, including MCMC methods

def val_list(key: str, d: list) -> list:
    '''
    Given a key and a list of dicts, this extracts all values with a certain key into a list to return.

    Inputs:
        key (str): the key for which values are to be extracted
        d (list of dicts): the list of dicts for which all values of a certain key are needed

    Outputs:
        l (list): list of all values in the list of dictionaries with the given key
    '''
    l = []
    for i in d:
        try:
            l.append(float(i[key]))
        except ValueError:
            l.append(i[key])
        except KeyError:
            l.append(None)
    return l

def angle_list(k1: str, k2: list, d: list) -> list:
    '''
    Given two keys and a list of dicts, this extracts all values with a certain key into a list to return.
    This is specifically designed to return angle values, which, in this code, come with two keys: the lobe (blue or red), and the type of angle

    Inputs:
        k1 (str): the first key for which values are to be extracted, in this case, "Position angle", "Opening angle", "Opening corrected",
            "Position angle error", "Opening angle error", or "Opening corrected error"
        k2 (str): the second key for which values are to be extracted, in this case "Blue" or "Red"
        d (list of dicts): the list of dicts for which all values of a certain key are needed

    Outputs:
        l (list of floats): list of all values in the list of dictionaries with the given keys
    '''
    l = []
    for i in d:
        for k in k2:
            try:
                angle = float(i[k1][k])
                l.append(angle)
            except KeyError:
                l.append(None)
            except TypeError:
                l.append(None)
    return l

def removeValues(removeList, nameList, *args):
    '''
    This function returns all lists, with all values with the keys corresponding to removeList removed

    Inputs:
        removeList (list): a list of all "keys" to have removed. Since there aren't lists, nameList functions as the list of keys.
            Thus, nameList must have the same length as all *args
        nameList (list): a list of all "keys" corresponding to the lists in *args. nameList must have the same length as all *args
        *args (list(s)): any list(s) to have values removed. Every list must be the same length as nameList

    Outputs:
        returnLists (list of numpy arrays): all lists with removeList values removed, packaged into a list in the same order as *args
    '''
    removeIndices = []
    newNames = []
    for i in range(len(nameList)):
        for j in removeList:
            if nameList[i]==j:
                removeIndices.append(i)
    removeIndices.sort(reverse="True")
    returnLists = []
    # print(len(newNames))
    for i in args:
        iList = np.delete(i,removeIndices)
        returnLists.append(iList)
        # print(len(iList))
    return returnLists

def toNp(*args):
    '''
    This function returns all lists passed as numpy arrays.

    Inputs:
        *args (list(s)): any list(s) to be converted into numpy arrays

    Outputs:
        returnLists (list of numpy arrays): all lists converted into numpy arrays, packaged into a list in the same order as *args
    '''
    returnLists = []
    for i in args:
        newList = np.array(i, dtype=float)
        returnLists.append(newList)
    return returnLists

def isNaN(num):
    '''
    A simple check whether num is NaN or not. Returns a boolean.
    '''
    return num != num

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

def rician(x, b, loc, scale):
    '''
    A Rician function.

    Inputs:
        x (float or numpy array): x of the function
        b (float): corresponding width of the function
        loc (float): corresponding to x-shift of the function
        scale (float): corresponding to amplitude of the function

    Output:
        output (float or numpy array): resultant Rician function
    '''
    return stats.rice.pdf(x, b=b, loc=loc, scale=scale)

def movingAverage(xdata, ydata, window, type, llim=0, rlim=0):
    '''
    This function returns a simple moving average for a dataset.

    Inputs:
        xdata (list of floats): x-coordinates of data for which a moving average is to be found. The moving average "moves" through x.
        ydata (list of floats): y-coordinates of data for which a moving average is to be found. The ydata is what gets averaged.
        window (float): must be > 0. The interval size over which to average
        type (str): must be "lin" or "log". Accounts for either a linear or logarithmic x.
        llim (float): The leftmost value at which to start the moving average. Defaults to 0
        rlim (float): The rightmost value at which to end the moving average. Defaults to 0

    Outputs:
        xrange (list of floats): x-coordinates of moving average
        yrange (list of floats): y-coordinates of moving average
    '''
    density = 1
    if len(xdata) != len(ydata):
        raise Exception("Length mismatch of xdata and ydata")
    if type != "lin" and type != "log":
        raise Exception("Type should be 'lin' or 'log'")
    if window <=0:
        raise Exception("Choose a non-negative interval")
    if llim == 0:
        xmin = np.nanmin(xdata)
    else:
        xmin = llim
    if rlim == 0:
        xmax = np.nanmax(xdata)
    else:
        xmax = rlim
    if type == "lin":
        intervals = int(np.ceil(density*(xmax-xmin)/(window)))
        xrange = np.linspace(xmin,(xmin+intervals*window/density),(intervals+1))
        xints = xrange - (window/density)
        xints = np.append(xints,(xrange[-1]+(window/density)))
    elif type == "log":
        if window <= 2:
            raise Exception("Choose a larger interval")
        intervals = int(np.ceil(np.log(xmax-xmin)/np.log(window)))
        xrange = np.logspace(np.log(xmin)/np.log(window),np.log(xmax)/np.log(window),(intervals+1),base=window)
        xints = [xrange[0]-1]
        for i in range(len(xrange)-1):
            xints.append(np.sqrt(xrange[i]*xrange[i+1]))
        xints.append(xrange[-1]+1)
    yrange = np.zeros(len(xrange))
    for i in range(len(xints)-1):
        movingSum = 0
        points = 0
        for j in range(len(xdata)):
            try:
                xbound = xints[i+1]
            except IndexError:
                xbound = xints[-1]
            if (xdata[j] > xints[i]) and (xdata[j] < xbound):
                if not np.isnan(ydata[j]):
                    movingSum = movingSum + ydata[j]
                    points = points+1
        try:
            yrange[i] = float(movingSum) / float(points)
        except ZeroDivisionError:
            yrange[i] = 0
    return (xrange,yrange)

def listMethod(list1,list2,list3):
    '''
    This is an auxiliary method for list comprehension: given three lists of identical sizes, it returns all of list1 for which
    list2 and list3 have values at that index.

    Inputs:
        list1 (list): list whose values are to be returned
        list2 (list): one list for comparison: if any of these values are NaN, do not return the corresponding value for list1.
            For this reason, list1 and list2 must have the same length
        list3 (list): the other list for comparison: if any of these values are NaN, do not return the corresponding value for list1.
            For this reason, list1 and list2 must have the same length
    
    Outputs:
        newlist (list): All values of list1 for which the corresponding values of list2 and list3 are not NaN
    '''
    if len(list1) != len(list2) or len(list1) != len(list3):
        raise Exception("Inconsistent list sizes")
    newList = [list1[i] for i in range(len(list2)) if not (isNaN(list2[i]) or isNaN(list3[i]))]
    newList = np.asarray(newList)
    return newList

'''
The following functions are all built for MCMC fitting methods.
'''
#tanh function
def f1(theta, x):
    '''
    Returns a*tanh((x-b)/c), where theta = [a,b,c]
    '''
    # a: theta_0
    # b: Tbol_shift
    # c: Tbol_0
    a,b,c = theta
    a,b,c = float(a), float(b), float(c)
    #wait...how would this work in radians?
    f1 = a * np.tanh((x-b)/c)
    return f1

#(1-exp) function
def f2(theta, x):
    '''
    Returns a*(1-exp((b-x)/c)), where theta = [a,b,c]
    '''
    # a: theta_0
    # b: Tbol_shift
    # c: Tbol_0
    a,b,c = theta
    a,b,c = float(a), float(b), float(c)
    f2 = a * (1 - np.exp(-1*(x-b)/c))
    return f2

#power law function
def f3(theta, x):
    '''
    Returns a*(x^b)+c, where theta = [a,b,c]
    '''
    a,b,c = theta
    a,b,c = float(a), float(b), float(c)
    f3 = a * np.power(x, b) + c
    return f3

#split power law function
def f4(theta, x):
    '''
    Returns a broken power law function: a*((x/d)^b) where x<d, and a*((x/d)^c) where x>d, where theta = [a,b,c,d]
    '''
    # d: cutoff in K
    a,b,c,d = theta
    a,b,c,d = float(a), float(b), float(c), float(d)
    f4 = np.where(x < d, a*(scimath.powm1((x/d),b)+1), a*(scimath.powm1((x/d),c)+1))
    return f4

def log_prior_1(theta):
    '''
    The prior for function 1 (tanh) requires a to be between 0 and 180, and for b and c to be greater than 0.
    Again, theta = [a,b,c]
    '''
    a,b,c = theta
    a,b,c = float(a), float(b), float(c)
    if (0 < a < 180) and (b > 0) and (c > 0):
        return 0.0
    else:
        return -np.inf
    
def log_prior_2(theta):
    '''
    The prior for function 2 (1-exp) requires a to be between 0 and 180, and for c to be greater than 0.
    Again, theta = [a,b,c]
    '''
    a,b,c = theta
    a,b,c = float(a), float(b), float(c)
    if (0 < a < 180) and (c > 0):
        return 0.0
    else:
        return -np.inf
    
def log_prior_3(theta):
    '''
    The prior for function 3 (power law) requires a to be greater than 0.
    Again, theta = [a,b,c]
    '''
    a,b,c = theta
    a,b,c = float(a), float(b), float(c)
    if a > 0:
        return 0.0
    else:
        return -np.inf
    
def log_prior_4(theta):
    '''
    The prior for function 4 (broken power law) requires a and d to be greater than 0.
    Again, theta = [a,b,c,d]
    '''
    a,b,c,d = theta
    a,b,c,d = float(a), float(b), float(c), float(d)
    if a > 0 and d > 0:
        return 0.0
    else:
        return -np.inf

def log_likelihood(theta, x, y, sigma, fit):
    '''
    This function provides the likelihood of the fit chosen for dataset (x,y).
    Fit must be 1, 2, 3, or 4, corresponding to the function to use. (tanh: 1; 1-exp: 2; power: 3; broken power: 4)
    Theta corresponds to the theta for each function, and has the input parameters.
    '''
    x = np.asarray(x)
    y = np.asarray(y)
    if fit==1:
        y_model = f1(theta, x)
    elif fit==2:
        y_model = f2(theta, x)
    elif fit==3:
        y_model = f3(theta, x)
    elif fit==4:
        y_model = f4(theta, x)
    else:
        raise Exception("Provide valid fit index (tanh: 1; 1-exp: 2; power: 3; broken power: 4)")
    # return -0.5 * np.sum(((y-y_model)/sigma) ** 2) + np.sum(np.log(2*np.pi*(sigma**2)))
    return -0.5 * np.sum((y-y_model) ** 2)

def log_posterior_1(theta, x, y, sigma):
    '''
    Calculates the posterior (prior * likelihood) for fit function 1 (tanh)
    '''
    x = x.astype(np.float64)
    y = y.astype(np.float64)
    lp = log_prior_1(theta)
    if (math.isinf(lp) and lp < 0):
        return -np.inf
    return lp + log_likelihood(theta, x, y, sigma, 1)

def log_posterior_2(theta, x, y, sigma):
    '''
    Calculates the posterior (prior * likelihood) for fit function 2 (1-exp)
    '''
    x = x.astype(np.float64)
    y = y.astype(np.float64)
    lp = log_prior_2(theta)
    if (math.isinf(lp) and lp < 0):
        return -np.inf
    return lp + log_likelihood(theta, x, y, sigma, 2)

def log_posterior_3(theta, x, y, sigma):
    '''
    Calculates the posterior (prior * likelihood) for fit function 3 (power law)
    '''
    x = x.astype(np.float64)
    y = y.astype(np.float64)
    lp = log_prior_3(theta)
    if (math.isinf(lp) and lp < 0):
        return -np.inf
    return lp + log_likelihood(theta, x, y, sigma, 3)

def log_posterior_4(theta, x, y, sigma):
    '''
    Calculates the posterior (prior * likelihood) for fit function 4 (broken power law)
    '''
    x = x.astype(np.float64)
    # y = y.astype(np.float64)
    lp = log_prior_4(theta)
    if (math.isinf(lp) and lp < 0):
        return -np.inf
    return lp + log_likelihood(theta, x, y, sigma, 4)

def emcee_fit(x_data, y_data, error_data, guesses, nwalkers, nburn, nsteps, xfit, fit):
    '''
    Computes the best fit for the data with the MCMC method.

    Inputs:
        x_data (list of floats or numpy array): x-coordinates of the data to fit
        y_data (list of floats or numpy array): y-coordinates of the data to fit
        error_data (list of floats): error values of y_data
        guesses (list of floats): initial guesses for the parameters, based off which fitting function used.
            (Essentially, first guess for theta)
        nwalkers (int): number of MCMC walkers to use
        nburn (int): number of "burn-in" steps to use
        nsteps (int): number of MCMC steps to take per walker
        xfit (list of floats or numpy array): x-coordinates for which to find the fit. Usually a linspaced array across the x domain.
        fit (int): must be 1, 2, 3, or 4, corresponding to the function to use. (tanh: 1; 1-exp: 2; power: 3; broken power: 4)
    
    Outputs:
        mu (numpy array): mean of average best fit found
        sigma (numpy array): standard deviation of average best fit found
        sampler:
        parameters:
    '''
    error_data = error_data*5
    data = (x_data, y_data, error_data)
    ndim = len(guesses)
    p0 = [np.array(guesses) + 1e-7 * np.random.randn(ndim) for i in range(nwalkers)]
    if fit==1:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior_1, args=data)
    elif fit==2:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior_2, args=data)
    elif fit==3:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior_3, args=data)
    elif fit==4:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior_4, args=data)
    else:
        raise Exception("Invalid fit (try 1 for tanh, 2 for (1-exp), 3 for power, 4 for broken power)")
    p0, _, _ = sampler.run_mcmc(p0, nburn)
    sampler.reset()
    pos, prob, state = sampler.run_mcmc(p0, nsteps)
    emcee_trace = sampler.chain[:, nburn:, :].reshape(-1, ndim).T
    try:
        alpha, beta, gamma, delta = emcee_trace[:4]
    except ValueError:
        alpha, beta, gamma = emcee_trace[:3]
    samples = sampler.flatchain
    parameters = samples[np.argmax(sampler.flatlnprobability)]
    if fit==1:
        yfit = alpha[:,None] * np.tanh((xfit-beta[:,None])/gamma[:,None])
    elif fit==2:
        yfit = alpha[:,None] * (1 - np.exp(-1*(xfit-beta[:,None])/gamma[:,None]))
    elif fit==3:
        yfit = alpha[:,None] * np.power(xfit, beta[:,None]) + gamma[:,None]
    elif fit==4:
        yfit = np.where(xfit < delta[:,None], alpha[:,None] * (scimath.powm1((xfit/delta[:,None]), beta[:,None])+1), alpha[:,None] * (scimath.powm1((xfit/delta[:,None]), gamma[:,None])+1))
    else:
        raise Exception("Invalid fit (try 1 for tanh, 2 for (1-exp), 3 for power, 4 for broken power)")
    mu = yfit.mean(0)
    sig = 2 * yfit.std(0)
    return (mu, sig, samples, parameters)

figdir = "./figures/"
jsondir = "./json/"
jsonpath = jsondir + "labels.json"

#get list of labels ready
data = json.loads(open(jsonpath).read())
l_list = data["Data"]

movingAverages = False #plot moving averages?
excludeOutliers = True #exclude outlier sources?

#Data values to neglect (just names of the sources):
neglect = ["HOPS-85","HOPS-174","HOPS-179", "HOPS-369", "HOPS-385"]

#Some portions of the data analysis have different sources to neglect than others.
#oa_neglect are all the sources for which the position angles are fine, but not the opening angles
oa_neglect = ["HOPS-59-A","HOPS-166","HOPS-287","HOPS-321","HOPS-345","HOPS-359","HOPS-386-B","HOPS-405","HOPS-408"]
#oab_neglect are all the sources where only the opening angle for the blueshifted lobe is to be neglected
oab_neglect = ["HOPS-10", "HOPS-12-B", "HOPS-182-A", "HOPS-247"]
#oar_neglect are all the sources where only the opening angle for the redshifted lobe is to be neglected
oar_neglect = ["HOPS-29", "HOPS-94", "HOPS-95", "HOPS-123", "HOPS-198", "HOPS-203", "HOPS-281-A"]
#single_b are all the sources for which there is only a blueshifted lobe
single_b = ["HOPS-19", "HOPS-45", "HOPS-50", "HOPS-84-B", "HOPS-86-A", "HOPS-90", "HOPS-92-B", "HOPS-139", "HOPS-160", "HOPS-206", "HOPS-229", "HOPS-244", "HOPS-288-A-A", "HOPS-312-A", "HOPS-323-B", "HOPS-365", "HOPS-383", "HOPS-384-B", "HOPS-395", "HOPS-397"]
#single_b are all the sources for which there is only a redshifted lobe
single_r = ["HOPS-1", "HOPS-78-A", "HOPS-84-A", "HOPS-92-A", "HOPS-96", "HOPS-152", "HOPS-288-B", "HOPS-312-B", "HOPS-340", "HOPS-384-A", "HOPS-386-A", "HOPS-387-A", "HOPS-406"]

ofInterest = ["HOPS-11", "HOPS-68", "HOPS-87", "HOPS-124", "HOPS-169", "HOPS-172", "HOPS-178", "HOPS-206", "HOPS-372", "HOPS-399", "HOPS-405", "HOPS-406"]

#code for writing to CSV

csvpath = jsondir + "angles.csv"
hasCsv = False
if not hasCsv:
    csvheader = ["Name", "Disk PA", "e Disk PA", "PA (b)", "PA (r)", "e PA (b)", "e PA (r)", "OA (b)", "OA (r)", "e OA (b)", "e OA (r)", "IC OA (b)", "IC OA (r)", "e IC OA (b)", "e IC OA (r)"]

    with open(csvpath, "w+", newline='') as acsv:
        csv_writer = csv.writer(acsv)
        csv_writer.writerow(csvheader)
        for i in l_list:
            pab = i['Position angle']['Blue']
            par = i['Position angle']['Red']
            if not (pab==None and par==None):
                name = i['Name']
                try:
                    dpa = i['dPa']
                except KeyError:
                    dpa = None
                try:
                    edpa = i['e_pa (deg)']
                except KeyError:
                    edpa = None
                epab = i['Position angle error']['Blue']
                epar = i['Position angle error']['Red']
                oab = i['Opening angle']['Blue']
                oar = i['Opening angle']['Red']
                eoab = i['Opening angle error']['Blue']
                eoar = i['Opening angle error']['Red']
                coab = i['Opening corrected']['Blue']
                coar = i['Opening corrected']['Red']
                ecoab = i['Opening corrected error']['Blue']
                ecoar = i['Opening corrected error']['Red']
                csv_writer.writerow([name,dpa,edpa,pab,par,epab,epar,oab,oar,eoab,eoar,coab,coar,ecoab,ecoar])

#Define lists from json files
name_list = val_list("Name", l_list)
tbol_list = val_list("Tbol", l_list)
lbol_list = val_list("Lbol", l_list)
dpa_list = val_list("dPa", l_list)
edpa_list = val_list("e_pa (deg)", l_list)
class_list = val_list("Class", l_list)

#Some dPa < -90, some dPa > 90, it doesn't make a difference so put them all in this range
for i in range(len(dpa_list)):
    try:
        if dpa_list[i] == -99.0:
            dpa_list[i] = None
        else:
            if dpa_list[i] > 90:
                dpa_list[i] -= 180
            if dpa_list[i] < -90:
                dpa_list[i] += 180
    except TypeError:
        pass
# print(min(list(filter(None,pa_list))))
# print(max(list(filter(None,pa_list))))

dpa_list_shift = []
for i in dpa_list:
    try:
        if i > 0:
            dpa_list_shift.append(i-90)
        else:
            dpa_list_shift.append(i+90)
    except TypeError:
        dpa_list_shift.append(None)

out_blue_list = angle_list("Position angle", ["Blue"], l_list)
out_red_list = angle_list("Position angle", ["Red"], l_list)

#Compute differences BEFORE correcting for 90+ degree angles
#not that it's foolproof
difference = []
for i in range(len(out_blue_list)):
    try:
        difference.append(abs(out_blue_list[i] - out_red_list[i]))
    except TypeError:
        difference.append(None)

#Move all position angles to the range (-90, 90)
for i in range(len(out_blue_list)):
    try:
        if out_blue_list[i] > 90:
            out_blue_list[i] = out_blue_list[i] - 180
    except TypeError:
        pass
    try:
        if out_red_list[i] > 90:
            out_red_list[i] = out_red_list[i] - 180
    except TypeError:
        pass

out_blue_error = angle_list("Position angle error", ["Blue"], l_list)
out_red_error = angle_list("Position angle error", ["Red"], l_list)

open_blue_list = angle_list("Opening angle", ["Blue"], l_list)
open_red_list = angle_list("Opening angle", ["Red"], l_list)

rb_difference = []
for i in range(len(open_blue_list)):
    try:
        rb_difference.append(open_blue_list[i] - open_red_list[i])
    except TypeError:
        rb_difference.append(None)

open_blue_error = angle_list("Opening angle error", ["Blue"], l_list)
open_red_error = angle_list("Opening angle error", ["Red"], l_list)

c_blue_list = angle_list("Opening corrected", ["Blue"], l_list)
c_red_list = angle_list("Opening corrected", ["Red"], l_list)

c_rb_difference = []
for i in range(len(c_blue_list)):
    try:
        c_rb_difference.append(c_blue_list[i] - c_red_list[i])
    except TypeError:
        c_rb_difference.append(None)

c_blue_error = angle_list("Opening corrected error", ["Blue"], l_list)
c_red_error = angle_list("Opening corrected error", ["Red"], l_list)

#convert all to numpy arrays and remove any outliers
difference, rb_difference, c_rb_difference, tbol_list, lbol_list, open_blue_error, open_blue_list, open_red_error, open_red_list, c_blue_error, c_red_error, c_blue_list, c_red_list, out_blue_list, out_blue_error, out_red_list, out_red_error, dpa_list_shift = toNp(difference, rb_difference, c_rb_difference, tbol_list, lbol_list, open_blue_error, open_blue_list, open_red_error, open_red_list, c_blue_error, c_red_error, c_blue_list, c_red_list, out_blue_list, out_blue_error, out_red_list, out_red_error, dpa_list_shift)
if excludeOutliers:
    name_list, tbol_list, lbol_list, dpa_list, edpa_list, class_list, difference, rb_difference, c_rb_difference, tbol_list, lbol_list, open_blue_error, open_blue_list, open_red_error, open_red_list, c_blue_error, c_red_error, c_blue_list, c_red_list, out_blue_list, out_blue_error, out_red_list, out_red_error, dpa_list_shift = removeValues(neglect, name_list, name_list, tbol_list, lbol_list, dpa_list, edpa_list, class_list, difference, rb_difference, c_rb_difference, tbol_list, lbol_list, open_blue_error, open_blue_list, open_red_error, open_red_list, c_blue_error, c_red_error, c_blue_list, c_red_list, out_blue_list, out_blue_error, out_red_list, out_red_error, dpa_list_shift)

#Extracting MASSES data from a txt file

massespath = "../MASSES.txt" #location of MASSES txt file
masses_data = []
tbol_masses = []
open_masses = []
open_error_masses = []
with open (massespath, 'r') as f :
    reader = csv.DictReader(f, fieldnames=["Name", "Tbol", "Position angle", "Opening angle", "Confidence"] , delimiter = "\t")
    for row in reader:
        masses_data.append(row)
        try:
            tbol_masses.append(float(row["Tbol"]))
        except ValueError:
            pass
        rows = row["Opening angle"].split("Â±")
        rows[0] = rows[0].split("<")
        try:
            open_masses.append(float(rows[0][0]))
            try:
                open_error_masses.append(float(rows[1]))
            except IndexError:
                open_error_masses.append(np.nan)
            except ValueError:
                open_error_masses.append(np.nan)
        except ValueError:
            try:
                open_masses.append(float(rows[0][1]))
                open_error_masses.append(np.nan)
            except IndexError:
                pass

masses_data = masses_data[1:]

#Plot Tbol vs opening angle (original and corrected)
#Remove all angle values of 0, that is, no information given
#70K determines class cutoff

ylimtop = np.nanmax(open_masses) *1.1
tbolOfInterest = []
lbolOfInterest = []
oaBOfInterest = []
oaBICOfInterest = []
oaROfInterest = []
oaRICOfInterest = []

#if the issue is with OA only, turning this into NaN might work better. PA and everything else are still relevant.
for i in range(len(name_list)):
    #knock out oab, oar, and oa all at once
    if (name_list[i] in oab_neglect) or (name_list[i] in oa_neglect):
        open_blue_list[i] = float('NaN')
        open_blue_error[i] = float('NaN')
        c_blue_list[i] = float('NaN')
        c_blue_error[i] = float('NaN')
    if (name_list[i] in oar_neglect) or (name_list[i] in oa_neglect):
        open_red_list[i] = float('NaN')
        open_red_error[i] = float('NaN')
        c_red_list[i] = float('NaN')
        c_red_error[i] = float('NaN')
    if (name_list[i] in ofInterest):
        tbolOfInterest.append(tbol_list[i])
        lbolOfInterest.append(lbol_list[i])
        oaBOfInterest.append(open_blue_list[i])
        oaBICOfInterest.append(c_blue_list[i])
        oaROfInterest.append(open_red_list[i])
        oaRICOfInterest.append(c_red_list[i])

fig = plt.figure(figsize = (13,5))

xlim = 11
ax1 = plt.subplot(1,2,1)
ax1.semilogx()
#VANDAM
ax1.errorbar(tbol_list, open_blue_list, yerr=open_blue_error, fmt='None', ecolor="b", elinewidth=1)
ax1.errorbar(tbol_list, open_red_list, yerr=open_red_error, fmt='None', ecolor="r", elinewidth=1)
ax1.scatter(tbol_list, open_blue_list, c="blue", s=10, label="VANDAM, blue lobe")
ax1.scatter(tbol_list, open_red_list, c="red", s=10, label="VANDAM, red lobe")
#MASSES
ax1.errorbar(tbol_masses, open_masses, yerr=open_error_masses, fmt='None', ecolor="gray", elinewidth=1)
ax1.scatter(tbol_masses, open_masses, c="gray", s=10, label="MASSES data", marker="v")
ax1.set_xlim(left=xlim)
ax1.set_ylim(bottom=0,top=ylimtop)
ax1.set_title("$T_{bol}$ vs Outflow Opening Angle")
ax1.set_xlabel("$T_{bol}$ (K)", fontsize=12)
ax1.set_ylabel("Outflow Opening Angle ($^\circ$)", fontsize=12)
ax1.legend(loc="upper left")

ax2 = plt.subplot(1,2,2)
ax2.semilogx()
#VANDAM
ax2.errorbar(tbol_list, c_blue_list, yerr=c_blue_error, fmt='None', ecolor="b", elinewidth=1)
ax2.errorbar(tbol_list, c_red_list, yerr=c_red_error, fmt='None', ecolor="r", elinewidth=1)
ax2.scatter(tbol_list, c_blue_list, c="blue", s=10, label="VANDAM, blue lobe")
ax2.scatter(tbol_list, c_red_list, c="red", s=10, label="VANDAM, red lobe")
#MASSES
ax2.errorbar(tbol_masses, open_masses, yerr=open_error_masses, fmt='None', ecolor="gray", elinewidth=1)
ax2.scatter(tbol_masses, open_masses, c="gray", s=10, label="MASSES data", marker="v")
ax2.set_xlim(left=xlim)
ax2.set_ylim(bottom=0,top=ylimtop)
ax2.set_title("$T_{bol}$ vs Corrected Outflow Opening Angle")
ax2.set_xlabel("$T_{bol}$ (K)", fontsize=12)
ax2.set_ylabel("Corrected Outflow Opening Angle ($^\circ$)", fontsize=12)
ax2.legend(loc="upper left")

path = figdir + "TbolOpeningBR.png"
plt.savefig(path, dpi=1000)
plt.close()

#Define lists where blue/red values are averaged
open_list = []
for i in range(len(open_blue_list)):
    open_list.append(np.average([open_blue_list[i], open_red_list[i]]))
c_list = []
for i in range(len(c_blue_list)):
    c_list.append(np.average([c_blue_list[i], c_red_list[i]]))
open_error = []
for i in range(len(open_blue_error)):
    open_error.append(0.5*np.sqrt(open_blue_error[i]**2 + open_red_error[i]**2))
c_error = []
for i in range(len(c_blue_error)):
    c_error.append(0.5*np.sqrt(c_blue_error[i]**2 + c_red_error[i]**2))

#Tbol vs opening angle (original and corrected)
#Remove all angle values of 0, that is, no information given

fig = plt.figure(figsize = (13,5))

#Choose whether to use the MCMC or scipy method of fitting: if MCMC, initialize MCMC parameters
use_emcee = True
if use_emcee:
    nwalkers = 100
    nburn = 100
    nsteps = 500
use_scipy = False

#Process lists so that emcee can handle them:
if True:
    #t1 refers to un-corrected, t2 refers to corrected. they may be different because c is unknown for some values
    tbol_list_t1 = [tbol_list[i] for i in range(len(tbol_list)) if not (isNaN(tbol_list[i]) or isNaN(open_list[i]))]
    open_list_t = [open_list[i] for i in range(len(tbol_list)) if not (isNaN(tbol_list[i]) or isNaN(open_list[i]))]
    open_error_t = [open_error[i] for i in range(len(tbol_list)) if not (isNaN(tbol_list[i]) or isNaN(open_list[i]))]
    class_list_t1 = [class_list[i] for i in range(len(tbol_list)) if not (isNaN(tbol_list[i]) or isNaN(open_list[i]))]
    tbol_list_t1 = np.asarray(tbol_list_t1)
    open_list_t = np.asarray(open_list_t)
    open_error_t = np.asarray(open_error_t)
    class_list_t1 = np.asarray(class_list_t1)
    tbol_list_t2 = [tbol_list[i] for i in range(len(tbol_list)) if not (isNaN(tbol_list[i]) or isNaN(c_list[i]))]
    c_list_t = [c_list[i] for i in range(len(tbol_list)) if not (isNaN(tbol_list[i]) or isNaN(c_list[i]))]
    c_error_t = [c_error[i] for i in range(len(tbol_list)) if not (isNaN(tbol_list[i]) or isNaN(c_list[i]))]
    class_list_t2 = [class_list[i] for i in range(len(tbol_list)) if not (isNaN(tbol_list[i]) or isNaN(c_list[i]))]
    tbol_list_t2 = np.asarray(tbol_list_t2)
    c_list_t = np.asarray(c_list_t)
    c_error_t = np.asarray(c_error_t)
    class_list_t2 = np.asarray(class_list_t2)

ylimtop = np.nanmax(open_masses) *1.1

ax1 = plt.subplot(1,2,1)
ax1.semilogx()
ax1.errorbar(tbol_list, open_list, yerr=open_error, fmt='None', ecolor="k", elinewidth=1)
ax1.scatter(tbol_list, open_list, c="black", s=10, label="VANDAM data")
ax1.set_xlim(left=xlim)
ax1.set_ylim(bottom=0,top=ylimtop)
ax1.set_title("$T_{bol}$ vs Outflow Opening Angle, Averaged")
ax1.set_xlabel("$T_{bol}$ (K)", fontsize=12)
ax1.set_ylabel("Outflow Opening Angle ($^\circ$)", fontsize=12)

#MASSES
ax1.errorbar(tbol_masses, open_masses, yerr=open_error_masses, fmt='None', ecolor="gray", elinewidth=1)
ax1.scatter(tbol_masses, open_masses, c="gray", s=10, label="MASSES data", marker="v")

#Using scipy fit
if use_scipy:
    #Using parameters from Dunham et. al *
    p0_t = [65, 5, 42]
    popt_t, pcov_t = curve_fit(f1, tbol_list, open_list, p0=p0_t, sigma=open_error, nan_policy="omit")
    xrange = np.linspace(xlim, max(list(filter(None, tbol_list))))
    sigma_t = np.sqrt(np.diag(pcov_t)) * 1
    ax1.plot(xrange, f1(xrange, popt_t[0], popt_t[1], popt_t[2]), 'c-', label="tanh fit")
    ax1.fill_between(xrange, f1(xrange, popt_t[0]+sigma_t[0], popt_t[1]+sigma_t[1], popt_t[2]+sigma_t[2]), f1(xrange, popt_t[0]-sigma_t[0], popt_t[1]-sigma_t[1], popt_t[2]-sigma_t[2]), color="c", alpha=0.2)

    p0_e = [65, 10, 26]
    popt_e, pcov_e = curve_fit(f2, tbol_list, open_list, p0=p0_e, sigma=open_error, nan_policy="omit")
    xrange = np.linspace(xlim, max(list(filter(None, tbol_list))))
    sigma_e = np.sqrt(np.diag(pcov_e)) * 1
    ax1.plot(xrange, f2(xrange, popt_e[0], popt_e[1], popt_e[2]), 'm-', label="(1-exp) fit")
    ax1.fill_between(xrange, f2(xrange, popt_e[0]+sigma_e[0], popt_e[1]-sigma_e[1], popt_e[2]-sigma_e[2]), f2(xrange, popt_e[0]-sigma_e[0], popt_e[1]+sigma_e[1], popt_e[2]+sigma_e[2]), color="m", alpha=0.2)

#Using emcee
if use_emcee:
    #Starting guesses are from Dunham et. al
    xfit = np.linspace(xlim, max(list(filter(None, tbol_list))))
    
    yfit_t, sig_t, samples_t, parameters_t = emcee_fit(tbol_list_t1, open_list_t, open_error_t, np.array([65, 5, 42]), nwalkers, nburn, nsteps, xfit, 1)
    yfit_e, sig_e, samples_e, parameters_e = emcee_fit(tbol_list_t1, open_list_t, open_error_t, np.array([65, 10, 26]), nwalkers, nburn, nsteps, xfit, 2)
    ax1.plot(xfit, yfit_t, c="#143232", label="tanh fit (VANDAM)", linewidth=3, alpha=0.8)
    ax1.plot(xfit, yfit_e, c="#006161", label="(1-exp) fit (VANDAM)", linewidth=3, alpha=0.8)
    # ax1.fill_between(xfit, yfit_t - sig_t, yfit_t + sig_t, color='c', alpha=0.2)
    # ax1.fill_between(xfit, yfit_e - sig_e, yfit_e + sig_e, color='m', alpha=0.2)
    ax1.legend(loc="upper left")

    print("tanh fit:\nθ = " + str(round(parameters_t[0],2)) + " deg, Tbol_shift = " + str(round(parameters_t[1],2)) + "K, Tbol_0 = " + str(round(parameters_t[2],2)) + "K")
    print("(1-exp) fit:\nθ = " + str(round(parameters_e[0],2)) + " deg, Tbol_shift = " + str(round(parameters_e[1],2)) + "K, Tbol_0 = " + str(round(parameters_e[2],2)) + "K")
    
ax2 = plt.subplot(1,2,2)
ax2.semilogx()
ax2.errorbar(tbol_list, c_list, yerr=c_error, fmt='None', ecolor="k", elinewidth=1)
ax2.scatter(tbol_list, c_list, c="black", s=10, label="VANDAM data")
ax2.set_xlim(left=xlim)
ax2.set_ylim(bottom=0,top=ylimtop)
ax2.set_title("$T_{bol}$ vs Corrected Outflow Opening Angle, Averaged")
ax2.set_xlabel("$T_{bol}$ (K)", fontsize=12)
ax2.set_ylabel("Corrected Outflow Opening Angle ($^\circ$)", fontsize=12)

#MASSES
ax2.errorbar(tbol_masses, open_masses, yerr=open_error_masses, fmt='None', ecolor="gray", elinewidth=1)
ax2.scatter(tbol_masses, open_masses, c="gray", s=10, label="MASSES data", marker="v")

#Using scipy fit
if use_scipy:
    #Using parameters from Dunham et. al *
    p0_t = [65, 5, 42]
    popt_t, pcov_t = curve_fit(f1, tbol_list, c_list, p0=p0_t, sigma=c_error, nan_policy="omit")
    xrange = np.linspace(xlim, max(list(filter(None, tbol_list))))
    sigma_t = np.sqrt(np.diag(pcov_t)) * 1
    ax2.plot(xrange, f1(xrange, popt_t[0], popt_t[1], popt_t[2]), 'c-', label="tanh fit")
    ax2.fill_between(xrange, f1(xrange, popt_t[0]+sigma_t[0], popt_t[1]+sigma_t[1], popt_t[2]+sigma_t[2]), f1(xrange, popt_t[0]-sigma_t[0], popt_t[1]-sigma_t[1], popt_t[2]-sigma_t[2]), color="c", alpha=0.2)

    p0_e = [65, 10, 26]
    popt_e, pcov_e = curve_fit(f2, tbol_list, c_list, p0=p0_e, sigma=c_error, nan_policy="omit")
    xrange = np.linspace(xlim, max(list(filter(None, tbol_list))))
    print(popt_e)
    sigma_e = np.sqrt(np.diag(pcov_e)) * 1
    ax2.plot(xrange, f2(xrange, popt_e[0], popt_e[1], popt_e[2]), 'm-', label="(1-exp) fit")
    ax2.fill_between(xrange, f2(xrange, popt_e[0]+sigma_e[0], popt_e[1]-sigma_e[1], popt_e[2]-sigma_e[2]), f2(xrange, popt_e[0]-sigma_e[0], popt_e[1]+sigma_e[1], popt_e[2]+sigma_e[2]), color="m", alpha=0.2)

    ax2.legend(loc="upper left")

#Using emcee
if use_emcee:
    #Starting guesses are from Dunham et. al
    xfit = np.linspace(xlim, max(list(filter(None, tbol_list))))
    yfit_ct, sig_ct, samples_ct, parameters_ct = emcee_fit(tbol_list_t2, c_list_t, c_error_t, np.array([65, 5, 42]), nwalkers, nburn, nsteps, xfit, 1)
    yfit_ce, sig_ce, samples_ce, parameters_ce = emcee_fit(tbol_list_t2, c_list_t, c_error_t, np.array([65, 10, 26]), nwalkers, nburn, nsteps, xfit, 2)
    ax2.plot(xfit, yfit_ct, c="#143232", label="tanh fit (VANDAM)", linewidth=3, alpha=0.8)
    ax2.plot(xfit, yfit_ce, c="#006161", label="(1-exp) fit (VANDAM)", linewidth=3, alpha=0.8)
    # ax2.fill_between(xfit, yfit_ct - sig_ct, yfit_ct + sig_ct, color='c', alpha=0.2)
    # ax2.fill_between(xfit, yfit_ce - sig_ce, yfit_ce + sig_ce, color='m', alpha=0.2)
    ax2.legend(loc="upper left")

    print("tanh fit (with inclination correction):\nθ = " + str(round(parameters_ct[0],2)) + " deg, Tbol_shift = " + str(round(parameters_ct[1],2)) + "K, Tbol_0 = " + str(round(parameters_ct[2],2)) + "K")
    print("(1-exp) fit (with inclination correction):\nθ = " + str(round(parameters_ce[0],2)) + " deg, Tbol_shift = " + str(round(parameters_ce[1],2)) + "K, Tbol_0 = " + str(round(parameters_ce[2],2)) + "K")

path = figdir + "TbolOpeningTEFits.png"
plt.savefig(path, dpi=1000)
plt.close()

#Power law fitting

#repeat, fitting classes separately (run the previous cell first)
fig = plt.figure(figsize = (13,5))
ylimtop = np.nanmax(tbol_masses) * 1.0
ylimbottom = np.nanmin(tbol_masses) * 0.9
nwalkers = 100
nburn = 100
nsteps = 500
cutoff_1 = 70 #70K

if True:
    tbol_0 = np.ndarray.flatten(tbol_list_t1)[(np.where(class_list_t1=='0.0'))[0]]
    tbol_1 = np.ndarray.flatten(tbol_list_t1)[(np.where(class_list_t1=="I"))[0]]
    tbol_c0 = np.ndarray.flatten(tbol_list_t2)[(np.where(class_list_t2=='0.0'))[0]]
    tbol_c1 = np.ndarray.flatten(tbol_list_t2)[(np.where(class_list_t2=="I"))[0]]
    open_0 = np.ndarray.flatten(open_list_t)[(np.where(class_list_t1=='0.0'))[0]]
    open_1 = np.ndarray.flatten(open_list_t)[(np.where(class_list_t1=="I"))[0]]
    open_c0 = np.ndarray.flatten(c_list_t)[(np.where(class_list_t2=='0.0'))[0]]
    open_c1 = np.ndarray.flatten(c_list_t)[(np.where(class_list_t2=="I"))[0]]
    open_error_0 = np.ndarray.flatten(open_error_t)[(np.where(class_list_t1=='0.0'))[0]]
    open_error_1 = np.ndarray.flatten(open_error_t)[(np.where(class_list_t1=="I"))[0]]
    open_error_c0 = np.ndarray.flatten(c_error_t)[(np.where(class_list_t2=='0.0'))[0]]
    open_error_c1 = np.ndarray.flatten(c_error_t)[(np.where(class_list_t2=="I"))[0]]
    #there has got to be a better way to do this

xfit_0 = np.linspace(xlim, cutoff_1)
xfit_1 = np.linspace(cutoff_1, np.nanmax(tbol_1))

ax1 = plt.subplot(1,2,1)
ax1.loglog()
ax1.plot(np.linspace(cutoff_1,cutoff_1), np.linspace(ylimbottom,ylimtop), 'k--', linewidth=1, label=str(cutoff_1)+"K")
ax2.errorbar(tbol_0, open_0, yerr=open_error_0, fmt='None', ecolor="#004600", elinewidth=1)
ax2.errorbar(tbol_1, open_1, yerr=open_error_1, fmt='None', ecolor="#209252", elinewidth=1)
ax1.scatter(tbol_0, open_0, c="#004600", s=10, label="VANDAM, Class 0")
ax1.scatter(tbol_1, open_1, c="#209252", s=10, label="VANDAM, Class I")
ax1.set_xlim(left=xlim)
ax1.set_ylim(bottom=ylimbottom,top=ylimtop)
ax1.set_title("$T_{bol}$ vs Outflow Opening Angle, Averaged")
ax1.set_xlabel("$T_{bol}$ (K)", fontsize=12)
ax1.set_ylabel("Outflow Opening Angle ($^\circ$)", fontsize=12)

#MASSES
ax1.errorbar(tbol_masses, open_masses, yerr=open_error_masses, fmt='None', ecolor="gray", elinewidth=1)
ax1.scatter(tbol_masses, open_masses, c="gray", s=10, label="MASSES data", marker="v")

#Using emcee
if use_emcee:
    #Starting guesses are from Dunham et. al
    yfit_0, sig_0, samples_0, parameters_0 = emcee_fit(tbol_0, open_0, open_error_0, np.array([np.pow(10,0.72),0.55,0]), nwalkers, nburn, nsteps, xfit_0, 3)
    yfit_1, sig_1, samples_1, parameters_1 = emcee_fit(tbol_1, open_1, open_error_1, np.array([np.pow(10,1.92),-0.05,0]), nwalkers, nburn, nsteps, xfit_1, 3)
    ax1.plot(xfit_0, yfit_0, c="#004600", label="power fit, Class 0", linewidth=3, alpha=0.8)
    ax1.plot(xfit_1, yfit_1, c="#209252", label="power fit, Class I", linewidth=3, alpha=0.8)
    # ax1.fill_between(xfit_0, yfit_0 - sig_0, yfit_0 + sig_0, color='c', alpha=0.2)
    # ax1.fill_between(xfit_1, yfit_1 - sig_1, yfit_1 + sig_1, color='m', alpha=0.2)
    ax1.legend(loc="upper left")

    print("Class 0 power fit (ax^b+c):\na = " + str(round(parameters_0[0],2)) + ", b = " + str(round(parameters_0[1],2)) + ", c = " + str(round(parameters_0[2],2)) + " deg")
    print("Class I power fit (ax^b+c):\na = " + str(round(parameters_1[0],2)) + ", b = " + str(round(parameters_1[1],2)) + ", c = " + str(round(parameters_1[2],2)) + " deg")

ax2 = plt.subplot(1,2,2)
ax2.loglog()
ax2.plot(np.linspace(cutoff_1,cutoff_1), np.linspace(ylimbottom,ylimtop), 'k--', linewidth=1, label=str(cutoff_1)+"K")
ax2.errorbar(tbol_c0, open_c0, yerr=open_error_c0, fmt='None', ecolor="#004600", elinewidth=1)
ax2.errorbar(tbol_c1, open_c1, yerr=open_error_c1, fmt='None', ecolor="#209252", elinewidth=1)
ax2.scatter(tbol_c0, open_c0, c="#004600", s=10, label="VANDAM, Class 0")
ax2.scatter(tbol_c1, open_c1, c="#209252", s=10, label="VANDAM, Class I")
ax2.set_xlim(left=xlim)
ax2.set_ylim(bottom=ylimbottom,top=ylimtop)
ax2.set_title("$T_{bol}$ vs Corrected Outflow Opening Angle, Averaged")
ax2.set_xlabel("$T_{bol}$ (K)", fontsize=12)
ax2.set_ylabel("Corrected Outflow Opening Angle ($^\circ$)", fontsize=12)

#MASSES
ax2.errorbar(tbol_masses, open_masses, yerr=open_error_masses, fmt='None', ecolor="gray", elinewidth=1)
ax2.scatter(tbol_masses, open_masses, c="gray", s=10, label="MASSES data", marker="v")

#Using emcee
if use_emcee:
    #Starting guesses are from Dunham et. al
    yfit_c0, sig_c0, samples_c0, parameters_c0 = emcee_fit(tbol_c0, open_c0, open_error_c0, np.array([np.pow(10,0.72),0.55,0]), nwalkers, nburn, nsteps, xfit_0, 3)
    yfit_c1, sig_c1, samples_c1, parameters_c1 = emcee_fit(tbol_c1, open_c1, open_error_c1, np.array([np.pow(10,1.92),-0.05,0]), nwalkers, nburn, nsteps, xfit_1, 3)
    ax2.plot(xfit_0, yfit_c0, c="#004600", label="power fit, Class 0", linewidth=3, alpha=0.8)
    ax2.plot(xfit_1, yfit_c1, c="#209252", label="power fit, Class I", linewidth=3, alpha=0.8)
    # ax2.fill_between(xfit_0, yfit_c0 - sig_c0, yfit_c0 + sig_c0, color='c', alpha=0.2)
    # ax2.fill_between(xfit_1, yfit_c1 - sig_c1, yfit_c1 + sig_c1, color='m', alpha=0.2)
    ax2.legend(loc="upper left")

    print("Class 0 power fit (ax^b+c):\na = " + str(round(parameters_c0[0],2)) + ", b = " + str(round(parameters_c0[1],2)) + ", c = " + str(round(parameters_c0[2],2)) + " deg")
    print("Class I power fit (ax^b+c):\na = " + str(round(parameters_c1[0],2)) + ", b = " + str(round(parameters_c1[1],2)) + ", c = " + str(round(parameters_c1[2],2)) + " deg")

path = figdir + "TbolOpeningPowerFits.png"
plt.savefig(path, dpi=1000)
plt.close()

#broken power law
#repeat, fitting classes separately and treating cutoff as a parameter (run the previous cell first)

fig = plt.figure(figsize = (13,5))
xfit = np.linspace(xlim, np.nanmax(tbol_list_t1))

ax1 = plt.subplot(1,2,1)
ax1.loglog()
ax1.errorbar(tbol_list, open_list, yerr=open_error, fmt='None', ecolor="k", elinewidth=1)
ax1.scatter(tbol_list, open_list, c="black", s=10, label="VANDAM data")
ax1.set_xlim(left=xlim)
ax1.set_ylim(bottom=ylimbottom,top=ylimtop)
ax1.set_title("$T_{bol}$ vs Outflow Opening Angle, Averaged")
ax1.set_xlabel("$T_{bol}$ (K)", fontsize=12)
ax1.set_ylabel("Outflow Opening Angle ($^\circ$)", fontsize=12)

#MASSES
ax1.errorbar(tbol_masses, open_masses, yerr=open_error_masses, fmt='None', ecolor="gray", elinewidth=1)
ax1.scatter(tbol_masses, open_masses, c="gray", s=10, marker="v", label="MASSES data")

#Using emcee
if use_emcee:
    #Starting guesses are from Dunham et. al
    yfit_b, sig_b, samples_b, parameters_b = emcee_fit(tbol_list_t1, open_list_t, open_error_t, np.array([54.3,0.55,-0.05,70]), nwalkers, nburn, nsteps, xfit, 4)
    ax1.plot(np.linspace(parameters_b[3],parameters_b[3]), np.linspace(ylimbottom,ylimtop), 'k--', linewidth=1, label=str(round(parameters_b[3],2))+"K")
    ax1.plot(xfit, yfit_b, c="#136C3A", label="Broken power law fit (VANDAM)", linewidth=3, alpha=0.8)
    # ax1.fill_between(xfit, yfit_b - sig_b, yfit_b + sig_b, color='c', alpha=0.2)
    ax1.legend(loc="upper left")

    print("Broken power fit (a(x/d)^b if x<d; a(x/d)^c if x>d):\na = " + str(round(parameters_b[0],2)) + ", b = " + str(round(parameters_b[1],2)) + ", c = " + str(round(parameters_b[2],2)) + ", d = " + str(round(parameters_b[3],2)) + " K")

ax2 = plt.subplot(1,2,2)
ax2.loglog()
ax2.errorbar(tbol_list, c_list, yerr=c_error, fmt='None', ecolor="k", elinewidth=1)
ax2.scatter(tbol_list, c_list, c="black", s=10, label="VANDAM data")
ax2.set_xlim(left=xlim)
ax2.set_ylim(bottom=ylimbottom,top=ylimtop)
ax2.set_title("$T_{bol}$ vs Corrected Outflow Opening Angle, Averaged")
ax2.set_xlabel("$T_{bol}$ (K)", fontsize=12)
ax2.set_ylabel("Corrected Outflow Opening Angle ($^\circ$)", fontsize=12)

#MASSES
ax2.errorbar(tbol_masses, open_masses, yerr=open_error_masses, fmt='None', ecolor="gray", elinewidth=1)
ax2.scatter(tbol_masses, open_masses, c="gray", s=10, marker="v", label="MASSES data")

#Using emcee
if use_emcee:
    #Starting guesses are from Dunham et. al
    yfit_b, sig_b, samples_b, parameters_b = emcee_fit(tbol_list_t2, c_list_t, c_error_t, np.array([54.3,0.55,-0.05,70]), nwalkers, nburn, nsteps, xfit, 4)
    ax2.plot(np.linspace(parameters_b[3],parameters_b[3]), np.linspace(ylimbottom,ylimtop), 'k--', linewidth=1, label=str(round(parameters_b[3],2))+"K")
    ax2.plot(xfit, yfit_b, c="#136C3A", label="Broken power law fit (VANDAM)", linewidth=3, alpha=0.8)
    # ax1.fill_between(xfit, yfit_b - sig_b, yfit_b + sig_b, color='c', alpha=0.2)
    ax2.legend(loc="upper left")

    print("Broken power fit (a(x/d)^b if x<d; a(x/d)^c if x>d):\na = " + str(round(parameters_b[0],2)) + ", b = " + str(round(parameters_b[1],2)) + ", c = " + str(round(parameters_b[2],2)) + ", d = " + str(round(parameters_b[3],2)) + " K")

path = figdir + "TbolOpeningBrokenPowerFits.png"
plt.savefig(path, dpi=1000)
plt.close()

#Now for the Lbol comparison. MASSES data does not have Lbol information as far as I know
#Here, Lbol is plotted as a colorbar on top of the Tbol vs opening angle plot

fig = plt.figure(figsize = (14,5))
xfit = np.linspace(xlim, np.nanmax(tbol_list_t1))

ax1 = plt.subplot(1,2,1)
ax1.loglog()
data_tl = ax1.scatter(tbol_list, open_list, c=lbol_list, s=10, norm="log", cmap="plasma")
clb = plt.colorbar(data_tl,label="$L_{bol}$")

ax1.set_xlim(left=xlim)
ax1.set_ylim(bottom=ylimbottom,top=ylimtop)
ax1.set_title("$T_{bol}$ vs Outflow Opening Angle, Averaged")
ax1.set_xlabel("$T_{bol}$ (K)", fontsize=12)
ax1.set_ylabel("Outflow Opening Angle ($^\circ$)", fontsize=12)

ax2 = plt.subplot(1,2,2)
ax2.loglog()

data_tl = ax2.scatter(tbol_list, c_list, c=lbol_list, s=15, norm="log", cmap="plasma")
plt.colorbar(data_tl,label="$L_{bol}$")
ax2.set_xlim(left=xlim)
ax2.set_ylim(bottom=ylimbottom,top=ylimtop)
ax2.set_title("$T_{bol}$ vs Corrected Outflow Opening Angle, Averaged")
ax2.set_xlabel("$T_{bol}$ (K)", fontsize=12)
ax2.set_ylabel("Corrected Outflow Opening Angle ($^\circ$)", fontsize=12)

path = figdir + "TbolOpeningLbolScale.png"
plt.savefig(path, dpi=1000)
plt.close()

#Lbol vs opening angle (original and corrected), blueshifted and redshifted lobes

fig = plt.figure(figsize = (13,5))
ylimtop = 180*1.1

ax1 = plt.subplot(1,2,1)
ax1.semilogx()
ax1.errorbar(lbol_list, open_blue_list, yerr=open_blue_error, fmt='None', ecolor="b", elinewidth=1)
ax1.errorbar(lbol_list, open_red_list, yerr=open_red_error, fmt='None', ecolor="r", elinewidth=1)
ax1.scatter(lbol_list, open_blue_list, c="blue", s=10, label="VANDAM, blue lobe")
ax1.scatter(lbol_list, open_red_list, c="red", s=10, label="VANDAM, red lobe")
ax1.set_ylim(bottom=0,top=ylimtop)
ax1.set_title("$L_{bol}$ vs Outflow Opening Angle")
ax1.set_xlabel("$L_{bol}$ ($L_\odot$)", fontsize=12)
ax1.set_ylabel("Outflow Opening Angle ($^\circ$)", fontsize=12)
ax1.legend()

ax2 = plt.subplot(1,2,2)
ax2.semilogx()
ax2.errorbar(lbol_list, c_blue_list, yerr=c_blue_error, fmt='None', ecolor="b", elinewidth=1)
ax2.errorbar(lbol_list, c_red_list, yerr=c_red_error, fmt='None', ecolor="r", elinewidth=1)
ax2.scatter(lbol_list, c_blue_list, c="blue", s=10, label="VANDAM, blue lobe")
ax2.scatter(lbol_list, c_red_list, c="red", s=10, label="VANDAM, red lobe")
ax2.set_ylim(bottom=0,top=ylimtop)
ax2.set_title("$L_{bol}$ vs Corrected Outflow Opening Angle")
ax2.set_xlabel("$L_{bol}$ ($L_\odot$)", fontsize=12)
ax2.set_ylabel("Corrected Outflow Opening Angle ($^\circ$)", fontsize=12)
ax2.legend()

path = figdir + "LbolOpeningBR.png"
plt.savefig(path, dpi=1000)
plt.close()

#Lbol vs opening angle (original and corrected), both lobes averaged
#Remove all angle values of 0, that is, no information given

fig = plt.figure(figsize = (13,5))
ylimtop = np.nanmax(open_list) *1.1

ax1 = plt.subplot(1,2,1)
ax1.semilogx()
ax1.errorbar(lbol_list, open_list, yerr=open_error, fmt='None', ecolor="k", elinewidth=1)
ax1.scatter(lbol_list, open_list, c="black", s=10, label="VANDAM")
ax1.set_ylim(bottom=0,top=ylimtop)
ax1.set_title("$L_{bol}$ vs Outflow Opening Angle, Averaged")
ax1.set_xlabel("$L_{bol}$ ($L_\odot$)", fontsize=12)
ax1.set_ylabel("Outflow Opening Angle ($^\circ$)", fontsize=12)
ax1.legend()

ax2 = plt.subplot(1,2,2)
ax2.semilogx()
ax2.errorbar(lbol_list, c_list, yerr=c_error, fmt='None', ecolor="k", elinewidth=1)
ax2.scatter(lbol_list, c_list, c="black", s=10, label="VANDAM")
ax2.set_ylim(bottom=0,top=ylimtop)
ax2.set_title("$L_{bol}$ vs Corrected Outflow Opening Angle, Averaged")
ax2.set_xlabel("$L_{bol}$ ($L_\odot$)", fontsize=12)
ax2.set_ylabel("Corrected Outflow Opening Angle ($^\circ$)", fontsize=12)
ax2.legend()

path = figdir  + "LbolOpening.png"
plt.savefig(path, dpi=1000)
plt.close()

#Compute Perason's R for a number of datasets:
#Tbol vs. opening angle for both classes 0 and I
#Lbol vs. opening angle for both classes 0 and I
#Values closer to 1 indicate higher correlation, values closer to -1 indicate inverse correlation, and values close to 0 indicate no correlation

if True:
    #l1 refers to un-corrected, l2 refers to corrected. they may be different because c is unknown for some values
    lbol_list_l1 = listMethod(lbol_list,lbol_list,open_list)
    open_list_l = listMethod(open_list,lbol_list,open_list)
    class_list_l1 = listMethod(class_list,lbol_list,open_list)
    lbol_list_l2 = listMethod(lbol_list,lbol_list,c_list)
    c_list_l = listMethod(c_list,lbol_list,c_list)
    class_list_l2 = listMethod(class_list,lbol_list,c_list)

if True:
    lbol_0 = np.ndarray.flatten(lbol_list_l1)[(np.where(class_list_l1=='0.0'))[0]]
    lbol_1 = np.ndarray.flatten(lbol_list_l1)[(np.where(class_list_l1=="I"))[0]]
    lbol_c0 = np.ndarray.flatten(lbol_list_l2)[(np.where(class_list_l2=='0.0'))[0]]
    lbol_c1 = np.ndarray.flatten(lbol_list_l2)[(np.where(class_list_l2=="I"))[0]]
    open_0l = np.ndarray.flatten(open_list_l)[(np.where(class_list_l1=='0.0'))[0]]
    open_1l = np.ndarray.flatten(open_list_l)[(np.where(class_list_l1=="I"))[0]]
    open_c0l = np.ndarray.flatten(c_list_l)[(np.where(class_list_l2=='0.0'))[0]]
    open_c1l = np.ndarray.flatten(c_list_l)[(np.where(class_list_l2=="I"))[0]]

rt0 = stats.pearsonr(tbol_0,open_0)
print("Pearson's R for Tbol and Opening Angle, Class 0:\t" + str(float(rt0[0])) + "; p-value: " + str(float(rt0[1])))
rt1 = stats.pearsonr(tbol_1,open_1)
print("Pearson's R for Tbol and Opening Angle, Class I:\t" + str(float(rt1[0])) + "; p-value: " + str(float(rt1[1])))
rtc0 = stats.pearsonr(tbol_c0,open_c0)
print("Pearson's R for Tbol and Inclination-Corrected Opening Angle, Class 0:\t" + str(float(rtc0[0])) + "; p-value: " + str(float(rtc0[1])))
rtc1 = stats.pearsonr(tbol_c1,open_c1)
print("Pearson's R for Tbol and Inclination-Corrected Opening Angle, Class I:\t" + str(float(rtc1[0])) + "; p-value: " + str(float(rtc1[1])))
print("")
rl0 = stats.pearsonr(lbol_0,open_0l)
print("Pearson's R for Lbol and Opening Angle, Class 0:\t" + str(float(rl0[0])) + "; p-value: " + str(float(rl0[1])))
rl1 = stats.pearsonr(lbol_1,open_1l)
print("Pearson's R for Lbol and Opening Angle, Class I:\t" + str(float(rl1[0])) + "; p-value: " + str(float(rl1[1])))
rlc0 = stats.pearsonr(lbol_c0,open_c0l)
print("Pearson's R for Lbol and Inclination-Corrected Opening Angle, Class 0:\t" + str(float(rlc0[0])) + "; p-value: " + str(float(rlc0[1])))
rlc1 = stats.pearsonr(lbol_c1,open_c1l)
print("Pearson's R for Lbol and Inclination-Corrected Opening Angle, Class I:\t" + str(float(rlc1[0])) + "; p-value: " + str(float(rlc1[1])))

#What is the distribution of differences between outflow position angles?
#To do this, we calculate how "misaligned" the position angles of the lobes are. Most of them average 0
#If, for example, the blue lobe is positioned at 40 degrees and its red lobe at 30 degrees, we say the PA difference is 10 degrees

windowlog = np.e

tbol_filter = tbol_list.copy()
lbol_filter = lbol_list.copy()
difference_filter = difference.copy()
names_filter = name_list.copy()

names_filter, difference_filter, tbol_filter, lbol_filter = removeValues(single_b, names_filter, names_filter, difference_filter, tbol_filter, lbol_filter)
names_filter, difference_filter, tbol_filter, lbol_filter = removeValues(single_r, names_filter, names_filter, difference_filter, tbol_filter, lbol_filter)
names_filter_pa = names_filter.copy()
difference_filter_pa = difference_filter.copy()

#First as a histogram...
for i in range(len(difference_filter)):
    if difference_filter[i] > 90:
        difference_filter[i] = 180 - difference_filter[i]
bs = 3
pa_diff_bins = np.arange(-1.5,179.5,bs)
pa_hist = plt.hist(difference_filter,bins=pa_diff_bins)
pa_hist_n = plt.hist(difference_filter,bins=pa_diff_bins,color="#2A76B1",ec="#2A76B1",density=True)

#use a Rician distribution instead!
#all PA difference values are positive
target_func = rician
p0_i = [2, -10, 10]
pa_x = pa_hist[1]+(bs/2)
popt, pcov = curve_fit(target_func,(pa_hist_n[1]+(bs/2))[:-1],pa_hist_n[0],p0=p0_i, maxfev=100000)
print(popt)

xrange = np.linspace(0,180,180)
plt.plot(xrange, (pa_hist[0][0]/pa_hist_n[0][0])*target_func(xrange, *popt), 'k-', label="Rician fit: b=" + str(round((popt[0]),2)) + "; loc=" + str(round((popt[1]),2)) + "; scale=" + str(round((popt[2]),2)))
plt.legend()
plt.title("Difference Between Position Angles of Lobes")
plt.xlabel("Difference (degrees)")
plt.ylabel("Number of Sources")

path = figdir + "PADiff_hist.png"
plt.savefig(path, dpi=1000)
plt.close()

#Tbol and Lbol vs red/blue outflow position angle difference
#Now as a scatter plot
fig = plt.figure(figsize = (13,5))

#Compare PA difference with Tbol
xlimt = 14
ax1 = plt.subplot(1,2,1)
ax1.semilogx()
ax1.scatter(tbol_filter, difference_filter, c="black", s=15)
if movingAverages:
    (xma,yma) = movingAverage(tbol_filter, difference_filter, windowlog, "log", llim=xlimt)
    ax1.plot(xma,yma,label="Moving average")
ax1.set_title("$T_{bol}$ vs Lobe Position Angle Difference")
ax1.set_xlabel("$T_{bol}$ (K)")
ax1.set_ylabel("Lobe Position Angle Difference ($^\circ$)")
ax1.set_ylim(bottom=0,top=90)
ax1.set_xlim(left=xlimt)

#Now compare PA difference with Lbol
xliml = 0.2
ax2 = plt.subplot(1,2,2)
ax2.semilogx()
ax2.scatter(lbol_filter, difference_filter, c="black", s=15)
if movingAverages:
    (xma,yma) = movingAverage(lbol_filter, difference_filter, windowlog, "log", llim=xliml)
    ax2.plot(xma,yma,label="Moving average")
ax2.set_title("$L_{bol}$ vs Lobe Position Angle Difference")
ax2.set_xlabel("$L_{bol}$ (K)")
ax2.set_ylabel("Lobe Position Angle Difference ($^\circ$)")
ax2.set_ylim(bottom=0,top=90)
ax2.set_xlim(left=xliml)

path = figdir + "PADiff_scatter.png"
plt.savefig("paper_figures/PADiff_scatter.png", dpi=1000)
plt.close()

#Find spread of PA difference
spread = (np.sqrt(np.nanvar(difference_filter)))
paDiffOfInterest = []
print(spread)
print("Sources with lobes with a large discrepancy in position angle (> 25 degrees)")
for i in range(len(names_filter_pa)):
    if abs(difference_filter[i]) > 25:
        print(names_filter_pa[i] + " " + str(difference_filter[i]))
    if names_filter_pa[i] in ofInterest:
        paDiffOfInterest.append(difference_filter_pa[i])

#What is the distribution of differences between opening angles (blue - red)?
#Again, OA difference is usually around 0. But some sources have lobes that are far wider than the other.
#For example, if the blueshifted lobe has an opening angle of 40 degrees and the redshifted lobe has an opening angle of 20 degrees,
#we say the OA difference is 20 degrees. If the lobes were reversed, it would be -20 degrees
#(No reason for the sign convention; but a sign difference is important here.)
windowlog = np.e

names_filter = name_list.copy()
rb_difference_filter = rb_difference.copy()
c_rb_difference_filter = c_rb_difference.copy()
tbol_filter = tbol_list.copy()
lbol_filter = lbol_list.copy()
names_filter, rb_difference_filter, c_rb_difference_filter, tbol_filter, lbol_filter = removeValues(oab_neglect, names_filter, names_filter, rb_difference_filter, c_rb_difference_filter, tbol_filter, lbol_filter)
names_filter, rb_difference_filter, c_rb_difference_filter, tbol_filter, lbol_filter = removeValues(oar_neglect, names_filter, names_filter, rb_difference_filter, c_rb_difference_filter, tbol_filter, lbol_filter)
names_filter, rb_difference_filter, c_rb_difference_filter, tbol_filter, lbol_filter = removeValues(oa_neglect, names_filter, names_filter, rb_difference_filter, c_rb_difference_filter, tbol_filter, lbol_filter)
names_filter, rb_difference_filter, c_rb_difference_filter, tbol_filter, lbol_filter = removeValues(single_b, names_filter, names_filter, rb_difference_filter, c_rb_difference_filter, tbol_filter, lbol_filter)
names_filter, rb_difference_filter, c_rb_difference_filter, tbol_filter, lbol_filter = removeValues(single_r, names_filter, names_filter, rb_difference_filter, c_rb_difference_filter, tbol_filter, lbol_filter)
names_filter_oa = names_filter.copy()

#First as a histogram
bs = 5
oa_diff_bins = np.arange(-182.5,182.5,bs)
oa_hist = plt.hist(rb_difference,bins=oa_diff_bins)

target_func = gaussian
p0_i = [10, 0, 5, 0]
pa_x = pa_hist[1]+(bs/2)
popt, pcov = curve_fit(target_func,(oa_hist[1]+(bs/2))[:-1],oa_hist[0],p0=p0_i, maxfev=100000)

xrange = np.linspace(-180,180,360)
plt.plot(xrange, gaussian(xrange,popt[0],popt[1],popt[2],popt[3]), 'k-', label="Gaussian fit: μ=" + str(round((popt[1]),2)) + "; σ=" + str(round((popt[2]),2)))
plt.legend()
plt.title("Difference Between Opening Angles of Lobes (Blue - Red)")
plt.xlabel("Difference (degrees)")
plt.ylabel("Number of Sources")

path = figdir + "OADiff_hist.png"
plt.savefig(path, dpi=1000)
plt.close()

#Tbol and Lbol vs red/blue outflow position angle difference
#Now as a scatter plot
fig = plt.figure(figsize = (13,5))

#First, Tbol vs OA difference
ax1 = plt.subplot(1,2,1)
ax1.semilogx()
ax1.scatter(tbol_filter, rb_difference_filter, c="black", s=15)
if movingAverages:
    (xma,yma) = movingAverage(tbol_filter, rb_difference_filter, windowlog, "log", llim=xlimt)
    ax1.plot(xma,yma,label="Moving average")
ax1.set_title("$T_{bol}$ vs Opening Angle Difference (Blue - Red)")
ax1.set_xlabel("$T_{bol}$ (K)")
ax1.set_ylabel("Opening Angle Difference ($^\circ$)")
ax1.set_ylim(bottom=-180,top=180)
ax1.set_xlim(left=xlimt)

#Next, Lbol vs OA difference
ax2 = plt.subplot(1,2,2)
ax2.semilogx()
ax2.scatter(lbol_filter, rb_difference_filter, c="black", s=15)
if movingAverages:
    (xma,yma) = movingAverage(lbol_filter, rb_difference_filter, windowlog, "log", llim=xliml)
    ax2.plot(xma,yma,label="Moving average")
ax2.set_title("$L_{bol}$ vs Opening Angle Difference (Blue - Red)")
ax2.set_xlabel("$L_{bol}$ (K)")
ax2.set_ylabel("Opening Angle Difference ($^\circ$)")
ax2.set_ylim(bottom=-180,top=180)
ax2.set_xlim(left=xliml)
# ax2.legend()

path = figdir + "OADiff_scatter.png"
plt.savefig(path, dpi=1000)
plt.close()

#Compute spread of OA difference
spread = (np.sqrt(np.nanvar(rb_difference_filter)))
oaDiffOfInterest = []
coaDiffOfInterest = []
print(spread)
print("Sources with lobes with a large discrepancy in opening angle (> " + str(spread) + " degrees)")
for i in range(len(names_filter_oa)):
    if abs(rb_difference_filter[i]) > spread:
        print(names_filter_oa[i] + " " + str(rb_difference_filter[i]))
    if names_filter_oa[i] in ofInterest:
        oaDiffOfInterest.append(rb_difference_filter[i])
        coaDiffOfInterest.append(c_rb_difference_filter[i])

#What is the distribution of differences between opening angles (blue - red)?
#This produces the same plots as earlier, but now with the inclination-corrected opening angles
windowlog = np.e

#First as a histogram
oa_hist = plt.hist(c_rb_difference,bins=oa_diff_bins)

target_func = gaussian
p0_i = [10, 0, 5, 0]
pa_x = pa_hist[1]+(bs/2)
popt, pcov = curve_fit(target_func,(oa_hist[1]+(bs/2))[:-1],oa_hist[0],p0=p0_i, maxfev=100000)

xrange = np.linspace(-180,180,360)
plt.plot(xrange, gaussian(xrange,popt[0],popt[1],popt[2],popt[3]), 'k-', label="Gaussian fit: μ=" + str(round((popt[1]),2)) + "; σ=" + str(round((popt[2]),2)))
plt.legend()
plt.title("Difference Between Inclination-Corrected Opening Angles")
plt.xlabel("Difference (degrees)")
plt.ylabel("Number of Sources")

path = figdir + "OACDiff_hist.png"
plt.savefig(path, dpi=1000)
plt.close()

#Tbol and Lbol vs inclination-corrected OA difference
#Scatter plots
fig = plt.figure(figsize = (13,5))

#Tbol vs inclination-corrected OA difference
ax1 = plt.subplot(1,2,1)
ax1.semilogx()
ax1.scatter(tbol_filter, c_rb_difference_filter, c="black", s=15)
if movingAverages:
    (xma,yma) = movingAverage(tbol_filter, c_rb_difference_filter, windowlog, "log", llim=xlimt)
    ax1.plot(xma,yma,label="Moving average")
ax1.set_title("$T_{bol}$ vs Inclination-Corrected Opening Angle Difference")
ax1.set_xlabel("$T_{bol}$ (K)")
ax1.set_ylabel("Opening Angle Difference ($^\circ$)")
ax1.set_ylim(bottom=-180,top=180)
ax1.set_xlim(left=xlimt)

#Lbol vs inclination-corrected OA difference
ax2 = plt.subplot(1,2,2)
ax2.semilogx()
ax2.scatter(lbol_filter, c_rb_difference_filter, c="black", s=15)
if movingAverages:
    (xma,yma) = movingAverage(lbol_filter, c_rb_difference_filter, windowlog, "log", llim=xliml)
    ax2.plot(xma,yma,label="Moving average")
ax2.set_title("$L_{bol}$ vs Inclination-Corrected Opening Angle Difference")
ax2.set_xlabel("$L_{bol}$ (K)")
ax2.set_ylabel("Opening Angle Difference ($^\circ$)")
ax2.set_ylim(bottom=-180,top=180)
ax2.set_xlim(left=xliml)

# print(np.sqrt(np.nanvar(rb_difference_filter)))

path = figdir + "OACDiff_scatter.png"
plt.savefig(path, dpi=1000)
plt.close()

#Now, let's see how OA difference and PA difference compare
#Another set of scatter plots

fig = plt.figure(figsize = (13,5))
windowlin = 10

names_filter_diff = names_filter_pa.copy() #449
names_filter_diff, difference_filter = removeValues(oab_neglect, names_filter_diff, names_filter_diff, difference_filter)
names_filter_diff, difference_filter = removeValues(oar_neglect, names_filter_diff, names_filter_diff, difference_filter)
names_filter_diff, difference_filter = removeValues(oa_neglect, names_filter_diff, names_filter_diff, difference_filter)

ax1 = plt.subplot(1,2,1)
ax1.scatter(difference_filter, rb_difference_filter, c="k", s=15)
if movingAverages:
    (xma,yma) = movingAverage(difference_filter,rb_difference_filter,windowlin,"lin")
    plt.plot(xma,yma,label="Moving average")
ax1.set_xlim(left=0,right=90)
ax1.set_ylim(bottom=np.nanmin(rb_difference_filter)*1.1,top=np.nanmax(rb_difference_filter)*1.1)
ax1.set_xlabel("Position angle difference")
ax1.set_ylabel("Opening angle difference")
ax1.set_title("PA difference vs OA difference")

#Inclination-corrected
ax2 = plt.subplot(1,2,2)
ax2.scatter(difference_filter, c_rb_difference_filter, c="k", s=15)
if movingAverages:
    (xma,yma) = movingAverage(difference_filter,c_rb_difference_filter,windowlin,"lin")
    plt.plot(xma,yma,label="Moving average")
ax2.set_xlim(left=0,right=90)
ax2.set_ylim(bottom=np.nanmin(rb_difference_filter)*1.1,top=np.nanmax(rb_difference_filter)*1.1)
ax2.set_xlabel("Position angle difference")
ax2.set_ylabel("Corrected opening angle difference")
ax2.set_title("PA difference vs corrected OA difference")

path = figdir + "PADiff_OADiff_scatter.png"
plt.savefig(path, dpi=1000)
plt.close()

#From here on out, we want to analyze whether protostellar disk and protostellar outflows (their position angle)
#must always be perpendicular.

#Let's start by subtracting 90 from the disk position angle
#Then, plot it against the outflow position angle

ax1 = plt.subplot()
x = np.linspace(-90,90,100)
y = np.linspace(-90,90,100)
ax1.plot(x,y,c="black")
ax1.scatter(dpa_list_shift, out_blue_list, c="blue", s=10, label="Blue lobe")
ax1.scatter(dpa_list_shift, out_red_list, c="red", s=10, label="Red lobe")
ax1.set_title("Disk Position Angle-90 vs Outflow Position Angle")
ax1.set_xlabel("Disk Position Angle-90 ($^\circ$)")
ax1.set_ylabel("Outflow Position Angle ($^\circ$)")
ax1.legend()

#Now we subtract them each other in order to find how much it differs from 0
#in which case 0 means perfectly perpendicular

b = []
r = []
for i in range(len(dpa_list_shift)):
    try:
        b_i = dpa_list_shift[i]-out_blue_list[i]
        if b_i > 90:
            b_i -= 180
        if b_i < -90:
            b_i += 180
        b.append(b_i)
    except TypeError:
        b.append(None)
    try:
        r_i = dpa_list_shift[i]-out_red_list[i]
        if r_i > 90:
            r_i -= 180
        if r_i < -90:
            r_i += 180
        r.append(r_i)
    except TypeError:
        r.append(None)
#something I need to fix
ind = np.where(name_list=="HOPS-1")[0][0]
b[ind] = None
ind = np.where(name_list=="HOPS-139")[0][0]
r[ind] = None
b = np.array(b,dtype=np.float64)
r = np.array(r,dtype=np.float64)

ax1 = plt.subplot()
x = np.linspace(-90,90,100)
x1 = np.linspace(-90,-90,100)
y = np.linspace(0,0,100)
ax1.plot(x,y,c="black")
ax1.scatter(out_blue_list, b, c="blue", s=10, label="Blue lobe")
ax1.scatter(out_red_list, r, c="red", s=10, label="Red lobe")

ax1.set_xlabel("Outflow Position Angle ($^\circ$)")
ax1.set_ylabel("Difference ($^\circ$)")
ax1.legend()
plt.title("Perpendicularity of Outflow to Disk")

path = figdir + "PADiskDifference.png"
plt.savefig(path, dpi=1000)
plt.close()

#Find new spread of perpendicularity measures
#Again, closer to 0 means more perfectly perpendicular

r_filter = list(filter(None,r))
r_filter = list(filter(lambda x: x < 50 and x > -50,r_filter))
b_filter = list(filter(None,b))
b_filter = list(filter(lambda x: x < 50 and x > -50,b_filter))

r_spread = math.sqrt(np.var(r_filter))
b_spread = math.sqrt(np.var(b_filter))

print(r_spread)
print(b_spread)

#Display the same perpendicularity information, now as a histogram

br_avg = []

for i in range(len(b)):
    if b[i] == None:
        if r[i] == None:
            br_avg.append(None)
        else:
            br_avg.append(r[i])
    elif r[i] == None:
        br_avg.append(b[i])
    else:
        br_avg.append((b[i]+r[i])/2)
br_avg = np.asarray(br_avg, dtype=np.float32)

bs = 4
pod_diff_bins = np.arange(-90.5,90.5,bs)
pod_hist = plt.hist(br_avg,bins=pod_diff_bins)

target_func = gaussian
p0_i = [20, 0, 10, 0]
pod_x = pod_hist[1]+(bs/2)
popt, pcov = curve_fit(target_func,(pod_hist[1]+(bs/2))[:-1],pod_hist[0],p0=p0_i, maxfev=100000)

xrange = np.linspace(-90,90,360)
plt.plot(xrange, gaussian(xrange,popt[0],popt[1],popt[2],popt[3]), 'k-', label="Gaussian fit: μ=" + str(round((popt[1]),2)) + "; σ=" + str(round((popt[2]),2)))

plt.legend(loc="upper left")
plt.title("Perpendicularity of Outflow to Disk")
plt.xlabel("Difference ($^\circ$)")
plt.ylabel("Number of Sources")

path = figdir + "PADiskDifference_hist.png"
plt.savefig(path, dpi=1000)
plt.close()

#How might perpendicularity depend on other parameters?
#Difference between disk position and outflow, vs Tbol, Lbol, and opening angle
#This set of plots displays blueshifted and redshifted lobes separately

fig = plt.figure(figsize = (20,5))
windowlin = 10
windowlog = np.e

actualtbol = [tbol_list[i] for i in range(len(b)) if not np.isnan(b[i])] + [tbol_list[i] for i in range(len(r)) if not np.isnan(r[i])]
actuallbol = [lbol_list[i] for i in range(len(b)) if not np.isnan(b[i])] + [lbol_list[i] for i in range(len(r)) if not np.isnan(r[i])]
actualopen = [open_blue_list[i] for i in range(len(b)) if not np.isnan(b[i])] + [open_red_list[i] for i in range(len(r)) if not np.isnan(r[i])]

#Compare perpendicularity with Tbol
ax1 = plt.subplot(1,3,1)
ax1.semilogx()
x = np.linspace(min(actualtbol),max(actualtbol),100)
y = np.linspace(0,0,100)
ax1.plot(x,y,c="black")
ax1.scatter(tbol_list, b, c="blue", s=10, label="Blue lobe")
ax1.scatter(tbol_list, r, c="red", s=10, label="Red lobe")
if movingAverages:
    (xma,yma) = movingAverage(np.concatenate((tbol_list,tbol_list)),np.concatenate((b,r)), windowlog, "log", llim=min(actualtbol), rlim=max(actualtbol))
    ax1.plot(xma,yma,label="Moving average")
ax1.set_title("Difference between dPa-90$^\circ$ and Outflow Angle vs $T_{bol}$")
ax1.set_xlabel("$T_{bol}$ (K)")
ax1.set_ylabel("Difference ($^\circ$)")
ax1.legend()

#Compare perpendicularity with Lbol
ax2 = plt.subplot(1,3,2)
ax2.semilogx()
x = np.linspace(min(actuallbol),max(actuallbol),100)
y = np.linspace(0,0,100)
ax2.plot(x,y,c="black")
ax2.scatter(lbol_list, b, c="blue", s=10, label="Blue lobe")
ax2.scatter(lbol_list, r, c="red", s=10, label="Red lobe")
if movingAverages:
    (xma,yma) = movingAverage(np.concatenate((lbol_list,lbol_list)),np.concatenate((b,r)), windowlog, "log", llim=min(actuallbol), rlim=max(actuallbol))
    ax2.plot(xma,yma,label="Moving average")
ax2.set_title("Difference between dPa-90$^\circ$ and Outflow Angle vs $L_{bol}$")
ax2.set_xlabel("$L_{bol}$ ($L_\odot$)")
ax2.set_ylabel("Difference ($^\circ$)")
ax2.legend()

#Compare perpendicularity with opening angle
ax3 = plt.subplot(1,3,3)
# ax3.semilogy()
x = np.linspace(min(actualopen),max(actualopen),100)
y = np.linspace(0,0,100)
ax3.plot(x,y,c="black")
ax3.scatter(open_blue_list, b, c="blue", s=10, label="Blue lobe")
ax3.scatter(open_red_list, r, c="red", s=10, label="Red lobe")
if movingAverages:
    (xma,yma) = movingAverage(np.concatenate((open_blue_list,open_red_list)),np.concatenate((b,r)),windowlin,"lin")
    ax3.plot(xma,yma,label="Moving average")
ax3.set_title("Difference between dPa-90$^\circ$ and Outflow Angle vs Opening Angle")
ax3.set_xlabel("Opening Angle ($^\circ$)")
ax3.set_ylabel("Difference ($^\circ$)")
ax3.legend()

path = figdir + "PADiskDifference_scatters.png"
plt.savefig(path, dpi=1000)
plt.close()

#How might perpendicularity depend on other parameters?
#Difference between disk position and outflow, vs Tbol, Lbol, and opening angle
#Now with blueshifted and redshifted lobes averaged

fig = plt.figure(figsize = (20,5))
windowlin = 10

#Redo the filtering out of single-lobed sources, etc.
br_avg_filter = br_avg.copy() #481
names_filter_diff = name_list.copy() #481
names_filter_diff, br_avg_filter = removeValues(single_b, names_filter_diff, names_filter_diff, br_avg_filter)
names_filter_diff, br_avg_filter = removeValues(single_r, names_filter_diff, names_filter_diff, br_avg_filter)
#down to 449
names_filter_diff, br_avg_filter = removeValues(oab_neglect, names_filter_diff, names_filter_diff, br_avg_filter)
names_filter_diff, br_avg_filter = removeValues(oar_neglect, names_filter_diff, names_filter_diff, br_avg_filter)
names_filter_diff, br_avg_filter = removeValues(oa_neglect, names_filter_diff, names_filter_diff, br_avg_filter)
#and 429 again

#Compare average perpendicuarity with Tbol
ax1 = plt.subplot(1,3,1)
# ax1.semilogx()
x = np.linspace(np.nanmin(difference_filter),np.nanmax(difference_filter),100)
y = np.linspace(0,0,100)
ax1.plot(x,y,c="black")
ax1.scatter(difference_filter, br_avg_filter, c="k", s=10)
if movingAverages:
    (xma,yma) = movingAverage(difference_filter,br_avg_filter,windowlin,"lin")
    ax1.plot(xma,yma,label="Moving average")
ax1.set_title("Perpendicularity vs Position Angle Difference")
ax1.set_xlabel("PA Difference ($^\circ$)")
ax1.set_ylabel("PA/DPA Difference ($^\circ$)")
# ax1.legend()

#Compare average perpendicularity with Lbol
ax2 = plt.subplot(1,3,2)
# ax2.semilogx()
x = np.linspace(np.nanmin(rb_difference_filter),np.nanmax(rb_difference_filter),100)
y = np.linspace(0,0,100)
ax2.plot(x,y,c="black")
ax2.scatter(rb_difference_filter, br_avg_filter, c="k", s=10)
if movingAverages:
    (xma,yma) = movingAverage(rb_difference_filter,br_avg_filter,windowlin,"lin")
    ax2.plot(xma,yma,label="Moving average")
ax2.set_title("Perpendicularity vs Opening Angle Difference")
ax2.set_xlabel("OA Difference ($^\circ$)")
ax2.set_ylabel("PA/DPA Difference ($^\circ$)")
# ax2.legend()

#Compare average perpendicularity with opening angle
ax3 = plt.subplot(1,3,3)
# ax3.semilogx()
x = np.linspace(np.nanmin(c_rb_difference_filter),np.nanmax(c_rb_difference_filter),100)
y = np.linspace(0,0,100)
ax3.plot(x,y,c="black")
ax3.scatter(c_rb_difference_filter, br_avg_filter, c="k", s=10)
if movingAverages:
    (xma,yma) = movingAverage(c_rb_difference_filter,br_avg_filter,windowlin,"lin")
    ax3.plot(xma,yma,label="Moving average")
ax3.set_title("Perpendicularity vs Corrected Opening Angle Difference")
ax3.set_xlabel("Corrected OA Difference ($^\circ$)")
ax3.set_ylabel("PA/DPA Difference ($^\circ$)")
# ax3.legend()

path = figdir + "PADiskDifference_avg_scatters.png"
plt.savefig(path, dpi=1000)
plt.close()