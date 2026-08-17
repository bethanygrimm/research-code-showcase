import numpy as np 

# longM31 = 10.68*np.pi/180.0
# latM31 = 41.27*np.pi/180.0
LONGM31 = 121.17*np.pi/180.0
LATM31 = -21.57*np.pi/180.0

def angleFromM31(longd,latd):
    '''
    Description:
        This calculates the distance of an object from M31 in degrees, as viewed from Earth and projected onto the sky.
        If the longitude of the object exceeds the longitude of M31, it returns a positive value, and vice versa
        Could be generalized to be a general distance function

    Inputs:
        longd (float): the galactic longitude (l) of the object or scan, in degrees
        latd (float): the galactic latitude (b) of the object or scan, in degrees

    Returns:
        angleFromM31 (float): the distance of the object from M31 in degrees
    '''
    longr=longd*np.pi/180.0
    latr=latd*np.pi/180.0
    ldiff = longr - LONGM31
    sign = int(ldiff/abs(ldiff))
    return sign * np.arccos(np.sin(LATM31)*np.sin(latr) + np.cos(LATM31)*np.cos(latr)*np.cos(ldiff)) * 180.0/np.pi