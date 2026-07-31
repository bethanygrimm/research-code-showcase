import numpy as np

long0_gsr=87.8
lat0_gsr=1.7
Vapex_gsr=232.3

long0_lgsr=np.array([91., 93.0, 95.])
lat0_lgsr=np.array([-2., -4.0, -6.])
Vapex_lgsr=np.array([311., 316.0, 321.])

#longitude and latitude values in radians
long0r_lgsr=long0_lgsr*np.pi/180.0
lat0r_lgsr=lat0_lgsr*np.pi/180.0

def velConversion(longd,latd,VLSR):
    '''
    Description:
        this is the NED 2/17/2012 definition
        ned.ipac.caltech.edu/help/velc_help.html  see my Wright's cloud M31 notebook p. 87
        GIVEN VLSR FIND OTHER VELOCITIES USING THE NED DEFINITIONS

        code by Jay Lockman
        edited by Bethany Grimm
        VHEL from Robert Braun & Thilker

    Inputs:
        longd (float): the galactic longitude (l) of the object or scan, in degrees
        latd (float): the galactic latitude (b) of the object or scan, in degrees
        VLSR (float): the LSR velocity of the object or scan, in km/s

    Returns:
        VHEL (float): the heliocentric velocity of the object or scan, in km/s
        VGSR (float): the GSR velocity of the object or scan, in km/s
        VLGSR (float): the LGSR velocity of the object or scan, in km/s. Derived from VHEL.
        VLGSR1 (float): the LGSR velocity of the object or scan, in km/s. Derived from VGSR.
    '''
    longr=longd*np.pi/180.0
    latr=latd*np.pi/180.0
    cos_long = np.cos(longr)
    sin_long = np.sin(longr)
    sin_lat  = np.sin(latr)
    cos_lat  = np.cos(latr)
    VHEL = VLSR-9.0*cos_long*cos_lat-12.0*sin_long*cos_lat-7.0*sin_lat
    VGSR = VLSR + 220.0*sin_long*cos_lat
    VLGSR1 = VGSR - 62.*cos_long*cos_lat+40.*sin_long*cos_lat -35.*sin_lat
    VLGSR=VHEL+Vapex_lgsr[1]*(np.sin(latr)*np.sin(lat0r_lgsr[1]) + 
        np.cos(latr)*np.cos(lat0r_lgsr[1])*np.cos(longr-long0r_lgsr[1]))
    # print(longd,latd,VLSR,VHEL,VGSR,VLGSR,VLGSR1)
    return VHEL, VGSR, VLGSR, VLGSR1