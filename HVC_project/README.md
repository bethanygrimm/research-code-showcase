# HVC Project Code

## About the Project
This is the code I worked on during my summer undergraduate research appointment at Green Bank Observatory. I worked on a project **analyzing high-velocity clouds** (HVCs) between the galaxies M31 (Andromeda) and M33 (Triangulum). Data was collected using the Green Bank Telescope (GBT), which mapped **neutral hydrogen** abundances in the sky. The HVCs of interest between M31 and M33 exist in a line between the two galaxies and are not rotationally bound to either galaxy. The HVCs' origins are still unknown, so this code exists to **reduce GBT data** and **analyze the dynamics and parameters of these HVCs**.

Once the report is published, it will be linked here!

## Dependencies
Due to the nature of GBT data, the data required to run these files are **not** included in the repository. However, the folder `HVCData` includes data from various HVCs.
`.pro` files are executed with GBTIDL (a variant of IDL), and `.py` files are executed with Python.

If a `figures` directory is initialized in `HVC_project`, `cloudPlots.py` and `findTcal.py` can be run, and the resultant figures will show up there.

## Included Files
### GBTIDL
- `auxiliary.pro`: Auxiliary GBTIDL functions for `align.pro` and `edgeoffkeep_hanning.pro`
- `fluxcheck.pro`: Allows a user to go through calibration scans and manually select regions to be baselined and subtracted (in other words, this is how we subtract noise from a telescope scan). Then, the function extracts the longitude/latitude and maximum intensity of each scan and outputs it into a .txt file.
- `edgeoffkeep_hanning.pro`: This function is where the bulk of GBT data reduction takes place! It goes through all scans of specified files, "subtracts" reference data (taken to be the scans on the edge of the map), averages, and Hanning smooths the data.
- `align.pro`: When taking telescope observations, minute changes in velocity are inevitable. This function corrects for them by "aligning" spectra from different scans together.
### Python
- `angleFromM31.py`: This contains a function that calculates the distance of an object from M31 (in degrees, as projected on the plane of the sky)
- `velConversion.py`: This contains a function that converts Local Standard of Rest (LSR) velocity to heliocentric (HEL) velocity, Galactic Standard of Rest (GSR) velocity, and Local Galactic Standard of Rest (LGSR) velocity
- `findTcal.py`: Given a set of data from calibration scans (as output of `fluxcheck.pro`), this finds the *actual* calibration temperatures for the telescope.
- `cloudParameters.py`: Given a set of Gaussian fit logs from CARTA, this compiles log files into a machine-readable CSV. (The Gaussian fits are 1-D fits at different spatial points on the HVC, fit to spectral features in the velocity domain.)
- `cloudPlots.py`: Now knowing the parameters of the HVC from `cloudParameters.py`, this computes values of interest, such as cloud area (viewed on the plane of the sky), mass, etc. This creates many supplementary plots to visualize the HVCs and their properties.
