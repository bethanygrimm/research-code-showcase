# Star Formation Project Code

## About the Project
This is the code I worked on in my star formation group at UT Austin. The dataset included 330 protostars in the Orion nebula, ALMA data publically available [here](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi%3A10.7910%2FDVN%2FCWTUCI&version=&q=&fileTypeGroupFacet=%22FITS%22&fileAccess=&fileTag=&fileSortField=&fileSortOrder=&tagPresort=true&folderPresort=true). I assembled moment-0 maps of the protostars (`moment_map_script.py`), then calculated position angle and opening angle of each resolvable protostellar outflow (`position_angle_script.py`), and finally performed statistical analysis on the obtained angles (`angle_statistics.py`). This aims to answer the questions of 1\) whether protostellar outflows must always be perpendicular to their disks (generally they are, within a certain spread, though there are some interesting outliers) and 2\) if protostellar opening angle can be correlated with any other protstellar properties (yes! There is a positive correlation between opening angle and bolometric temperature, and hence opening angle and stellar age).

Any relevant publications will be linked here!

## Dependencies
Due to the size of ALMA data, the data required to run these files are **not** included in the repository. Relevant data of protostellar properties are located in `table6.txt`, `table8.txt`, and `tableE.txt`, in the folder `data`

If a `figures` directory is initialized in `star_formation`, `angle_statistics.py` can be run, and the resultant figures will show up there.

The `sample_results` shows sample plots from `moment_map_script.py`, `position_angle_script.py` (for source HH212MMS), and `angle_statistics.py`.

## Included Files
### Python
- `labels.py`: These functions handle all the json files storing the protostar data: it reads input from txt files, and updates each source with its angle information when needed.
- `moment_map_script.py`: This script displays moment-0 maps for up to four sources, given a user-provided frequency range. Completely optional to run.
- `position_angle_script.py`: This script takes the moment-0 map for a source, masks out noise, and performs iterative Gaussian fits to calculate the protostellar outflow position angle and opening angles. This one requires a lot of user input and finagling to lock in the fit parameters (and only handles one source at a time); it's better suited as a Jupyter notebook.
- `angle_statistics.py`: This script creates plots and performs fits in order to analyze the relationship between opening angle and other parameters. It also compares the protostellar disk position angle and outflow position angle to see if they are perpendicular, and identifies any sources with large discrepancies between lobes. After initializing a `figures` directory in `star_formation`, this one can be run as-is.