# TADF_DRC
These files can be used to obtain the decay rate constants for TADF decays using the bi-exponential or tri-exponential models
The different files are listed below, along with a brief description of them:
a.	‘2model.py’: This file contains the Python script to calculate decay rate constants and their associated uncertainties with the Bi-exponential model. The program requires the degassed decay datafile evaluated with ‘Transient_PL_Udur.py’, as well as its respective aerated PLQY and aerated-degassed ratio values.
b.	‘3model.py’: This file contains the Python script to calculate decay rate constants and their associated uncertainties with the Tri-exponential model. The program requires the degassed and aerated decays evaluated with ‘Transient_PL_Udur.py’, as well as their respective aerated PLQY and aerated-degassed ratio values.
c.	‘Example_data’: This folder includes example data files ‘P190AD_evaluated.txt’ and ‘P190DD_evaluated.txt’, which contain the 190 mg/ml sample aerated and degassed decays, respectively. These can be used for code testing. The PLQY and aerated-degassed ratio for these measurements are 0.2518 and 1.38, respectively.


The following files were also used in this work, but are not included in the repository as they were provided externally:
a.	‘EX15.TXT’: This file contains the integration times and delay times utilized for time-resolved measurements. Provided by OEM group.
b.	‘Transient_PL_Udur.py’: This file contains the code utilized to calculate the time-resolved decay from the individual spectra obtained with the iCCD camera. Provided by OEM group.
