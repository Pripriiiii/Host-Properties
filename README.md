# Host-Properties
This repository stores and updates datasets and relevant notebooks for host galaxy property computations.

## Repository Contents

### Jupyter Notebooks

- **Host_Properties**:  
  This notebook covers essential procedures for computing host properties (LOGMASS, (U-R) colour, dDLR).

### Python Cripts

- **Photometry**:
  Separated script for computing global photometry only.

- **Store_phot**:
  Separated script for saving photometry data and writing it into text files used for sedfitting.

- **UR_Colour**:
  Separated script for computing (U-R) colour only.


### Dataset

- **DEBASS_allHost**:
Contains RA&DEC for SN&Host; and host galaxy property measurements of DEBASS sample. (Currently 321 hosts. Note that dDLR is in processing.)

- **FN_Host_Properties**:
Contains RA&DEC for SN&Host; and host galaxy property measurements of the Foundation sample.

- **In_Both**:
Contains DEBASS SNe Ia that image cut-outs are available in both PS1 and SDSS surveys.

- **Only_PS1**:
Contains DEBASS SNe Ia that image cut-outs are available only in the PS1 survey.

- **FN_PS1_ONLY**:
Contains FN SNe Ia that image cut-outs are available only in the PS1 survey.

- **FN_PS1_SDSS**:
Contains FN SNe Ia that image cut-outs are available in both PS1 and SDSS surveys.

- **Model Prediction**:
Contains extracted data of the Childress2014 Mass Distribution Model.

- **DEBASS_shape_low_threshold, FN_shape_low_threshold2**:
Contains shape parameters with a low threshold for DEBASS and FN samples, which are used to measure dDLR. Note that we find larger values for dDLR than expected.

- **alpha_lyr_stis_005.ascii**, **sr.dat**, **sux.dat**:
Bessell filter curves in SALT2 are used for host (U-R) colour computation.


## How to Use

Clone this repository to your local machine:
   ```bash
   git clone https://github.com/Pripriiiii/DEBASS-Host.git
