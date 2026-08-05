# PySMACKS
## Overview

`PySMACKS` is a Python framework for analysis of single-molecule FRET (smFRET) data, from the raw microscopy movies to the corrected FRET histogram and kinetic analysis.

The initial motivation for `PySMACKS` was to port the [SMACKS framework](https://www.singlemolecule.uni-freiburg.de/software/smacks), used within our group, into a modern Python framework. Additional features have then been integrated 
in order to provide a complete analytical pipeline within a comprehensive framework. The choice of Python was motivated by the vast community of Python developers in smFRET and biophysics fields, 
and the availability of state-of-the-art computational packages, especially in machine-learning and deep-learning methods.

`PySMACKS` includes features that are necessary for our ongoing research projects, and for shich we could not find existing solutions in the smFRET fields. For imstance, the ability to correct for the microscope stage drift is critical for long
smFRET experiments, hence the development and integration of this feature into `PySMACKS`.

`PySMACKS` is modular in its architecture and is meant to follow the FAIR data principles. Especially, the choice of the [OpenFRET format](https://pypi.org/project/openfret) for the storage of the traces datasets is motivated by 
an interoperability purpose, with the idea that this open format makes it easy to translate to another open format (for instance netCDF) without loss of metadata. 
Open science is a core value of the development team of `PySMACKS` and we hope that the smFRET community will benefit from this framework and build upon it.

`PySMACKS` is provided as a Python package and is freely available from [https://pypi.org/project/PySMACKS/](https://pypi.org/project/PySMACKS/).

## Table of Contents

  - [Key Features](#key-features)
  - [Installation](#installation)
  - [Workflow](#workflow)
    - [Channels registration](#channels-registration)
	- [Traces extraction](#traces-extraction)
	- [Traces filtering](#trace-filtering)
	- [FRET histogram](#stoichiometry-efficiency-fret-histogram)
	- [Kinetic analysis](#kinetic-analysis)
  - [Citation](#citation)
    - [Persistent identifiers](#persistent-identifiers)
  - [License](#license)
  - [Contact](#contact)
  
## Key Features

In the current version (1.0.0):

- [Channels registration](registration.md)
    - including Optical Flow or Affine transformation correction methods
    - QC plot of the registration results
- [Traces extraction](extraction.md)
    - including ALEX- and non-ALEX- smFRET data
    - Microscope stage drift correction
    - local background correction over time using the median of surrounding pixels, or advanced strategy (Total Variation, Minimum of Total Variation)
    - Storage of traces within OpenFRET JSON files
- [Traces viewer and filtering](filtering.md)
    - interactive plot
    - manual labelling tool for downstream analysis of highlighted sections and/or ML-DL curation of training datasets
- [Stoichiometry-Efficiency FRET histogram](FRET_hist.md)
    - automatic calculation of FRET correction factors
- [Kinetic analysis](kinetics.md)
    - Trace-by-Trace HMM
    - Ensemble HMM
- Metadata viewer

## Installation

`PySMACKS` is available as a Python package from [https://pypi.org/project/PySMACKS/](https://pypi.org/project/PySMACKS/)

The suggested installation method is to use pip within a virtual environment:

  1. Set a virtual environment
  
    ```
    python -m venv NAME_OF_YOUR_VENV
	```
	
	And activate it
	
	```
    NAME_OF_YOUR_VENV\Scripts\activate
	```
	
  2. Install the PySMACKS package using pip

    ```
    pip install PySMACKS
    ```
	
  3. Launch the PySMACKS GUI using the command
  
    ```
    PySMACKS_GUI
    ```
	
### Dependencies

| Libraries | Licence | Note  |
|----------|----------|----------|
| [numpy](https://pypi.org/project/numpy/)  |   BSD-3-Clause    |         |
| [matplotlib](https://pypi.org/project/matplotlib/)  |    Python Software Foundation License    |       |
| [pandas](https://pypi.org/project/pandas/)  |    BSD 3-Clause    |       |
| [PyQt5](https://pypi.org/project/PyQt5/)  |   GPL v3     |       |
| [scipy](https://pypi.org/project/scipy/)  |   BSD License     |       |
| [scikit-image](https://pypi.org/project/scikit-image/)  |    BSD License    |       |
| [big-fish](https://pypi.org/project/big-fish/)  |   BSD 3-Clause     |   [paper](https://doi.org/10.1101/2021.07.20.453024)    |
| [hmmlearn](https://pypi.org/project/hmmlearn/)  |    BSD License    |       |
| [opencv-python](https://pypi.org/project/opencv-python/)  |    Apache 2.0    |       |
| [openfret](https://pypi.org/project/openfret/)  |   MIT     |    [GitHub repo](https://github.com/simol-lab/OpenFRET)   |
| [pillow](https://pypi.org/project/pillow/)  |    MIT-CMU    |       |
| [ruptures](https://pypi.org/project/ruptures/)  |   BSD-2-Clause     |   [paper](https://doi.org/10.1016/j.sigpro.2019.107299) * not used in the current version of `PySMACKS`*   |
| [seaborn](https://pypi.org/project/seaborn/)  |   BSD License     |       |
| [trackpy](https://pypi.org/project/trackpy/)  |        |  [GitHub repo](https://github.com/soft-matter/trackpy)     |
| [tvdcondat2013](https://pypi.org/project/TVDCondat2013/)  |   CeCILL FREE SOFTWARE     |   [GitHub repo](https://github.com/bgailleton/TVD_Condat2013)    |
| [uuid7](https://pypi.org/project/uuid7/)  |    MIT    |       |

## Workflow

The analysis pipeline can be divided into two separate sections:

  - [Traces extraction](#traces-extraction)
  - [Traces analysis](#traces-analysis)
  
`PySMACKS` is meant to be modular and interoperable: raw movies can be processed by `PySMACKS` to extract the single-molecule traces, and analyzed within `PySMACKS`; 
or these traces can be exported to JSON files and processed using another framework.

Basically, the data management within `PySMACKS` follows the logic: `TIFF --> JSON (OpenFRET format)`.
  
### Traces extraction

#### Input data

In the current version (1.0.0), `PySMACKS` only handles TIFF stacks as data input for the traces extraction. Especially, the two channels (donor and acceptor) should be in separate TIFF files.

Additionally, for memory efficiency, the full movies should be split into multiple submovies, (usually a submovie < 500 MB), to avoid saturation of the RAM.

#### ALEX smFRET

`PySMACKS` supports both ALEX and non-ALEX smFRET data; you can specify it during the traces extraction process. 

If your data are ALEX, it is expected that the acceptor submovies frames alternate between acceptor-acceptor frame (ALEX) and donor-acceptor frame (FRET), with the first frame of the movie being ALEX.

#### Channels registration

You can correct the chromatic aberrations and channels misalignments using optical beads TIFF movies (donor and acceptor separate files) and the [registration module](registration.md) of `PySMACKS`.

Currently, two methods can be used for the registration:

  - Optical flow
  - Affine transformation
  
The results of the registration step are stored within numpy files (*.npy)

#### Extraction step

The individual traces are [detected](extraction.md) using a LoG filter and a tracking algorithm.

The drift of the microscope stage can be corrected at this stage.

The individual traces are stored within a JSON file (compressed as .json.zip) following the [OpenFRET format](https://pypi.org/project/openfret) (source repo: [https://github.com/simol-lab/OpenFRET](https://github.com/simol-lab/OpenFRET))

All metadata are stored within the OpenFRET structure, including the local background intensity over time for each individual traces.

### Traces analysis

#### Input data

For the trace analysis (filtering, FRET histogram, kinetics), the input data is the JSON file (compressed as .json.zip) following the [OpenFRET format](https://pypi.org/project/openfret).

The outputs of the analysis pipeline are saved within the same .json file through its metadata.

#### Trace filtering

`PySMACKS` contains a [trace viewer](filtering.md), which you can use to visualize traces one-by-one and manually filter them. You can also highlight relevant sections of these traces using an interactive tool.

Especially, this labelling interactive tool can be used to curate datasets for ML/DL training purposes.

*NB: the development of a fully integrated deep learning approach for traces filtering in under development and will be implemented in future releases of PySMACKS.*

#### Stoichiometry-Efficiency FRET histogram

If your data are ALEX-smFRET, after filtering and labelling of FRET, Donor-Only and Acceptor-Only sections of individual traces, you can plot the stoichiometry-efficiency (SE) FRET histogram using this [module](FRET_hist.md).

Especially, the calculation of the FRET corrections factors (alpha, delta, beta and gamma) is fully automated with intermediate quality check by the user.

A corrected FRET efficiency histogram can be computed.

#### Kinetic analysis

The [last module](kinetics.md) of the pipeline is dedicated to the extraction of the kinetic rates governing states transitions, using a Hidden-Markov Model (HMM) framework.

In details, this HMM approach relies on the [SMACKS framework](https://www.singlemolecule.uni-freiburg.de/software/smacks) (original study: [10.1016/j.bpj.2016.08.023](https://doi.org/10.1016/j.bpj.2016.08.023)), initially 
developed in IGOR language, and ported to Python within `PySMACKS`.

Users can conduct HMM inferences on a Trace-by-Trace basis (one individual HMM model optimized for each individual traces), or on an Ensemble way.

The outputs are rates matrix (TbT or Ensemble), related to the discrete states transition rates of the single-molecules.

## Citation

Molcrette, B., & Schmid, S. (2026). PySMACKS (Version 1.0.0b2) [Computer software]. https://github.com/molcretb/PySMACKS

### Persistent identifiers

  - [RRID:SCR_028788](https://scicrunch.org/resolver/RRID:SCR_028788)
  - [biotools:pysmacks](https://bio.tools/pysmacks)
  - [swh:1:snp:55dd878eb83b5de7a64e0f104e9db2dacf2c6669;origin=https://github.com/molcretb/PySMACKS](https://archive.softwareheritage.org/swh:1:snp:55dd878eb83b5de7a64e0f104e9db2dacf2c6669;origin=https://github.com/molcretb/PySMACKS)

## License

Distributed under the [MIT License](https://opensource.org/license/MIT).

## Contact

Please use the [GitHub issues](https://github.com/molcretb/PySMACKS/issues) to report any problem with `PySMACKS`, thanks.

Other requests: send email to `bastien.molcrette@unibas.ch`