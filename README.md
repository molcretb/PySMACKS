# PySMACKS
Python version of the SMACKS framework for smFRET analysis

Modules that are ready for test:

* Traces_extractor_pipeline.ipynb : JP notebook used to run the traces extraction from the raw TIFF movies, including chromatic aberration and drift corrections
* PySMACKS_FRET_corrections.ipnyb : JP notebook used to process the FRET corrections $\alpha$ - $\delta$ - $\beta$ - $\gamma$, and plot the corrected FRET efficiency histogram
* GUI_traces_viewer.py : A GUI traces viewer to manually check the extracted traces, filter them using SNR threshold and view the corresponding spot from the raw TIFF movies; will serve as a basis for the quality checks
* GUI_traces_annotater.py : similar to GUI_traces_viewer.py, with additional function to select sub-regions of individual traces for downstream analysis, as well as some user-friendly feature, such as zooming with the mouse roll
