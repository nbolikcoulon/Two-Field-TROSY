# Two-Field-TROSY
This repository contains scripts to calculate two-dimensional NMR spectra for a two-field TROSY type of pulse sequence, applied on 13C-19F spin pairs. For more information on the pulse sequence, see: Nicolas Bolik-Coulon, Philippe Pelupessy, Guillaume Bouvignies, Fabien Ferrage, J. Magn. Reson. Open, 2020, 100007

The scripts use python 3, numpy, spicy and the nmrglue packages. Launch using the command line: python simulatesequence.py

The following architectures must be created in the folder containing the scripts to save the simulated spectra:
- Spectra/25ns/3FTyr/Noise
- Spectra/100ns/3FTyr/Noise
- Spectra/25ns/4FPhe/Noise
- Spectra/100ns/4FPhe/Noise

The following architectures must be created in the folder containing the scripts to save the figures for the cross-sections:
- Fig/25ns/3FTyr/13C
- Fig/25ns/3FTyr/19F
- Fig/100ns/3FTyr/13C
- Fig/100ns/3FTyr/19F
- Fig/25ns/4FPhe/13C
- Fig/25ns/4FPhe/19F
- Fig/100ns/4FPhe/13C
- Fig/100ns/4FPhe/19F
