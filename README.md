# Noise-Performance-Modeling
MATLAB code that plots the power spectral density of the energy jitter, frequency jitter, and phase jitter using dynamical methods. <br>
This package reproduces Fig. 4(a-c) from the below [paper](http://www.photonics.umbc.edu/publications/PdfPapers/PAJ286.pdf): <br>
S. Wang, T. F. Carruthers, and C. R. Menyuk, "Efficiently modeling the noise performance of short-pulse lasers with a computational implementation of dynamical methods," J. Opt. Soc. Am. B 35, 2521-2531 (2018).

<b>Summary:</b> This program calculates the power spectral density (PSD) in the energy, phase, and frequency domain using dynamical methods.
The results are visualized on line graphs after loading the file directory for a short-pulse parameter dataframe or processing user created parameters.

<b>How to: </b> Open FigureGenerator.m, change the values of the parameters in the section "%% Object setup parameters" and run. In the current version, the parameters are set in the module "Parameters_SESAM300_paper2.m"

When generating figures with entirely new parameters, it will take several minutes and generate three figures:

Figure 1: Energy Jitter Power Spectral Density;
Figure 2: Frequency Jitter Power Spectral Density;
Figure 3: Phase Jitter Power Spectral Density;

These figures are identical to Figure 4 of the aforementioned article.

A typo: In the above paper, in Table 3, the value of w_A should be 100 pJ.
