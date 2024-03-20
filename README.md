# VIE_FFT
Fast VIE-FFT method for Modeling Optical Response of Nanostructures

## data folder
data folder: permittivity of materials

result folder: run c/c++ code to generate results into the folder

vie_fft: c/c++ source codes

## c/c++ files
revise Model.h and ReadFile.h (comment/uncomment command lines) to switch the two cases (nanoparticle and solar cell)

## matlab files
Postprocess1_NP and Postprocess2_SolarCell MATLAB files are used to merge and postprocess data in the result folder for nanoparticle and solar cell cases, respectively.

Mie_LPR_Multilayer MATLAB file is used to calculate Mie series for extinction and absorption cross sections of spherical nanoparticle.

Model_Grid MATLAB file is used to generate staircase grids. The grid file can be obtained by uncommenting command lines of GetMaterial() function in Model.h file.

## compiler configuration
The C++ Compiler configuration for 64-bit version of Windows is given by

Additional Include Directories:  Intel Installed Path\MKL\Include; Intel Installed Path\MKL\Include\fftw

Additional Library Directories:  Intel Installed Path\lib\intel64; Intel Installed Path\mkl\em64t\lib

Additional Dependencies:  mkl_core.lib mkl_intel_lp64.lib mkl_sequential.lib
