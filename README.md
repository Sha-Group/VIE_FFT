# VIE_FFT
A Volume Integral Equation-Fast Fourier Transform Method for Modeling Optical Response of Nanostructures

## Publication
1.  Wei E.I. Sha, Wallace C.H. Choy, Yongpin P. Chen, and Weng Cho Chew, “Optical Design of Organic Solar Cell with Hybrid Plasmonic System,” Optics Express, vol. 19, no. 17, pp. 15908-15918, Aug. 2011.

2.  Chengnian Huang and Wei E.I. Sha, “A Parallel Block Preconditioner-Based VIE-FFT Algorithm for Modeling the Electromagnetic Response From Nanostructures,” IEEE Transactions on Antennas and Propagation, vol. 72, no. 1, pp. 1051-1056, Jan. 2024.


## Data folder
/data folder: permittivity of materials

/result folder: run c/c++ code to generate results into the folder

/vie_fft folder: c/c++ source codes

## C/C++ Files
Revise Model.h and ReadFile.h (comment/uncomment command lines) to switch the two cases (nanoparticle and solar cell)

## Matlab Files
Postprocess1_NP and Postprocess2_SolarCell MATLAB files are used to merge and postprocess data in the result folder for nanoparticle and solar cell cases, respectively.

Mie_LPR_Multilayer MATLAB file is used to calculate Mie series for extinction and absorption cross sections of spherical nanoparticle.

Model_Grid MATLAB file is used to generate staircase grids. The grid file can be obtained by uncommenting command lines of GetMaterial() function in Model.h file.

## Compiler Configuration
The C++ Compiler configuration for 64-bit version of Windows is given by

Additional Include Directories:  Intel Installed Path\MKL\Include; Intel Installed Path\MKL\Include\fftw

Additional Library Directories:  Intel Installed Path\lib\intel64; Intel Installed Path\mkl\em64t\lib

Additional Dependencies:  mkl_core.lib mkl_intel_lp64.lib mkl_sequential.lib
