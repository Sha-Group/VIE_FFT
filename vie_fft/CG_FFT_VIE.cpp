// CG_FFT_VIE.cpp : Defines the entry point for the console application.

// http://www.zjuisee.zju.edu.cn/weisha/SourceForge/sourceforge.html  (**please read all the users' remarks carefully at the website**)

// Please keep all the comments when redistributing the program
// BICGSTAB-FFT Volume Integral Equation Program for Nano-Optics Applications (metallic nanopaticles, solar cells, etc)
// Rectangular grids (roof-top basis functions) are adopted (No near-field calculation)
// FFT algorithm relies on Intel MKL (mkl_core.lib mkl_intel_lp64.lib mkl_sequential.lib)
// The program can be naturally parallel (each frequency is independent)
// All the data will be combined and postprocessed by MATLAB program (one case is for nanoparticle, another one is for solar cell.)
// Developer:Wei E.I. Sha
// Email: weisha@zju.edu.cn
// Updated at: Jan. 2019
// please help to cite our works:
// W.E.I. Sha, W.C.H. Choy, Y.P. Chen, and W.C. Chew, "Optical Design of Organic Solar Cell with Hybrid Plasmonic System," Optics Express, vol. 19, no. 17, pp. 15908-15918, 2011.
// W.E.I. Sha, W.C.H. Choy, Y.G. Liu, and W.C. Chew, "Near-Field Multiple Scattering Effects of Plasmonic Nanospheres Embedded into Thin-Film Organic Solar Cells," Applied Physics Letters, vol. 99, no. 11, pp. 113304, 2011.

#include "Model.h"
#include "GreenFunction.h"
#include "CG_Iteration.h"
#include "NearFar.h"
#include "ReadFile.h"

void Initial()
{
	Energy_Flag=0;           //energy flag (no use)
	PPW_x=10;                //points per wavelength along x
	PPW_y=10;                //points per wavelength along y
	PPW_z=10;                //points per wavelength along z
}

int main(int argc, char* argv[])
{
	clock_t start, finish;     //time variables
    start = clock();           //time begins

	ReadParameter();           //read parameters

	//char *str1 = argv[1];
	//int i=atoi(str1);        //get frequency number by command line (no use)

	int b_fre=1;               //frequency starts
	int e_fre=FRE_N+1;         //frequency ends

	for (int i=b_fre; i<e_fre; i++)
	{
		Initial ();  //initialize

		do
		{
			Energy_Flag=1;             //no adaptive grids
			Parameter_Input(i);        //input parameters
			//PreModel();              //premodel (no use)
			GetMaterial();             //material identifier
			GetFlag();                 //current identifier
			Green_Function_Compute();  //get green function
			Incident_Wave();           //get incident field
			FFT_Program();             //FFT main program (Intel MKL)
			FarField1();               //far-field calculation1 (extinction cross section)
			FarField2();               //far-field calculation2 (radar cross section or scattering pattern)
			NearField();               //near-field calculation
		}
		while(Energy_Flag==0);         //energy change is small (no use)
	}

    finish = clock();          //timing ends
	double  duration = (double)(finish - start)/CLOCKS_PER_SEC;  //total time
	printf("Total CPU time: %.1f sec\n",duration);  //show calculation time

	return 0;
}
