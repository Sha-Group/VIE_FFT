#ifndef READFILE
#define READFILE

const int FRE_N=41;     //***total frequency points***
const int WBegin=400;   //***smallest wavelength of interest (nm)***
const int WInterval=10; //***wavelength distance (nm)***
const int WScale=100;   //***scaling factor for wavelengths (numerical tricks)***

const char file_m1[256]={"../data/PEDOT_epr.txt"};     //***PEDOT for spacer***
const char file_m2[256]={"../data/Au_epr.txt"};        //***Metal for spheres***
const char file_m3[256]={"../data/P3HT_PCBM_epr.txt"}; //***P3HT:PCBM for active layer***

double Inc_Wavelength[FRE_N];  //incident wavelengths

double m1_r[FRE_N]={0.0};  //real part of permittivity (PEDOT)
double m1_i[FRE_N]={0.0};  //imaginary part of permittivity (PEDOT)
double m2_r[FRE_N]={0.0};  //real part of permittivity (Metal)
double m2_i[FRE_N]={0.0};  //imaginary part of permittivity (Metal)
double m3_r[FRE_N]={0.0};  //real part of permittivity (P3HT:PCBM)
double m3_i[FRE_N]={0.0};  //imaginary part of permittivity (P3HT:PCBM)

//read all the parameters (permittivities and wavelengths)
void ReadParameter()
{
	FILE* file_r1=NULL;  //file1
	FILE* file_r2=NULL;  //file2
	FILE* file_r3=NULL;  //file3

	file_r1=fopen(file_m1,"r");  //PEDOT data
	file_r2=fopen(file_m2,"r");  //Metal data
	file_r3=fopen(file_m3,"r");  //P3HT:PCBM data

	//read parameters by frequency (dispersive media)
	for (int i=0; i<FRE_N; i++)
	{
		fscanf(file_r1,"%lf %lf\n",&m1_r[i],&m1_i[i]);
		fscanf(file_r2,"%lf %lf\n",&m2_r[i],&m2_i[i]);
		fscanf(file_r3,"%lf %lf\n",&m3_r[i],&m3_i[i]);
		Inc_Wavelength[i]=double(WBegin+i*WInterval)/double(WScale);  //incident wavelengths with scaled factors
	}

	fclose(file_r1);
	fclose(file_r2);
	fclose(file_r3);
}

#endif

