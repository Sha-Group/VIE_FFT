#ifndef FFT_SW
#define FFT_SW

#include "Model.h"
//#include ".\CG_Iteration.h"

//Ref. IEEE trans. AP. Vol. 37, No. 5, May 1989

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

//reset current Jx Jy Jz
void Current_Reset(MKL_Complex16*** Current, const int Flag, const int* Dim)
{
	//index
	int index[3];
	int i,j,k;

	//loop
	for (i=0; i<Dim[0]; i++)
	{
		for (j=0; j<Dim[1]; j++)
		{
			for (k=0; k<Dim[2]; k++)
			{
				index[0]=int(Flag_X[i][j][k]);
				index[1]=int(Flag_Y[i][j][k]);
				index[2]=int(Flag_Z[i][j][k]);

				Current[i][j][k].real=Current[i][j][k].real*double(index[Flag]);  //reset current (real part)
				Current[i][j][k].imag=Current[i][j][k].imag*double(index[Flag]);  //reset current (imaginary part)
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

//Self interaction
void L_DC_Dq(MKL_Complex16*** temp_xyz, MKL_Complex16*** current_xyz,\
		     MKL_Complex16*** const Green, const double delta_xyz, const int Index_Source, const int* Dim)
{

	int i,j,k;

	//complex
	using namespace std;
	complex<double> V1;
	complex<double> V2;
	complex<double> V3;
	complex<double> V4;
	complex<double> V5;
	complex<double> V6;

	int index1[3];   //dimension
	int index2[3];   //dimension

	//index
	index1[0]=2*Dim[0];
	index1[1]=2*Dim[1];
	index1[2]=2*Dim[2];

	//FFT Initial
	DFTI_DESCRIPTOR_HANDLE my_desc1_handle;
	long status, l[3];

	l[0] = index1[0]; l[1] = index1[1]; l[2]=index1[2];  //FFT dimention
	status = DftiCreateDescriptor(&my_desc1_handle, DFTI_DOUBLE, DFTI_COMPLEX, 3, l);  //create descriptor
	status = DftiCommitDescriptor(my_desc1_handle);  //commit descriptor

	//FFT forward
	status = DftiComputeForward( my_desc1_handle, current_xyz[0][0]);

	//normalized
	V5=complex<double>(0,k0*Z0*delta_xyz/index1[0]/index1[1]/index1[2]);
	V6=complex<double>(0,Z0/k0/delta_xyz/index1[0]/index1[1]/index1[2]);

	//multiplication in frequency domain
	for (i=0; i<index1[0]; i++)
	{
		for (j=0; j<index1[1]; j++)
		{
			for (k=0; k<index1[2]; k++)
			{
				//index
				index2[0]=i;
				index2[1]=j;
				index2[2]=k;
				//current in frequency domain
				V1=complex<double>(current_xyz[i][j][k].real,current_xyz[i][j][k].imag);
				//green function in frequency domain
				V2=complex<double>(Green[i][j][k].real,Green[i][j][k].imag);
				//source
				V3=complex<double>(1.0-cos(2*PI*index2[Index_Source]/double(index1[Index_Source])),\
					               sin(2*PI*index2[Index_Source]/double(index1[Index_Source])));
				//field
				V4=complex<double>(cos(2*PI*index2[Index_Source]/double(index1[Index_Source]))-1.0,\
					               sin(2*PI*index2[Index_Source]/double(index1[Index_Source])));
				V1=V1*V2*(V5+V6*V3*V4);
				current_xyz[i][j][k].real=real(V1);
				current_xyz[i][j][k].imag=imag(V1);
			}
		}
	}

	//inverse FFT
	status = DftiComputeBackward(my_desc1_handle, current_xyz[0][0]); 

//	Current_Reset(current_xyz,Index_Source);

	//summation
	for (i=0; i<Dim[0]; i++)
	{
		for (j=0; j<Dim[1]; j++)
		{
			for (k=0; k<Dim[2]; k++)
			{
				temp_xyz[i][j][k].real=temp_xyz[i][j][k].real+current_xyz[i][j][k].real;
				temp_xyz[i][j][k].imag=temp_xyz[i][j][k].imag+current_xyz[i][j][k].imag;
			}
		}
	}

	//FFT end
	status = DftiFreeDescriptor(&my_desc1_handle);
}


//Self interaction (conjugate)
void L_DC_Dq_Conjugate(MKL_Complex16*** temp_xyz, MKL_Complex16*** current_xyz,\
		               MKL_Complex16*** const Green, const double delta_xyz, const int Index_Source, const int* Dim)
{
	int i,j,k;

	//complex
	using namespace std;
	complex<double> V1;
	complex<double> V2;
	complex<double> V3;
	complex<double> V4;
	complex<double> V5;
	complex<double> V6;

	int index1[3];   //dimension
	int index2[3];   //dimension

	//dimension
	index1[0]=2*Dim[0];
	index1[1]=2*Dim[1];
	index1[2]=2*Dim[2];

	//FFT Initial
	DFTI_DESCRIPTOR_HANDLE my_desc1_handle;
	long status, l[3];

	l[0] = index1[0]; l[1] = index1[1]; l[2]=index1[2];  //FFT dimention
	status = DftiCreateDescriptor(&my_desc1_handle, DFTI_DOUBLE, DFTI_COMPLEX, 3, l);
	status = DftiCommitDescriptor(my_desc1_handle);


	//forward fft
	status = DftiComputeForward( my_desc1_handle, current_xyz[0][0]);


	//normalized
	V5=complex<double>(0,-k0*Z0*delta_xyz/index1[0]/index1[1]/index1[2]);
	V6=complex<double>(0,-Z0/k0/delta_xyz/index1[0]/index1[1]/index1[2]);

	//multiplication in frequency domain
	for (i=0; i<index1[0]; i++)
	{
		for (j=0; j<index1[1]; j++)
		{
			for (k=0; k<index1[2]; k++)
			{
				//loop index
				index2[0]=i;
				index2[1]=j;
				index2[2]=k;
				//current in frequency domain
				V1=complex<double>(current_xyz[i][j][k].real,current_xyz[i][j][k].imag);
				//green function in frequency domain
				V2=complex<double>(Green[i][j][k].real,-Green[i][j][k].imag);
				//field
				V3=complex<double>(1.0-cos(2*PI*index2[Index_Source]/double(index1[Index_Source])),\
					               -sin(2*PI*index2[Index_Source]/double(index1[Index_Source])));
				//source
				V4=complex<double>(cos(2*PI*index2[Index_Source]/double(index1[Index_Source]))-1.0,\
					              -sin(2*PI*index2[Index_Source]/double(index1[Index_Source])));
				V1=V1*V2*(V5+V6*V3*V4);

				current_xyz[i][j][k].real=real(V1);
				current_xyz[i][j][k].imag=imag(V1);
			}
		}
	}


	//inverse FFT
	status = DftiComputeBackward(my_desc1_handle, current_xyz[0][0]); 
//	Current_Reset(current_xyz,Index_Source);

	//summation
	for (i=0; i<Dim[0]; i++)
	{
		for (j=0; j<Dim[1]; j++)
		{
			for (k=0; k<Dim[2]; k++)
			{
				temp_xyz[i][j][k].real=temp_xyz[i][j][k].real+current_xyz[i][j][k].real;
				temp_xyz[i][j][k].imag=temp_xyz[i][j][k].imag+current_xyz[i][j][k].imag;
			}
		}
	}

	//FFT ends
	status = DftiFreeDescriptor(&my_desc1_handle);
}
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

//FFT for charge
void L_Dq(MKL_Complex16*** temp_xyz, MKL_Complex16*** current_xyz,\
		  MKL_Complex16*** const Green, const int Index_Source,\
		  const int Index_Field, const double delta_xyz_source, const int* Dim)
{
	int i,j,k;
	int index1[3];   //dimension
	int index2[3];   //dimension

	//dimension index
	index1[0]=2*Dim[0];
	index1[1]=2*Dim[1];
	index1[2]=2*Dim[2];

	//complex var
	using namespace std;
	complex<double> V1;
	complex<double> V2;
	complex<double> V3;
	complex<double> V4;
	complex<double> V5;

	//FFT initial
	DFTI_DESCRIPTOR_HANDLE my_desc1_handle;
	long status, l[3];

	l[0] = index1[0]; l[1] = index1[1]; l[2]=index1[2];  //FFT dimention
	status = DftiCreateDescriptor(&my_desc1_handle, DFTI_DOUBLE, DFTI_COMPLEX, 3, l);
	status = DftiCommitDescriptor(my_desc1_handle);

	//forward FFT
	status = DftiComputeForward( my_desc1_handle, current_xyz[0][0]);
	//normalized
	V5=complex<double>(0,Z0/k0/delta_xyz_source/index1[0]/index1[1]/index1[2]);

	//multiplication in frequency domain
	for (i=0; i<index1[0]; i++)
	{
		for (j=0; j<index1[1]; j++)
		{
			for (k=0; k<index1[2]; k++)
			{
				//loop index
				index2[0]=i;
				index2[1]=j;
				index2[2]=k;

				//(exp(j*2*pi*(s_n-1)/(2*N))-1)*(1-exp(-j*2*pi*(s_m-1)/(2*M)));  %  M source£¬N field
				//current
				V1=complex<double>(current_xyz[i][j][k].real,current_xyz[i][j][k].imag);
				//green function
				V2=complex<double>(Green[i][j][k].real,Green[i][j][k].imag);
				//source
				V3=complex<double>(1.0-cos(2*PI*index2[Index_Source]/double(index1[Index_Source])),\
					               sin(2*PI*index2[Index_Source]/double(index1[Index_Source])));
				//field
				V4=complex<double>(cos(2*PI*index2[Index_Field]/double(index1[Index_Field]))-1.0,\
					               sin(2*PI*index2[Index_Field]/double(index1[Index_Field])));
				V1=V5*V1*V2*V3*V4;
				current_xyz[i][j][k].real=real(V1);
				current_xyz[i][j][k].imag=imag(V1);
			}
		}
	}

	//inverse FFT
	status = DftiComputeBackward(my_desc1_handle, current_xyz[0][0]);  

//	Current_Reset(current_xyz,Index_Source);
//	Current_Reset(current_xyz,Index_Field);

	//summation
	for (i=0; i<Dim[0]; i++)
	{
		for (j=0; j<Dim[1]; j++)
		{
			for (k=0; k<Dim[2]; k++)
			{
				temp_xyz[i][j][k].real=temp_xyz[i][j][k].real+current_xyz[i][j][k].real;
				temp_xyz[i][j][k].imag=temp_xyz[i][j][k].imag+current_xyz[i][j][k].imag;
			}
		}
	}

	//FFT ends
	status = DftiFreeDescriptor(&my_desc1_handle);
}

//FFT for charge (conjugate)
void L_Dq_Conjugate(MKL_Complex16*** temp_xyz, MKL_Complex16*** current_xyz,\
		            MKL_Complex16*** const Green, const int Index_Source,\
					const int Index_Field, const double delta_xyz_field, const int* Dim)
{
	int i,j,k;
	int index1[3];   //dimension
	int index2[3];   //dimention

	//dimension index
	index1[0]=2*Dim[0];
	index1[1]=2*Dim[1];
	index1[2]=2*Dim[2];

	//complex
	using namespace std;
	complex<double> V1;
	complex<double> V2;
	complex<double> V3;
	complex<double> V4;
	complex<double> V5;

	//FFT initial
	DFTI_DESCRIPTOR_HANDLE my_desc1_handle;
	long status, l[3];

	l[0] = index1[0]; l[1] = index1[1]; l[2]=index1[2];  //FFT dimension
	status = DftiCreateDescriptor(&my_desc1_handle, DFTI_DOUBLE, DFTI_COMPLEX, 3, l);
	status = DftiCommitDescriptor(my_desc1_handle);

	//forward FFT fir current
	status = DftiComputeForward( my_desc1_handle, current_xyz[0][0]);
	//normalized
	V5=complex<double>(0,-Z0/k0/delta_xyz_field/index1[0]/index1[1]/index1[2]);


	//multiplication in frequency domain
	for (i=0; i<index1[0]; i++)
	{
		for (j=0; j<index1[1]; j++)
		{
			for (k=0; k<index1[2]; k++)
			{
				//loop index
				index2[0]=i;
				index2[1]=j;
				index2[2]=k;

				//(exp(-j*2*pi*(s_n-1)/(2*N))-1)*(1-exp(-j*2*pi*(s_m-1)/(2*M)));  %  M source£¬N field
				//current in frequency domain
				V1=complex<double>(current_xyz[i][j][k].real,current_xyz[i][j][k].imag);
				//green function in frequency domain
				V2=complex<double>(Green[i][j][k].real,-Green[i][j][k].imag);
				//field
				V3=complex<double>(1.0-cos(2*PI*index2[Index_Field]/double(index1[Index_Field])),\
					               -sin(2*PI*index2[Index_Field]/double(index1[Index_Field])));
				//source
				V4=complex<double>(cos(2*PI*index2[Index_Source]/double(index1[Index_Source]))-1.0,\
					              -sin(2*PI*index2[Index_Source]/double(index1[Index_Source])));
				V1=V5*V1*V2*V3*V4;
				current_xyz[i][j][k].real=real(V1);
				current_xyz[i][j][k].imag=imag(V1);
			}
		}
	}

	//inverse FFT
	status = DftiComputeBackward(my_desc1_handle, current_xyz[0][0]); 
//	Current_Reset(current_xyz,Index_Source);
//	Current_Reset(current_xyz,Index_Field);

	//summation
	for (i=0; i<Dim[0]; i++)
	{
		for (j=0; j<Dim[1]; j++)
		{
			for (k=0; k<Dim[2]; k++)
			{
				temp_xyz[i][j][k].real=temp_xyz[i][j][k].real+current_xyz[i][j][k].real;
				temp_xyz[i][j][k].imag=temp_xyz[i][j][k].imag+current_xyz[i][j][k].imag;
			}
		}
	}

	//FFT ends
	status = DftiFreeDescriptor(&my_desc1_handle);
}



/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

//Jx (non-FFT terms)
void LXX_DO(MKL_Complex16*** temp_x, MKL_Complex16*** current_x, const int* Dim)
{
	//loop var
	int i,j,k;

	//complex var
	using namespace std;
	complex<double> V1;
	complex<double> V2;
	complex<double> V3;
	complex<double> V4;
	complex<double> V5;

	//coef
	V4=complex<double>(0.0,-delta_x/2.0*Z0/k0);
	V5=complex<double>(1.0,0.0);

	//summation
	for (i=0; i<Dim[0]-1; i++)
	{
		for (j=0; j<Dim[1]; j++)
		{
			for (k=0; k<Dim[2]; k++)
			{
				V1=complex<double>(current_x[i][j][k].real,current_x[i][j][k].imag);
				V2=complex<double>(Real_EPR[Material[i+1][j][k]]-1.0,Imag_EPR[Material[i+1][j][k]]);
				V3=complex<double>(Real_EPR[Material[i][j][k]]-1.0,Imag_EPR[Material[i][j][k]]);
				V1=V4*(V5/V2+V5/V3)*V1;
				temp_x[i][j][k].real=temp_x[i][j][k].real+real(V1);
				temp_x[i][j][k].imag=temp_x[i][j][k].imag+imag(V1);
			}
		}
	}
}

//Jx (non-FFT terms, conjugate)
void LXX_DO_Conjugate(MKL_Complex16*** temp_x, MKL_Complex16*** current_x, const int* Dim)
{
	//loop var
	int i,j,k;

	//complex var
	using namespace std;
	complex<double> V1;
	complex<double> V2;
	complex<double> V3;
	complex<double> V4;
	complex<double> V5;

	//coef
	V4=complex<double>(0,delta_x/2.0*Z0/k0);
	V5=complex<double>(1.0,0.0);

	//summation
	for (i=0; i<Dim[0]-1; i++)
	{
		for (j=0; j<Dim[1]; j++)
		{
			for (k=0; k<Dim[2]; k++)
			{
				V1=complex<double>(current_x[i][j][k].real,current_x[i][j][k].imag);
				V2=complex<double>(Real_EPR[Material[i+1][j][k]]-1.0,-Imag_EPR[Material[i+1][j][k]]);
				V3=complex<double>(Real_EPR[Material[i][j][k]]-1.0,-Imag_EPR[Material[i][j][k]]);
				V1=V4*(V5/V2+V5/V3)*V1;
				temp_x[i][j][k].real=temp_x[i][j][k].real+real(V1);
				temp_x[i][j][k].imag=temp_x[i][j][k].imag+imag(V1);
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

//Jy (non-FFT terms)
void LYY_DO(MKL_Complex16*** temp_y, MKL_Complex16*** current_y, const int* Dim)
{
	//loop var
	int i,j,k;

	//complex var
	using namespace std;
	complex<double> V1;
	complex<double> V2;
	complex<double> V3;
	complex<double> V4;
	complex<double> V5;

	//coef
	V4=complex<double>(0.0,-delta_y/2.0*Z0/k0);
	V5=complex<double>(1.0,0.0);

	//summation
	for (i=0; i<Dim[0]; i++)
	{
		for (j=0; j<Dim[1]-1; j++)
		{
			for (k=0; k<Dim[2]; k++)
			{
				V1=complex<double>(current_y[i][j][k].real,current_y[i][j][k].imag);
				V2=complex<double>(Real_EPR[Material[i][j+1][k]]-1.0,Imag_EPR[Material[i][j+1][k]]);
				V3=complex<double>(Real_EPR[Material[i][j][k]]-1.0,Imag_EPR[Material[i][j][k]]);
				V1=V4*(V5/V2+V5/V3)*V1;
				temp_y[i][j][k].real=temp_y[i][j][k].real+real(V1);
				temp_y[i][j][k].imag=temp_y[i][j][k].imag+imag(V1);
			}
		}
	}
}

//Jy (non-FFT terms, conjugate)
void LYY_DO_Conjugate(MKL_Complex16*** temp_y, MKL_Complex16*** current_y, const int* Dim)
{
	//loop var
	int i,j,k;

	//complex var
	using namespace std;
	complex<double> V1;
	complex<double> V2;
	complex<double> V3;
	complex<double> V4;
	complex<double> V5;

	//coef
	V4=complex<double>(0,delta_y/2.0*Z0/k0);
	V5=complex<double>(1.0,0.0);

	//summation
	for (i=0; i<Dim[0]; i++)
	{
		for (j=0; j<Dim[1]-1; j++)
		{
			for (k=0; k<Dim[2]; k++)
			{
				V1=complex<double>(current_y[i][j][k].real,current_y[i][j][k].imag);
				V2=complex<double>(Real_EPR[Material[i][j+1][k]]-1.0,-Imag_EPR[Material[i][j+1][k]]);
				V3=complex<double>(Real_EPR[Material[i][j][k]]-1.0,-Imag_EPR[Material[i][j][k]]);
				V1=V4*(V5/V2+V5/V3)*V1;
				temp_y[i][j][k].real=temp_y[i][j][k].real+real(V1);
				temp_y[i][j][k].imag=temp_y[i][j][k].imag+imag(V1);
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

//Jz (non-FFT terms)
void LZZ_DO(MKL_Complex16*** temp_z, MKL_Complex16*** current_z, const int* Dim)
{
	//loop var
	int i,j,k;

	//complex var
	using namespace std;
	complex<double> V1;
	complex<double> V2;
	complex<double> V3;
	complex<double> V4;
	complex<double> V5;

	//coef
	V4=complex<double>(0.0,-delta_z/2.0*Z0/k0);
	V5=complex<double>(1.0,0.0);

	//summation
	for (i=0; i<Dim[0]; i++)
	{
		for (j=0; j<Dim[1]; j++)
		{
			for (k=0; k<Dim[2]-1; k++)
			{
				V1=complex<double>(current_z[i][j][k].real,current_z[i][j][k].imag);
				V2=complex<double>(Real_EPR[Material[i][j][k+1]]-1.0,Imag_EPR[Material[i][j][k+1]]);
				V3=complex<double>(Real_EPR[Material[i][j][k]]-1.0,Imag_EPR[Material[i][j][k]]);
				V1=V4*(V5/V2+V5/V3)*V1;
				temp_z[i][j][k].real=temp_z[i][j][k].real+real(V1);
				temp_z[i][j][k].imag=temp_z[i][j][k].imag+imag(V1);
			}
		}
	}
}

//Jz (non-FFT terms, conjugate)
void LZZ_DO_Conjugate(MKL_Complex16*** temp_z, MKL_Complex16*** current_z, const int* Dim)
{
	//loop var
	int i,j,k;

	//complex var
	using namespace std;
	complex<double> V1;
	complex<double> V2;
	complex<double> V3;
	complex<double> V4;
	complex<double> V5;

	//coef
	V4=complex<double>(0,delta_z/2.0*Z0/k0);
	V5=complex<double>(1.0,0.0);

	//summation
	for (i=0; i<Dim[0]; i++)
	{
		for (j=0; j<Dim[1]; j++)
		{
			for (k=0; k<Dim[2]-1; k++)
			{
				V1=complex<double>(current_z[i][j][k].real,current_z[i][j][k].imag);
				V2=complex<double>(Real_EPR[Material[i][j][k+1]]-1.0,-Imag_EPR[Material[i][j][k+1]]);
				V3=complex<double>(Real_EPR[Material[i][j][k]]-1.0,-Imag_EPR[Material[i][j][k]]);
				V1=V4*(V5/V2+V5/V3)*V1;
				temp_z[i][j][k].real=temp_z[i][j][k].real+real(V1);
				temp_z[i][j][k].imag=temp_z[i][j][k].imag+imag(V1);
			}
		}
	}
}

#endif
