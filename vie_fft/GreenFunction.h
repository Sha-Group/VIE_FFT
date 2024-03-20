#ifndef GREENFUNCTION
#define GREENFUNCTION

#include "Model.h"


//green function calculation
void Green_Function_Compute()
{
	int i,j,k,ii,jj,kk;

	//complex
	using namespace std;
	complex<double> V1;
	complex<double> V2;
	complex<double> V3;

	//real
	double V4;
	double V5;
//	double V6;

	double Delta_Volume=delta_x*delta_y*delta_z;  //cell volume
	double Diameter=2*pow(Delta_Volume/(4.0/3.0*PI),1.0/3.0);  //diameter of equivalent sphere
	double Parameter=0.5*k0*Diameter;  //paramter

	//numerical quadrature (no use)
	double vertex[8][3]={{delta_x/2.0,delta_y/2.0,delta_z/2.0},
	                     {delta_x/2.0,delta_y/2.0,-delta_z/2.0},
						 {delta_x/2.0,-delta_y/2.0,delta_z/2.0},
	                     {delta_x/2.0,-delta_y/2.0,-delta_z/2.0},
	                     {-delta_x/2.0,delta_y/2.0,delta_z/2.0},
	                     {-delta_x/2.0,delta_y/2.0,-delta_z/2.0},
	                     {-delta_x/2.0,-delta_y/2.0,delta_z/2.0},
	                     {-delta_x/2.0,-delta_y/2.0,-delta_z/2.0}};  //corner
	
	double cubic_center[6][3]={{delta_x/2.0,0.0,0.0},
	                           {-delta_x/2.0,0.0,0.0},
						       {0.0,delta_y/2.0,0.0},
	                           {0.0,-delta_y/2.0,0.0},
						       {0.0,0.0,delta_z/2.0},
	                           {0.0,0.0,-delta_z/2.0}};  //face center

	double mid_center[6][3]={{delta_x/4.0,0.0,0.0},
	                         {-delta_x/4.0,0.0,0.0},
						     {0.0,delta_y/4.0,0.0},
	                         {0.0,-delta_y/4.0,0.0},
						     {0.0,0.0,delta_z/4.0},
	                         {0.0,0.0,-delta_z/4.0}};  //mid points between face and center point


	//analytical solution for integral
	for (i=0; i<=Nx; i++)
	{
		for (j=0; j<=Ny; j++)
		{
			for (k=0; k<=Nz; k++)
			{
				if ((i+j+k)!=0)  //non singular points
				{
					//exp(-i*k0*norm(cc))*(sinh(-0.5*i*k0*x_eq)/(-0.5*i*k0*x_eq)
					//-cosh(-0.5*i*k0*x_eq))/(1/3*pi*k0^2*x_eq^2*norm(cc))*4/3*pi*(0.5*x_eq)^3

					V4=sqrt(pow(i*delta_x,2.0)+pow(j*delta_y,2.0)+pow(k*delta_z,2.0));
					V1=complex<double>(0.0,-k0*V4);
					V2=complex<double>(0.0,-Parameter);
					V3=exp(V1)*(sinh(V2)/V2-cosh(V2));
					V5=1.0/3.0*PI*k0*k0*pow(Diameter,2.0)*V4;
					Green_Function[i][j][k].real=real(V3)/V5*Delta_Volume;
					Green_Function[i][j][k].imag=imag(V3)/V5*Delta_Volume;

				}
				else //sigularity treatment
				{
					//((1+0.5*i*k0*x_eq)*exp(-0.5*i*k0*x_eq)-1)/(1/6*pi*k0^2*x_eq^3)
					//*4/3*pi*(0.5*x_eq)^3
					V1=complex<double>(1.0,+Parameter);
					V2=complex<double>(0.0,-Parameter);
					V3=complex<double>(1.0,0.0);
					V3=V1*exp(V2)-V3;
					V4=1.0/6.0*PI*k0*k0*pow(Diameter,3.0);
					V5=Delta_Volume;
					Green_Function[i][j][k].real=real(V3)/V4*V5;
					Green_Function[i][j][k].imag=imag(V3)/V4*V5;
				}
			}
		}
	}


	//periodic extension
	for (i=0; i<NNx; i++)
	{
		for (j=0; j<NNy; j++)
		{
			for (k=0; k<NNz; k++)
			{
				if (i>=Nx)  //outer part
				{
					ii=NNx-i;
				}
				else
				{
					ii=i;
				}

				if (j>=Ny)  //outer part
				{
					jj=NNy-j;
				}
				else
				{
					jj=j;
				}

				if (k>=Nz)  //outer part
				{
					kk=NNz-k;
				}
				else
				{
					kk=k;
				}

				Green_Function[i][j][k].real=Green_Function[ii][jj][kk].real;
				Green_Function[i][j][k].imag=Green_Function[ii][jj][kk].imag;
			}
		}
	}

	//initialize fft function
	DFTI_DESCRIPTOR_HANDLE my_desc1_handle;
	long status, l[3];

	l[0] = NNx; l[1] = NNy; l[2]=NNz;  //FFT dimentions
	status = DftiCreateDescriptor(&my_desc1_handle, DFTI_DOUBLE, DFTI_COMPLEX, 3, l);  //create Descriptor
	status = DftiCommitDescriptor(my_desc1_handle);  //commit Descriptor
	status = DftiComputeForward(my_desc1_handle, Green_Function[0][0]);  //forward FFT
	status = DftiFreeDescriptor(&my_desc1_handle);  //FFT end

}

//Incident wave
void Incident_Wave()
{
	//dummy variables
	int i,j,k;

	//complex
	using namespace std;
	complex<double> V1;
	complex<double> V2;
	complex<double> V3;
	complex<double> V4;

	//coordinates
	double xm;
	double ym;
	double zm;

	//wavenumbers
	double k_x=k0*sin(PI*Theta/180.0)*cos(PI*Phi/180.0);
	double k_y=k0*sin(PI*Theta/180.0)*sin(PI*Phi/180.0);
	double k_z=k0*cos(PI*Theta/180.0);

	//Ex_inc
	for (i=0; i<Nx; i++)
	{
		for (j=0; j<Ny; j++)
		{
			for (k=0; k<Nz; k++)
			{
				//analytical expression of Ex_inc
				//exp(-j*(k_x*x_m+k_y*y_m+k_z*z_m))
				//*(exp(-j*k_x*deltax/2)-exp(j*k_x*deltax/2))/(-j*k_x)

				//Ex_inc
				xm=(i+1)*delta_x;
				ym=(j+0.5)*delta_y;
				zm=(k+0.5)*delta_z;

				V1=complex<double>(0.0,-(k_x*xm+k_y*ym+k_z*zm));
				V2=complex<double>(0.0,-(k_x*delta_x/2.0));
				V3=complex<double>(0.0,(k_x*delta_x/2.0));
				V4=complex<double>(0.0,-(k_x));

				if (fabs(k_x)<EPS)
				{
					V4=exp(V1)*delta_x;
				}
				else
				{
					V4=exp(V1)*(exp(V2)-exp(V3))/(V4);
				}

				Ex_inc[i][j][k].real=real(V4)*(cos(PI*Theta/180.0)*cos(PI*Phi/180.0)*\
					                 cos(PI*Psi/180.0)-sin(PI*Phi/180.0)*sin(PI*Psi/180.0));
				Ex_inc[i][j][k].imag=imag(V4)*(cos(PI*Theta/180.0)*cos(PI*Phi/180.0)*\
					                 cos(PI*Psi/180.0)-sin(PI*Phi/180.0)*sin(PI*Psi/180.0));
				//Ey_inc
				xm=(i+0.5)*delta_x;
				ym=(j+1)*delta_y;
				zm=(k+0.5)*delta_z;

				V1=complex<double>(0.0,-(k_x*xm+k_y*ym+k_z*zm));
				V2=complex<double>(0.0,-(k_y*delta_y/2.0));
				V3=complex<double>(0.0,(k_y*delta_y/2.0));
				V4=complex<double>(0.0,-(k_y));

				if (fabs(k_y)<EPS)
				{
					V4=exp(V1)*delta_y;
				}
				else
				{
					V4=exp(V1)*(exp(V2)-exp(V3))/(V4);
				}

				Ey_inc[i][j][k].real=real(V4)*(cos(PI*Theta/180.0)*sin(PI*Phi/180.0)*\
					                 cos(PI*Psi/180.0)+cos(PI*Phi/180.0)*sin(PI*Psi/180.0));
				Ey_inc[i][j][k].imag=imag(V4)*(cos(PI*Theta/180.0)*sin(PI*Phi/180.0)*\
					                 cos(PI*Psi/180.0)+cos(PI*Phi/180.0)*sin(PI*Psi/180.0));

				//Ez_inc
				xm=(i+0.5)*delta_x;
				ym=(j+0.5)*delta_y;
				zm=(k+1)*delta_z;

				V1=complex<double>(0.0,-(k_x*xm+k_y*ym+k_z*zm));
				V2=complex<double>(0.0,-(k_z*delta_z/2.0));
				V3=complex<double>(0.0,(k_z*delta_z/2.0));
				V4=complex<double>(0.0,-(k_z));

				if (fabs(k_z)<EPS)
				{
					V4=exp(V1)*delta_z;
				}
				else
				{
					V4=exp(V1)*(exp(V2)-exp(V3))/(V4);
				}

				Ez_inc[i][j][k].real=-real(V4)*sin(PI*Theta/180.0)*cos(PI*Psi/180.0);
				Ez_inc[i][j][k].imag=-imag(V4)*sin(PI*Theta/180.0)*cos(PI*Psi/180.0);
			}//k
		}//j
	}//i
}
#endif


////高斯求和
//double sum_r=0.0;
//double sum_i=0.0;

////顶点
//for (ll=0; ll<8; ll++)
//{
//	V4=sqrt(pow(i*delta_x+vertex[ll][0],2.0)+
//		    pow(j*delta_y+vertex[ll][1],2.0)+
//			pow(k*delta_z+vertex[ll][2],2.0));
//	sum_r=sum_r+5*cos(k0*V4)/(4*PI*V4);
//	sum_i=sum_i-5*sin(k0*V4)/(4*PI*V4);
//}

////六个面中点
//for (ll=0; ll<6; ll++)
//{
//	V4=sqrt(pow(i*delta_x+cubic_center[ll][0],2.0)+
//		    pow(j*delta_y+cubic_center[ll][1],2.0)+
//			pow(k*delta_z+cubic_center[ll][2],2.0));

//	sum_r=sum_r+8*cos(k0*V4)/(4*PI*V4);
//	sum_i=sum_i-8*sin(k0*V4)/(4*PI*V4);
//}

////六个面和中心的中点
//for (ll=0; ll<6; ll++)
//{
//	V4=sqrt(pow(i*delta_x+mid_center[ll][0],2.0)+
//		    pow(j*delta_y+mid_center[ll][1],2.0)+
//			pow(k*delta_z+mid_center[ll][2],2.0));

//	sum_r=sum_r+128*cos(k0*V4)/(4*PI*V4);
//	sum_i=sum_i-128*sin(k0*V4)/(4*PI*V4);
//}

////立方体中心
//V4=sqrt(pow(i*delta_x,2.0)+
//		pow(j*delta_y,2.0)+
//	    pow(k*delta_z,2.0));
//sum_r=sum_r+(-496)*cos(k0*V4)/(4*PI*V4);
//sum_i=sum_i-(-496)*sin(k0*V4)/(4*PI*V4);

////高斯求和
//sum_r=sum_r*Delta_Volume/360;
//sum_i=sum_i*Delta_Volume/360;

//Green_Function[i][j][k].real=sum_r;
//Green_Function[i][j][k].imag=sum_i;