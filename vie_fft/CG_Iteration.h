#ifndef CG_ITERATION
#define CG_ITERATION


#include "fft.h"
//CG or BICG in Prof. Jianming Jin's FEM book
//BICG-Stab in http://www.netlib.org/templates/matlab/
//CORS by Y.F. Jing in PIER 99, 427-451, 2009

//iterative variables
MKL_Complex16*** Temp_Current;  //expanded current (zero padding)

MKL_Complex16*** Apx;			//impedance matrix1(x component)
MKL_Complex16*** Apy;			//impedance matrix1(y component)
MKL_Complex16*** Apz;			//impedance matrix1(z component)

MKL_Complex16*** Aqx;			//impedance matrix2(x component)
MKL_Complex16*** Aqy;			//impedance matrix2(y component)
MKL_Complex16*** Aqz;			//impedance matrix2(z component)

MKL_Complex16*** q1x;			//modified gradient direction(x)
MKL_Complex16*** q1y;			//modified gradient direction(y)
MKL_Complex16*** q1z;			//modified gradient direction(z)

MKL_Complex16*** p1x;			//gradient direction(x)
MKL_Complex16*** p1y;			//gradient direction(y)
MKL_Complex16*** p1z;			//gradient direction(z)

MKL_Complex16*** r0x;           //residual(x)
MKL_Complex16*** r0y;           //residual(y)
MKL_Complex16*** r0z;           //residual(z)

MKL_Complex16*** q1x_c;			//modified gradient direction(x conj)
MKL_Complex16*** q1y_c;			//modified gradient direction(y conj)
MKL_Complex16*** q1z_c;			//modified gradient direction(z conj)

MKL_Complex16*** p1x_c;			//gradient direction(x conj)
MKL_Complex16*** p1y_c;			//gradient direction(y conj)
MKL_Complex16*** p1z_c;			//gradient direction(z conj)

MKL_Complex16*** r0x_c;         //residual1(x conj)
MKL_Complex16*** r0y_c;         //residual1(y conj)
MKL_Complex16*** r0z_c;         //residual1(z conj)

unsigned char*** Pre_Flag;      //preconditioner flag (no use)
const double CG_THR=0.001;      //Threshold  ******
const int Iteration_Type=3;     //iterative solvers（1:CG;2:BICG;3:BICGSTAB;4:CORS;5：CORSTAB）

//expand currents
void Current_Expansion(MKL_Complex16*** Temp_Current_V, MKL_Complex16*** const Current, int Flag, const int* Dim)
{
    //dummy variables
	int i,j,k;
	//total data size
	int Total=8*Dim[0]*Dim[1]*Dim[2];
	//index
	int index[3][3];
	//Jx dimension
	index[0][0]=Dim[0]-1;index[0][1]=Dim[1];index[0][2]=Dim[2];
	//Jy dimension
	index[1][0]=Dim[0];index[1][1]=Dim[1]-1;index[1][2]=Dim[2];
	//Jz dimension
	index[2][0]=Dim[0];index[2][1]=Dim[1];index[2][2]=Dim[2]-1;

	//set zero
	memset(Temp_Current_V[0][0],0,Total*sizeof(MKL_Complex16));

	//loop
	for (i=0; i<Dim[0]; i++)
	{
		for (j=0; j<Dim[1]; j++)
		{
			for (k=0; k<Dim[2]; k++)
			{
				if (i<index[Flag][0] && j<index[Flag][1] && k<index[Flag][2])
				{
					Temp_Current_V[i][j][k].real=Current[i][j][k].real;
					Temp_Current_V[i][j][k].imag=Current[i][j][k].imag;
				}
			}
		}
	}
}

//conjugate x component by FFT(L_11*, L_12*, L_13*)
void X_FFT_Conjugate(MKL_Complex16*** Field, MKL_Complex16*** const Sourx, MKL_Complex16*** const Soury,\
		             MKL_Complex16*** const Sourz, MKL_Complex16*** temp, MKL_Complex16*** const Green, const int* Dim)
{
	Current_Expansion(temp, Sourx, 0, Dim);                     //expand Jx
	LXX_DO_Conjugate(Field, temp, Dim);                         //Jx (non-FFT terms, conjugate)
	L_DC_Dq_Conjugate(Field, temp, Green, delta_x, 0, Dim);     //Jx (FFT terms, conjugate)

	Current_Expansion(temp, Soury, 1, Dim);                     //expand Jy
	L_Dq_Conjugate(Field, temp, Green, 1, 0, delta_x, Dim);     //Jy (FFT terms, conjugate)

	Current_Expansion(temp, Sourz, 2, Dim);                     //expand Jz
	L_Dq_Conjugate(Field, temp, Green, 2, 0, delta_x, Dim);     //Jz (FFT terms, conjugate))
}

//conjugate y component by FFT(L_21*, L_22*, L_23*)
void Y_FFT_Conjugate(MKL_Complex16*** Field, MKL_Complex16*** const Sourx, MKL_Complex16*** const Soury,\
		             MKL_Complex16*** const Sourz, MKL_Complex16*** temp, MKL_Complex16*** const Green, const int* Dim)
{
	Current_Expansion(temp, Soury, 1, Dim);                     //expand Jy
	LYY_DO_Conjugate(Field, temp, Dim);                         //Jy (non-FFT terms, conjugate)
	L_DC_Dq_Conjugate(Field, temp, Green, delta_y, 1, Dim);     //Jy (FFT terms, conjugate)

	Current_Expansion(temp, Sourx, 0, Dim);                     //expand Jx
	L_Dq_Conjugate(Field, temp, Green, 0, 1, delta_y, Dim);     //Jx (FFT terms, conjugate)

	Current_Expansion(temp, Sourz, 2, Dim);                     //expand Jz
	L_Dq_Conjugate(Field, temp, Green, 2, 1, delta_y, Dim);     //Jz (FFT terms, conjugate)
}

//conjugate z component by FFT(L_31*, L_32*, L_33*)
void Z_FFT_Conjugate(MKL_Complex16*** Field,MKL_Complex16*** const Sourx, MKL_Complex16*** const Soury,\
		             MKL_Complex16*** const Sourz, MKL_Complex16*** temp, MKL_Complex16*** const Green, const int* Dim)
{
	Current_Expansion(temp, Sourz, 2, Dim);                     //expand Jz
	LZZ_DO_Conjugate(Field, temp, Dim);                         //Jz (non-FFT terms, conjugate)
	L_DC_Dq_Conjugate(Field, temp, Green, delta_z, 2, Dim);     //Jz (FFT terms, conjugate)

	Current_Expansion(temp, Sourx, 0, Dim);                     //expand Jx
	L_Dq_Conjugate(Field, temp, Green, 0, 2, delta_z, Dim);     //Jx (FFT terms, conjugate)

	Current_Expansion(temp, Soury, 1, Dim);                     //expand Jy
	L_Dq_Conjugate(Field, temp, Green, 1, 2, delta_z, Dim);     //Jy (FFT terms, conjugate)
}

//x component by FFT (L_11, L_12, L_13)
void X_FFT(MKL_Complex16*** Field, MKL_Complex16*** const Sourx, MKL_Complex16*** const Soury,\
		   MKL_Complex16*** const Sourz, MKL_Complex16*** temp, MKL_Complex16*** const Green, const int* Dim)
{
	Current_Expansion(temp, Sourx, 0, Dim);                //expand Jx
	LXX_DO(Field, temp, Dim);                              //Jx (non-FFT terms)
	L_DC_Dq(Field, temp, Green, delta_x, 0, Dim);          //Jx (FFT terms)

	Current_Expansion(temp, Soury, 1, Dim);                //expand Jy
	L_Dq(Field, temp, Green, 1, 0, delta_y, Dim);          //Jy (FFT terms)

	Current_Expansion(temp, Sourz, 2, Dim);                //expand Jz
	L_Dq(Field, temp, Green, 2, 0, delta_z, Dim);          //Jz (FFT terms)
}

//y component by FFT(L_21, L_22, L_23)
void Y_FFT(MKL_Complex16*** Field, MKL_Complex16*** const Sourx, MKL_Complex16*** const Soury,\
		   MKL_Complex16*** const Sourz, MKL_Complex16*** temp, MKL_Complex16*** const Green, const int* Dim)
{
	Current_Expansion(temp, Soury, 1, Dim);                //expand Jy
	LYY_DO(Field, temp, Dim);                              //Jy (non-FFT terms)
	L_DC_Dq(Field, temp, Green, delta_y, 1, Dim);          //Jy (FFT terms)

	Current_Expansion(temp, Sourx, 0, Dim);                //expand Jx
	L_Dq(Field, temp, Green, 0, 1, delta_x, Dim);          //Jx (FFT terms)

	Current_Expansion(temp, Sourz, 2, Dim);                //expand Jz
	L_Dq(Field, temp, Green, 2, 1, delta_z, Dim);          //Jz (FFT terms)
}

//z component by FFT(L_31, L_32, L_33)
void Z_FFT(MKL_Complex16*** Field, MKL_Complex16*** const Sourx, MKL_Complex16*** const Soury,\
		   MKL_Complex16*** const Sourz, MKL_Complex16*** temp, MKL_Complex16*** const Green, const int* Dim)
{
	Current_Expansion(temp, Sourz, 2, Dim);                //expand Jz
	LZZ_DO(Field, temp, Dim);                              //Jz (non-FFT terms)
	L_DC_Dq(Field, temp, Green, delta_z, 2, Dim);          //Jz (FFT terms)

	Current_Expansion(temp, Sourx, 0, Dim);                //expand Jx
	L_Dq(Field, temp, Green, 0, 2, delta_x, Dim);          //Jx (FFT terms)

	Current_Expansion(temp, Soury, 1, Dim);                //expand Jy
	L_Dq(Field, temp, Green, 1, 2, delta_y, Dim);          //Jy (FFT terms)
}

//conjugate gradient-FFT
void CG_FFT()
{
	int Total_inc=Nx*Ny*Nz;   //total size
    int Dim[3]={Nx,Ny,Nz};    //dimension

	//allocate memory Temp_Current
	Temp_Current=Memory3D(Temp_Current,NNx,NNy,NNz);
	//allocate memory r0
	r0x=Memory3D(r0x,Nx,Ny,Nz);
	r0y=Memory3D(r0y,Nx,Ny,Nz);
	r0z=Memory3D(r0z,Nx,Ny,Nz);
	//allocate memory q1
	q1x=Memory3D(q1x,Nx,Ny,Nz);
	q1y=Memory3D(q1y,Nx,Ny,Nz);
	q1z=Memory3D(q1z,Nx,Ny,Nz);
	//allocate memory p1
	p1x=Memory3D(p1x,Nx,Ny,Nz);
	p1y=Memory3D(p1y,Nx,Ny,Nz);
	p1z=Memory3D(p1z,Nx,Ny,Nz);
	//allocate memory Ap
	Apx=Memory3D(Apx,Nx,Ny,Nz);
	Apy=Memory3D(Apy,Nx,Ny,Nz);
	Apz=Memory3D(Apz,Nx,Ny,Nz);

    //temp variables
    MKL_Complex16 temp,temp1,temp2;

	//reset E_inc 
	Current_Reset(Ex_inc, 0, Dim);
	Current_Reset(Ey_inc, 1, Dim);
	Current_Reset(Ez_inc, 2, Dim);
	//r0=E_inc
    cblas_zcopy(Total_inc, Ex_inc[0][0], 1, r0x[0][0], 1);
    cblas_zcopy(Total_inc, Ey_inc[0][0], 1, r0y[0][0], 1);
    cblas_zcopy(Total_inc, Ez_inc[0][0], 1, r0z[0][0], 1);
	//r0 reverse
    temp.real=-1.0;
    temp.imag=0.0;
    cblas_zscal(Total_inc, &temp, r0x[0][0], 1);
    cblas_zscal(Total_inc, &temp, r0y[0][0], 1);
    cblas_zscal(Total_inc, &temp, r0z[0][0], 1);

	//p1 set to 0
	memset(p1x[0][0],0,Total_inc*sizeof(MKL_Complex16));
	memset(p1y[0][0],0,Total_inc*sizeof(MKL_Complex16));
	memset(p1z[0][0],0,Total_inc*sizeof(MKL_Complex16));
	//p1 (p1=A*E_inc)
	X_FFT_Conjugate(p1x, Ex_inc, Ey_inc, Ez_inc, Temp_Current, Green_Function, Dim);
	Y_FFT_Conjugate(p1y, Ex_inc, Ey_inc, Ez_inc, Temp_Current, Green_Function, Dim);
	Z_FFT_Conjugate(p1z, Ex_inc, Ey_inc, Ez_inc, Temp_Current, Green_Function, Dim);
	//reset p1
	Current_Reset(p1x, 0, Dim);
	Current_Reset(p1y, 1, Dim);
	Current_Reset(p1z, 2, Dim);

	//variables
	double v1,v2,alpha,beta,v4,Error;
    //norm v3=||p1||
	cblas_zdotc_sub(Total_inc, p1x[0][0], 1, p1x[0][0], 1, &temp);
	cblas_zdotc_sub(Total_inc, p1y[0][0], 1, p1y[0][0], 1, &temp1);
	cblas_zdotc_sub(Total_inc, p1z[0][0], 1, p1z[0][0], 1, &temp2);
	double v3=temp.real+temp1.real+temp2.real;
    //norm vx=||E_inc||
	cblas_zdotc_sub(Total_inc, Ex_inc[0][0], 1, Ex_inc[0][0], 1, &temp);
	cblas_zdotc_sub(Total_inc, Ey_inc[0][0], 1, Ey_inc[0][0], 1, &temp1);
	cblas_zdotc_sub(Total_inc, Ez_inc[0][0], 1, Ez_inc[0][0], 1, &temp2);
    double vx=temp.real+temp1.real+temp2.real;

	//iterative steps
	for (int m=0; m<Total_inc; m++)
	{

		/*gradient*/
		//v1=v3
		v1=v3;

		//Ap initial
		memset(Apx[0][0],0,Total_inc*sizeof(MKL_Complex16));
		memset(Apy[0][0],0,Total_inc*sizeof(MKL_Complex16));
		memset(Apz[0][0],0,Total_inc*sizeof(MKL_Complex16));
		//FFT(A*p1)
		X_FFT(Apx, p1x, p1y, p1z, Temp_Current, Green_Function, Dim);
		Y_FFT(Apy, p1x, p1y, p1z, Temp_Current, Green_Function, Dim);
		Z_FFT(Apz, p1x, p1y, p1z, Temp_Current, Green_Function, Dim);
		//reset Ap
		Current_Reset(Apx, 0, Dim);
		Current_Reset(Apy, 1, Dim);
		Current_Reset(Apz, 2, Dim);
		//v2=||Ap||
        cblas_zdotc_sub(Total_inc, Apx[0][0], 1, Apx[0][0], 1, &temp);
        cblas_zdotc_sub(Total_inc, Apy[0][0], 1, Apy[0][0], 1, &temp1);
        cblas_zdotc_sub(Total_inc, Apz[0][0], 1, Apz[0][0], 1, &temp2);
        v2=temp.real+temp1.real+temp2.real;

		//step set by alpha=v1/v2
		alpha=v1/v2;
		//update (J=J+alpha*p1)
		temp.real=alpha;
		temp.imag=0.;
        cblas_zaxpy(Total_inc, &temp, p1x[0][0], 1, Current_X[0][0], 1);
        cblas_zaxpy(Total_inc, &temp, p1y[0][0], 1, Current_Y[0][0], 1);
        cblas_zaxpy(Total_inc, &temp, p1z[0][0], 1, Current_Z[0][0], 1);
		//update (r0=r0+alpha*Ap)
        cblas_zaxpy(Total_inc, &temp, Apx[0][0], 1, r0x[0][0], 1);
        cblas_zaxpy(Total_inc, &temp, Apy[0][0], 1, r0y[0][0], 1);
        cblas_zaxpy(Total_inc, &temp, Apz[0][0], 1, r0z[0][0], 1);
		//v4=||r0||
        cblas_zdotc_sub(Total_inc, r0x[0][0], 1, r0x[0][0], 1, &temp);
        cblas_zdotc_sub(Total_inc, r0y[0][0], 1, r0y[0][0], 1, &temp1);
        cblas_zdotc_sub(Total_inc, r0z[0][0], 1, r0z[0][0], 1, &temp2);
        v4=temp.real+temp1.real+temp2.real;

		/*conjugate*/

		//q1 initial
		memset(q1x[0][0],0,Total_inc*sizeof(MKL_Complex16));
		memset(q1y[0][0],0,Total_inc*sizeof(MKL_Complex16));
		memset(q1z[0][0],0,Total_inc*sizeof(MKL_Complex16));
		//q1 by FFT(q1=A'*r0)
		X_FFT_Conjugate(q1x, r0x, r0y, r0z, Temp_Current, Green_Function, Dim);
		Y_FFT_Conjugate(q1y, r0x, r0y, r0z, Temp_Current, Green_Function, Dim);
		Z_FFT_Conjugate(q1z, r0x, r0y, r0z, Temp_Current, Green_Function, Dim);
		//q1 reset
		Current_Reset(q1x, 0, Dim);
		Current_Reset(q1y, 1, Dim);
		Current_Reset(q1z, 2, Dim);
		//v3=||q1||
        cblas_zdotc_sub(Total_inc, q1x[0][0], 1, q1x[0][0], 1, &temp);
        cblas_zdotc_sub(Total_inc, q1y[0][0], 1, q1y[0][0], 1, &temp1);
        cblas_zdotc_sub(Total_inc, q1z[0][0], 1, q1z[0][0], 1, &temp2);
        v3=temp.real+temp1.real+temp2.real;

		//gradient step
		beta=v3/v1;
		//find conjugate direction (p1=beta*p1-q1)
        temp.real=beta;
		temp.imag=0.;
        cblas_zscal(Total_inc, &temp, p1x[0][0], 1);
        cblas_zscal(Total_inc, &temp, p1y[0][0], 1);
        cblas_zscal(Total_inc, &temp, p1z[0][0], 1);
        temp.real=-1.0;
		temp.imag=0.;
        cblas_zaxpy(Total_inc, &temp, q1x[0][0], 1, p1x[0][0], 1);
        cblas_zaxpy(Total_inc, &temp, q1y[0][0], 1, p1y[0][0], 1);
        cblas_zaxpy(Total_inc, &temp, q1z[0][0], 1, p1z[0][0], 1);

		/*truction condition*/
		Error=v4/vx;
		if (Error<CG_THR)
		{
		    printf("CG Error: %lf\n",Error);
			printf("Steps: %d\n\n",m);
			break;
		}

		printf("CG Error: %lf\n",Error);
	}


	//release
	RMemory3D(Apx);
	RMemory3D(Apy);
	RMemory3D(Apz);
	RMemory3D(p1x);
	RMemory3D(p1y);
	RMemory3D(p1z);
	RMemory3D(q1x);
	RMemory3D(q1y);
	RMemory3D(q1z);
	RMemory3D(r0x);
	RMemory3D(r0y);
	RMemory3D(r0z);
	RMemory3D(Temp_Current);

	RMemory3D(Ex_inc);
	RMemory3D(Ey_inc);
	RMemory3D(Ez_inc);
	RMemory3D(Green_Function);
	RMemory3D(Flag_X);
	RMemory3D(Flag_Y);
	RMemory3D(Flag_Z);
}

//BICG-FFT
void BICG_FFT()
{
	//total size
	int Total_inc=Nx*Ny*Nz;
    int Dim[3]={Nx,Ny,Nz};    //dimensions

	//variables
	MKL_Complex16 v1,v2,vx,alpha,beta,s1,s2,s3;
    //temp variables
    MKL_Complex16 temp,temp1,temp2;
	//error
	double Error;
	//complex
	using namespace std;
	complex<double> V1;
	complex<double> V2;
	complex<double> V3;

	//allocate memory Temp_Current
	Temp_Current=Memory3D(Temp_Current,NNx,NNy,NNz);
	//allocate memory r0
	r0x=Memory3D(r0x,Nx,Ny,Nz);
	r0y=Memory3D(r0y,Nx,Ny,Nz);
	r0z=Memory3D(r0z,Nx,Ny,Nz);
	//allocate memory p1
	p1x=Memory3D(p1x,Nx,Ny,Nz);
	p1y=Memory3D(p1y,Nx,Ny,Nz);
	p1z=Memory3D(p1z,Nx,Ny,Nz);
	//allocate memory Ap
	Apx=Memory3D(Apx,Nx,Ny,Nz);
	Apy=Memory3D(Apy,Nx,Ny,Nz);
	Apz=Memory3D(Apz,Nx,Ny,Nz);

	//E_in rest
	Current_Reset(Ex_inc, 0, Dim);
	Current_Reset(Ey_inc, 1, Dim);
	Current_Reset(Ez_inc, 2, Dim);
	//r0=E_inc
    cblas_zcopy(Total_inc, Ex_inc[0][0], 1, r0x[0][0], 1);
    cblas_zcopy(Total_inc, Ey_inc[0][0], 1, r0y[0][0], 1);
    cblas_zcopy(Total_inc, Ez_inc[0][0], 1, r0z[0][0], 1);
	//v1=||r0||
    cblas_zdotu_sub(Total_inc, r0x[0][0], 1, r0x[0][0], 1, &temp);
    cblas_zdotu_sub(Total_inc, r0y[0][0], 1, r0y[0][0], 1, &temp1);
    cblas_zdotu_sub(Total_inc, r0z[0][0], 1, r0z[0][0], 1, &temp2);
	v1.real=temp.real+temp1.real+temp2.real;
	v1.imag=temp.imag+temp1.imag+temp2.imag;

	//p1 set to 0
	memset(p1x[0][0],0,Total_inc*sizeof(MKL_Complex16));
	memset(p1y[0][0],0,Total_inc*sizeof(MKL_Complex16));
	memset(p1z[0][0],0,Total_inc*sizeof(MKL_Complex16));
	//p1=r0
    cblas_zcopy(Total_inc, r0x[0][0], 1, p1x[0][0], 1);
    cblas_zcopy(Total_inc, r0y[0][0], 1, p1y[0][0], 1);
    cblas_zcopy(Total_inc, r0z[0][0], 1, p1z[0][0], 1);
	//p1=p1/v1
    V1=complex<double> (1.0,0.0);
    V2=complex<double> (v1.real,v1.imag);
    V3=V1/V2;
    temp.real=real(V3);
    temp.imag=imag(V3);
    cblas_zscal(Total_inc, &temp, p1x[0][0], 1);
    cblas_zscal(Total_inc, &temp, p1y[0][0], 1);
    cblas_zscal(Total_inc, &temp, p1z[0][0], 1);

	//convergence cond
	vx=v1;
	int converge_time=0;  //times

	//iterative steps
	for (int m=0; m<Total_inc; m++)
	{

		/*gradient*/

		//Ap initial
		memset(Apx[0][0],0,Total_inc*sizeof(MKL_Complex16));
		memset(Apy[0][0],0,Total_inc*sizeof(MKL_Complex16));
		memset(Apz[0][0],0,Total_inc*sizeof(MKL_Complex16));
		//Ap by FFT(A*p1)
		X_FFT(Apx, p1x, p1y, p1z, Temp_Current, Green_Function, Dim);
		Y_FFT(Apy, p1x, p1y, p1z, Temp_Current, Green_Function, Dim);
		Z_FFT(Apz, p1x, p1y, p1z, Temp_Current, Green_Function, Dim);
		//Ap reset
		Current_Reset(Apx, 0, Dim);
		Current_Reset(Apy, 1, Dim);
		Current_Reset(Apz, 2, Dim);
		//v1=Ap.*p1
        cblas_zdotu_sub(Total_inc, Apx[0][0], 1, p1x[0][0], 1, &temp);
        cblas_zdotu_sub(Total_inc, Apy[0][0], 1, p1y[0][0], 1, &temp1);
        cblas_zdotu_sub(Total_inc, Apz[0][0], 1, p1z[0][0], 1, &temp2);
        v1.real=temp.real+temp1.real+temp2.real;
        v1.imag=temp.imag+temp1.imag+temp2.imag;

		//gradient step (alpha=1/v1)
        V1=complex<double> (1.0,0.0);
        V2=complex<double> (v1.real,v1.imag);
        V3=V1/V2;
        alpha.real=real(V3);
        alpha.imag=imag(V3);

		//update J=J+alpha*p1
        cblas_zaxpy(Total_inc, &alpha, p1x[0][0], 1, Current_X[0][0], 1);
        cblas_zaxpy(Total_inc, &alpha, p1y[0][0], 1, Current_Y[0][0], 1);
        cblas_zaxpy(Total_inc, &alpha, p1z[0][0], 1, Current_Z[0][0], 1);
		//reverse
		alpha.real=-alpha.real;
		alpha.imag=-alpha.imag;

		//update r0=r0-alpha*Ap
        cblas_zaxpy(Total_inc, &alpha, Apx[0][0], 1, r0x[0][0], 1);
        cblas_zaxpy(Total_inc, &alpha, Apy[0][0], 1, r0y[0][0], 1);
        cblas_zaxpy(Total_inc, &alpha, Apz[0][0], 1, r0z[0][0], 1);
		//v2=||r0||
        cblas_zdotu_sub(Total_inc, r0x[0][0], 1, r0x[0][0], 1, &temp);
        cblas_zdotu_sub(Total_inc, r0y[0][0], 1, r0y[0][0], 1, &temp1);
        cblas_zdotu_sub(Total_inc, r0z[0][0], 1, r0z[0][0], 1, &temp2);
        v2.real=temp.real+temp1.real+temp2.real;
        v2.imag=temp.imag+temp1.imag+temp2.imag;

		//gradient step (beta=1/v2)
        V1=complex<double> (1.0,0.0);
        V2=complex<double> (v2.real,v2.imag);
        V3=V1/V2;
        beta.real=real(V3);
        beta.imag=imag(V3);

		//find conjuagte direction p1=p1+beta*r0
        cblas_zaxpy(Total_inc, &beta, r0x[0][0], 1, p1x[0][0], 1);
        cblas_zaxpy(Total_inc, &beta, r0y[0][0], 1, p1y[0][0], 1);
        cblas_zaxpy(Total_inc, &beta, r0z[0][0], 1, p1z[0][0], 1);

		/*error trunction*/
		//error
		Error=(v2.real*v2.real+v2.imag*v2.imag)/(vx.real*vx.real+vx.imag*vx.imag);
		Error=sqrt(Error);

		if (Error<CG_THR)
		{
            printf("CG Error: %lf\n",Error);
			printf("Steps: %d\n\n",m);
			break;
		}

		printf("CG Error: %lf\n",Error);
	}


	//release
	RMemory3D(Apx);
	RMemory3D(Apy);
	RMemory3D(Apz);
	RMemory3D(p1x);
	RMemory3D(p1y);
	RMemory3D(p1z);
	RMemory3D(r0x);
	RMemory3D(r0y);
	RMemory3D(r0z);
	RMemory3D(Temp_Current);

	RMemory3D(Ex_inc);
	RMemory3D(Ey_inc);
	RMemory3D(Ez_inc);
	RMemory3D(Green_Function);
	RMemory3D(Flag_X);
	RMemory3D(Flag_Y);
	RMemory3D(Flag_Z);
}

//BICG-Stab-FFT
void BICGSTAB_FFT()
{
	int Total_inc=Nx*Ny*Nz;   //total size
    int Dim[3]={Nx,Ny,Nz};    //dimensions

	//allocate memory Temp_Current
	Temp_Current=Memory3D(Temp_Current,NNx,NNy,NNz);
	//allocate memory r0
	r0x=Memory3D(r0x,Nx,Ny,Nz);
	r0y=Memory3D(r0y,Nx,Ny,Nz);
	r0z=Memory3D(r0z,Nx,Ny,Nz);
	//allocate memory q1
	q1x=Memory3D(q1x,Nx,Ny,Nz);
	q1y=Memory3D(q1y,Nx,Ny,Nz);
	q1z=Memory3D(q1z,Nx,Ny,Nz);
	//allocate memory p1
	p1x=Memory3D(p1x,Nx,Ny,Nz);
	p1y=Memory3D(p1y,Nx,Ny,Nz);
	p1z=Memory3D(p1z,Nx,Ny,Nz);
	//allocate memory Ap
	Apx=Memory3D(Apx,Nx,Ny,Nz);
	Apy=Memory3D(Apy,Nx,Ny,Nz);
	Apz=Memory3D(Apz,Nx,Ny,Nz);
	//allocate memory Aq
	Aqx=Memory3D(Aqx,Nx,Ny,Nz);
	Aqy=Memory3D(Aqy,Nx,Ny,Nz);
	Aqz=Memory3D(Aqz,Nx,Ny,Nz);

	//E_in reset
	Current_Reset(Ex_inc, 0, Dim);
	Current_Reset(Ey_inc, 1, Dim);
	Current_Reset(Ez_inc, 2, Dim);

	//r0;  r = b - A*x; x=0;
    cblas_zcopy(Total_inc, Ex_inc[0][0], 1, r0x[0][0], 1);
    cblas_zcopy(Total_inc, Ey_inc[0][0], 1, r0y[0][0], 1);
    cblas_zcopy(Total_inc, Ez_inc[0][0], 1, r0z[0][0], 1);

	//variables
	MKL_Complex16 bnrm2,alpha,omega,nega_omega,rho,rho_1,beta,temp,temp1,temp2;
	double Error;
	//||b||=||E_inc||
	cblas_zdotc_sub(Total_inc, Ex_inc[0][0], 1, Ex_inc[0][0], 1, &temp);
	cblas_zdotc_sub(Total_inc, Ey_inc[0][0], 1, Ey_inc[0][0], 1, &temp1);
	cblas_zdotc_sub(Total_inc, Ez_inc[0][0], 1, Ez_inc[0][0], 1, &temp2);
	bnrm2.real=temp.real+temp1.real+temp2.real;
	bnrm2.imag=temp.imag+temp1.imag+temp2.imag;

	using namespace std;
	complex<double> V1;
	complex<double> V2;
	complex<double> V3;
	complex<double> V4;
	complex<double> V5;

	//stablized variables
	omega.real=1.0;
    omega.imag=0.0;

	//iterative steps
	for (int m=0; m<Total_inc; m++)
	{
		//rho   = (b'*r0);
		cblas_zdotc_sub(Total_inc, Ex_inc[0][0], 1, r0x[0][0], 1, &temp);
        cblas_zdotc_sub(Total_inc, Ey_inc[0][0], 1, r0y[0][0], 1, &temp1);
        cblas_zdotc_sub(Total_inc, Ez_inc[0][0], 1, r0z[0][0], 1, &temp2);
		rho.real=temp.real+temp1.real+temp2.real;
		rho.imag=temp.imag+temp1.imag+temp2.imag;

        if (m > 0)
        {
            //beta  = ( rho/rho_1 )*( alpha/omega );
            V1=complex<double> (rho.real,rho.imag);
            V2=complex<double> (rho_1.real,rho_1.imag);
            V3=complex<double> (alpha.real,alpha.imag);
            V4=complex<double> (omega.real,omega.imag);
            V5=(V1/V2)*(V3/V4);
            beta.real=real(V5);
            beta.imag=imag(V5);

            //p1 = r0 + beta*( p1 - omega*Ap );
			temp.real=-omega.real;
			temp.imag=-omega.imag;
            cblas_zaxpy(Total_inc, &temp, Apx[0][0], 1, p1x[0][0], 1);
            cblas_zaxpy(Total_inc, &temp, Apy[0][0], 1, p1y[0][0], 1);
            cblas_zaxpy(Total_inc, &temp, Apz[0][0], 1, p1z[0][0], 1);

            cblas_zscal(Total_inc, &beta, p1x[0][0], 1);
            cblas_zscal(Total_inc, &beta, p1y[0][0], 1);
            cblas_zscal(Total_inc, &beta, p1z[0][0], 1);

			temp.real=1.0;
			temp.imag=0.0;
			cblas_zaxpy(Total_inc, &temp, r0x[0][0], 1, p1x[0][0], 1);
			cblas_zaxpy(Total_inc, &temp, r0y[0][0], 1, p1y[0][0], 1);
			cblas_zaxpy(Total_inc, &temp, r0z[0][0], 1, p1z[0][0], 1);
        }
        else
        {
            //p1 = r0;
            cblas_zcopy(Total_inc, r0x[0][0], 1, p1x[0][0], 1);
            cblas_zcopy(Total_inc, r0y[0][0], 1, p1y[0][0], 1);
            cblas_zcopy(Total_inc, r0z[0][0], 1, p1z[0][0], 1);
        }

		/*gradient*/

		//Ap initial
		memset(Apx[0][0],0,Total_inc*sizeof(MKL_Complex16));
		memset(Apy[0][0],0,Total_inc*sizeof(MKL_Complex16));
		memset(Apz[0][0],0,Total_inc*sizeof(MKL_Complex16));
		//Ap by FFT(Ap=A*p1)
		X_FFT(Apx, p1x, p1y, p1z, Temp_Current, Green_Function, Dim);
		Y_FFT(Apy, p1x, p1y, p1z, Temp_Current, Green_Function, Dim);
		Z_FFT(Apz, p1x, p1y, p1z, Temp_Current, Green_Function, Dim);
		//Ap reset
		Current_Reset(Apx, 0, Dim);
		Current_Reset(Apy, 1, Dim);
		Current_Reset(Apz, 2, Dim);

		//alpha = rho / (b'*Ap);
		cblas_zdotc_sub(Total_inc, Ex_inc[0][0], 1, Apx[0][0], 1, &temp);
	    cblas_zdotc_sub(Total_inc, Ey_inc[0][0], 1, Apy[0][0], 1, &temp1);
	    cblas_zdotc_sub(Total_inc, Ez_inc[0][0], 1, Apz[0][0], 1, &temp2);
		temp.real=temp.real+temp1.real+temp2.real;
		temp.imag=temp.imag+temp1.imag+temp2.imag;

		V2=complex<double> (temp.real,temp.imag);
		V1=complex<double> (rho.real,rho.imag);
		V3=V1/V2;
		alpha.real=real(V3);
		alpha.imag=imag(V3);

		//q1 = r0 - alpha*Ap;
		cblas_zcopy(Total_inc, r0x[0][0], 1, q1x[0][0], 1);
		cblas_zcopy(Total_inc, r0y[0][0], 1, q1y[0][0], 1);
		cblas_zcopy(Total_inc, r0z[0][0], 1, q1z[0][0], 1);
		temp.real=-alpha.real;
		temp.imag=-alpha.imag;
		cblas_zaxpy(Total_inc, &temp, Apx[0][0], 1, q1x[0][0], 1);
		cblas_zaxpy(Total_inc, &temp, Apy[0][0], 1, q1y[0][0], 1);
		cblas_zaxpy(Total_inc, &temp, Apz[0][0], 1, q1z[0][0], 1);

		/* gradient again*/

		//Aq initial
		memset(Aqx[0][0],0,Total_inc*sizeof(MKL_Complex16));
		memset(Aqy[0][0],0,Total_inc*sizeof(MKL_Complex16));
		memset(Aqz[0][0],0,Total_inc*sizeof(MKL_Complex16));
		//Aq by FFT(Aq=A*q1)
		X_FFT(Aqx, q1x, q1y, q1z, Temp_Current, Green_Function, Dim);
		Y_FFT(Aqy, q1x, q1y, q1z, Temp_Current, Green_Function, Dim);
		Z_FFT(Aqz, q1x, q1y, q1z, Temp_Current, Green_Function, Dim);
		//Aq reset
		Current_Reset(Aqx, 0, Dim);
		Current_Reset(Aqy, 1, Dim);
		Current_Reset(Aqz, 2, Dim);

		//omega = ( Aq'*q1) / ( Aq'*Aq );
		cblas_zdotc_sub(Total_inc, Aqx[0][0], 1, q1x[0][0], 1, &temp);
		cblas_zdotc_sub(Total_inc, Aqy[0][0], 1, q1y[0][0], 1, &temp1);
		cblas_zdotc_sub(Total_inc, Aqz[0][0], 1, q1z[0][0], 1, &temp2);
		temp.real=temp.real+temp1.real+temp2.real;
		temp.imag=temp.imag+temp1.imag+temp2.imag;
		V1=complex<double> (temp.real,temp.imag);

		cblas_zdotc_sub(Total_inc, Aqx[0][0], 1, Aqx[0][0], 1, &temp);
		cblas_zdotc_sub(Total_inc, Aqy[0][0], 1, Aqy[0][0], 1, &temp1);
		cblas_zdotc_sub(Total_inc, Aqz[0][0], 1, Aqz[0][0], 1, &temp2);
		temp.real=temp.real+temp1.real+temp2.real;
		temp.imag=temp.imag+temp1.imag+temp2.imag;
		V2=complex<double> (temp.real,temp.imag);

		V3=V1/V2;
		omega.real=real(V3);
		omega.imag=imag(V3);

		//x = x + alpha*p1 + omega*q1;
		cblas_zaxpy(Total_inc, &alpha, p1x[0][0], 1, Current_X[0][0], 1);
		cblas_zaxpy(Total_inc, &alpha, p1y[0][0], 1, Current_Y[0][0], 1);
		cblas_zaxpy(Total_inc, &alpha, p1z[0][0], 1, Current_Z[0][0], 1);
		cblas_zaxpy(Total_inc, &omega, q1x[0][0], 1, Current_X[0][0], 1);
		cblas_zaxpy(Total_inc, &omega, q1y[0][0], 1, Current_Y[0][0], 1);
		cblas_zaxpy(Total_inc, &omega, q1z[0][0], 1, Current_Z[0][0], 1);

		//r0 = q1 - omega*Aq;
		cblas_zcopy(Total_inc, q1x[0][0], 1, r0x[0][0], 1);
		cblas_zcopy(Total_inc, q1y[0][0], 1, r0y[0][0], 1);
		cblas_zcopy(Total_inc, q1z[0][0], 1, r0z[0][0], 1);
		temp.real=-omega.real;
		temp.imag=-omega.imag;
		cblas_zaxpy(Total_inc, &temp, Aqx[0][0], 1, r0x[0][0], 1);
		cblas_zaxpy(Total_inc, &temp, Aqy[0][0], 1, r0y[0][0], 1);
		cblas_zaxpy(Total_inc, &temp, Aqz[0][0], 1, r0z[0][0], 1);

		//error = ||r||/||b||;  % check convergence
		cblas_zdotc_sub(Total_inc, r0x[0][0], 1, r0x[0][0], 1, &temp);
		cblas_zdotc_sub(Total_inc, r0y[0][0], 1, r0y[0][0], 1, &temp1);
		cblas_zdotc_sub(Total_inc, r0z[0][0], 1, r0z[0][0], 1, &temp2);
		temp.real=temp.real+temp1.real+temp2.real;
		temp.imag=temp.imag+temp1.imag+temp2.imag;

		/*trunction*/
		Error=sqrt(temp.real/bnrm2.real);
		if (Error<CG_THR)
		{
            printf("CG Error: %lf\n",Error);
			printf("Steps: %d\n\n",m);
			break;
		}

        //rho_1 = rho;
        rho_1.real=rho.real;
        rho_1.imag=rho.imag;

        //show error
		printf("CG Error: %lf\n",Error);
	}


	//release error
	RMemory3D(Apx);
	RMemory3D(Apy);
	RMemory3D(Apz);
	RMemory3D(Aqx);
	RMemory3D(Aqy);
	RMemory3D(Aqz);
	RMemory3D(p1x);
	RMemory3D(p1y);
	RMemory3D(p1z);
	RMemory3D(q1x);
	RMemory3D(q1y);
	RMemory3D(q1z);
	RMemory3D(r0x);
	RMemory3D(r0y);
	RMemory3D(r0z);
	RMemory3D(Temp_Current);

	RMemory3D(Ex_inc);
	RMemory3D(Ey_inc);
	RMemory3D(Ez_inc);
	RMemory3D(Green_Function);
	RMemory3D(Flag_X);
	RMemory3D(Flag_Y);
	RMemory3D(Flag_Z);
}


//CORS-FFT
void CORS_FFT()
{
	int Total_inc=Nx*Ny*Nz;   //total size
    int Dim[3]={Nx,Ny,Nz};    //dimentions

	//allocate memory Temp_Current
	Temp_Current=Memory3D(Temp_Current,NNx,NNy,NNz);
	//allocate memory r0
	r0x=Memory3D(r0x,Nx,Ny,Nz);
	r0y=Memory3D(r0y,Nx,Ny,Nz);
	r0z=Memory3D(r0z,Nx,Ny,Nz);
	//allocate memory q1
	q1x=Memory3D(q1x,Nx,Ny,Nz);
	q1y=Memory3D(q1y,Nx,Ny,Nz);
	q1z=Memory3D(q1z,Nx,Ny,Nz);
	//allocate memory p1
	p1x=Memory3D(p1x,Nx,Ny,Nz);
	p1y=Memory3D(p1y,Nx,Ny,Nz);
	p1z=Memory3D(p1z,Nx,Ny,Nz);
	//allocate memory Ap
	Apx=Memory3D(Apx,Nx,Ny,Nz);
	Apy=Memory3D(Apy,Nx,Ny,Nz);
	Apz=Memory3D(Apz,Nx,Ny,Nz);
	//allocate memory Aq
	Aqx=Memory3D(Aqx,Nx,Ny,Nz);
	Aqy=Memory3D(Aqy,Nx,Ny,Nz);
	Aqz=Memory3D(Aqz,Nx,Ny,Nz);
	//allocate memory r0_c
	r0x_c=Memory3D(r0x_c,Nx,Ny,Nz);
	r0y_c=Memory3D(r0y_c,Nx,Ny,Nz);
	r0z_c=Memory3D(r0z_c,Nx,Ny,Nz);
	//allocate memory q1_c
	q1x_c=Memory3D(q1x_c,Nx,Ny,Nz);
	q1y_c=Memory3D(q1y_c,Nx,Ny,Nz);
	q1z_c=Memory3D(q1z_c,Nx,Ny,Nz);
	//allocate memory p1_c
	p1x_c=Memory3D(p1x_c,Nx,Ny,Nz);
	p1y_c=Memory3D(p1y_c,Nx,Ny,Nz);
	p1z_c=Memory3D(p1z_c,Nx,Ny,Nz);

	//E_in reset
	Current_Reset(Ex_inc, 0, Dim);
	Current_Reset(Ey_inc, 1, Dim);
	Current_Reset(Ez_inc, 2, Dim);

	//r0;  r = b - A*x; x=0;
    cblas_zcopy(Total_inc, Ex_inc[0][0], 1, r0x[0][0], 1);
    cblas_zcopy(Total_inc, Ey_inc[0][0], 1, r0y[0][0], 1);
    cblas_zcopy(Total_inc, Ez_inc[0][0], 1, r0z[0][0], 1);

	//variables
	MKL_Complex16 bnrm2,alpha,rho,rho_1,beta,temp,temp1,temp2;
	double Error;
	//||b||=||E_inc||
	cblas_zdotc_sub(Total_inc, Ex_inc[0][0], 1, Ex_inc[0][0], 1, &temp);
	cblas_zdotc_sub(Total_inc, Ey_inc[0][0], 1, Ey_inc[0][0], 1, &temp1);
	cblas_zdotc_sub(Total_inc, Ez_inc[0][0], 1, Ez_inc[0][0], 1, &temp2);
	bnrm2.real=temp.real+temp1.real+temp2.real;
	bnrm2.imag=temp.imag+temp1.imag+temp2.imag;

	using namespace std;
	complex<double> V1;
	complex<double> V2;
	complex<double> V3;
	complex<double> V4;
	complex<double> V5;

	//iterative steps
	for (int m=0; m<Total_inc; m++)
	{
        /*gradient*/

		//Ap initial
		memset(Apx[0][0],0,Total_inc*sizeof(MKL_Complex16));
		memset(Apy[0][0],0,Total_inc*sizeof(MKL_Complex16));
		memset(Apz[0][0],0,Total_inc*sizeof(MKL_Complex16));
		//Ap by FFT(Ap=A*r0)
		X_FFT(Apx, r0x, r0y, r0z, Temp_Current, Green_Function, Dim);
		Y_FFT(Apy, r0x, r0y, r0z, Temp_Current, Green_Function, Dim);
		Z_FFT(Apz, r0x, r0y, r0z, Temp_Current, Green_Function, Dim);
		//Ap rest
		Current_Reset(Apx, 0, Dim);
		Current_Reset(Apy, 1, Dim);
		Current_Reset(Apz, 2, Dim);
		//rho = (b'*Ap);
		cblas_zdotc_sub(Total_inc, Ex_inc[0][0], 1, Apx[0][0], 1, &temp);
        cblas_zdotc_sub(Total_inc, Ey_inc[0][0], 1, Apy[0][0], 1, &temp1);
        cblas_zdotc_sub(Total_inc, Ez_inc[0][0], 1, Apz[0][0], 1, &temp2);
		rho.real=temp.real+temp1.real+temp2.real;
		rho.imag=temp.imag+temp1.imag+temp2.imag;

        if (m > 0)
        {
            //beta = rho/rho_1
            V1=complex<double> (rho.real,rho.imag);
            V2=complex<double> (rho_1.real,rho_1.imag);
            V3=(V1/V2);
            beta.real=real(V3);
            beta.imag=imag(V3);

            //p1 = r0 + beta*p1_c;
			temp.real=beta.real;
			temp.imag=beta.imag;
            cblas_zcopy(Total_inc, r0x[0][0], 1, p1x[0][0], 1);
            cblas_zcopy(Total_inc, r0y[0][0], 1, p1y[0][0], 1);
            cblas_zcopy(Total_inc, r0z[0][0], 1, p1z[0][0], 1);
			cblas_zaxpy(Total_inc, &temp, p1x_c[0][0], 1, p1x[0][0], 1);
			cblas_zaxpy(Total_inc, &temp, p1y_c[0][0], 1, p1y[0][0], 1);
			cblas_zaxpy(Total_inc, &temp, p1z_c[0][0], 1, p1z[0][0], 1);

            //q1_c = Ap + beta*r0_c;
            cblas_zcopy(Total_inc, Apx[0][0], 1, q1x_c[0][0], 1);
            cblas_zcopy(Total_inc, Apy[0][0], 1, q1y_c[0][0], 1);
            cblas_zcopy(Total_inc, Apz[0][0], 1, q1z_c[0][0], 1);
			cblas_zaxpy(Total_inc, &temp, r0x_c[0][0], 1, q1x_c[0][0], 1);
			cblas_zaxpy(Total_inc, &temp, r0y_c[0][0], 1, q1y_c[0][0], 1);
			cblas_zaxpy(Total_inc, &temp, r0z_c[0][0], 1, q1z_c[0][0], 1);

            //q1 = q1_c + beta*( r0_c + beta*q1 );
            cblas_zscal(Total_inc, &temp, q1x[0][0], 1);
            cblas_zscal(Total_inc, &temp, q1y[0][0], 1);
            cblas_zscal(Total_inc, &temp, q1z[0][0], 1);

			temp1.real=1.0;
			temp1.imag=0.0;
			cblas_zaxpy(Total_inc, &temp1, r0x_c[0][0], 1, q1x[0][0], 1);
			cblas_zaxpy(Total_inc, &temp1, r0y_c[0][0], 1, q1y[0][0], 1);
			cblas_zaxpy(Total_inc, &temp1, r0z_c[0][0], 1, q1z[0][0], 1);

            cblas_zscal(Total_inc, &temp, q1x[0][0], 1);
            cblas_zscal(Total_inc, &temp, q1y[0][0], 1);
            cblas_zscal(Total_inc, &temp, q1z[0][0], 1);

			cblas_zaxpy(Total_inc, &temp1, q1x_c[0][0], 1, q1x[0][0], 1);
			cblas_zaxpy(Total_inc, &temp1, q1y_c[0][0], 1, q1y[0][0], 1);
			cblas_zaxpy(Total_inc, &temp1, q1z_c[0][0], 1, q1z[0][0], 1);
        }
        else
        {
            //p1 = r0;
            cblas_zcopy(Total_inc, r0x[0][0], 1, p1x[0][0], 1);
            cblas_zcopy(Total_inc, r0y[0][0], 1, p1y[0][0], 1);
            cblas_zcopy(Total_inc, r0z[0][0], 1, p1z[0][0], 1);

            //q1_c = Ap;
            cblas_zcopy(Total_inc, Apx[0][0], 1, q1x_c[0][0], 1);
            cblas_zcopy(Total_inc, Apy[0][0], 1, q1y_c[0][0], 1);
            cblas_zcopy(Total_inc, Apz[0][0], 1, q1z_c[0][0], 1);

            //q1 = Ap;
            cblas_zcopy(Total_inc, Apx[0][0], 1, q1x[0][0], 1);
            cblas_zcopy(Total_inc, Apy[0][0], 1, q1y[0][0], 1);
            cblas_zcopy(Total_inc, Apz[0][0], 1, q1z[0][0], 1);
        }

		/*gradient again*/

		//Aq initial
		memset(Aqx[0][0],0,Total_inc*sizeof(MKL_Complex16));
		memset(Aqy[0][0],0,Total_inc*sizeof(MKL_Complex16));
		memset(Aqz[0][0],0,Total_inc*sizeof(MKL_Complex16));
		//Aq by FFT(Aq=A*q1)
		X_FFT(Aqx, q1x, q1y, q1z, Temp_Current, Green_Function, Dim);
		Y_FFT(Aqy, q1x, q1y, q1z, Temp_Current, Green_Function, Dim);
		Z_FFT(Aqz, q1x, q1y, q1z, Temp_Current, Green_Function, Dim);
		//Aq reset
		Current_Reset(Aqx, 0, Dim);
		Current_Reset(Aqy, 1, Dim);
		Current_Reset(Aqz, 2, Dim);

		//alpha = rho / ( b'*Aq );
		V1=complex<double> (rho.real,rho.imag);

		cblas_zdotc_sub(Total_inc, Ex_inc[0][0], 1, Aqx[0][0], 1, &temp);
		cblas_zdotc_sub(Total_inc, Ey_inc[0][0], 1, Aqy[0][0], 1, &temp1);
		cblas_zdotc_sub(Total_inc, Ez_inc[0][0], 1, Aqz[0][0], 1, &temp2);
		temp.real=temp.real+temp1.real+temp2.real;
		temp.imag=temp.imag+temp1.imag+temp2.imag;
		V2=complex<double> (temp.real,temp.imag);

		V3=V1/V2;
		alpha.real=real(V3);
		alpha.imag=imag(V3);

		//p1_c = p1 - alpha*q1;
		cblas_zcopy(Total_inc, p1x[0][0], 1, p1x_c[0][0], 1);
		cblas_zcopy(Total_inc, p1y[0][0], 1, p1y_c[0][0], 1);
		cblas_zcopy(Total_inc, p1z[0][0], 1, p1z_c[0][0], 1);
		temp.real=-alpha.real;
		temp.imag=-alpha.imag;
		cblas_zaxpy(Total_inc, &temp, q1x[0][0], 1, p1x_c[0][0], 1);
		cblas_zaxpy(Total_inc, &temp, q1y[0][0], 1, p1y_c[0][0], 1);
		cblas_zaxpy(Total_inc, &temp, q1z[0][0], 1, p1z_c[0][0], 1);

		//r0_c = q1_c - alpha*Aq;
		cblas_zcopy(Total_inc, q1x_c[0][0], 1, r0x_c[0][0], 1);
		cblas_zcopy(Total_inc, q1y_c[0][0], 1, r0y_c[0][0], 1);
		cblas_zcopy(Total_inc, q1z_c[0][0], 1, r0z_c[0][0], 1);

		cblas_zaxpy(Total_inc, &temp, Aqx[0][0], 1, r0x_c[0][0], 1);
		cblas_zaxpy(Total_inc, &temp, Aqy[0][0], 1, r0y_c[0][0], 1);
		cblas_zaxpy(Total_inc, &temp, Aqz[0][0], 1, r0z_c[0][0], 1);

		//x = x + alpha*(2*p1-alpha*q1);
        temp.real=2*alpha.real;
		temp.imag=2*alpha.imag;
		cblas_zaxpy(Total_inc, &temp, p1x[0][0], 1, Current_X[0][0], 1);
		cblas_zaxpy(Total_inc, &temp, p1y[0][0], 1, Current_Y[0][0], 1);
		cblas_zaxpy(Total_inc, &temp, p1z[0][0], 1, Current_Z[0][0], 1);

        V1=complex<double> (alpha.real,alpha.imag);
        V2=complex<double> (-1.0,0.0);
        V3=V1*V1*V2;
        temp.real=real(V3);
		temp.imag=imag(V3);
		cblas_zaxpy(Total_inc, &temp, q1x[0][0], 1, Current_X[0][0], 1);
		cblas_zaxpy(Total_inc, &temp, q1y[0][0], 1, Current_Y[0][0], 1);
		cblas_zaxpy(Total_inc, &temp, q1z[0][0], 1, Current_Z[0][0], 1);

		//r0 = r0-alpha*(2*q1_c-alpha*Aq);
        temp.real=-2.0*alpha.real;
		temp.imag=-2.0*alpha.imag;
		cblas_zaxpy(Total_inc, &temp, q1x_c[0][0], 1, r0x[0][0], 1);
		cblas_zaxpy(Total_inc, &temp, q1y_c[0][0], 1, r0y[0][0], 1);
		cblas_zaxpy(Total_inc, &temp, q1z_c[0][0], 1, r0z[0][0], 1);

        V1=complex<double> (alpha.real,alpha.imag);
        V3=V1*V1;
        temp.real=real(V3);
		temp.imag=imag(V3);
		cblas_zaxpy(Total_inc, &temp, Aqx[0][0], 1, r0x[0][0], 1);
		cblas_zaxpy(Total_inc, &temp, Aqy[0][0], 1, r0y[0][0], 1);
		cblas_zaxpy(Total_inc, &temp, Aqz[0][0], 1, r0z[0][0], 1);

		//error = ||r0||/||b||;  % check convergence
		cblas_zdotc_sub(Total_inc, r0x[0][0], 1, r0x[0][0], 1, &temp);
		cblas_zdotc_sub(Total_inc, r0y[0][0], 1, r0y[0][0], 1, &temp1);
		cblas_zdotc_sub(Total_inc, r0z[0][0], 1, r0z[0][0], 1, &temp2);
		temp.real=temp.real+temp1.real+temp2.real;
		temp.imag=temp.imag+temp1.imag+temp2.imag;

		/*trunction*/
		Error=sqrt(temp.real/bnrm2.real);
		if (Error<CG_THR)
		{
            printf("CG Error: %lf\n",Error);
			printf("Steps: %d\n\n",m);
			break;
		}

        //rho_1 = rho;
        rho_1.real=rho.real;
        rho_1.imag=rho.imag;

        //show error
		printf("CG Error: %lf\n",Error);
	}


	//释放内存
	RMemory3D(Apx);
	RMemory3D(Apy);
	RMemory3D(Apz);
	RMemory3D(Aqx);
	RMemory3D(Aqy);
	RMemory3D(Aqz);
	RMemory3D(p1x);
	RMemory3D(p1y);
	RMemory3D(p1z);
	RMemory3D(q1x);
	RMemory3D(q1y);
	RMemory3D(q1z);
	RMemory3D(r0x);
	RMemory3D(r0y);
	RMemory3D(r0z);
    RMemory3D(p1x_c);
	RMemory3D(p1y_c);
	RMemory3D(p1z_c);
	RMemory3D(q1x_c);
	RMemory3D(q1y_c);
	RMemory3D(q1z_c);
	RMemory3D(r0x_c);
	RMemory3D(r0y_c);
	RMemory3D(r0z_c);
	RMemory3D(Temp_Current);

	RMemory3D(Ex_inc);
	RMemory3D(Ey_inc);
	RMemory3D(Ez_inc);
	RMemory3D(Green_Function);
	RMemory3D(Flag_X);
	RMemory3D(Flag_Y);
	RMemory3D(Flag_Z);
}

//BICORSSTAB-FFT
void BICORSTAB_FFT()
{
	int Total_inc=Nx*Ny*Nz;   //total sizes
    int Dim[3]={Nx,Ny,Nz};    //dimensions

	//allocate memory Temp_Current
	Temp_Current=Memory3D(Temp_Current,NNx,NNy,NNz);
	//allocate memory r0
	r0x=Memory3D(r0x,Nx,Ny,Nz);
	r0y=Memory3D(r0y,Nx,Ny,Nz);
	r0z=Memory3D(r0z,Nx,Ny,Nz);
	//allocate memory q1
	q1x=Memory3D(q1x,Nx,Ny,Nz);
	q1y=Memory3D(q1y,Nx,Ny,Nz);
	q1z=Memory3D(q1z,Nx,Ny,Nz);
	//allocate memory p1
	p1x=Memory3D(p1x,Nx,Ny,Nz);
	p1y=Memory3D(p1y,Nx,Ny,Nz);
	p1z=Memory3D(p1z,Nx,Ny,Nz);
	//allocate memory Ap
	Apx=Memory3D(Apx,Nx,Ny,Nz);
	Apy=Memory3D(Apy,Nx,Ny,Nz);
	Apz=Memory3D(Apz,Nx,Ny,Nz);
	//allocate memory Aq
	Aqx=Memory3D(Aqx,Nx,Ny,Nz);
	Aqy=Memory3D(Aqy,Nx,Ny,Nz);
	Aqz=Memory3D(Aqz,Nx,Ny,Nz);
	//allocate memory q1
	q1x_c=Memory3D(q1x_c,Nx,Ny,Nz);
	q1y_c=Memory3D(q1y_c,Nx,Ny,Nz);
	q1z_c=Memory3D(q1z_c,Nx,Ny,Nz);
	//allocate memory p1
	p1x_c=Memory3D(p1x_c,Nx,Ny,Nz);
	p1y_c=Memory3D(p1y_c,Nx,Ny,Nz);
	p1z_c=Memory3D(p1z_c,Nx,Ny,Nz);

	//E_in reset
	Current_Reset(Ex_inc, 0, Dim);
	Current_Reset(Ey_inc, 1, Dim);
	Current_Reset(Ez_inc, 2, Dim);

	//r0;  r = b - A*x; x=0;
	cblas_zcopy(Total_inc, Ex_inc[0][0], 1, r0x[0][0], 1);
	cblas_zcopy(Total_inc, Ey_inc[0][0], 1, r0y[0][0], 1);
	cblas_zcopy(Total_inc, Ez_inc[0][0], 1, r0z[0][0], 1);

	//variables
	MKL_Complex16 bnrm2,alpha,omega,rho,rho_1,beta,temp,temp1,temp2;
	double Error;
	//||b||=||E_inc||
	cblas_zdotc_sub(Total_inc, Ex_inc[0][0], 1, Ex_inc[0][0], 1, &temp);
	cblas_zdotc_sub(Total_inc, Ey_inc[0][0], 1, Ey_inc[0][0], 1, &temp1);
	cblas_zdotc_sub(Total_inc, Ez_inc[0][0], 1, Ez_inc[0][0], 1, &temp2);
	bnrm2.real=temp.real+temp1.real+temp2.real;
	bnrm2.imag=temp.imag+temp1.imag+temp2.imag;

	using namespace std;
	complex<double> V1;
	complex<double> V2;
	complex<double> V3;
	complex<double> V4;
	complex<double> V5;

	//iterative steps
	for (int m=0; m<Total_inc; m++)
	{
		/*gradient*/

		//Ap initial
		memset(Apx[0][0],0,Total_inc*sizeof(MKL_Complex16));
		memset(Apy[0][0],0,Total_inc*sizeof(MKL_Complex16));
		memset(Apz[0][0],0,Total_inc*sizeof(MKL_Complex16));
		//Ap by FFT(Ap=A*r0)
		X_FFT(Apx, r0x, r0y, r0z, Temp_Current, Green_Function, Dim);
		Y_FFT(Apy, r0x, r0y, r0z, Temp_Current, Green_Function, Dim);
		Z_FFT(Apz, r0x, r0y, r0z, Temp_Current, Green_Function, Dim);
		//Ap reset
		Current_Reset(Apx, 0, Dim);
		Current_Reset(Apy, 1, Dim);
		Current_Reset(Apz, 2, Dim);

		//rho = (b'*Ap);
		cblas_zdotc_sub(Total_inc, Ex_inc[0][0], 1, Apx[0][0], 1, &temp);
		cblas_zdotc_sub(Total_inc, Ey_inc[0][0], 1, Apy[0][0], 1, &temp1);
		cblas_zdotc_sub(Total_inc, Ez_inc[0][0], 1, Apz[0][0], 1, &temp2);
		rho.real=temp.real+temp1.real+temp2.real;
		rho.imag=temp.imag+temp1.imag+temp2.imag;

		if (m > 0)
		{
			//beta  = ( rho/rho_1 )*( alpha/omega );
			V1=complex<double> (rho.real,rho.imag);
			V2=complex<double> (rho_1.real,rho_1.imag);
			V3=complex<double> (alpha.real,alpha.imag);
			V4=complex<double> (omega.real,omega.imag);
			V5=(V1/V2)*(V3/V4);
			beta.real=real(V5);
			beta.imag=imag(V5);

			//p1 = r0 + beta*( p1 - omega*q1 );
			temp.real=-omega.real;
			temp.imag=-omega.imag;
			cblas_zaxpy(Total_inc, &temp, q1x[0][0], 1, p1x[0][0], 1);
			cblas_zaxpy(Total_inc, &temp, q1y[0][0], 1, p1y[0][0], 1);
			cblas_zaxpy(Total_inc, &temp, q1z[0][0], 1, p1z[0][0], 1);

			cblas_zscal(Total_inc, &beta, p1x[0][0], 1);
			cblas_zscal(Total_inc, &beta, p1y[0][0], 1);
			cblas_zscal(Total_inc, &beta, p1z[0][0], 1);

			temp.real=1.0;
			temp.imag=0.0;
			cblas_zaxpy(Total_inc, &temp, r0x[0][0], 1, p1x[0][0], 1);
			cblas_zaxpy(Total_inc, &temp, r0y[0][0], 1, p1y[0][0], 1);
			cblas_zaxpy(Total_inc, &temp, r0z[0][0], 1, p1z[0][0], 1);

			//q1 = Ap + beta*( q1 - omega*Aq );
			temp.real=-omega.real;
			temp.imag=-omega.imag;
			cblas_zaxpy(Total_inc, &temp, Aqx[0][0], 1, q1x[0][0], 1);
			cblas_zaxpy(Total_inc, &temp, Aqy[0][0], 1, q1y[0][0], 1);
			cblas_zaxpy(Total_inc, &temp, Aqz[0][0], 1, q1z[0][0], 1);

			cblas_zscal(Total_inc, &beta, q1x[0][0], 1);
			cblas_zscal(Total_inc, &beta, q1y[0][0], 1);
			cblas_zscal(Total_inc, &beta, q1z[0][0], 1);

			temp.real=1.0;
			temp.imag=0.0;
			cblas_zaxpy(Total_inc, &temp, Apx[0][0], 1, q1x[0][0], 1);
			cblas_zaxpy(Total_inc, &temp, Apy[0][0], 1, q1y[0][0], 1);
			cblas_zaxpy(Total_inc, &temp, Apz[0][0], 1, q1z[0][0], 1);
		}
		else
		{
			//p1 = r0;
			cblas_zcopy(Total_inc, r0x[0][0], 1, p1x[0][0], 1);
			cblas_zcopy(Total_inc, r0y[0][0], 1, p1y[0][0], 1);
			cblas_zcopy(Total_inc, r0z[0][0], 1, p1z[0][0], 1);
			//q1 = Ap;
			cblas_zcopy(Total_inc, Apx[0][0], 1, q1x[0][0], 1);
			cblas_zcopy(Total_inc, Apy[0][0], 1, q1y[0][0], 1);
			cblas_zcopy(Total_inc, Apz[0][0], 1, q1z[0][0], 1);
		}

		/*gradient again*/

		//Aq initial
		memset(Aqx[0][0],0,Total_inc*sizeof(MKL_Complex16));
		memset(Aqy[0][0],0,Total_inc*sizeof(MKL_Complex16));
		memset(Aqz[0][0],0,Total_inc*sizeof(MKL_Complex16));
		//Aq by FFT(Aq=A*q1)
		X_FFT(Aqx, q1x, q1y, q1z, Temp_Current, Green_Function, Dim);
		Y_FFT(Aqy, q1x, q1y, q1z, Temp_Current, Green_Function, Dim);
		Z_FFT(Aqz, q1x, q1y, q1z, Temp_Current, Green_Function, Dim);
		//Aq rest
		Current_Reset(Aqx, 0, Dim);
		Current_Reset(Aqy, 1, Dim);
		Current_Reset(Aqz, 2, Dim);

		//alpha = rho / ( b'*Aq );
		V1=complex<double> (rho.real,rho.imag);

		cblas_zdotc_sub(Total_inc, Ex_inc[0][0], 1, Aqx[0][0], 1, &temp);
		cblas_zdotc_sub(Total_inc, Ey_inc[0][0], 1, Aqy[0][0], 1, &temp1);
		cblas_zdotc_sub(Total_inc, Ez_inc[0][0], 1, Aqz[0][0], 1, &temp2);
		temp.real=temp.real+temp1.real+temp2.real;
		temp.imag=temp.imag+temp1.imag+temp2.imag;
		V2=complex<double> (temp.real,temp.imag);

		V3=V1/V2;
		alpha.real=real(V3);
		alpha.imag=imag(V3);

		//p1_c = r0 - alpha*q1;
		cblas_zcopy(Total_inc, r0x[0][0], 1, p1x_c[0][0], 1);
		cblas_zcopy(Total_inc, r0y[0][0], 1, p1y_c[0][0], 1);
		cblas_zcopy(Total_inc, r0z[0][0], 1, p1z_c[0][0], 1);
		temp.real=-alpha.real;
		temp.imag=-alpha.imag;
		cblas_zaxpy(Total_inc, &temp, q1x[0][0], 1, p1x_c[0][0], 1);
		cblas_zaxpy(Total_inc, &temp, q1y[0][0], 1, p1y_c[0][0], 1);
		cblas_zaxpy(Total_inc, &temp, q1z[0][0], 1, p1z_c[0][0], 1);

		//q1_c = Ap - alpha*Aq;
		cblas_zcopy(Total_inc, Apx[0][0], 1, q1x_c[0][0], 1);
		cblas_zcopy(Total_inc, Apy[0][0], 1, q1y_c[0][0], 1);
		cblas_zcopy(Total_inc, Apz[0][0], 1, q1z_c[0][0], 1);
		temp.real=-alpha.real;
		temp.imag=-alpha.imag;
		cblas_zaxpy(Total_inc, &temp, Aqx[0][0], 1, q1x_c[0][0], 1);
		cblas_zaxpy(Total_inc, &temp, Aqy[0][0], 1, q1y_c[0][0], 1);
		cblas_zaxpy(Total_inc, &temp, Aqz[0][0], 1, q1z_c[0][0], 1);

		//omega = ( q1_c'*p1_c) / ( q1_c'*q1_c );
		cblas_zdotc_sub(Total_inc, q1x_c[0][0], 1, p1x_c[0][0], 1, &temp);
		cblas_zdotc_sub(Total_inc, q1y_c[0][0], 1, p1y_c[0][0], 1, &temp1);
		cblas_zdotc_sub(Total_inc, q1z_c[0][0], 1, p1z_c[0][0], 1, &temp2);
		temp.real=temp.real+temp1.real+temp2.real;
		temp.imag=temp.imag+temp1.imag+temp2.imag;
		V1=complex<double> (temp.real,temp.imag);

		cblas_zdotc_sub(Total_inc, q1x_c[0][0], 1, q1x_c[0][0], 1, &temp);
		cblas_zdotc_sub(Total_inc, q1y_c[0][0], 1, q1y_c[0][0], 1, &temp1);
		cblas_zdotc_sub(Total_inc, q1z_c[0][0], 1, q1z_c[0][0], 1, &temp2);
		temp.real=temp.real+temp1.real+temp2.real;
		temp.imag=temp.imag+temp1.imag+temp2.imag;
		V2=complex<double> (temp.real,temp.imag);

		V3=V1/V2;
		omega.real=real(V3);
		omega.imag=imag(V3);

		//x = x + alpha*p1 + omega*p1_c;
		cblas_zaxpy(Total_inc, &alpha, p1x[0][0], 1, Current_X[0][0], 1);
		cblas_zaxpy(Total_inc, &alpha, p1y[0][0], 1, Current_Y[0][0], 1);
		cblas_zaxpy(Total_inc, &alpha, p1z[0][0], 1, Current_Z[0][0], 1);
		cblas_zaxpy(Total_inc, &omega, p1x_c[0][0], 1, Current_X[0][0], 1);
		cblas_zaxpy(Total_inc, &omega, p1y_c[0][0], 1, Current_Y[0][0], 1);
		cblas_zaxpy(Total_inc, &omega, p1z_c[0][0], 1, Current_Z[0][0], 1);

		//r0 = p1_c - omega*q1_c;
		cblas_zcopy(Total_inc, p1x_c[0][0], 1, r0x[0][0], 1);
		cblas_zcopy(Total_inc, p1y_c[0][0], 1, r0y[0][0], 1);
		cblas_zcopy(Total_inc, p1z_c[0][0], 1, r0z[0][0], 1);
		temp.real=-omega.real;
		temp.imag=-omega.imag;
		cblas_zaxpy(Total_inc, &temp, q1x_c[0][0], 1, r0x[0][0], 1);
		cblas_zaxpy(Total_inc, &temp, q1y_c[0][0], 1, r0y[0][0], 1);
		cblas_zaxpy(Total_inc, &temp, q1z_c[0][0], 1, r0z[0][0], 1);

		//error = ||r0||/||b||  % check convergence
		cblas_zdotc_sub(Total_inc, r0x[0][0], 1, r0x[0][0], 1, &temp);
		cblas_zdotc_sub(Total_inc, r0y[0][0], 1, r0y[0][0], 1, &temp1);
		cblas_zdotc_sub(Total_inc, r0z[0][0], 1, r0z[0][0], 1, &temp2);
		temp.real=temp.real+temp1.real+temp2.real;
		temp.imag=temp.imag+temp1.imag+temp2.imag;

		/*trunction*/
		Error=sqrt(temp.real/bnrm2.real);
		if (Error<CG_THR)
		{
			printf("迭代步数: %d\n\n",m);
			break;
		}

		//rho_1 = rho;
		rho_1.real=rho.real;
		rho_1.imag=rho.imag;

		//show error
		printf("CG Error: %lf\n",Error);
	}


	//release
	RMemory3D(Apx);
	RMemory3D(Apy);
	RMemory3D(Apz);
	RMemory3D(Aqx);
	RMemory3D(Aqy);
	RMemory3D(Aqz);
	RMemory3D(p1x);
	RMemory3D(p1y);
	RMemory3D(p1z);
	RMemory3D(q1x);
	RMemory3D(q1y);
	RMemory3D(q1z);
	RMemory3D(r0x);
	RMemory3D(r0y);
	RMemory3D(r0z);
	RMemory3D(p1x_c);
	RMemory3D(p1y_c);
	RMemory3D(p1z_c);
	RMemory3D(q1x_c);
	RMemory3D(q1y_c);
	RMemory3D(q1z_c);
	RMemory3D(Temp_Current);

	RMemory3D(Ex_inc);
	RMemory3D(Ey_inc);
	RMemory3D(Ez_inc);
	RMemory3D(Green_Function);
	RMemory3D(Flag_X);
	RMemory3D(Flag_Y);
	RMemory3D(Flag_Z);
}

//VIE-FFT program
void FFT_Program()
{
	//select different iterative solvers
	if (Iteration_Type==1)
	{
		CG_FFT();              //CG_FFT
	}
	else if (Iteration_Type==2)
	{
		BICG_FFT();            //BICG_FFT
	}
	else if (Iteration_Type==3)
	{
		BICGSTAB_FFT();        //BICGSTAB_FFT
	}
	else if (Iteration_Type==4)
	{
		CORS_FFT();            //CORS_FFT
	}
	else
	{
		BICORSTAB_FFT();       //BICORSTAB_FFT
	}
}

#endif
