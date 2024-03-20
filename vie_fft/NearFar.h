#ifndef NEARFAR
#define NEARFAR

#include "Model.h"

//for RCS data
double RCS_Theta[181];
double RCS_Phi[181];

//flag for enegy convergence (no use)
int Energy_Flag=0;
//forward scattering (extinction cross section) 
double Forward_RCS;

//extinction cross section
void FarField1()
{
	int i,j,k;           //dummy variables
//	int angle;           //angles
	double theta_s;      //spherical angle for scattered field (theta)
	double phi_s;        //spherical angle for scattered field (phi)
	double v1,v2,v3,v4;  //angles variables
	double xm,ym,zm;     //position variables
	double vx,vy,vz;     //merged angle variables
	double mv1;          //current variables
	double Forward_RCS_TM, Forward_RCS_TE;  //forward scattering (TM,TE)

	//complex
	using namespace std;
	complex<double> V1;
	complex<double> V2;
	complex<double> Sum_x(0.0,0.0);
	complex<double> Sum_y(0.0,0.0);
	complex<double> Sum_z(0.0,0.0);

	//cell volume
	double Volume=delta_x*delta_y*delta_z;

    //forward scattered directions
	theta_s=Theta;
	phi_s=Phi;

	//angle variables
	v1=cos(PI*theta_s/180.0);
	v2=sin(PI*theta_s/180.0);
	v3=cos(PI*phi_s/180.0);
	v4=sin(PI*phi_s/180.0);

	//merged angle variables
	vx=v2*v3;
	vy=v2*v4;
	vz=v1;

	//summation variables set to 0
	Sum_x=complex<double>(0.0,0.0);
	Sum_y=complex<double>(0.0,0.0);
	Sum_z=complex<double>(0.0,0.0);

	//all current contributions
	for (i=0; i<Nx; i++)
	{
		for (j=0; j<Ny; j++)
		{
			for (k=0; k<Nz; k++)
			{
				//positions of currents
				xm=(i+0.5)*delta_x;
				ym=(j+0.5)*delta_y;
				zm=(k+0.5)*delta_z;

				V1=complex<double>(0.0,k0*(xm*vx+ym*vy+zm*vz));  //phase
				V1=exp(V1);

				//Jx
				V2=complex<double>(Current_X[i][j][k].real,Current_X[i][j][k].imag); 
				Sum_x=Sum_x+V1*V2;  

				//Jy
				V2=complex<double>(Current_Y[i][j][k].real,Current_Y[i][j][k].imag); 
				Sum_y=Sum_y+V1*V2; 

				//JZ
				V2=complex<double>(Current_Z[i][j][k].real,Current_Z[i][j][k].imag);  
				Sum_z=Sum_z+V1*V2;  
			}
		}
	}
	//x projection to theta component
	V2=complex<double>(0.0,k0*0.5*delta_x*vx);
	Sum_x=exp(V2)*Sum_x*Volume;
	//y projection to theta component
	V2=complex<double>(0.0,k0*0.5*delta_y*vy);
	Sum_y=exp(V2)*Sum_y*Volume;
	//z projection to theta component
	V2=complex<double>(0.0,k0*0.5*delta_z*vz);
	Sum_z=exp(V2)*Sum_z*Volume;

	//all current for forward scattering (TM imaginary parts)
	mv1=real(Sum_x)*v1*v3+real(Sum_y)*v1*v4-real(Sum_z)*v2;
	Forward_RCS_TM=Z0*(mv1);

	//all current for forward scattering (TE imaginary parts)
	mv1=real(Sum_x)*v4-real(Sum_y)*v3;
	Forward_RCS_TE=Z0*(mv1);

	//combination
	Forward_RCS=cos(PI*Psi/180.0)*Forward_RCS_TM-sin(PI*Psi/180.0)*Forward_RCS_TE;

}

//near-field calculation for active layer absorption
//metal absorption is excluded
void NearField()
{
	int i,j,k;  //dummy variables

	using namespace std;
	complex<double> V1;
	complex<double> V2;
	complex<double> V3;
	complex<double> V4;
	complex<double> V5;
	complex<double> V6;
	complex<double> V7;
	V5=complex<double>(0.0,-1.0/2.0*Z0/k0);
	V6=complex<double>(1.0,0.0);
	V7=complex<double>(2.0,0.0);

	//files for near-fields
    char string_field_l[256];  //yoz plane
	char string_field_h[256];  //zox plane
	char string_field_v[256];  //xoy plane
	char string_fre[256];      //frequency point
    char string_metal[256];    //metal absorption
    char string_active[256];   //active layer absorption
    char string_grid[256];     //grid

	//file paths
	strcpy(string_field_l, "../result/LField");
	strcpy(string_field_h, "../result/HField");
	strcpy(string_field_v, "../result/VField");
	sprintf(string_fre,"%d",Fre_Point);
	strcat(string_field_l,string_fre);
	strcat(string_field_h,string_fre);
	strcat(string_field_v,string_fre);
    strcat(string_field_l,".txt");
	strcat(string_field_h,".txt");
	strcat(string_field_v,".txt");

	//near-field files
    FILE* file0;
	file0=fopen(string_field_l,"w");

	FILE* file1;
	file1=fopen(string_field_h,"w");

	FILE* file2;
	file2=fopen(string_field_v,"w");

	//absorption and grid information
	strcpy(string_metal, "../result/Ene_Metal");
	strcpy(string_active, "../result/Ene_Active");
	strcpy(string_grid, "../result/Grid");
	strcat(string_metal,string_fre);
	strcat(string_active,string_fre);
	strcat(string_grid,string_fre);
    strcat(string_metal,".txt");
	strcat(string_active,".txt");
	strcat(string_grid,".txt");

	FILE* file3;
	file3=fopen(string_metal,"w");

	FILE* file4;
	file4=fopen(string_active,"w");

	FILE* file5;
	file5=fopen(string_grid,"w");

	double Ene_Metal=0.0;                                    //metal absorption
	double Ene_Active=0.0;                                   //active layer absorption (excluding metals)
	double small_volume=delta_x*delta_y*delta_z;             //cell volume
	double Absorb_Volume=0.0;                                //total volume of active material

	//all currents
	for (i=0; i<Nx-1; i++)
	{
		for (j=0; j<Ny-1; j++)
		{
			for (k=0; k<Nz-1; k++)
			{
				//Ex
				double Abs_Current_X=sqrt(Current_X[i][j][k].real*Current_X[i][j][k].real+
									      Current_X[i][j][k].imag*Current_X[i][j][k].imag);
				V2=complex<double>(Real_EPR[Material[i+1][j][k]]-1.0,Imag_EPR[Material[i+1][j][k]]);
				V3=complex<double>(Real_EPR[Material[i][j][k]]-1.0,Imag_EPR[Material[i][j][k]]);
				V2=(V6/V2+V6/V3);
				V1=V2*V5;
				V4=V7/V2+V6;
				double Abs_Ex=Abs_Current_X*sqrt(real(V1)*real(V1)+imag(V1)*imag(V1));

                if (Material[i+1][j][k]==2 && Material[i][j][k]==2)  //metals
                {
                    Ene_Metal=Ene_Metal+k0/Z0*imag(-V4)*Abs_Ex*Abs_Ex*small_volume;  //absorption of metals
                }
                if (Material[i+1][j][k]==3 && Material[i][j][k]==3)  //active materials
                {
                    Ene_Active=Ene_Active+k0/Z0*imag(-V4)*Abs_Ex*Abs_Ex*small_volume;  //absorption of active materials
					Absorb_Volume=Absorb_Volume+small_volume;
                }

				//Ey
				double Abs_Current_Y=sqrt(Current_Y[i][j][k].real*Current_Y[i][j][k].real+
									      Current_Y[i][j][k].imag*Current_Y[i][j][k].imag);
				V2=complex<double>(Real_EPR[Material[i][j+1][k]]-1.0,Imag_EPR[Material[i][j+1][k]]);
				V3=complex<double>(Real_EPR[Material[i][j][k]]-1.0,Imag_EPR[Material[i][j][k]]);
				V2=(V6/V2+V6/V3);
				V1=V2*V5;
				V4=V7/V2+V6;
				double Abs_Ey=Abs_Current_Y*sqrt(real(V1)*real(V1)+imag(V1)*imag(V1));

                if (Material[i][j+1][k]==2 && Material[i][j][k]==2)  //metals
                {
                    Ene_Metal=Ene_Metal+k0/Z0*imag(-V4)*Abs_Ey*Abs_Ey*small_volume;  //absorption of metals
                }
                if (Material[i][j+1][k]==3 && Material[i][j][k]==3)  //active materials
                {
                    Ene_Active=Ene_Active+k0/Z0*imag(-V4)*Abs_Ey*Abs_Ey*small_volume;  //absorption of active materials
					Absorb_Volume=Absorb_Volume+small_volume;
                }

				//Ez
				double Abs_Current_Z=sqrt(Current_Z[i][j][k].real*Current_Z[i][j][k].real+
	    							      Current_Z[i][j][k].imag*Current_Z[i][j][k].imag);
				V2=complex<double>(Real_EPR[Material[i][j][k+1]]-1.0,Imag_EPR[Material[i][j][k+1]]);
				V3=complex<double>(Real_EPR[Material[i][j][k]]-1.0,Imag_EPR[Material[i][j][k]]);
				V2=(V6/V2+V6/V3);
				V1=V2*V5;
				V4=V7/V2+V6;
				double Abs_Ez=Abs_Current_Z*sqrt(real(V1)*real(V1)+imag(V1)*imag(V1));

                if (Material[i][j][k+1]==2 && Material[i][j][k]==2)  //metals
                {
                    Ene_Metal=Ene_Metal+k0/Z0*imag(-V4)*Abs_Ez*Abs_Ez*small_volume;  //absorption of metals
                }
                if (Material[i][j][k+1]==3 && Material[i][j][k]==3)  //active materials
                {
                    Ene_Active=Ene_Active+k0/Z0*imag(-V4)*Abs_Ez*Abs_Ez*small_volume;  //absorption of active materials
					Absorb_Volume=Absorb_Volume+small_volume;
                }

				if (i==int(Nx/2.0))  //yoz plane(V)
				{
					fprintf(file2, "%.15lf %.15lf %.15lf %.15lf %.15lf %.15lf\n", Current_X[i][j][k].real, Current_X[i][j][k].imag,\
					Current_Y[i][j][k].real, Current_Y[i][j][k].imag, Current_Z[i][j][k].real, Current_Z[i][j][k].imag);
				}

				if (j==int(Ny/2.0))  //zox plane(H)
				{
					fprintf(file1, "%.15lf %.15lf %.15lf %.15lf %.15lf %.15lf\n", Current_X[i][j][k].real, Current_X[i][j][k].imag,\
					Current_Y[i][j][k].real, Current_Y[i][j][k].imag, Current_Z[i][j][k].real, Current_Z[i][j][k].imag);
				}

                if (k==int(Nz/2.0))  //xoy plane(L)
				{
					fprintf(file0, "%.15lf %.15lf %.15lf %.15lf %.15lf %.15lf\n", Current_X[i][j][k].real, Current_X[i][j][k].imag,\
					Current_Y[i][j][k].real, Current_Y[i][j][k].imag, Current_Z[i][j][k].real, Current_Z[i][j][k].imag);
				}
			}
		}
	}

	//output results

	//metal extinction and absorption cross section as a function of wavelength
    fprintf(file3, "%d %.15lf %.15lf\n",WBegin+(Fre_Point)*WInterval,
			(Forward_RCS/Cross_Section),
			(Ene_Metal*Z0/Cross_Section));
	//absorption of active layer per unit volume
    fprintf(file4, "%d %.15lf\n",WBegin+(Fre_Point)*WInterval,
			3.0*Ene_Active/Absorb_Volume); //output results
	//modeling information
	fprintf(file5, "%d %d %d %d %.15lf %.15lf %.15lf\n", Fre_Point, Nx, Ny, Nz, delta_x, delta_y, delta_z);

    fclose(file0);
	fclose(file1);
	fclose(file2);
	fclose(file3);
	fclose(file4);
	fclose(file5);

	//release current
	RMemory3D(Current_X);
	RMemory3D(Current_Y);
	RMemory3D(Current_Z);
	//realse material
	RMemory3D(Material);
}


//far-field scattering pattern
void FarField2()
{
	int i,j,k;           //dummy variables
	int angle;           //angles
	double theta_s;      //spherical angle for scattered field (theta)
	double phi_s;        //spherical angle for scattered field (phi)
	double v1,v2,v3,v4;  //angles variables
	double xm,ym,zm;     //position variables
	double vx,vy,vz;     //merged angle variables
	double mv1,mv2;      //absolute values of current variables

	//complex
	using namespace std;
	complex<double> V1;
	complex<double> V2;
	complex<double> Sum_x(0.0,0.0);
	complex<double> Sum_y(0.0,0.0);
	complex<double> Sum_z(0.0,0.0);

    //E plane
	phi_s=Phi;
	//cell volume
	double Volume=delta_x*delta_y*delta_z;

	//theta angles
	for (angle=0; angle<=180; angle++)  //scattering angle
	{
		//spherical angles (theta)
		theta_s=angle;

		//angles variables
		v1=cos(PI*theta_s/180.0);
		v2=sin(PI*theta_s/180.0);
		v3=cos(PI*phi_s/180.0);
		v4=sin(PI*phi_s/180.0);
		//merged angle variables
		vx=v2*v3;
		vy=v2*v4;
		vz=v1;
		//summation set to 0
		Sum_x=complex<double>(0.0,0.0);
		Sum_y=complex<double>(0.0,0.0);
		Sum_z=complex<double>(0.0,0.0);
		//all currents
		for (i=0; i<Nx; i++)
		{
			for (j=0; j<Ny; j++)
			{
				for (k=0; k<Nz; k++)
				{
					//position
					xm=(i+0.5)*delta_x;
					ym=(j+0.5)*delta_y;
					zm=(k+0.5)*delta_z;

					V1=complex<double>(0.0,k0*(xm*vx+ym*vy+zm*vz));  //phase
					V1=exp(V1);

					//Jx
					V2=complex<double>(Current_X[i][j][k].real,Current_X[i][j][k].imag);  
					Sum_x=Sum_x+V1*V2;  

					//Jy
					V2=complex<double>(Current_Y[i][j][k].real,Current_Y[i][j][k].imag);  
					Sum_y=Sum_y+V1*V2;  

					//Jz
					V2=complex<double>(Current_Z[i][j][k].real,Current_Z[i][j][k].imag); 
					Sum_z=Sum_z+V1*V2; 
				}
			}
		}
		//x projection to theta component
		V2=complex<double>(0.0,k0*0.5*delta_x*vx);
		Sum_x=exp(V2)*Sum_x*Volume;
		//y projection to theta component
		V2=complex<double>(0.0,k0*0.5*delta_y*vy);
		Sum_y=exp(V2)*Sum_y*Volume;
		//z projection to theta component
		V2=complex<double>(0.0,k0*0.5*delta_z*vz);
		Sum_z=exp(V2)*Sum_z*Volume;
		//total contribution
		mv1=real(Sum_x)*v1*v3+real(Sum_y)*v1*v4-real(Sum_z)*v2;
		mv2=imag(Sum_x)*v1*v3+imag(Sum_y)*v1*v4-imag(Sum_z)*v2;
		//RCS results for theta components
		RCS_Theta[angle]=k0*k0*Z0*Z0*(mv1*mv1+mv2*mv2)/(4.0*PI);
		//log10 scale
		RCS_Theta[angle]=10.0*log10(RCS_Theta[angle]/(2*PI/k0)/(2*PI/k0));
	}

	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	//H plane
	phi_s=Phi+90.;

	//theta angles
	for (angle=0; angle<=180; angle++)  //scattering angle
	{
		//spherical angles (theta)
		theta_s=angle;

		//angles variables
		v1=cos(PI*theta_s/180.0);
		v2=sin(PI*theta_s/180.0);
		v3=cos(PI*phi_s/180.0);
		v4=sin(PI*phi_s/180.0);
		//merged angle variables
		vx=v2*v3;
		vy=v2*v4;
		vz=v1;
		//summation set to 0
		Sum_x=complex<double>(0.0,0.0);
		Sum_y=complex<double>(0.0,0.0);
		Sum_z=complex<double>(0.0,0.0);
		//all currents
		for (i=0; i<Nx; i++)
		{
			for (j=0; j<Ny; j++)
			{
				for (k=0; k<Nz; k++)
				{
					//positions
					xm=(i+0.5)*delta_x;
					ym=(j+0.5)*delta_y;
					zm=(k+0.5)*delta_z;

					V1=complex<double>(0.0,k0*(xm*vx+ym*vy+zm*vz));  //phase
					V1=exp(V1);

					//Jx
					V2=complex<double>(Current_X[i][j][k].real,Current_X[i][j][k].imag); 
					Sum_x=Sum_x+V1*V2;  

					//Jy
					V2=complex<double>(Current_Y[i][j][k].real,Current_Y[i][j][k].imag); 
					Sum_y=Sum_y+V1*V2;  

				}
			}
		}
		//x projection to phi component
		V2=complex<double>(0.0,k0*0.5*delta_x*vx);
		Sum_x=exp(V2)*Sum_x*Volume;
		//y projection to phi component
		V2=complex<double>(0.0,k0*0.5*delta_y*vy);
		Sum_y=exp(V2)*Sum_y*Volume;
		//total contribution
		mv1=real(Sum_x)*v4-real(Sum_y)*v3;
		mv2=imag(Sum_x)*v4-imag(Sum_y)*v3;
		//RCS results
		RCS_Phi[angle]=k0*k0*Z0*Z0*(mv1*mv1+mv2*mv2)/(4.0*PI);
		//log10 scale
		RCS_Phi[angle]=10*log10(RCS_Phi[angle]/(2*PI/k0)/(2*PI/k0));
	}

	//far field files
	char string_rcs[256];  //path
	char string_fre[256];  //frequency
	//path
	strcpy(string_rcs, "../result/RCS");
	sprintf(string_fre,"%d",Fre_Point);
	strcat(string_rcs,string_fre);
	strcat(string_rcs,".txt");

	//write to file
	FILE* file;
	file=fopen(string_rcs,"w");  //open file

	for (i=0; i<=180; i++)
	{
		fprintf(file,"%.15lf %.15lf\n",RCS_Theta[i],RCS_Phi[i]);  //RCS at Theta and Phi plane
	}

	fclose(file);  //close file
}

#endif
