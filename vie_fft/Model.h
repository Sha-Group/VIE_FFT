#ifndef MODEL
#define MODEL
#include "common.h"
#include "ReadFile.h"
#include <vector>

//Object size [0,Box_L]x[0,Box_W]x[0,Box_H]
//current index starts from 0. Current unknowns Jx,Jy,Jz along the x,y,z directions are N_{x,y,z}-1, respectively.
//charge index starts from 0. Charge unknowns are equal to grid sizes.
//material index starts from 0, simple averaged scheme is used epr_i(J(i))=0.5*(epr_{i+1}+epr_{i})

//global variables
MKL_Complex16*** Current_X;       //Current along x
MKL_Complex16*** Current_Y;       //Current along y
MKL_Complex16*** Current_Z;       //Current along z
MKL_Complex16*** Green_Function;  //scalar Green function
unsigned char*** Material;        //material index
MKL_Complex16*** Ex_inc;          //incident Ex field
MKL_Complex16*** Ey_inc;          //incident Ey field
MKL_Complex16*** Ez_inc;          //incident Ez field
unsigned char*** Flag_X;          //for polarization current identifier (Jx)
unsigned char*** Flag_Y;          //for polarization current identifier (Jy)
unsigned char*** Flag_Z;          //for polarization current identifier (Jz)

double k0;                        //wavenumber in free space
double Frequency;                 //frequency

double delta_x;                   //spatial step along x
double delta_y;                   //spatial step along y
double delta_z;                   //spatial step along z

int Nx;                           //Grid number along x
int Ny;                           //Grid number along y
int Nz;                           //Grid number along z
int NNx;                          //Unknowns along x for FFT
int NNy;                          //Unknowns along y for FFT
int NNz;                          //Unknowns along z for FFT

//***points per direction x y z (control staircase error)
int St_Nx=30;
int St_Ny=30;
int St_Nz=30;

double Theta=0;            //***spherical angle (wrt z axis)
double Phi=0;              //***spherical angle (wrt x axis)
double Psi=0;              //***TM(Hx,Hy)-0 degree,TE(Ex,Ey)-90 degree

//______________________________________________________________________________________//
//______________________________________________________________________________________//

//Metallic sphere
double Dia=60/double(WScale);        //***sphere diameter
double Cross_Section=PI*Dia*Dia/4.0; //***geometric cross section of sphere
double Box_L=Dia;				     //***length of bouncing box
double Box_W=Dia;				     //***width of bouncing box
double Box_H=Dia;				     //***height of bouncing box

//______________________________________________________________________________________//
//______________________________________________________________________________________//

////Solar Cell
//double Dia=20/double(WScale);       //***sphere diameter
//double Cross_Section=PI*Dia*Dia/4.0;//***geometric cross section of sphere
//double Inter_void=Dia;              //***void size
//double Inter_Sphx=0.2*Dia;          //***distance between adjacent spheres
//const int Num_Sphx=5;               //***sphere number
//
//double Box_L=Num_Sphx*Dia+(Num_Sphx-1)*Inter_Sphx+Inter_void;  //***length of bouncing box
//double Box_W=Dia+Inter_void;          //***width of bouncing box
//double Thick1=25/double(WScale);      //***spacer layer thickness
//double Thick2=50/double(WScale);      //***active layer thickness
//double Box_H=(Thick1+Thick2);         //***height of bouncing box

//______________________________________________________________________________________//
//______________________________________________________________________________________//

int PPW_x;     //points per wavelength (x)
int PPW_y;     //points per wavelength (y)
int PPW_z;     //points per wavelength (z)

const int NUM_DIE=4;                          //***material types（including dummy material）
double Real_EPR[NUM_DIE]={1.1,2.0,2.0,2.0};   //***Index[0] is for dummy material (real part of permittivity) for avoiding zero denominator 1/(\epsilon_r-1)
double Imag_EPR[NUM_DIE]={0.0,0.0,0.0,0.0};   //***Index[0] is for dummy material (imaginary part of permittivity)

int Fre_Point;  //current frequency point

//read parameters and allocate memory
void Parameter_Input(int fre_point)
{
	int i;

	Fre_Point=fre_point-1;  //current frequency point

	//assign materials
	Real_EPR[1]=m1_r[Fre_Point];
	Imag_EPR[1]=m1_i[Fre_Point];
	Real_EPR[2]=m2_r[Fre_Point];
	Imag_EPR[2]=m2_i[Fre_Point];
	Real_EPR[3]=m3_r[Fre_Point];
	Imag_EPR[3]=m3_i[Fre_Point];
//	Real_EPR[2]=4.0;
//	Imag_EPR[2]=0.0;

	//assign paramters
	double lambda0=Inc_Wavelength[Fre_Point];   //free space wavelength
	k0=2*PI/lambda0;                            //free space wavenumber
	Frequency=V0/(lambda0);                     //incident frequency

	//finding smallest grid size
	using namespace std;
	complex<double> k_d;

	double lambda_small=SIGMAX;  //minimum dielectric wavelength 
	double real_kd;              //real part Re(k_d)
	double imag_kd;              //imaginary part Im(k_d)

	for (i=0; i<NUM_DIE; i++)
	{
		k_d=complex<double>(Real_EPR[i],Imag_EPR[i]);  //complex permittivity
		k_d=k0*sqrt(k_d);  //dielectric wavenumber

		real_kd=fabs(real(k_d));
		imag_kd=fabs(imag(k_d));
		real_kd=2*PI/real_kd;      //dielectric wavelength
        imag_kd=1.0/(imag_kd+EPS); //decay length

		//find minimum wavelength（skin depth & oscillation）
		if (real_kd<lambda_small)
		{
			lambda_small=real_kd;
		}
		if (imag_kd<lambda_small)
		{
			lambda_small=imag_kd;
		}
	}

	//spatial steps
	delta_x=lambda_small/PPW_x;
	delta_y=lambda_small/PPW_y;
	delta_z=lambda_small/PPW_z;

	//staircase error control (points per direction)
	double stair_x=Dia/St_Nx;
	double stair_y=Dia/St_Ny;
	double stair_z=Dia/St_Nz;

	//skin depth & oscillation & staircase error
	if (delta_x>stair_x)
	{
		delta_x=stair_x;
	}
	if (delta_y>stair_y)
	{
		delta_y=stair_y;
	}
	if (delta_z>stair_z)
	{
		delta_z=stair_z;
	}

	//grid numbers
	Nx=int(Box_L/delta_x)+2;
	Ny=int(Box_W/delta_y)+2;
	Nz=int(Box_H/delta_z)+2;

	//unknowns
    NNx=2*Nx;
	NNy=2*Ny;
	NNz=2*Nz;

	//allocate memory
	Current_X=Memory3D(Current_X,Nx,Ny,Nz);
	Current_Y=Memory3D(Current_Y,Nx,Ny,Nz);
	Current_Z=Memory3D(Current_Z,Nx,Ny,Nz);
    Green_Function=Memory3D(Green_Function,NNx,NNy,NNz);
	Material=Memory3D(Material,Nx,Ny,Nz);
	Ex_inc=Memory3D(Ex_inc,Nx,Ny,Nz);
	Ey_inc=Memory3D(Ey_inc,Nx,Ny,Nz);
	Ez_inc=Memory3D(Ez_inc,Nx,Ny,Nz);
	Flag_X=Memory3D(Flag_X,Nx,Ny,Nz);
	Flag_Y=Memory3D(Flag_Y,Nx,Ny,Nz);
	Flag_Z=Memory3D(Flag_Z,Nx,Ny,Nz);

	//initial values
	int Total=NNx*NNy*NNz;
	int Total_inc=Nx*Ny*Nz;
	//current, green function and incident fields (0)
	memset(Current_X[0][0],0,Total_inc*sizeof(MKL_Complex16));
	memset(Current_Y[0][0],0,Total_inc*sizeof(MKL_Complex16));
	memset(Current_Z[0][0],0,Total_inc*sizeof(MKL_Complex16));
	memset(Green_Function[0][0],0,Total*sizeof(MKL_Complex16));
	memset(Ex_inc[0][0],0,Total_inc*sizeof(MKL_Complex16));
	memset(Ey_inc[0][0],0,Total_inc*sizeof(MKL_Complex16));
	memset(Ez_inc[0][0],0,Total_inc*sizeof(MKL_Complex16));
	//material settings
	memset(Material[0][0],0,Total_inc*sizeof(unsigned char));
	//current identifier
	memset(Flag_X[0][0],0,Total_inc*sizeof(unsigned char));
	memset(Flag_Y[0][0],0,Total_inc*sizeof(unsigned char));
	memset(Flag_Z[0][0],0,Total_inc*sizeof(unsigned char));

	//show information on screen
	printf("frequency points：%d\n\n",Fre_Point);
	printf("wavelength: %lf nm\n",lambda0*WScale);
	printf("incident spherical angles: %lf %lf\n",Theta,Phi);
	printf("polarization angle: %lf\n",Psi);
	printf("grid sizes: %d X %d X %d\n",Nx,Ny,Nz);
	//printf("离散步长是: %lf X %lf X %lf\n",delta_x,delta_y,delta_z);
	printf("\n");
	printf("\n");

}

//material assignment（return material index）
int Model_Identifier(const double posx,const double posy,const double posz)
{

//______________________________________________________________________________________//
//______________________________________________________________________________________//

	//one sphere
	double x_center=Box_L/2.0;
	double y_center=Box_W/2.0;
	double z_center=Box_H/2.0;

    double v1=sqrt(pow((posx-x_center),2.0)+pow((posy-y_center),2.0)+pow((posz-z_center),2.0));  //distance from the sphere center

	if (v1<=0.5*Dia) //inner part of the sphere
	{
		return 2;    //metal
	}
	else
	{
		return 0;   //air
	}

//______________________________________________________________________________________//
//______________________________________________________________________________________//

//	//two-layer organic solar cells
//    double x_center0=Box_L/2.0;
//	double y_center0=Box_W/2.0;
//	double z_center0=Box_H/2.0;
//
//	double x_center1=Box_L/2.0;
//	double y_center1=Box_W/2.0;
//	double z_center1=Thick1/2;
//	double z_side1=Thick1/2;
//
//    double x_center2=Box_L/2.0;
//	double y_center2=Box_W/2.0;
//	double z_center2=Thick1+Thick2/2;
//	double z_side2=Thick2/2;
//
//	//sohere center（5 spheres）
//	double x_dis=Dia+Inter_Sphx;
//	double sp_cx[Num_Sphx]={x_center0-2*x_dis,x_center0-x_dis,x_center0,\
//		                    x_center0+x_dis,x_center0+2*x_dis};
//	double sp_cy[Num_Sphx]={y_center0,y_center0,y_center0,y_center0,y_center0};
////	double sp_cz[Num_Sphx]={Thick1-Dia/2,Dia/2,Thick1-Dia/2,Dia/2,Thick1-Dia/2};
//	double sp_cz[Num_Sphx]={Thick1/2,Thick1/2,Thick1/2,Thick1/2,Thick1/2};
//
//	if (fabs(posx-x_center1)<x_center1 && fabs(posy-y_center1)<y_center1 &&\
//		fabs(posz-z_center1)<z_side1) //within PEDOT
//	{
//		double radius_square=(0.5*Dia)*(0.5*Dia);  //radius square
//		int sp_flag=0;                             //flag 
//
//		for (int m=0; m<Num_Sphx; m++)
//		{
//			double r_dis=pow((posx-sp_cx[m]),2.0)+pow((posy-sp_cy[m]),2.0)+\
//				         pow((posz-sp_cz[m]),2.0);
//
//			//within the sphere
//			if (float(r_dis)<float(radius_square))
//			{
//				sp_flag=1;   //flag set for metals
//				return 2;    //metals
//			}
//		}
//
//		if (sp_flag==0) //out of the metal
//		{
//			return 1;  //within PEDOT
//		}
//
//	}
//	else if (fabs(posx-x_center2)<x_center2 && fabs(posy-y_center2)<y_center2 &&\
//		     fabs(posz-z_center2)<z_side2)
//	{
//		return 3;  //within active layer
//	}
//
//	else
//	{
//		return 0;  //air
//	}

}



//assign materials
void GetMaterial()
{
	int i,j,k;

	//model files to check grids
//	FILE* file;
//	file=fopen("../result/Model.txt","w");

	for (i=0; i<Nx; i++)
	{
		for (j=0; j<Ny; j++)
		{
			for (k=0; k<Nz; k++)
			{
				double posx=i*delta_x+0.5*delta_x;
				double posy=j*delta_y+0.5*delta_y;
				double posz=k*delta_z+0.5*delta_z;

				int flag=Model_Identifier(posx, posy, posz);  //material identifier
				Material[i][j][k]=(unsigned char)flag;        //assign material

				//fprintf(file,"%d %d %d %d\n",i,j,k,flag);   //save to file
			}
		}
	}
//	fclose(file);
//	exit(0);
}

//get identifiers for currents
void GetFlag()
{
	//dummy variables
	int i,j,k;

	//position variables
	double posx;
	double posy;
	double posz;

	for (i=0; i<Nx; i++)
	{
		for (j=0; j<Ny; j++)
		{
			for (k=0; k<Nz; k++)
			{
				//x rooftop basis
				posx=(i+1)*delta_x;
				posy=(j+0.5)*delta_y;
				posz=(k+0.5)*delta_z;

				//Flag_X for Jx
				int flagx1=Model_Identifier(posx, posy, posz);
				int flagx2=Model_Identifier(posx-0.5*delta_x, posy, posz);
				int flagx3=Model_Identifier(posx+0.5*delta_x, posy, posz);

				if (flagx1!=0 && flagx2!=0 && flagx3!=0)
				{
					Flag_X[i][j][k]=(unsigned char)1;  //inner
				}

				//y rooftop basis
				posx=(i+0.5)*delta_x;
				posy=(j+1)*delta_y;
				posz=(k+0.5)*delta_z;

				//Flag_Y for Jy
				int flagy1=Model_Identifier(posx, posy, posz);
				int flagy2=Model_Identifier(posx, posy-0.5*delta_y, posz);
				int flagy3=Model_Identifier(posx, posy+0.5*delta_y, posz);

				if (flagy1!=0 && flagy2!=0 && flagy3!=0)
				{
					Flag_Y[i][j][k]=(unsigned char)1;  //inner
				}

				//z rooftop basis
				posx=(i+0.5)*delta_x;
				posy=(j+0.5)*delta_y;
				posz=(k+1)*delta_z;

				//Flag_Z for Jz
				int flagz1=Model_Identifier(posx, posy, posz);
				int flagz2=Model_Identifier(posx, posy, posz-0.5*delta_z);
				int flagz3=Model_Identifier(posx, posy, posz+0.5*delta_z);

				if (flagz1!=0 && flagz2!=0 && flagz3!=0)
				{
					Flag_Z[i][j][k]=(unsigned char)1;  //inner
				}
			}
		}
	}
}

#endif


/*
	//偶极子辐射
	//一层结构有机太阳能
    double x_center0=Box_L/2.0;
	double y_center0=Box_W/2.0;
	double z_center0=Box_H/2.0;

	//球中心
//	double x_dis=Dia+Inter_Sphx;
	double sp_cx[Num_Sphx]={x_center0};
	double sp_cy[Num_Sphx]={y_center0};
	double sp_cz[Num_Sphx]={z_center0};

	if (fabs(posx-x_center0)<x_center0 && fabs(posy-y_center0)<y_center0 &&\
		fabs(posz-z_center0)<z_center0) //长方体内部1
	{
		double radius_square=(0.5*Dia)*(0.5*Dia);  //半径平方
		int sp_flag=0;                             //球内部标志

		for (int m=0; m<Num_Sphx; m++)
		{
			double r_dis=pow((posx-sp_cx[m]),2.0)+pow((posy-sp_cy[m]),2.0)+\
				         pow((posz-sp_cz[m]),2.0);

			//球内部
			if (float(r_dis)<float(radius_square))
			{
				sp_flag=1;   //金属
				return 2;    //金属内部
			}
		}

		if (sp_flag==0)
		{
			return 1;  //介质内部
		}
	}

	else
	{
		return 0;  //外部
	}
*/


/*
	//六角
	int m;  //循环变量

	//Box中心
	double x_center=Box_L/2.0;
	double y_center=Box_W/2.0;
	double z_center=Box_H/2.0;

	//7个球中心
	double cx[7]={-Dia/2.0,Dia/2.0,-Dia,0.0,Dia,-Dia/2.0,Dia/2.0};
	double cy[7]={sqrt(3.0)/2.0*Dia,sqrt(3.0)/2.0*Dia,0.0,0.0,0.0,\
		         -sqrt(3.0)/2.0*Dia,-sqrt(3.0)/2.0*Dia};
	//球中心移位
	for (m=0; m<7; m++)
	{
		cx[m]=cx[m]+x_center;
		cy[m]=cy[m]+y_center;
	}

	//六角形顶点坐标（循环）
	double tri_cx[7]={-Dia/2.0,-Dia,-Dia/2.0,Dia/2.0,Dia,Dia/2.0,-Dia/2.0};
	double tri_cy[7]={sqrt(3.0)/2.0*Dia,0.0,-sqrt(3.0)/2.0*Dia,-sqrt(3.0)/2.0*Dia,\
		              0.0,sqrt(3.0)/2.0*Dia,sqrt(3.0)/2.0*Dia};
	//六角形中心移位
	for (m=0; m<7; m++)
	{
		tri_cx[m]=tri_cx[m]+x_center;
		tri_cy[m]=tri_cy[m]+y_center;
	}

	//判断是否六角形内部
	double sum_angle=0.0;
	for (m=0; m<6; m++)
	{
		double x1=tri_cx[m]-posx;
		double x2=tri_cx[m+1]-posx;
		double y1=tri_cy[m]-posy;
		double y2=tri_cy[m+1]-posy;

		double abs1=sqrt(x1*x1+y1*y1);
		double abs2=sqrt(x2*x2+y2*y2);
		double angle;

		if (abs1>EPS && abs2>EPS)
		{
			angle=(x1*x2+y1*y2)/abs1/abs2;
		}
		else
		{
			angle=0.0;
		}

		sum_angle=sum_angle+acos(angle);
	}

	int flag=-1;  //标志

	//在六角形内部
	if (fabs(sum_angle-2*PI)<EPS)
	{
		//球循环
		for (m=0; m<7; m++)
		{
			//某个球内部
			if ((pow((posx-cx[m]),2.0)+pow((posy-cy[m]),2.0))
				<=pow((Dia/2.0),2.0))
			{
				flag=0;  //介质
			}
		}

		if (flag==-1)  //不再任何一个球内部
		{
			flag=2;  //金属
		}
	}

    //不是金属
    if (!(flag==2))
	{
		flag=0;  //介质
	}

	//超出box边界
	if (fabs(posx-x_center)>x_center || fabs(posy-y_center)>y_center || fabs(posz-z_center)>z_center)  //v1<r
	{
		flag=0;  //外部(非金属也非介质)
	}

	return flag;

*/


////有机太阳能建模2
//double Cross_Section=PI*Dia*Dia/4.0;
//double Inter_void=4*Dia;   //***控制传输层两边大小
//double Inter_Sphx=0.2*Dia; //***纳米球间隔
//const int Num_Sphx=1;      //***纳米球个数
//
//double Box_L=Dia+Inter_void;          //***包围介质的Box长x
//double Box_W=Dia+Inter_void;          //***包围介质的Box宽y
//double Box_H=Dia+Inter_void;          //***包围介质的Box高z
//
////不规则金属颗粒
//double Box_L=3*Dia;              //***包围介质的Box长x
//double Box_W=(sqrt(3.0)+1.0)*Dia;//***包围介质的Box宽y
//double Box_H=Dia;                //***包围介质的Box高z

////三角分形天线
//const int F_Level=3;                 //***层数
//const int N_Fractal=1+3+9;           //***多少个空心小单元
//double Dia=12.5/double(WScale);      //***最小尺寸
//double Box_L=8*Dia;                  //***包围介质的Box长x
//double Box_W=Box_L*sqrt(3.0)/2.0;    //***包围介质的Box宽y
//double Box_H=30/double(WScale);	     //***包围介质的Box高z
//double Cross_Section=pow(0.75,F_Level)*Box_L*Box_W/2.0;  //截面
//
////三角形单元的参数
//struct element
//{
//	int level;
//	double v1[2];
//	double v2[2];
//	double v3[2];
//}ele;
//
//using namespace std;
//vector <element> tri;
//vector <element>::iterator tri_iter;
//
//double tri_pos[N_Fractal+1][4][2];  //三角形四个循环顶点坐标x,y

//	//判断高度范围
//	if (fabs(posz-Box_H/2.0)<Box_H/2.0)
//	{
//		int sp_flag=0;  //标志
//
//		for (int n=0; n<N_Fractal; n++)  //所有三角形
//		{
//			double sum_angle=0.0;
//			for (int m=0; m<3; m++)
//			{
//				double x1=tri_pos[n][m][0]-posx;
//				double x2=tri_pos[n][m+1][0]-posx;
//				double y1=tri_pos[n][m][1]-posy;
//				double y2=tri_pos[n][m+1][1]-posy;
//
//				double abs1=sqrt(x1*x1+y1*y1);
//				double abs2=sqrt(x2*x2+y2*y2);
//				double angle;
//
//				if (abs1>EPS && abs2>EPS)  //防止分母0
//				{
//					angle=(x1*x2+y1*y2)/abs1/abs2;
//				}
//				else  //距离某个顶点很近
//				{
//					sp_flag=1;
//					return 0;
//				}
//				//求内角和
//				sum_angle=sum_angle+acos(angle);
//			}
//
//			//在六角形内部
//			if (fabs(sum_angle-2*PI)<EPS)
//			{
//				sp_flag=1;
//				return 0;
//			}
//		}
//
//		if (sp_flag==0)  //不在空三角形
//		{
//			double esum_angle=0.0;
//			for (int p=0; p<3; p++)
//			{
//				double xx1=tri_pos[N_Fractal][p][0]-posx;
//				double xx2=tri_pos[N_Fractal][p+1][0]-posx;
//				double yy1=tri_pos[N_Fractal][p][1]-posy;
//				double yy2=tri_pos[N_Fractal][p+1][1]-posy;
//
//				double e_abs1=sqrt(xx1*xx1+yy1*yy1);
//				double e_abs2=sqrt(xx2*xx2+yy2*yy2);
//				double e_angle;
//
//				if (e_abs1>EPS && e_abs2>EPS)  //防止分母0
//				{
//					e_angle=(xx1*xx2+yy1*yy2)/e_abs1/e_abs2;
//				}
//				else  //距离某个顶点很近
//				{
//					return 2;
//				}
//				//求内角和
//				esum_angle=esum_angle+acos(e_angle);
//			}
//			//在大三角形内部
//			if (fabs(esum_angle-2*PI)<EPS)
//			{
//				return 2;
//			}
//			else
//			{
//				return 0;
//			}
//		}  //不在空三角形
//	}  //高度判断
//	else
//	{
//		return 0;
//	}

////fractal antennas
//void PreModel()
//{
////	using namespace std;
////
////	//中心
////	double tri_centerx=Box_L/2.0;
////	double tri_centery=Box_W/2.0;
////
////	//0层是等边三角形
////	ele.level=0;
////	ele.v1[0]=0;
////	ele.v1[1]=0;
////	ele.v2[0]=Box_L;
////	ele.v2[1]=0;
////	ele.v3[0]=tri_centerx;
////	ele.v3[1]=Box_W;
////	//存入数组
////	tri.push_back(ele);
////	//中间挖掉三角形数量
////	int empnum_tri=0;
////
////	//循环
////	do
////	{
////		element ele_r;
////		ele_r=tri.back(); //获得最后元素
////		tri.pop_back();   //弹出最后元素
////
////		if (ele_r.level!=3)
////		{
////			//分成三个非空三角形
////			ele.level=ele_r.level+1;
////			ele.v1[0]=ele_r.v1[0];
////			ele.v1[1]=ele_r.v1[1];
////			ele.v2[0]=0.5*(ele_r.v1[0]+ele_r.v2[0]);
////			ele.v2[1]=ele_r.v2[1];
////			ele.v3[0]=0.75*ele_r.v1[0]+0.25*ele_r.v2[0];
////			ele.v3[1]=0.5*(ele_r.v2[1]+ele_r.v3[1]);
////			tri.push_back(ele);
////
////			ele.level=ele_r.level+1;
////			ele.v1[0]=0.5*(ele_r.v1[0]+ele_r.v2[0]);
////			ele.v1[1]=ele_r.v1[1];
////			ele.v2[0]=ele_r.v2[0];
////			ele.v2[1]=ele_r.v2[1];
////			ele.v3[0]=0.25*ele_r.v1[0]+0.75*ele_r.v2[0];
////			ele.v3[1]=0.5*(ele_r.v2[1]+ele_r.v3[1]);
////			tri.push_back(ele);
////
////			ele.level=ele_r.level+1;
////			ele.v1[0]=0.75*ele_r.v1[0]+0.25*ele_r.v2[0];
////			ele.v1[1]=0.5*(ele_r.v1[1]+ele_r.v3[1]);
////			ele.v2[0]=0.75*ele_r.v2[0]+0.25*ele_r.v1[0];
////			ele.v2[1]=0.5*(ele_r.v2[1]+ele_r.v3[1]);
////			ele.v3[0]=ele_r.v3[0];
////			ele.v3[1]=ele_r.v3[1];
////			tri.push_back(ele);
////
////			//中心挖掉三角形坐标
////			tri_pos[empnum_tri][0][0]=0.5*(ele_r.v1[0]+ele_r.v2[0]);
////			tri_pos[empnum_tri][0][1]=0.5*(ele_r.v1[1]+ele_r.v2[1]);
////			tri_pos[empnum_tri][1][0]=0.5*(ele_r.v2[0]+ele_r.v3[0]);
////			tri_pos[empnum_tri][1][1]=0.5*(ele_r.v2[1]+ele_r.v3[1]);
////			tri_pos[empnum_tri][2][0]=0.5*(ele_r.v3[0]+ele_r.v1[0]);
////			tri_pos[empnum_tri][2][1]=0.5*(ele_r.v3[1]+ele_r.v1[1]);
////			tri_pos[empnum_tri][3][0]=tri_pos[empnum_tri][0][0];
////			tri_pos[empnum_tri][3][1]=tri_pos[empnum_tri][0][1];
////
////			empnum_tri++;  //三角形个数增加
////		}
////	} while (empnum_tri!=N_Fractal);
////
////	//最后一个大三角形
////	tri_pos[N_Fractal][0][0]=0;
////	tri_pos[N_Fractal][0][1]=0;
////	tri_pos[N_Fractal][1][0]=Box_L;
////	tri_pos[N_Fractal][1][1]=0;
////	tri_pos[N_Fractal][2][0]=tri_centerx;
////	tri_pos[N_Fractal][2][1]=Box_W;
////	tri_pos[N_Fractal][3][0]=0;
////	tri_pos[N_Fractal][3][1]=0;
//}