#include "mpi.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
using namespace std;

int main(int argc, char* argv[])
{
    char string[255];
    char s[255];

    int m;
    int myid,numprocs;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);

    remove("../result/Ene_Metal.txt");  //删除能量文件
	remove("../result/Grid.txt");       //删除网格文件
	remove("../result/Ene_Active.txt"); //删除能量文件

	MPI_Barrier(MPI_COMM_WORLD);

	for (m=1; m<=41; m++)  //频率
	{
		if (m%numprocs==myid)
		{
            strcpy(string, "./CG_FFT_VIE ");
            sprintf(s,"%d\n",m);
            strcat(string,s);
            system(string);
		}
	}

    MPI_Barrier(MPI_COMM_WORLD);

	char string_fre[256];      //扫频点
    char string_metal[256];    //金属吸收消失截面文件路径
    char string_active[256];   //有机材料吸收能量文件路径
    char string_grid[256];     //网格剖分文件路径
    //文件指针
	FILE* file0;
	FILE* file1;
	FILE* file2;
	FILE* file3;
	FILE* file4;
	FILE* file5;

    //组合删除文件
    if (myid==1)  //标号1处理器
    {
        //最终拼合文件
        file3=fopen("../result/Ene_Metal.txt","w");
        file4=fopen("../result/Ene_Active.txt","w");
        file5=fopen("../result/Grid.txt","w");

        for (m=0; m<41; m++)
        {
            //读入数据文件
            strcpy(string_metal, "../result/Ene_Metal");
            strcpy(string_active, "../result/Ene_Active");
            strcpy(string_grid, "../result/Grid");
            sprintf(string_fre,"%d",m);
            strcat(string_metal,string_fre);
            strcat(string_active,string_fre);
            strcat(string_grid,string_fre);
            strcat(string_metal,".txt");
            strcat(string_active,".txt");
            strcat(string_grid,".txt");

            file0=fopen(string_metal,"r");
            file1=fopen(string_active,"r");
            file2=fopen(string_grid,"r");

            int wavelength1,wavelength2;
            double extinct;
            double absorb1,absorb2;
            int Fre_Point,Nx,Ny,Nz;
            double delta_x,delta_y,delta_z;

            //读入金属吸收截面，消失截面
            fscanf(file0,"%d %lf %lf\n",&wavelength1,&extinct,&absorb1);
            //读入吸收材料的单位体积的吸收能量
            fscanf(file1,"%d %lf\n",&wavelength2,&absorb2);
            //读入建模的信息
            fscanf(file2,"%d %d %d %d %lf %lf %lf\n",&Fre_Point,&Nx,&Ny,&Nz,&delta_x,&delta_y,&delta_z);

            //写入金属吸收截面，消失截面
            fprintf(file3,"%d %.15lf %.15lf\n",wavelength1,extinct,absorb1);
            //写入吸收材料的单位体积的吸收能量
            fprintf(file4,"%d %.15lf\n",wavelength2,absorb2);
            //写入建模的信息
            fprintf(file5, "%d %d %d %d %.15lf %.15lf %.15lf\n",Fre_Point,Nx,Ny,Nz,delta_x,delta_y,delta_z);

            //关闭删除读入文件
            fclose(file0);
            fclose(file1);
            fclose(file2);
            remove(string_metal);
            remove(string_active);
            remove(string_grid);
        }
        //关闭写入文件
        fclose(file3);
        fclose(file4);
        fclose(file5);
    }

    MPI_Finalize();

    cout << "Finished" << endl;
    return 0;
}
