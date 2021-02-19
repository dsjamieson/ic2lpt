#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "real.h"
#include <omp.h>

#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))
#define pow4(x) ((x)*(x)*(x)*(x))

#define myreal float

FILE *fd;
#pragma omp threadprivate(fd)

long totNpart;
long myNpart;
double g_boxsize;
double g_Omega_m;
double g_Omega_l;
double g_redshift;
double a_exp;
double g_hubble;
double part_mass;
int g_files;


struct io_header_1
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header1;



void ShowHeader()
{
  int i=0, n;
  for(i=0;i<6;i++)
    {
      header1.npart[i] =(header1.npart[i]);
      header1.npartTotal[i] = (header1.npartTotal[i]);
      header1.mass[i]=(header1.mass[i]);
    }

  header1.BoxSize=(header1.BoxSize);
}



void  write_num_bytes(FILE *f, int nbytes)
{
	if (sizeof(int)==4)
		fwrite(&nbytes,sizeof(int),1,f);
	else
	{
		fprintf(stderr,"size of int is NOT 4!!!");
		exit(1);
	}
}

void write_header(struct io_header_1 * h)
{
(*h).npart[0]=0;
(*h).npart[1]=myNpart;
(*h).npart[2]=0;
(*h).npart[3]=0;
(*h).npart[4]=0;
(*h).npart[5]=0;

(*h).npartTotal[0]=0;
(*h).npartTotal[1]=totNpart;
(*h).npartTotal[2]=0;
(*h).npartTotal[3]=0;
(*h).npartTotal[4]=0;
(*h).npartTotal[5]=0;

(*h).mass[0]=0;
(*h).mass[1]=part_mass;
(*h).mass[2]=0;
(*h).mass[3]=0;
(*h).mass[4]=0;
(*h).mass[5]=0;
	
(*h).redshift=g_redshift;
(*h).time=1./(1.+(*h).redshift);
(*h).num_files=g_files;
(*h).BoxSize=g_boxsize;
(*h).Omega0=g_Omega_m;
(*h).OmegaLambda=g_Omega_l;
(*h).HubbleParam=g_hubble;


}


void write_snapshot(char *fname, int file){
	      char   *buf;
	      myreal buffer,sqrta;
          long i;

	      buf   = (char *)  calloc(256,   sizeof(char));
	
	      if(g_files>1)
		  sprintf(buf,"gadget_%s.%d", fname, file);
	      else
		  sprintf(buf,"gadget_%s.dat", fname);
	      fd=fopen(buf,"w");
	      if(fd==NULL)
	      {
		  printf("can't create file `%s`\n",buf);
		  exit(0);
	      }
	     
	      write_header(&header1);
	      write_num_bytes(fd,sizeof(header1));
	      fwrite(&header1, sizeof(struct io_header_1), 1, fd);
	      write_num_bytes(fd,sizeof(header1));
	      
	      ShowHeader(); 
	      write_num_bytes(fd,myNpart*3*sizeof(myreal));
	
	      for(i=0;i<(header1).npart[1];i++){ 
		  buffer=(myreal) xx[i+myNpart*file]*1000.;
		  fwrite(&buffer, sizeof(myreal),1 , fd);
		  buffer=(myreal) yy[i+myNpart*file]*1000.;
		  fwrite(&buffer, sizeof(myreal),1 , fd);
		  buffer=(myreal) zz[i+myNpart*file]*1000.;
		  fwrite(&buffer, sizeof(myreal),1 , fd);
	      }
	      write_num_bytes(fd,myNpart*3*sizeof(myreal));
	      
	      write_num_bytes(fd,myNpart*3*sizeof(myreal));
      
	      sqrta=sqrt((header1).time);
	      for(i=0;i<(header1).npart[1];i++){ 
		  buffer=vx[i+myNpart*file]/sqrta/g_hubble;
		  fwrite(&buffer, sizeof(myreal),1 , fd);
		  buffer=vy[i+myNpart*file]/sqrta/g_hubble;
		  fwrite(&buffer, sizeof(myreal),1 , fd);
		  buffer=vz[i+myNpart*file]/sqrta/g_hubble;
		  fwrite(&buffer, sizeof(myreal),1 , fd);
	      }
	      write_num_bytes(fd,myNpart*3*sizeof(myreal));
	      
	      write_num_bytes(fd,myNpart*sizeof(unsigned int));    
	      unsigned int u;
	      for(i=0;i<(header1).npart[1];i++){ 
		  u=i+file*myNpart;
		  fwrite(&u, sizeof(unsigned int),1 , fd);
	      }
	      write_num_bytes(fd,myNpart*sizeof(unsigned int));
	       
	    fclose(fd);

	return;
}


void write_gadget(int np,double b,double z_i, double h,double omega_m, double omega_lambda, char*outprefix,int no_files){

    int file;
    g_Omega_m=omega_m;
    g_Omega_l=omega_lambda;
    g_hubble=h;
    g_redshift=z_i;
//    totNpart=pow3(np);
    totNpart=((long)np)*np*np;
    g_boxsize=b*1000.;
    g_files=no_files;


    double G=6.672e-8;
    double SOLAR_MASS=1.989e33;
    double H0=3.2407789e-18;
    double PI=3.14159265358979323846;
    double CM_PER_MPC = 3.085678e24;

    double  rhocrit = 3.*H0*H0/(8.*PI*G); // in gr/cm^3
    rhocrit = pow3(CM_PER_MPC)*rhocrit/ SOLAR_MASS; // en Msun/Mpc^3

    part_mass = g_Omega_m*rhocrit/1.e10*pow3(g_boxsize)/totNpart/pow3(1000.); // 27.7424215e10

    a_exp=1./(1.+g_redshift);


    if (totNpart%no_files!=0){
	printf("totNpart is not divideable by no_files!=0 !!!\n");
	exit(1);
    } 

    myNpart=totNpart/no_files;
#pragma omp parallel private(file)
 { 
#pragma omp for 
    for (file=0;file<no_files; file++){
	write_snapshot(outprefix,file);
    }
 } 

}
