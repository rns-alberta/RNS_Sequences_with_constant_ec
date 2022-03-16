/************************************************************************** 
*                         NSSS.c
* 
* Computes constant neutron star sequences with constant baryonic mass
* 
**************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h> 

#include "consts.h"
#include "struct.h"

#include "nrutil.h"
#include "equil.h"
#include "equil_util.h"
#include "findmodel.h"
#include "surface.h"
#include "stableorbit.h"
#include "interpol.h"


/* Main; where it all starts and ends */

int main(int argc, char **argv)     /* Number of command line arguments, Command line arguments */
{ NeutronStar star;
  EOS eos;
  NeutronStar star1;
  EOS eos1;  
  int i, ierr;
  double
    e_min, e_max, e_max_mass,
    e_center=1e15,                     /* central en. density */
    B,                            /* Quark Bag Constant */
    K=3.0,                        /* Second parameter in "quark" eos */
    spin_freq=100,                  /* Spin Frequency */
    Gamma_P=0.0;                      /* Gamma for polytropic EOS */  
                
  int j;

  int a = 0, numseq=2;
  int spin_lim = 0;
  float e_c[4], M_0[4];
  float M0, Mtot,Mstat, Rstat, Radius, freq,freqK,energy_value, temp_energy, ratio_r = 1.0, ej,Volp;
  float maxmass, maxradius;   // Mass and radius of the maximum mass neutron star for an EOS
  float T, W;
  long int angmom;
  //double Kfreq, Kfreq_j;

  //andreas
  double lumi;
  double u_lumi;
  double poten[4];
  int call;
  double ratio_ch=0.005;
  double energy_min=0.2;//10^{15}$ g/$cm^3
  double e_ch=0.05;//10^{15}$ g/$cm^3
  FILE *out;

  //andreas

  FILE *fpointer;

  char eos_file[80] = "no EOS file specified";   /* EOS file name */
  char eos_type[80] = "tab";                     /* EOS type (poly or tab) */
  char data_dir[80] = "junk";                    /* Data output directory */
  char filename[100] = "NS_data_";


  //andreas

  /* READ IN THE COMMAND LINE OPTIONS */
  for(i=1;i<argc;i++) 
    if(argv[i][0]=='-'){
      switch(argv[i][1]){

      case 'f':
	/* IF A TABULATED EOS WAS CHOSEN, CHOOSE THE
	   NAME OF THE FILE */
	sscanf(argv[i+1],"%s",eos_file);
	break;
	
      case 'e':
	/* CHOOSE THE CENTRAL ENERGY DENSITY OF THE 
	   NEUTRON STAR (IN g/cm^3) */
	sscanf(argv[i+1],"%lf",&e_min);
	if(strcmp(eos_type,"poly")!=0)
	  e_min *= C*C*KSCALE;
	break;

      }
    }
 
  strncat(filename, eos_file, 65);
  strncat(filename, ".txt", 65);
  printf("%s\n", filename);
  fpointer = fopen(filename, "a");


  strncat(filename, eos_file, 65);
  strncat(filename, "_table.txt", 65);
  printf("%s\n", filename);
  out = fopen(filename, "a"); 

  /* PRINT THE HEADER */
  if(strcmp(eos_type,"tab")==0)
    printf("%s,  MDIVxSDIV=%dx%d\n",eos_file,MDIV,SDIV);
  if(strcmp(eos_type,"quark")==0)
    printf("Quark star with B=%f, MDIVxSDIV=%dx%d\n",B/1.602e33/KSCALE,MDIV,SDIV);

  /* SetUpStar loads in the eos and sets up the grid */
  /* Source code for SetUpStar can be found in findmodel.c */

  ierr = SetUpStar(eos_file, eos_type, data_dir, Gamma_P, B, K,
		   &eos, &star);



  //printf("The star infrastructure has been set up! \n");

  e_center = e_min;
  temp_energy = e_center;

  // Computing the star with the maximum mass and its corresponding radius
  ierr = MakeSphere(&eos, &star, e_max_mass);
  rns(1.0, e_max_mass, &eos, &star); 
  maxmass = star.Mass/MSUN;
  maxradius = star.R_e*1e-5;

 
  // Computing one neutron star
  if(1){

    
    ratio_r=1.00;
    double MaxM;
    //Computing the non-rotating spherical neutron star
    ierr = SetUpStar(eos_file, eos_type, data_dir, Gamma_P, B, K,&eos, &star);
    ierr = MakeSphere(&eos, &star, e_center); 
    rns(1.00, e_center, &eos, &star); 
    MaxM = star.Mass/MSUN;
    

    
    while(e_center>energy_min){

      //Computing the non-rotating spherical neutron star
      ierr = SetUpStar(eos_file, eos_type, data_dir, Gamma_P, B, K,&eos, &star);
      ierr = MakeSphere(&eos, &star, e_center); 
      rns(1.00, e_center, &eos, &star); 
      Mstat = star.Mass/MSUN;
      Rstat = star.R_e*1e-5;   

      // Computing the one star with given value of ratio_r

      while(star.Omega/(2.0*PI)<star.Omega_K/(2.0*PI)){

	// Compute a spherical star with the correct central energy density
	//ierr = MakeSphere(&eos, &star, e_center);
	//Mstat = star.Mass/MSUN;
	//Rstat = star.R_e*1e-5;
	// Computing the star with the maximum mass and its corresponding radius

        /*
	if(ratio_r < 0.7){
	  rns(0.7, e_center, &eos, &star);
	}*/
	rns(ratio_r, e_center, &eos, &star); 


	//andreas
	call = Surface(&eos,&star,&poten[0],&poten[1],&poten[2],&poten[3]);
	//fprintf(fpointer, "%7g %7g %7g %8g %4g %6g %6g %g %g %g\n", 
	//star.e_center, star.Mass/MSUN, star.Mass_0/MSUN,  star.R_e*1e-5, ratio_r, star.Omega/(2.0*PI), star.Omega_K/(2.0*PI), star.ang_mom, T, W);
	T = (0.5*star.ang_mom*star.Omega)/(C*C);
	W = star.Mp + T - star.Mass;

	Mtot=star.Mass/MSUN;
	M0=star.Mass_0/MSUN;
	Radius=star.R_e*1e-5;
	freq=star.Omega/(2.0*PI);
	freqK=star.Omega_K/(2.0*PI);
	angmom=star.ang_mom;
	Volp=star.Vp;
	//printf(" %g\n",star.ang_mom);
         
/*
	double ratioloop=ratio_r;
	double angmomnext,ecnext;
	if(Mtot>MaxM+0.0001){
          temp_energy = e_center;
          eos1=eos;
          star1=star;	  
	  while(1){   
	    ratioloop = ratioloop - 0.001;

	    for(j=0;j<3;j++){
	      ej = temp_energy - 0.001*j;

	      ierr = MakeSphere(&eos1, &star1, ej);
	      //rns(ratio_r, ej, &eos1, &star1); 
	      if(ratioloop < 0.7)
		rns(0.7, ej, &eos1, &star1);
  
	      rns(ratioloop, ej, &eos1, &star1); 

	      e_c[j] = ej;
	      M_0[j] = star1.Mass_0/MSUN;
	    }

	    energy_value = polyinter(M0, e_c, M_0);

	    temp_energy = energy_value;
	    ierr = MakeSphere(&eos1, &star1, energy_value);

	    if(ratioloop < 0.7)
	      rns(0.7, ej, &eos1, &star1);
    
	    rns(ratioloop, ej, &eos1, &star1); 

	    if( (star1.Omega/(2.0*PI)) > star1.Omega_K/(2.0*PI)){
             //angmomnext=0; 
             break;
	    }
	    if(isnan(star1.Mass/MSUN)){
	      printf("Mass is NAN\n");
              //angmomnext=0; 
	      break;
	    }

	    if((round(M0*100.0)/100.0) == (round(star1.Mass_0/MSUN * 100.0)/100.0)){
	      angmomnext=star1.ang_mom;
	      ecnext=ej;
	      break;
	    }


	  }
	}*/
	//printf(" %g\n",star1.ang_mom);


          printf("%g %.5f  %.5f  %.5f %.5f %.3f %.5f %.3f %.5f\n",
            star.e_center, star.Mass/MSUN, star.Mass_0/MSUN, Mstat, star.R_e*1e-5, ratio_r, Rstat, star.Omega/(2.0*PI), star.Omega_K/(2.0*PI));
           fprintf(fpointer, "%7g %7g %7g %6g %8g %4g %6g %7g %6g  %g  %g  %g %g\n", 
             star.e_center, star.Mass/MSUN, star.Mass_0/MSUN, Mstat, star.R_e*1e-5, ratio_r, Rstat, star.Omega/(2.0*PI), star.Omega_K/(2.0*PI), star.ang_mom, T/MSUN, W/MSUN,star.r_surf[MDIV]/star.r_surf[1]);
             //fprintf(fpointer, "%7g %7g %7g %6g %f %8g %4g %6g %7g %6g  %g  %g  %g %f %f %f\n", 
             //star.e_center, star.Mass/MSUN, star.Mass_0/MSUN, Mstat, 0.0, star.R_e*1e-5, ratio_r, Rstat, star.Omega/(2.0*PI), star.Omega_K/(2.0*PI), star.ang_mom, T,W, 0.0,0.0,0.0 );
	if(ratio_r==1.00){
             fprintf(out, "%s %g %s %g %s %g %s \n","&",star.e_center,"&",star.Mass/MSUN,"&",star.R_e*1e-5,"\\\\");
	}
	/*if(ratio_r>0.995){
	  ratio_ch=0.0001;
	}
	else{
	  ratio_ch=0.01;  
	}*/
	ratio_r-=ratio_ch; 
      }
      //break;
      ratio_r=1.00;
      e_center-=e_ch;
    }

    //andreas
  }


  fclose(out);//andreas
  fclose(fpointer);
  return 0;
}









