/************************************************************************** 
*                         SpinUp.c
* 
* Computes zero spin stars
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
#include "interpol.h"
#include "spinup_func.h"
#include "stableorbit.h"
#include "surface.h"

/* Main; where it all starts and ends */

int main(int argc, char **argv)     /* Number of command line arguments, Command line arguments */
{ NeutronStar star;
  EOS eos;    
  int i;
  int ierr;
  double
  r_min=0.8, r_max=1.0, e_max_mass,
  e_center=1e15,                   /* central en. density */
  B,                               /* Quark Bag Constant */
  K=3.0,                           /* Second parameter in "quark" eos */
  spin_freq=100,                   /* Spin Frequency */
  Gamma_P=0.0;                     /* Gamma for polytropic EOS */  

  int spin_lim = 0;
  double Mstat, Rstat, ratio_r = 1.0;
  float maxmass, maxradius;   // Mass and radius of the maximum mass neutron star for an EOS
  float T, W;

  FILE *fpointer;

  char eos_file[80] = "no EOS file specified";   /* EOS file name */
  char eos_type[80] = "tab";                     /* EOS type (poly or tab) */
  char data_dir[80] = "junk";                    /* Data output directory */
  char filename[100] = "./";

  //New variables
  int numstar=2; /* Number of stars computed */
  double rstep;

  /* READ IN THE COMMAND LINE OPTIONS */
  for(i=1;i<argc;i++) 
    if(argv[i][0]=='-'){
      switch(argv[i][1]){

      case 'q':
	/* CHOOSE THE EOS TYPE: EITHER "tab" or "poly" or "quark"
	   (default is tab) */
	sscanf(argv[i+1],"%s",eos_type);
	break;  

      case 'b':
	sscanf(argv[i+1],"%lf",&B);
	B *= 1.602e33*KSCALE;
	break;       

      case 'f':
	/* IF A TABULATED EOS WAS CHOSEN, CHOOSE THE
	   NAME OF THE FILE */
	sscanf(argv[i+1],"%s",eos_file);
	break;

      case 'd':
	/* CHOOSE THE NAME OF THE OUTPUT DIRECTORY */
	sscanf(argv[i+1],"%s",data_dir);
	break;

      case 'e':
	/* CHOOSE THE CENTRAL ENERGY DENSITY OF THE 
	   NEUTRON STAR (IN g/cm^3) */
	sscanf(argv[i+1],"%lf",&e_center);
	if(strcmp(eos_type,"poly")!=0)
	  e_center *= C*C*KSCALE;
	break;

      case 'n':  
	/* CHOOSE THE NUMBER OF SEQUENCES 
	   PRODUCED */
	sscanf(argv[i+1],"%d",&numstar);
	break;

      case 'm':  
	/* CHOOSE CENTRAL ENERGY DENSITY FOR THE
	   MAXIMUM MASS NEUTRON STAR FOR THE EOS */
	sscanf(argv[i+1],"%lf",&e_max_mass);
	break;

      case 'o':
	/* Name of Output file */
	sscanf(argv[i+1],"%s",filename);
	break;

     case 's':
	/* CHOOSE THE SPIN FREQUENCY (HZ) */
	sscanf(argv[i+1],"%lf",&spin_freq);
	printf("spin=%g\n",spin_freq);
	break;

     case 'r':
        /* CHOOSE THE RATIO OF r_p and r_e */
        sscanf(argv[i+1],"%lf",&ratio_r);
        break;

     case 't':
        /* CHOOSE TO COMPUTE SEQUENCE UP TO 
	  A CERTAIN VALUE */
        sscanf(argv[i+1],"%d",&spin_lim);
        break;

     case 'h': 
	fprintf(stderr,"\nQuick help:\n\n");
	fprintf(stderr,"  -q EOS type (tab)\n"); 
	fprintf(stderr,"     tab   : tabulated \n");
        fprintf(stderr,"     quark : simple quark model \n"); 
	fprintf(stderr,"  -b bag constant in MeV/fm^3 for quark models\n");
	fprintf(stderr,"  -f EOS file \n");
	fprintf(stderr,"  -d directory output goes to \n");
	fprintf(stderr,"  -e lowest central energy density to be used, in gr/cm^3\n");
	fprintf(stderr,"  -n number of stars [2] \n");
	fprintf(stderr,"  -h this menu\n\n");
	exit(1);
	break;  
      }
    }
 
  printf("opening for output %s\n", filename);
  fpointer = fopen(filename, "a");



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
  
  //Opens files for creating plots
  FILE *fpointer_plot;
  fpointer_plot = fopen("plot_Mvsloge.txt", "w");
  FILE *fpointer_Macc;
  fpointer_Macc = fopen("plot_Macc.txt", "w");

  //Set the goal values of M and Omega
  MakeSphere(&eos, &star, e_center);
  rns(ratio_r, e_center, &eos, &star);

  //Sets initial values for star mass and spin
  double M_star = 2.35, Omega_star = 700;
  printf("\nInitial Mass = %f Msun, Spin = %6.5e Hz\n", M_star, Omega_star);
 
  double *M_nospin, *r_nospin;
  M_nospin = dvector (0, numstar-1);
  r_nospin = dvector (0, numstar-1);
  int stop = 0;

  //Initiate step size
  if (numstar == 1) {
    rstep = r_max;
  }
  else { 
    rstep = (r_max-r_min)/ (numstar-1.0);
  }
 
  ratio_r = r_max;
  
  printf("Searching for initial mass \n");

  ierr = MakeSphere(&eos, &star, e_center);

  loop_M(&eos, &star, ratio_r, e_center, numstar, Mstat, Rstat, M_nospin, r_nospin, rstep, M_star, stop, fpointer_plot);



  /*Find and interpolate value of M */
      
  double r_pt, start[2];
  int count;
  r_pt = find_M(&eos, &star, e_center, M_star, M_nospin, r_nospin, numstar, 
		Mstat, Rstat, fpointer_plot);

  printf("Initial mass has been found \n \n");

  /*Find and interpolate value of Spin */
  float temp_energy;
  double *Omega_spin, *e_spin, *r_spin;
  Omega_spin = dvector (0, 100);
  e_spin = dvector (0, 100);
  r_spin = dvector (0, 100);
  stop = 0;
  temp_energy = e_center;

  // Computing the star with the maximum mass and its corresponding radius
  ierr = MakeSphere(&eos, &star, e_max_mass);
  rns(1.0, e_max_mass, &eos, &star); 
  maxmass = star.Mass/MSUN;
  maxradius = star.R_e*1e-5;

  // Computing sequences with the same M0
  printf("Searching for initial spin \n");
  count = loop_spin(&eos, &star, ratio_r, e_center, temp_energy, maxmass, Omega_star, Omega_spin, r_spin, e_spin, r_pt, Mstat, Rstat, T, W, stop, fpointer_plot);

  //Find and interpolate value of M0
  find_spin(&eos, &star, ratio_r, count, numstar, M_star, Omega_star, Omega_spin, e_spin, r_spin, Mstat, Rstat, fpointer_plot, start);

  e_center = start[0];
  ratio_r = start[1];
  printf("Initial spin has been found \n");

  printf("Starting point found \n \n");

  Surface(&eos, &star);
  orbit(&eos, &star);
  //printf("orbit freq at isco = %f \n", star.orbitP);

  double Spin_final;

  printf("Searching for the equilibrium radius where the orbital frequency = spin frequency\n");
  double r_equil;
  r_equil = calc_omega(&eos, &star);
  
  double M0tot = 0, M0prev = star.Mass_0/MSUN, M0acc, Mtot = 0, Mprev = star.Mass/MSUN, Macc;

  printf("----------------------------------\n");

  while (ratio_r < 0.998) {

  //Find goal values of M0 and J
  MakeSphere(&eos, &star, e_center);
  rns(ratio_r, e_center, &eos, &star);

  double l;
  double delta_m = 1e-4;
  double M0_star, J_star;
  int check_step=0;

  l = find_l(&star, r_equil);
  printf("l = %6.5e \n", l);

  J_star = star.ang_mom - (delta_m*MSUN * l);  //cgs 

  //Option to make delta_m smaller close to spin = 0
  if (J_star < 3e48) {
    check_step += 1;
  }
  if (check_step ==1) {
    delta_m *= 1e-1;
    printf("stepsize too big\n");
  }

  M0_star = star.Mass_0/MSUN - delta_m;  //dimensionless
  J_star = star.ang_mom - (delta_m*MSUN * l);  //cgs 
  //M0_star = 3.2;
  //J_star = 2e49;
  printf("M0_star = %f, J_star = %6.5e\n", M0_star, J_star);
 
  printf("Searching for M0_star\n");
  double *M0_nospin, *r_nospin;
  M0_nospin = dvector (0, numstar-1);
  r_nospin = dvector (0, numstar-1);
  int stop = 0, direction = 0;
  //If direction = 0, spin is decreasing. If direction = 1, spin is increasing.

  //Initiate step size
  if (numstar == 1) {
    rstep = r_max;
  }
  else { 
    rstep = (r_max-r_min)/ (numstar-1.0);
  }
 
  if (direction ==1) ratio_r = r_max;

  ierr = MakeSphere(&eos, &star, e_center);

  loop_M0(&eos, &star, ratio_r, e_center, numstar, Mstat, Rstat, M0_nospin, r_nospin, rstep, M0_star, stop, fpointer_plot, direction);

    /***********************************************************************************/

  /*Interpolate value of M0 */
      
  double r_pt;
  int count;
  r_pt = find_M0(&eos, &star, e_center, M0_star, M0_nospin, r_nospin, numstar, Mstat, Rstat, fpointer_plot);
 
  printf("M0_star has been found \n \n");
  printf("Searching for J_star \n");

  float temp_energy;
  double *J_spin, *e_spin, *r_spin;
  J_spin = dvector (0, 1e4);
  e_spin = dvector (0, 1e4);
  r_spin = dvector (0, 1e4);
  stop = 0;
  temp_energy = e_center;

  // Computing the star with the maximum mass and its corresponding radius
  ierr = MakeSphere(&eos, &star, e_max_mass);
  rns(1.0, e_max_mass, &eos, &star); 
  maxmass = star.Mass/MSUN;
  maxradius = star.R_e*1e-5;

  // Computing sequences with the same M0
  count = loop_J(&eos, &star, ratio_r, e_center, temp_energy, maxmass, J_star, J_spin, r_spin, e_spin, r_pt, Mstat, Rstat, T, W, stop, fpointer_plot);

  /*Interpolate value of M0 */
  r_pt = find_J(&eos, &star, ratio_r, count, numstar, M0_star, J_star, J_spin, e_spin, r_spin, Mstat, Rstat, fpointer_plot);

  M0acc = M0prev - star.Mass_0/MSUN;
  M0tot += M0acc;
  Macc = Mprev - star.Mass/MSUN;
  Mtot += Macc;
  fprintf(fpointer_Macc, "%f, %f, %f\n", 
	  star.Omega/(2.0*PI), Mtot, Macc, M0tot, M0acc);
  M0prev = star.Mass_0/MSUN;
  Mprev = star.Mass/MSUN;
  printf("Total mass accreted = %f, Mass accreted in this step = %f\n", Mtot, Macc);
  printf("Total baryonic mass accreted = %f, Baryonic mass accreted in this step = %f\n", M0tot, M0acc);

  e_center = star.e_center;
  ratio_r = r_pt;
  Spin_final = star.Omega/(2.0*PI);
  //printf("e_center = %f\n", e_center);
  //printf("ratio_r = %f\n", ratio_r);
  printf("Spin frequency = %f\n", Spin_final);
  
  Surface(&eos, &star);
  orbit(&eos, &star);
  printf("Orbit frequency at isco = %f \n", star.orbitP);
  if (star.orbitP < Spin_final) {
    printf("orbital freq < spin freq/n");
    break;
  }
  
  printf("----------------------------------------------\n");
  }

  fclose(fpointer_plot);
  fclose(fpointer_Macc);

  fclose(fpointer); 
  return 0;
}









