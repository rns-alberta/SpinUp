/************************************************************************** 
*                         SPINUP_func.c
* 
* Computes zero spin stars
* 
**************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "spinup_func.h"
#include "struct.h"
#include "consts.h"
#include "findmodel.h"
#include "equil_util.h"
#include "nrutil.h"
#include "equil.h"
#include "interpol.h"

void loop_M(EOS *eos, NeutronStar *star, double ratio_r, double e_center, int numstar, double Mstat, double Rstat, double *M_nospin, double *r_nospin, double rstep, double M_star, int stop, FILE *fpointer_plot) {
  int k = 0;

  printf("e_c \tMass \t Mass_0\t  StatM\t  Radius   R-ratio StatR   Spin\t       K freq      Ang mom\n");
  printf("e15 \tMsun \t Msun\t  Msun\t  km\t   --  \t   km \t   Hz \t       Hz \t   g cm^2 /s \n");

  for (k=0; k<=numstar-1; k++) {
    rns(ratio_r, e_center, eos, star);

    printf("%.5f %.5f  %.5f  %.5f %.5f %.5f %.5f %.5e %.5e %.5e\n",
	   star->e_center, star->Mass/MSUN, star->Mass_0/MSUN, Mstat, 
	   star->R_e*1e-5, ratio_r, Rstat, star->Omega/(2.0*PI), 
	   star->Omega_K/(2.0*PI), star->ang_mom);
    /*fprintf(fpointer, "%7g %7g %7g %6g %4g %8g %4g %6g %7g %6g  %g  %g  %g  %g  %g\n", 
	     star->e_center, star->Mass/MSUN, star->Mass_0/MSUN, Mstat, maxmass, 
	      star->R_e*1e-5, ratio_r, Rstat, star->Omega/(2.0*PI), 
	      star->Omega_K/(2.0*PI), star->ang_mom, T, W, maxradius, 
	      maxmass/maxradius);*/
    
    fprintf(fpointer_plot, "%f, %f \n", 
             star->e_center, star->Mass/MSUN); 

    M_nospin[k] = star->Mass/MSUN;
    r_nospin[k] = ratio_r;
    ratio_r -= rstep;
    if(star->Mass/MSUN > M_star)
      stop += 1;
    if(stop>1 && k>2)
      break;
    }
    }

double find_M(EOS *eos, NeutronStar *star, double e_center, double M_star, double *M_nospin, double *r_nospin, int numstar, double Mstat, double Rstat, FILE *fpointer_plot) {
  double dif, dift, r_pt;
  int j, ns, new;
  dif = M_star - M_nospin[0];

  for (j=0; j<numstar; j++) {
    dift = M_star - M_nospin[j];
    //printf("dif = %f, dift = %f \n", dif, dift);
    if (dift*dif <0) {
      ns = j;
      break;
    }
    dif = dift;
  }
  ns -=1;
  //printf("ns = %d \n", ns);

  printf("Finding point immediately before M0 = %f \ninitial point: M0 = %f, r = %f \n", M_star, M_nospin[ns], r_nospin[ns]);

  new=IMIN(IMAX(ns - 1,0),numstar-4);
  r_pt = interpolate(&M_nospin[new], &r_nospin[new], M_star);
  //printf("r_pt = %f \n", r_pt);
  rns(r_pt, e_center, eos, star);
  printf("The interpolated point has \n%.5f %.5f  %.5f  %.5f %.5f %.5f %.5f %.5e %.5e %.5e \n",
	   star->e_center, star->Mass/MSUN, star->Mass_0/MSUN, Mstat, 
	   star->R_e*1e-5, r_pt, Rstat, star->Omega/(2.0*PI), 
	   star->Omega_K/(2.0*PI), star->ang_mom);
  fprintf(fpointer_plot, "%f, %f \n", 
           star->e_center, star->Mass/MSUN);

  printf("The percent error is %6.5e \n", (M_star - star->Mass/MSUN)/M_star*100);
  return r_pt;  
  }

double loop_spin (EOS *eos, NeutronStar *star, double ratio_r, double e_center, double temp_energy, double maxmass, double Omega_star, double *Omega_spin, double *r_spin, double *e_spin, double r_pt, double Mstat, double Rstat, double T, double W, int stop, FILE *fpointer_plot) {
  printf("e_c \tMass \t Mass_0\t  StatM\t  Radius   R-ratio StatR    Spin\tK freq      Ang mom\n");
  printf("e15 \tMsun \t Msun\t  Msun\t  km\t   --  \t   km \t    Hz\t        Hz \t    g cm^2 /s \n");

  int j, ierr;
  int count;
  double M, ej, energy_value;
  float e_c[4], M_0[4];
  count = 0;
  
  //printf("r_pt = %f\n", r_pt);
  ratio_r = r_pt;
  temp_energy = e_center;
  //printf("temp_energy = %g \n",temp_energy);
  
  ierr = MakeSphere(eos, star, e_center);
  rns(ratio_r, e_center, eos, star); 
    
  Mstat = star->Mass/MSUN;
  Rstat = star->R_e*1e-5;

  T = (0.5*star->ang_mom*star->Omega)/(C*C);
  W = star->Mp + T - star->Mass;
  
  printf("%.5f %.5f  %.5f  %.5f %.5f %.5f %.5f %.5e %.5e %.5e\n",
	   star->e_center, star->Mass/MSUN, star->Mass_0/MSUN, Mstat, 
	   star->R_e*1e-5, ratio_r, Rstat, star->Omega/(2.0*PI), 
	   star->Omega_K/(2.0*PI), star->ang_mom);

  /*fprintf(fpointer, "%7g %7g %7g %6g %4g %8g %4g %6g %7g %6g  %g  %g %g  %g  %g\n", star.e_center, star.Mass/MSUN, star.Mass_0/MSUN, Mstat, maxmass, star.R_e*1e-5, ratio_r, Rstat, star.Omega/(2.0*PI), star.Omega_K/(2.0*PI), star.ang_mom, T, W, maxradius, maxmass/maxradius);*/

  /*fprintf(fpointer_plot, "%f, %f \n", 
    star->e_center, star->Mass/MSUN);*/

  Omega_spin[count] = star->Omega/(2.0*PI);
  e_spin[count] = star->e_center;
  r_spin[count] = ratio_r;
  //printf("J = %.5e, e = %f, count = %d \n", J_spin[count], e_spin[count], count);
   
  count += 1;

  M = star->Mass/MSUN;
  
  //printf("Omega_star = %f, Omega = %f\n", Omega_star, star->Omega/(2.0*PI));
  if (Omega_star - star->Omega/(2.0*PI) > 0) {
    while(1){   
      //printf("---------------------------------------------------------------------------\n");
      ratio_r = ratio_r - 5e-3;

      //printf("ratio_r = %f\n", ratio_r);

      for(j=0;j<3;j++){

	ej = temp_energy - 1e-2*j;
	//printf("ej = %f\n", ej);

	ierr = MakeSphere(eos, star, ej);
	  
	if(ratio_r < 0.7)
	  rns(0.7, ej, eos, star);
  
	  rns(ratio_r, ej, eos, star); 

	  //printf("%g %.5f  %.5f  %.5f %.5f %.3f %.4f %.4f \n",
	  //    ej, star.Mass/MSUN, star.Mass_0/MSUN, Mstat, star.R_e*1e-5, ratio_r, Rstat, star.Omega/(2.0*PI));

	  e_c[j] = ej;
	  M_0[j] = star->Mass/MSUN;
	}
	//printf("-----------------------------------------------\n");
	//printf("M_0 = %g\n",M0);
	//printf("energies:  %g, %g, %g\n",e_c[0], e_c[1], e_c[2]);
	//printf("M_0's: %g, %g, %g\n",M_0[0], M_0[1], M_0[2]);

	energy_value = polyinter(M, e_c, M_0);

	//printf("---------------------------------------------------------------------------\n");
	temp_energy = energy_value;
	ierr = MakeSphere(eos, star, energy_value);
	//rns(ratio_r, energy_value, &eos, &star); 

	if(ratio_r < 0.7)
	  rns(0.7, ej, eos, star);
    
	rns(ratio_r, ej, eos, star); 

	printf("%.5f %.5f  %.5f  %.5f %.5f %.5f %.5f %.5e %.5e %.5e\n",
	   energy_value, star->Mass/MSUN, star->Mass_0/MSUN, Mstat, 
	   star->R_e*1e-5, ratio_r, Rstat, star->Omega/(2.0*PI), 
	   star->Omega_K/(2.0*PI), star->ang_mom);

	/*fprintf(fpointer_plot, "%f, %f \n", 
	  energy_value, star->Mass/MSUN);*/

	Omega_spin[count] = star->Omega/(2.0*PI);
	e_spin[count] = energy_value;
	r_spin[count] = ratio_r;
	//printf("J = %.5e, e = %f, count = %d \n", J_spin[count], e_spin[count], count);
	if(star->Omega/(2.0*PI) > Omega_star)
	  stop += 1;
	if(stop>0 && count>2)
	  break;
	count += 1;

	if( (star->Omega/(2.0*PI)) > (star->Omega_K/(2.0*PI))) {
	  printf("Max reached\n");
	  break;
	}
	if(isnan(star->Mass/MSUN)){
	  printf("Mass is NAN\n");
	  break;
	}
      }
      }

    else {
	while(1){   
	//printf("---------------------------------------------------------------------------\n");
	  ratio_r = ratio_r + 5e-3;

	//printf("ratio_r = %f\n", ratio_r);

	for(j=0;j<3;j++){
	  ej = temp_energy + 1e-2*j;
	  //printf("ej = %f\n", ej);

	  ierr = MakeSphere(eos, star, ej);
	   
	  if(ratio_r < 0.7)
	    rns(0.7, ej, eos, star);
  
	  rns(ratio_r, ej, eos, star); 

	  //printf("%g %.5f  %.5f  %.5f %.5f %.3f %.4f %.4f \n",
	  //ej, star->Mass/MSUN, star->Mass_0/MSUN, Mstat, star->R_e*1e-5, ratio_r, Rstat, star->Omega/(2.0*PI));

	  e_c[j] = ej;
	  M_0[j] = star->Mass/MSUN;
	  //printf("M_0 = %f\n", star->Mass/MSUN);
	}
	//printf("-----------------------------------------------\n");
	//printf("M_0 = %g\n",M0);
	//printf("energies:  %g, %g, %g\n",e_c[0], e_c[1], e_c[2]);
	//printf("M_0's: %g, %g, %g\n",M_0[0], M_0[1], M_0[2]);

	energy_value = polyinter(M, e_c, M_0);

	//printf("---------------------------------------------------------------------------\n");
	temp_energy = energy_value;
	//printf("new energy = %f\n", energy_value);
	
	ierr = MakeSphere(eos, star, energy_value);
	//rns(ratio_r, energy_value, &eos, &star); 


	if(ratio_r < 0.7)
	  rns(0.7, ej, eos, star);
    
	rns(ratio_r, ej, eos, star); 

	printf("%.5f %.5f  %.5f  %.5f %.5f %.5f %.5f %.5e %.5e %.5e\n",
	   energy_value, star->Mass/MSUN, star->Mass_0/MSUN, Mstat, 
	   star->R_e*1e-5, ratio_r, Rstat, star->Omega/(2.0*PI), 
	   star->Omega_K/(2.0*PI), star->ang_mom);

	/*fprintf(fpointer_plot, "%f, %f \n", 
	  energy_value, star->Mass/MSUN);*/

	Omega_spin[count] = star->Omega/(2.0*PI);
	e_spin[count] = energy_value;
	r_spin[count] = ratio_r;
	//printf("J = %.5e, e = %f, count = %d \n", J_spin[count], e_spin[count], count);

	if(star->Omega/(2.0*PI)< Omega_star)
	  stop += 1;
	if(stop>0 && count>2)
	  break;
	count += 1;

	if( (star->Omega/(2.0*PI)) > (star->Omega_K/(2.0*PI))) {
	  printf("Max reached\n");
	  break;
	}
	if(isnan(star->Mass/MSUN)){
	  printf("Mass is NAN\n");
	  break;
	}
      }

	T = (0.5*star->ang_mom*star->Omega)/(C*C);
	W = star->Mp + T - star->Mass;
	//printf("Mp = %g T = %g Mass = %g  J = %g W = %g T/W = %g \n", star.Mp, T, star.Mass, star.ang_mom, W, T/W);
    //printf("M0 = %g \t Mass_0 = %g\n", M0, star.Mass_0/MSUN);
	/*if((round(M0*100.0)/100.0) == (round(star->Mass_0/MSUN * 100.0)/100.0))
	  fprintf(fpointer, "%7g %7g %7g %6g %4g %8g %4g %6g %7g %6g  %g  %g  %g  %g  %g\n", energy_value, star->Mass/MSUN, star->Mass_0/MSUN, Mstat, maxmass, star-.R_e*1e-5, ratio_r, Rstat, star->Omega/(2.0*PI), star->Omega_K/(2.0*PI), star->ang_mom, T, W, maxradius, maxmass/maxradius);*/

      }
  return count;
}

void find_spin (EOS *eos, NeutronStar *star, double ratio_r, int count, int numstar, double M_star, double Omega_star, double *Omega_spin, double *e_spin, double *r_spin, double Mstat, double Rstat, FILE *fpointer_plot, double start[]) {
  count += 1;
  int ns = 0;

  printf("Finding point immediately before Spin = %.5e \n", Omega_star);

  double dif, dift, e_pt, r_pt;
  int j=0, new;
  dif = Omega_star - Omega_spin[0];

  while (1) {
    dift = Omega_star - Omega_spin[j];
    //printf("j = %d, dif = %.5e, dift = %.5e \n", j, dif, dift);
    if (dift*dif <0) {
      ns = j;
      //printf("ns = %d \n", ns);
      break;
    }
    if (Omega_spin[j] == 0) {
      printf("index not found\n");
      break;
    }
    dif = dift;
    j += 1;
  }
  ns -=1;
  //printf("ns = %d \n", ns);

  printf("initial point: J = %.5e, e = %f \n", Omega_spin[ns], e_spin[ns]);

  new=IMIN(IMAX(ns - 1,0),count-4);
  e_pt = interpolate(&Omega_spin[new], &e_spin[new], Omega_star);
  r_pt = interpolate(&e_spin[new], &r_spin[new], e_pt);
  //printf("r_pt = %f \n", r_pt);
      
  ratio_r = r_pt;
  star->e_center = e_pt;
  MakeSphere(eos, star, e_pt);
  rns(r_pt, e_pt, eos, star);
  printf("The interpolated point has \n%.5f %.5f  %.5f  %.5f %.5f %.5f %.5f %.5f %.5f %.5e \n",
	   star->e_center, star->Mass/MSUN, star->Mass_0/MSUN, Mstat, 
	   star->R_e*1e-5, ratio_r, Rstat, star->Omega/(2.0*PI), 
	   star->Omega_K/(2.0*PI), star->ang_mom); 

  /*fprintf(fpointer_plot, "%f, %f \n", 
    star->e_center, star->Mass/MSUN);*/

  printf("percent error in M0 =  %6.5e, J = %6.5e \n", (M_star - star->Mass/MSUN)/M_star*100, (Omega_star - star->Omega/(2.0*PI))/Omega_star*100);

  //printf("e_pt = %f \n", e_pt);
  start[0] = e_pt;
  start[1] = r_pt;
}

double calc_omega(EOS *eos, NeutronStar *star){
  int initial_i;
  initial_i = (SDIV - 1.0)/2 +1;
  int num = SDIV;

  double r0s, r1s, rho0r, rho1r, ome0r, ome1r, gam0r, gam1r;
  double OMEGA_P, omega_P, vel_P, r_is_P, rho_P;
  double OMEGAP[num], omega[num], velP[num], r_is[num], rho[num];
  double one[num], two[num], three[num];
  double first[num], second[num], third[num];
  double s_used[num];
  int i;
  double Omega_star = 700, Omega_final, r_final, dif, dift; 
  dif = -1e4;

  for(i=initial_i; i<(SDIV-1); i++){

    /* Convert the derivatives from ds to dr values */
    r0s    = star->r_e * star->metric.s_gp[i] / (1.0 - star->metric.s_gp[i]);
    r1s    = star->r_e * pow(1.0 - star->metric.s_gp[i], -2.0);
    rho0r  = star->metric.rho[i][1];
    rho1r  = star->metric.rho_s[i][1] / r1s;
    ome0r  = star->metric.omega[i][1];
    ome1r  = star->metric.omega_s[i][1] / r1s;
    gam0r  = star->metric.gama[i][1];
    gam1r  = star->metric.gama_s[i][1] / r1s;
    //printf("derv = %6.5e\n", star->metric.gama_s[i][1]);

    r_is[i]   = r0s;
    rho[i]    = rho0r;
    omega[i]  = ome0r;
    s_used[i] = star->metric.s_gp[i];

    one[i]    = SQ(ome1r*r0s*exp(-rho0r));
    two[i]    = (2.0/r0s)*(rho1r + gam1r);
    three[i]  = -(SQ(rho1r) - SQ(gam1r)); 
    //printf("gam1r = %6.5e\n", gam1r);
    first[i]  = r0s / (2.0 - r0s*(rho1r - gam1r));
    second[i] = ome1r * r0s * exp(-rho0r);
    third[i]  = one[i] + two[i] + three[i];

    velP[i]   = first[i]*(second[i] + pow(third[i], 0.5));
    //printf("one = %6.5e\n", one[i]);
    //printf("two = %6.5e\n", two[i]);
    //printf("three = %6.5e\n", three[i]);
    //printf("third = %6.5e\n", third[i]);
    //printf("velP = %6.5e\n", velP[i]);

    OMEGAP[i] = omega[i] + (velP[i] / r_is[i])*exp(rho[i]);
    //printf("\t\t\t\t s=%6.5e \t\tOmega=%6.5e\n", s_used[i],OMEGAP[i]);
    //printf("\t\t\t\t\t\ts_used=%6.5e\n", s_used[i]);
    OMEGAP[i] *= C/(sqrt(KAPPA))/(2.0*PI);
    r_is[i] *= pow(KAPPA, 0.5);
    //printf("r = %6.5e Omega = %6.5e\n", r_is[i], OMEGAP[i]);

    dift = Omega_star - OMEGAP[i];
    //printf("dif = %f, dift = %f \n", dif, dift);
     if (dift*dif <0) {
       if (abs(dif)<abs(dift)) {
	 Omega_final = OMEGAP[i-1];
	 r_final = r_is[i-1];
       }
       else {
	 Omega_final = OMEGAP[i];
	 r_final = r_is[i];
       }
       printf("Alfven radius = %6.5e with frequency = %6.5e\n", r_final, Omega_final);
       break;
    }
    dif = dift;

  }
  return r_final;
}

double find_l ( NeutronStar *star, double r_equil) {
  double r_e = star->r_e, l, r_re, r_isco;
  int index;
  double s_gp_point, v_plus_point, rho_point, gama_point, r;
  /*If r_type = 0, the accretion radius is at the isco. If r_type = 1, 
    the accretion radius is at the equlibrium radius (where 
    spin frequency = orbital frequency).*/
  int r_type = 1;
  
  if (r_type == 0) {
    r_isco = 9 * star->Mass/MSUN; //km
    //printf("r_isco =%6.5e\n", r_isco);
    r_re = r_isco/r_e;
  }
  else {
    r_re = r_equil/r_e;
  }
  //printf("r_re = %6.5e\n", r_re);

  s_gp_point = r_re / (1+ r_re);
  //printf("s_gp_point = %.10f \n", s_gp_point);

  index = s_gp_point *(SDIV - 1) +1;
  if (index > SDIV-3) index = 198;
  //printf("index = %d \n", index);
  
  //printf("s_gp = %.10f\n", star->metric.s_gp[index]);
  v_plus_point = interpolate(&star->metric.s_gp[index], &star->v_plus[index], s_gp_point);
  rho_point = interpolate_matrix(&star->metric.s_gp[index], &star->metric.rho[index], s_gp_point);
  gama_point = interpolate_matrix(&star->metric.s_gp[index], &star->metric.gama[index], s_gp_point);
  
  //printf("v_plus_point = %f \n", v_plus_point);
  //printf("rho_point = %6.5e \n", rho_point); 
  //printf("gama_point = %6.5e \n", gama_point);

  r = r_e*s_gp_point /(1- s_gp_point);
  //printf("r_e = %f, r = %f \n", r_e, r);
  
  l = v_plus_point *r *exp((gama_point - rho_point)/2) /pow((1 - v_plus_point*v_plus_point), 0.5); 
  l *= sqrt(KAPPA)*C; // l is in cgs units 
  //printf("l = %6.5e \n", l);
  return l;
}

void loop_M0(EOS *eos, NeutronStar *star, double ratio_r, double e_center, int numstar, double Mstat, double Rstat, double *M0_nospin, double *r_nospin, double rstep, double M0_star, int stop, FILE *fpointer_plot, int direction) {
  int k = 0, check_step=0;

printf("e_c \tMass \t Mass_0\t  StatM\t  Radius   R-ratio StatR   Spin\t       K freq      Ang mom\n");
printf("e15 \tMsun \t Msun\t  Msun\t  km\t   --  \t   km \t   Hz \t       Hz \t   g cm^2 /s \n");

  for (k=0; k<=numstar-1; k++) {
    rns(ratio_r, e_center, eos, star);

    printf("%.5f %.5f  %.5f  %.5f %.5f %.5f %.5f %.5e %.5e %.5e\n",
	   star->e_center, star->Mass/MSUN, star->Mass_0/MSUN, Mstat, 
	   star->R_e*1e-5, ratio_r, Rstat, star->Omega/(2.0*PI), 
	   star->Omega_K/(2.0*PI), star->ang_mom);

    /*fprintf(fpointer, "%7g %7g %7g %6g %4g %8g %4g %6g %7g %6g  %g  %g  %g  %g  %g\n", 
	     star->e_center, star->Mass/MSUN, star->Mass_0/MSUN, Mstat, maxmass, 
	      star->R_e*1e-5, ratio_r, Rstat, star->Omega/(2.0*PI), 
	      star->Omega_K/(2.0*PI), star->ang_mom, T, W, maxradius, 
	      maxmass/maxradius);*/
    
    fprintf(fpointer_plot, "%f, %f \n", 
             star->e_center, star->Mass/MSUN); 

    M0_nospin[k] = star->Mass_0/MSUN;
    r_nospin[k] = ratio_r;

    //Option to make step size smaller close to ratio_r = 1
    if (ratio_r >= 0.99) check_step+=1;
    if (check_step==1) rstep*=0.1;
    //printf("check_step = %d\n", check_step);
    //printf("rstep = %f\n", rstep);

    if (direction ==1) ratio_r -= rstep;
    else ratio_r += rstep;

    //printf("k = %d, stop = %d\n", k, stop);
    if(star->Mass_0/MSUN < M0_star){
      stop += 1;
    }
    if(stop>1 && k>2)
      break;
  }
}


double find_M0(EOS *eos, NeutronStar *star, double e_center, double M0_star, double *M0_nospin, double *r_nospin, int numstar, double Mstat, double Rstat, FILE *fpointer_plot) {
  double dif, dift, r_pt;
  int j, ns, new;

  dif = M0_star - M0_nospin[0];

  for (j=0; j<numstar; j++) {
    dift = M0_star - M0_nospin[j];
    //printf("dif = %f, dift = %f \n", dif, dift);
    if (dift*dif <0) {
      ns = j;
      break;
    }
    dif = dift;
  }
  ns -=1;
  //printf("ns = %d \n", ns);

  printf("Finding point immediately before M0 = %f \ninitial point: M0 = %f, r = %f \n", M0_star, M0_nospin[ns], r_nospin[ns]);

  new=IMIN(IMAX(ns - 1,0),numstar-4);
  r_pt = interpolate(&M0_nospin[new], &r_nospin[new], M0_star);
  //printf("r_pt = %f \n", r_pt);
  rns(r_pt, e_center, eos, star);
  printf("The interpolated point has \n%.5f %.5f  %.5f  %.5f %.5f %.5f %.5f %.5e %.5e %.5e\n",
	   star->e_center, star->Mass/MSUN, star->Mass_0/MSUN, Mstat, 
	   star->R_e*1e-5, r_pt, Rstat, star->Omega/(2.0*PI), 
	   star->Omega_K/(2.0*PI), star->ang_mom);

  fprintf(fpointer_plot, "%f, %f \n", 
           star->e_center, star->Mass/MSUN);

  printf("percent error is %6.5e \n", (M0_star - star->Mass_0/MSUN)/M0_star*100);
  return r_pt;  
  }

double loop_J (EOS *eos, NeutronStar *star, double ratio_r, double e_center, double temp_energy, double maxmass, double J_star, double *J_spin, double *r_spin, double *e_spin, double r_pt, double Mstat, double Rstat, double T, double W, int stop, FILE *fpointer_plot) {
printf("e_c \tMass \t Mass_0\t  StatM\t  Radius   R-ratio StatR   Spin\t       K freq      Ang mom\n");
printf("e15 \tMsun \t Msun\t  Msun\t  km\t   --  \t   km \t   Hz \t       Hz \t   g cm^2 /s \n");

  int j, ierr;
  int count, check_step;
  double M0, ej, energy_value, rstep = 5e-3;
  float e_c[4], M_0[4];
  count = 0;
  
  ratio_r = r_pt;
  temp_energy = e_center;
  //printf("Energy center = %g \n",e_center);
  
  ierr = MakeSphere(eos, star, e_center);
  rns(ratio_r, e_center, eos, star); 
    
  Mstat = star->Mass/MSUN;
  Rstat = star->R_e*1e-5;

  T = (0.5*star->ang_mom*star->Omega)/(C*C);
  W = star->Mp + T - star->Mass;

  /*fprintf(fpointer, "%7g %7g %7g %6g %4g %8g %4g %6g %7g %6g  %g  %g %g  %g  %g\n", star.e_center, star.Mass/MSUN, star.Mass_0/MSUN, Mstat, maxmass, star.R_e*1e-5, ratio_r, Rstat, star.Omega/(2.0*PI), star.Omega_K/(2.0*PI), star.ang_mom, T, W, maxradius, maxmass/maxradius);*/

  printf("%.5f %.5f  %.5f  %.5f %.5f %.5f %.5f %.5e %.5e %.5e\n",
	   star->e_center, star->Mass/MSUN, star->Mass_0/MSUN, Mstat, 
	   star->R_e*1e-5, ratio_r, Rstat, star->Omega/(2.0*PI), 
	 star->Omega_K/(2.0*PI), star->ang_mom);

  fprintf(fpointer_plot, "%f, %f \n", 
	    star->e_center, star->Mass/MSUN);

  J_spin[count] = star->ang_mom;
  e_spin[count] = star->e_center;
  r_spin[count] = ratio_r;
   
  count += 1;

  M0 = star->Mass_0/MSUN;
     
  if (J_star - star->ang_mom > 0) {
    while(1){   

      //Option to make step size smaller close to ratio_r = 1
      if (ratio_r >= 0.99) check_step+=1;
      if (check_step==1) rstep *= 0.1;
      //printf("check_step = %d\n", check_step);
      //printf("rstep = %f\n", rstep);

      ratio_r -= rstep;

      //printf("ratio_r = %f\n", ratio_r);

      for(j=0;j<3;j++){
	ej = temp_energy - 1e-2*j;
       
	//printf("ej = %f\n", ej);

	ierr = MakeSphere(eos, star, ej);
	  
	if(ratio_r < 0.7)
	  rns(0.7, ej, eos, star);
  
	  rns(ratio_r, ej, eos, star); 

	  /*printf("%g %.5f  %.5f  %.5f %.5f %.3f %.4f %.4f \n",
	    ej, star.Mass/MSUN, star.Mass_0/MSUN, Mstat, star.R_e*1e-5, 
	    ratio_r, Rstat, star.Omega/(2.0*PI));*/

	  e_c[j] = ej;
	  M_0[j] = star->Mass_0/MSUN;
	}
	//printf("-----------------------------------------------\n");
	//printf("M_0 = %g\n",M0);
	//printf("energies:  %g, %g, %g\n",e_c[0], e_c[1], e_c[2]);
	//printf("M_0's: %g, %g, %g\n",M_0[0], M_0[1], M_0[2]);

	energy_value = polyinter(M0, e_c, M_0);

	//printf("---------------------------------------------------------------------------\n");
	temp_energy = energy_value;
	ierr = MakeSphere(eos, star, energy_value);
	//rns(ratio_r, energy_value, &eos, &star); 

	if(ratio_r < 0.7)
	  rns(0.7, ej, eos, star);
    
	rns(ratio_r, ej, eos, star); 

        printf("%.5f %.5f  %.5f  %.5f %.5f %.5f %.5f %.5e %.5e %.5e\n",
	   energy_value, star->Mass/MSUN, star->Mass_0/MSUN, Mstat, 
	   star->R_e*1e-5, ratio_r, Rstat, star->Omega/(2.0*PI), 
	   star->Omega_K/(2.0*PI), star->ang_mom);

	fprintf(fpointer_plot, "%f, %f \n", 
	       energy_value, star->Mass/MSUN);

	J_spin[count] = star->ang_mom;
	e_spin[count] = energy_value;
	r_spin[count] = ratio_r;
	//printf("J = %.5e, e = %f, count = %d \n", J_spin[count], e_spin[count], count);
	if(star->ang_mom > J_star)
	  stop += 1;
	if(stop>0 && count>2)
	  break;
	count += 1;

	if( (star->Omega/(2.0*PI)) > (star->Omega_K/(2.0*PI))) {
	  printf("Max reached\n");
	  break;
	}
	if(isnan(star->Mass/MSUN)){
	  printf("Mass is NAN\n");
	  break;
	}
      }
      }

    else {
	while(1){   
	
	  //Option to make step size smaller close to ratio_r = 1
	  if (ratio_r >= 0.99) check_step+=1;
	  if (check_step==1) rstep *= 0.1;
	  //printf("check_step = %d\n", check_step);
	  //printf("rstep = %f\n", rstep);

	  ratio_r += rstep;

	  //printf("ratio_r = %f\n", ratio_r);

	  for(j=0;j<3;j++){
	    ej = temp_energy + 1e-2*j;
	    //printf("ej = %f\n", ej);

	    ierr = MakeSphere(eos, star, ej);
	   
	    if(ratio_r < 0.7)
	      rns(0.7, ej, eos, star);
  
	    rns(ratio_r, ej, eos, star); 

	    //printf("%g %.5f  %.5f  %.5f %.5f %.3f %.4f %.4f \n",
	    //    ej, star.Mass/MSUN, star.Mass_0/MSUN, Mstat, star.R_e*1e-5, ratio_r, Rstat, star.Omega/(2.0*PI));

	    e_c[j] = ej;
	    M_0[j] = star->Mass_0/MSUN;
	  }
	  //printf("-----------------------------------------------\n");
	  //printf("M_0 = %g\n",M0);
	  //printf("energies:  %g, %g, %g\n",e_c[0], e_c[1], e_c[2]);
	  //printf("M_0's: %g, %g, %g\n",M_0[0], M_0[1], M_0[2]);

	  energy_value = polyinter(M0, e_c, M_0);

	  //printf("---------------------------------------------------------------------------\n");
	  temp_energy = energy_value;
	  ierr = MakeSphere(eos, star, energy_value);
	  //rns(ratio_r, energy_value, &eos, &star); 


	  if(ratio_r < 0.7)
	    rns(0.7, ej, eos, star);
    
	  rns(ratio_r, ej, eos, star); 

	  printf("%.5f %.5f  %.5f  %.5f %.5f %.5f %.5f %.5e %.5e %.5e\n",
		 energy_value, star->Mass/MSUN, star->Mass_0/MSUN, Mstat, 
		 star->R_e*1e-5, ratio_r, Rstat, star->Omega/(2.0*PI), 
		 star->Omega_K/(2.0*PI), star->ang_mom);

	  fprintf(fpointer_plot, "%f, %f \n", 
		  energy_value, star->Mass/MSUN);

	  J_spin[count] = star->ang_mom;
	  e_spin[count] = energy_value;
	  r_spin[count] = ratio_r;
	  //printf("J = %.5e, e = %f, count = %d \n", J_spin[count], e_spin[count], count);

	  if(star->ang_mom < J_star)
	    stop += 1;
	  if(stop>0 && count>2)
	    break;
	  count += 1;

	  if( (star->Omega/(2.0*PI)) > (star->Omega_K/(2.0*PI))) {
	    printf("Max reached\n");
	    break;
	  }
	  if(isnan(star->Mass/MSUN)){
	    printf("Mass is NAN\n");
	    break;
	  }
	}

	T = (0.5*star->ang_mom*star->Omega)/(C*C);
	W = star->Mp + T - star->Mass;
	//printf("Mp = %g T = %g Mass = %g  J = %g W = %g T/W = %g \n", star.Mp, T, star.Mass, star.ang_mom, W, T/W);
    //printf("M0 = %g \t Mass_0 = %g\n", M0, star.Mass_0/MSUN);
	/*if((round(M0*100.0)/100.0) == (round(star->Mass_0/MSUN * 100.0)/100.0))
	  fprintf(fpointer, "%7g %7g %7g %6g %4g %8g %4g %6g %7g %6g  %g  %g  %g  %g  %g\n", energy_value, star->Mass/MSUN, star->Mass_0/MSUN, Mstat, maxmass, star-.R_e*1e-5, ratio_r, Rstat, star->Omega/(2.0*PI), star->Omega_K/(2.0*PI), star->ang_mom, T, W, maxradius, maxmass/maxradius);*/

      }
  return count;
}

double find_J (EOS *eos, NeutronStar *star, double ratio_r, int count, int numstar, double M0_star, double J_star, double *J_spin, double *e_spin, double *r_spin, double Mstat, double Rstat, FILE *fpointer_plot) {
  count += 1;
  int ns = 0;

  printf("Finding point immediately before J = %.5e \n", J_star);

  double dif, dift, r_pt, e_pt;
  int j=0, new;
  dif = J_star - J_spin[0];

  while (1) {
    dift = J_star - J_spin[j];
    //printf("j = %d, dif = %.5e, dift = %.5e \n", j, dif, dift);
    if (dift*dif <0) {
      ns = j;
      //printf("ns = %d \n", ns);
      break;
    }
    if (J_spin[j] == 0) {
      printf("index not found\n");
      break;
    }
    dif = dift;
    j += 1;
  }
  ns -=1;
  //printf("ns = %d \n", ns);

  printf("initial point: J = %.5e, e = %f \n", J_spin[ns], e_spin[ns]);

  new=IMIN(IMAX(ns - 1,0),count-4);
  e_pt = interpolate(&J_spin[new], &e_spin[new], J_star);
  r_pt = interpolate(&e_spin[new], &r_spin[new], e_pt);
  //printf("r_pt = %f \n", r_pt);
      
  ratio_r = r_pt;
  star->e_center = e_pt;
  MakeSphere(eos, star, e_pt);
  rns(r_pt, e_pt, eos, star);
  printf("The interpolated point has \n%.5f %.5f  %.5f  %.5f %.5f %.5f %.5f %.5e %.5e %.5e\n",
	   star->e_center, star->Mass/MSUN, star->Mass_0/MSUN, Mstat, 
	   star->R_e*1e-5, ratio_r, Rstat, star->Omega/(2.0*PI), 
	   star->Omega_K/(2.0*PI), star->ang_mom); 

  fprintf(fpointer_plot, "%f, %f \n", 
    star->e_center, star->Mass/MSUN);

  printf("percent error in M0 =  %6.5e, J = %6.5e \n", fabs((M0_star - star->Mass_0/MSUN)/M0_star*100), fabs((J_star - star->ang_mom)/J_star*100));

  return r_pt;
}


