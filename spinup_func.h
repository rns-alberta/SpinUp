#include "struct.h"

double find_l (NeutronStar *star, double r_equil);

void loop_M0 (EOS *eos, NeutronStar *star, double ratio_r, double e_center, int numstar, double Mstat, double Rstat, double *M0_nospin, double *r_nospin, double rstep, double M0_star, int stop, FILE *fpointer_plot, int direction);

double find_M0(EOS *eos, NeutronStar *star, double e_center, double M0_star, double *M0_nospin, double *r_nospin, int numstar, double Mstat, double Rstat, FILE *fpointer_plot);

double loop_J(EOS *eos, NeutronStar *star, double ratio_r, double e_center, double temp_energy, double maxmass, double J_star, double *J_spin, double *r_spin, double *e_spin, double r_pt, double Mstat, double Rstat, double T, double W, int stop, FILE *fpointer_plot);

double find_J (EOS *eos, NeutronStar *star, double ratio_r, int count, int numstar, double M0_star, double J_star, double *J_spin, double *e_spin, double *r_spin, double Mstat, double Rstat, FILE *fpointer_plot);

void loop_M (EOS *eos, NeutronStar *star, double ratio_r, double e_center, int numstar, double Mstat, double Rstat, double *M_nospin, double *r_nospin, double rstep, double M_star, int stop, FILE *fpointer_plot);

double find_M(EOS *eos, NeutronStar *star, double e_center, double M_star, double *M_nospin, double *r_nospin, int numstar, double Mstat, double Rstat, FILE *fpointer_plot);

double loop_spin(EOS *eos, NeutronStar *star, double ratio_r, double e_center, double temp_energy, double maxmass, double Omega_star, double *Omega_spin, double *r_spin, double *e_spin, double r_pt, double Mstat, double Rstat, double T, double W, int stop, FILE *fpointer_plot);

void find_spin (EOS *eos, NeutronStar *star, double ratio_r, int count, int numstar, double M0_star, double Omega_star, double *Omega_spin, double *e_spin, double *r_spin, double Mstat, double Rstat, FILE *fpointer_plot, double start[]);

double calc_omega(EOS *eos, NeutronStar *star);










