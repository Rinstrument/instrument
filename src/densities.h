#ifndef DENSITIES
#define DENSITIES

double lognorm_dens_dx(double x, double mean, double sd);
double logtruncnorm_dens_dx(double x, double mean, double sd, double a, double b);
double log_normd_dx(double x, double mean, double sd);
double log_unifd_dx(double x, double l, double u);

#endif