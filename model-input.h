#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <float.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>


int N_M;

double complex TDModelY(double x, double eta, double r, double Ma);
double complex TDModelZ(double x, double eta, double r, double Ma);
double complex TDModelW(double x, double eta, double r, double Ma);
double TDpsi(double x, double eta, double r, double Ma);
double TDsinpsi(double x, double eta, double r, double Ma);
double TDcospsi(double x, double eta, double r, double Ma);
double TDtheta(double x, double eta, double r, double Ma);
double TDth12(double x, double eta, double r, double Ma);
double TDth13(double x, double eta, double r, double Ma);
double TDth23(double x, double eta, double r, double Ma);
double TDdCP(double x, double eta, double r, double Ma);
double TDdm21(double x, double eta, double r, double Ma);
double TDdm31(double x, double eta, double r, double Ma);
double MODEL_init();
double MtoS(double osc_para[6], double M_para[]);
double model_restriction(double model []);
