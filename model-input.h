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
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
//#include "FASE_GLoBES.h"

int N_M;
/******************************************the user-defined relation*/
/*******************************************************************/

double MODEL_init();
/********************************************diagonaliser ********/
//int STAN_OSC(double complex M[], double output[6]);
//int ModelTO( double OSC_PARAMS[6],double M_para[]);
/*******************************************************************/

double MtoS(double osc_para[6], double M_para[]);
double model_restriction(double model []);
