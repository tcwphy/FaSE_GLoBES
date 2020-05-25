#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_math.h>
#include <globes/globes.h>   /* GLoBES library */
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include "FASE_GLoBES.h"
/* Macros */
int N_M;


double MODEL_init(int N) /*initialise the model-input*/
{
    if (N>6){printf("the number of parameters cannot be larger than 6.\n");}
    N_M=N;
    return 0;
}

double MtoS(double osc_para[6], double M_para[])
/*revert model parameters to
 the standard oscillation parameters */
{
    
    /*theta12*/        osc_para[0]=;
    /*theta13*/        osc_para[1]=;
    /*theta23*/        osc_para[2]=;
    /*deltaCP*/        osc_para[3]=;
    /*Delta_m^2_{21}*/ osc_para[4]=;
    /*Delta_m^2_{31}*/ osc_para[5]=;


    return 0;
}

double model_restriction(double model [])
{
    
    return 0;
}

