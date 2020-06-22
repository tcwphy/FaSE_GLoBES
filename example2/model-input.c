/*--------------------------------------------------------------------
 This example is using the warped flavor symmetry model
 as an example. More details can be checked in arXiv:1906.02208[hep-ph]
 The user using this model, should cite the reference:
 
 @article{Chen:2015jta,
 author = "Chen, Peng and Ding, Gui-Jun and Rojas, Alma. D. and Vaquera-Araujo, C.A. and Valle, J.W.F.",
 title = "{Warped flavor symmetry predictions for neutrino physics}",
 eprint = "1509.06683",
 archivePrefix = "arXiv",
 primaryClass = "hep-ph",
 reportNumber = "IFIC-15-XX",
 doi = "10.1007/JHEP01(2016)007",
 journal = "JHEP",
 volume = "01",
 pages = "007",
 year = "2016"
 }
 -----------------------------------------------------------------------*/
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
#include "FaSE_GLoBES.h"
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
{// set up for normal ordering according to arXiv:1509.06683
    /*example: tri-direct*/
    
    double theta_nu   = M_para[0];
    double phi_nu     = M_para[1];
    double sdm        = M_para[2];
    double ldm        = M_para[3];
    
    osc_para[4]=sdm;
    osc_para[5]=ldm;
    
    double sin2the= sin(2*theta_nu);
    double cosphi = cos(phi_nu);
    double s12= sqrt(1/(2-sin2the*cosphi));
    double s13= sqrt((1+sin2the*cosphi)/3);
    double s23= sqrt((1-sin2the*sin(M_PI/6-phi_nu))/(2-sin2the*cosphi));
    double JCP= -1*cos(2*theta_nu)/(6*sqrt(3));
    
    double theta12=asin(s12);
    double theta13=asin(s13);
    double theta23=asin(s23);
    
    double sindCP=JCP/(cos(theta12)*s12*cos(theta23)*s23*cos(theta13)*cos(theta13)*s13);
    double deltaCP=asin(sindCP);
  
    
    osc_para[0]=theta12;
    osc_para[1]=theta13;
    osc_para[2]=theta23;
    osc_para[3]=deltaCP;
    return 0;
}

double model_restriction(double M_para [])
{// set up for normal ordering according to arXiv:1509.06683
    
    double theta_nu   = M_para[0];
    double phi_nu     = M_para[1];
    double sdm        = M_para[2];
    double ldm        = M_para[3];
    if(theta_nu<0) {return 1;} /*for NO*/
    if(theta_nu>.5*M_PI) {return 1;} /*for NO*/
    if(ldm<0) {return 1;} /*for NO*/
    
    return 0;
}

