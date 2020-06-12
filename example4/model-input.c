/*--------------------------------------------------------------------
This example is using a S4 modular flavour symmetry model
as an example. More details can be checked in arXiv:1906.02208[hep-ph]
The user using this model, should cite the reference:

 @article{deMedeirosVarzielas:2019cyj,
 author = "de Medeiros Varzielas, Ivo and King, Stephen F. and Zhou, Ye-Ling",
 title = "{Multiple modular symmetries as the origin of flavor}",
 eprint = "1906.02208",
 archivePrefix = "arXiv",
 primaryClass = "hep-ph",
 doi = "10.1103/PhysRevD.101.055033",
 journal = "Phys. Rev. D",
 volume = "101",
 number = "5",
 pages = "055033",
 year
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
#include "FASE_GLoBES.h"
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
    
    double mu1          = M_para[0]; //
    double mu2          = M_para[1];
    double theta_R      = M_para[2];
    double alpha1       = M_para[3];
    double alpha2       = M_para[4];

    double complex cR=cos(theta_R)*cos(alpha1)+I*cos(theta_R)*sin(alpha1);
    double complex sR=sin(theta_R)*cos(alpha2)+I*sin(theta_R)*sin(alpha2);
    double complex sR_con=sin(theta_R)*cos(alpha2)-I*sin(theta_R)*sin(alpha2);
    double complex cR_con=cos(theta_R)*cos(alpha1)-I*cos(theta_R)*sin(alpha1);
    double complex cRs=cR*cR;
    double complex sRs=sR*sR;
    double complex cRs_con=cR_con*cR_con;
    double complex sRs_con=sR_con*sR_con;
    
    double complex A=mu1*cRs+mu2*sRs_con;
    double complex B=mu1*sRs+mu2*cRs_con;
    double complex C=mu1*cR*sR-mu2*cR_con*sR_con;
    
    double complex omega= -.5+I*.5*sqrt(3);
    double complex omega_s=omega*omega;
    
    double complex M[]={A,A*(-2*omega_s)-C,-2*A*omega+C,-2*omega_s*A-C,4*omega*A+B+4*omega_s*C,\
        4*A-B+I*2*sqrt(3)*C,-2*omega*A+C,4*A-B+2*I*sqrt(3)*C,4*omega_s*A+B-4*omega*C};
    
       STAN_OSC(M,osc_para);
 
    
    
    return 0;
}

double model_restriction(double M_para [])
{
    double mu1          = M_para[0]; //
    double mu2          = M_para[1];
    double theta_R      = M_para[2];
    double alpha1       = M_para[3];
    double alpha2       = M_para[4];
    
    return 0;
}

