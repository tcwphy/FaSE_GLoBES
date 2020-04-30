#if HAVE_CONFIG_H   /* config.h should come before any other includes */
#  include "config.h"
#endif


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


#include <globes/globes.h>   /* GLoBES library */
#include "FASE_GLoBES.h"
//#include "model-input_diag.h"
#include "model-input.h"

/***************************************************************************
 *                            M A I N   P R O G R A M                      *
 ***************************************************************************/

int main(int argc, char *argv[])
{
    /* Initialize libglobes */
    glbInit(argv[0]);
    glbInitExperiment("exp/MOMENT_FIX_FLUX_150KM_addATM_NC.glb",&glb_experiment_list[0],&glb_num_of_exps);
    /*Initialize FASE*/
    MODEL_init(4);
    
    /* Register non-standard probability engine. This has to be done
     * before any calls to glbAllocParams() or glbAllocProjections() */
    glb_init_probability_engine();
    glbRegisterProbabilityEngine(6,
                                 &FASE_glb_probability_matrix,
                                 &FASE_glb_set_oscillation_parameters,
                                 &FASE_glb_get_oscillation_parameters,
                                 NULL);
    /*define the true value*/
    float degree   = M_PI/180;
    float x_true,eta_true,r_true,ma_true;
    x_true=-3.65029; eta_true=1.13067*M_PI; r_true=0.511325; ma_true=3.71199e-3;
    
    /* Initialize parameter and projection vector(s) */
    glb_params true_values = glbAllocParams();
    glb_params test_values = glbAllocParams();
    glb_params input_errors = glbAllocParams();
    glb_params centers = glbAllocParams();
    
    /*define model parameter: true_values( x, eta, r, ma)*/
    glbDefineParams(true_values,x_true,eta_true,r_true,ma_true,0,0);
    glbSetDensityParams(true_values,1.0,GLB_ALL);
    glbCopyParams(true_values,test_values);
    
    /*assign the parameter set*/
    PARA=STAN; /* if the user is interested in the sensitivity on Standard Oscillation parameters*/ /*Mark4*/
    PARA=MODEL; /* if the user is interested in the sensitivity on model parameters*/
    /*true value setup for GLoBES*/
    glbSetOscillationParameters(true_values);
    glbSetRates();
    
    /*projection setup*/
    glb_projection projection = glbAllocProjection();
    glb_projection free = glbAllocProjection();
    glbDefineProjection(free,GLB_FREE,GLB_FREE,GLB_FREE,GLB_FREE,GLB_FIXED,GLB_FIXED);
    glbSetDensityProjectionFlag(free, GLB_FIXED, GLB_ALL);
    glbDefineProjection(projection,GLB_FIXED,GLB_FIXED,GLB_FREE,GLB_FREE,GLB_FIXED,GLB_FIXED);
    glbSetDensityProjectionFlag(projection, GLB_FIXED, GLB_ALL);

    /*prior setup*/
    glbRegisterPriorFunction(FASE_prior_OSC,NULL,NULL,NULL);
    
    UPPER_prior[0]=36.27*degree;  LOWER_prior[0]=31.61*degree;   Central_prior[0]=33.82*degree;
    UPPER_prior[1]=8.99*degree; LOWER_prior[1]=8.22*degree;  Central_prior[1]=8.61*degree;
    UPPER_prior[2]=52.4*degree;  LOWER_prior[2]=40.3*degree;   Central_prior[2]=49.6*degree;
    UPPER_prior[3]=392*degree;  LOWER_prior[3]=125*degree;   Central_prior[3]=215*degree;
    UPPER_prior[4]=8.01e-5;     LOWER_prior[4]=6.79e-5;      Central_prior[4]=7.39e-5;
    UPPER_prior[5]=2.625e-3;     LOWER_prior[5]=2.427e-3;      Central_prior[5]=2.525e-3;
    
    int i;
    for (i=0;i<6;i++) {UPPER_prior[i]=fabs(UPPER_prior[i]-Central_prior[i])/3;
        LOWER_prior[i]=fabs(LOWER_prior[i]-Central_prior[i])/3;}
    /* centers and input_errors will not be used, but need to be given in for central values
    and prior width. Otherwise, GLoBES will complain. */
    for (i=0; i<6; i++) glbSetOscParams(centers,0,i); glbSetDensityParams(centers,1.0,GLB_ALL); glbCopyParams(centers,input_errors);
    glbSetCentralValues(centers); glbSetInputErrors(input_errors);
    
    /*for chi^2 minimum*/
    glbSetProjection(free);
    float res0=glbChiNP(true_values,NULL,GLB_ALL);

    /*two loops for chi^2(x,eta)*/
    float x,eta,r,ma,dx,deta,lower_x,upper_x,lower_eta,upper_eta;
    FILE* File=fopen("data/constraint_x_eta_test.dat", "w");
    lower_x=-5; upper_x=-3; lower_eta=0.8*M_PI; upper_eta=1.5*M_PI;
    dx=(upper_x-lower_x)/100; deta=(upper_eta-lower_eta)/100;
    glbSetProjection(projection);
    for (x=lower_x;x<=upper_x;x=x+dx){
    for (eta=lower_eta;eta<=upper_eta;eta=eta+deta){
     glbSetOscParams(test_values,x,0); glbSetOscParams(test_values,eta,1);
        float res=glbChiNP(test_values,NULL,GLB_ALL);
            fprintf(File,"%f %f %f\n",x,eta/M_PI,res-res0);
         } fprintf(File,"\n");
    }
 

    /* Destroy parameter and projection vector(s) */
    glbFreeParams(true_values);
    glbFreeParams(test_values); 
    glbFreeParams(input_errors); 
    glbFreeProjection(free);
    glbFreeProjection(projection);
    exit(0);
}

