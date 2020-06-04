#if HAVE_CONFIG_H   /* config.h should come before any other includes */
#  include "config.h"
#endif


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <float.h>

/*
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
*/
#include <gsl/gsl_cdf.h>
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
    glbInitExperiment("exp/DUNE_GLoBES.glb",&glb_experiment_list[0],&glb_num_of_exps);
//    glbInitExperiment("exp/MOMENT_FIX_FLUX_150KM_addATM_NC.glb",&glb_experiment_list[0],&glb_num_of_exps);
    /*Initialize FASE*/
    MODEL_init(5);
    
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
    double th12_true=33.82*degree;
    double th13_true=8.61*degree;
    double th23_true=49.6*degree;
    double dCP_true=215*degree;
    double DM21_true=7.39e-5;
    double DM31_true=2.525e-3;
    /* Initialize parameter and projection vector(s) */
    glb_params true_values = glbAllocParams();
    glb_params test_values = glbAllocParams();
    glb_params input_errors = glbAllocParams();
    glb_params centers = glbAllocParams();
    double osc_para[6];
    double M_para[5];
    M_para[0]=0; M_para[1]=0; M_para[2]=dCP_true; M_para[3]=DM21_true; M_para[4]=DM31_true;
    MtoS(osc_para, M_para);
    /*define model parameter: true_values( x, eta, r, ma)*/
    glbDefineParams(test_values,th13_true,th23_true,DM21_true,DM31_true,-1,0);
    glbSetDensityParams(true_values,1.0,GLB_ALL);
    int i;
    
    glbDefineParams(true_values,th12_true,th13_true,th23_true,dCP_true,DM21_true,DM31_true);
    glbSetDensityParams(true_values,1.0,GLB_ALL);

    
    /*assign the parameter set*/
    PARA=MODEL; /* if the user is interested in the sensitivity on model parameters*/
    PARA=STAN; /* if the user is interested in the sensitivity on Standard Oscillation parameters*/
    /*true value setup for GLoBES*/
    glbSetOscillationParameters(true_values);
    glbSetRates();
    
    /*projection setup*/
    glb_projection free = glbAllocProjection();
    glbDefineProjection(free,GLB_FREE,GLB_FREE,GLB_FREE,GLB_FREE,GLB_FIXED,GLB_FIXED);
    glbSetDensityProjectionFlag(free, GLB_FREE, GLB_ALL);

    /*prior setup*/
    glbRegisterPriorFunction(FASE_prior_OSC,NULL,NULL,NULL);
    UPPER_prior[0]=36.27*degree;  LOWER_prior[0]=31.61*degree;   Central_prior[0]=33.82*degree;
    UPPER_prior[1]=8.99*degree;   LOWER_prior[1]=8.22*degree;    Central_prior[1]=8.61*degree;
    UPPER_prior[2]=52.4*degree;   LOWER_prior[2]=40.3*degree;    Central_prior[2]=49.6*degree;
    UPPER_prior[3]=392*degree;    LOWER_prior[3]=125*degree;     Central_prior[3]=215*degree;
    UPPER_prior[4]=8.01e-5;       LOWER_prior[4]=6.79e-5;        Central_prior[4]=7.39e-5;
    UPPER_prior[5]=2.625e-3;      LOWER_prior[5]=2.427e-3;       Central_prior[5]=2.525e-3;
    
    
    for (i=0;i<6;i++) {UPPER_prior[i]=fabs(UPPER_prior[i]-Central_prior[i])/3;
        LOWER_prior[i]=fabs(LOWER_prior[i]-Central_prior[i])/3;}

//        for (i=0;i<6;i++) {UPPER_prior[i]=0; LOWER_prior[i]=0;}
    /* centers and input_errors will not be used, but need to be given in for central values
    and prior width. Otherwise, GLoBES will complain. */
    for (i=0; i<6; i++) glbSetOscParams(centers,0.1,i); glbSetDensityParams(centers,1.0,GLB_ALL); glbCopyParams(centers,input_errors);
    glbSetCentralValues(centers); glbSetInputErrors(input_errors);
    
    
    /*two loops for chi^2(x,eta)*/
    float th23,dCP,dth23,ddCP,lower_th23,upper_th23,lower_dCP,upper_dCP;
    FILE* File=fopen("data/TM1_th23_dCP(DUNE).dat", "w");
    lower_th23=40; upper_th23=52.5; lower_dCP=125; upper_dCP=395;
    dth23=(upper_th23-lower_th23)/100; ddCP=(upper_dCP-lower_dCP)/100;
    int dof=2;
    glbSetProjection(free);
    for (th23=lower_th23;th23<=upper_th23;th23=th23+dth23){
    for (dCP=lower_dCP;dCP<=upper_dCP;dCP=dCP+ddCP){
    glbSetOscParams(true_values,th23*degree,GLB_THETA_23);
    glbSetOscParams(true_values,dCP*degree,GLB_DELTA_CP);
    PARA=STAN; /* if the user is interested in the sensitivity on Standard Oscillation parameters*/
    glbSetOscillationParameters(true_values);
    glbSetRates();
    PARA=MODEL; /* if the user is interested in the sensitivity on Standard Oscillation parameters*/
        
    glbSetOscParams(test_values,-1,4);
        float res=glbChiNP(test_values,NULL,GLB_ALL);
        glbSetOscParams(test_values,1,4);
        float res2=glbChiNP(test_values,NULL,GLB_ALL);
        if(res2<res) res=res2;
        
        fprintf(File,"%f %f %f\n",th23,dCP,gsl_cdf_chisq_Qinv(gsl_cdf_chisq_Q(fabs(res),dof),1));
         } fprintf(File,"\n");
    }
 

    /* Destroy parameter and projection vector(s) */
    glbFreeParams(true_values);
    glbFreeParams(test_values); 
    glbFreeParams(input_errors); 
    glbFreeProjection(free);
    exit(0);
}

