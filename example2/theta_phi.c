#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <float.h>
#include <globes/globes.h>   /* GLoBES library */
#include "FaSE_GLoBES.h"
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
    double th12_true=33.82*degree;
    double th13_true=8.61*degree;
    double th23_true=49.6*degree;
    double dCP_true=215*degree;
    double DM21_true=7.39e-5;
    double DM31_true=2.525e-3;
    double phi_test=1*M_PI;
    double theta_test=0.19*M_PI;
    /* Initialize parameter and projection vector(s) */
    glb_params true_values = glbAllocParams();
    glb_params test_values = glbAllocParams();
    glb_params input_errors = glbAllocParams();
    glb_params centers = glbAllocParams();
    
    /*define model parameter: true_values( x, eta, r, ma)*/
    glbDefineParams(true_values,th12_true,th13_true,th23_true,dCP_true,DM21_true,DM31_true);
    glbSetDensityParams(true_values,1.0,GLB_ALL);

    glbDefineParams(test_values,theta_test,phi_test,DM21_true,DM31_true,0,0);
    glbSetDensityParams(true_values,1.0,GLB_ALL);
    
    /*assign the parameter set*/
    PARA=STAN; /* if the user is interested in the sensitivity on Standard Oscillation parameters*/
    PARA=MODEL; /* if the user is interested in the sensitivity on model parameters*/
    
    /*true value setup for GLoBES*/
    PARA=STAN;
    glbSetOscillationParameters(true_values);
    glbSetRates();
    PARA=MODEL;
    /*projection setup*/
    glb_projection projection = glbAllocProjection();
    glb_projection free = glbAllocProjection();
    glbDefineProjection(free,GLB_FREE,GLB_FREE,GLB_FREE,GLB_FREE,GLB_FIXED,GLB_FIXED);
    glbSetDensityProjectionFlag(free, GLB_FREE, GLB_ALL);
    glbDefineProjection(projection,GLB_FIXED,GLB_FIXED,GLB_FREE,GLB_FREE,GLB_FIXED,GLB_FIXED);
    glbSetDensityProjectionFlag(projection, GLB_FREE, GLB_ALL);

    /*prior setup*/
    glbRegisterPriorFunction(FASE_prior_OSC,NULL,NULL,NULL);
    
    UPPER_prior[0]=36.27*degree;  LOWER_prior[0]=31.61*degree;   Central_prior[0]=33.82*degree;
    UPPER_prior[1]=8.99*degree; LOWER_prior[1]=8.22*degree;  Central_prior[1]=8.61*degree;
    UPPER_prior[2]=52.4*degree;  LOWER_prior[2]=40.3*degree;   Central_prior[2]=49.6*degree;
    UPPER_prior[3]=392*degree;  LOWER_prior[3]=125*degree;   Central_prior[3]=215*degree;
    UPPER_prior[4]=8.01e-5;     LOWER_prior[4]=6.79e-5;      Central_prior[4]=7.39e-5;
    UPPER_prior[5]=2.625e-3;     LOWER_prior[5]=2.427e-3;      Central_prior[5]=2.525e-3;
    
    int i;
    for (i=0;i<6;i++) {UPPER_prior[i]=0; LOWER_prior[i]=0;}

    /* centers and input_errors will not be used, but need to be given in for central values
    and prior width. Otherwise, GLoBES will complain. */
    for (i=0; i<6; i++) glbSetOscParams(centers,0,i); glbSetDensityParams(centers,1.0,GLB_ALL); glbCopyParams(centers,input_errors);
    glbSetCentralValues(centers); glbSetInputErrors(input_errors);
    
    glbSetProjection(free);
    float res0=glbChiNP(test_values,NULL,GLB_ALL); /*chi-square minimum*/

    
    float x,y,r,ma,dx,dy,lower_x,upper_x,lower_y,upper_y;
    FILE* File=fopen("data/theta_phi(MOMENT).dat", "w");
    lower_x=0.16*M_PI; upper_x=0.25*M_PI; lower_y=0.85*M_PI; upper_y=1.1*M_PI;
    dx=(upper_x-lower_x)/100; dy=(upper_y-lower_y)/100;
    glbSetProjection(projection);
    
    for (x=lower_x;x<=upper_x;x=x+dx){
    for (y=lower_y;y<=upper_y;y=y+dy){
     glbSetOscParams(test_values,x,0); glbSetOscParams(test_values,y,1);
    float res=glbChiNP(test_values,NULL,GLB_ALL);
            fprintf(File,"%f %f %f\n",x/M_PI,y/M_PI,res-res0);
    } fprintf(File,"\n");
    }

    /* Destroy parameter and projection vector(s) */
    glbFreeParams(true_values);
    glbFreeParams(test_values); 
    glbFreeParams(input_errors); 
    glbFreeParams(centers);
    glbFreeProjection(free);
    glbFreeProjection(projection);
    exit(0);
}

