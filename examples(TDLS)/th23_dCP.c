#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <float.h>
#include <gsl/gsl_cdf.h>
#include <globes/globes.h>   /* GLoBES library */
#include "FASE_GLoBES.h"
#include "model-input.h"

/***************************************************************************
 *                            M A I N   P R O G R A M                      *
 ***************************************************************************/

int main(int argc, char *argv[])
{
    /* Initialize libglobes */
    glbInit(argv[0]);
    glbInitExperiment("exp/MOMENT_150KM_addATM_NC.glb",&glb_experiment_list[0],&glb_num_of_exps);

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
    double M_para[4],OSC_para[6];
    M_para[0]=x_true; M_para[1]=eta_true; M_para[2]=r_true; M_para[3]=ma_true;
    
    MtoS(OSC_para, M_para);
    double th12_true=OSC_para[0];
    double th13_true=OSC_para[1];
    double th23_true=OSC_para[2];
    double dCP_true=OSC_para[3];
    double DM21_true=OSC_para[4];
    double DM31_true=OSC_para[5];

    /* Initialize parameter and projection vector(s) */
    glb_params true_values = glbAllocParams();
    glb_params test_values = glbAllocParams();
    glb_params input_errors = glbAllocParams();
    glb_params centers = glbAllocParams();
    glb_params out = glbAllocParams();
    glb_params out2 = glbAllocParams();
    /*define model parameter: true_values( x, eta, r, ma)*/
    glbDefineParams(test_values,x_true,eta_true,r_true,ma_true,0,0);
    glbSetDensityParams(test_values,1.0,GLB_ALL);
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
    glbDefineProjection(free,GLB_FIXED,GLB_FIXED,GLB_FREE,GLB_FREE,GLB_FIXED,GLB_FIXED);
    glbSetDensityProjectionFlag(free, GLB_FREE, GLB_ALL);
    
    /*prior setup*/
    glbRegisterPriorFunction(FASE_prior_OSC,NULL,NULL,NULL);
    UPPER_prior[0]=36.27*degree;  LOWER_prior[0]=31.61*degree;   Central_prior[0]=33.82*degree;
    UPPER_prior[1]=8.99*degree; LOWER_prior[1]=8.22*degree;  Central_prior[1]=8.61*degree;
    UPPER_prior[2]=52.4*degree;  LOWER_prior[2]=40.3*degree;   Central_prior[2]=49.7*degree;
    UPPER_prior[3]=392*degree;  LOWER_prior[3]=125*degree;   Central_prior[3]=216*degree;
    UPPER_prior[4]=8.01e-5;     LOWER_prior[4]=6.79e-5;      Central_prior[4]=7.39e-5;
    UPPER_prior[5]=2.625e-3;     LOWER_prior[5]=2.427e-3;      Central_prior[5]=2.525e-3;
    

    for (i=0;i<6;i++) {UPPER_prior[i]=fabs(UPPER_prior[i]-Central_prior[i])/3;
        LOWER_prior[i]=fabs(LOWER_prior[i]-Central_prior[i])/3;}
    /* centers and input_errors will not be used, but need to be given in for central values
    and prior width. Otherwise, GLoBES will complain. */
    for (i=0; i<6; i++) glbSetOscParams(centers,0,i); glbSetDensityParams(centers,1.0,GLB_ALL);
    glbSetDensityParams(centers,1.0,GLB_ALL);
    glbSetDensityParams(input_errors,0.05,GLB_ALL);
    glbSetCentralValues(centers); glbSetInputErrors(input_errors);
    
    /*two loops for chi^2(x,eta)*/
    float th23,dCP,dth23,ddCP,lower_th23,upper_th23,lower_dCP,upper_dCP;
    FILE* File=fopen("data/th23_dCP(MOMENT).dat", "w");
    lower_th23=40; upper_th23=52.5; lower_dCP=125; upper_dCP=395;
    dth23=(upper_th23-lower_th23)/50; ddCP=(upper_dCP-lower_dCP)/50;
    int dof=4;
    glbSetProjection(free);
    
    PARA=MODEL;
    float res0=glbChiNP(test_values,out,GLB_ALL);
    glbCopyParams(out,out2);
    for (th23=lower_th23;th23<=upper_th23;th23=th23+dth23){
    for (dCP=lower_dCP;dCP<=upper_dCP;dCP=dCP+ddCP){
    glbSetOscParams(true_values,th23*degree,GLB_THETA_23);
    glbSetOscParams(true_values,dCP*degree,GLB_DELTA_CP);
    PARA=STAN; /* if the user is interested in the sensitivity on Standard Oscillation parameters*/
    glbSetOscillationParameters(true_values);
    glbSetRates();
    PARA=MODEL; /* if the user is interested in the sensitivity on Standard Oscillation parameters*/
        float res=glbChiNP(test_values,out,GLB_ALL);
        float x_out   = glbGetOscParams(out,0);
        float eta_out = glbGetOscParams(out,1);
        float r_out   = glbGetOscParams(out,2);
        float ma_out  = glbGetOscParams(out,3);
            fprintf(File,"%f %f %f %f %g %g %g %g\n",th23,dCP,gsl_cdf_chisq_Qinv(gsl_cdf_chisq_Q(fabs(res),dof),1),gsl_cdf_chisq_Qinv(gsl_cdf_chisq_Q(fabs(res-res0),dof),1),x_out,eta_out/M_PI,r_out,ma_out);
         } fprintf(File,"\n");
    }
 
    /* Destroy parameter and projection vector(s) */
    glbFreeParams(true_values);
    glbFreeParams(test_values); 
    glbFreeParams(input_errors);
    glbFreeParams(centers);
    glbFreeParams(out);
    glbFreeParams(out2);
    glbFreeProjection(free);
    exit(0);
}

