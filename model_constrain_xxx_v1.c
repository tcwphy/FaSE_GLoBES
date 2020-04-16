/* GLoBES -- General LOng Baseline Experiment Simulator
 * (C) 2002 - 2007,  The GLoBES Team
 *
 * GLoBES is mainly intended for academic purposes. Proper
 * credit must be given if you use GLoBES or parts of it. Please
 * read the section 'Credit' in the README file.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * Example: Non-Standard-Interactions and user-defined priors
 * Compile with ``make example6''
 *
 * This example is similar to Chapter 4 of hep-ph/0502147
 */
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
    glbInitExperiment("exp/MOMENT_150KM_addATM_NC.glb",&glb_experiment_list[0],&glb_num_of_exps);
    
     int i;
//     int M_N=4;
    MODEL_init(4);

     float degree   = M_PI/180;
     float x,eta,r,ma,x_true,eta_true,r_true,ma_true,dx,deta;
     float lower_x,upper_x,lower_eta,upper_eta;
     float theta12,theta23,theta13,deltacp,sdm,ldm,Dm21,Dm31;
     float theta12_MODEL,theta13_MODEL,theta23_MODEL,delta_MODEL,sdm_MODEL,ldm_MODEL;

     FILE* File=fopen("data/constraint_x_eta.dat", "w");

    lower_x=-9; upper_x=-3; lower_eta=0.8*M_PI; upper_eta=2.0*M_PI;
    dx=(upper_x-lower_x)/500; deta=(upper_eta-lower_eta)/500;


    UPPER_prior[0]=36.27*degree;  LOWER_prior[0]=31.61*degree;   Central_prior[0]=33.82*degree;
    UPPER_prior[1]=8.99*degree; LOWER_prior[1]=8.22*degree;  Central_prior[1]=8.61*degree;
    UPPER_prior[2]=52.4*degree;  LOWER_prior[2]=40.3*degree;   Central_prior[2]=49.6*degree;
    UPPER_prior[3]=392*degree;  LOWER_prior[3]=125*degree;   Central_prior[3]=215*degree;
    UPPER_prior[4]=8.01e-5;     LOWER_prior[4]=6.79e-5;      Central_prior[4]=7.39e-5;
    UPPER_prior[5]=2.625e-3;     LOWER_prior[5]=2.47e-3;      Central_prior[5]=2.525e-3;
    
    for (i=0;i<6;i++) {UPPER_prior[i]=fabs(UPPER_prior[i]-Central_prior[i])/3;
                       LOWER_prior[i]=fabs(LOWER_prior[i]-Central_prior[i])/3;}
    
    x_true=-3.65029; eta_true=1.13067*M_PI; r_true=0.511325; ma_true=0.00371199;
    
    PARA=STAN; /* if the user is interested in the sensitivity on Standard Oscillation parameters*/ /*Mark4*/
    PARA=MODEL; /* if the user is interested in the sensitivity on model parameters*/
    
    /* Register non-standard probability engine. This has to be done
     * before any calls to glbAllocParams() or glbAllocProjections() */
    glb_init_probability_engine();
    glbRegisterProbabilityEngine(6,
                                 &FASE_glb_probability_matrix,
                                 &FASE_glb_set_oscillation_parameters,
                                 &FASE_glb_get_oscillation_parameters,
                                 NULL);
    
    /* Initialize parameter and projection vector(s) */
    glb_params true_values = glbAllocParams();
    glb_params test_values = glbAllocParams();
    glb_params input_errors = glbAllocParams();
    glb_params out = glbAllocParams();
    glb_params centers = glbAllocParams();
    glb_projection projection = glbAllocProjection();
    glb_projection free = glbAllocProjection();    
    
    /*define model parameter: true_values( x, eta, r, ma)*/
    glbDefineParams(true_values,x_true,eta_true,r_true,ma_true,0,0);
    glbDefineParams(test_values,x_true,eta_true,r_true,ma_true,0,0);
    glbDefineParams(centers,33.82*degree,8.61*degree,49.7*degree,217*degree,7.39e-5,2.525e-3);
    
    glbSetDensityParams(true_values,1.0,GLB_ALL);
    glbSetDensityParams(test_values,1.0,GLB_ALL);
    glbSetDensityParams(centers,1.0,GLB_ALL);
 
    glbDefineProjection(free,GLB_FREE,GLB_FREE,GLB_FREE,GLB_FREE,GLB_FIXED,GLB_FIXED);
    glbSetDensityProjectionFlag(free, GLB_FREE, GLB_ALL);

    glbDefineProjection(projection,GLB_FIXED,GLB_FIXED,GLB_FREE,GLB_FREE,GLB_FIXED,GLB_FIXED);
    glbSetDensityProjectionFlag(projection, GLB_FREE, GLB_ALL);


    glbDefineParams(input_errors,0.78*degree,0.127*degree,1.88*degree,38.5*degree,0.203e-5,0.032e-3);//nufit 4.0
    glbSetDensityParams(input_errors,0.05,GLB_ALL);
    glbSetCentralValues(centers);
    glbSetInputErrors(input_errors);
 

    glbRegisterPriorFunction(FASE_prior_OSC,NULL,NULL,NULL);
    
    PARA=MODEL;
    glbSetOscillationParameters(true_values);
    glbSetRates();
    glbSetProjection(free);
    float res0=glbChiNP(true_values,NULL,GLB_ALL);
    glbSetProjection(projection);
    

    /*if the user gives the true values in oscillation parameters
    PARA=STAN;
    glbSetOscillationParameters(true_values);
    glbSetRates(); PARA=MODEL;
    glbSetProjection(free);
    float res0=glbChiNP(true_values,NULL,GLB_ALL);
    glbSetProjection(projection);
    */

    
    for (x=lower_x;x<=upper_x;x=x+dx){ glbSetOscParams(test_values,x,0);
    for (eta=lower_eta;eta<=upper_eta;eta=eta+deta){ glbSetOscParams(test_values,eta,1);
        float res=glbChiNP(test_values,out,GLB_ALL);
        float r_out   = glbGetOscParams(out,2);
        float ma_out  = glbGetOscParams(out,3);
//        printf("%f %f %f %f %f\n",x,eta/M_PI,res-res0,r_out,ma_out);
        fprintf(File,"%f %f %f %f %f\n",x,eta/M_PI,res-res0,r_out,ma_out);
         } fprintf(File,"\n");}
    

    
    /* Destroy parameter and projection vector(s) */
    glbFreeParams(true_values);
    glbFreeParams(test_values); 
    glbFreeParams(input_errors); 
    glbFreeParams(out);
    glbFreeProjection(free);
    glbFreeProjection(projection);
    exit(0);
}

