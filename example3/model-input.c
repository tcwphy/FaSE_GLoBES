/*--------------------------------------------------------------------
 This example is using TM1 sum rule as an example.
 More details can be checked in arXiv:0812.0436[hep-ph]
 The user using this model, should cite the reference:
 
 @article{Albright:2008rp,
 author = "Albright, Carl H. and Rodejohann, Werner",
 title = "{Comparing Trimaximal Mixing and Its Variants with Deviations from Tri-bimaximal Mixing}",
 eprint = "0812.0436",
 archivePrefix = "arXiv",
 primaryClass = "hep-ph",
 reportNumber = "FERMILAB-PUB-08-535-T",
 doi = "10.1140/epjc/s10052-009-1074-3",
 journal = "Eur. Phys. J. C",
 volume = "62",
 pages = "599--608",
 year = "2009"
 }
 
 @article{Xing:2006ms,
 author = "Xing, Zhi-zhong and Zhou, Shun",
 title = "{Tri-bimaximal Neutrino Mixing and Flavor-dependent Resonant Leptogenesis}",
 eprint = "hep-ph/0607302",
 archivePrefix = "arXiv",
 doi = "10.1016/j.physletb.2007.08.009",
 journal = "Phys. Lett. B",
 volume = "653",
 pages = "278--287",
 year = "2007"
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
{// set up for normal ordering according to arXiv:1904.06660v3
    /*example: atmospheric sum rule TM1 */
    
    double th13          = M_para[0]; //
    double th23          = M_para[1];
    double sdm        = M_para[2];
    double ldm        = M_para[3];
    double Sign_sind        = M_para[4];
    
    double s13=sin(th13);
    double ssq13=s13*s13;
    double ct223=cos(2*th23)/sin(2*th23);
    double s12=sqrt(1-ssq13)/(sqrt(3)*cos(th13));
    double th12=asin(s12);
    
    double cosd=-ct223*(1-5*ssq13)/(2*sqrt(2)*s13*sqrt(1-3*ssq13));
    double dCP=acos(cosd);
    if(Sign_sind<0) dCP=-dCP;
    
    
    osc_para[0]=th12;
    osc_para[1]=th13;
    osc_para[2]=th23;
    osc_para[3]=dCP;
    osc_para[4]=sdm;
    osc_para[5]=ldm;
    return 0;
}

double model_restriction(double M_para [])
{// set up for normal ordering according to arXiv:1904.06660v3
    double th13          = M_para[0];
    double th23          = M_para[1];
    double s13=sin(th13);
    double ssq13=s13*s13;
    double ct223=cos(2*th23)/sin(2*th23);
    double s12=sqrt(1-ssq13)/(sqrt(3)*cos(th13));
    double th12=asin(s12);
    
    double cosd=-ct223*(1-5*ssq13)/(2*sqrt(2)*s13*sqrt(1-3*ssq13));
    
    if ((s12>1)||(s12<0)||(cosd>1)||(cosd<0)) return 1;
    
    return 0;
}

