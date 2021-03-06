/*--------------------------------------------------------------------
 This example is using littlest seesaw model for tri-direct approach
 as an example. More details can be checked in arXiv:1304.6264[hep-ph],
 arXiv:1512.07531[hep-ph] and arXiv:1607.05276[hep-ph].
 The user using this model, should cite the references:
 
 @article{King:2016yvg,
 author = "King, Stephen F. and Luhn, Christoph",
 title = "{Littlest Seesaw model from S$_{4} \times$ U(1)}",
 eprint = "1607.05276",
 archivePrefix = "arXiv",
 primaryClass = "hep-ph",
 reportNumber = "QFET-2016-12, SI-HEP-2016-20",
 doi = "10.1007/JHEP09(2016)023",
 journal = "JHEP",
 volume = "09",
 pages = "023",
 year = "2016"
 }
 
 @article{King:2015dvf,
 author = "King, Stephen F.",
 title = "{Littlest Seesaw}",
 eprint = "1512.07531",
 archivePrefix = "arXiv",
 primaryClass = "hep-ph",
 doi = "10.1007/JHEP02(2016)085",
 journal = "JHEP",
 volume = "02",
 pages = "085",
 year = "2016"
 }
 
 @article{King:2013iva,
 author = "King, Stephen F.",
 title = "{Minimal predictive see-saw model with normal neutrino mass hierarchy}",
 eprint = "1304.6264",
 archivePrefix = "arXiv",
 primaryClass = "hep-ph",
 doi = "10.1007/JHEP07(2013)137",
 journal = "JHEP",
 volume = "07",
 pages = "137",
 year = "2013"
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

/*Define the formula (5) in the working note provided by Ding.*/
double complex TDModelY(double x, double eta, double r, double Ma){
    double Ms=Ma*r;
    
    double complex result;
    result = (5*SQR(x) + 2*x + 2)*(Ma + cexp(I*eta)*Ms)/(SQR(x) + x +1)*1/2;
    
    return result;
    
} //Checked!

/*Define the formula (5) in the working note provided by Ding.*/
double complex TDModelZ(double x, double eta, double r, double Ma){
    double Ms=Ma*r;
    
    double complex result;
    result = (-1)*sqrt(5*SQR(x) + 2*x + 2)/(SQR(x) + x + 1) * 1/2 * ((x+2)*Ma - x*(2*x+1) * cexp(I*eta)*Ms);
    
    return result;
} //Checked!

/*Define the formula (5) in the working note provided by Ding.*/
double complex TDModelW(double x, double eta, double r, double Ma){
    double Ms=Ma*r;
    
    double complex result;
    result = 1/(SQR(x) + x + 1) * 1/2 * ( SQR(x+2)*Ma + SQR(x)*SQR(2*x+1) * cexp(I*eta)*Ms);
    
    return result;
} //Checked!

/*Define the formula (7) in the working note provided by Ding.*/
double TDpsi(double x, double eta, double r, double Ma){
    double complex y, z, w;
    y = TDModelY(x, eta, r, Ma);
    z = TDModelZ(x, eta, r, Ma);
    w = TDModelW(x, eta, r, Ma);
    
    double complex tmp;
    double sinus, result;
    tmp = conj(y)*z + w*conj(z) ;
    sinus = cimag( tmp ) / cabs( tmp );
    
    result = asin(sinus);
    
    return result;
} //Checked!

/*Define the formula (7) in the working note provided by Ding.*/
double TDtheta(double x, double eta, double r, double Ma){
    double complex y, z, w;
    y = TDModelY(x, eta, r, Ma);
    z = TDModelZ(x, eta, r, Ma);
    w = TDModelW(x, eta, r, Ma);
    
    double norm, denorm;
    double sinus, result;
    norm = 2 * cabs( conj(y)*z + w*conj(z) ); //Cross checked!
    denorm = sqrt( SQR( SQR_ABS(w) - SQR_ABS(y) ) + 4*SQR_ABS( conj(y) * z + w * conj(z) ) );
    
    sinus =  norm/denorm ;
    
    result = asin(sinus)/2;
    
    return result;
} //checked!

/*Now we come to define the neutrino mixing angles */
double TDth12(double x, double eta, double r, double Ma){
    double angle, sinus, result;
    angle = TDtheta(x, eta, r, Ma);
    sinus = 1 - 3*SQR(x)/  (3*SQR(x) + 2 * (SQR(x) + x + 1) * SQR(cos(angle)) ) ;
    result = asin( sqrt(sinus) );
    
    return result;
    
} //checked!

double TDth13(double x, double eta, double r, double Ma){
    double angle, sinus, result;
    angle = TDtheta(x, eta, r, Ma);
    sinus = 2*( SQR(x) + x + 1 ) * SQR(sin(angle))  /  (5*SQR(x) + 2*x + 2) ;
    result = asin( sqrt(sinus) );
    
    return result;
    
} //checked!

double TDth23(double x, double eta, double r, double Ma){
    double angle1, angle2, sinus, result; //angle1 is theta in Eq. (7), angle2 is psi in Eq. (7)
    angle1 = TDtheta(x, eta, r, Ma);
    angle2 = TDpsi(x, eta, r, Ma);
    double norm, denorm;
    norm = sin( 2*angle1 ) * sin( angle2 ) * x * sqrt( 3*(5*SQR(x) + 2*x + 2) );
    denorm = 3*SQR(x) + 2 * ( SQR(x) + x + 1 ) * SQR( cos(angle1) );
    
    sinus = 0.5+0.5*norm/denorm;
    result = asin( sqrt(sinus) );
    
    return result;
    
} //checked!

double TDcosdCP(double x, double eta, double r, double Ma)
{
    double two_theta23=2*TDth23(x,eta,r,Ma);
    double theta13=TDth13(x,eta,r,Ma);
    double cot223=cos(two_theta23)/sin(two_theta23);
    double sin13=sin(theta13);
    double cos13_square=cos(theta13)*cos(theta13);
    double cosdelta;
    cosdelta=cot223*(3*SQR(x)-(4*SQR(x)+x+1)*cos13_square);
    cosdelta=cosdelta/(sqrt(3)*fabs(x)*sin13*sqrt((5*x*x+2*x+2)*cos13_square-3*x*x));
    
    return cosdelta;
    
}


double TDdCP(double x, double eta, double r, double Ma){
    double judgeAngle = TDpsi(x, eta, r, Ma);
    double theta23 = TDth23(x, eta, r, Ma);
    double theta13 = TDth13(x, eta, r, Ma);
    double rightSide, norm, denorm, result, sinCP;
    norm = SQR(cos(theta13)/sin(theta13)) * SQR(cos(2*theta23)) * SQR( SQR(x) + x + 1 );
    denorm = 3*SQR(x)*( 3*SQR(x)*SQR(tan(theta13)) - 2*(SQR(x) + x + 1) );
    rightSide = sqrt(1+norm/denorm)/sin(2*theta23);
    
    int sign;
    
    if (x*cos(judgeAngle) > 0 ){
        sign = +1;
        sinCP = sign * rightSide;
    }
    else{
        sign = -1;
        sinCP = sign * rightSide;
    }
    
    //  return result = asin(sinCP);
    if (TDcosdCP(x,eta,r,Ma)<0) {result = M_PI-asin(sinCP)-2*M_PI;}
    else    result = asin(sinCP);
    return result;
} //checked!

double TDdm21(double x, double eta, double r, double Ma){
    double complex y, z, w;
    y = TDModelY(x, eta, r, Ma);
    z = TDModelZ(x, eta, r, Ma);
    w = TDModelW(x, eta, r, Ma);
    
    double part1, part2, result;
    part1 = SQR_ABS(y) + SQR_ABS(w) + 2*SQR_ABS(z);
    part2 = sqrt( SQR(SQR_ABS(w)-SQR_ABS(y)) + 4*SQR_ABS(conj(y)*z+w*conj(z)) );
    result = 0.5 * ( part1 - part2);
    
    return result;
    
} //checked!

double TDdm31(double x, double eta, double r, double Ma){
    double complex y, z, w;
    y = TDModelY(x, eta, r, Ma);
    z = TDModelZ(x, eta, r, Ma);
    w = TDModelW(x, eta, r, Ma);
    
    double part1, part2, result;
    part1 = SQR_ABS(y) + SQR_ABS(w) + 2*SQR_ABS(z);
    part2 = sqrt( SQR(SQR_ABS(w)-SQR_ABS(y)) + 4*SQR_ABS(conj(y)*z+w*conj(z)) );
    result = 0.5 * ( part1 + part2);
    
    return result;
    
} //checked!



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
    /*example: tri-direct*/
    
    double x   = M_para[0];
    double eta = M_para[1];
    double r   = M_para[2];
    double ma  = M_para[3];
    
    /*theta12*/        osc_para[0]=TDth12(x,eta,r, ma);
    /*theta13*/        osc_para[1]=TDth13(x,eta,r, ma);
    /*theta23*/        osc_para[2]=TDth23(x,eta,r, ma);
    /*deltaCP*/        osc_para[3]= TDdCP(x,eta,r, ma);
    /*Delta_m^2_{21}*/ osc_para[4]=TDdm21(x,eta,r, ma);
    /*Delta_m^2_{31}*/ osc_para[5]=TDdm31(x,eta,r, ma);

    /*
    double ms=ma*r;
    
    double complex Mass_Matrix[] = {ma+ms*(cos(eta) + I*sin(eta)), ma*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))+x*ms*(cos(eta) + I*sin(eta)), ma*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))+x*ms*(cos(eta) + I*sin(eta)), ma*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))+x*ms*(cos(eta) + I*sin(eta)), ma*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))+x*x*ms*(cos(eta) + I*sin(eta)), ma+x*x*ms*(cos(eta) + I*sin(eta)), ma*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))+x*ms*(cos(eta) + I*sin(eta)), ma+x*x*ms*(cos(eta) + I*sin(eta)), ma*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))+x*x*ms*(cos(eta) + I*sin(eta))};
    
    double out[6];
    STAN_OSC(Mass_Matrix,osc_para);
    */
    return 0;
}

double model_restriction(double model [])
{
    
    double x=model[0];
    double eta=model[1];
    double r=model[2];
    double ma=model[3];
    
    if(ma<0) {return 1;}
    if(r<0.) {return 1;}
    return 0;
}

