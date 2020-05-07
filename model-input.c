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
//#include <gsl/gsl_statistics.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include "FASE_GLoBES.h"
/* Macros */
//#define SQR(x)          ((x)*(x))                       /*x^2*/
//#define SQR_ABS(x)      (SQR(creal(x))+SQR(cimag(x)))   /*|x|^2*/
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

/*
int STAN_OSC(double complex M[], double out[6])
{
    //https://github.com/daviddoria/Examples/blob/master/c%2B%2B/GSL/Test1/Test1.cpp
    //http://linux.math.tifr.res.in/manuals/html/gsl-ref-html/gsl-ref_14.html
    
    
    double complex MMsquare[] = { 1.0 + 0.0*I, 10.0 + 0.0*I, 3.0 + 0.0*I,
        0.0 + 0.0*I, 2.0 + 0.0*I, 100.0 + 0.0*I,
        1.0 + 0.0*I, 5.0 + 0.0*I, 3.0 + 0.0*I };
    
    
    double complex beta = 0.0 + 0.0*I;
    double complex alpha = 1.0 + 0.0*I;
    
    cblas_zgemm (CblasRowMajor, CblasConjTrans, CblasNoTrans, 3, 3, 3, &alpha, M, 3, M, 3, &beta, MMsquare, 3); //eq 4.5 page 107 Giunti
    
    
    
    
    gsl_matrix_complex *MM = gsl_matrix_complex_alloc(3, 3);
    int i,j;
    for(i=0;i<3;i++)
    {
        for(j=0;j<3;j++)
        {
            //            gsl_matrix_complex_set(MM, i, j, gsl_complex_rect(creal(M[i*3+j]),cimag(M[i*3+j])));
            gsl_matrix_complex_set(MM, i, j, gsl_complex_rect(creal(MMsquare[i*3+j]),cimag(MMsquare[i*3+j])));
        }
    }
    
    
    gsl_eigen_hermv_workspace * w = gsl_eigen_hermv_alloc(3);
    
    gsl_vector *eval = gsl_vector_alloc(3);
    gsl_matrix_complex  *evec = gsl_matrix_complex_alloc(3, 3);
    
    gsl_eigen_hermv (MM, eval, evec, w);
    gsl_eigen_hermv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);
    
    double complex U[] = { 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I,
        0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I,
        0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I };

    int n = 0;
    for(i=0;i<3;i++)
    { for(j=0;j<3;j++) {
        gsl_vector_complex_view outcol = gsl_matrix_complex_column (evec, i);
        gsl_complex z = gsl_vector_complex_get (&outcol.vector, j);
        U[j*3+n] = GSL_REAL(z) + GSL_IMAG(z)*I;} n+=1;
    }
    
    double s13 = cabs(U[2]);       double the13=asin(s13);
    double t12 = cabs(U[1])/cabs(U[0]);           double the12=atan(t12);
    double t23 = cabs(U[5])/cabs(U[8]);           double the23=atan(t23);
    
    
    
    double s12=sin(the12);
    double c12=cos(the12);
    double s23=sin(the23);
    double c13=cos(the13);
    double dCP=carg(U[1]*conj(U[2])*conj(U[4])*U[5]+s12*s12*s13*s13*c13*c13*s23*s23);

    
    double m1=   sqrt(gsl_vector_get(eval,0));
    double m2=   sqrt(gsl_vector_get(eval,1));
    double m3=   sqrt(gsl_vector_get(eval,2));
    
    
    out[0]=the12; out[1]=the13; out[2]=the23;
    out[3]=dCP;   out[4]=m2*m2-m1*m1;
    out[5]=m3*m3-m1*m1;
    
    return 0;
}



int ModelTO( double OSC_PARAMS[6],double M_para[])
{
    
    
    //for the 2 heavy sterile case keep these zero.
    double alpha = 0.0;
    double mc = 0.0;
    double gamma = 0.0;
    
    
    
    double x=M_para[0];
    double eta=M_para[1];
    double r=M_para[2];
    double ma=M_para[3];
    
    double ms=ma*r;
    
    double complex Mass_Matrix[] = {ma+ms*(cos(eta) + I*sin(eta)), ma*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))+x*ms*(cos(eta) + I*sin(eta)), ma*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))+x*ms*(cos(eta) + I*sin(eta)), ma*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))+x*ms*(cos(eta) + I*sin(eta)), ma*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))+x*x*ms*(cos(eta) + I*sin(eta)), ma+x*x*ms*(cos(eta) + I*sin(eta)), ma*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))+x*ms*(cos(eta) + I*sin(eta)), ma+x*x*ms*(cos(eta) + I*sin(eta)), ma*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))+x*x*ms*(cos(eta) + I*sin(eta))};
    
    double out[6];
    STAN_OSC(Mass_Matrix,OSC_PARAMS);
    
    return 0;
}
 */
double MtoS(double osc_para[6], double M_para[])
/*revert model parameters to
 the standard oscillation parameters */
{
    /*example: tri-direct*/
    
    double x   = M_para[0];
    double eta = M_para[1];
    double r   = M_para[2];
    double ma  = M_para[3];
    
    /*theta12*/        //osc_para[0]=TDth12(x,eta,r, ma);
    /*theta13*/        //osc_para[1]=TDth13(x,eta,r, ma);
    /*theta23*/        //osc_para[2]=TDth23(x,eta,r, ma);
    /*deltaCP*/        //osc_para[3]= TDdCP(x,eta,r, ma);
    /*Delta_m^2_{21}*/ //osc_para[4]=TDdm21(x,eta,r, ma);
    /*Delta_m^2_{31}*/ //osc_para[5]=TDdm31(x,eta,r, ma);

    
    double ms=ma*r;
    
    double complex Mass_Matrix[] = {ma+ms*(cos(eta) + I*sin(eta)), ma*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))+x*ms*(cos(eta) + I*sin(eta)), ma*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))+x*ms*(cos(eta) + I*sin(eta)), ma*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))+x*ms*(cos(eta) + I*sin(eta)), ma*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))+x*x*ms*(cos(eta) + I*sin(eta)), ma+x*x*ms*(cos(eta) + I*sin(eta)), ma*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))+x*ms*(cos(eta) + I*sin(eta)), ma+x*x*ms*(cos(eta) + I*sin(eta)), ma*(cos(6.6666e-1*M_PI) + I*sin(6.6666e-1*M_PI))+x*x*ms*(cos(eta) + I*sin(eta))};
    
    double out[6];
    STAN_OSC(Mass_Matrix,osc_para);
    
    return 0;
}

double model_restriction(double model [])
{
    
    double x=model[0];
    double eta=model[1];
    double r=model[2];
    double ma=model[3];
    
    if(ma<0) {return 1;}
    if(r<0) {return 1;}
    
    return 0;
}

