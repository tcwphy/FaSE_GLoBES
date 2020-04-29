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
/* Macros */
#define SQR(x)          ((x)*(x))                       /*x^2*/
#define SQR_ABS(x)      (SQR(creal(x))+SQR(cimag(x)))   /*|x|^2*/

int N_M;
typedef struct OSC_PARAMS { double theta12; double theta13; double theta23; double delta; double alpha21; double alpha31; double Dm21; double Dm31; double m3; double m2; double m1; } OSC_PARAMS;


double MODEL_init(int N) /*initialise the model-input*/
{
    if (N>6){printf("the number of parameters cannot be larger than 6.\n");}
    N_M=N;

    return 0;
}



int STAN_OSC(double complex M[], double out[6])
{
/*https://github.com/daviddoria/Examples/blob/master/c%2B%2B/GSL/Test1/Test1.cpp
 http://linux.math.tifr.res.in/manuals/html/gsl-ref-html/gsl-ref_14.html
 */

    double complex MMsquare[] = { 1.0 + 0.0*I, 10.0 + 0.0*I, 3.0 + 0.0*I,
        0.0 + 0.0*I, 2.0 + 0.0*I, 100.0 + 0.0*I,
        1.0 + 0.0*I, 5.0 + 0.0*I, 3.0 + 0.0*I };
    
    
    double complex beta = 0.0 + 0.0*I;
    double complex alpha = 1.0 + 0.0*I;
    
    cblas_zgemm (CblasRowMajor, CblasConjTrans, CblasNoTrans, 3, 3, 3, &alpha, M, 3, M, 3, &beta, MMsquare, 3); /*eq 4.5 page 107 Giunti*/

    
    
    
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

//    double complex U[];
  /*
    gsl_complex M00 = gsl_matrix_complex_get(MM,0,0);
    gsl_complex M01 = gsl_matrix_complex_get(MM,0,1);
    gsl_complex M02 = gsl_matrix_complex_get(MM,0,2);
    */
//    printf(" gsl %g %g %g %g %g %g \n",GSL_REAL(M00),GSL_IMAG(M00),GSL_REAL(M01),GSL_IMAG(M01),GSL_REAL(M02),GSL_IMAG(M02));
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
//    double dCP=carg(-U[0]*conj(U[2])*conj(U[3])*U[5]-c12*c12*s13*s13*c13*c13*s23*s23);
    double dCP=carg(U[1]*conj(U[2])*conj(U[4])*U[5]+s12*s12*s13*s13*c13*c13*s23*s23);
//    double dCP=carg(conj(U[1])*U[2]*U[4]*conj(U[5])+s12*s12*s13*s13*c13*c13*s23*s23);
    /* check if it is correct: I think the phase defined by gsl is opposite to the PMNS.
     i.e. the matrix by gsl in dual to PMNS */
    
/*    double absJ=cabs(U[1]*conj(U[2])*conj(U[4])*U[5]+s12*s12*s13*s13*c13*c13*s23*s23);
    double sind = cimag(U[1]*conj(U[2])*conj(U[4])*U[5]+s12*s12*s13*s13*c13*c13*s23*s23)/absJ;
    double cosd = creal(U[1]*conj(U[2])*conj(U[4])*U[5]+s12*s12*s13*s13*c13*c13*s23*s23)/absJ;
    dCP=acos(cosd);
    if(sind<0) dCP=2*M_PI-dCP;
  */
    /*
     gsl_complex M00 = gsl_matrix_complex_get(evec,0,0);
     gsl_complex M01 = gsl_matrix_complex_get(evec,0,1);
     gsl_complex M02 = gsl_matrix_complex_get(evec,0,2);
     */
    //    printf(" gsl %g %g %g %g %g %g \n",GSL_REAL(M00),GSL_IMAG(M00),GSL_REAL(M01),GSL_IMAG(M01),GSL_REAL(M02),GSL_IMAG(M02));
    
    
//    printf("dCP %g %g \n",sind,cosd);
    double m1=   sqrt(gsl_vector_get(eval,0));
    double m2=   sqrt(gsl_vector_get(eval,1));
    double m3=   sqrt(gsl_vector_get(eval,2));

/*    double m1=   gsl_vector_get(eval,0);
    double m2=   gsl_vector_get(eval,1);
    double m3=   gsl_vector_get(eval,2);*/

    
//    printf("my %g %g %g\n",m1,m2,m3);
    
    out[0]=the12; out[1]=the13; out[2]=the23;
    out[3]=dCP;   out[4]=m2*m2-m1*m1;
    out[5]=m3*m3-m1*m1;
    
//    printf(" out %g %g %g %g %g %g\n",out[0]*180/M_PI,out[1]*180/M_PI,out[2]*180/M_PI,out[3]*180/M_PI,out[4],out[5]);
    return 0;
}



int ModelTO( double OSC_PARAMS[6],double M_para[])
{
    
    
    //for the 2 heavy sterile case keep these zero.
    double alpha = 0.0;
    double mc = 0.0;
    double gamma = 0.0;
    
    /*
    double N=M_para[0];
    double ma=M_para[1];
    double mb=M_para[2];
    double beta=M_para[3];
    */
    /*     double complex Mass_Matrix[] = { mb*(cos(beta) + I*sin(beta)), N*mb*(cos(beta) + I*sin(beta)), (N -2)*mb*(cos(beta) + I*sin(beta)),
     N*mb*(cos(beta) + I*sin(beta)), ma*(cos(alpha) + I*sin(alpha)) + N*N*mb*(cos(beta) + I*sin(beta)), ma*(cos(alpha) + I*sin(alpha)) + N*(N-2)*mb*(cos(beta) + I*sin(beta)),
     (N-2)*mb*(cos(beta) + I*sin(beta)), ma*(cos(alpha) + I*sin(alpha)) + N*(N-2)*mb*(cos(beta) + I*sin(beta)), ma*(cos(alpha) + I*sin(alpha)) + (N-2)*(N-2)*mb*(cos(beta) + I*sin(beta)) + mc*(cos(gamma)*I*sin(gamma)) };
     */
    /*
    double complex Mass_Matrix[] = { 1,0,0,0,1,0,0,0,1};
    */
    
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


double MtoS(double osc_para[6], double M_para[])
/*revert model parameters to
 the standard oscillation parameters */
{
    
/*  double x=M_para[0];
    double eta=M_para[1];
    double r=M_para[2];
    double ma=M_para[3];
    printf("x %g; eta %g pi; r %g; ma %g\n",x,eta/M_PI,r,ma);
    ModelTO(osc_para,M_para);*/
  
    double x=M_para[0];
    double eta=M_para[1];
    double r=M_para[2];
    double ma=M_para[3];
    
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
