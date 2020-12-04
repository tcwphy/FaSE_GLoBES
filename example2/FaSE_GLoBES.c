/*Probability engine for a given flavour symmetry,
 which needs to be given in model-input.c. This code
 using output oscillation parameter values from
 model-input.c, compute the probability spectra and
 the prior value. This code is based on the default
 GLoBES probability engine.
 ---------------------------------------------- */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <float.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <globes/globes.h>   /* GLoBES library */
#include"FaSE_GLoBES.h"
#include"model-input.h"

extern char **glb_param_names;

typedef struct {
    double *osc_params;
    size_t length;
} glb_osc_type;

typedef struct {
    double *density_params;
    size_t length;
} glb_density_type;

struct glb_params_type {
    glb_osc_type *osc;
    glb_density_type *density;
    int iterations;
};

typedef struct {
    int *osc_params;
    size_t length;
} glb_osc_proj_type;

typedef struct {
    int *density_params;
    size_t length;
} glb_density_proj_type;

struct glb_projection_type {
    glb_osc_proj_type *osc;
    glb_density_proj_type *density;
};

//typedef struct OSC_PARAMS { double theta12; double theta13; double theta23; double delta; double alpha21; double alpha31; double Dm21; double Dm31; double m3; double m2; double m1; } OSC_PARAMS;

/* Interface for non-standard implementations of the probability engine */
int glb_oscp;
glb_probability_matrix_function glb_hook_probability_matrix;
glb_set_oscillation_parameters_function glb_hook_set_oscillation_parameters;
glb_get_oscillation_parameters_function glb_hook_get_oscillation_parameters;
void *glb_probability_user_data=NULL;
/* Internal temporary variables */
gsl_matrix_complex *U=NULL; /* The vacuum mixing matrix                           */
gsl_matrix_complex *H=NULL; /* Neutrino Hamiltonian                               */
gsl_matrix_complex *Q=NULL; /* Eigenvectors of Hamiltonian (= eff. mixing matrix) */
gsl_vector *lambda=NULL;    /* Eigenvalues of Hamiltonian                         */
gsl_matrix_complex *S=NULL; /* The neutrino S-matrix                              */
gsl_matrix_complex *H0_template=NULL;  /* Used in the construction of the vac. Hamiltonian */
gsl_matrix_complex *S1=NULL, *T0=NULL; /* Temporary matrix storage                         */


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
    out[3]=dCP;   out[4]=gsl_vector_get(eval,1)-gsl_vector_get(eval,0);
    out[5]=gsl_vector_get(eval,2)-gsl_vector_get(eval,0);
    gsl_vector_free (eval);
    gsl_matrix_complex_free (evec);
    gsl_matrix_complex_free (MM);
    gsl_eigen_hermv_free(w);
    return 0;
}


int STAN_OSC_U(double complex U[], double out[4])
{
    
    double s13 = cabs(U[2]);       double the13=asin(s13);
    double t12 = cabs(U[1])/cabs(U[0]);           double the12=atan(t12);
    double t23 = cabs(U[5])/cabs(U[8]);           double the23=atan(t23);
    
    double s12=sin(the12);
    double c12=cos(the12);
    double s23=sin(the23);
    double c13=cos(the13);
    double dCP=carg(U[1]*conj(U[2])*conj(U[4])*U[5]+s12*s12*s13*s13*c13*c13*s23*s23);
    
    out[0]=the12; out[1]=the13; out[2]=the23;
    out[3]=dCP;
    
    return 0;
}



inline double square(double x)
{
    return x*x;
}

/***************************************************************************
 * Function glb_free_probability_engine                                    *
 ***************************************************************************
 * Destroys internal data structures of the probability engine.            *
 ***************************************************************************/
int glb_free_probability_engine()
{
    if (T0!=NULL)     { gsl_matrix_complex_free(T0);  T0 = NULL; }
    if (S1!=NULL)     { gsl_matrix_complex_free(S1);  S1 = NULL; }
    if (H0_template!=NULL) { gsl_matrix_complex_free(H0_template);  H0_template = NULL; }
    
    if (S!=NULL)      { gsl_matrix_complex_free(S);   S = NULL; }
    if (lambda!=NULL) { gsl_vector_free(lambda);      lambda = NULL; }
    if (Q!=NULL)      { gsl_matrix_complex_free(Q);   Q = NULL; }
    if (H!=NULL)      { gsl_matrix_complex_free(H);   H = NULL; }
    if (U!=NULL)      { gsl_matrix_complex_free(U);   U = NULL; }
    
    return 0;
}
/***************************************************************************
 * Function glb_init_probability_engine                                    *
 ***************************************************************************
 * Allocates internal data structures for the probability engine.          *
 ***************************************************************************/
int glb_init_probability_engine()
{
    glb_free_probability_engine();
    
    U = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
    H = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
    Q = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
    lambda = gsl_vector_alloc(GLB_NU_FLAVOURS);
    S = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
    
    H0_template = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
    S1 = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
    T0 = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
    
    return 0;
}

int zheevc3(double complex A[3][3], double w[3])
// ----------------------------------------------------------------------------
// Calculates the eigenvalues of a hermitian 3x3 matrix A using Cardano's
// analytical algorithm.
// Only the diagonal and upper triangular parts of A are accessed. The access
// is read-only.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The hermitian input matrix
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error
// ----------------------------------------------------------------------------
{
    double m, c1, c0;
    
    // Determine coefficients of characteristic poynomial. We write
    //       | a   d   f  |
    //  A =  | d*  b   e  |
    //       | f*  e*  c  |
    double complex de = A[0][1] * A[1][2];                            // d * e
    double dd = SQR_ABS(A[0][1]);                                  // d * conj(d)
    double ee = SQR_ABS(A[1][2]);                                  // e * conj(e)
    double ff = SQR_ABS(A[0][2]);                                  // f * conj(f)
    m  = creal(A[0][0]) + creal(A[1][1]) + creal(A[2][2]);
    c1 = (creal(A[0][0])*creal(A[1][1])  // a*b + a*c + b*c - d*conj(d) - e*conj(e) - f*conj(f)
          + creal(A[0][0])*creal(A[2][2])
          + creal(A[1][1])*creal(A[2][2]))
    - (dd + ee + ff);
    c0 = creal(A[2][2])*dd + creal(A[0][0])*ee + creal(A[1][1])*ff
    - creal(A[0][0])*creal(A[1][1])*creal(A[2][2])
    - 2.0 * (creal(A[0][2])*creal(de) + cimag(A[0][2])*cimag(de));
    // c*d*conj(d) + a*e*conj(e) + b*f*conj(f) - a*b*c - 2*Re(conj(f)*d*e)
    
    double p, sqrt_p, q, c, s, phi;
    p = SQR(m) - 3.0*c1;
    q = m*(p - (3.0/2.0)*c1) - (27.0/2.0)*c0;
    sqrt_p = sqrt(fabs(p));
    
    phi = 27.0 * ( 0.25*SQR(c1)*(p - c1) + c0*(q + 27.0/4.0*c0));
    phi = (1.0/3.0) * atan2(sqrt(fabs(phi)), q);
    
    c = sqrt_p*cos(phi);
    s = (1.0/M_SQRT3)*sqrt_p*sin(phi);
    
    w[1]  = (1.0/3.0)*(m - c);
    w[2]  = w[1] + s;
    w[0]  = w[1] + c;
    w[1] -= s;
    
    return 0;
}


void zhetrd3(double complex A[3][3], double complex Q[3][3],
             double d[3], double e[2])
// ----------------------------------------------------------------------------
// Reduces a hermitian 3x3 matrix to real tridiagonal form by applying
// (unitary) Householder transformations:
//            [ d[0]  e[0]       ]
//    A = Q . [ e[0]  d[1]  e[1] ] . Q^T
//            [       e[1]  d[2] ]
// The function accesses only the diagonal and upper triangular parts of
// A. The access is read-only.
// ---------------------------------------------------------------------------
{
    const int n = 3;
    double complex u[n], q[n];
    double complex omega, f;
    double K, h, g;
    int i,j;
    // Initialize Q to the identitity matrix
#ifndef EVALS_ONLY
    for (i=0; i < n; i++)
    {
        Q[i][i] = 1.0;
        for (j=0; j < i; j++)
            Q[i][j] = Q[j][i] = 0.0;
    }
#endif
    
    // Bring first row and column to the desired form
    h = SQR_ABS(A[0][1]) + SQR_ABS(A[0][2]);
    if (creal(A[0][1]) > 0)
        g = -sqrt(h);
    else
        g = sqrt(h);
    e[0] = g;
    f    = g * A[0][1];
    u[1] = conj(A[0][1]) - g;
    u[2] = conj(A[0][2]);
    
    omega = h - f;
    if (creal(omega) > 0.0)
    {
        omega = 0.5 * (1.0 + conj(omega)/omega) / creal(omega);
        K = 0.0;
        for (i=1; i < n; i++)
        {
            f    = conj(A[1][i]) * u[1] + A[i][2] * u[2];
            q[i] = omega * f;                  // p
            K   += creal(conj(u[i]) * f);      // u* A u
        }
        K *= 0.5 * SQR_ABS(omega);
        
        for (i=1; i < n; i++)
            q[i] = q[i] - K * u[i];
        
        d[0] = creal(A[0][0]);
        d[1] = creal(A[1][1]) - 2.0*creal(q[1]*conj(u[1]));
        d[2] = creal(A[2][2]) - 2.0*creal(q[2]*conj(u[2]));
        
        // Store inverse Householder transformation in Q
#ifndef EVALS_ONLY
        for (j=1; j < n; j++)
        {
            f = omega * conj(u[j]);
            for (i=1; i < n; i++)
                Q[i][j] = Q[i][j] - f*u[i];
        }
#endif
        
        // Calculate updated A[1][2] and store it in f
        f = A[1][2] - q[1]*conj(u[2]) - u[1]*conj(q[2]);
    }
    else
    {
        for (i=0; i < n; i++)
            d[i] = creal(A[i][i]);
        f = A[1][2];
    }
    
    // Make (23) element real
    e[1] = cabs(f);
#ifndef EVALS_ONLY
    if (e[1] != 0.0)
    {
        f = conj(f) / e[1];
        for (i=1; i < n; i++)
            Q[i][n-1] = Q[i][n-1] * f;
    }
#endif
}



int zheevq3(double complex A[3][3], double complex Q[3][3], double w[3])
// ----------------------------------------------------------------------------
// Calculates the eigenvalues and normalized eigenvectors of a hermitian 3x3
// matrix A using the QL algorithm with implicit shifts, preceded by a
// Householder reduction to real tridiagonal form.
// The function accesses only the diagonal and upper triangular parts of A.
// The access is read-only.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The hermitian input matrix
//   Q: Storage buffer for eigenvectors
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error (no convergence)
// ----------------------------------------------------------------------------
// Dependencies:
//   zhetrd3()
// ----------------------------------------------------------------------------
{
    const int n = 3;
    double e[3];                 // The third element is used only as temporary workspace
    double g, r, p, f, b, s, c;  // Intermediate storage
    double complex t;
    int nIter;
    int m,l,i,k;
    
    // Transform A to real tridiagonal form by the Householder method
    zhetrd3(A, Q, w, e);
    
    // Calculate eigensystem of the remaining real symmetric tridiagonal matrix
    // with the QL method
    //
    // Loop over all off-diagonal elements
    for (l=0; l < n-1; l++)
    {
        nIter = 0;
        while (1)
        {
            // Check for convergence and exit iteration loop if off-diagonal
            // element e(l) is zero
            for (m=l; m <= n-2; m++)
            {
                g = fabs(w[m])+fabs(w[m+1]);
                if (fabs(e[m]) + g == g)
                    break;
            }
            if (m == l)
                break;
            
            if (nIter++ >= 30)
                return -1;
            
            // Calculate g = d_m - k
            g = (w[l+1] - w[l]) / (e[l] + e[l]);
            r = sqrt(SQR(g) + 1.0);
            if (g > 0)
                g = w[m] - w[l] + e[l]/(g + r);
            else
                g = w[m] - w[l] + e[l]/(g - r);
            
            s = c = 1.0;
            p = 0.0;
            for (i=m-1; i >= l; i--)
            {
                f = s * e[i];
                b = c * e[i];
                if (fabs(f) > fabs(g))
                {
                    c      = g / f;
                    r      = sqrt(SQR(c) + 1.0);
                    e[i+1] = f * r;
                    c     *= (s = 1.0/r);
                }
                else
                {
                    s      = f / g;
                    r      = sqrt(SQR(s) + 1.0);
                    e[i+1] = g * r;
                    s     *= (c = 1.0/r);
                }
                
                g = w[i+1] - p;
                r = (w[i] - g)*s + 2.0*c*b;
                p = s * r;
                w[i+1] = g + p;
                g = c*r - b;
                
                // Form eigenvectors
#ifndef EVALS_ONLY
                for (k=0; k < n; k++)
                {
                    t = Q[k][i+1];
                    Q[k][i+1] = s*Q[k][i] + c*t;
                    Q[k][i]   = c*Q[k][i] - s*t;
                }
#endif
            }
            w[l] -= p;
            e[l]  = g;
            e[m]  = 0.0;
        }
    }
    
    return 0;
}

int zheevh3(double complex A[3][3], double complex Q[3][3], double w[3])
// ----------------------------------------------------------------------------
// Calculates the eigenvalues and normalized eigenvectors of a hermitian 3x3
// matrix A using Cardano's method for the eigenvalues and an analytical
// method based on vector cross products for the eigenvectors. However,
// if conditions are such that a large error in the results is to be
// expected, the routine falls back to using the slower, but more
// accurate QL algorithm. Only the diagonal and upper triangular parts of A need
// to contain meaningful values. Access to A is read-only.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The hermitian input matrix
//   Q: Storage buffer for eigenvectors
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error
// ----------------------------------------------------------------------------
// Dependencies:
//   zheevc3(), zhetrd3(), zheevq3()
// ----------------------------------------------------------------------------
// Version history:
//   v1.1: Simplified fallback condition --> speed-up
//   v1.0: First released version
// ----------------------------------------------------------------------------
{
#ifndef EVALS_ONLY
    double norm;          // Squared norm or inverse norm of current eigenvector
    //  double n0, n1;        // Norm of first and second columns of A
    double error;         // Estimated maximum roundoff error
    double t, u;          // Intermediate storage
    int j;                // Loop counter
#endif
    
    // Calculate eigenvalues
    zheevc3(A, w);
    
#ifndef EVALS_ONLY
    //  n0 = SQR(creal(A[0][0])) + SQR_ABS(A[0][1]) + SQR_ABS(A[0][2]);
    //  n1 = SQR_ABS(A[0][1]) + SQR(creal(A[1][1])) + SQR_ABS(A[1][2]);
    
    t = fabs(w[0]);
    if ((u=fabs(w[1])) > t)
        t = u;
    if ((u=fabs(w[2])) > t)
        t = u;
    if (t < 1.0)
        u = t;
    else
        u = SQR(t);
    error = 256.0 * DBL_EPSILON * SQR(u);
    //  error = 256.0 * DBL_EPSILON * (n0 + u) * (n1 + u);
    
    Q[0][1] = A[0][1]*A[1][2] - A[0][2]*creal(A[1][1]);
    Q[1][1] = A[0][2]*conj(A[0][1]) - A[1][2]*creal(A[0][0]);
    Q[2][1] = SQR_ABS(A[0][1]);
    
    // Calculate first eigenvector by the formula
    //   v[0] = conj( (A - w[0]).e1 x (A - w[0]).e2 )
    Q[0][0] = Q[0][1] + A[0][2]*w[0];
    Q[1][0] = Q[1][1] + A[1][2]*w[0];
    Q[2][0] = (creal(A[0][0]) - w[0]) * (creal(A[1][1]) - w[0]) - Q[2][1];
    norm    = SQR_ABS(Q[0][0]) + SQR_ABS(Q[1][0]) + SQR(creal(Q[2][0]));
    
    // If vectors are nearly linearly dependent, or if there might have
    // been large cancellations in the calculation of A(I,I) - W(1), fall
    // back to QL algorithm
    // Note that this simultaneously ensures that multiple eigenvalues do
    // not cause problems: If W(1) = W(2), then A - W(1) * I has rank 1,
    // i.e. all columns of A - W(1) * I are linearly dependent.
    if (norm <= error)
        return zheevq3(A, Q, w);
    else                      // This is the standard branch
    {
        norm = sqrt(1.0 / norm);
        for (j=0; j < 3; j++)
            Q[j][0] = Q[j][0] * norm;
    }
    
    // Calculate second eigenvector by the formula
    //   v[1] = conj( (A - w[1]).e1 x (A - w[1]).e2 )
    Q[0][1]  = Q[0][1] + A[0][2]*w[1];
    Q[1][1]  = Q[1][1] + A[1][2]*w[1];
    Q[2][1]  = (creal(A[0][0]) - w[1]) * (creal(A[1][1]) - w[1]) - creal(Q[2][1]);
    norm     = SQR_ABS(Q[0][1]) + SQR_ABS(Q[1][1]) + SQR(creal(Q[2][1]));
    if (norm <= error)
        return zheevq3(A, Q, w);
    else
    {
        norm = sqrt(1.0 / norm);
        for (j=0; j < 3; j++)
            Q[j][1] = Q[j][1] * norm;
    }
    
    // Calculate third eigenvector according to
    //   v[2] = conj(v[0] x v[1])
    Q[0][2] = conj(Q[1][0]*Q[2][1] - Q[2][0]*Q[1][1]);
    Q[1][2] = conj(Q[2][0]*Q[0][1] - Q[0][0]*Q[2][1]);
    Q[2][2] = conj(Q[0][0]*Q[1][1] - Q[1][0]*Q[0][1]);
#endif
    
    return 0;
}




int FASE_glb_set_oscillation_parameters(glb_params p, void *user_data)  /*Mark2*/
{
    double complex (*_U)[GLB_NU_FLAVOURS]
    = (double complex (*)[GLB_NU_FLAVOURS]) gsl_matrix_complex_ptr(U, 0, 0);
    int i;
    /*printf("set\n");*/
    /* Copy parameters */
    /*PARA_in=p->osc->osc_params[5];*/
    /*if (PARA_in==CSDN)*/
    if (PARA==MODEL)
        
    {
        double osc_para[6];
        double M_para[N_M];
        
        for(i=0;i<N_M;i++) M_para[i] = p->osc->osc_params[i];
        
        //        printf("in MODEL x %g eta %g r %g ma %g\n",x_in,eta_in,r_in,ma_in);
        //      find_CSND(N_in,ma_in,mb_in,eta_in,OSC_PARAMS_CSDN);
        
        MtoS(osc_para, M_para);
        
        th12  = osc_para[0];
        th13  = osc_para[1];
        th23  = osc_para[2];
        delta = osc_para[3];
        mq[0] = fabs(osc_para[5]);
        mq[1] = fabs(osc_para[5])+osc_para[4];
        mq[2] = fabs(osc_para[5])+osc_para[5];
//        printf("in CSDN %g %g %g %g %g %g\n",th12*180/M_PI,th13*180/M_PI,th23*180/M_PI,delta*180/M_PI,mq[1]-mq[0],mq[2]-mq[0]);
        /*printf("in CSDN %g %g %g %g %g %g\n",th12*180/M_PI,th13*180/M_PI,th23*180/M_PI,delta*180/M_PI,mq[1]-mq[0],mq[2]-mq[0]);*/
        
    }
    else if (PARA==STAN)
        
    {
        th12  = p->osc->osc_params[0];
        th13  = p->osc->osc_params[1];
        th23  = p->osc->osc_params[2];
        delta = p->osc->osc_params[3];
        mq[0] = fabs(p->osc->osc_params[5]);
        mq[1] = fabs(p->osc->osc_params[5]) + p->osc->osc_params[4];
        mq[2] = fabs(p->osc->osc_params[5]) + p->osc->osc_params[5];
        //        printf("in PMNS %g %g %g %g %g %g\n",th12,th13,th23,delta,mq[1]-mq[0],mq[2]-mq[0]);
        
    }
    
    //        printf("in MODEL x %g eta %g r %g ma %g\n",x_in,eta_in,r_in,ma_in);
    //        printf("in PMNS %g %g %g %g %g %g\n",th12,th13,th23,delta,mq[1]-mq[0],mq[2]-mq[0]);
    //        printf("%g %g %g %g %g %g %g %g %g %g\n",x_in,eta_in,r_in,ma_in,th12,th13,th23,delta,mq[1]-mq[0],mq[2]-mq[0]);
    /* Compute vacuum mixing matrix */
    _U[0][0] = cos(th12)*cos(th13);
    _U[0][1] = sin(th12)*cos(th13);
    _U[0][2] = sin(th13) * cexp(-I * delta);
    
    _U[1][0] = -sin(th12)*cos(th23) - cos(th12)*sin(th23)*sin(th13) * cexp(I*delta);
    _U[1][1] =  cos(th12)*cos(th23) - sin(th12)*sin(th23)*sin(th13) * cexp(I*delta);
    _U[1][2] =  sin(th23)*cos(th13);
    
    _U[2][0] =  sin(th12)*sin(th23) - cos(th12)*cos(th23)*sin(th13) * cexp(I*delta);
    _U[2][1] = -cos(th12)*sin(th23) - sin(th12)*cos(th23)*sin(th13) * cexp(I*delta);
    _U[2][2] =  cos(th23)*cos(th13);
    
    /* Calculate energy independent matrix H0 * E */
    gsl_matrix_complex_set_zero(H0_template);
    gsl_matrix_complex_set_zero(H);
    for (i=0; i < GLB_NU_FLAVOURS; i++)
        gsl_matrix_complex_set(H0_template, i, i, gsl_complex_rect(0.5*mq[i], 0.0));
    
    gsl_matrix_complex *T = gsl_matrix_complex_alloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
    gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, GSL_COMPLEX_ONE, H0_template, U, /* T=H0.U^\dagger */
                   GSL_COMPLEX_ZERO, T);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, U, T,             /* H0=U.T */
                   GSL_COMPLEX_ZERO, H0_template);
    gsl_matrix_complex_free(T);
    return 0;
}

int FASE_glb_get_oscillation_parameters(glb_params p, void *user_data) /*Mark2*/
{
    /*printf("get\n");*/
    if (PARA==STAN) glbDefineParams(p, th12, th13, th23, delta, mq[1] - mq[0], mq[2] - mq[0]); /*to read the corresponding oscaillation parameters*/
    else if (PARA==MODEL) glbDefineParams(p, x_in, eta_in, r_in, ma_in, 0, 0);
    
    return 0;
}


int glb_hamiltonian_cd(double E, double V, int cp_sign)
{
    double inv_E = 1.0 / E;
    double complex (*_H)[GLB_NU_FLAVOURS]
    = (double complex (*)[GLB_NU_FLAVOURS]) gsl_matrix_complex_ptr(H, 0, 0);
    double complex (*_H0_template)[GLB_NU_FLAVOURS]
    = (double complex (*)[GLB_NU_FLAVOURS]) gsl_matrix_complex_ptr(H0_template, 0, 0);
    int i, j;
    
    if (cp_sign > 0)
    {
        for (i=0; i < GLB_NU_FLAVOURS; i++)
            for (j=0; j < GLB_NU_FLAVOURS; j++)
                _H[i][j] = _H0_template[i][j] * inv_E;
    }
    else
    {
        for (i=0; i < GLB_NU_FLAVOURS; i++)
            for (j=0; j < GLB_NU_FLAVOURS; j++)
                _H[i][j] = conj(_H0_template[i][j] * inv_E); /* delta_CP -> -delta_CP */
    }
    
    _H[0][0] = _H[0][0] + cp_sign*V;
    return 0;
}


int glb_S_matrix_cd(double E, double L, double V, int cp_sign)
{
    /* Introduce some abbreviations */
    double complex (*_S)[3]  = (double complex (*)[3]) gsl_matrix_complex_ptr(S,0,0);
    double complex (*_Q)[3]  = (double complex (*)[3]) gsl_matrix_complex_ptr(Q,0,0);
    double complex (*_T0)[3] = (double complex (*)[3]) gsl_matrix_complex_ptr(T0,0,0);
    double *_lambda = gsl_vector_ptr(lambda,0);
    int status;
    int i, j, k;
    double inv_E;
    
    if (fabs(V) < V_THRESHOLD)                       /* Vacuum */
    {
        /* Use vacuum mixing angles and masses */
        inv_E = 0.5/E;
        for (i=0; i < GLB_NU_FLAVOURS; i++)
            _lambda[i] = mq[i] * inv_E;
        
        if (cp_sign > 0)
            gsl_matrix_complex_memcpy(Q, U);
        else
        {
            double complex (*_U)[3]  = (double complex (*)[3]) gsl_matrix_complex_ptr(U,0,0);
            for (i=0; i < GLB_NU_FLAVOURS; i++)
                for (j=0; j < GLB_NU_FLAVOURS; j++)
                    _Q[i][j] = conj(_U[i][j]);
        }
    }
    else                                             /* Matter */
    {
        double complex (*_H)[3] = (double complex (*)[3]) gsl_matrix_complex_ptr(H,0,0);
        
        /* Calculate neutrino Hamiltonian */
        if ((status=glb_hamiltonian_cd(E, V, cp_sign)) != 0)
            return status;
        
        /* Calculate eigenvalues of Hamiltonian */
        if ((status=zheevh3(_H, _Q, _lambda)) != 0)
            return status;
    }
    
    /* Calculate S-Matrix in mass basis in matter ... */
    double phase;
    gsl_matrix_complex_set_zero(S);
    for (i=0; i < GLB_NU_FLAVOURS; i++)
    {
        phase    = -L * _lambda[i];
        _S[i][i] = cos(phase) + I*sin(phase);
    }
    
    /* ... and transform it to the flavour basis */
    gsl_matrix_complex_set_zero(T0);
    double complex *p = &_T0[0][0];
    
    
    for (i=0; i < GLB_NU_FLAVOURS; i++)              /* T0 = S.Q^\dagger */
        for (j=0; j < GLB_NU_FLAVOURS; j++)
        {
            for (k=0; k < GLB_NU_FLAVOURS; k++)
            {
                *p += ( creal(_S[i][k])*creal(_Q[j][k])+cimag(_S[i][k])*cimag(_Q[j][k]) )
                + I * ( cimag(_S[i][k])*creal(_Q[j][k])-creal(_S[i][k])*cimag(_Q[j][k]) );
            }
            p++;
        }
    gsl_matrix_complex_set_zero(S);
    p = &_S[0][0];
    for (i=0; i < GLB_NU_FLAVOURS; i++)              /* S = Q.T0 */
        for (j=0; j < GLB_NU_FLAVOURS; j++)
        {
            for (k=0; k < GLB_NU_FLAVOURS; k++)
            {
                *p += ( creal(_Q[i][k])*creal(_T0[k][j])-cimag(_Q[i][k])*cimag(_T0[k][j]) )
                + I * ( cimag(_Q[i][k])*creal(_T0[k][j])+creal(_Q[i][k])*cimag(_T0[k][j]) );
            }
            p++;
        }
    
    return 0;
}

int glb_filtered_probability_matrix_cd(double P[3][3], double E, double L, double V,
                                       double sigma, int cp_sign)
{
    /* Introduce some abbreviations */
    double complex (*_Q)[3]  = (double complex (*)[3]) gsl_matrix_complex_ptr(Q,0,0);
    double complex (*_T0)[3] = (double complex (*)[3]) gsl_matrix_complex_ptr(T0,0,0);
    double *_lambda = gsl_vector_ptr(lambda,0);
    int status;
    int i, j, k, l;
    
    if (fabs(V) < V_THRESHOLD)                       /* Vacuum */
    {
        /* Use vacuum mixing angles and masses */
        double inv_E = 0.5/E;
        for (i=0; i < GLB_NU_FLAVOURS; i++)
            _lambda[i] = mq[i] * inv_E;
        
        if (cp_sign > 0)
            gsl_matrix_complex_memcpy(Q, U);
        else
        {
            double complex (*_U)[3]  = (double complex (*)[3]) gsl_matrix_complex_ptr(U,0,0);
            for (i=0; i < GLB_NU_FLAVOURS; i++)
                for (j=0; j < GLB_NU_FLAVOURS; j++)
                    _Q[i][j] = conj(_U[i][j]);
        }
    }
    else                                             /* Matter */
    {
        double complex (*_H)[3] = (double complex (*)[3]) gsl_matrix_complex_ptr(H,0,0);
        
        /* Calculate neutrino Hamiltonian */
        if ((status=glb_hamiltonian_cd(E, V, cp_sign)) != 0)
            return status;
        
        /* Calculate eigenvalues of Hamiltonian */
        if ((status=zheevh3(_H, _Q, _lambda)) != 0)
            return status;
    }
    
    // Calculate probability matrix (see GLoBES manual for a discussion of the algorithm)
    double phase, filter_factor;
    double t = -0.5/1.0e-18 * SQR(sigma) / SQR(E);
    gsl_matrix_complex_set_zero(T0);
    for (i=0; i < GLB_NU_FLAVOURS; i++)
        for (j=i+1; j < GLB_NU_FLAVOURS; j++)
        {
            phase         = -L * (_lambda[i] - _lambda[j]);
            filter_factor = exp(t * SQR(phase));
            _T0[i][j]     = filter_factor * (cos(phase) + I*sin(phase));
        }
    
    for (k=0; k < GLB_NU_FLAVOURS; k++)
        for (l=0; l < GLB_NU_FLAVOURS; l++)
        {
            P[k][l] = 0.0;
            for (i=0; i < GLB_NU_FLAVOURS; i++)
            {
                /*complex t = conj(_Q[k][i]) * _Q[l][i];*/ /*not very sure about this*/
                t = conj(_Q[k][i]) * _Q[l][i];
                for (j=i+1; j < GLB_NU_FLAVOURS; j++)
                    P[k][l] += 2.0 * creal(_Q[k][j] * conj(_Q[l][j]) * t * _T0[i][j]);
                P[k][l] += SQR_ABS(_Q[k][i]) * SQR_ABS(_Q[l][i]);
            }
        }
    
    return 0;
}



int FASE_glb_probability_matrix(double P[3][3], int cp_sign, double E,
                              int psteps, const double *length, const double *density,
                              double filter_sigma, void *user_data) 
{
    int status;
    int i, j;
    
    /* Convert energy to eV */
    E *= 1.0e9;
    
    if (filter_sigma > 0.0)                     /* With low-pass filter */
    {
        if (psteps == 1)
            glb_filtered_probability_matrix_cd(P, E, GLB_KM_TO_EV(length[0]), density[0]*GLB_V_FACTOR*GLB_Ne_MANTLE,
                                               filter_sigma, cp_sign);
        else
            return -1;
    }
    else                                        /* Without low-pass filter */
    {
        if (psteps > 1)
        {
            gsl_matrix_complex_set_identity(S1);                                 /* S1 = 1 */
            for (i=0; i < psteps; i++)
            {
                status = glb_S_matrix_cd(E, GLB_KM_TO_EV(length[i]), density[i]*GLB_V_FACTOR*GLB_Ne_MANTLE, cp_sign);
                if (status != 0)
                    return status;
                gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, S, S1, /* T0 = S.S1 */
                               GSL_COMPLEX_ZERO, T0);
                gsl_matrix_complex_memcpy(S1, T0);                                 /* S1 = T0 */
            }
            gsl_matrix_complex_memcpy(S, S1);                                    /* S = S1 */
        }
        else
        {
            status = glb_S_matrix_cd(E, GLB_KM_TO_EV(length[0]), density[0]*GLB_V_FACTOR*GLB_Ne_MANTLE, cp_sign);
            if (status != 0)
                return status;
        }
        
        double complex (*_S)[3] = (double complex (*)[3]) gsl_matrix_complex_ptr(S,0,0);
        for (i=0; i < GLB_NU_FLAVOURS; i++)
            for (j=0; j < GLB_NU_FLAVOURS; j++)
                P[j][i] = SQR_ABS(_S[i][j]);
    }
    
    return 0;
}



double FASE_prior_OSC(const glb_params in, void* user_data) /*the prior is always used in the standard parameterisation*/
{
    glb_projection p = glbAllocProjection();
    glb_params input_errors = glbAllocParams();
    glbGetInputErrors(input_errors);
    int i;
    double central_STAN[6];
    double pv = 0.0;
    double fitvalue,centralvalue,inputerror;
    double osc_para[6];
    double M_para[N_M];
//         FILE* File_in=fopen("data/for_prior.dat", "w");
    /* Add oscillation parameter priors */
    if(PARA==MODEL){
/*
        x   = glbGetOscParams(in,0);
        eta = glbGetOscParams(in,1);
        r   = glbGetOscParams(in,2);
        ma  = glbGetOscParams(in,3);*/
        
//        for(i=0;i<N_M;i++) M_para[i] = in->osc->osc_params[i];
        for(i=0;i<N_M;i++) M_para[i] = glbGetOscParams(in,i);
//        printf("in OSC %g %g %g %g\n",M_para[0],M_para[1],M_para[2],M_para[3]);
 
        MtoS(osc_para, M_para);
//         printf("in OSC %g %g %g %g %g %g\n",osc_para[0],osc_para[1],osc_para[2],osc_para[3],osc_para[4],osc_para[5]);
//        printf("")
            if (model_restriction(M_para)==1) return 1e8;
        
       }
    else{for(i=0;i<6;i++) osc_para[i] = glbGetOscParams(in,i);}

    for(i=0;i<6;i++){
        fitvalue=osc_para[i]; centralvalue=Central_prior[i];
        if(i==3){
            if(fitvalue>2*M_PI){int run; run=fitvalue/2/M_PI; fitvalue=fitvalue-run*2*M_PI;}
            else if(fitvalue<0){int run; run=-fitvalue/2/M_PI; run=run+1; fitvalue=fitvalue+run*2*M_PI;}}
         if(fitvalue>centralvalue) {inputerror=UPPER_prior[i];} else {inputerror=LOWER_prior[i];}
           if(inputerror>1e-12) { pv+=square((centralvalue-fitvalue)/inputerror);
//        fprintf(File_in,"%i %g %g %g %g\n",i,fitvalue,centralvalue,inputerror,pv);
//        printf("in\n");
        }
    
    }
//    printf("in M pv %g\n",pv);
    /* Add matter parameter priors */
    for(i=0;i<glb_num_of_exps;i++){
        if(glbGetDensityProjectionFlag(p,i)==GLB_FREE)
        {   fitvalue=glbGetDensityParams(in,i);
            centralvalue=1.0;
            inputerror=glbGetDensityParams(input_errors,i);
            if(inputerror>1e-12)
                pv+=square((centralvalue-fitvalue)/inputerror);
        }}
    glbFreeProjection(p);
    glbFreeParams(input_errors);
//    printf("in M pv %g\n",pv);
    return pv;
}

double FASE_prior_model(const glb_params in, void* user_data)
{
    glb_params input_errors = glbAllocParams();
    glbGetInputErrors(input_errors);

    glb_projection p = glbAllocProjection();
    int i;
    double pv = 0.0;
    double fitvalue,centralvalue,inputerror;
    
    /* Add model parameter priors */
    if(PARA==STAN){ printf("The input parameters do NOT fit to the setup for prior!! \n"); return 1e6;}
    
    for(i=0;i<N_M;i++){
        fitvalue=in->osc->osc_params[i]; centralvalue=Central_prior[i];
        if(fitvalue>centralvalue) {inputerror=UPPER_prior[i];}
        else {inputerror=LOWER_prior[i];}
         if(inputerror>1e-12) { pv+=square((centralvalue-fitvalue)/inputerror);} }
    
    /* Add matter parameter priors */
    for(i=0;i<glb_num_of_exps;i++){
        if(glbGetDensityProjectionFlag(p,i)==GLB_FREE)
        {   fitvalue=glbGetDensityParams(in,i);
            centralvalue=1.0;
            inputerror=glbGetDensityParams(input_errors,i);
            if(inputerror>1e-12)
                pv+=square((centralvalue-fitvalue)/inputerror); } }
    glbFreeProjection(p);
    glbFreeParams(input_errors);
    return pv;
}
