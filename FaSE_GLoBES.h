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
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include "model-input.h"

/*#define GLB_NU_FLAVOURS 3*/
#define GLB_SIGMA_E 6        /* Index of non-standard parameter sigma_E */
#define GLB_ETA 0
#define GLB_MA 1
#define GLB_MB 2
#define STAN 0
#define MODEL -100
#define DUNE_nu 0
#define DUNE_anu 1
#define DUNE_app 0
#define DUNE_disapp 1
#define T2HK 2
#define RAM 100
#define    top_bin_number 10000
#define SQR(x)      ((x)*(x))                        /* x^2   */
#define SQR_ABS(x)  (SQR(creal(x)) + SQR(cimag(x)))  /* |x|^2 */
#define GLB_NU_FLAVOURS 3
#define GLB_V_FACTOR        7.5e-14   /* Conversion factor for matter potentials */
//#define GLB_Ne_MANTLE       0.497      /* Effective electron numbers for calculation */
#define GLB_Ne_MANTLE       0.5        /* Effective electron numbers for calculation */
#define GLB_Ne_CORE         0.468      /*   of MSW potentials                        */
#define V_THRESHOLD 0.001*GLB_V_FACTOR*GLB_Ne_MANTLE  /* The minimum matter potential below */
/* which vacuum algorithms are used   */
#define M_SQRT3  1.73205080756887729352744634151      /* sqrt(3) */

/* Interface for non-standard implementations of the probability engine */
double th12;
double th13;
double th23;
double deltacp;
double sdm;
double ldm;
double sigma_E;
double delta;
int N_in;
static double mq[3];
int    PARA,PARA_in,ran,ran_in,n_ram,FIT,PRIOR;
double x_in,eta_in,r_in,ma_in;

double UPPER_prior[6],LOWER_prior[6],Central_prior[6];
int STAN_OSC(double complex M[], double output[6]);
int STAN_OSC_U(double complex U[], double out[6]);
inline double square(double x);
int sinsq(double complex U[], double * a);
int glb_free_probability_engine();
int glb_init_probability_engine();
int zheevc3(double complex A[3][3], double w[3]);
void zhetrd3(double complex A[3][3], double complex Q[3][3],
             double d[3], double e[2]);
int zheevq3(double complex A[3][3], double complex Q[3][3], double w[3]);
int FASE_glb_set_oscillation_parameters(glb_params p, void *user_data);
int FASE_glb_get_oscillation_parameters(glb_params p, void *user_data);
int glb_hamiltonian_cd(double E, double V, int cp_sign);
int glb_S_matrix_cd(double E, double L, double V, int cp_sign);
int glb_filtered_probability_matrix_cd(double P[3][3], double E, double L, double V,
                                       double sigma, int cp_sign);
int FASE_glb_probability_matrix(double P[3][3], int cp_sign, double E,
                              int psteps, const double *length, const double *density,
                              double filter_sigma, void *user_data);
int DUNE_random_true_n(int m);
int extra_random_true_n(int m);
int random_truerate(int m, int exp, int rule);
inline double glb_likelihood2(double true_rate, double fit_rate);
inline double glb_prior2(double x, double center, double sigma);
double glbChiSpectrumTilt2(int exp, int rule, int n_params, double *x, double *errors,
                           void *user_data);
double glbChiSpectrumTilt3(int exp, int rule, int n_params, double *x, double *errors,
                           void *user_data);
double FASE_prior_OSC(const glb_params in, void* user_data);
double FASE_prior_model(const glb_params in, void* user_data);
