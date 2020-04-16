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
 Some simple output functions for GLoBES examples
 */


//#ifndef TRI_DIRECT_PE_H
//#define TRI_DIRECT_PE_H


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


/*#define GLB_NU_FLAVOURS 3*/
//#define GLB_SIGMA_E 6
//#define GLB_ETA 0
//#define GLB_MA 1
//#define GLB_MB 2
//#define STAN 0
//#define MODEL -100
//#define SQR(x)      ((x)*(x))                        /* x^2   */
//#define SQR_ABS(x)  (SQR(creal(x)) + SQR(cimag(x)))  /* |x|^2 */
//#define GLB_NU_FLAVOURS 3
//#define GLB_V_FACTOR        7.5e-14   /* Conversion factor for matter potentials */

//#define GLB_Ne_MANTLE       0.5        /* Effective electron numbers for calculation */
//#define GLB_Ne_CORE         0.468      /*   of MSW potentials                        */
//#define V_THRESHOLD 0.001*GLB_V_FACTOR*GLB_Ne_MANTLE  /* The minimum matter potential below */
/* which vacuum algorithms are used   */
//#define M_SQRT3  1.73205080756887729352744634151      /* sqrt(3) */


/* Interface for non-standard implementations of the probability engine */
typedef struct OSC_PARAMS { double theta12; double theta13; double theta23; double delta; double alpha21; double alpha31; double Dm21; double Dm31; double m3; double m2; double m1; } OSC_PARAMS;

int N_M;
double MODEL_init();
int STAN_OSC(double complex M[], OSC_PARAMS output[6]);
int ModelTO( double OSC_PARAMS[6],double M_para[]);
double MtoS(double osc_para[6], double M_para[]);
double model_restriction(double model []);
