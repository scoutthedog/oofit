#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "levmar.h"
#include "oofit.h"

// This file constructs implimentations of the levmar functions
// that can be used for fitting dose-response data
// Levenberg-Marquandt algorithm

// range of possible doses, mostly for reference / testing
// double dose[19] = {-10,-9.53,-9,-8.53,-8,-7.53,-7,-6.53,-6,-5.53,-5,-4.53,-4,-3.53,-3,-2.53,-2,-1.53,-1};

static void drc1(double *p, double *x, int m, int n, void *adata) {
// pass the dose array as a void pointer
// convert the void pointer to a double pointer
double * ptr = (double *) adata;
register int i;
for (i=0; i<n; ++i) {
    // get the value of the thing that the pointer is pointing to as d
    // double d = *ptr;
    double d = *ptr;
    x[i] = 100/(1+(pow(10, (p[0]-d)*p[1])));
    ++ptr; // increment the pointer so it points to the next item in the array
    // this is a really bad/weird way of doing this but I wrote it when I
    // was still learning how to use pointers ... don't judge me
    }
}
static void drc2(double *p, double *x, int m, int n, void *adata) {
double *ptr = (double *) adata;
register int i;
for (i=0; i<n; ++i) {
    double d = *ptr;
    x[i] = p[2] + ((100-p[2])/(1+(pow(10,(p[0]-d)*p[1]))));
    ++ptr;
    }
}
static void drc3(double *p, double *x, int m, int n, void *adata) {
double *ptr = (double *) adata;
register int i;
for (i=0; i<n; ++i) {
    double d = *ptr;
    x[i] = p[2]/(1+(pow(10, (p[0]-d)*p[1])));
    ++ptr;
    }
}
static void drc4(double *p, double *x, int m, int n, void *adata) {
double *ptr = (double *) adata;
register int i;
for (i=0; i <n; ++i) {
    double d = *ptr;
    x[i] = p[2]+((p[3]-p[2])/(1+(pow(10, (p[0]-d)*p[1]))));
    ++ptr;
    }
}
static double ec50est(double * res, double * dose) {
    // approximates the EC50/IC50 value by finding the 
    // response that is closes to 50%, then returning
    // the corresponding dose. I think this method works
    // well enough for data where finding the EC50 value
    // is reasonably possible
    int i;
    int n = sizeof(*dose);
    double d50[n];
    for (i=0; i<n; ++i) {
      d50[i] = abs(*(res + i) -  50.0);
      //printf("%.7g \n", d50[i]); debugging
    }
    double big = 100;
    for(i = 0; i < n; ++i) {
        if(big > d50[i])
        big = d50[i];
    }
    for (i=0; i<n; ++i) {
        if (d50[i] == big)
            return *(dose + i); 
        }
    printf("Warning: reached unreachable code");
    return 1;
}
static int slope(double * res) {
    // used for finding an initial value for the hillslope
    // Basically just tries to figure out if the curve is mostly
    // positive or negative in slope and returns 1 or -1
    int i;
    int n = sizeof(*res);
    double inc;
    for (i=0; i<n-1; ++i) {
        // find the differences between each of the points along the curve
        // sum them up
        inc += (*(res + i + 1)) - (*(res + i));
     }
    // if the sum is positive, return 1
    if (inc > 0) {
        return 1;
     }
    // if negative return -1
    if (inc < 0) {
        return -1;
    } else {
        return 1;
     }
}

double * fit1 (double * res, double * dose, double *lb, double *ub, double * ret_arr) {
     // This is the two parameter "contrained" version of fitting.
    // keeps ymin set to 0 and ymax set to 100
    // you can specify your own lower and upper bounds
    int m=2;
    int n = sizeof(res)-1;
    double p[m];
     double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
    opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
    opts[4]=LM_DIFF_DELTA; // for finite differences
    p[0] = ec50est(res, dose);
    p[1] = slope(res);
    int ret = dlevmar_bc_dif(drc1, p, res, m, n, lb, ub, NULL, 1000, opts, info, NULL, NULL, (void *) dose);
    ret_arr[0] = p[0]; ret_arr[1] = p[1]; ret_arr[2] = 0;
    ret_arr[3] = 100; ret_arr[4] = 1; ret_arr[5] = (double) ret;
    return ret_arr;
}
double * fit1d (double * res, double * dose, double * ret_arr) {
    // defualt value version 
    double lb[2] = {-DBL_MAX, -DBL_MAX};
    double ub[2] = {DBL_MAX, DBL_MAX};
    return fit1(res, dose, lb, ub, ret_arr);
}
double * fit2 (double * res, double *dose, double *lb, double *ub, double * ret_arr) {
    int m=3;
    int n = sizeof(res)-1;
    double p[m];
    double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
    opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
    opts[4]=LM_DIFF_DELTA; // for finite differences
    p[0] = ec50est(res, dose);
    p[1] = slope(res);
     p[2] = 0.0;
    int ret = dlevmar_bc_dif(drc2, p, res, m, n, lb, ub, NULL, 1000, opts, info, NULL, NULL, (void *) dose);
    ret_arr[0] = p[0]; ret_arr[1] = p[1]; ret_arr[2] = p[2];
    ret_arr[3] = 100; ret_arr[4] = 2; ret_arr[5] = (double) ret;
    return ret_arr;
}
double * fit2d (double * res, double * dose, double * ret_arr) {
      double lb[3] = {-DBL_MAX, -DBL_MAX, 0};
    double ub[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
    return fit2(res, dose, lb, ub, ret_arr);
}
double * fit3 (double * res, double *dose, double *lb, double *ub, double * ret_arr) {
    int m=3;
    int n = sizeof(res)-1;
    double p[m];
    double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
    opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
    opts[4]=LM_DIFF_DELTA; // for finite differences
    p[0] = ec50est(res, dose);
    p[1] = slope(res);
     p[2] = 100.0;
    int ret = dlevmar_bc_dif(drc3, p, res, m, n, lb, ub, NULL, 1000, opts, info, NULL, NULL, (void *) dose);
    ret_arr[0] = p[0]; ret_arr[1] = p[1]; ret_arr[2] = 0;
    ret_arr[3] = p[2]; ret_arr[4] = 3; ret_arr[5] = (double) ret;
    return ret_arr;
}
double * fit3d (double * res, double * dose, double * ret_arr) {
    double lb[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
    double ub[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
    return fit3(res, dose, lb, ub, ret_arr);
}
double * fit4 (double * res, double *dose, double *lb, double *ub, double * ret_arr) {
    int m=4;
    int n = sizeof(res)-1;
    double p[m];
    double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
    opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
    opts[4]=LM_DIFF_DELTA; // for finite differences
    p[0] = ec50est(res, dose);
    p[1] = slope(dose);
    p[2] = 0.0;
    p[3] = 100.0;
    int ret = dlevmar_bc_dif(drc4, p, res, m, n, lb, ub, NULL, 1000, opts, info, NULL, NULL, (void *) dose);
    ret_arr[0] = p[0]; ret_arr[1] = p[1]; ret_arr[2] = p[2];
    ret_arr[3] = p[3]; ret_arr[4] = 4; ret_arr[5] = (double) ret;
    return ret_arr;

}
double * fit4d (double * res, double * dose, double * ret_arr) {
    double lb[4] = {-DBL_MAX, -DBL_MAX, 0, -DBL_MAX};
    double ub[4] = {DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX};
    return fit4(res, dose, lb, ub, ret_arr);
}
void pfit (double * ret_arr) {
    printf("EasyFit---\nEC50      : %.7g\nhillslope : %.7g\nymin      : %.7g\nymax      : %.7g\nequation  : %.7g\nconv      : %.7g\n",
		  ret_arr[0], ret_arr[1], ret_arr[2], ret_arr[3], ret_arr[4], ret_arr[5]);
}
/*
int main() {
    double res_arr[7] = {20.232,30.2,40.2,50.2,70.1,95.2,100};
    double dose_arr[7] = {-7.52,-7,-6.52,-6,-5.52,-5,-4};
    double ret_arr[6];
    fit1d(res_arr, dose_arr, ret_arr);
    pfit(ret_arr);
    fit2d(res_arr, dose_arr, ret_arr);
    pfit(ret_arr);
    fit3d(res_arr, dose_arr, ret_arr);
    pfit(ret_arr);
    fit4d(res_arr, dose_arr, ret_arr);
    pfit(ret_arr);
  
    exit(0);
}
*/
