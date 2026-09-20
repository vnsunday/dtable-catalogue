#include <stdlib.h>
#include <stdio.h>
#include <cmath>

#define PI  3.14159265359

double gauss_function(double x, double m, double sigma) {
    return std::exp( - (x - m) * (x - m) / (2 * sigma * sigma)   )  /(sigma * sqrt (2* PI)) ;
}

double bernoulli_function(int k, double p) {
    return 0.0;
}

/*============================================================
 * Output: 
 *      (mean1, mean2): Confidence interval
 *============================================================*/
double mean_estimation(double x[], int n, double sigma, double& mean1, double& mean2) {
    return 0.0;
}

double statistic_pick_intervals(double x[], int n) {
    return 0.0;
}

int main() {
    return 0;
}