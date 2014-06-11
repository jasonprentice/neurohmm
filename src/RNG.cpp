//
//  RNG.cpp
//  
//
//  Created by Jason Prentice on 6/11/14.
//
//

#include <vector>
#include <ctime>

#include "RNG.h"

using namespace std;

RNG::RNG() {
    // Initialize mersenne twister RNG
    rng_pr = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng_pr, time(NULL));
}
RNG::~RNG() {
    gsl_rng_free(rng_pr);
}
int RNG::discrete(const vector<double>& p) const {
    double u = gsl_rng_uniform(rng_pr);
    double c = p[0];
    int ix=0;
    while (c<u) {
        ix++;
        c += p[ix];
    }
    return ix;
}
bool RNG::bernoulli(double p) const {
    return (gsl_rng_uniform(rng_pr) < p);
}
double RNG::uniform() const {
    return (gsl_rng_uniform(rng_pr));
}

vector<int> RNG::randperm(int nmax) const {
    vector<int> nvals (nmax);
    for (int i=0; i<nmax; i++) {
        nvals[i] = i;
    }
    for (int i=0; i<nmax; i++) {
        // select random integer ix between i and nmax-1
        // swap i with ix
        int ix = i + gsl_rng_uniform_int(rng_pr, nmax-i);
        int tmp = nvals[i];
        nvals[i] = nvals[ix];
        nvals[ix] = tmp;
    }
    return nvals;
}
