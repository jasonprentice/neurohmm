//
//  RNG.h
//  
//
//  Created by Jason Prentice on 6/11/14.
//
//  Wrapper around the Gnu Scientific Library RNG

#ifndef ____RNG__
#define ____RNG__

#include <vector>
#include <gsl/gsl_rng.h>

// ************ RNG ***************
class RNG {
public:
    RNG();
    ~RNG();
    int                     discrete( const std::vector<double>& ) const;
    bool                    bernoulli( double ) const;
    double                  uniform() const;
    std::vector<int>        randperm( int ) const;
private:
    gsl_rng *               rng_pr;
};


#endif /* defined(____RNG__) */
