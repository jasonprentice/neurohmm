//
//  BasinModel.h
//  
//
//  Created by Jason Prentice on 11/22/12.
//
//

#ifndef ____BasinModel__
#define ____BasinModel__

#include <vector>

#include "myMatrix.h"
#include "paramsStruct.h"
#include "RNG.h"
#include "EMBasins.h"


class BasinData;

class BasinModel
{
public:
                            BasinModel( int N, int basin_num, RNG* rng )
                                : N(N), basin_num(basin_num), rng(rng) {};
    void                    reset_stats();

    template <class BasinT>
    void                    increment_stats( const State<typename EMBasins<BasinT>::BasinData>& );
    void                    normalize_stats();
    inline double           get_norm() const {
                                return norm; };
protected:
    std::vector<double>     stats;
    double                  norm;
    int                     N;
    int                     basin_num;
    RNG *                   rng;
};


#endif /* defined(____BasinModel__) */
