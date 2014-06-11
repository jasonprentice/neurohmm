//
//  IndependentBasin.h
//  
//
//  Created by Jason Prentice on 6/11/14.
//
//

#ifndef ____IndependentBasin__
#define ____IndependentBasin__

#include <vector>

#include "BasinModel.h"
#include "myMatrix.h"

class paramsStruct;

// ************************** IndependentBasin ************************

class IndependentBasin : public BasinModel {
public:
    IndependentBasin( int, int, RNG*, double );
    static std::vector<int> get_active_constraints( const State& );
    void                    doMLE( double );
    double                  P_state( const State& ) const;
    std::vector<char>       sample() const;
    paramsStruct            get_params();
private:
    myMatrix<double>        m;
    std::vector<char>       above_thresh_bool;
    std::vector<int>        above_thresh_list;
    double                  prefactor;
    
    void                    update_thresh_list();
};
// ***************************************************************


#endif /* defined(____IndependentBasin__) */
