//
//  Autocorr.h
//  
//
//  Created by Jason Prentice on 6/11/14.
//
//

#ifndef ____Autocorr__
#define ____Autocorr__

#include <vector>
#include "HMM.h"

// ************ Autocorr ***********
template <class BasinT>
class Autocorr : public HMM<BasinT>
{
public:
    Autocorr(std::vector<std::vector<double> >& st, double binsize, int nbasins);
    ~Autocorr();
    
    std::vector<double> train(int niter);
    std::vector<int> viterbi();
    std::vector<double> get_basin_trans();
    
private:
    
    void update_forward();
    void update_backward();
    void update_P();
    void update_w();
    void update_basin_trans_indep();
    
    std::vector<double> trans_at_t(int);
    
    std::vector<double*> basin_trans;
    double logli();
};

#endif /* defined(____Autocorr__) */
