//
//  BasinModel.cpp
//  
//
//  Created by Jason Prentice on 11/22/12.
//
//


#include "RNG.h"
#include "EMBasins.h"
#include "BasinModel.h"

using namespace std;


void BasinModel::reset_stats() {
    norm = 0;
    for (vector<double>::iterator it=stats.begin(); it!=stats.end(); ++it) {
        *it = 0;
    }
    return;
}

template <class BasinT>
void BasinModel::increment_stats(const State<typename EMBasins<BasinT>::BasinData>& this_state) {
//    double wt = this_state.freq * this_state.P[basin_num];
    double wt = this_state.aux_data.weight[basin_num];
    norm += wt;
    for (vector<int>::const_iterator it=this_state.aux_data.active_constraints.begin(); it!=this_state.aux_data.active_constraints.end(); ++it) {
        stats[*it] += wt;
    }
    return;
}

void BasinModel::normalize_stats() {
    for (vector<double>::iterator it=stats.begin(); it!=stats.end(); ++it) {
        *it /= norm;
    }
    return;
}

//int BasinModel::nparams() const {
//    return stats.size();
//}
//vector<double> BasinModel::get_params() const {
//    return stats;
//}
