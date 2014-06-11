//
//  IndependentBasin.cpp
//  
//
//  Created by Jason Prentice on 6/11/14.
//
//


#include "IndependentBasin.h"

using namespace std;


// IndependentBasin


IndependentBasin::IndependentBasin(int N, int basin_num, RNG* rng, double min) : BasinModel(N,basin_num,rng), above_thresh_bool(N,0), prefactor(1)
{
    stats.assign(N, 0);
    for (vector<double>::iterator it = stats.begin(); it != stats.end(); ++it) {
        //        double u = 0.1*((double) rand() / (double) RAND_MAX) + min;
        double u = 0.1*rng->uniform() + min;
        //        double u = 0.1*((double) rand(time) / (double) RAND_MAX) + 0.45;
        (*it) = u;
    }
    m.assign(stats,N,1);
    update_thresh_list();
    
}

vector<int> IndependentBasin::get_active_constraints(const State& this_state) {
    return this_state.on_neurons;
}


void IndependentBasin::doMLE(double alpha) {
    m.assign(stats, N, 1);
    update_thresh_list();
    return;
}

void IndependentBasin::update_thresh_list() {
    // We want to compute prod_i [m_i^sigma_i * (1-m_i)^(1-sigma_i)].
    // We want to take advantage of the fact that most sigma_i=0, and not do
    // the product over all i for every state: instead write it as
    // prod_i (1-m_i) prod_i (m_i/(1-m_i))^sigma_i.
    // However, if m_i ~ 1 this risks underflow so we only use
    // this trick for m_i < thresh.
    
    double thresh = 0.9;
    prefactor = 1;
    above_thresh_list.clear();
    above_thresh_bool.assign(N,0);
    for (int i=0; i<N; i++) {
        if (m.at(i) < thresh) {
            prefactor *= (1 - m.at(i));
        } else {
            above_thresh_list.push_back(i);
            above_thresh_bool[i] = 1;
        }
    }
    
    return;
}

double IndependentBasin::P_state(const State& this_state) const {
    
    double P = prefactor;
    for (vector<int>::const_iterator neuron_iter = this_state.on_neurons.begin(); neuron_iter != this_state.on_neurons.end(); ++neuron_iter) {
        if (above_thresh_bool[*neuron_iter] == 0) {
            P *= (m.at(*neuron_iter) / (1-m.at(*neuron_iter)));
        }
    }
    
    for (vector<int>::const_iterator thresh_iter=above_thresh_list.begin(); thresh_iter!=above_thresh_list.end(); ++thresh_iter) {
        if (this_state.word[*thresh_iter]==1) {
            P *= m.at(*thresh_iter);
        } else {
            P *= (1 - m.at(*thresh_iter));
        }
    }
    
    return P;
}

vector<char> IndependentBasin::sample() const {
    vector<char> this_sample (N);
    for (int i=0; i<N; i++) {
        this_sample[i] = (rng->bernoulli(m.at(i))) ? 1 : 0;
    }
    return this_sample;
}

paramsStruct IndependentBasin::get_params() {
    
    //    myMatrix m (stats, N, 1);
    paramsStruct params;
    params.addField("m", m);
    return params;
    
}
//int IndependentBasin::nparams() const {
//    return N;
//}
//
//vector<double> IndependentBasin::get_params() const {
//    return m;
//}
//
