#include <cmath>
#include <algorithm>
#include <iostream>

#include "RNG.h"
#include "BasinModel.h"

using namespace std;

template <class BasinT>
EMBasins<BasinT>::EMBasins(int N, int nbasins) : w(nbasins), N(N), nbasins(nbasins) {
    rng = new RNG();
    //srand(0);
}


template <class BasinT>
EMBasins<BasinT>::EMBasins(const vector<vector<double> >& st, vector<double> unobserved_l, vector<double> unobserved_u,double binsize, int nbasins)
               : rng ( new RNG() ),
                 data (st,unobserved_l,unobserved_u,binsize),
                 N (data.get_N()),
                 T (data.get_T()),
                 nbasins ( nbasins ),
                 w ( nbasins ) {}//, m ( nbasins * N ) {}


template <class BasinT>
EMBasins<BasinT>::~EMBasins() {
    delete rng;
}

template <class BasinT>
vector<double> EMBasins<BasinT>::train_model(int niter, bool ret_train_logli) {
    
    bindesc_t desc = ret_train_logli ? train : test;
    initParams();
    
    vector<double> train_logli (niter);
    for (int i=0; i<niter; i++) {
        //        cout << "Iteration " << i << endl;
        printf("Iteration %d\n", i);
        Estep();
        Mstep();
        
        //        logli[i] = update_P();
        //        test_logli[i] = update_P_test();
        train_logli[i] = logli(desc);
    }
    return train_logli;
}

template <class BasinT>
void EMBasins<BasinT>::initParams() {
    w.assign(nbasins,1/(double)nbasins);
    basins.clear();
    // Initialize each basin model
    for (int i=0; i<nbasins; i++) {
        basins.emplace_back(N,i,rng,0.45);
    }
    
    return;
}

template <class BasinT>
void EMBasins<BasinT>::Estep() {
    update_P();
    return;
}


template <class BasinT>
void EMBasins<BasinT>::Mstep() {
    for (int j=0; j<nbasins; j++) {
        basins[j].reset_stats();
    }
    
    for_each(data.begin(train), data.end(train), [this] (const State<BasinData>& state) {
        for (int j=0; j<nbasins; j++) {
            basins[j].increment_stats(state);
        } } );
             
    for (int j=0; j<nbasins; j++) {
        basins[j].normalize_stats();
    }
    
    double alpha = 0.002;
    //        if (i >= niter/2) {
    //            alpha = 0.002 + (1-0.002)*exp(-(double) (i-niter/2) * (10.0/(((double)(niter/2)-1))));
    //        }
    //        cout << alpha << endl;
    for (int j=0; j<nbasins; j++) {
        basins[j].doMLE(alpha);
    }
    update_w();
    
    
    return;
}

template <class BasinT>
void EMBasins<BasinT>::update_w() {
    for (int i=0; i<nbasins; i++) {
        w[i] = basins[i].get_norm() / T;
    }
    
    return;
}


template <class BasinT>
void EMBasins<BasinT>::update_P(bindesc_t desc) {
    for_each(data.begin(desc), data.end(desc), [this] (State<BasinData>& this_state) {
        double Z = 0;
        for (int i=0; i<nbasins; i++) {
            this_state.aux_data.P[i] = basins[i].P_state(this_state);
            Z += w[i] * this_state.aux_data.P[i];
        }
        for (int i=0; i<nbasins; i++) {
            this_state.aux_data.weight[i] = this_state.freq * w[i] * this_state.aux_data.P[i] / Z;
        }
        this_state.aux_data.pred_prob = Z;
    });
    return;
}

/*
template<class BasinT>
State<BasinData> EMBasins<BasinT>::state_obs(int obs, int t) const {
    // obs = 0:     test only
    // obs = 1:     train only
    // obs = 2:     all
    
    if (state_list[t] && (obs==1 || obs==2)) {
        return *state_list[t];
    } else if (!state_list[t] && (obs==0 || obs==2)) {
        
        State<BasinData> this_state;
        this_state.word.assign(this->N, 0);
        for (int n=0; n<this->N; n++) {
            if ((this->raster[t])[n] == '1') {
                this_state.on_neurons.push_back(n);
                this_state.word[n] = 1;
            }
        }
        
        vector<double> emiss (this->nbasins, 1);
        for (int k=0; k<this->nbasins; k++) {
            emiss[k] = (this->basins)[k].P_state(this_state);
        }
        this_state.P = emiss;
        this_state.identifier = 0;
        return this_state;
    }
    return State<BasinData>();
}

template<class BasinT>
vector<double> EMBasins<BasinT>::emiss_obs(int obs, int t) const {
    
    State this_state = state_obs(obs, t);
    if (this_state.identifier == -1) {
        return vector<double> (this->nbasins,1);
    } else {
        return this_state.P;
    }
}
*/

template <class BasinT>
double EMBasins<BasinT>::logli(bindesc_t desc) const {
    
    if (desc == test) {
        update_P(test);
    }
    
    double logli = 0;
    int nsamp = 0;
    for_each(data.begin(desc), data.end(desc), [&] (const State<BasinData>& this_state) {
        double delta = log(this_state.aux_data.pred_prob) - logli;
        nsamp++;
        logli += delta / nsamp;
    });
    
    /*
    for (int t=0; t<T; t++) {
        State this_state = state_obs(obs, t);
        if (this_state.identifier != -1) {
            if (this_state.pred_prob == -1) {
                this_state.pred_prob = 0;
                for (int k=0; k<nbasins; k++) {
                    this_state.pred_prob += w[k] * this_state.P[k];
                }
            }
            double delta = log(this_state.pred_prob) - logli;
            nsamp++;
            logli += delta / nsamp;
        }
    }
    */
    return logli;
}


template <class BasinT>
vector<double> EMBasins<BasinT>::all_prob( bindesc_t desc ) const {
    vector<double> prob (data.nstates(desc));
    transform(data.begin(desc), data.end(desc), prob.begin(), [] (const State<BasinData>& this_state) {
            return this_state.aux_data.pred_prob;
    });
    return prob;
}

template <class BasinT>
vector<double> EMBasins<BasinT>::P( bindesc_t desc ) const {
    vector<double> P (data.nstates(desc) * nbasins, 0);
    unsigned long pos = 0;
    for (auto it=data.begin(desc); it != data.end(desc); ++it) {
        for (int i=0; i<nbasins; i++) {
            P[pos*nbasins + i] = (w[i] * it->aux_data.P[i]) / it->aux_data.pred_prob;
        }
        pos++;
    }
    return P;
}

template <class BasinT>
vector<paramsStruct> EMBasins<BasinT>::basin_params() const {
    //    int nparams = basins[0].nparams();
    vector<paramsStruct> params (nbasins);
    
    for (int i=0; i<nbasins; i++) {
        params[i] = basins[i].get_params();
    }
    
    return params;
}

template <class BasinT>
vector<char> EMBasins<BasinT>::sample(int nsamples) {
    vector<char> samples (N*nsamples);
    for (int i=0; i<nsamples; i++) {
        int basin_ind = rng->discrete(w);
        vector<char> this_sample = basins[basin_ind].sample();
        for (int n=0; n<N; n++) {
            samples[i*N+n] = this_sample[n];
        }
    }
    return samples;
}



/*
 template <class BasinT>
 vector<double> EMBasins<BasinT>::test_prob() const {
 vector<double> prob (test_states.size(), 0);
 int pos = 0;
 for (const_state_iter it=test_states.begin(); it != test_states.end(); ++it) {
 const State& this_state = it->second;
 prob[pos] = this_state.pred_prob;
 pos++;
 }
 return prob;
 
 }
 */


/*
 template <class BasinT>
 vector<unsigned long> EMBasins<BasinT>::test_hist() const {
 vector<unsigned long> hist (test_states.size(), 0);
 int pos = 0;
 for (const_state_iter it=test_states.begin(); it != test_states.end(); ++it) {
 const State& this_state = it->second;
 hist[pos] = this_state.freq;
 pos++;
 }
 return hist;
 }
 */
/*
 template <class BasinT>
 vector<double> EMBasins<BasinT>::crossval(int niter, int k) {
 int blocksize = floor(raster.size() / k);
 // Generate random permutation of time bins
 vector<int> tperm = rng->randperm(raster.size());
 vector<double> all_logli (k*niter);
 for (int i=0; i<k; i++) {
 //    populate train_states and test_states
 train_states.clear();
 test_states.clear();
 for (int t=0; t<tperm.size(); t++) {
 string this_str = raster[tperm[t]];
 State this_state = train_states[this_str];
 this_state.freq = 1;
 map<string, State>& curr_map = (t < i*blocksize || t >= (i+1)*blocksize)
 ? train_states : test_states;
 pair<state_iter, bool> ins = curr_map.insert(pair<string,State> (this_str,this_state));
 if (!ins.second) {
 (((ins.first)->second).freq)++;
 }
 }
 
 //    train (returns logli of test set)
 vector<double> logli = train(niter);
 for (int j=0; j<niter; j++) {
 all_logli[i*niter + j] = logli[j];
 }
 
 }
 return all_logli;
 }
 */

