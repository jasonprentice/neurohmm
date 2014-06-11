#include <ctime>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <iostream>

#include "RNG.h"
#include "BasinModel.h"

using namespace std;

template <class BasinT>
EMBasins<BasinT>::EMBasins(int N, int nbasins) : w(nbasins), N(N), nbasins(nbasins) {
    rng = new RNG();
    srand(time(NULL));
    //srand(0);
}


template <class BasinT>
EMBasins<BasinT>::EMBasins(vector<vector<double> >& st, vector<double> unobserved_l, vector<double> unobserved_u,double binsize, int nbasins) : w(nbasins), N( st.size() ), nbasins(nbasins), nsamples(0) {
    
    rng = new RNG();
    srand(time(NULL));
    
    //N = st.size();
    //srand(time(NULL));
    //srand(0);
    
    for (vector<double>::iterator it = unobserved_l.begin(); it != unobserved_l.end(); ++it) {
        *it = floor(*it / binsize);
    }
    for (vector<double>::iterator it = unobserved_u.begin(); it != unobserved_u.end(); ++it) {
        *it = floor(*it / binsize);
    }
    
    
    // Build state structure from spike times in st:
    cout << "Building state histogram..." << endl;
    // Dump all the spikes into a vector and sort by spike time
    
    vector<Spike> all_spikes = sort_spikes(st, binsize);
    T = all_spikes.back().bin;
    
    vector<string> words (T,"");
    
    string silent_str (N,'0');
    
    // Add silent state with frequency of zero
    State this_state;
    string this_str = silent_str;
    this_state.freq = 0;
    this_state.P.assign(nbasins, 0);
    this_state.weight.assign(nbasins, 0);
    this_state.word.assign(N,0);
    train_states.insert(pair<string,State> (silent_str,this_state));
    //    test_states.insert(pair<string,State> (silent_str,this_state));
    
    int curr_bin = 0;
    for (vector<Spike>::iterator it=all_spikes.begin(); it!=all_spikes.end(); ++it) {
        
        int next_bin = it->bin;
        int next_cell = it->neuron_ind;
        
        bool bin_observed = true;
        for (int n=0; n<unobserved_l.size(); n++) {
            if (curr_bin >= unobserved_l[n] && curr_bin < unobserved_u[n]) {
                bin_observed = false;
                break;
            }
        }
        
        if (next_bin > curr_bin) {
            
            raster.push_back(this_str);
            if (bin_observed) {
                // Add new state; if it's already been discovered increment its frequency
                this_state.active_constraints = BasinT::get_active_constraints(this_state);
                
                pair<state_iter, bool> ins = train_states.insert(pair<string,State> (this_str,this_state));
                if (!ins.second) {
                    (((ins.first)->second).freq)++;
                }
                
            }
            // Update probabilities of time bins [curr_bin, next_bin)
            words[curr_bin] = this_str;
            for (int n=curr_bin+1; n<next_bin; n++) {
                words[n] = silent_str;
            }
            
            // All states between curr_bin and next_bin (exclusive) are silent; update frequency of silent state accordingly
            for (int t=curr_bin+1; t<next_bin; t++) {
                bool t_observed = true;
                for (int n=0; n<unobserved_l.size(); n++) {
                    if (t >= unobserved_l[n] && t < unobserved_u[n]) {
                        t_observed = false;
                        break;
                    }
                }
                if (t_observed) {
                    (train_states[silent_str].freq)++;
                }
                raster.push_back(silent_str);
            }
            
            
            // Reset state and jump to next bin
            this_str = silent_str;
            this_state.freq = 1;
            this_state.on_neurons.clear();
            this_state.P.assign(nbasins,0);
            this_state.weight.assign(nbasins,0);
            this_state.word.assign(N,0);
            
            curr_bin = next_bin;
        }
        
        // Add next_cell to this_state
        if (this_state.word[next_cell] == 0) {  // Don't want to count a cell twice in one bin
            this_str[next_cell] = '1';
            this_state.on_neurons.push_back(next_cell);
            this_state.word[next_cell] = 1;
        }
        
    }
    if (train_states[silent_str].freq == 0) {
        train_states.erase(silent_str);
    }
    
    // Now train_states contains all states found in the data together with their frequencies.
    
    
    int identifier = 0;
    for (state_iter it = train_states.begin(); it != train_states.end(); ++it) {
        nsamples += (it->second).freq;
        (it->second).identifier = identifier;
        identifier++;
    }
    
    state_list.assign(T, NULL);
    for (int t=0; t<T; t++) {
        bool t_observed = true;
        for (int n=0; n<unobserved_l.size(); n++) {
            if (t >= unobserved_l[n] && t < unobserved_u[n]) {
                t_observed = false;
                break;
            }
        }
        if (t_observed) {
            State& this_state = train_states.at(words[t]);
            state_list[t] = &this_state;
        } else {
            state_list[t] = 0;
        }
    }
    
}

template <class BasinT>
State EMBasins<BasinT>::build_state(vector<char> word) const {
    State this_state;
    this_state.P.assign(nbasins,0);
    this_state.weight.assign(nbasins,0);
    this_state.word = word;
    for (int i=0; i<N; i++) {
        if (word[i]) {
            this_state.on_neurons.push_back(i);
        }
    }
    this_state.active_constraints = BasinT::get_active_constraints(this_state);
    return this_state;
}


template <class BasinT>
EMBasins<BasinT>::~EMBasins() {
    delete rng;
}

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
template <class BasinT>
vector<double> EMBasins<BasinT>::train(int niter, bool ret_train_logli) {
    
    initParams();
    
    vector<double> train_logli (niter);
    for (int i=0; i<niter; i++) {
        //        cout << "Iteration " << i << endl;
        printf("Iteration %d\n", i);
        Estep();
        Mstep();
        
        //        logli[i] = update_P();
        //        test_logli[i] = update_P_test();
        train_logli[i] = logli(ret_train_logli);
    }
    return train_logli;
}

template <class BasinT>
void EMBasins<BasinT>::initParams() {
    w.assign(nbasins,1/(double)nbasins);
    basins.clear();
    // Initialize each basin model
    for (int i=0; i<nbasins; i++) {
        basins.push_back(BasinT(N,i,rng,0.45));
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
    
    for (state_iter it = train_states.begin(); it!=train_states.end(); ++it) {
        for (int j=0; j<nbasins; j++) {
            basins[j].increment_stats(it->second);
        }
    }
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
        w[i] = basins[i].get_norm() / nsamples;
    }
    
    return;
}


template <class BasinT>
void EMBasins<BasinT>::update_P() {
    for (state_iter it=train_states.begin(); it != train_states.end(); ++it) {
        State& this_state = it->second;
        double Z = 0;
        for (int i=0; i<nbasins; i++) {
            this_state.P[i] = basins[i].P_state(this_state);
            Z += w[i] * this_state.P[i];
        }
        for (int i=0; i<nbasins; i++) {
            this_state.weight[i] = this_state.freq * w[i] * this_state.P[i] / Z;
        }
        this_state.pred_prob = Z;
    }
    return;
}


template<class BasinT>
State EMBasins<BasinT>::state_obs(int obs, int t) const {
    // obs = 0:     test only
    // obs = 1:     train only
    // obs = 2:     all
    
    if (state_list[t] && (obs==1 || obs==2)) {
        return *state_list[t];
    } else if (!state_list[t] && (obs==0 || obs==2)) {
        
        State this_state;
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
    return State();
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


template <class BasinT>
double EMBasins<BasinT>::logli(bool obs) const {
    
    double logli = 0;
    int nsamp = 0;
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
    
    return logli;
}


template <class BasinT>
vector<unsigned long> EMBasins<BasinT>::state_hist() const {
    vector<unsigned long> hist (train_states.size(), 0);
    int pos = 0;
    for (const_state_iter it=train_states.begin(); it != train_states.end(); ++it) {
        const State& this_state = it->second;
        hist[pos] = this_state.freq;
        pos++;
    }
    return hist;
}

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
template <class BasinT>
vector<double> EMBasins<BasinT>::all_prob() const {
    vector<double> prob (train_states.size(), 0);
    int pos = 0;
    for (const_state_iter it=train_states.begin(); it != train_states.end(); ++it) {
        const State& this_state = it->second;
        prob[pos] = this_state.pred_prob;
        pos++;
    }
    return prob;
    
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

template <class BasinT>
vector<double> EMBasins<BasinT>::P() const {
    vector<double> P (train_states.size() * nbasins, 0);
    unsigned long pos = 0;
    for (const_state_iter it=train_states.begin(); it != train_states.end(); ++it) {
        const State& this_state = it->second;
        
        for (int i=0; i<nbasins; i++) {
            P[pos*nbasins + i] = (w[i] * this_state.P[i]) / this_state.pred_prob;
        }
        pos++;
    }
    return P;
}

template <class BasinT>
vector<paramsStruct> EMBasins<BasinT>::basin_params() {
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

template <class BasinT>
vector<Spike> EMBasins<BasinT>::sort_spikes(const vector<vector<double> >& st, double binsize) const {
    vector<Spike> all_spikes;
    for (int i=0; i<N; i++) {
        for (vector<double>::const_iterator it = st[i].begin(); it != st[i].end(); ++it) {
            Spike this_spike;
            double t = (*it);
            this_spike.bin = floor(t/binsize);
            this_spike.neuron_ind = i;
            all_spikes.push_back(this_spike);
        }
    }
    sort(all_spikes.begin(), all_spikes.end(), SpikeComparison());
    return all_spikes;
}




template <class BasinT>
vector<char> EMBasins<BasinT>::word_list() {
    vector<char> out (train_states.size() * N);
    vector<char>::iterator out_it = out.begin();
    for (state_iter it = train_states.begin(); it != train_states.end(); ++it) {
        vector<char> word = (it->second).word;
        for (vector<char>::iterator w_it = word.begin(); w_it!=word.end(); ++w_it) {
            *out_it++ = *w_it;
        }
    }
    return out;
}




