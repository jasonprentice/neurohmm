//
//  EMBasins.cpp
//  
//
//  Created by Jason Prentice on 11/13/12.
//
//

#include "EMBasins.h"
#include "BasinModel.h"
#include "TreeBasin.h"

#include "matrix.h"
#include "mex.h"

#include <queue>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <exception>
#include <cstdlib>
#include <ctime>

#include <thread>


//#include <gperftools/profiler.h>

// Selects which basin model to use
typedef IndependentBasin BasinType;

template <typename T>
void writeOutputMatrix(int pos, vector<T> value, int N, int M, mxArray**& plhs) {
    mxArray* out_matrix = mxCreateDoubleMatrix(N,M,mxREAL);
    double* pr = mxGetPr(out_matrix);
    for (typename vector<T>::iterator it=value.begin(); it!=value.end(); ++it) {
        *pr++ = (double) *it;
    }
    plhs[pos] = out_matrix;
    return;
}

void writeOutputStruct(int pos, vector<paramsStruct>& value, mxArray**& plhs) {
//    for (int n=0; n<value.size(); n++) {
    mxArray* out_struct = mxCreateStructMatrix(value.size(),1,value[0].get_nfields(), value[0].fieldNamesArray());
    int n = 0;
    for (vector<paramsStruct>::iterator it = value.begin(); it != value.end(); ++it) {
        for (int i=0; i < it->get_nfields(); i++) {
            vector<double>* data = it->getFieldData(i);
            
            int N = it->getFieldN(i);
            int M = it->getFieldM(i);
            mxArray* field_data = mxCreateDoubleMatrix(N, M, mxREAL);
            double* pr = mxGetPr(field_data);
            double* iter_pr = pr;
            for (vector<double>::iterator data_it=data->begin(); data_it != data->end(); ++data_it) {
                *iter_pr++ = *data_it;
            }            
            mxSetFieldByNumber(out_struct, n, i, field_data);
        }
        n++;
    }
    plhs[pos] = out_struct;
    return;
}


vector<double> mpow(const vector<double>& matrix, int n, int k) {

    if (k==1) {
        return matrix;
    } else if (k==0) {
        vector<double> ident (n*n,0);
        for (int i=0; i<n; i++) {
            ident[n*i + i] = 1;
        }
        return ident;
    } else {
        vector<double> ml = mpow(matrix, n, floor(k/2));
        vector<double> mr = mpow(matrix, n, k - floor(k/2));
        vector<double> mout (n*n);
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                for (int k=0; k<n; k++) {
                    mout[n*i+j] += ml[n*k + j] * mr[n*i + k];
                }
            }
        }
        return mout;        
    }

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
 // [freq,w,m,P,logli,prob] = EMBasins(st, unobserved_edges, binsize, nbasins, niter)

    cout << "Reading inputs..." << endl;
    int N = mxGetNumberOfElements(prhs[0]);
    vector<vector<double> > st (N);    
    for (int i=0; i<N; i++) {
        mxArray* elem = mxGetCell(prhs[0], i);
        double* elem_pr = mxGetPr(elem);
        int nspikes = mxGetNumberOfElements(elem);
        for (int n=0; n<nspikes; n++) {
            st[i].push_back(elem_pr[n]);
        }
    }
    
    int n_unobserved_blocks = mxGetM(prhs[1]);
    vector<double> unobserved_edges_low (n_unobserved_blocks);
    vector<double> unobserved_edges_high (n_unobserved_blocks);
    if (n_unobserved_blocks > 0 && mxGetN(prhs[1]) != 2) {
        cerr << "Unobserved edges must be empty or have two columns." << endl;
    } else {
        double* unobserved_edges_pr = mxGetPr(prhs[1]);
        for (int n=0; n<n_unobserved_blocks; n++) {
            unobserved_edges_low[n] = unobserved_edges_pr[n];
            unobserved_edges_high[n] = unobserved_edges_pr[n_unobserved_blocks + n];
        }
    }
    
    double binsize = *mxGetPr(prhs[2]);
    int nbasins = (int) *mxGetPr(prhs[3]);
    int niter = (int) *mxGetPr(prhs[4]);

    // ProfilerStart("./EMBasins.prof");
   /*
    // Autocorrelation model
    Autocorr<BasinType> basin_obj(st, binsize, nbasins);
    vector<double> logli = basin_obj.train(niter);
    cout << "Viterbi..." << endl;
    vector<int> alpha = basin_obj.viterbi();
//    cout << "P...." << endl;
//    vector<double> P = basin_obj.get_P();
//    cout << "Pred prob..." << endl;
//    pair<vector<double>, vector<double> > tmp = basin_obj.pred_prob();
//    vector<double> pred_prob = tmp.first;
//    vector<double> hist = tmp.second;
    //    vector<unsigned long> hist = basin_obj.state_hist();
//    int T = floor(P.size() / nbasins);
    int T = alpha.size();
    
    cout << "Params..." << endl;
    vector<paramsStruct> params = basin_obj.basin_params();
    
    writeOutputMatrix(0, logli, niter, 1, plhs);
//    writeOutputMatrix(1, basin_obj.get_trans(), nbasins, nbasins, plhs);
    //    writeOutputMatrix(2, P, nbasins, T, plhs);
//    writeOutputMatrix(2, basin_obj.emiss_prob(), nbasins, T, plhs);
    //    cout << "Microstates..." << endl;
    //    writeOutputMatrix(2, basin_obj.state_v_time(), 1, T, plhs);
    writeOutputMatrix(1, alpha, T, 1, plhs);
//    writeOutputMatrix(4, pred_prob, 1, pred_prob.size(), plhs);
//    writeOutputMatrix(5, hist, 1, hist.size(), plhs);
    writeOutputStruct(2, params, plhs);
    writeOutputMatrix(3, basin_obj.get_forward(), nbasins,T,plhs);
    writeOutputMatrix(4, basin_obj.get_backward(), nbasins,T,plhs);
    writeOutputMatrix(5, basin_obj.get_basin_trans(), 4*N,nbasins,plhs);
    writeOutputMatrix(6, basin_obj.w, nbasins,1,plhs);
    writeOutputMatrix(7, basin_obj.P_indep(), nbasins, T, plhs);
//    cout << "Samples..." << endl;
//    writeOutputMatrix(3, basin_obj.sample(100000), N,100000, plhs);
    */

   
    
    // Hidden Markov model
    cout << "Constructing HMM object" << endl;
    HMM<BasinType> basin_obj(st, unobserved_edges_low, unobserved_edges_high, binsize, nbasins);
    bool ret_train_logli = true;
    cout << "Training." << endl;
    vector<double> logli = basin_obj.train(niter, ret_train_logli);

    cout << "Viterbi..." << endl;
    vector<int> alpha = basin_obj.viterbi(2);
    cout << "P...." << endl;
    vector<double> P = basin_obj.get_P();
    cout << "Pred prob..." << endl;
    pair<vector<double>, vector<double> > tmp = basin_obj.pred_prob();
    vector<double> pred_prob = tmp.first;
    vector<double> hist = tmp.second;
//    vector<unsigned long> hist = basin_obj.state_hist();
    int T = floor(P.size() / nbasins);

    cout << "Params..." << endl;
    vector<paramsStruct> params = basin_obj.basin_params();
    
    writeOutputMatrix(0, logli, niter, 1, plhs);
    writeOutputMatrix(1, basin_obj.get_trans(), nbasins, nbasins, plhs);
    writeOutputMatrix(2, P, nbasins, T, plhs);
//    writeOutputMatrix(2, basin_obj.emiss_prob(), nbasins, T, plhs);
//    cout << "Microstates..." << endl;
//    writeOutputMatrix(2, basin_obj.state_v_time(), 1, T, plhs);
    writeOutputMatrix(3, alpha, T, 1, plhs);
    writeOutputMatrix(4, pred_prob, 1, pred_prob.size(), plhs);
    writeOutputMatrix(5, hist, 1, hist.size(), plhs);
    writeOutputStruct(6, params, plhs);
    //cout << "Samples..." << endl;
    writeOutputMatrix(7, basin_obj.sample(100000), N,100000, plhs);
//    writeOutputMatrix(7, basin_obj.word_list(), N, hist.size(), plhs);
//    writeOutputMatrix(6, basin_obj.stationary_prob(), 1,nbasins, plhs);
    
    
    /*
    // Mixture model
    cout << "Initializing EM..." << endl;
    EMBasins<BasinType> basin_obj(st, unobserved_edges_low, unobserved_edges_high, binsize, nbasins);
        
    cout << "Training model..." << endl;
    bool ret_train_logli = false;
    vector<double> logli = basin_obj.train(niter, ret_train_logli);
    //vector<double> test_logli = basin_obj.test_logli;
    
//    cout << "Testing..." << endl;
//    vector<double> P_test = basin_obj.test(st_test,binsize);
    
    vector<paramsStruct> params = basin_obj.basin_params();
    int nstates = basin_obj.nstates();
    cout << nstates << " states." << endl;
    
//    cout << "Getting samples..." << endl;
    int nsamples = 100000;
    vector<char> samples = basin_obj.sample(nsamples);


    cout << "Writing outputs..." << endl;    
    
    writeOutputMatrix(0, basin_obj.w, nbasins,1, plhs);    
    writeOutputStruct(1, params, plhs);
//    writeOutputMatrix(2, basin_obj.word_list(), N, nstates, plhs);
    writeOutputMatrix(2, samples, N, nsamples, plhs);
    writeOutputMatrix(3, basin_obj.state_hist(), nstates, 1, plhs);
    writeOutputMatrix(4, basin_obj.P(), nbasins, nstates, plhs);
    writeOutputMatrix(5, basin_obj.all_prob(), nstates, 1, plhs);
    writeOutputMatrix(6, logli, niter, 1, plhs);
//    writeOutputMatrix(6, P_test, nbasins, P_test.size()/nbasins, plhs);
    */
    
    
    /*
    // k-fold cross-validation
    int kfolds = 10;
    vector<double> logli = basin_obj.crossval(niter, kfolds);
    
    writeOutputMatrix(0, logli, niter, kfolds, plhs);
    */
    
  //  ProfilerStop();
    return;
}


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

template <class BasinT>
EMBasins<BasinT>::EMBasins(int N, int nbasins) : N(N), nbasins(nbasins), w(nbasins) {
    rng = new RNG();
    srand(time(NULL));
    //srand(0);
}


template <class BasinT>
EMBasins<BasinT>::EMBasins(vector<vector<double> >& st, vector<double> unobserved_l, vector<double> unobserved_u,double binsize, int nbasins) : N( st.size() ), nbasins(nbasins), nsamples(0), w(nbasins) {
    
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

};

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
        cout << "Iteration " << i << endl;

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



// ************* HMM **********************

template <class BasinT>
HMM<BasinT>::HMM(vector<vector<double> >& st, vector<double> unobserved_l, vector<double> unobserved_u, double binsize, int nbasins) : EMBasins<BasinT> (st, unobserved_l, unobserved_u, binsize, nbasins),
                            forward(this->T * nbasins,0),
                            backward(this->T * nbasins,0),
                            trans(this->T*nbasins,0),
                            w0 (nbasins),
                            barrier (4) {}

template <class BasinT>
vector<int> HMM<BasinT>::state_v_time() const {
    vector<int> states (this->T,-1);
    for (int t=0; t<this->T; t++) {
        if (this->state_list[t]) {
            states[t] = this->state_list[t]->identifier;
        }
    }
    return states;
}
                      
template <class BasinT>
vector<double> HMM<BasinT>::emiss_prob() const {
    vector<double> emiss (this->nbasins * this->T);
    for (int t=0; t<this->T; t++) {
//        State& this_state = this->train_states[words[t]];
//        State& this_state = *(state_list[t]);
        for (int i=0; i<this->nbasins; i++) {
//            emiss[t*this->nbasins + i] = this->basins[i].P_state(this_state);
            if (this->state_list[t]) {
                emiss[t*this->nbasins + i] = (this->state_list[t]->P)[i];
            } else {
                emiss[t*this->nbasins + i] = 1;
            }
        }
    }
    return emiss;
}

template <class BasinT>
vector<double> HMM<BasinT>::get_forward() const {
    return forward;
}

template <class BasinT>
vector<double> HMM<BasinT>::get_backward() const {
    return backward;
}

template <class BasinT>
vector<double> HMM<BasinT>::train(int niter, bool ret_train_logli) {
    this->initParams();

    vector<double> train_logli;
    
   
    
    // Create 3 threads: one to do update_forward, then update_trans; one to do update_backward; and one to do logli. The main thread will update emission probabilities and basin parameters
    
    // 1. Main thread updates emission probabilities.
    // 2. Worker threads update forward, backward, and logli in parallel (no race condition as long as emission probabilities and trans are not written to)
    // 3. Main thread does M-step on basin params while worker thread updates trans
    
    // Still need some mechanism to ensure each worker thread runs once per main loop iteration
    // Could create/destroy threads inside loop, but probably expensive?
    
    forward_updated = backward_updated = emiss_updated = false;
    
    
    thread forward_trans_thread(HMM<BasinT>::forward_trans_thread_fun, this);
    thread backward_thread(HMM<BasinT>::backward_thread_fun, this);
    thread logli_thread(HMM<BasinT>::logli_thread_fun, this, &train_logli, ret_train_logli);

    cout << "Beginning train loop." << endl;
    for (int i=0; i<niter; i++) {
        //cout << "Iteration " << i << endl;
        mexPrintf("Iteration %d\n", i);
        mexEvalString("drawnow");
        
      
        update_emiss();

        // Signal all worker threads
        unique_lock<mutex> lck (this->emiss_flag_mtx);
        emiss_updated = true;
        cv_emiss.notify_all();
        lck.unlock();
        
        // Wait for forward thread
        unique_lock<mutex> lck_fwd(this->fwd_flag_mtx);
        while (!this->forward_updated) {
            this->cv_fwd.wait(lck_fwd);
        }
        lck_fwd.unlock();
        
        // Wait for backward_thread
        unique_lock<mutex> lck_bkwd(this->bkwd_flag_mtx);
        while (!this->backward_updated) {
            this->cv_bkwd.wait(lck_bkwd);
        }
        lck_bkwd.unlock();
        
        
        update_P();
        
//        this->Estep();
        Mstep();

        if (this->barrier.wait()) {
            forward_updated = backward_updated = emiss_updated = false;
        };
        
//        forward_trans_thread.join();
//        backward_thread.join();
//        logli_thread.join();
        
//        train_logli[i] = logli(ret_train_logli);
    }
    return train_logli;
}

template <class BasinT>
void HMM<BasinT>::forward_trans_thread_fun(HMM<BasinT>* instance) {
    while (1) {
        
        // Wait for signal from main thread.
        unique_lock<mutex> lck (instance->emiss_flag_mtx);
        while (!instance->emiss_updated) {
            instance->cv_emiss.wait(lck);
        }
        lck.unlock();
        
        instance->update_forward();

        // Signal main thread
        unique_lock<mutex> lck_fwd (instance->fwd_flag_mtx);
        instance->forward_updated = true;
        instance->cv_fwd.notify_all();
        lck_fwd.unlock();

        
        // Wait for backward_thread
        unique_lock<mutex> lck_bkwd(instance->bkwd_flag_mtx);
        while (!instance->backward_updated) {
            instance->cv_bkwd.wait(lck_bkwd);
        }
        lck_bkwd.unlock();
        
        instance->update_trans();
        
        if (instance->barrier.wait()) {
            instance->forward_updated = instance->backward_updated = instance->emiss_updated = false;
        };

//    return;
    }
}

template <class BasinT>
void HMM<BasinT>::logli_thread_fun(HMM<BasinT>* instance, vector<double>* train_logli, bool ret_train_logli) {
    while (1) {
        // Wait for signal from main thread.
        unique_lock<mutex> lck (instance->emiss_flag_mtx);
        while (!instance->emiss_updated) {
            instance->cv_emiss.wait(lck);
        }
        lck.unlock();
        
        train_logli->push_back( instance->logli(ret_train_logli) );
        
        if (instance->barrier.wait()) {
            instance->forward_updated = instance->backward_updated = instance->emiss_updated = false;
        };

//    return;
    }
}

template <class BasinT>
void HMM<BasinT>::backward_thread_fun(HMM<BasinT>* instance) {
    while (1) {
        // Wait for signal from main thread.
        unique_lock<mutex> lck (instance->emiss_flag_mtx);
        while (!instance->emiss_updated) {
            instance->cv_emiss.wait(lck);
        }
        lck.unlock();
        
        instance->update_backward();

        // Signal forward_trans_thread and main thread
        unique_lock<mutex> bkwd_lck (instance->bkwd_flag_mtx);
        instance->backward_updated = true;
        instance->cv_bkwd.notify_all();
        bkwd_lck.unlock();
        
        if (instance->barrier.wait()) {
            instance->forward_updated = instance->backward_updated = instance->emiss_updated = false;
        };

   // return;
    }
}

template <class BasinT>
void HMM<BasinT>::initParams() {
    this->basins.clear();
    // Initialize each basin model
    for (int i=0; i<this->nbasins; i++) {
        this->basins.push_back(BasinT(this->N,i,this->rng,0.45));
    }

    
    w0.assign(this->nbasins, 1/(double)this->nbasins);
    trans.assign(this->nbasins*this->nbasins, 1/(double)this->nbasins);
    return;
}

template <class BasinT>
void HMM<BasinT>::Estep() {

    for (state_iter it=this->train_states.begin(); it != this->train_states.end(); ++it) {
        State& this_state = it->second;
        for (int i=0; i<this->nbasins; i++) {
            this_state.P[i] = this->basins[i].P_state(this_state);
        }
    }

    update_forward();
    update_backward();
    update_P();
    
    return;
}

template <class BasinT>
void HMM<BasinT>::Mstep() {

    for (int j=0; j<this->nbasins; j++) {
        this->basins[j].reset_stats();
    }
    
    for (state_iter it = this->train_states.begin(); it!=this->train_states.end(); ++it) {
        for (int j=0; j<this->nbasins; j++) {
            this->basins[j].increment_stats(it->second);
        }
    }
    for (int j=0; j<this->nbasins; j++) {
        this->basins[j].normalize_stats();
    }

    //        double alpha = (i<niter/2) ? 1 - (double)i/(niter/2) : 0;
    //        double alpha = 0.00002;
    double alpha = 0.002;
    //        if (i >= niter/2) {
    //            alpha = 0.002 + (1-0.002)*exp(-(double) (i-niter/2) * (10.0/(((double)(niter/2)-1))));
    //        }
    //        cout << alpha << endl;
    for (int j=0; j<this->nbasins; j++) {
        this->basins[j].doMLE(alpha);
    }
    
//    update_trans();
    return;
}

template <class BasinT>
void HMM<BasinT>::update_emiss() {
    for (state_iter it=this->train_states.begin(); it != this->train_states.end(); ++it) {
        State& this_state = it->second;
        for (int i=0; i<this->nbasins; i++) {
            this_state.P[i] = this->basins[i].P_state(this_state);
        }
    }

    return;
}

template <class BasinT>
void HMM<BasinT>::update_forward() {
    double norm = 0;
    for (int n=0; n<this->nbasins; n++) {
        if (this->state_list[this->T-1]) {
            forward[(this->T-1)*this->nbasins+n] = (this->state_list[this->T-1]->P)[n];
        } else {
            forward[(this->T-1)*this->nbasins+n] = 1;
        }
        norm += forward[(this->T-1)*this->nbasins+n];
    }
    
    for (int n=0; n<this->nbasins; n++) {
        forward[(this->T-1)*this->nbasins+n] /= norm;
    }
    
    
    for (int t=this->T-2; t>=0; t--) {
        double norm = 0;
        for (int n=0; n<this->nbasins; n++) {
            if (this->state_list[t]) {
                forward[t*this->nbasins + n] = (this->state_list[t]->P)[n];
            } else {
                forward[t*this->nbasins + n] = 1;
            }
            double tmp = 0;
            for (int m=0; m<this->nbasins; m++) {
                tmp += trans[n*this->nbasins + m] * forward[(t+1)*this->nbasins + m];
            }
            forward[t*this->nbasins + n] *= tmp;            
            norm += forward[t*this->nbasins + n];
        }
        
        for (int n=0; n<this->nbasins; n++) {
            forward[t*this->nbasins + n] /= norm;
        }
         
    }
    
    return;
}


template <class BasinT>
void HMM<BasinT>::update_backward() {    
    for (int n=0; n<this->nbasins; n++) {
        backward[n] = w0[n];
    }
 
    for (int t=1; t<this->T; t++) {
        double norm = 0;
        for (int n=0; n<this->nbasins; n++) {
            backward[this->nbasins*t + n] = 0;
            for (int m=0; m<this->nbasins; m++) {
                if (this->state_list[t-1]) {
                    backward[this->nbasins*t + n] += (this->state_list[t-1]->P)[m] * trans[m*this->nbasins+n] * backward[(t-1)*this->nbasins + m];
                } else {
                    backward[this->nbasins*t + n] += trans[m*this->nbasins+n] * backward[(t-1)*this->nbasins + m];
                }
            }
            norm += backward[this->nbasins*t + n];
        }
        for (int n=0; n<this->nbasins; n++) {
            backward[this->nbasins*t+n] /= norm;
        }

    }
    
    return;
}

template <class BasinT>
void HMM<BasinT>::update_trans() {
    // Update w0
    double norm=0;
    for (int n=0; n<this->nbasins; n++) {
        w0[n] *= forward[n];
        norm += w0[n];
    }

    for (int n=0; n<this->nbasins; n++) {
        w0[n] /= norm;
    }
    
    // Update trans
    vector<double> num (this->nbasins*this->nbasins,0);
    vector<double> prob (this->nbasins*this->nbasins);
    for (int t=1; t<this->T; t++) {
//        State& this_state = this->train_states.at(words[t-1]);

        double norm = 0;
        for (int n=0; n<this->nbasins; n++) {
            double tmp;
            if (this->state_list[t-1]) {
                tmp = (this->state_list[t-1]->P)[n] * backward[(t-1)*this->nbasins+n];
            } else {
                tmp = backward[(t-1)*this->nbasins+n];
            }
            for (int m=0; m<this->nbasins; m++) {
                prob[n*this->nbasins+m] = tmp * trans[n*this->nbasins + m] * forward[t*this->nbasins + m];
                norm += prob[n*this->nbasins + m];
            }
        }
        for (int n=0; n<this->nbasins*this->nbasins; n++) {
            prob[n] /= norm;
        }
        for (int n=0; n<this->nbasins; n++) {
            for (int m=0; m<this->nbasins; m++) {
                double delta_num = prob[n*this->nbasins+m] - num[n*this->nbasins+m];
                num[n*this->nbasins+m] += delta_num / t;
            }
        }
    }
    for (int n=0; n<this->nbasins; n++) {
        double norm = 0;
        for (int m=0; m<this->nbasins; m++) {
            norm += num[n*this->nbasins+m];
        }
        for (int m=0; m<this->nbasins; m++) {
            trans[n*this->nbasins+m] = num[n*this->nbasins+m]/norm;
        }
    }
    return;
}

template <class BasinT>
void HMM<BasinT>::update_P() {

    for (state_iter it=this->train_states.begin(); it != this->train_states.end(); ++it) {
        State& this_state = it->second;
        this_state.weight.assign(this->nbasins,0);
    }
    
    vector<double> denom  (this->nbasins,0);
    int nsamp = 0;
    vector<double>  this_P (this->nbasins,0);
    for (int t=0; t<this->T; t++) {
//        State& this_state = this->train_states.at(words[t]);
        if (this->state_list[t]) {
            State& this_state = *(this->state_list[t]);
            
            double norm = 0;
            for (int i=0; i<this->nbasins; i++) {
                this_P[i] = forward[t*this->nbasins + i] * backward[t*this->nbasins + i];
                norm += this_P[i];
            }
            for (int i=0; i<this->nbasins; i++) {
                this_P[i] /= norm;
                
                //double delta = this_P[i] - this_state.weight[i];
               // this_state.weight[i] += delta / (i+1);
                this_state.weight[i] += this_P[i];
                //denom[i] += this_P[i];
//                denom[i] += (delta_denom / ((t/tskip)+1));
                double delta_denom = this_P[i] - denom[i];
                denom[i] += (delta_denom / (nsamp+1));
            }
            nsamp++;
        }
    }
    
    for (state_iter it=this->train_states.begin(); it != this->train_states.end(); ++it) {
        State& this_state = it->second;
        for (int i=0; i<this->nbasins; i++) {
//            this_state.weight[i] /= (ceil(T/tskip)*denom[i]);
            this_state.weight[i] /= (nsamp*denom[i]);
//            this_state.weight[i] /= (denom[i]);
//            this_state.P[i] = this->basins[i].P_state(this_state);
        }
    }

    return;
}

template <class BasinT>
double HMM<BasinT>::logli(bool obs) const {
    
    vector<int> alpha = viterbi(obs);
    vector<double> emiss = this->emiss_obs(obs, 0);

    double logli = log(w0[alpha[0]] * emiss[alpha[0]]);
    double nsamp = this->state_list[0]==0 ? 0 : 1;
    for (int t=1; t<this->T; t++) {
        emiss = this->emiss_obs(obs, t);
        double delta = log(trans[alpha[t-1]*this->nbasins + alpha[t]]) + log(emiss[alpha[t]]) - logli;;
        if (this->state_list[t]) nsamp++;
        logli += delta / t;
    }
    return logli * (this->T/nsamp);
}


template <class BasinT>
vector<double> HMM<BasinT>::get_P() const {
    vector<double> P (this->T*this->nbasins);
    for (int t=0; t<this->T; t++) {
        double norm = 0;
        for (int n=0; n<this->nbasins; n++) {
            P[t*this->nbasins+n] = forward[t*this->nbasins+n]*backward[t*this->nbasins+n];
            norm += P[t*this->nbasins+n];
        }
        for (int n=0; n<this->nbasins; n++) {
            P[t*this->nbasins+n] /= norm;
        }
    }
    return P;
}

template <class BasinT>
vector<double> HMM<BasinT>::P_indep() const {
    vector<double> P (this->T*this->nbasins);
    for (int t=0; t<this->T; t++) {
        double norm = 0;        
        for (int n=0; n<this->nbasins; n++) {
            P[t*this->nbasins+n] = this->w[n] * this->state_list[t]->P[n];
            norm += P[t*this->nbasins+n];
        }
        for (int n=0; n<this->nbasins; n++) {
            P[t*this->nbasins+n] /= norm;
        }

    }
    return P;
}

template <class BasinT>
vector<double> HMM<BasinT>::get_trans() const{
    return trans;
}



template <class BasinT>
vector<int> HMM<BasinT>::viterbi(int obs) const {

    vector<int> alpha_max (this->T,0);
    vector<int> argmax (this->T*this->nbasins, 0);

    vector<double> max (this->nbasins, 1);
    vector<double> emiss (this->nbasins);
    for (int t = this->T-1; t>=0; t--) {
//        State& this_state = this->train_states.at(words[t]);
        emiss = this->emiss_obs(obs,t);
        double norm = 0;
        for (int n=0; n<this->nbasins; n++) {
            double this_max = 0;
            int this_arg = 0;
            for (int m=0; m<this->nbasins; m++) {
                double tmp = emiss[m] * trans[n*this->nbasins+m] * max[m];
//                if (state_list[t]) {
//                    tmp = (state_list[t]->P)[m] * trans[n*this->nbasins+m] * max[m];
//                } else {
//                    tmp = trans[n*this->nbasins+m] * max[m];
//                }
                
                if (tmp > this_max) {
                    this_max = tmp;
                    this_arg = m;
                }
            }
            max[n] = this_max;
            norm += this_max;
            argmax[t*this->nbasins + n] = this_arg;
        }
        for (int n=0; n<this->nbasins; n++) {
            max[n] /= norm;
        }
    }
//    State& this_state = this->train_states.at(words[0]);
    double this_max = 0;
    int this_arg = 0;
    emiss = this->emiss_obs(obs,0);
    for (int m=0; m<this->nbasins; m++) {
        double tmp = emiss[m] * w0[m] * max[m];
        if (tmp > this_max) {
            this_max = tmp;
            this_arg = m;
        }
    }
    alpha_max[0] = this_arg;
    for (int t=1; t<this->T; t++) {
        alpha_max[t] = argmax[t*this->nbasins + alpha_max[t-1]];
    }
    
    return alpha_max;
}

template <class BasinT>
vector<double> HMM<BasinT>::stationary_prob() const {
    // Stationary basin probability
    vector<double> trans_pow = mpow(trans, this->nbasins, 1000);
    vector<double> w (this->nbasins,0);
    for (int i=0; i<this->nbasins; i++) {
        for (int j=0; j<this->nbasins; j++) {
            w[i] += trans_pow[(this->nbasins)*j+i]*w0[j];
        }
    }
    return w;
}


template <class BasinT>
vector<char> HMM<BasinT>::sample(int nsamples) const {
    vector<char> samples (this->N*nsamples);
    int basin_ind = (this->rng)->discrete(w0);
    vector<char> this_sample = (this->basins)[basin_ind].sample();
    for (int n=0; n<this->N; n++) {
        samples[n] = this_sample[n];
    }
    for (int t=1; t<nsamples; t++) {
        vector<double> this_trans (this->nbasins);
        for (int k=0; k<this->nbasins; k++) {
            this_trans[k] = trans[basin_ind*this->nbasins + k];
        }
        basin_ind = (this->rng)->discrete(this_trans);
        this_sample = (this->basins)[basin_ind].sample();
        for (int n=0; n<this->N; n++) {
            samples[t*this->N + n] = this_sample[n];
        }
    }
    return samples;
}

template <class BasinT>
pair<vector<double>, vector<double> > HMM<BasinT>::pred_prob() const {
    
//    this->test_states.clear();
    map<string,State> test_states;
    for (int t=0; t<this->T; t++) {
        State this_state = this->state_obs(2,t);
        vector<char> this_word = this_state.word;
        string this_str (this_word.size(), '0');
        for (int i=0; i<this_word.size(); i++) {
            this_str[i] = this_word[i];
        }

        pair<map<string, State>::iterator, bool> ins = test_states.insert(pair<string,State> (this_str, this_state));
        State& inserted_state = (ins.first)->second;
        inserted_state.freq++;
        
        /*
        if (this->state_list[t]) {
            State this_state = *(this->state_list[t]);
            this_state.freq = 0;
            vector<char> this_word = this_state.word;
            string this_str (this_word.size(), '0');
            for (int i=0; i<this_word.size(); i++) {
                this_str[i] = this_word[i];
            }
            pair<map<string, State>::iterator, bool> ins = this->test_states.insert(pair<string,State> (this_str, this_state));
            State& inserted_state = (ins.first)->second;
            inserted_state.freq++;
        }
         */
    }
    
    vector<double> w = stationary_prob();
    vector<double> prob (test_states.size(), 0);
    vector<double> freq (test_states.size(), 0);
    int ix = 0;
    for (map<string, State>::iterator it = test_states.begin(); it != test_states.end(); ++it) {
        State& this_state = it->second;
        for (int i=0; i<this->nbasins; i++) {
            prob[ix] += w[i] * this_state.P[i];
            freq[ix] = (it->second).freq;
        }
        ix++;
    }

    return pair<vector<double>, vector<double> > (prob, freq);
    
    
}


template <class BasinT>
Autocorr<BasinT>::Autocorr(vector<vector<double> >& st, double binsize, int nbasins) : HMM<BasinT>(st,vector<double>(),vector<double>(),binsize,nbasins), basin_trans (nbasins * st.size()) {

    for (vector<double*>::iterator it = basin_trans.begin(); it != basin_trans.end(); ++it) {
        *it = new double[4];
    }
}


template <class BasinT>
Autocorr<BasinT>::~Autocorr() {
    for (vector<double*>::iterator it = basin_trans.begin(); it != basin_trans.end(); ++it) {
        delete[] *it;
    }
}


template <class BasinT>
void Autocorr<BasinT>::update_basin_trans_indep() {
    vector<vector<double> > basin_trans_num (this->nbasins * this->N, vector<double> (4,0));
    vector<vector<double> > basin_trans_den (this->nbasins * this->N, vector<double> (2,0));
    vector<double> denom  (this->nbasins,0);
    
    for (int t=1; t<this->T; t++) {
        vector<double> P_joint(this->nbasins * this->nbasins, 0);
        double norm = 0;
        for (int n=0; n<this->nbasins; n++) {
            for (int m=0; m<this->nbasins; m++) {
                P_joint[n*this->nbasins+m] = this->state_list[t]->P[m] * this->state_list[t-1]->P[n] * this->w[n] * this->w[m];
                norm += P_joint[n*this->nbasins + m];
            }
        }
        for (int n=0; n<this->nbasins*this->nbasins; n++) {
            P_joint[n] /= norm;
        }
        for (int i=0; i<this->nbasins; i++) {
            for (int n=0; n<this->N; n++) {
                char sp_prev = (this->state_list[t-1]->word)[n];
                char sp_this = (this->state_list[t]->word)[n];

//                for (int a=0; a<4; a++) {
//                    double delta_num = -basin_trans_num[this->N*i+n][a];
//                    delta_num += (a == (sp_this + 2*sp_prev)) ? P_joint[i*this->nbasins+i] : 0;
//                    basin_trans_num[this->N*i + n][a] += delta_num / (t+1);
//                }
//
//                for (int a=0; a<2; a++) {
//                    double delta_den = -basin_trans_den[this->N*i+n][a];
//                    delta_den += (a == sp_prev) ? P_joint[i*this->nbasins+i] : 0;
//                    basin_trans_den[this->N*i + n][a] += delta_den / (t+1);
//                }
                
                basin_trans_num[this->N*i + n][sp_this + 2*sp_prev] += P_joint[i*this->nbasins+i];
                basin_trans_den[this->N*i + n][sp_prev] += P_joint[i*this->nbasins+i];
            }
        }
        
    }
    
    for (int i=0; i<this->nbasins * this->N; i++) {
        if (basin_trans_den[i][0] > 0) {
            basin_trans[i][0] = basin_trans_num[i][0] / basin_trans_den[i][0];
            // basin_trans[i][1] = 1 - basin_trans[i][0];
            basin_trans[i][1] = basin_trans_num[i][1] / basin_trans_den[i][0];
        } else {
            basin_trans[i][0] = 0.5;
            basin_trans[i][1] = 0.5;
        }
        if (basin_trans_den[i][1] > 0) {
            basin_trans[i][2] = basin_trans_num[i][2] / basin_trans_den[i][1];
            basin_trans[i][3] = basin_trans_num[i][3] / basin_trans_den[i][1];
        } else {
            basin_trans[i][2] = 0.5;
            basin_trans[i][3] = 0.5;
   
        }
        // basin_trans[i][3] = 1 - basin_trans[i][2];
    }

    
    return;
}



template <class BasinT>
void Autocorr<BasinT>::update_P() {
    
    for (state_iter it=this->train_states.begin(); it != this->train_states.end(); ++it) {
        State& this_state = it->second;
        this_state.weight.assign(this->nbasins,0);
    }
    
    vector<vector<double> > basin_trans_num (this->nbasins * this->N, vector<double> (4,0));
    vector<vector<double> > basin_trans_den (this->nbasins * this->N, vector<double> (2,0));
    vector<double> denom  (this->nbasins,0);
    
    for (int t=1; t<this->T; t++) {
        //        State& this_state = this->train_states.at(words[t]);
        State& this_state = *(this->state_list[t]);
        
        vector<double> this_P  (this->nbasins,0);
        vector<double> this_trans = trans_at_t(t);

        vector<double> P_joint(this->nbasins * this->nbasins, 0);
        double norm = 0;
        for (int n=0; n<this->nbasins; n++) {
            for (int m=0; m<this->nbasins; m++) {
                P_joint[n*this->nbasins+m] = this->backward[t*this->nbasins + m] * this_trans[n*this->nbasins + m] * this->forward[(t-1)*this->nbasins + n];
                norm += P_joint[n*this->nbasins + m];
            }
        }
        for (int n=0; n<this->nbasins*this->nbasins; n++) {
            P_joint[n] /= norm;
        }
        

        for (int i=0; i<this->nbasins; i++) {
            for (int j=0; j<this->nbasins; j++) {
                if (j != i) {
                    this_P[i] += P_joint[j*this->nbasins + i];
                }
            }
        }
        
        for (int i=0; i<this->nbasins; i++) {
            
            double delta_denom = this_P[i] - denom[i];

            this_state.weight[i] += this_P[i];

            denom[i] += (delta_denom / (t+1));
            
            for (int n=0; n<this->N; n++) {
                char sp_prev = (this->state_list[t-1]->word)[n];
                char sp_this = this_state.word[n];
                
                for (int a=0; a<4; a++) {
                    double delta_num = -basin_trans_num[this->N*i+n][a];
                    delta_num += (a == (sp_this + 2*sp_prev)) ? P_joint[i*this->nbasins+i] : 0;
                    basin_trans_num[this->N*i + n][a] += delta_num / (t+1);
                }
                
                for (int a=0; a<2; a++) {
                    double delta_den = -basin_trans_den[this->N*i+n][a];
                    delta_den += (a == sp_prev) ? P_joint[i*this->nbasins+i] : 0;
                    basin_trans_den[this->N*i + n][a] += delta_den / (t+1);
                }
//                basin_trans_num[this->N*i + n][sp_this + 2*sp_prev] += P_joint[i*this->nbasins+i];
//                basin_trans_den[this->N*i + n][sp_prev] += P_joint[i*this->nbasins+i];
            }
        }
    }
    
    for (int i=0; i<this->nbasins * this->N; i++) {
        if (basin_trans_den[i][0] > 0) {
            basin_trans[i][0] = basin_trans_num[i][0] / basin_trans_den[i][0];
            // basin_trans[i][1] = 1 - basin_trans[i][0];
            basin_trans[i][1] = basin_trans_num[i][1] / basin_trans_den[i][0];
        } else {
            basin_trans[i][0] = 0.5;
            basin_trans[i][1] = 0.5;
        }
        if (basin_trans_den[i][1] > 0) {
            basin_trans[i][2] = basin_trans_num[i][2] / basin_trans_den[i][1];
            basin_trans[i][3] = basin_trans_num[i][3] / basin_trans_den[i][1];
        } else {
            basin_trans[i][2] = 0.5;
            basin_trans[i][3] = 0.5;
            
        }
        // basin_trans[i][3] = 1 - basin_trans[i][2];
    }
    
    for (state_iter it=this->train_states.begin(); it != this->train_states.end(); ++it) {
        State& this_state = it->second;
        for (int i=0; i<this->nbasins; i++) {
            this_state.weight[i] /= (this->T * denom[i]);
            //            this_state.weight[i] /= (denom[i]);
            this_state.P[i] = this->basins[i].P_state(this_state);
        }
    }
    
    return;
}

template <class BasinT>
void Autocorr<BasinT>::update_w() {
    vector<double> next_w (this->nbasins,0);
    for (int t=0; t<this->T; t++) {
        for (int i=0; i<this->nbasins; i++) {
            double delta_w = this->forward[t*this->nbasins + i] * this->backward[t*this->nbasins + i] - next_w[i];
            next_w[i] += delta_w / (t+1);
        }
    }
    double norm = 0;
    for (int i=0; i<this->nbasins; i++) {
        norm += next_w[i];
    }
    for (int i=0; i<this->nbasins; i++) {
        next_w[i] /= norm;
    }
    this->w = next_w;
    return;
}


template <class BasinT>
vector<double> Autocorr<BasinT>::trans_at_t(int t) {
    vector<double> trans (this->nbasins * this->nbasins, 0);
    for (int a=0; a<this->nbasins; a++) {
        for (int b=0; b<this->nbasins; b++) {
            trans[this->nbasins*a + b] = this->w[b];
            if (a==b) {
                for (int n=0; n<this->N; n++) {
                    char sp_prev = this->state_list[t-1]->word[n];
                    char sp_this = this->state_list[t]->word[n];
                    trans[this->nbasins*a + b] *= basin_trans[this->N*a + n][sp_this + 2*sp_prev];
                }
            } else {
                trans[this->nbasins*a + b] *= this->state_list[t]->P[b];
            }
        }
    }
    return trans;
}

template <class BasinT>
void Autocorr<BasinT>::update_forward() {
    for (int n=0; n<this->nbasins; n++) {
        this->forward[n] = this->w[n] * (this->state_list[0]->P)[n];
    }
    for (int t=1; t<this->T; t++) {
        //        State& this_state = this->train_states.at(words[t-1]);
        vector<double> this_trans = trans_at_t(t);
        double norm = 0;
        for (int n=0; n<this->nbasins; n++) {
            this->forward[this->nbasins*t + n] = 0;
            for (int m=0; m<this->nbasins; m++) {
                this->forward[this->nbasins*t + n] += this_trans[m*this->nbasins+n] * this->forward[(t-1)*this->nbasins + m];
            }
            norm += this->forward[this->nbasins*t + n];
        }
        for (int n=0; n<this->nbasins; n++) {
            this->forward[this->nbasins*t+n] /= norm;
        }
        
    }
    return;

    
}

template <class BasinT>
void Autocorr<BasinT>::update_backward() {

    int tmax = this->T - 1;
    for (int n=0; n<this->nbasins; n++) {
        //        forward[(T-1)*this->nbasins+n] = final_state.P[n];
        this->backward[tmax*this->nbasins+n] = 1;
    }
    
    for (int t=tmax-1; t>=0; t--) {
        //        State& this_state = this->train_states[words[t]];
        double norm = 0;
        vector<double> this_trans = trans_at_t(t+1);
        
        for (int n=0; n<this->nbasins; n++) {
            //            forward[t*this->nbasins + n] = this_state.P[n];
            this->backward[t*this->nbasins + n] = 0;
            
            for (int m=0; m<this->nbasins; m++) {
                this->backward[t*this->nbasins + n] += this_trans[n*this->nbasins + m] * this->backward[(t+1)*this->nbasins + m];
            }
            norm += this->backward[t*this->nbasins + n];
        }
        
        for (int n=0; n<this->nbasins; n++) {
            this->backward[t*this->nbasins + n] /= norm;
        }
        
    }
    return;
}

template <class BasinT>
vector<double> Autocorr<BasinT>::train(int niter) {
    /*
    (this->w).assign(this->nbasins, 1/(double)this->nbasins);
    this->basins.clear();
    // Initialize each basin model
    for (int i=0; i<this->nbasins; i++) {
        this->basins.push_back(BasinT(this->N,i,this->rng));
    }
    
    // Initialize emission probabilities
    for (state_iter it=this->train_states.begin(); it != this->train_states.end(); ++it) {
        State& this_state = it->second;
        for (int i=0; i<this->nbasins; i++) {
            this_state.P[i] = this->basins[i].P_state(this_state);
        }
    }
    
    for (int i=0; i<this->nbasins * this->N; i++) {
//        basin_trans[i][0] = 0.1*((double) rand() / (double) RAND_MAX) + 0.45;
        basin_trans[i][0] = 0.5;
        basin_trans[i][1] = 1 - basin_trans[i][0];
//        basin_trans[i][2] = 0.1*((double) rand() / (double) RAND_MAX) + 0.45;
        basin_trans[i][2] = 0.5;
        basin_trans[i][3] = 1 - basin_trans[i][2];
    }
     */
    int uncorr_iter = 20;
    uncorr_iter = (uncorr_iter < niter) ? uncorr_iter : niter;
//    vector<double> train_logli_begin = this->EMBasins<BasinT>::train(uncorr_iter);
    vector<double> train_logli_begin = this->HMM<BasinT>::train(uncorr_iter);
//    for (int i=0; i<this->nbasins * this->N; i++) {
//     //        basin_trans[i][0] = 0.1*((double) rand() / (double) RAND_MAX) + 0.45;
//         basin_trans[i][0] = 0.5;
//         basin_trans[i][1] = 1 - basin_trans[i][0];
//         //        basin_trans[i][2] = 0.1*((double) rand() / (double) RAND_MAX) + 0.45;
//         basin_trans[i][2] = 0.5;
//         basin_trans[i][3] = 1 - basin_trans[i][2];
//    }

    //this->EMBasins<BasinT>::update_P();
    
    update_basin_trans_indep();
//    
    update_forward();
    update_backward();
    update_P();
    update_w();
    
    vector<double> train_logli (niter);
    for (int i=0; i<uncorr_iter; i++) {
        train_logli[i] = train_logli_begin[i];
    }
    for (int i=uncorr_iter; i<niter; i++) {
        cout << "Iteration " << i << endl;
        
        // E step
        
        for (int j=0; j<this->nbasins; j++) {
            this->basins[j].reset_stats();
        }
        
        for (state_iter it = this->train_states.begin(); it!=this->train_states.end(); ++it) {
            for (int j=0; j<this->nbasins; j++) {
                this->basins[j].increment_stats(it->second);
            }
        }
        for (int j=0; j<this->nbasins; j++) {
            this->basins[j].normalize_stats();
        }
        
        // M step
        
        //        double alpha = (i<niter/2) ? 1 - (double)i/(niter/2) : 0;
        double alpha = 0.002;
        //        if (i >= niter/2) {
        //            alpha = 0.002 + (1-0.002)*exp(-(double) (i-niter/2) * (10.0/(((double)(niter/2)-1))));
        //        }
        //        cout << alpha << endl;
        for (int j=0; j<this->nbasins; j++) {
            
            this->basins[j].doMLE(alpha);
        }

        cout << "Forward..." << endl;
        update_forward();
        cout << "Backward..." << endl;
        update_backward();
        cout << "P..." << endl;
        update_P();
        update_w();
        cout << "logli..." << endl;
        train_logli[i] = logli();
    }
    
    

    return train_logli;
    
}


template <class BasinT>
vector<int> Autocorr<BasinT>::viterbi() {
    
    vector<int> alpha_max (this->T,0);
    vector<int> argmax (this->T*this->nbasins, 0);
    
    vector<double> max (this->nbasins, 1);
    for (int t=(this->T-1); t>0; t--) {
        double norm = 0;
        vector<double> this_trans = trans_at_t(t);
        for (int n=0; n<this->nbasins; n++) {
            double this_max = 0;
            int this_arg = 0;
            for (int m=0; m<this->nbasins; m++) {
                double tmp = this_trans[n*this->nbasins+m] * max[m];
                
                if (tmp > this_max) {
                    this_max = tmp;
                    this_arg = m;
                }
            }
            max[n] = this_max;
            norm += this_max;
            argmax[t*this->nbasins + n] = this_arg;
        }
        for (int n=0; n<this->nbasins; n++) {
            max[n] /= norm;
        }
    }
    //    State& this_state = this->train_states.at(words[0]);
    double this_max = 0;
    int this_arg = 0;
    for (int m=0; m<this->nbasins; m++) {
        double tmp = this->w[m] * (this->state_list[0]->P)[m] * max[m];
        if (tmp > this_max) {
            this_max = tmp;
            this_arg = m;
        }
    }
    alpha_max[0] = this_arg;
    for (int t=1; t<this->T; t++) {
        alpha_max[t] = argmax[t*this->nbasins + alpha_max[t-1]];
    }
    
    return alpha_max;
}

template <class BasinT>
double Autocorr<BasinT>::logli() {
    
    vector<int> alpha = viterbi();
    //    State& init_state = this->train_states.at(words[0]);
    cout << "Viterbi done." << endl;
    
    double logli = log(this->w[alpha[0]] * (this->state_list[0]->P)[alpha[0]]);
    
    
    for (int t=1; t<this->T; t++) {
        //        State& this_state = this->train_states.at(words[t]);
        vector<double> this_trans = trans_at_t(t);
        double delta = log(this_trans[alpha[t-1]*this->nbasins + alpha[t]]) - logli;
        
        logli += delta / t;
    }
    return logli;
}

template <class BasinT>
vector<double> Autocorr<BasinT>::get_basin_trans() {
    
    vector<double> basin_trans_out (4 * this->nbasins * this->N);
    for (int n=0; n<this->N; n++) {
        for (int a=0; a<this->nbasins; a++) {
            for (int i=0; i<4; i++) {
                basin_trans_out[4*this->N*a + 4*n + i] = basin_trans[this->N*a + n][i];
            }
        }
    }
    return basin_trans_out;
}

