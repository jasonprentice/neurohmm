//
//  HMM.cpp
//  
//
//  Created by Jason Prentice on 6/11/14.
//
//

#include <thread>

#include "RNG.h"
#include "BasinModel.h"
#include "TreeBasin.h"
#include "IndependentBasin.h"


using namespace std;

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


// ************* HMM **********************

template <class BasinT>
HMM<BasinT>::HMM(vector<vector<double> >& st, vector<double> unobserved_l, vector<double> unobserved_u, double binsize, int nbasins)
    : EMBasins<BasinT> (st, unobserved_l, unobserved_u, binsize, nbasins),
      forward(this->T * nbasins,0),
      backward(this->T * nbasins,0),
      trans(this->T*nbasins,0),
      w0 (nbasins) {}



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
    
    
    
    
    cout << "Beginning train loop." << endl;
    
    for (int i=0; i<niter; i++) {
        //cout << "Iteration " << i << endl;
        printf("Iteration %d\n", i);
      //  mexEvalString("drawnow");
        
        forward_updated = backward_updated = emiss_updated = false;
        
        // Create 3 threads: one to do update_forward, then update_trans; one to do update_backward; and one to do logli. The main thread will update emission probabilities and basin parameters
        
        // 1. Main thread updates emission probabilities.
        // 2. Worker threads update forward, backward, and logli in parallel (no race condition as long as emission probabilities and trans are not written to)
        // 3. Main thread does M-step on basin params while worker thread updates trans
        
        
        thread forward_trans_thread (HMM<BasinT>::forward_trans_thread_fun, this);
        thread backward_thread (HMM<BasinT>::backward_thread_fun, this);
        thread logli_thread (HMM<BasinT>::logli_thread_fun, this, &train_logli, ret_train_logli);
        
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
        
        
        forward_trans_thread.join();
        backward_thread.join();
        logli_thread.join();
        
    }
    return train_logli;
}

template <class BasinT>
void HMM<BasinT>::forward_trans_thread_fun(HMM<BasinT>* instance) {
    
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
    
    return;
}

template <class BasinT>
void HMM<BasinT>::logli_thread_fun(HMM<BasinT>* instance, vector<double>* train_logli, bool ret_train_logli) {
    
    // Wait for signal from main thread.
    unique_lock<mutex> lck (instance->emiss_flag_mtx);
    while (!instance->emiss_updated) {
        instance->cv_emiss.wait(lck);
    }
    lck.unlock();
    
    train_logli->push_back( instance->logli(ret_train_logli) );
    
    return;
}

template <class BasinT>
void HMM<BasinT>::backward_thread_fun(HMM<BasinT>* instance) {
    
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
    
    return;
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
        State this_state = this->state_obs(0,t);
        if (this_state.identifier != -1) {
            
            vector<char> this_word = this_state.word;
            string this_str (this_word.size(), '0');
            for (int i=0; i<this_word.size(); i++) {
                this_str[i] = this_word[i];
            }
            
            pair<state_iter, bool> ins = test_states.insert(pair<string,State> (this_str, this_state));
            State& inserted_state = (ins.first)->second;
            inserted_state.freq++;
            
        }
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
    
    return pred_prob_helper(test_states);
    
    
}

template <class BasinT>
pair<vector<double>, vector<double> > HMM<BasinT>::pred_prob_helper(const map<string,State>& test_states) const {
    vector<double> w = stationary_prob();
    vector<double> prob (test_states.size(), 0);
    vector<double> freq (test_states.size(), 0);
    int ix = 0;
    for (auto it = test_states.begin(); it != test_states.end(); ++it) {
        State this_state = it->second;
        for (int i=0; i<this->nbasins; i++) {
            prob[ix] += w[i] * this_state.P[i];
            freq[ix] = (it->second).freq;
        }
        ix++;
    }
    
    return pair<vector<double>, vector<double> > (prob,freq);
}


template <class BasinT>
pair<vector<double>, vector<double> > HMM<BasinT>::sample_pred_prob(const vector<char>& sample) const {
    int nsamples = sample.size() / this->N;
    
    map<string,State> sample_states;
    vector<char> word (this->N);
    for (int t=0; t<nsamples; t++) {
        for (int i=0; i<this->N; i++) {
            word[i] = sample[this->N * t + i];
        }
        
        string this_str (this->N, '0');
        for (int i=0; i<this->N; i++) {
            this_str[i] = word[i];
        }
        
        State this_state = this->build_state(word);
        
        pair<state_iter, bool> ins = sample_states.insert(pair<string,State> (this_str, this_state));
        State& inserted_state = (ins.first)->second;
        inserted_state.freq++;
    }
    
    for (state_iter it = sample_states.begin(); it != sample_states.end(); ++it) {
        State& this_state = it->second;
        for (int i=0; i<this->nbasins; i++) {
            this_state.P[i] = this->basins[i].P_state(this_state);
        }
    }
    
    return pred_prob_helper(sample_states);
}

