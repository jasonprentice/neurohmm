//
//  Autocorr.cpp
//  
//
//  Created by Jason Prentice on 6/11/14.
//
//

#include "Autocorr.h"


using namespace std;

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

