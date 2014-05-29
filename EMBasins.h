//
//  EMBasins.h
//  
//
//  Created by Jason Prentice on 11/13/12.
//
//

#ifndef ____EMBasins__
#define ____EMBasins__

#include <gsl/gsl_rng.h>

#include <vector>
#include <string>
#include <map>

#include <mutex>
#include <condition_variable>

using namespace std;


// ************ RNG ***************
class RNG
{
public:
    RNG();
    ~RNG();
    int discrete(const vector<double>&) const;
    bool bernoulli(double) const;
    double uniform() const;
    vector<int> randperm(int) const;
private:
    gsl_rng* rng_pr;
};

struct State
{
    State() : freq(0), pred_prob(-1), identifier(-1) {};
    
    vector<int> active_constraints;
    vector<int> on_neurons;
    vector<double> P;
    vector<double> weight;
    double freq;
    double pred_prob;

    vector<char> word;
    
    int identifier;
};
// *********************************
typedef map<string,State>::iterator state_iter;
typedef map<string,State>::const_iterator const_state_iter;

// ************ Spike ***************
struct Spike
{
    int bin;
    int neuron_ind;
};
// *********************************
// ************ SpikeComparison ***************
class SpikeComparison
{
public:
    bool operator() (const Spike& lhs, const Spike& rhs) const
    {
        return (lhs.bin < rhs.bin);
    }
};
// *********************************

//typedef priority_queue<Spike, vector<Spike>, SpikeComparison> SpikeHeap;

// ************ EMBasins ***************
class paramsStruct;

template <class BasinT>
class EMBasins
{
public:
    EMBasins(int N, int nbasins);
    EMBasins(vector<vector<double> >& st, vector<double>, vector<double>, double binsize, int nbasins);
    ~EMBasins();
    
    vector<double> train(int, bool);
//    vector<double> crossval(int niter, int k);      // k-fold cross-validation
//    vector<double> test(const vector<vector<double> >& st, double binsize);
    
    int nstates() const {return train_states.size();};
    vector<unsigned long> state_hist() const;
//    vector<unsigned long> test_hist() const;
    vector<double> all_prob() const;
//    vector<double> test_prob() const;
    vector<double> P() const;    
    vector<paramsStruct> basin_params();
    vector<char> sample(int);
    vector<char> word_list();
    
    vector<double> w;       // 1 x nbasins
    vector<double> m;       // N x nbasins
    
    //vector<double> test_logli;
    
protected:
    int nbasins, N, T;
    double nsamples;
    RNG* rng;
    vector<BasinT> basins;
    
    vector<string> raster;
    map<string, State> train_states;
    vector<State*> state_list;
    
    vector<Spike> sort_spikes(const vector<vector<double> >&, double) const;
   
    State state_obs(int,int) const;
    vector<double> emiss_obs(int,int) const;
    double logli(bool) const;

private:
    void update_w();

    void update_P();
    
    void initParams();
    void Estep();
    void Mstep();

};

// *********************************

// ************ HMM ***************
template <class BasinT>
class HMM : public EMBasins<BasinT>
{
public:
    HMM(vector<vector<double> >& st, vector<double>, vector<double>, double binsize, int nbasins);
    
    vector<double> train(int, bool);
    vector<int> viterbi(int) const;
    vector<double> emiss_prob() const;
    vector<double> get_forward() const;
    vector<double> get_backward() const;
    vector<double> get_P() const;
    vector<double> get_trans() const;
    vector<double> stationary_prob() const;
    pair<vector<double>, vector<double> > pred_prob() const;
    vector<int> state_v_time() const;
    
    vector<char> sample(int) const;
    vector<double> P_indep() const;
protected:
    vector<double> forward;         // Forward filtering distribution
    vector<double> backward;        // Backward filtering distribution
    
    void update_forward();
    void update_backward();
    void update_P();

    double logli(bool) const;
private:
    vector<double> w0;
    vector<double> trans;           // State transition probability matrix
    
    void update_emiss();
    void update_trans();
    
    static void forward_trans_thread_fun(HMM<BasinT>*);
    static void logli_thread_fun(HMM<BasinT>*,vector<double>*, bool);
    static void backward_thread_fun(HMM<BasinT>*);
    

    
    void initParams();
    void Estep();
    void Mstep();
    
    bool emiss_updated, backward_updated, forward_updated;
    condition_variable cv_emiss, cv_bkwd, cv_fwd;
    mutex emiss_flag_mtx, bkwd_flag_mtx, fwd_flag_mtx;
    

};
// *********************************


// ************ Autocorr ***********
template <class BasinT>
class Autocorr : public HMM<BasinT>
{
public:
    Autocorr(vector<vector<double> >& st, double binsize, int nbasins);
    ~Autocorr();
    
    vector<double> train(int niter);
    vector<int> viterbi();
    vector<double> get_basin_trans();
    
private:
    
    void update_forward();
    void update_backward();
    void update_P();
    void update_w();
    void update_basin_trans_indep();
    
    vector<double> trans_at_t(int);
    
    vector<double*> basin_trans;
    double logli();
};




#endif /* defined(____EMBasins__) */
