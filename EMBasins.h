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


// ************ RNG ***************
class RNG {
public:
                            RNG();
                            ~RNG();
    int                     discrete( const std::vector<double>& ) const;
    bool                    bernoulli( double ) const;
    double                  uniform() const;
    std::vector<int>             randperm( int ) const;
private:
    gsl_rng *               rng_pr;
};

struct State {
    std::vector<double>          P;
    std::vector<double>          weight;
    std::vector<int>             active_constraints;
    std::vector<int>             on_neurons;
    double                  freq;
    double                  pred_prob;
    std::vector<char>            word;
    int                     identifier;
    
    
                            State() :
                                freq(0), pred_prob(-1), identifier(-1) {};
    
    
};
// *********************************
typedef std::map<std::string,State>::iterator state_iter;
typedef std::map<std::string,State>::const_iterator const_state_iter;

// ************ Spike ***************
struct Spike {
    int                     bin;
    int                     neuron_ind;
};
// *********************************
// ************ SpikeComparison ***************
class SpikeComparison {
public:
    bool                    operator() ( const Spike& lhs, const Spike& rhs ) const {
        return (lhs.bin < rhs.bin);
    }
};
// *********************************

//typedef priority_queue<Spike, std::vector<Spike>, SpikeComparison> SpikeHeap;


// ************ EMBasins ***************
class paramsStruct;

template <class BasinT>
class EMBasins {
public:
    std::vector<double>          w;       // 1 x nbasins
    std::vector<double>          m;       // N x nbasins
    
                            EMBasins( int, int );
                            EMBasins( std::vector<std::vector<double> >&, std::vector<double>, std::vector<double>, double, int );
                            ~EMBasins();
    std::vector<double>          train( int, bool );
    std::vector<char>            word_list();
    std::vector<unsigned long>   state_hist() const;
    std::vector<double>          all_prob() const;
    std::vector<double>          P() const;
    std::vector<paramsStruct>    basin_params();
    std::vector<char>            sample( int );
    int                     nstates() const {
                                return train_states.size(); };
protected:
    int                     N, nbasins, T;
    double                  nsamples;
    RNG *                   rng;
    std::vector<BasinT>          basins;
    
    std::vector<std::string>          raster;
    std::map<std::string, State>      train_states;
    std::vector<State*>          state_list;
    
    std::vector<Spike>           sort_spikes( const std::vector<std::vector<double> >&, double ) const;
    State                   build_state( std::vector<char> ) const;
    State                   state_obs( int, int ) const;
    std::vector<double>          emiss_obs( int, int ) const;
    double                  logli( bool ) const;

private:
    void                    update_w();
    void                    update_P();
    void                    initParams();
    void                    Estep();
    void                    Mstep();

};

// *********************************

// ************ HMM ***************
template <class BasinT>
class HMM : public EMBasins<BasinT> {
public:
                            HMM( std::vector<std::vector<double> >&, std::vector<double>, std::vector<double>, double, int );
    
    std::vector<double>          train( int, bool );
    std::vector<int>             viterbi( int ) const;
    std::vector<double>          emiss_prob() const;
    std::vector<double>          get_forward() const;
    std::vector<double>          get_backward() const;
    std::vector<double>          get_P() const;
    std::vector<double>          get_trans() const;
    std::vector<double>          stationary_prob() const;
    std::pair<std::vector<double>,
         std::vector<double> >   pred_prob() const;
    std::pair<std::vector<double>,
         std::vector<double> >   sample_pred_prob( const std::vector<char>& ) const;
    std::vector<int>             state_v_time() const;
    
    std::vector<char>            sample( int ) const;
    std::vector<double>          P_indep() const;
protected:
    std::vector<double>          forward;         // Forward filtering distribution
    std::vector<double>          backward;        // Backward filtering distribution
    
    void                    update_forward();
    void                    update_backward();
    void                    update_P();

    double                  logli( bool ) const;
private:
    std::vector<double>     trans;           // State transition probability matrix
    std::vector<double>     w0;
    bool                    emiss_updated,
                            backward_updated,
                            forward_updated;
    std::condition_variable cv_emiss,
                            cv_bkwd,
                            cv_fwd;
    std::mutex              emiss_flag_mtx,
                            bkwd_flag_mtx,
                            fwd_flag_mtx;

    static void             forward_trans_thread_fun( HMM<BasinT> * );
    static void             logli_thread_fun( HMM<BasinT> *, std::vector<double> *, bool );
    static void             backward_thread_fun( HMM<BasinT> * );
    std::pair<std::vector<double>,
         std::vector<double> >   pred_prob_helper( const std::map<std::string,State>& ) const;
    void                    update_emiss();
    void                    update_trans();
    void                    initParams();
    void                    Estep();
    void                    Mstep();
    
    

};
// *********************************


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




#endif /* defined(____EMBasins__) */
