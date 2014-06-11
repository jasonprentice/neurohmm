//
//  EMBasins.h
//  
//
//  Created by Jason Prentice on 11/13/12.
//
//

#ifndef ____EMBasins__
#define ____EMBasins__


#include <vector>
#include <string>
#include <map>


struct State {
    std::vector<double>          P;
    std::vector<double>          weight;
    std::vector<int>             active_constraints;
    std::vector<int>             on_neurons;
    double                       freq;
    double                       pred_prob;
    std::vector<char>            word;
    int                          identifier;
    
                                State()
                                    : freq(0), pred_prob(-1), identifier(-1) {};
    
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
class RNG;

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

#include "EMBasins.cpp"




#endif /* defined(____EMBasins__) */
