//
//  SpikeData.h
//  
//
//  Created by Jason Prentice on 6/13/14.
//
//

#ifndef ____SpikeData__
#define ____SpikeData__

#include <vector>
#include <map>

template <typename T>
class State {
public:
    int                          freq;
    T                            aux_data;
    
                                 State()
                                    : freq(0), identifier(-1) {};
                                 State( const std::vector<char>& );
    
    inline std::vector<int>      get_on_neurons() const { return on_neurons; };
    inline std::vector<char>     get_word() const { return word; };
    inline int                   get_identifier() const { return identifier; };
private:
    std::vector<int>             on_neurons;
    std::vector<char>            word;
    int                          identifier;
    
    State&                       operator=( const State& );
};


enum bindesc_t { train, test, unobserved };


template <typename U>
class SpikeData {
public:
    typedef std::pair<bindesc_t, const State<U>* > bintag_t;
    typedef typename std::vector<State<U> >::iterator state_iter_t;
    typedef typename std::vector<State<U> >::const_iterator const_state_iter_t;
    
                                 SpikeData( std::vector<std::vector<double> >&, std::vector<double>, std::vector<double>, double );

    inline state_iter_t          begin( bindesc_t desc ) {
                                    return states[desc].begin(); };
    inline const_state_iter_t    begin( bindesc_t desc ) const {
                                    return states[desc].begin(); };
    inline state_iter_t          end( bindesc_t desc ) {
                                    return states[desc].end(); };
    inline const_state_iter_t    end( bindesc_t desc ) const {
                                    return states[desc].end(); };
    
    inline int                   get_T() const { return T; };
    inline int                   get_N() const { return N; };
    std::vector<char>            word_list( bindesc_t ) const;
    std::vector<int>             state_hist( bindesc_t ) const;
    inline int                   nstates( bindesc_t desc ) const {
                                    return states[desc].size(); };
private:
    typedef std::map<bindesc_t, std::vector<State<U> > > state_map_t;
    
    int                          N, T;
    state_map_t                  states;
    std::vector<bintag_t>        state_list;
    
    typedef std::pair<bindesc_t, std::string> key_t;
    void                         insert_word(int, const std::vector<double>&, const std::vector<double>&,
                                             std::string, const std::vector<char>&,
                                             std::map<key_t, State<U> >&, std::vector<std::string>&);
    struct Spike {
        int                     bin;
        int                     neuron_ind;
    };
    struct SpikeComparison {
        inline bool             operator() ( const Spike& lhs, const Spike& rhs ) const {
            return (lhs.bin < rhs.bin); };
    };
    std::vector<Spike>           sort_spikes( const std::vector<std::vector<double> >&, double ) const;
};


#endif /* defined(____SpikeData__) */
