//
//  SpikeData.cpp
//  
//
//  Created by Jason Prentice on 6/13/14.
//
//

#include <algorithm>

template <typename T>
State<T>::State( const std::vector<char>& word ) : word (word), freq (0), identifier (-1) {
    for (int i=0; i<word.size(); i++) {
        if (word[i] == 1) {
            on_neurons.push_back(i);
        }
    }
    aux_data = T(this);
}

template <typename U>
void SpikeData<U>::insert_word(const std::vector<char>& this_word,
                               std::string this_str, int curr_bin,
                               const std::vector<double>& unobserved_l,
                               const std::vector<double>& unobserved_u,
                               std::map<key_t, int>& state_index
                               std::vector<key_t>& keys) {
    
    bindesc_t desc = train;
    for (int n=0; n<unobserved_l.size(); n++) {
        if (curr_bin >= unobserved_l[n] && curr_bin < unobserved_u[n]) {
            desc = test;
            break;
        }
    }
    
    key_t key (desc, this_str);
    keys[curr_bin] = key;
    
    auto elem = std::pair<key_t,int> (key, states[desc].size());
    auto ins = state_index.insert( elem );
    if ( !ins.second ) {
        int index = ins.first->second;
        states[desc][index].freq++;
    } else {
        states[desc].emplace_back(this_word);
    }
    
    return;
}

template <typename U>
SpikeData<U>::SpikeData(const std::vector<std::vector<double> >& st, std::vector<double> unobserved_l, std::vector<double> unobserved_u, double binsize ) : N (st.size()) {
    
    std::transform(unobserved_l.begin(), unobserved_l.end(), unobserved_l.begin(), [binsize] (t) {
        return std::floor(t/binsize); });
    std::transform(unobserved_u.begin(), unobserved_u.end(), unobserved_u.begin(), [binsize] (t) {
        return std::floor(t/binsize); });
    
    // Build state structure from spike times in st:
    cout << "Building state histogram..." << endl;
    
    std::vector<Spike> all_spikes = sort_spikes(st, binsize);
    T = all_spikes.back().bin;
    
    std::vector<key_t> keys (T);
    std::map<key_t, int> state_index;
    
    std::string silent_str (N,'0');

    std::string this_str (silent_str);
    vector<char> this_word (N,0);
    int curr_bin = 0;
    for (auto it=all_spikes.begin(); it!=all_spikes.end(); ++it) {
        int next_bin = it->bin;
        int next_cell = it->neuron_ind;
        
        if (next_bin > curr_bin) {
            insert_word(this_word, this_str, curr_bin, unobserved_l, unobserved_u, state_index, keys);
            this_str = silent_str;
            this_word.assign(N,0);
            // All states between curr_bin and next_bin (exclusive) are silent
            for (int t=curr_bin+1; t<next_bin; t++) {
                insert_word(this_word, silent_str, t, unobserved_l, unobserved_u, state_index, keys);
            }
            curr_bin = next_bin;
        }
        // Add next_cell to this_state
        this_word[next_cell] = 1;
        this_str[next_cell] = '1';
    }
    
    int identifier = 0;
    for (auto it = state_index.begin(); it != state_index.end(); ++it) {
        int index = it->second;
        bindesc_t desc = (it->first).first;
        states[desc][index].identifier = identifier;
        identifier++;
    }
    
    for (int t=0; t<T; t++) {
        binddesc_t desc = keys[t].first;
        int index = state_index[keys[t]];
        state_list.emplace_back(desc, &(states[desc].at(index)));
    }
}


template <typename U>
vector<unsigned long> SpikeData<U>::state_hist( bindesc_t desc ) const {
    vector<unsigned long> hist (states[desc].size(), 0);
    int pos = 0;
    for (auto it=states[desc].begin(); it != states[desc].end(); ++it) {
        const State& this_state = it->second;
        hist[pos] = this_state.freq;
        pos++;
    }
    return hist;
}


template <typename U>
vector<char> SpikeData<U>::word_list( bindesc_t desc ) {
    vector<char> out (states[desc].size() * N);
    vector<char>::iterator out_it = out.begin();
    for (auto it = states[desc].begin(); it != states[desc].end(); ++it) {
        vector<char> word = (it->second).get_word();
        for (vector<char>::iterator w_it = word.begin(); w_it!=word.end(); ++w_it) {
            *out_it++ = *w_it;
        }
    }
    return out;
}




