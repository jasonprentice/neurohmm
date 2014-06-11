//
//  HMM.h
//  
//
//  Created by Jason Prentice on 6/11/14.
//
//

#ifndef ____HMM__
#define ____HMM__

#include <vector>
#include <string>
#include <map>
#include <mutex>
#include <condition_variable>

#include "EMBasins.h"

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

#include "HMM.cpp"

#endif /* defined(____HMM__) */
