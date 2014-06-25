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
#include "SpikeData.h"

class paramsStruct;
class BasinModel;
class RNG;

template <class BasinT>
class EMBasins {
public:
    struct BasinData {
        std::vector<double>          P;
        std::vector<double>          weight;
        std::vector<int>             active_constraints;
        double                       pred_prob;
        
        BasinData( const State<BasinData>& this_state )
        : active_constraints (BasinT::get_active_constraints(this_state)), pred_prob (-1) {};
    };

                                EMBasins( int, int );
                                EMBasins( const std::vector<std::vector<double> >&, std::vector<double>, std::vector<double>, double, int );
                                ~EMBasins();
    std::vector<double>          train_model( int, bool );
    inline std::vector<double>   get_w() const
                                        { return w; };
/*    inline std::vector<double>   get_m() const
                                        { return m; }; */
    std::vector<double>          all_prob( bindesc_t ) const;
    std::vector<double>          P( bindesc_t ) const;
    std::vector<paramsStruct>    basin_params() const;
    std::vector<char>            sample( int );
protected:

//    friend BasinModel;
    
    RNG *                   rng;
    std::vector<BasinT>     basins;
    SpikeData<BasinData>    data;
    int                     N, T, nbasins;

/*    State<BasinData>        state_obs( int, int ) const;
    std::vector<double>     emiss_obs( int, int ) const;*/
    double                  logli( bindesc_t ) const;
private:
    std::vector<double>          w;       // 1 x nbasins
//    std::vector<double>          m;       // N x nbasins

    
                            EMBasins( const EMBasins& );
    EMBasins&               operator=( const EMBasins& );

    void                    update_w();
    void                    update_P( bindesc_t );
    void                    initParams();
    void                    Estep();
    void                    Mstep();
};


#include "EMBasins.cpp"




#endif /* defined(____EMBasins__) */
