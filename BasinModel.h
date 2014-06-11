//
//  BasinModel.h
//  
//
//  Created by Jason Prentice on 11/22/12.
//
//

#ifndef ____BasinModel__
#define ____BasinModel__

#include <iostream>
#include <vector>
#include <map>
#include <string>


struct State;       // Defined in EMBasins.h
class RNG;

// *********************** myMatrix ****************************
template <class T>
class myMatrix {
public:
                            myMatrix();
                            myMatrix( std::vector<T>&, int, int );
    void                    assign( std::vector<T>&, int, int );
    int                     get_N();
    int                     get_M();
    const T&                at( const int, const int ) const;     // Subscripted index
    const T&                at( const int ) const;            // Linear index
    std::vector<T> *        data();
private:
    std::vector<T>          matrix_data;
    int                     N, M;        // N = #rows, M = #columns
};

// myMatrix definition

template <class T>
myMatrix<T>::myMatrix() : N(0),M(0) {};

template <class T>
myMatrix<T>::myMatrix(std::vector<T>& _data, int _N, int _M) : N(_N), M(_M) {
    if (_data.size() != _N*_M) {
        std::cerr << "Matrix dimensions must agree." << std::endl;
        N = 0; M = 0;
    } else {
        matrix_data = _data;
    }
}

template <class T>
void myMatrix<T>::assign(std::vector<T>& _data, int _N, int _M) {
    if (_data.size() != _N*_M) {
        std::cerr << "Matrix dimensions must agree." << std::endl;
        N = 0; M = 0;
    } else {
        N = _N; M = _M;
        matrix_data = _data;
    }
    return;
}

template <class T>
int myMatrix<T>::get_N() { return N; }

template <class T>
int myMatrix<T>::get_M() { return M; }

template <class T>
const T& myMatrix<T>::at(const int i, const int j) const {
//    if (i < 0 || i >= N || j < 0 || j >= M) {
//        std::cerr << "Index exceeds matrix dimensions." << std::endl;
//    } else {
        return matrix_data[j*N + i];
//    }
}
template <class T>
const T& myMatrix<T>::at(const int i) const {
//    if (i<0 || i>= N*M) {
//        std::cerr << "Index exceeds matrix dimensions." << std::endl;
//    } else {
    
        return matrix_data[i];
//    }
}

template <class T>
std::vector<T>* myMatrix<T>::data() { return &matrix_data; }
// *****************************************************************

// *********************** paramsStruct ****************************
class paramsStruct {            // Variable-sized structure holding double matrices
public:
                            paramsStruct();
    int                     get_nfields();
    void                    addField( std::string, myMatrix<double>& );
    const char**            fieldNamesArray();
    std::vector<double>*         getFieldData( int );
    int                     getFieldN( int );
    int                     getFieldM( int );
    const char*             getFieldName( int );
    
private:
    int                     nfields;
    std::map<std::string,
        myMatrix<double>* > fields;
    std::vector<const char*> fieldNames;
    
};
// *****************************************************************

// ************************** BasinModel ************************
class BasinModel
{
public:
                            BasinModel( int N, int basin_num, RNG* rng )
                                : N(N), basin_num(basin_num), rng(rng) {};
    void                    reset_stats();
    void                    increment_stats( const State& );
    void                    normalize_stats();
    inline double           get_norm() const {
                                return norm; };
protected:
    std::vector<double>     stats;
    double                  norm;
    int                     N;
    int                     basin_num;
    RNG *                   rng;
};
// ***************************************************************

// ************************** IndependentBasin ************************

class IndependentBasin : public BasinModel {
public:
                            IndependentBasin( int, int, RNG*, double );
    static std::vector<int> get_active_constraints( const State& );
    void                    doMLE( double );
    double                  P_state( const State& ) const;
    std::vector<char>       sample() const;
    paramsStruct            get_params();
private:
    myMatrix<double>        m;
    std::vector<char>       above_thresh_bool;
    std::vector<int>        above_thresh_list;
    double                  prefactor;
    
    void                    update_thresh_list();
};
// ***************************************************************


#endif /* defined(____BasinModel__) */
