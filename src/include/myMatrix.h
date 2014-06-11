//
//  myMatrix.h
//  
//
//  Created by Jason Prentice on 6/11/14.
//
//

#ifndef _myMatrix_h
#define _myMatrix_h

#include <vector>
#include <iostream>

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


#endif
