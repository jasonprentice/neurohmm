//
//  FitBasinModel.cpp
//  
//
//  Created by Jason Prentice on 11/13/12.
//
//

#include <cmath>
#include <iostream>
#include <vector>

#include "matrix.h"
#include "mex.h"

#include "paramsStruct.h"
#include "HMM.h"
#include "EMBasins.h"

//#define MIXTURE
#define MARKOV

using namespace std;

class TreeBasin;
class IndependentBasin;

typedef TreeBasin BasinType;
const bool ret_train_logli = false;

template <typename T>
void writeOutputMatrix(int pos, vector<T> value, int N, int M, mxArray**& plhs) {
    mxArray* out_matrix = mxCreateDoubleMatrix(N,M,mxREAL);
    double* pr = mxGetPr(out_matrix);
    for (typename vector<T>::iterator it=value.begin(); it!=value.end(); ++it) {
        *pr++ = (double) *it;
    }
    plhs[pos] = out_matrix;
    return;
}

void writeOutputStruct(int pos, vector<paramsStruct>& value, mxArray**& plhs) {
//    for (int n=0; n<value.size(); n++) {
    mxArray* out_struct = mxCreateStructMatrix(value.size(),1,value[0].get_nfields(), value[0].fieldNamesArray());
    int n = 0;
    for (vector<paramsStruct>::iterator it = value.begin(); it != value.end(); ++it) {
        for (int i=0; i < it->get_nfields(); i++) {
            vector<double>* data = it->getFieldData(i);
            
            int N = it->getFieldN(i);
            int M = it->getFieldM(i);
            mxArray* field_data = mxCreateDoubleMatrix(N, M, mxREAL);
            double* pr = mxGetPr(field_data);
            double* iter_pr = pr;
            for (vector<double>::iterator data_it=data->begin(); data_it != data->end(); ++data_it) {
                *iter_pr++ = *data_it;
            }            
            mxSetFieldByNumber(out_struct, n, i, field_data);
        }
        n++;
    }
    plhs[pos] = out_struct;
    return;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
 // [freq,w,m,P,logli,prob] = EMBasins(st, unobserved_edges, binsize, nbasins, niter)

  
    
    cout << "Reading inputs..." << endl;
    int N = mxGetNumberOfElements(prhs[0]);
    vector<vector<double> > st (N);    
    for (int i=0; i<N; i++) {
        mxArray* elem = mxGetCell(prhs[0], i);
        double* elem_pr = mxGetPr(elem);
        int nspikes = mxGetNumberOfElements(elem);
        for (int n=0; n<nspikes; n++) {
            st[i].push_back(elem_pr[n]);
        }
    }
    
    int n_unobserved_blocks = mxGetM(prhs[1]);
    vector<double> unobserved_edges_low (n_unobserved_blocks);
    vector<double> unobserved_edges_high (n_unobserved_blocks);
    if (n_unobserved_blocks > 0 && mxGetN(prhs[1]) != 2) {
        cerr << "Unobserved edges must be empty or have two columns." << endl;
    } else {
        double* unobserved_edges_pr = mxGetPr(prhs[1]);
        for (int n=0; n<n_unobserved_blocks; n++) {
            unobserved_edges_low[n] = unobserved_edges_pr[n];
            unobserved_edges_high[n] = unobserved_edges_pr[n_unobserved_blocks + n];
        }
    }
    
    double binsize = *mxGetPr(prhs[2]);
    int nbasins = (int) *mxGetPr(prhs[3]);
    int niter = (int) *mxGetPr(prhs[4]);

    // ProfilerStart("./EMBasins.prof");
  
#ifdef MARKOV
            // Hidden Markov model
            cout << "Constructing HMM object" << endl;
            HMM<BasinType> basin_obj(st, unobserved_edges_low, unobserved_edges_high, binsize, nbasins);
            vector<double> logli = basin_obj.train(niter, ret_train_logli);

            cout << "Viterbi..." << endl;
            vector<int> alpha = basin_obj.viterbi(2);
            cout << "P...." << endl;
            vector<double> P = basin_obj.get_P();
            cout << "Pred prob..." << endl;
            pair<vector<double>, vector<double> > tmp = basin_obj.pred_prob();
            vector<double> pred_prob = tmp.first;
            vector<double> hist = tmp.second;
        //    vector<unsigned long> hist = basin_obj.state_hist();
            int T = floor(P.size() / nbasins);

            cout << "Params..." << endl;
            vector<paramsStruct> params = basin_obj.basin_params();
            
            writeOutputMatrix(0, logli, niter, 1, plhs);
            writeOutputMatrix(1, basin_obj.get_trans(), nbasins, nbasins, plhs);
            writeOutputMatrix(2, P, nbasins, T, plhs);
//            writeOutputMatrix(2, basin_obj.emiss_prob(), nbasins, T, plhs);
        //    cout << "Microstates..." << endl;
            writeOutputMatrix(3, alpha, T, 1, plhs);
            writeOutputMatrix(4, pred_prob, 1, pred_prob.size(), plhs);
            writeOutputMatrix(5, hist, 1, hist.size(), plhs);
            writeOutputStruct(6, params, plhs);
//            writeOutputMatrix(7, basin_obj.state_v_time(), 1, T, plhs);
            //cout << "Samples..." << endl;
    
            int nsamples = 100000;
            vector<char> sample = basin_obj.sample(nsamples);
            writeOutputMatrix(7, sample, N,nsamples, plhs);
            
            pair<vector<double>, vector<double> > tmp_samp = basin_obj.sample_pred_prob(sample);
            vector<double> samp_prob = tmp_samp.first;
            vector<double> samp_hist = tmp_samp.second;
            writeOutputMatrix(8, samp_prob, 1, samp_prob.size(), plhs);
            writeOutputMatrix(9, samp_hist, 1, samp_prob.size(), plhs);
            
    
           // writeOutputMatrix(7, basin_obj.word_list(), N, hist.size(), plhs);
        //    writeOutputMatrix(6, basin_obj.stationary_prob(), 1,nbasins, plhs);
#endif
#ifdef MIXTURE

            // Mixture model
            cout << "Initializing EM..." << endl;
            EMBasins<BasinType> basin_obj(st, unobserved_edges_low, unobserved_edges_high, binsize, nbasins);
                
            cout << "Training model..." << endl;

            vector<double> logli = basin_obj.train(niter, ret_train_logli);
            //vector<double> test_logli = basin_obj.test_logli;
            
        //    cout << "Testing..." << endl;
        //    vector<double> P_test = basin_obj.test(st_test,binsize);
            
            vector<paramsStruct> params = basin_obj.basin_params();
            int nstates = basin_obj.nstates();
            cout << nstates << " states." << endl;
            
        //    cout << "Getting samples..." << endl;
            int nsamples = 100000;
            vector<char> samples = basin_obj.sample(nsamples);


            cout << "Writing outputs..." << endl;    
            
            writeOutputMatrix(0, basin_obj.w, nbasins,1, plhs);    
            writeOutputStruct(1, params, plhs);
        //    writeOutputMatrix(2, basin_obj.word_list(), N, nstates, plhs);
            writeOutputMatrix(2, samples, N, nsamples, plhs);
            writeOutputMatrix(3, basin_obj.state_hist(), nstates, 1, plhs);
            writeOutputMatrix(4, basin_obj.P(), nbasins, nstates, plhs);
            writeOutputMatrix(5, basin_obj.all_prob(), nstates, 1, plhs);
            writeOutputMatrix(6, logli, niter, 1, plhs);
        //    writeOutputMatrix(6, P_test, nbasins, P_test.size()/nbasins, plhs);
#endif
    
    /*
     // Autocorrelation model
     Autocorr<BasinType> basin_obj(st, binsize, nbasins);
     vector<double> logli = basin_obj.train(niter);
     cout << "Viterbi..." << endl;
     vector<int> alpha = basin_obj.viterbi();
     //    cout << "P...." << endl;
     //    vector<double> P = basin_obj.get_P();
     //    cout << "Pred prob..." << endl;
     //    pair<vector<double>, vector<double> > tmp = basin_obj.pred_prob();
     //    vector<double> pred_prob = tmp.first;
     //    vector<double> hist = tmp.second;
     //    vector<unsigned long> hist = basin_obj.state_hist();
     //    int T = floor(P.size() / nbasins);
     int T = alpha.size();
     
     cout << "Params..." << endl;
     vector<paramsStruct> params = basin_obj.basin_params();
     
     writeOutputMatrix(0, logli, niter, 1, plhs);
     //    writeOutputMatrix(1, basin_obj.get_trans(), nbasins, nbasins, plhs);
     //    writeOutputMatrix(2, P, nbasins, T, plhs);
     //    writeOutputMatrix(2, basin_obj.emiss_prob(), nbasins, T, plhs);
     //    cout << "Microstates..." << endl;
     //    writeOutputMatrix(2, basin_obj.state_v_time(), 1, T, plhs);
     writeOutputMatrix(1, alpha, T, 1, plhs);
     //    writeOutputMatrix(4, pred_prob, 1, pred_prob.size(), plhs);
     //    writeOutputMatrix(5, hist, 1, hist.size(), plhs);
     writeOutputStruct(2, params, plhs);
     writeOutputMatrix(3, basin_obj.get_forward(), nbasins,T,plhs);
     writeOutputMatrix(4, basin_obj.get_backward(), nbasins,T,plhs);
     writeOutputMatrix(5, basin_obj.get_basin_trans(), 4*N,nbasins,plhs);
     writeOutputMatrix(6, basin_obj.w, nbasins,1,plhs);
     writeOutputMatrix(7, basin_obj.P_indep(), nbasins, T, plhs);
     //    cout << "Samples..." << endl;
     //    writeOutputMatrix(3, basin_obj.sample(100000), N,100000, plhs);
     */

    /*
    // k-fold cross-validation
    int kfolds = 10;
    vector<double> logli = basin_obj.crossval(niter, kfolds);
    
    writeOutputMatrix(0, logli, niter, kfolds, plhs);
    */
    
  //  ProfilerStop();
    return;
}


